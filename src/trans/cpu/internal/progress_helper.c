#define _GNU_SOURCE
#include <stdio.h>
#include <mpi.h>
#include <pthread.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <sched.h>
#include <unistd.h>
#include <stdbool.h>
#include <stdint.h>
#include <assert.h>
#include <stdatomic.h>
#include <sys/syscall.h>
#include "progress_helper.h"

#ifdef SYS_gettid
#define gettid() ((pid_t)syscall(SYS_gettid))
#else
#error "SYS_gettid unavailable on this system"
#endif

/* User configuration */
static bool verbose = false;  /* be quiet or verbose */
int read_counters = 1;  /* Non-zero to activate the counters */
char* storage_array_filename_template = "/tmp/app.rank";  /* Template for the counters filename */
double read_interval = 1e-1;  /* interval to check the counters 100 ms */
int32_t report_all = 0;  /* report all values, aka. ignore the cutoff_amount */
int32_t cutoff_amount = 4*1024;  /* if sum of counters below this threshold do not report */
int32_t storage_array_length = 1024*1024;  /* the size of the counter's storage array. Bigger
                                            * has less frequent flush but uses more memory. */

#define PROGRESS_THREAD_ENV "PROGRESS_THREAD_ENV"
#define PROGRESS_THREAD_COUNTERS "PROGRESS_THREAD_COUNTERS"
#define PROGRESS_THREAD_FILE "PROGRESS_THREAD_FILE"
#define MAX_REQSET_SIZE 20

typedef struct files_to_read_s {
    char* name;
    FILE* fp;
    int64_t last_read;
} files_to_read_t;

files_to_read_t files_to_read[] = {
    { .name = "/sys/class/infiniband/mlx0_1/ports/1/counters/port_xmit_data", .fp = NULL, .last_read = -1 },
    { .name = "/sys/class/infiniband/mlx0_1/ports/1/counters/port_rcv_data", .fp = NULL, .last_read = -1 },
};
int files_to_read_size = 0;

/* NO NEED TO CHANGE ANYTHING BELOW THIS */

typedef struct pt_reqset_s {
    struct pt_reqset_s *prev;
    struct pt_reqset_s *next;
    int                 count;
    int                 flags;
    MPI_Request*        array_of_requests;
    MPI_Status*         array_of_statuses;
    /* for management of independently completed requests */
    int                 reported_completion;
    __sig_atomic_t                 detected_completion;
    int*                idx_completed_reqs;
} pt_reqset_t;

typedef struct mpi_helper_s {
    MPI_Comm            comm;
    MPI_Request         sync_req;
    pthread_t           thread;
    pthread_mutex_t     cq_mutex;
    pt_reqset_t        *cq;
    int                 recvint;
    int                 read_files;
    int                 nothing_happened;  /* keep track of how active theprogress thread is:
                                            * 0 is active, as the value increases the thread will
                                            * start yielding resources via nano_sleep. This value
                                            * will be capped by max_slack.
                                            */
    int                 max_slack;          /* when the thread is not busy how long could it go to sleep */
    int                 be_quiet;           /* temporarily force the progress thread to almost stop */
    pthread_mutex_t     quiet_mutex;        /* block on this if quiet mode is requested */
} mpi_helper_t;

static mpi_helper_t* mpi_helper = NULL;
static double startup_wtime;
static int32_t* storage_array = NULL;
static int32_t storage_array_count = 0;
FILE* storage_array_fp = NULL;

#define LINELEN 32

static int flush_storage_array(void)
{
    printf("save %d counters\n", storage_array_count);
    fwrite(storage_array, storage_array_count, sizeof(int32_t), storage_array_fp);
    storage_array_count = 0;  /* start from the begining */
    return 0;
}

static int open_files(files_to_read_t* ftr, int count)
{
    char line[LINELEN];
    int cnt = 0;
    for(int i = 0; i < count; i++) {
        assert(NULL == ftr[i].fp);
        ftr[i].fp = fopen(ftr[i].name, "r");
        if(NULL == ftr[i].fp) {
            /* this file cannot be opened, skip it */
            continue;
        }
        cnt++;
        
        /* Get some sane starting points */
        (void)fread(line, 1, LINELEN, ftr[i].fp);
        ftr[i].last_read = atol(line);     /* get the new counter */
    }
    return cnt;
}

static int close_files(files_to_read_t* ftr, int count)
{
    int cnt = 0;
    for(int i = 0; i < count; i++) {
        if(NULL == ftr[i].fp) continue;
        fclose(ftr[i].fp);
        ftr[i].fp = NULL;
        cnt++;
    }
    return cnt;
}

/**
 * Read the counter files and report the values if above a predefined threshold.
 * Return the number of values reported, or 0 if their sum was below the threshold.
 */
static int read_files(files_to_read_t* ftr, int count)
{
    int32_t cnt = 0, *amounts, sum = 0;
    char line[LINELEN];
    int64_t last_value;
    size_t len;
    
    double elapsed = MPI_Wtime() - startup_wtime;

    if((storage_array_count+1+count) > storage_array_length)
        flush_storage_array();
    storage_array[storage_array_count] = (int)(elapsed * 1e6);
    amounts = &storage_array[storage_array_count+1];
    for(int i = 0; i < count; i++) {
        if(NULL == ftr[i].fp) {
            amounts[i] = 0;
            continue;
        }
        fseek(ftr[i].fp, 0, SEEK_SET);  /* rollback to begining */
        len = fread(line, 1, LINELEN, ftr[i].fp);
        if(0 >= len) {
            printf("read incorrect data %ld\n", len);
            amounts[i] = 0;
            continue;
        }
        last_value = ftr[i].last_read;                    /* save the old value */
        ftr[i].last_read = atol(line);                    /* get the new counter */
        amounts[i] = (int)ftr[i].last_read - last_value;  /* compute the diff */
        sum += amounts[i];
        cnt++;
    }

    if(report_all || (sum > cutoff_amount)) {
        storage_array_count += (count + 1);  /* move onto the next section */
        return (count + 1);
    }
    return 0;
}

void* mpi_helper_thread_routine( void* args )
{
    mpi_helper_t* helper = (mpi_helper_t*)args;
    cpu_set_t *binding_mask = NULL;
    int flag, rank;
    double last = MPI_Wtime(), now;
    pt_reqset_t* active_rsets = NULL;

    startup_wtime = last;
    MPI_Comm_rank(helper->comm, &rank);

    /* We want to bind the progress thread away from all the compute threads.
     * Check if the environment variable was provided, and if yes follow the required
     * binding
     */
    if( NULL != getenv(PROGRESS_THREAD_ENV) ) {
        int cpusetsize = 256, rank, idx, i, setsize;
        cpu_set_t *index_mask = CPU_ALLOC(cpusetsize);
        binding_mask = CPU_ALLOC(cpusetsize);
        setsize = CPU_ALLOC_SIZE(cpusetsize);
        MPI_Comm_rank(mpi_helper->comm, &rank);

        char *binding_env = strdup(getenv(PROGRESS_THREAD_ENV)),  /* save a copy to alter */
             *saved_base = binding_env;
        /* Mark all the allowed cores */
        CPU_ZERO_S(setsize, index_mask);
        while(NULL != binding_env) {
            char* token = strsep(&binding_env, ",");
            int idx = atoi(token);
            CPU_SET_S(idx, setsize, index_mask);
            if(verbose)
                printf("rank %d: allow %d-indexed core into the progress binding mask\n", rank, idx);
        }
        free(saved_base);

        /* Get this process binding and trim it for the orogress thread */
        sched_getaffinity(0  /* this process */, setsize, binding_mask);
        for( i = idx = 0; i < cpusetsize; i++ ) {
            if( CPU_ISSET_S(i, setsize, binding_mask) ) {
                if( !CPU_ISSET_S(idx, setsize, index_mask) ) {
                    CPU_CLR_S(i, setsize, binding_mask);
                }
                idx++;
            }
        }
        if(verbose) {
            for( i = 0; i < cpusetsize; i++ ) {
                if( CPU_ISSET_S(i, setsize, binding_mask) ) {
                    printf("rank %d: progress thread will be bound to core %d\n", rank, i);
                }
            }
        }
        CPU_FREE(index_mask);
        /* And now bind thyself */
        pid_t tid = gettid();
        sched_setaffinity(tid, setsize, binding_mask);
    }

    if(helper->read_files) {
        /* And now prepare the output file */
        char *filename_env = getenv(PROGRESS_THREAD_FILE);
        char *filename;
        if( NULL != strchr(filename_env, '%') ) {
            fprintf(stderr, "The PROGRESS_THREAD_FILE environment cannot contain %% on rank %d\n", rank);
            helper->read_files = 0;
        } else {
	        asprintf(&filename, "%s.rank%d",
                     (NULL == filename_env) ? storage_array_filename_template : filename_env, rank);
            storage_array_fp = fopen(filename, "w");
            if(NULL == storage_array_fp) {
                helper->read_files = 0;
	            fprintf( stderr, "Cannot open the output file %s. Bail out!\n", filename);
            } else {
                storage_array = (int32_t*)malloc(storage_array_length*sizeof(int32_t));
                storage_array_count = 0;

                /* one rank per node is reading the IB counter */
                (void)open_files(files_to_read, sizeof(files_to_read)/sizeof(files_to_read_t));
            }
            free(filename);
        }
    }

    do {
        pt_reqset_t* cmd;
        
        if( NULL != helper->cq ) {  /* absorb one new reqsets from the application*/
            pthread_mutex_lock(&helper->cq_mutex);
            cmd = helper->cq;
            if( cmd->next == cmd ) {  /* single command */
                helper->cq = NULL;
            } else {  /* pick the first one, put the rest back into the cq */
                cmd->next->prev = cmd->prev;
                cmd->prev->next = cmd->next;
                helper->cq = cmd->next;
                cmd->next = cmd;
                cmd->prev = cmd;
            }
            pthread_mutex_unlock(&helper->cq_mutex);

            /* Only persistent requests are started */
            if( !(cmd->flags & PT_REQSET_NON_PERSISTENT) )
                MPI_Startall(cmd->count, cmd->array_of_requests);

            if( NULL == active_rsets ) {
                active_rsets = cmd;
            } else {
                cmd->next = active_rsets;
                cmd->prev = active_rsets->prev;
                active_rsets->prev->next = cmd;
                active_rsets->prev = cmd;
            }
            assert(cmd->flags & PT_REQSET_ACTIVE);
            helper->nothing_happened = 0;  /* we now have work to do */
        }
        /* Check the status of active reqsets */
        if( NULL != (cmd = active_rsets) ) {
            do {
                pt_reqset_t *tmpitem = cmd->next;
                if (cmd->flags & PT_REQSET_CHECK_ANY) {  /* report independent completion */
                    int outcount = 0;
                    MPI_Testsome(cmd->count, cmd->array_of_requests, &outcount,
                                 &cmd->idx_completed_reqs[cmd->detected_completion],
                                 (NULL == cmd->array_of_statuses) ? MPI_STATUSES_IGNORE : &cmd->array_of_statuses[cmd->detected_completion]);
                    if( 0 != outcount ) {
                        atomic_thread_fence(memory_order_release);
                        atomic_store_explicit((_Atomic int*) &cmd->detected_completion, cmd->detected_completion + outcount, memory_order_relaxed);
                    }
                    /* Everything done ? */
                    flag = (cmd->detected_completion == cmd->count);  
                } else {  /* report all completions once*/
                    MPI_Testall(cmd->count, cmd->array_of_requests, &flag,
                                (NULL == cmd->array_of_statuses) ? MPI_STATUSES_IGNORE : cmd->array_of_statuses);
                    if( flag ) {
                        atomic_thread_fence(memory_order_release);
                        atomic_store_explicit((_Atomic int*) &cmd->detected_completion, cmd->count, memory_order_relaxed);
                    }
                }
                if( flag ) {
                    /* This active reqset is completed. Mark it and go to the next */
                    if( cmd->next == cmd ) {  /* single active reqset */
                        assert(cmd == active_rsets);
                        active_rsets = NULL;
                        tmpitem = NULL;
                    } else {
                        cmd->prev->next = cmd->next;
                        cmd->next->prev = cmd->prev;
                        if( cmd == active_rsets ) {
                            active_rsets = cmd->next;
                        }
                    }
                    cmd->next = NULL;
                    cmd->prev = NULL;
                    /* No need for protection here, nobody else should alter the reqset flags. The current flags
                     * should be active, so using ^should turn off active and turn on completed. */
                    cmd->flags ^= (PT_REQSET_ACTIVE | PT_REQSET_COMPLETED);
                    assert(PT_REQSET_COMPLETED == (cmd->flags & ((PT_REQSET_ACTIVE | PT_REQSET_COMPLETED))));
                    /* drop the current [completed] cmd */
                }  /* otherwise move onto the next cmd */
                cmd = tmpitem;
            } while (cmd != active_rsets);
            flag = 0;  /* make sure we dont allow the progress thread to completeand quit */
        } else {
            /* No MPI progress this iteration so let's force MPI to do something */
            MPI_Test(&helper->sync_req, &flag, MPI_STATUS_IGNORE);
            helper->nothing_happened++;
            if( helper->nothing_happened > 50 ) {
                /* Try to slow down if there is nothing to do */
                struct timespec ts = {.tv_nsec = 1000 * helper->nothing_happened};
                if( ts.tv_nsec > helper->max_slack ) ts.tv_nsec = helper->max_slack;
                nanosleep(&ts, NULL);
            }
        }
        if(helper->read_files) {
            now = MPI_Wtime();
            if( (now - last) < read_interval )
                continue;
            read_files(files_to_read, sizeof(files_to_read)/sizeof(files_to_read_t));
            last = now;
        }
        if( helper->be_quiet ) { /* we were asked to pause */
            /* lock the mutex. As the main thread will own the mutex in quiet mode, this
             * thread will go to sleep and will only be awaken when the main thread will
             * release the mutex, aka in unpause. At that moment is should release the 
             * newly acquired mutex and reset the be_quiet to 0.
             */
            pthread_mutex_lock(&helper->quiet_mutex);
            pthread_mutex_unlock(&helper->quiet_mutex);
            helper->be_quiet = 0;
        }
    } while( 0 == flag );

    if(helper->read_files) {
        /* Flush and close the output file */
        flush_storage_array();
        fclose(storage_array_fp);storage_array_fp = NULL;
        
        (void)close_files(files_to_read, sizeof(files_to_read)/sizeof(files_to_read_t));
        free(storage_array);
        storage_array = NULL;
    }
    CPU_FREE(binding_mask);

    return NULL;
}

int start_MPI_helper(void)
{
    int rc, my_rank;
    MPI_Comm node_comm;

    mpi_helper = (mpi_helper_t*)malloc(sizeof(mpi_helper_t));
    rc = MPI_Comm_dup(MPI_COMM_WORLD, &mpi_helper->comm);
    
    mpi_helper->read_files = 0;  /* default: not reading the counters */

    if( read_counters ) {
        if (NULL != getenv(PROGRESS_THREAD_FILE) ) {
            int my_node_rank;
            MPI_Comm_split_type(mpi_helper->comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &node_comm);
            MPI_Comm_rank(node_comm, &my_node_rank);
            mpi_helper->read_files = (my_node_rank == 0);
            MPI_Comm_free(&node_comm);
        }
    }
    
    rc = MPI_Comm_rank(mpi_helper->comm, &my_rank);
    rc = MPI_Irecv(&mpi_helper->recvint, 1, MPI_INT, my_rank, 0, mpi_helper->comm,
                   &mpi_helper->sync_req);

    mpi_helper->cq = NULL;
    mpi_helper->nothing_happened = 0;
    mpi_helper->max_slack = 100000;  /* in nano-seconds */
    mpi_helper->be_quiet = 0;  /* go do your job */
    pthread_mutex_init(&mpi_helper->cq_mutex, NULL);
    pthread_mutex_init(&mpi_helper->quiet_mutex, NULL);
    rc = pthread_create(&mpi_helper->thread, NULL, mpi_helper_thread_routine, mpi_helper);
    return rc;
}

int pause_MPI_helper(void)
{
    int rc = pthread_mutex_trylock(&mpi_helper->quiet_mutex);
    if(0 != rc ) {
        fprintf(stderr, "The quiescence mutex is already blocked. Ignoring the pause command\n");
        return -1;
    }
    mpi_helper->be_quiet = 1;
    return 0;
}

int unpause_MPI_helper(int wait)
{
    pthread_mutex_unlock(&mpi_helper->quiet_mutex);
    /* This will eventually release the progress thread, who will then set the be_quiet back to 0
     * to signal that he is back to work.
     */
    if( wait ) {
        while( mpi_helper->be_quiet ) {
            struct timespec ts = {.tv_nsec = 1000};
            /* do something */
            nanosleep(&ts, NULL);
        }
    }
    return mpi_helper->be_quiet;
}


int stop_MPI_helper(void)
{
    int rc = 0, my_rank;
    void* ret;

    rc = MPI_Comm_rank(mpi_helper->comm, &my_rank);
    if (1 == mpi_helper->be_quiet ) {  /* the progress thread is quiet, wake it up. */
        unpause_MPI_helper(1 /* wait until awake */);
    }
    /* send the message to complete the helper thread */
    rc = MPI_Send(&rc, 1, MPI_INT, my_rank, 0, mpi_helper->comm);
    /* wait until the helper thread completes */
    pthread_join(mpi_helper->thread, &ret);

    assert(MPI_REQUEST_NULL == mpi_helper->sync_req);

    MPI_Comm_free(&mpi_helper->comm);
    free(mpi_helper);
    mpi_helper = NULL;

    return 0;
}

static int pt_reqset_array_active = 0;
static pt_reqset_t pt_reqset_array[MAX_REQSET_SIZE];

int pt_reqset_register(int count, MPI_Request* array_of_requests, int flags, int* gid)
{
    pt_reqset_t* rset;

    /* Find a storage place for the reqset */
    *gid = -1;  /* insane default value */
    if( pt_reqset_array_active == MAX_REQSET_SIZE ) {  /* no more room */
        fprintf(stderr, "Maximum number of entries (%d) in the reqset array reached! Recompile with a larger number\n",
                MAX_REQSET_SIZE);
        return MPI_ERR_NO_SPACE;
    }
    for( int idx = 0; idx < MAX_REQSET_SIZE; idx++ ) {
        if( NULL == pt_reqset_array[idx].array_of_requests ) {
            *gid = idx;
            break;
        }
    }
    if( *gid == -1 ) {
        fprintf(stderr, "Unable to find a location for the new requests set. Bailing out !");
        return MPI_ERR_ARG;
    }
    rset = &pt_reqset_array[*gid];
    rset->count = count;
    rset->flags = flags;
    rset->array_of_requests = array_of_requests;  /* will be overwritten if not allowed to use the user array */
    rset->idx_completed_reqs = NULL;
    if( !(flags & PT_REQSET_USE_REQUEST_ARRAY) ) {
        rset->array_of_requests = (MPI_Request*)malloc(count * sizeof(MPI_Request));
        memcpy(rset->array_of_requests, array_of_requests, count * sizeof(MPI_Request));
    }
    rset->array_of_statuses = NULL;
    if( !(rset->flags & PT_REQSET_STATUSES_IGNORE) ) {
        rset->array_of_statuses = (MPI_Status*)malloc(count * sizeof(MPI_Status));
    }
    assert(0 == (flags & PT_REQSET_ACTIVE));

    /* Intent to check completion of requests indepdently ? */
    if( flags & PT_REQSET_CHECK_ANY ) {
        rset->idx_completed_reqs = (int *)malloc(count * sizeof(int));
    }
    return MPI_SUCCESS;
}

#define SET_AND_CHECK_REQUEST_SET_ID(ID, RSET, CODE) \
do { \
    if( (ID) >= MAX_REQSET_SIZE ) { \
        fprintf(stderr, "GID (%d) larger than the recorded requests sets (%d)\n", (ID), MAX_REQSET_SIZE ); \
        CODE; \
    } \
    (RSET) = &pt_reqset_array[(ID)]; \
    if( NULL == (RSET) ) { \
        fprintf(stderr, "GID (%d) is not recorded for any registered requests set\n", (ID)); \
        CODE; \
    } \
} while (0)

#define CHECK_ACTIVE_REQSET(RSET, GID, MSG, CODE)                                                                              \
    do {                                                                                                                       \
        if (0 == ((RSET)->flags & (PT_REQSET_ACTIVE | PT_REQSET_COMPLETED)) ) {                                                \
            fprintf(stderr, "Request set %d is neither active nor completed. It is illegal to "                                \
            "check its completion in %s\n", (GID), (MSG));                                                                     \
            CODE;                                                                                                              \
        }                                                                                                                      \
    } while (0)

int pt_reqset_unregister(int* gid)
{
    pt_reqset_t* rset;

    SET_AND_CHECK_REQUEST_SET_ID(*gid, rset, return MPI_ERR_ARG);
    if( rset->flags & PT_REQSET_ACTIVE ) {
        fprintf(stderr, "An active request set (gid = %d) cannot be unregistered\n", *gid);
        return MPI_ERR_ARG;
    }
    *gid = -1;
    free(rset->array_of_requests);
    rset->array_of_requests = NULL;
    if( NULL != rset->array_of_statuses ) {
        free(rset->array_of_statuses);
        rset->array_of_statuses = NULL;
    }
    
    return MPI_SUCCESS;
}

int pt_reqset_start(int gid)
{
    pt_reqset_t* rset;
    
    SET_AND_CHECK_REQUEST_SET_ID(gid, rset, return MPI_ERR_ARG);
    if( rset->flags & PT_REQSET_ACTIVE ) {
        fprintf(stderr, "An active request set (gid = %d) cannot be restarted before completion\n", gid);
        return MPI_ERR_ARG;
    }
    /* There is no need for an atomic operation here */
    rset->flags |= PT_REQSET_ACTIVE;
    rset->reported_completion = 0;
    rset->detected_completion = 0;

    /* signal the progress thread that a new command is there to be processed */
    pthread_mutex_lock(&mpi_helper->cq_mutex);
    pt_reqset_t* pt_command_queue = mpi_helper->cq;
    if( NULL == pt_command_queue ) {
        mpi_helper->cq = pt_command_queue = rset;
    } else {
        rset->prev = pt_command_queue->prev;
        rset->prev->next = rset;
    }
    pt_command_queue->prev = rset;
    rset->next = mpi_helper->cq;
    pthread_mutex_unlock(&mpi_helper->cq_mutex);

    /* Check the status of the progress thread and make sure it is running,
     * or nobody will be there to service our requests.
     */
    if( 1 == mpi_helper->be_quiet ) {
        pthread_mutex_unlock(&mpi_helper->quiet_mutex);
    }
    return MPI_SUCCESS;
}

static int pt_reqset_copy_status(pt_reqset_t* rset, MPI_Status* statuses)
{
    rset->flags ^= PT_REQSET_COMPLETED;
    if( (MPI_STATUSES_IGNORE != statuses) && !(rset->flags & PT_REQSET_STATUSES_IGNORE) ) {
        assert( NULL != rset->array_of_statuses);
        /* copy the statuses */
        memcpy(statuses, rset->array_of_statuses, rset->count * sizeof(MPI_Status));
    }
    if( rset->flags & PT_REQSET_NON_PERSISTENT ) {
        int gid = (int)(((uintptr_t)((char*)rset - (char*)&pt_reqset_array[0])) / sizeof(pt_reqset_array[0]));
        /* remove the gid requests set */
        return pt_reqset_unregister(&gid);
    }
    return MPI_SUCCESS;
}

int pt_reqset_test(int gid, int* flag, MPI_Status* statuses)
{
    pt_reqset_t* rset;

    SET_AND_CHECK_REQUEST_SET_ID(gid, rset, {return MPI_ERR_ARG;});
    CHECK_ACTIVE_REQSET(rset, gid, "pt_reqset_test", {return MPI_ERR_ARG;});

    *flag = 0;
    if( rset->flags & PT_REQSET_COMPLETED ) {
        *flag = 1;
        return pt_reqset_copy_status(rset, statuses);
    }
    return MPI_SUCCESS;
}

int pt_reqset_wait(int gid, MPI_Status* statuses)
{
    pt_reqset_t* rset;

    SET_AND_CHECK_REQUEST_SET_ID(gid, rset, {return MPI_ERR_ARG;});
    CHECK_ACTIVE_REQSET(rset, gid, "pt_reqset_test", {return MPI_ERR_ARG;});

    while( !(rset->flags & PT_REQSET_COMPLETED) ) {
        struct timespec ts = {.tv_nsec = 1000};
        /* do something */
        nanosleep(&ts, NULL);
    }
    if( rset->flags & PT_REQSET_COMPLETED ) {
        return pt_reqset_copy_status(rset, statuses);
    }
    return MPI_SUCCESS;
}

int pt_reqset_waitany(int gid, int* idx, MPI_Status* status)
{
    pt_reqset_t* rset;

    SET_AND_CHECK_REQUEST_SET_ID(gid, rset, {return MPI_ERR_ARG;});
    CHECK_ACTIVE_REQSET(rset, gid, "pt_reqset_test", {return MPI_ERR_ARG;});
    if( (rset->detected_completion == rset->reported_completion) &&
        !(rset->flags & PT_REQSET_COMPLETED) ) {
        while (rset->detected_completion == rset->reported_completion) {
            struct timespec ts = {.tv_nsec = 1000};
            /* do something */
            nanosleep(&ts, NULL);
        }
    }
    if( rset->detected_completion != rset->reported_completion) {
        atomic_thread_fence(memory_order_acquire); /* make sure idx_completed_reqs and array_of_statuses are sound */
        *idx = rset->idx_completed_reqs[rset->reported_completion];
        if( (MPI_STATUS_IGNORE != status) && !(rset->flags & PT_REQSET_STATUSES_IGNORE) ) {
            *status = rset->array_of_statuses[rset->reported_completion];
        }
        rset->reported_completion++;  /* Move to the next completion */
	if( rset->reported_completion == rset->count ) {
            assert(rset->reported_completion == rset->detected_completion);
            pt_reqset_copy_status(rset, MPI_STATUSES_IGNORE  /* ignore the status but unregister the reqset */);
	}
        return MPI_SUCCESS;
    }
    /* we should never reach this state, as reporting the last request should mark the reqset
     * as inactive, and any further waitany shall fail early.
     */
    return MPI_ERR_ARG;
}

int pt_reqset_testany(int gid, int *idx, MPI_Status *status)
{
    pt_reqset_t *rset;

    SET_AND_CHECK_REQUEST_SET_ID(gid, rset, {return MPI_ERR_ARG;});
    CHECK_ACTIVE_REQSET(rset, gid, "pt_reqset_test", {return MPI_ERR_ARG;});

    if (rset->detected_completion != rset->reported_completion) {
        atomic_thread_fence(memory_order_acquire); /* make sure idx_completed_reqs and array_of_statuses are sound */
        *idx = rset->idx_completed_reqs[rset->reported_completion];
        if( (MPI_STATUS_IGNORE != status) && !(rset->flags & PT_REQSET_STATUSES_IGNORE) ) {
            *status = rset->array_of_statuses[rset->reported_completion];
        }
        rset->reported_completion++; /* Move to the next completion */
        return MPI_SUCCESS;
    }
    if (rset->flags & PT_REQSET_COMPLETED) {
        /* we are completing the reqset */
        *idx = rset->count;
        return pt_reqset_copy_status(rset, MPI_STATUSES_IGNORE /* ignore the status but unregister the reqset */);
    }
    *idx = -1;  /* no completed requests */
    return MPI_SUCCESS;
}
/**
 * Wrapper functions for the Fortran interface
 */
int pt_reqset_register_f(int count, int* array_of_requests, int flags, int *gid)
{
    MPI_Request *myreqs;
    myreqs = malloc(count * sizeof(MPI_Request));
    for(int i = 0; i < count; i++ )
        myreqs[i] = MPI_Request_f2c(array_of_requests[i]);

    return pt_reqset_register(count, myreqs, flags | PT_REQSET_USE_REQUEST_ARRAY, gid);
}

void pt_reqset_test_f(int* gid, int* flag, int* c_err)
{
    pt_reqset_t* rset;
    SET_AND_CHECK_REQUEST_SET_ID(*gid, rset, {*c_err = MPI_ERR_ARG; return;});
    CHECK_ACTIVE_REQSET(rset, *gid, "pt_reqset_test", {*c_err = MPI_ERR_ARG; return;});
    if( rset->flags & PT_REQSET_COMPLETED ) {
    /*
        if( !(rset->flags & PT_REQSET_STATUSES_IGNORE) && (MPI_STATUSES_IGNORE != (MPI_Status*)statuses_f) ) {
            assert( NULL != rset->array_of_statuses);
            for( int i = 0; i < rset->count; i++ ) {
                MPI_Status_c2f(rset->array_of_statuses, &statuses_f[i * (sizeof(MPI_Status) / sizeof(int))]);
            }
        }*/
        if( rset->flags & PT_REQSET_NON_PERSISTENT ) {
            /* remove the gid requests set */
            *c_err = pt_reqset_unregister(gid);
            return;
        }
        rset->flags ^= PT_REQSET_COMPLETED;
    }
    *c_err = MPI_SUCCESS;
}

void pt_reqset_testany_f(int* gid, int* idx, int* status_f, int* c_err)
{
    MPI_Status status;
    *c_err = pt_reqset_testany(*gid, idx, &status);
    if( MPI_STATUS_IGNORE != (MPI_Status*)status_f ) {
        MPI_Status_c2f(&status, status_f);
    }
}

void pt_reqset_wait_f(int* gid, int* c_err)
{
    pt_reqset_t* rset;
    int flag;

    SET_AND_CHECK_REQUEST_SET_ID(*gid, rset, {*c_err = MPI_ERR_ARG; return;});
    CHECK_ACTIVE_REQSET(rset, *gid, "pt_reqset_test", {*c_err = MPI_ERR_ARG; return;});

    while( !(rset->flags & PT_REQSET_COMPLETED) ) {
        struct timespec ts = {.tv_nsec = 1000};
        /* do something */
        nanosleep(&ts, NULL);
    }
    pt_reqset_test_f(gid, &flag,c_err);
}

void pt_reqset_waitany_f(int gid, int* idx, int* c_err)
{
    MPI_Status status;
    *c_err = pt_reqset_waitany(gid, idx, &status);
    /*
    if( MPI_STATUS_IGNORE != (MPI_Status*)status_f ) {
        MPI_Status_c2f(&status, status_f);
	}*/
}
