#include <mpi.h>
#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
 
typedef struct mpi_helper_s {
  MPI_Comm     comm;
  MPI_Request  sync_req;
  pthread_t    thread;
  int          recvint;
} mpi_helper_t;
 
static mpi_helper_t* mpi_helper = NULL;
 
void* mpi_helper_thread_routine( void* args )
{
  mpi_helper_t* helper = (mpi_helper_t*)args;
 
  MPI_Wait(&helper->sync_req, MPI_STATUS_IGNORE);
  /* nothing else to do, just return */
  return NULL;
}
 
void start_MPI_helper_(void)
{
  int rc, my_rank;
 
  mpi_helper = (mpi_helper_t*)malloc(sizeof(mpi_helper_t));
  rc = MPI_Comm_dup(MPI_COMM_WORLD, &mpi_helper->comm);
  rc = MPI_Comm_rank(mpi_helper->comm, &my_rank);
  rc = MPI_Irecv(&mpi_helper->recvint, 1, MPI_INT, my_rank, 0, mpi_helper->comm,
                  &mpi_helper->sync_req);
  rc = pthread_create(&mpi_helper->thread, NULL, mpi_helper_thread_routine, mpi_helper);
}
 
void stop_MPI_helper_(void)
{
  int rc = 0, my_rank;
  void* ret;
 
  rc = MPI_Comm_rank(mpi_helper->comm, &my_rank);
  /* send the message to complete the helper thread */
  rc = MPI_Send(&rc, 1, MPI_INT, my_rank, 0, mpi_helper->comm);
  /* wait until the helper thread completes */
  pthread_join(mpi_helper->thread, &ret);
  MPI_Comm_free(&mpi_helper->comm);
  free(mpi_helper);
  mpi_helper = NULL;
}
