#ifndef PROGRESS_HELPER_HEADER_HAS_BEEN_INCLUDED
#define PROGRESS_HELPER_HEADER_HAS_BEEN_INCLUDED

#include <mpi.h>

typedef enum pt_request_set_flags_e {
    PT_REQSET_NON_PERSISTENT    = (1 << 0),
    PT_REQSET_STATUSES_IGNORE   = (1 << 1),
    PT_REQSET_USE_REQUEST_ARRAY = (1 << 2),
    PT_REQSET_CHECK_ANY         = (1 << 3),
    /* Not to be set by the callee */
    PT_REQSET_ACTIVE            = (1 << 29),
    PT_REQSET_COMPLETED         = (1 << 30),
} pt_request_set_flags_e;

/**
 * Wait for one of the reqset requests to complete. Return it's index in the registered reqset as well as
 * it's status (if not MPI_STATUS_IGNORE). If all completions have been reported, which means that the reqset
 * is completed and unregistered, return an index equal to the number of requests in the reqset.
 */
int pt_reqset_waitany(int gid, int *idx, MPI_Status *status);
/**
 * Test is one of the reqset requests has completed. Return -1 in idx if no request has completed, the index
 * of the request in the registered reqset as well as it's status (if not MPI_STATUS_IGNORE) if completed requests
 * exists. If all completions have been reported, which means that the reqset
 * is completed and unregistered, return an index equal to the number of requests in the reqset.
 */
int pt_reqset_testany(int gid, int *idx, MPI_Status *status);
/**
 * Wait until all requests of the reqset have completed. Copy the status if requested.
 */
int pt_reqset_wait(int gid, MPI_Status *statuses);
/**
 * Check if all the requests in the reqset have completed. If yes, set the flag to 1, unregister the reqset
 * and return the statuses if requested.
 */
int pt_reqset_test(int gid, int* flag, MPI_Status* statuses);
/**
 * Start the reqset corresponding to the provided gid.
 */
int pt_reqset_start(int gid);
/**
 * Unregister an inactive reqset. If the reqset is active this operation will fail.
 */
int pt_reqset_unregister(int* gid);
/**
 * Register an array of request as a reqset. The reqset will remain inactive until a pt_reqset_start is called.
 */
int pt_reqset_register(int count, MPI_Request* array_of_requests, int flags, int* gid);

#endif  /* PROGRESS_HELPER_HEADER_HAS_BEEN_INCLUDED */
