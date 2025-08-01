#include <mpi.h>
#include "mpix_harmonize.h"

int MPIX_Harmonize_f2c(MPI_Fint fcomm, MPI_Fint *flag)
{
	return MPIX_Harmonize(PMPI_Comm_f2c(fcomm), flag);
}

