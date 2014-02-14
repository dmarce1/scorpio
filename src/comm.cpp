#include "comm.h"
#include <mpi.h>

int MPI_rank() {
	int proc;
	MPI_Comm_rank(MPI_COMM_WORLD, &proc);
	return proc;
}

int MPI_size() {
	int MPIsz;
	MPI_Comm_size(MPI_COMM_WORLD, &MPIsz);
	return MPIsz;
}
