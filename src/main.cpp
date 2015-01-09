#include "program.h"
#include <mpi.h>
#include <fenv.h>
#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>

extern "C" {
void backtrace_symbols_fd2(void * const *buffer, int size, int fd);
}

void handler(int sig) {
	void *array[15];
	size_t size;
	size = backtrace(array, 15);
	fprintf(stderr, "Error: signal %d:\n", sig);
	backtrace_symbols_fd(array, size, 2);
	exit(1);
}

int main(int argc, char* argv[]) {
	int rc;
	MPI_Init(NULL, NULL);
	//signal(SIGSEGV, handler);
	rc = Program().run(argc, argv);
	MPI_Finalize();
	return rc;
}

