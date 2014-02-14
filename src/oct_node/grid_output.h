
#ifndef GRID_OUTPUT_H_
#define GRID_OUTPUT_H_

#include "../real.h"
#include "../defs.h"

typedef double OReal;
#define DB_OREAL DB_DOUBLE

struct grid_output_t {
	OReal* x;
	OReal* y;
	OReal* z;
	OReal** ele;
	int* nodelist;
	int ni;
	int pi;
	int ei;
};

#endif /* GRID_OUTPUT_H_ */
