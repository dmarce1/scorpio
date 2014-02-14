#include "assert.h"

#ifndef NDEBUG
#include <stdio.h>
#include <stdlib.h>

void __assert(bool b, const char* f, int l) {
	if (!b) {
		printf("ASSERTION FAILED ON LINE %i IN %s\n", l, f);
		abort();
	}
}

#endif
