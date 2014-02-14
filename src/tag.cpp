#include "tag.h"
#include <stdio.h>
#include "assert.h"
#include <stdlib.h>

int tag_gen(int fid, int gid, int index, int subcycle) {
/*	if (index < 0 || index >= 8) {
		printf("%i-----------\n", index);
	}*/
	assert( index >= 0);
	assert( index < 26);
	int tag = (index + 26 * (fid + TAG_MAX * gid)) * (SUBCYCLEMAX + 1) + subcycle;
	if (tag < 0) {
		printf("Tag out of range %i %i %i %i %i\n", fid, gid, index, subcycle, tag);
		abort();
	}
//	printf("%i\n", tag);
	return tag;
}

