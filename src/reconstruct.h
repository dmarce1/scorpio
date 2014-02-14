#ifndef RECONSTRUCT_H_
#define RECONSTRUCT_H_

#include "defs.h"

#include "array3d.h"

#ifdef USE_HYDRO_GRID
class Reconstruct {
protected:
	void vanleer(State[GNX], State[GNX], State[GNX]) const;
	void minmod(State[GNX], State[GNX], State[GNX]) const;
	void ppm(State[GNX], State[GNX], State[GNX]) const;
public:
	Reconstruct();
	virtual ~Reconstruct();
	void operator()(State[GNX], State[GNX], State[GNX]) const;
};

#endif

#endif /* RECONSTRUCT_H_ */
