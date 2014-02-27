#ifndef RECONSTRUCT_H_
#define RECONSTRUCT_H_

#include "defs.h"

#include "array3d.h"

class Reconstruct {
public:
	Reconstruct();
	virtual ~Reconstruct();
	void operator()(State[GNX], State[GNX], State[GNX]) const;
};


#endif /* RECONSTRUCT_H_ */
