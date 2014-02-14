#include "oct_face.h"
//#include "oct_node.h"

OctFace invert(OctFace f) {
	static const OctFace op_face[OCT_NSIB] = { XU, XL, YU, YL, ZU, ZL };
	return op_face[f];
}

