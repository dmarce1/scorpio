/*
 * tag.h
 *
 *  Created on: Feb 25, 2013
 *      Author: dmarce1
 */

#ifndef TAG_H_
#define TAG_H_

#define TAG_RELAX  0
#define TAG_PHI    1
#define TAG_CHILD  2
#define TAG_COARSE 3
#define TAG_VDOWN  4
#define TAG_DPHI   5
#define TAG_FORCE  6
#define TAG_MOVE   7
#define TAG_FLUX   8
#define TAG_CF     9
#define TAG_C2P    10
#define TAG_ERROR  11
#define TAG_PHIFACE 12
#define TAG_FMM_MOMENT 13
#define TAG_FMM_4FORCE 14
#define TAG_FMM_BOUND 15
#define TAG_FMM_EXPANSION 16
#define TAG_DTUP 17
#define TAG_DTDOWN 18
#define TAG_MAX    19
#define SUBCYCLEMAX 16

int tag_gen( int fid, int gid, int face, int subcycle = 0 );

#endif /* TAG_H_ */
