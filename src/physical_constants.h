/*
 * physical_constants.h
 *
 *  Created on: May 2, 2013
 *      Author: dmarce1
 */

#ifndef PHYSICAL_CONSTANTS_H_
#define PHYSICAL_CONSTANTS_H_

#include "real.h"

class PhysicalConstants {
public:
	static Real G, A, B, amu, kb, sigma, c;
	static void set_cgs();
	PhysicalConstants();
	virtual ~PhysicalConstants();
};

#endif /* PHYSICAL_CONSTANTS_H_ */
