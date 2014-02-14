/*
 * physical_constants.cpp
 *
 *  Created on: May 2, 2013
 *      Author: dmarce1
 */

#include "physical_constants.h"

Real PhysicalConstants::A = 1.0;
Real PhysicalConstants::B = 1.0;
Real PhysicalConstants::G = 1.0;
Real PhysicalConstants::amu = 1.0;
Real PhysicalConstants::kb = 1.0;
Real PhysicalConstants::sigma = 1.0;
Real PhysicalConstants::c = 1.0;

void PhysicalConstants::set_cgs() {
	G = 6.67259e-8;
	A = 6.0023e+22;
	B = 2.0 * 9.7393e+5;
	amu = 1.6605402e-24;
	kb = 1.380658e-16;
	sigma = 5.67051e-5;
	c = 2.99792458e+10;
}

PhysicalConstants::PhysicalConstants() {
	// TODO Auto-generated constructor stub

}

PhysicalConstants::~PhysicalConstants() {
	// TODO Auto-generated destructor stub
}

