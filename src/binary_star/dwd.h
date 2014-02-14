/*
 * dwd.h
 *
 *  Created on: Oct 24, 2013
 *      Author: dmarce1
 */

#ifndef DWD1_H_
#define DWD1_H_


typedef struct {
	double m1, m2, q, rho1, rho2, a, r1, r2, fill_factor, x1, x2, omega, P;
} binary_parameters_t;


double density_at( binary_parameters_t* b, double x, double y, double z, int* star );

void binary_parameters_compute( binary_parameters_t* b );



#endif /* DWD_H_ */
