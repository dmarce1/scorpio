/*
 * poisson.h
 *
 *  Created on: Mar 17, 2013
 *      Author: dmarce1
 */

#ifndef POISSON_H_
#define POISSON_H_

#ifndef PNX
#define PNX 0
#endif

#define LMAX 4

#include "../multigrid/multigrid.h"

class Poisson: public MultiGrid {
private:
	static Real com[4];
	static Real** i_poles;
	static Real** r_poles;
	static Real compute_phi(Real, Real, Real);
	static void pole_init();
	static void poles_clear();
	static void poles_reduce();
	void poles_compute();
	void compute_local_physical_boundaries();
	void accumulate_com();
	static void compute_com();


protected:
	static void compute_physical_boundaries();
	static _3Vec get_center_of_mass();
public:
	static Real source_sum() {
		return r_poles[0][0];
	}
	static void run();
	Poisson();
	virtual ~Poisson();
};

#endif /* POISSON_H_ */
