#ifndef __SEDOV__H
#define __SEDOV__H

extern "C" {
void sedov_solution( double t, int N, const double* , double E, double rho,  double* dout, double* eout, double* vout );
}


#endif
