#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

const double pi = 3.1415926535897932384e0, eulercon = 0.577215664901532861e0, a2rad = pi / 180.0e0, rad2a = 180.0e0
		/ pi;

// physical constants
const double g = 6.6742867e-8, h = 6.6260689633e-27, hbar = 0.5e0 * h / pi, qe = 4.8032042712e-10,
		avo = 6.0221417930e23, clight = 2.99792458e10, kerg = 1.380650424e-16, ev2erg = 1.60217648740e-12, kev = kerg
				/ ev2erg, amu = 1.66053878283e-24, mn = 1.67492721184e-24, mp = 1.67262163783e-24,
		me = 9.1093821545e-28, rbohr = hbar * hbar / (me * qe * qe), fine = qe * qe / (hbar * clight), hion =
				13.605698140e0;

const double ssol = 5.67051e-5, asol = 4.0e0 * ssol / clight, weinlam = h * clight / (kerg * 4.965114232e0), weinfre =
		2.821439372e0 * kerg / h, rhonuc = 2.342e14;

// astronomical constants
const double msol = 1.9892e33, rsol = 6.95997e10, lsol = 3.8268e33, mearth = 5.9764e27, rearth = 6.37e8, ly =
		9.460528e17, pc = 3.261633e0 * ly, au = 1.495978921e13, secyer = 3.1558149984e7;

#include "dwd.h"

const double A = 6.0228e+22;
const double B = 2.0*9.81011e+5;

const double C0 = pi * g * B * B / A / 2.0;

double dydr(double y, double x, double r) {
	double a, b, c;
	a = C0 * x * x * sqrt(1.0 + x * x);
	if (r != 0.0) {
		b = 2.0 * y / r;
	} else {
		b = 0.0;
	}
	c = y * y / x / (1.0 + x * x);
	return -(a + b + c);
}

double dxdr(double y, double x, double r) {
	return y;
}
double dmdr(double y, double x, double r) {
	return 4.0 * pi * r * r * B * x * x * x;
}

double WD(double rho0, double* radius) {
	double x, y, r, rho;
	double dr = 1.0e+2;
	double dx1, dx2, dx3, dy1, dy2, dy3, dx, dy, m, dm1, dm2, dm3, dm, enth0,
			enth1;
	const int imax = 100000;
	const int imin = 50000;
	int i;
	bool redo = true;
	while( redo ) {
		x = pow(rho0 / B, 1.0 / 3.0);
		y = 0.0;
		m = 0.0;
		for ( i = 0; i < imax; i++) {
			rho = B * pow(x, 3);
	//	printf("%e %e %e\n", r, rho, m);
			r = double(i) * dr;
			dx1 = dxdr(y, x, r) * dr;
			dy1 = dydr(y, x, r) * dr;
			dm1 = dmdr(y, x, r) * dr;
			dx2 = dxdr(y + dy1, x + dx1, r + dr) * dr;
			dy2 = dydr(y + dy1, x + dx1, r + dr) * dr;
			dm2 = dmdr(y + dy1, x + dx1, r + dr) * dr;
			dx3 = dxdr(y + (dy1 + dy2) / 4.0, x + (dx1 + dx2) / 4.0, r + dr / 2.0)
					* dr;
			dy3 = dydr(y + (dy1 + dy2) / 4.0, x + (dx1 + dx2) / 4.0, r + dr / 2.0)
					* dr;
			dm3 = dmdr(y + (dy1 + dy2) / 4.0, x + (dx1 + dx2) / 4.0, r + dr / 2.0)
					* dr;
			dx = (dx1 + dx3 + 4.0 * dx2) / 6.0;
			*radius = r;
				x += dx;
			if (x <= 5.0e-4 * pow(rho0 / B, 1.0 / 3.0) ) {
			//	printf( "%e\n", m / 1.98e+33);
				if( i >= imin ) {
					redo = false;
				} else {
					dr /= 2.0;
				}
		//		printf( "%i %e\n", i, dr );
				break;
			}
			y += dy;
			m += dm;
			dy = (dy1 + dy3 + 4.0 * dy2) / 6.0;
			dm = (dm1 + dm3 + 4.0 * dm2) / 6.0;
		}
		if( i == imax ) {
			dr *= 2.0;
		}
	}
	return m / msol;

}

double find_l1( double q ) {

	double y, f, df, f1, f2, fr, df1, df2, dfr;
	const double qp1 = 1.0 + q;
	const double qp1inv = 1.0 / qp1;

	y = 0.0;
	do {
		f1 = -    pow(y + q*qp1inv,-2);
		f2 = +q * pow(y -   qp1inv,-2);
		fr = qp1*y;
		f = f1 + f2 + fr;
		df1 = +2.0*    pow(y + q*qp1inv,-3);
		df2 = -2.0*q * pow(y -   qp1inv,-3);
		dfr = qp1;
		df = df1 + df2 + dfr;
		y -= f / df;
	//	printf( "%e %e %e\n", y, df, f );
	} while( fabs( f ) > 1.0e-10 );
	return y;

}

double rand1() {
	return (double) rand() / (double) RAND_MAX;
}


double find_V( double q ) {

	const double qp1 = 1.0 + q;
	const double qp1inv = 1.0 / qp1;
	double x, y, z;
	double fx, fy, fz, r1inv, r2inv, phi, phi_l1;
	double h = 5.0e-3;
	int in;
	double r1inv3, r2inv3;
	const double l1_x = find_l1( q );
	r1inv = 1.0/sqrt( pow(l1_x+q*qp1inv,2)  );
	r2inv = 1.0/sqrt( pow(l1_x-  qp1inv,2) );
	phi_l1 = -1.0 *r1inv - q*r2inv  - 0.5*qp1*(l1_x*l1_x);
	in = 0;
	double dx;

	for( x = l1_x; x < 1.0+l1_x; x += h ) {
	for( y = h/2.0; y < 0.5; y += h ) {
	for( z = h/2.0; z < 0.5; z += h ) {
		r1inv = 1.0/sqrt( pow(x+q*qp1inv,2) + y*y + z*z );
		r2inv = 1.0/sqrt( pow(x-  qp1inv,2) + y*y + z*z );
		phi = -r1inv - q*r2inv  - 0.5*qp1*(x*x+y*y);
		if( phi < phi_l1 ) {
			r1inv3 = r1inv*r1inv*r1inv;
			r2inv3 = r2inv*r2inv*r2inv;
			dx = x - qp1inv;
			fx = -(x+q*qp1inv)*r1inv3 - q*dx*r2inv3 + qp1*x;
			fy = -y*r1inv3 - y*q*r2inv3 + qp1*y;
			fz = -z*r1inv3 - z*q*r2inv3;
			if( fx*dx + fy*y + fz*z <= 0.0 ) {
				in++;
			}
		}
	}}}
	return 4.0*in*h*h*h;
}


double wd_radius( double mass, double* rho0 ) {
	double rho_min, rho_max, rho_mid;
	double test_mass;
	rho_min = 1.0e+3;
	rho_max = 1.0e+12;
	double r;
	do {
		rho_mid = sqrt(rho_min*rho_max);
		test_mass = WD(rho_mid,&r);
		if( test_mass > mass ) {
			rho_max = rho_mid;
		} else {
			rho_min = rho_mid;
		}
	//	printf( "%e %e %e %e %e\n", rho_min, rho_mid, rho_max, test_mass, mass );
	} while( log(rho_max/rho_min) > 1.0e-9 );
	*rho0 = rho_mid;
	return r;
}

double binary_separation( double accretor_mass, double donor_mass, double fill_factor = 1.0) {
	double rho0_donor, rho0_accretor;
	double q = donor_mass / accretor_mass;
	double donor_radius = wd_radius( donor_mass, &rho0_donor );
	double normalized_roche_volume = find_V( q )*fill_factor;
	double roche_radius = pow(normalized_roche_volume / (4.0/3.0*M_PI),1.0/3.0);
	double separation = donor_radius / roche_radius;
	return separation;
}


void binary_parameters_compute( binary_parameters_t* b ) {
	b->q = b->m2 / b->m1;
	b->r1 = wd_radius( b->m1, &b->rho1 );
	b->r2 = wd_radius( b->m2, &b->rho2 );
	double normalized_roche_volume = find_V( b->q )*b->fill_factor;
	double roche_radius = pow(normalized_roche_volume / (4.0/3.0*M_PI),1.0/3.0);
	b->a = b->r2 / roche_radius;
	b->x1 = -b->q/(b->q+1.0)*b->a;
	b->x2 = +1.0 /(b->q+1.0)*b->a;
	b->omega = sqrt(g*msol*(1.0+b->q)/pow(b->a,3));
	b->P = 2.0*M_PI / b->omega;
}



static double y_dot(double x, double y, double r) {
	double a, b, c, x2p1, x3;
	const double c0 = M_PI * g * B * B / (2.0 * A);
	x2p1 = x * x + 1.0;
	if (r != 0.0) {
		x3 = x * x * x;
		a = 1.0 / (r * x * x2p1);
		b = y * (2.0 * x + 2.0 * x3 + r * y);
		c = c0 * x3 * r * sqrt(x2p1) * x2p1;
	} else {
		b = x2p1;
		a = c0 * x * x;
		c = 0.0;
	}
	return -a * (b + c);
}

static double x_dot(double x, double y, double r) {
	return y;
}



double wd_density(double rho0, double r0, double dr) {
	double x, y, r, dx1, dy1, dx2, dy2;
	int i, imax;
	x = pow(rho0 / B, 1.0 / 3.0);
	y = 0.0;
	r = 0.0;
	imax = int(r0 / dr + 0.5);
	for (i = 0; i < imax; i++) {
		dy1 = y_dot(x, y, r) * dr;
		dx1 = x_dot(x, y, r) * dr;
		dy2 = y_dot(x + dx1, y + dy1, r + dr) * dr;
		dx2 = x_dot(x + dx1, y + dy1, r + dr) * dr;
		y += (dy1 + dy2) / 2.0;
		x += (dx1 + dx2) / 2.0;
		r += dr;
		if (x < 0.0) {
			x = 0.0;
			break;
		}
	}
	//	printf("%e %e %e %i\n",dr, r, PhysCon::B*x*x*x, imax);
	return B * x * x * x;
}


double density_at( binary_parameters_t* b, double x, double y, double z, int* star ) {
	double r12, r22, rho;
	r12 = (pow( x - b->x1, 2 ) + y*y + z*z);
	r22 = (pow( x - b->x2, 2 ) + y*y + z*z);
	if( r12 < b->r1*b->r1  ) {
		rho = wd_density( b->rho1, sqrt(r12), b->r1/100.0 );
		*star = +1;
	} else if( r22 < b->r2*b->r2 ) {
		rho = wd_density( b->rho2, sqrt(r22),  b->r2/100.0 );
		*star = -1;
	} else {
		rho = 0.0;
		*star = 0;
	}
	return rho;
}
