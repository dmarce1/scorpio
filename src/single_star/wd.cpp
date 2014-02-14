#include <math.h>
#include <stdio.h>
#include "../physical_constants.h"


static double y_dot(double x, double y, double r) {
	double a, b, c, x2p1, x3;
	const double c0 = M_PI * PhysicalConstants::G * PhysicalConstants::B * PhysicalConstants::B / (2.0 * PhysicalConstants::A);
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

void wd(double rho0, double* rptr, double* mptr) {
	double x, y, rho, r, dr, dx1, dy1, dx2, dy2, mass;
	int N, i;
	N = 100;
	x = pow(rho0 / PhysicalConstants::B, 1.0 / 3.0);
	y = 0.0;
	r = 0.0;
	mass = 0.0;
	dr = sqrt(fabs(x / y_dot(x, y, r))) / (double) N;
	for (i = 0; x > 0.0; i++) {
		rho = PhysicalConstants::B * x * x * x;
		//printf( "%e %e\n", r, rho );
		dy1 = y_dot(x, y, r) * dr;
		dx1 = x_dot(x, y, r) * dr;
		dy2 = y_dot(x + dx1, y + dy1, r + dr) * dr;
		dx2 = x_dot(x + dx1, y + dy1, r + dr) * dr;
		mass += rho * dr * 4.0 * M_PI * r * r;
		y += (dy1 + dy2) / 2.0;
		x += (dx1 + dx2) / 2.0;
		r += dr;
	}
	mass /= 1.99e+33;
//	r /= 6.96e+10;
	*mptr = mass;
	*rptr = r;
}

static double wd_density(double rho0, double r0, double dr) {
	double x, y, r, dx1, dy1, dx2, dy2;
	int i, imax;
	x = pow(rho0 / PhysicalConstants::B, 1.0 / 3.0);
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
	return PhysicalConstants::B * x * x * x;
}

void wd_rho0(double mass, double* rho0ptr, double* rptr) {
	double dmax, dmin, dmid, m, r;
	dmax = 1.0e+13;
	dmin = 1.0e+2;
	while (log(dmax / dmin) > 1.0e-6) {
		dmid = (dmax + dmin) / 2.0;
		wd(dmid, &r, &m);
		if (m > mass) {
			dmax = dmid;
		} else {
			dmin = dmid;
		}
	}
	*rho0ptr = dmid;
	*rptr = r;
}
