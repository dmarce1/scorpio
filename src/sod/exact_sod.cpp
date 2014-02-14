/**************************************************************************************************************/
/**************************************************************************************************************/
/**************************************************************************************************************/
/**************************************************************************************************************/
/**************************************************************************************************************/
/**************************************************************************************************************/
/**************************************************************************************************************/
/**************************************************************************************************************/
/**************************************************************************************************************/
/* This C code was adopted from the FORTRAN code exact_sod.f at http://www.itam.nsc.ru/OLD2/flowlib/SRC/sod.f */
/**************************************************************************************************************/
/**************************************************************************************************************/
/**************************************************************************************************************/
/**************************************************************************************************************/
/**************************************************************************************************************/
/**************************************************************************************************************/
/**************************************************************************************************************/
/**************************************************************************************************************/
/**************************************************************************************************************/
/**************************************************************************************************************/
/**************************************************************************************************************/
/**************************************************************************************************************/

#include "exact_sod.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static double pl, pr, rhol, rhor, cl, cr, Gamma;
static double func(double pm);
static double rtbis(double x1, double x2, double xacc);

#ifdef TEST

void main( void ) {
	double x;
	double t= 0.245;
	sod_init_t in;
	sod_state_t out;
	in.rhol = 1.0;
	in.rhor = 0.125;
	in.pl = 1.0;
	in.pr = 0.1;
	in.Gamma = 1.4;
	for( x = -1.0; x <= 1.0; x+= 0.01 ) {
		exact_sod( &out, &in, x, t );
		printf( "%e %e %e %e\n", x, out.rho, out.v, out.p );
	}

}

#endif

void exact_sod(sod_state_t* out, sod_init_t* in, double x, double t) {
	Gamma = in->gamma;
	const double mu2 = (Gamma - 1.) / (Gamma + 1.);
	const double xmax = 1.0;
	const int numcells = 500;
	double pm, pressure, rhoml, vs, vt;
	double vm, density, velocity, rhomr;
	int i;
	pl = in->pl;
	pr = in->pr;
	rhol = in->rhol;
	rhor = in->rhor;
	cl = sqrt(Gamma * pl / rhol);
	cr = sqrt(Gamma * pr / rhor);
	pm = rtbis(pr, pl, 1.E-16);
	rhoml = rhol * pow((pm / pl), (1. / Gamma));
	vm = 2.0 * cl / (Gamma - 1.) * (1.0 - pow((pm / pl), ((Gamma - 1.0) / (2.0 * Gamma))));
	rhomr = rhor * ((pm + mu2 * pr) / (pr + mu2 * pm));
	vs = vm / (1. - rhor / rhomr);
	vt = cl - vm / (1. - mu2);

	if (x <= -cl * t) {
		density = rhol;
	} else if (x <= -vt * t) {
		density = rhol * pow((-mu2 * (x / (cl * t)) + (1. - mu2)), (2. / (Gamma - 1.)));
	} else if (x <= vm * t) {
		density = rhoml;
	} else if (x <= vs * t) {
		density = rhomr;
	} else {
		density = rhor;
	}
	if (x <= -cl * t) {
		pressure = pl;
	} else if (x <= -vt * t) {
		pressure = pl * pow((-mu2 * (x / (cl * t)) + (1. - mu2)), (2. * Gamma / (Gamma - 1.)));
	} else if (x <= vs * t) {
		pressure = pm;
	} else {
		pressure = pr;
	}
	if (x <= -cl * t) {
		velocity = 0.0;
	} else if (x <= -vt * t) {
		velocity = (1 - mu2) * (x / t + cl);
	} else if (x <= vs * t) {
		velocity = vm;
	} else {
		velocity = 0.0;
	}

	out->rho = density;
	out->v = velocity;
	out->p = pressure;
}

static double func(double pm) {
	double rc;
	const double mu2 = (Gamma - 1.) / (Gamma + 1.);
	rc = -2 * cl * (1 - pow((pm / pl), ((-1 + Gamma) / (2 * Gamma)))) / (cr * (-1 + Gamma)) + (-1 + pm / pr) * sqrt(((1 - mu2) / (Gamma * (mu2 + pm / pr))));
	return rc;
}

static double rtbis(double x1, double x2, double xacc) {
	double rc;
	int done;
	const int JMAX = 100;
	int j;
	double dx, f, fmid, xmid;
	fmid = func(x2);
	f = func(x1);
	if (f * fmid >= 0.) {
		printf("root must be bracketed in rtbis\n ");
		abort();
	}
	if (f < 0.) {
		rc = x1;
		dx = x2 - x1;
	} else {
		rc = x2;
		dx = x1 - x2;
	}
	done = 0;
	for (j = 1; j <= JMAX; j++) {
		dx = dx * 5.E-1;
		xmid = rc + dx;
		fmid = func(xmid);
		if (fmid <= 0.) {
			rc = xmid;
		}
		if (fabs(dx) < xacc || fmid == 0.) {
			done = 1;
			break;
		}
	}
	if (done == 0) {
		printf("too many bisections in rtbis\n");
		abort();
	}
	return rc;
}
