#ifndef __LEGENDRE__
#define __LEGENDRE__

#include <stdio.h>
#include "complex.h"
#include "real.h"
#include "vector.h"

class AssociatedLegendrePolynomial {
private:
	Real** P;
	int lmax;
public:
	AssociatedLegendrePolynomial(int _lmax) {
		lmax = _lmax;
		P = new Real*[lmax + 1];
		for (int l = 0; l <= lmax; l++) {
			P[l] = new Real[2 * l + 1];
			P[l] += l;
		}
	}
	~AssociatedLegendrePolynomial() {
		for (int l = 0; l <= lmax; l++) {
			P[l] -= l;
			delete[] P[l];
		}
		delete[] P;
	}
	void generate(Real x);
	Real get(int l, int m) {
		return P[l][m];
	}
};

class SphericalHarmonic {
private:
	AssociatedLegendrePolynomial P;
	Complex** Y;
	int lmax;
public:
	SphericalHarmonic(int _lmax) :
			P(_lmax) {
		lmax = _lmax;
		Y = new Complex*[lmax + 1];
		for (int l = 0; l <= lmax; l++) {
			Y[l] = new Complex[2 * l + 1];
			Y[l] += l;
		}
	}
	~SphericalHarmonic() {
		for (int l = 0; l <= lmax; l++) {
			Y[l] -= l;
			delete[] Y[l];
		}
		delete[] Y;
	}
	void generate(const _3Vec& v) {
		generate(v[0], v[1], v[2]);
	}
	void generate(Real x, Real y, Real z) {
		Real R = sqrt(x * x + y * y);
		Real r = sqrt(R * R + z * z);
		Real p, psi;
		P.generate(z / r);
		for (int l = 0; l <= lmax; l++) {
			for (int m = -l; m <= l; m++) {
				p = P.get(l, m);
				if (R != 0.0) {
					psi = atan2(y, x);
					Y[l][m].set_real(p * cos(Real(m) * psi));
					Y[l][m].set_imag(p * sin(Real(m) * psi));
				} else {
					Y[l][m].set_real(p);
					Y[l][m].set_imag(p);
				}
			}
		}
	}
	Complex operator()(int l, int m) const {
		if (fabs(m) >= l || l > lmax) {
			//	printf( "%e\n", Y[l][m].real());
			return Y[l][m];// * sqrt((2.0 * l + 1.0) / (4.0 * M_PI));
		} else {
			return Complex(0.0);
		}
	}
};

#endif
