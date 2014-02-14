#include <malloc.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "legendre.h"
#include <unistd.h>

Real factorial(int i) {
	if (i == 0) {
		return 1.0;
	} else {
		return factorial(i - 1) * Real(i);
	}
}

void AssociatedLegendrePolynomial::generate(Real x) {
	static Real** cons;
	static bool initialized = false;
	if (initialized == false) {
		cons = new Real*[lmax + 1];
		for (int l = 0; l <= lmax; l++) {
			cons[l] = new Real[2 * l + 1];
			cons[l] += l;
			for (int m = -l; m <= l; m++) {
				cons[l][m] = factorial(l + m) / factorial(l - m);
			}
		}
		initialized = true;
	}
	Real l0, m0, tmp;
	int mp1;
	P[0][0] = 1.0;
	if (lmax > 0) {
		if (x > -1.0 && x < 1.0) {
			tmp = 1.0 / sqrt(1.0 - x * x);
			P[1][0] = x;
			for (int l = 1; l < lmax; l++) {
				l0 = (Real) l;
				P[l + 1][0] = ((2.0 * l0 + 1.0) * x * P[l][0] - l0 * P[l - 1][0]) / (l0 + 1.0);
			}
			for (int l = 1; l <= lmax; l++) {
				l0 = (Real) l;
				for (int m = 0; m < l; m++) {
					m0 = (Real) m;
					mp1 = m + 1;
					P[l][mp1] = ((l0 - m0) * x * P[l][m] - (l0 + m0) * P[l - 1][m]) * tmp;
					P[l][-mp1] = (cons[l][-mp1] * P[l][mp1]) * (2 - (mp1 & 0x1));
				}
			}
			for (int l = 1; l <= lmax; l++) {
				for (int m = -l; m <= l; m++) {
					P[l][m] *= sqrt(cons[l][-m]);
				}
			}
		} else {
			for (int l = 1; l <= lmax; l++) {
				P[l][0] = pow(x, l);
				for (int m = -l; m <= l; m++) {
					if (m != 0) {
						P[l][m] = 0.0;
					}
				}
			}
		}
	}
}

