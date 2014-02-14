#ifndef REAL123_H_
#define REAL123_H_

#include <math.h>


typedef double Real;

Real max(Real, Real);
Real min(Real, Real);
Real max(Real, Real,Real);
Real min(Real, Real,Real);
Real sgn(Real);
Real abs(Real);
int ceiling(Real);
Real minmod(Real, Real);
Real minmod(Real, Real, Real);
Real minmod_theta(Real, Real, Real);
int nint(Real);

inline int max(int a, int b) {
	return a > b ? a : b;
}

inline int min(int a, int b) {
	return a < b ? a : b;
}

inline Real max(Real a, Real b) {
	return a > b ? a : b;
}

inline Real max(Real a, Real b, Real c) {
	return max(a, max(b, c));
}

inline Real min(Real a, Real b, Real c) {
	return min(a, min(b, c));
}

inline Real abs(Real a) {
	return fabs(a);
}

inline int nint(Real a) {
	return (int) (a + 0.5);
}

inline Real min(Real a, Real b) {
	return a > b ? b : a;
}

inline int ceiling(Real a) {
	int b = int(a);
	if (Real(int(a)) != a) {
		b++;
	}
	return b;
}

inline Real sgn(Real a) {
	if (a > 0.0) {
		return 1.0;
	} else if (a < 0.0) {
		return -1.0;
	} else {
		return 0.0;
	}
}

inline Real minmod(Real a, Real b) {
	return 0.5 * (sgn(a) + sgn(b)) * min(fabs(a), fabs(b));
}

inline Real minmod_theta(Real a, Real b, Real theta) {
	return minmod(theta * a, theta * b, 0.5 * (a + b));
}

inline Real minmod(Real a, Real b, Real c) {
	return minmod(a, minmod(b, c));
}

#endif /* REAL_H_ */
