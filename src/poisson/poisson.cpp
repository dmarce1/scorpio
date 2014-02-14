#include "poisson.h"
#include "../indexer3d.h"
#include "../legendre.h"

void Poisson::accumulate_com() {
	const int DX = (PNX / 2 - 1);
	const Real dv = pow(get_dx(), 3);
	Real dm;
	for (ChildIndex ci = 0; ci < 8; ci++) {
		if (get_child(ci) == NULL) {
			for (int k = 1 + ci.get_z() * DX; k < PNX / 2 + ci.get_z() * DX; k++) {
				for (int j = 1 + ci.get_y() * DX; j < PNX / 2 + ci.get_y() * DX; j++) {
					for (int i = 1 + ci.get_x() * DX; i < PNX / 2 + ci.get_x() * DX; i++) {
						dm = get_source(i, j, k) * dv;
						com[0] += dm;
						com[1] += Poisson::xc(i)*dm;
						com[2] += Poisson::yc(j)*dm;
						com[3] += Poisson::zc(k)*dm;
					}
				}
			}
		}
	}
}

Real Poisson::compute_phi(Real x, Real y, Real z) {
	static AssociatedLegendrePolynomial P(LMAX);
	Real real = 0.0;
	Real r, theta, phi, prpow, mphi, rpow;
	x -= com[1];
	y -= com[2];
	z -= com[3];
	r = sqrt(x * x + y * y + z * z);
	theta = acos(z / r);
	phi = atan2(y, x);
	P.generate(cos(theta));

	rpow = 1.0;
	for (int l = 0; l <= LMAX; l++) {
		rpow /= r;
		for (int m = -l; m <= l; m++) {
			prpow = P.get(l, m) * rpow;
			mphi = Real(m) * phi;
			real -= prpow * (r_poles[l][m] * cos(mphi) - i_poles[l][m] * sin(mphi));
		}
	}
	return real;

}

void Poisson::compute_local_physical_boundaries() {
	Vector<int, 3> ub, lb;
	for (int f = 0; f < OCT_NSIB; f++) {
		if (is_phys_bound(f)) {
			lb = 0;
			ub = PNX - 1;
			lb[f / 2] = ub[f / 2] = (f % 2) * (PNX - 1);
			for (Indexer3d i(lb, ub); !i.end(); i++) {
				set_phi(i[0], i[1], i[2], Poisson::compute_phi(MultiGrid::xc(i[0]), MultiGrid::yc(i[1]), MultiGrid::zc(i[2])));
			}
		}
	}
}

void Poisson::poles_compute() {
	static AssociatedLegendrePolynomial P(LMAX);
	int l, m, xlb, xub, ylb, yub, zlb, zub;
	Real theta, phi, r, x, y, z, prpow, mphi;
	Real dv = pow(get_dx(), 3) / (4.0 * M_PI);
	Real x2, z2, y2;
	for (ChildIndex ci = 0; ci < 8; ci++) {
		if (get_child(ci) == NULL) {
			xlb = 1 + ci.get_x() * (PNX / 2 - 1);
			ylb = 1 + ci.get_y() * (PNX / 2 - 1);
			zlb = 1 + ci.get_z() * (PNX / 2 - 1);
			xub = xlb + (PNX / 2 - 1) - 1;
			yub = ylb + (PNX / 2 - 1) - 1;
			zub = zlb + (PNX / 2 - 1) - 1;
			for (int j = ylb; j <= yub; j++) {
				y = MultiGrid::yc(j) - com[2];
				y2 = y * y;
				for (int i = xlb; i <= xub; i++) {
					x = MultiGrid::xc(i) - com[1];
					x2 = x * x;
					phi = atan2(y, x);
					for (int k = zlb; k <= zub; k++) {
						z = MultiGrid::zc(k) - com[3];
						z2 = z * z;
						r = sqrt(x2 + y2 + z2);
						theta = acos(z / r);
						P.generate(cos(theta));
						for (l = 0; l <= LMAX; l++) {
							for (m = -l; m <= l; m++) {
								prpow = P.get(l, m) * pow(r, l) * dv * get_source(i, j, k);
								mphi = Real(m) * phi;
								r_poles[l][m] += prpow * cos(mphi);
								i_poles[l][m] -= prpow * sin(mphi);
							}
						}
					}
				}
			}
		}
	}
}

Poisson::Poisson() {
// TODO Auto-generated constructor stub

}

Poisson::~Poisson() {
// TODO Auto-generated destructor stub
}

