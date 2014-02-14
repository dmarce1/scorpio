#include "../defs.h"
#ifdef HYDRO_GRAV_GRID
#include "wd.h"
#include "../state/state.h"
#include "single_star.h"


#define WDM 0.2

SingleStar::SingleStar() {
	// TODO Auto-generated constructor stub

}

SingleStar::~SingleStar() {
	// TODO Auto-generated destructor stub
}

Real SingleStar::radius(int i, int j, int k) {
	const Real x0 = 0.0e+10;
	const Real y0 = 0.0e+8;
	const Real z0 = 0.0e+9;
	/*const Real x0 = 0.0;
	 const Real y0 = 0.0;
	 const Real z0 = 0.0;*/
	return sqrt(pow(xc(i) - x0, 2) + pow(yc(j) - y0, 2) + pow(zc(k) - z0, 2));
}

static Real wd_d0, wd_r0;
static bool wd_init = false;

void SingleStar::set_refine_flags() {
	ChildIndex c;
	if (!wd_init) {
		wd_init = true;
		wd_rho0(WDM, &wd_d0, &wd_r0);
		State::rho_floor = wd_d0 * 1.0e-15;
		printf("%e %e\n", wd_d0, wd_r0);
		State::rho_floor = wd_d0 * 1.0e-12;
		//State::nohydro = true;
	}
	if (get_level() < 2) {
		for (int i = 0; i < OCT_NCHILD; i++) {
			set_refine_flag(i, true);
		}
	} else if (get_level() < get_max_level_allowed()) {
		for (int k = BW; k < GNX - BW; k++) {
			c.set_z(2 * k / GNX);
			for (int j = BW; j < GNX - BW; j++) {
				c.set_y(2 * j / GNX);
				for (int i = BW; i < GNX - BW; i++) {
					c.set_x(2 * i / GNX);
					if (!get_refine_flag(c)) {
						if ((*this)(i, j, k).rho() > wd_d0 * 1.0e-6) {
							//					if (get_shadow_error(i, j, k).mag() > 0.00001 * (1 << 2 * get_level())) {
							set_refine_flag(c, true);
						} else if (get_time() == 0.0 && radius(i, j, k) < get_dx() * 2.0) {
							set_refine_flag(c, true);
						}
					}
				}
			}
		}
	}
}

void SingleStar::initialize() {
    if (!wd_init) {
        wd_init = true;
        wd_rho0(WDM, &wd_d0, &wd_r0);
        State::rho_floor = wd_d0 * 1.0e-12;
        printf("%e %e\n", wd_d0, wd_r0);
        //State::nohydro = true;
    }
	State U;
	Real r, r0, d, e, tau;
	for (int k = BW; k < GNX - BW; k++) {
		for (int j = BW; j < GNX - BW; j++) {
			for (int i = BW; i < GNX - BW; i++) {
				U = Vector<Real, STATE_NF>(0.0);
				r = radius(i, j, k);
				if (r < wd_r0) {
					d = wd_density(wd_d0, r, r / 100.0);
					//	printf( "%e\n", d);
				} else {
					d = 0.0;
				}
				d = max(d, 1.0e-12 * wd_d0);
				U.set_rho(d);
				U.set_et(U.ed());
				if (r < 0.66666 * wd_r0) {
					U.set_frac(0, d);
				} else {
					U.set_frac(1, d);
				}
				(*this)(i, j, k) = U;
			}
		}
	}
}

#endif
