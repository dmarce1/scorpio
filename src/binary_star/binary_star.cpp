#include "lane_emden.h"
#include "binary_star.h"

#include "dwd.h"

#ifdef HYDRO_GRAV_GRID

#ifdef BINARY_STAR

binary_parameters_t BinaryStar::bparam;
bool BinaryStar::bparam_init = false;

void BinaryStar::compute_flow_off() {
	const Real da = get_dx() * get_dx();
	int k, j, i;
	static State this_DFO;
	this_DFO = Vector<Real, STATE_NF>(0.0);
	for (k = BW; k < GNX - BW; k++) {
		for (j = BW; j < GNX - BW; j++) {
			for (i = BW; i < GNX - BW; i++) {
				if (!zone_is_refined(i, j, k)) {
					if (is_phys_bound(XL) && i == BW) {
						this_DFO -= (get_flux(0, i, j, k)) * da;
					}
					if (is_phys_bound(XU) && i == GNX - BW - 1) {
						this_DFO += (get_flux(0, i + 1, j, k)) * da;
					}
					if (is_phys_bound(YL) && j == BW) {
						this_DFO -= (get_flux(1, i, j, k)) * da;
					}
					if (is_phys_bound(YU) && j == GNX - BW - 1) {
						this_DFO += (get_flux(1, i, j + 1, k)) * da;
					}
					if (is_phys_bound(ZL) && k == BW) {
						this_DFO -= (get_flux(2, i, j, k)) * da;
					}
					if (is_phys_bound(ZU) && k == GNX - BW - 1) {
						this_DFO += (get_flux(2, i, j, k + 1)) * da;
					}
				}
			}
		}
	}
	DFO += this_DFO;
}

BinaryStar::BinaryStar() {
}

BinaryStar::~BinaryStar() {
	// TODO Auto-generated destructor stub
}

Real BinaryStar::radius(int i, int j, int k) {
	return sqrt(xc(i) * xc(i) + yc(j) * yc(j) + zc(k) * zc(k));
}
void BinaryStar::set_refine_flags() {
	double refine_adjustment = 1;
	ChildIndex c;
	if (get_level() < 1) {
		for (int i = 0; i < OCT_NCHILD; i++) {
			set_refine_flag(i, true);
		}
	} else if (get_level() < get_max_level_allowed()) {
		Real mass_min, dxmin, vmin, this_mass;
		dxmin = dynamic_cast<BinaryStar*>(get_root())->get_dx()
				/ Real(1 << OctNode::get_max_level_allowed());
		vmin = dxmin * dxmin * dxmin;
#ifdef REFINE_ACC_MORE
		mass_min = refine_floor * vmin * refine_adjustment * 8.0;
#else
		mass_min = refine_floor * vmin * refine_adjustment;
#endif
		for (int k = BW; k < GNX - BW; k++) {
			c.set_z(2 * k / GNX);
			for (int j = BW; j < GNX - BW; j++) {
				c.set_y(2 * j / GNX);
				for (int i = BW; i < GNX - BW; i++) {
					c.set_x(2 * i / GNX);
#ifdef REFINE_ACC_MORE
					if ((get_level() == get_max_level_allowed() - 1
							&& (*this)(i, j, k).frac(0)
									> refine_floor * refine_adjustment)
							|| get_level() < get_max_level_allowed() - 1) {
#else
						if (get_level() < get_max_level_allowed()) {
#endif
						if (!get_refine_flag(c)) {
							//              set_refine_flag(c, true);
							Real ra = (X(i, j, k) - bparam.x1).mag();
							Real rd = (X(i, j, k) - bparam.x2).mag();
							if ((*this)(i, j, k).rho()
									> refine_floor * refine_adjustment) {
								set_refine_flag(c, true);
							} else if (get_time() == 0.0
									&& (ra < get_dx() || rd < get_dx())) {
								set_refine_flag(c, true);
							} else/* if (get_time() != 0.0)*/{
								this_mass = pow(get_dx(), 3)
										* (*this)(i, j, k).rho();
								if (this_mass > mass_min) {
									set_refine_flag(c, true);
								}
							}
						}
					}
				}
			}
		}
	}
}
/*
 void BinaryStar::set_refine_flags() {
 ChildIndex c;
 if (get_level() < 1) {
 for (int i = 0; i < OCT_NCHILD; i++) {
 set_refine_flag(i, true);
 }
 } else if (get_level() < get_max_level_allowed()) {
 Real mass_min, dxmin, vmin, this_mass;
 dxmin = dynamic_cast<BinaryStar*>(get_root())->get_dx() / Real(1 << OctNode::get_max_level_allowed());
 vmin = dxmin * dxmin * dxmin;
 mass_min = refine_floor * vmin;
 for (int k = BW; k < GNX - BW; k++) {
 c.set_z(2 * k / GNX);
 for (int j = BW; j < GNX - BW; j++) {
 c.set_y(2 * j / GNX);
 for (int i = BW; i < GNX - BW; i++) {
 c.set_x(2 * i / GNX);
 if ((get_level() == get_max_level_allowed() - 1 && (*this)(i, j, k).frac(0) > refine_floor) || get_level() < get_max_level_allowed() - 1) {
 if (!get_refine_flag(c)) {
 //		set_refine_flag(c, true);
 Real ra = (X(i, j, k) - a0).mag();
 Real rd = (X(i, j, k) - d0).mag();
 if ((*this)(i, j, k).rho() > refine_floor) {
 set_refine_flag(c, true);
 } else if (get_time() == 0.0 && (ra < get_dx() || rd < get_dx())) {
 set_refine_flag(c, true);
 } else{
 this_mass = pow(get_dx(), 3) * (*this)(i, j, k).rho();
 if (this_mass > mass_min) {
 set_refine_flag(c, true);
 }
 }
 }
 }
 }
 }
 }
 }
 }
 */

void BinaryStar::initialize() {
	if (!bparam_init) {
		bparam_init = true;
		bparam.fill_factor = 0.97;
		binary_parameters_compute(&bparam);
		State::rho_floor = 1.0e-12 * bparam.rho1;
		refine_floor = 1.0e-4 * bparam.rho2;
		dynamic_cast<HydroGrid*>(get_root())->HydroGrid::mult_dx(
				bparam.a * 5.0);
#ifndef USE_FMM
		dynamic_cast<MultiGrid*>(get_root())->MultiGrid::mult_dx(
				bparam.a * 5.0);
#endif
		State::set_omega(bparam.omega);
	}
	for (int k = BW - 1; k < GNX - BW + 1; k++) {
		for (int j = BW - 1; j < GNX - BW + 1; j++) {
			for (int i = BW - 1; i < GNX - BW + 1; i++) {
				int id;
				State U = Vector<Real, STATE_NF>(0.0);
				Real R2 = (xc(i) * xc(i) + yc(j) * yc(j));
				Real rho = density_at(&bparam, xc(i), yc(j), zc(k), &id);
				rho = max(rho, State::rho_floor);
				Real tau = pow(State::ei_floor, 1.0 / State::gamma);
				U.set_rho(rho);
				U.set_et(U.ed());
				U.set_tau(tau);
				U.set_sx(0.0);
				U.set_sy(bparam.omega * R2 * U.rho());
				U.set_sz(0.0);
				if (id == 1) {
					U.set_frac(0, U.rho());
					U.set_frac(1, 0.0);
				} else if (id == -1) {
					U.set_frac(1, U.rho());
					U.set_frac(0, 0.0);
				} else {
					U.set_frac(1, 0.0);
					U.set_frac(0, 0.0);
				}
				(*this)(i, j, k) = U;
			}
		}
	} /*Real K, E1, E2, period, rho_c1, rho_c2;
	 a0 = d0 = 0.0;
	 if (scf_code) {
	 const Real a = 0.075;
	 const Real R2 = 0.0075;
	 const Real n = 1.5;
	 rho_c2 = q_ratio * M1 * (pow(3.65375 / R2, 3.0) / (2.71406 * 4.0 * M_PI));
	 rho_c1 = rho_c2 / pow(q_ratio, 2.0 * n / (3.0 - n));

	 a0[0] = a * q_ratio / (q_ratio + 1.0);
	 d0[0] = -a / (q_ratio + 1.0);
	 K = pow(R2 / 3.65375, 2) * 4.0 * M_PI / (2.5) * pow(rho_c2, 1.0 / 3.0);
	 Ka = Kd = (5.0 / 8.0) * K;
	 polyK = K;
	 E1 = 1.0 / (sqrt((4.0 * M_PI * pow(rho_c1, 1.0 - 1.0 / n)) / ((n + 1.0) * K)));
	 E2 = 1.0 / (sqrt((4.0 * M_PI * pow(rho_c2, 1.0 - 1.0 / n)) / ((n + 1.0) * K)));
	 //	printf("%e %e %e %e\n", rho_c1, rho_c2, E1, E2);
	 period = sqrt(pow(a, 3) / (q_ratio + 1.0) / M1) * 2.0 * M_PI;
	 } else {
	 rho_c1 = rho_c2 = 1.0;
	 a0[0] = 0.025;
	 d0[0] = -0.025;
	 E1 = E2 = 0.0015;
	 K = pow(E1, 2) * 4.0 * M_PI / (2.5);
	 period = sqrt(pow((a0 - d0).mag(), 3) / 2.303394e-07) * 2.0 * M_PI;
	 }
	 State U;
	 Real d, e, gamma, tau;
	 gamma = 5.0 / 3.0;
	 Real ra, rd, r0, ek;
	 _3Vec v;
	 const Real Omega = 2.0 * M_PI / period;
	 State::set_omega(Omega);
	 Real f1, f2;
	 if (get_level() == 0) {
	 printf("Period = %e Omega = %e\n", period, Omega);
	 }
	 Real d_floor = 10.0 * State::rho_floor;
	 for (int k = BW - 1; k < GNX - BW + 1; k++) {
	 for (int j = BW - 1; j < GNX - BW + 1; j++) {
	 for (int i = BW - 1; i < GNX - BW + 1; i++) {
	 U = Vector<Real, STATE_NF>(0.0);
	 ra = (X(i, j, k) - a0).mag();
	 rd = (X(i, j, k) - d0).mag();
	 d = +rho_c1 * pow(lane_emden(ra / E1), 1.5);
	 d += rho_c2 * pow(lane_emden(rd / E2), 1.5);
	 d = max(d, d_floor);
	 if (ra < rd) {
	 f1 = d - d_floor / 2.0;
	 f2 = d_floor / 2.0;
	 } else {
	 f2 = d - d_floor / 2.0;
	 f1 = d_floor / 2.0;
	 }
	 if (State::cylindrical) {
	 Real R2 = (HydroGrid::xc(i) * HydroGrid::xc(i) + HydroGrid::yc(j) * HydroGrid::yc(j));
	 v[0] = 0.0;
	 v[1] = R2 * State::get_omega();
	 } else {
	 v[0] = -HydroGrid::yc(j) * State::get_omega();
	 v[1] = +HydroGrid::xc(i) * State::get_omega();

	 }
	 v[2] = 0.0;
	 e = K * pow(d, gamma) / (gamma - 1.0);
	 tau = pow(e, 1.0 / gamma);
	 U.set_rho(d);
	 U.set_frac(0, f1);
	 U.set_frac(1, f2);
	 U.set_et(e);
	 U.set_tau(tau);
	 U.set_sx(d * v[0]);
	 U.set_sy(d * v[1]);
	 U.set_sz(d * v[2]);
	 (*this)(i, j, k) = U;
	 }
	 }
	 }
	 */
}

#endif
#endif
