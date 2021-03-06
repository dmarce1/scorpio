#include "lane_emden.h"
#include "binary_star.h"
#include "../indexer3d.h"

#include "dwd.h"

#ifdef HYDRO_GRAV_GRID

#ifdef BINARY_STAR

binary_parameters_t BinaryStar::bparam;
bool BinaryStar::bparam_init = false;
Real BinaryStar::dtheta = 0.0;

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

void BinaryStar::physical_boundary(int dir) {
    Vector<int, 3> lb, ub;
    for (int dir = 0; dir < 3; dir++) {
        if (is_phys_bound(2 * dir + 0)) {
            lb = BW;
            ub = GNX - BW - 1;
            lb[dir] = 0;
            ub[dir] = BW - 1;
            for (Indexer3d i(lb, ub); !i.end(); i++) {
                Vector<int, 3> k = i;
                k[dir] = BW;
                _3Vec x = X(i[0], i[1], i[2]);
                U(i[0], i[1], i[2]) = U(k[0], k[1], k[2]);
                switch (dir) {
                case 0:
                    U(i[0], i[1], i[2]).set_sx(min(U(i[0], i[1], i[2]).sx(), 0.0));
                    break;
                case 1:
                    U(i[0], i[1], i[2]).set_sy(min(U(i[0], i[1], i[2]).sy(), 0.0));
                    break;
                case 2:
                    U(i[0], i[1], i[2]).set_sz(min(U(i[0], i[1], i[2]).sz(), 0.0));
                    break;
                }
                U(i[0], i[1], i[2]).set_lz_from_cartesian(x);
            }
        }
        if (is_phys_bound(2 * dir + 1)) {
            lb = BW;
            ub = GNX - BW - 1;
            lb[dir] = GNX - BW;
            ub[dir] = GNX - 1;
            for (Indexer3d i(lb, ub); !i.end(); i++) {
                Vector<int, 3> k = i;
                k[dir] = GNX - BW - 1;
                _3Vec x = X(i[0], i[1], i[2]);
                U(i[0], i[1], i[2]) = U(k[0], k[1], k[2]);
                switch (dir) {
                case 0:
                    U(i[0], i[1], i[2]).set_sx(max(U(i[0], i[1], i[2]).sx(), 0.0));
                    break;
                case 1:
                    U(i[0], i[1], i[2]).set_sy(max(U(i[0], i[1], i[2]).sy(), 0.0));
                    break;
                case 2:
                    U(i[0], i[1], i[2]).set_sz(max(U(i[0], i[1], i[2]).sz(), 0.0));
                    break;
                }
                U(i[0], i[1], i[2]).set_lz_from_cartesian(x);
            }
        }
    }
    Hydro::inc_instruction_pointer();
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
        dxmin = dynamic_cast<BinaryStar*>(get_root())->get_dx() / Real(1 << OctNode::get_max_level_allowed());
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
                            _3Vec tmp;
                            tmp = 0.0;
                            tmp[0] = bparam.x1;
                            Real ra = (X(i, j, k) - tmp).mag();
                            tmp = 0.0;
                            tmp[0] = bparam.x2;
                            Real rd = (X(i, j, k) - bparam.x2).mag();
                            if ((*this)(i, j, k).rho() > refine_floor * refine_adjustment) {
                                set_refine_flag(c, true);
                            } else if (get_time() == 0.0 && (ra < 2.0 * get_dx() || rd < 2.0 * get_dx())) {
                                set_refine_flag(c, true);
                            } else if (get_time() != 0.0) {
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

void BinaryStar::initialize() {
    if (!bparam_init) {
        bparam_init = true;
        bparam.fill_factor = 1.00;
        binary_parameters_compute(&bparam);
        State::rho_floor = min(1.0e-15 * bparam.rho1, 1.0e-13 * bparam.rho2);
        refine_floor = min(1.0e-6 * bparam.rho1, 1.0e-4 * bparam.rho2);
        dynamic_cast<Hydro*>(get_root())->Hydro::mult_dx(bparam.a * 128.0 * 0.96 / 16.0);
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
                if (id == 1) {
                    U.set_frac(0, rho);
                    U.set_frac(1, 0.0);
                } else if (id == -1) {
                    U.set_frac(1, rho);
                    U.set_frac(0, 0.0);
                } else {
                    U.set_frac(1, 0.0);
                    U.set_frac(0, 0.0);
                }
                U.set_et(U.ed());
                U.set_tau(tau);
                U.set_sx(0.0);
                if (rho > State::rho_floor) {
                    U.set_sy(bparam.omega * R2 * U.rho());
                } else {
                    U.set_sy(0.0);
                }
                U.set_sz(0.0);

                (*this)(i, j, k) = U;
            }
        }
    }

}

#endif
#endif
