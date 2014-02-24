/*
 * hydro_grav_grid.cpp
 *
 *  Created on: Mar 20, 2013
 *      Author: dmarce1
 */

#include "hydro_grav_grid.h"
#ifdef HYDRO_GRAV_GRID

void HydroGravGrid::write(FILE* fp) const {
    HydroGrid::write(fp);
    MultiGrid::write(fp);
}

void HydroGravGrid::read(FILE* fp) {
    HydroGrid::read(fp);
    MultiGrid::read(fp);
}

HydroGravGrid::HydroGravGrid() {
    // TODO Auto-generated constructor stub

}

HydroGravGrid::~HydroGravGrid() {
    // TODO Auto-generated destructor stub
}

void HydroGravGrid::deallocate_arrays() {
    HydroGrid::deallocate_arrays();
    MultiGrid::deallocate_arrays();
}

void HydroGravGrid::allocate_arrays() {
    HydroGrid::allocate_arrays();
    MultiGrid::allocate_arrays();
}

void HydroGravGrid::set_refine_flags() {
    HydroGrid::set_refine_flags();
}

HydroGravGrid* HydroGravGrid::new_octnode() const {
    return new HydroGravGrid;
}

Real HydroGravGrid::xf(int i) const {
    return HydroGrid::xf(i);
}

Real HydroGravGrid::yf(int i) const {
    return HydroGrid::yf(i);
}

Real HydroGravGrid::zf(int i) const {
    return HydroGrid::zf(i);
}

Real HydroGravGrid::xc(int i) const {
    return HydroGrid::xc(i);
}

Real HydroGravGrid::yc(int i) const {
    return HydroGrid::yc(i);
}

Real HydroGravGrid::zc(int i) const {
    return HydroGrid::zc(i);
}

Real HydroGravGrid::get_dx() const {
    return HydroGrid::get_dx();
}

int HydroGravGrid::nvar_output() const {
    if (using_shadow()) {
        return 2 * HydroGrid::nvar_output() + 6;
    } else {
        return HydroGrid::nvar_output() + 6;
    }
}

Real HydroGravGrid::get_output_point(int i, int j, int k, int l) const {
    const int o = BW - 1;
    int n;
    const int m = HydroGrid::nvar_output();
    if (using_shadow()) {
        n = 2 * m;
    } else {
        n = m;
    }
    if (l < m) {
        return HydroGrid::get_output_point(i, j, k, l);
    } else if (l < n) {
        return get_shadow_error(i, j, k)[l - n];
    } else if (l == n) {
        return (*this)(i, j, k).pg(HydroGrid::X(i, j, k));
    } else if (l == n + 1) {
        return get_phi(i - o, j - o, k - o);
    } else if (l == n + 2) {
        return 0.5 * (get_fx(i - o + 1, j - o, k - o) + get_fx(i - o, j - o, k - o));
    } else if (l == n + 3) {
        return 0.5 * (get_fy(i - o, j + 1 - o, k - o) + get_fy(i - o, j - o, k - o));
    } else if (l == n + 4) {
        return 0.5 * (get_fz(i - o, j - o, k + 1 - o) + get_fz(i - o, j - o, k - o));
    } else {
        return get_dphi(i - o, j - o, k - o);
    }
}

const char* HydroGravGrid::output_field_names(int l) const {
    static char str[256];
    int n;
    const int m = HydroGrid::nvar_output();
    if (using_shadow()) {
        n = 2 * m;
    } else {
        n = m;
    }
    if (l < m) {
        return HydroGrid::output_field_names(l);
    } else if (l < n) {
        sprintf(str, "%s_error", HydroGrid::output_field_names(l - m));
        return str;
    } else if (l == n) {
        return "p";
    } else if (l == n + 1) {
        return "phi";
    } else if (l == n + 2) {
        return "gx";
    } else if (l == n + 3) {
        return "gy";
    } else if (l == n + 4) {
        return "gz";
    } else {
        return "dphi";
    }
}

void HydroGravGrid::init() {
    HydroGrid::init();
    MultiGrid::init();
}

void HydroGravGrid::create_child(const ChildIndex& c) {
    Real phi;
    const int o = BW - 1;
    HydroGravGrid* g;
    HydroGrid::create_child(c);
    MultiGrid::create_multigrid_child(c);
    if (proc() == MPI_rank()) {
        g = dynamic_cast<HydroGravGrid*>(get_child(c));
        for (int k = BW; k < GNX - BW; k++) {
            for (int j = BW; j < GNX - BW; j++) {
                for (int i = BW; i < GNX - BW; i++) {
                    phi = (*g)(i, j, k).pot();
                    phi /= (*g)(i, j, k).rho();
                    g->set_phi(i - o, j - o, k - o, phi);
                }
            }
        }
    }
}

void HydroGravGrid::initialize() {

}

void HydroGravGrid::compute_dudt(int dir) {
    HydroGrid::compute_dudt(dir);
    if (dir == 0) {
        _3Vec x;
        Real d;
        Vector<Real, STATE_NF> D, D0;
        for (int k = BW; k < GNX - BW; k++) {
            for (int j = BW; j < GNX - BW; j++) {
                for (int i = BW; i < GNX - BW; i++) {
                    x = HydroGrid::X(i, j, k);
                    d = (*this)(i, j, k).rho();
                    D = 0.0;
                    D0 = get_dudt(i, j, k);
                    Real R = sqrt(x[0] * x[0] + x[1] * x[1]);
                    D[State::lz_index] = d * (x[0] * gy(i, j, k) - x[1] * gx(i, j, k));
                    D[State::sx_index] = d * gx(i, j, k);
                    D[State::sy_index] = d * gy(i, j, k);
                    D[State::sz_index] = d * gz(i, j, k);
                    D[State::et_index] = D0[State::pot_index] - (*this)(i, j, k).pot() / d * D0[State::d_index];
                    D[State::pot_index] = -D0[State::pot_index];
                    add_to_dudt(i, j, k, D);
                }
            }
        }
    }
}

#endif
