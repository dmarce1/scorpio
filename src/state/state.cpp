#include "../defs.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "../oct_node/oct_face.h"
#include "state.h"
Real State::gamma = 5.0 / 3.0;
Real State::rho_floor = 1.0e-12;
Real State::ei_floor = 1.0e-20;
const int State::d_index = 0;
const int State::sx_index = 1;
const int State::sy_index = 2;
const int State::sz_index = 3;
const int State::et_index = 4;
const int State::tau_index = 5;
const int State::frac_index = 6;
bool State::cylindrical = true;
_3Vec State::drift_vel = 0.0;

Vector<Real, STATE_NF> State::x_scalar_flux_coeff(const _3Vec& X) {
    Vector<Real, STATE_NF> c = 0.0;
    if (cylindrical) {
        c[sx_index] = X[0] / sqrt(X[0] * X[0] + X[1] * X[1]);
        c[sy_index] = 0.0;
    } else {
        c[sx_index] = 0.0;
    }
    return c;
}

Vector<Real, STATE_NF> State::y_scalar_flux_coeff(const _3Vec& X) {
    Vector<Real, STATE_NF> c = 0.0;
    if (cylindrical) {
        c[sx_index] = X[1] / sqrt(X[0] * X[0] + X[1] * X[1]);
        c[sy_index] = 0.0;
    } else {
        c[sy_index] = 0.0;
    }
    return c;
}

Vector<Real, STATE_NF> State::z_scalar_flux_coeff(const _3Vec& X) {
    Vector<Real, STATE_NF> c = 0.0;
    c[sz_index] = 0.0;
    return c;
}

Real State::scalar_flux(const _3Vec& X) const {
    return pg(X);
}

Vector<Real, STATE_NF> State::x_vector_flux(const _3Vec& X) const {
    Vector<Real, STATE_NF> flux;
    const Real v = vx(X);
    flux = (*this) * (v - drift_vel[0]);
    if (State::cylindrical) {
        flux[sy_index] -= X[1] * pg(X);
    } else {
        flux[sx_index] += pg(X);
    }
    flux[et_index] += v * pg(X);
    return flux;
}

Vector<Real, STATE_NF> State::y_vector_flux(const _3Vec& X) const {
    Vector<Real, STATE_NF> flux;
    const Real v = vy(X);
    flux = (*this) * (v - drift_vel[1]);
    if (State::cylindrical) {
        flux[sy_index] += X[0] * pg(X);
    } else {
        flux[sy_index] += pg(X);
    }
    flux[et_index] += v * pg(X);
    return flux;
}

Vector<Real, STATE_NF> State::z_vector_flux(const _3Vec& X) const {
    Vector<Real, STATE_NF> flux;
    const Real v = vz();
    flux = (*this) * (v - drift_vel[2]);
    flux[sz_index] += pg(X);
    flux[et_index] += v * pg(X);
    return flux;
}

Vector<Real, STATE_NF> State::source(const _3Vec& X, Real t) const {
    Vector<Real, STATE_NF> s = 0.0;
    if (!cylindrical) {
        s[sx_index] += omega * (*this)[sy_index];
        s[sy_index] -= omega * (*this)[sx_index];
    } else {
        Real R = sqrt(X[0] * X[0] + X[1] * X[1]);
        s[sx_index] += pow((*this)[sy_index], 2) / (rho() * pow(R, 3));
    }
#ifdef DRIVING
    const Real period = 2.0 * M_PI / omega;
    if (t < DRIVING_TIME * period) {
        if (cylindrical) {
            s[sy_index] -= (*this)[sy_index] * DRIVING / period;
        } else {
            printf("Error in state - cartesian not yet supported for drivingt\n");
            abort();
        }
    }
#endif
    return s;
}

const char* State::field_name(int i) {
    assert(i >= 0);assert(i < STATE_NF);
    switch (i) {
    case State::d_index:
        return "d";
    case et_index:
        return "et";
    case sx_index:
        if (cylindrical) {
            return "sr";
        } else {
            return "sx";
        }
    case sy_index:
        if (cylindrical) {
            return "lz";
        } else {
            return "sy";
        }
    case tau_index:
        return "tau";
    case sz_index:
        return "sz";
    case pot_index:
        return "pot";
    default:
        static char* str = NULL;
        if (str != NULL) {
            free(str);
        }
        asprintf(&str, "frac%i", i - frac_index + 1);
        return str;
    }
}

void State::enforce_dual_energy_formalism(const _3Vec& X, const State& n1, const State& n2, const State& n3, const State& n4,
        const State& n5, const State& n6) {
    Real max_et;
    max_et = max(et(), n1.et());
    max_et = max(max_et, n2.et());
    max_et = max(max_et, n3.et());
    max_et = max(max_et, n4.et());
    max_et = max(max_et, n5.et());
    max_et = max(max_et, n6.et());
    if (et() - ek(X) > 0.1 * max_et) {
        (*this)[tau_index] = pow(max(et() - ek(X) - ed(), ei_floor), 1.0 / gamma);
    }
}

void State::floor(const _3Vec& X) {
    Real rho1, rho2;
    Real de = pot();
    (*this)[d_index] = max((*this)[d_index], rho_floor);
    (*this)[tau_index] = max(pow(ei_floor, 1.0 / gamma), (*this)[tau_index]);
    if (NFRAC > 1) {
        Real tot;
        for (int i = 0; i < NFRAC; i++) {
            (*this)[frac_index + i] /= rho();
        }
        tot = 0.0;
        for (int i = 0; i < NFRAC; i++) {
            (*this)[frac_index + i] = min(1.0, (*this)[frac_index + i]);
            (*this)[frac_index + i] = max(0.0, (*this)[frac_index + i]);
            tot += (*this)[frac_index + i];
        }
        if (tot <= 0.0) {
            for (int i = 0; i < NFRAC; i++) {
                (*this)[frac_index + i] = (1.0 / Real(NFRAC)) * rho();
            }
        } else {
            for (int i = 0; i < NFRAC; i++) {
                (*this)[frac_index + i] *= rho() / tot;
            }
        }
    }
    de -= pot();
    (*this)[et_index] += de;
}

State::State() :
        Vector<Real, STATE_NF>() {
}

State::State(const Vector<Real, STATE_NF>& v) :
        Vector<Real, STATE_NF>(v) {
    return;
}

Real State::rho() const {
    return (*this)[d_index];
}

Real State::sx() const {
    return (*this)[sx_index];
}

Real State::sy() const {
    return (*this)[sy_index];
}

Real State::sz() const {
    return (*this)[sz_index];
}

Real State::et() const {
    return (*this)[et_index];
}

Real State::vx(const _3Vec& X) const {
    Real v_x;
    if (cylindrical) {
        Real sr, st, R;
        R = sqrt(X[0] * X[0] + X[1] * X[1]);
        sr = sx();
        st = sy() / R;
        v_x = (X[0] * sr - X[1] * st) / R / rho();
    } else {
        v_x = sx() / rho();
    }
    return v_x + X[1] * omega;
}

Real State::vy(const _3Vec& X) const {
    Real v_y;
    if (cylindrical) {
        Real sr, st, R;
        R = sqrt(X[0] * X[0] + X[1] * X[1]);
        sr = sx();
        st = sy() / R;
        v_y = (X[1] * sr + X[0] * st) / R / rho();
    } else {
        v_y = sy() / rho();
    }
    return v_y - X[0] * omega;
}

Real State::vz() const {
    return sz() / rho();
}

Real State::ek(const _3Vec& X) const {
    const Real x = vx(X);
    const Real y = vy(X);
    const Real z = vz();
    return 0.5 * rho() * (x * x + y * y + z * z);
}

Real State::hd() const {
    const Real x = pow(rho() / PhysicalConstants::B, 1.0 / 3.0);
    Real h;
    if (x < 0.01) {
        h = 4.0 * PhysicalConstants::A / PhysicalConstants::B * x * x;
    } else {
        h = 8.0 * PhysicalConstants::A / PhysicalConstants::B * (sqrt(x * x + 1.0) - 1.0);
    }
    return max(h, 0.0);
}

Real State::ed() const {
    return max(hd() * rho() - pd(), 0.0);
}

Real State::pd() const {
    const Real x = pow(rho() / PhysicalConstants::B, 1.0 / 3.0);
    Real p;
    if (x < 0.01) {
        p = 1.6 * PhysicalConstants::A * pow(x, 5);
    } else {
        p = PhysicalConstants::A * (x * (2.0 * x * x - 3.0) * sqrt(x * x + 1.0) + 3.0 * asinh(x));
    }
    return max(p, 0.0);
}

Real State::ei(const _3Vec& X) const {
    if (et() - ek(X) > 0.001 * et()) {
        return et() - ek(X) - ed();
    } else {
        return pow(max((*this)[tau_index], 0.0), gamma);
    }
}

Real State::pg(const _3Vec& X) const {
    return max((gamma - 1.0) * ei(X) + pd(), 0.0);
}

Real State::cs(const _3Vec& X) const {
    Real x, dp_depsilon, dp_drho, cs2;
    x = pow(rho() / PhysicalConstants::B, 1.0 / 3.0);
    dp_drho = max(((8.0 * PhysicalConstants::A) / (3.0 * PhysicalConstants::B)) * x * x / sqrt(x * x + 1.0) + gamma * (gamma - 1.0) * ei(X) / rho(), 0.0);
    dp_depsilon = (gamma - 1.0) * rho();
    cs2 = (pg(X) / (rho() * rho())) * dp_depsilon + dp_drho;
    return sqrt(cs2);
}

Real State::max_abs_x_eigen(const _3Vec& X) const {
    Real s = fabs(vx(X)) + cs(X);
    return s;
}

Real State::max_abs_y_eigen(const _3Vec& X) const {
    Real s = fabs(vy(X)) + cs(X);
    return s;
}

Real State::max_abs_z_eigen(const _3Vec& X) const {
    Real s = fabs(vz()) + cs(X);
    return s;
}

void State::set_rho(Real a, int frac) {
    (*this)[d_index + frac] = a;
}

void State::set_sx(Real a) {
    (*this)[sx_index] = a;
}

void State::set_sy(Real a) {
    (*this)[sy_index] = a;
}

void State::set_sz(Real a) {
    (*this)[sz_index] = a;
}

void State::set_et(Real a) {
    (*this)[et_index] = a;
}

void State::set_tau(Real a) {
    (*this)[tau_index] = a;
}

void State::to_prim(const _3Vec& X) {
//	(*this).to_cartesian(X);
    for (int i = 1; i < STATE_NF; i++) {
        (*this)[i] /= rho();
    }
    if (cylindrical) {
        (*this)[sy_index] -= omega * (X[0] * X[0] + X[1] * X[1]);
    } else {
        (*this)[sx_index] += X[1] * omega;
        (*this)[sy_index] -= X[0] * omega;
    }
}

void State::from_prim(const _3Vec& X) {
    if (cylindrical) {
        (*this)[sy_index] += omega * (X[0] * X[0] + X[1] * X[1]);
    } else {
        (*this)[sx_index] -= X[1] * omega;
        (*this)[sy_index] += X[0] * omega;
    }
    for (int i = 1; i < STATE_NF; i++) {
        (*this)[i] *= rho();
    }
//	(*this).from_cartesian(X);
}

Real State::omega = 0.0;
