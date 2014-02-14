#ifndef EULER1230_STATE_H_
#define EULER1230_STATE_H_

#include "../defs.h"

#include "../real.h"
#include "state.h"
#include "../vector.h"
#include "../oct_node/oct_face.h"
#include <stdio.h>
#include "../physical_constants.h"

#define NFRAC 2
#define STATE_NF (7+NFRAC)

class State: public Vector<Real, STATE_NF> {
protected:
    static _3Vec drift_vel;
    static Real omega;
public:
    void to_cartesian_momenta(const _3Vec X) {
        State U0 = *this;
        (*this)[State::sx_index] = U0.vx(X) * rho();
        (*this)[State::sy_index] = U0.vy(X) * rho();
    }
    void to_cylindrical_momenta(const _3Vec X) {
        Real R2 = X[0] * X[0] + X[1] * X[1];
        Real cos0 = X[0] / sqrt(R2);
        Real sin0 = X[1] / sqrt(R2);
        Real sx0 = (*this)[State::sx_index];
        Real sy0 = (*this)[State::sy_index];
        (*this)[State::sx_index] = sx0 * cos0 + sy0 * sin0;
        (*this)[State::sy_index] = sy0 * X[0] - sx0 * X[1];
    }
    static void set_gamma(Real g) {
        gamma = g;
    }
    Real abar() const {
        const Real a_He = 4.0;
        const Real a_CO = 1.0 / ((0.5 / 12.0) + (0.5 / 16.0));
        Real f_He = frac(1) / (frac(0) + frac(1));
        Real f_CO = frac(0) / (frac(0) + frac(1));
        return 1.0 / (f_He / a_He + f_CO / a_CO);
    }
    Real cv() const {
        Real a = abar();
        Real z = 0.5 * a;
        return ((z + 1.0) * PhysicalConstants::kb) / (a * PhysicalConstants::amu * (gamma - 1.0));
    }
    Real temp(const _3Vec x) const {
        return ei(x) / rho() / cv();
    }
    static void set_drift_vel(const _3Vec& v) {
        drift_vel = v;
    }
    static bool cylindrical;
    static Real get_omega() {
        return omega;
    }
    static Real gamma;
    static Real ei_floor;
    static Real rho_floor;
    static const int d_index;
    static const int sx_index;
    static const int sy_index;
    static const int sz_index;
    static const int et_index;
    static const int tau_index;
    static const int frac_index;
    void set_frac(int f, Real a) {
        (*this)[frac_index + f] = a;
    }
    Real frac(int i) const {
        return (*this)[frac_index + i];
    }
    State();
    State(const Vector<Real, STATE_NF>&);
    ~State() {
        return;
    }
    Vector<Real, STATE_NF> source(const _3Vec& X, Real t) const;
    Real max_abs_x_eigen(const _3Vec& X) const;
    Real max_abs_y_eigen(const _3Vec& X) const;
    Real max_abs_z_eigen(const _3Vec& X) const;
    Vector<Real, STATE_NF> x_vector_flux(const _3Vec& X) const;
    Vector<Real, STATE_NF> y_vector_flux(const _3Vec& X) const;
    Vector<Real, STATE_NF> z_vector_flux(const _3Vec& X) const;
    Real scalar_flux(const _3Vec& X) const;
    static Vector<Real, STATE_NF> x_scalar_flux_coeff(const _3Vec& X);
    static Vector<Real, STATE_NF> y_scalar_flux_coeff(const _3Vec& X);
    static Vector<Real, STATE_NF> z_scalar_flux_coeff(const _3Vec& X);
    void to_prim(const _3Vec& X);
    void from_prim(const _3Vec& X);
    static const char* field_name(int i);
    void floor(const _3Vec&);
    Real rho() const;
    Real tau() const;
    Real sz() const;
    Real sx() const;
    Real sy() const;
    Real vx(const _3Vec&) const;
    Real vy(const _3Vec&) const;
    Real vz() const;
    Real ek(const _3Vec&) const;
    Real et() const;
    Real ed() const;
    Real pd() const;
    Real hd() const;
    Real ei(const _3Vec& X) const;
    Real cs(const _3Vec& X) const;
    Real pg(const _3Vec& X) const;
    void set_rho(Real, int = 0);
    void set_sx(Real);
    void set_sy(Real);
    void set_sz(Real);
    void set_et(Real);
    void set_tau(Real);
    void reflect_on_z();

    void enforce_dual_energy_formalism(const _3Vec& X, const State& n1, const State& n2, const State& n3, const State& n4, const State& n5, const State& n6);

    static void set_omega(Real o) {
        omega = o;
    }
    Real rot_pot(const _3Vec& X) const {
        return rho() * rot_phi(X);
    }
    Real rot_phi(const _3Vec& X) const {
        const Real R2 = X[0] * X[0] + X[1] * X[1];
        return -0.5 * R2 * omega * omega;
    }

    static bool low_order_variable(int i) {
        return i == pot_index;
    }
    static bool smooth_variable(int i) {
        return i == pot_index;
    }
    static const int pot_index = STATE_NF - 1;
    Real pot() const {
        return (*this)[pot_index];
    }
    void set_pot(Real a) {
        (*this)[pot_index] = a;
    }
    State& operator=(const State& s) {
        this->Vector<Real, STATE_NF>::operator=(s);
        return *this;
    }
    State& operator=(const Vector<Real, STATE_NF>& s) {
        this->Vector<Real, STATE_NF>::operator=(s);
        return *this;
    }
    State(const State& s) {
        *this = s;
    }
    Real conserved_energy(const _3Vec X) {
        Real egas, egrav, erot;
        egas = et();
        if (cylindrical) {
            egrav = 0.5 * (pot() - rot_pot(X));
            erot = rot_pot(X);
            return egas + egrav + erot;
        } else {
            egrav = 0.5 * pot();
            return egas + egrav;
        }
    }
    void to_con(const _3Vec& X) {
        (*this)[et_index] += 0.5 * ((*this)[pot_index] + rot_pot(X));
    }
    void from_con(const _3Vec& X) {
        (*this)[et_index] -= 0.5 * ((*this)[pot_index] + rot_pot(X));
    }
    Real phi_eff() const {
        return pot() / rho();
    }
    Real phi(const _3Vec& x) const {
        return (pot() - rot_pot(x)) / rho();
    }
    Real virial(const _3Vec& x) const {
        Real ek = (vx(x) * vx(x) + vy(x) * vy(x) + vz() * vz()) * 0.5 * rho() + 1.5 * pd();
        Real ep = 0.5 * (pot() - rot_pot(x));
        return 2.0 * ek + ep;
    }
};

#endif /* EULER_STATE_H_ */
