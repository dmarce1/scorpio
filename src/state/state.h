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

#define STATE_NF 10

class State: public Vector<Real, STATE_NF> {
public:
    static Real gamma;
    static Real ei_floor;
    static Real rho_floor;
    static const int d_index = 0;
    static const int sx_index = 1;
    static const int sy_index = 2;
    static const int lz_index = 3;
    static const int sz_index = 4;
    static const int et_index = 5;
    static const int tau_index = 6;
    static const int pot_index = 7;
    static const int frac_index = STATE_NF - NFRAC;
    static Real omega;
    static Real omega0;
    static Real omega_dot;
    static _3Vec com_correction;
public:
    static void set_gamma(Real g);
    static Real get_omega();
    void set_frac(int f, Real a);
    Real frac(int i) const;
    State();
    State(const Vector<Real, STATE_NF>&);
    ~State();
    Vector<Real, STATE_NF> source(const _3Vec& X, Real t) const;
    Real max_abs_x_eigen(const _3Vec& X) const;
    Real max_abs_y_eigen(const _3Vec& X) const;
    Real max_abs_z_eigen(const _3Vec& X) const;
    Vector<Real, STATE_NF> x_flux(const _3Vec& X) const;
    Vector<Real, STATE_NF> y_flux(const _3Vec& X) const;
    Vector<Real, STATE_NF> z_flux(const _3Vec& X) const;
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
    void set_lz(Real);
    void set_tau(Real);
    void set_sx(Real);
    void set_sy(Real);
    void set_sz(Real);
    void set_et(Real);
    Real lz() const;
    Real set_lz_from_cartesian(const _3Vec& X);
    void reflect_on_z();
    void enforce_dual_energy_formalism(const _3Vec& X, const State& n1, const State& n2, const State& n3, const State& n4, const State& n5, const State& n6);
    static void set_omega(Real o);
    Real rot_pot(const _3Vec& X) const;
    Real rot_phi(const _3Vec& X) const;
    static bool low_order_variable(int i);
    static bool smooth_variable(int i);
    Real pot() const;
    void set_pot(Real a);
    State& operator=(const State& s);
    State& operator=(const Vector<Real, STATE_NF>& s);
    State(const State& s);
    Real conserved_energy(const _3Vec X);
    void to_con(const _3Vec& X);
    void from_con(const _3Vec& X);
    Real phi_eff() const;
    Real phi(const _3Vec& x) const;
    Real virial(const _3Vec& x) const;
};

#endif /* EULER_STATE_H_ */
