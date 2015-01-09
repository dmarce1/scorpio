#include "../defs.h"

#include "state.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "../oct_node/oct_face.h"
#include "state.h"
Real State::gamma = 5.0 / 3.0;
Real State::rho_floor = 1.0e-12;
Real State::ei_floor = 1.0e-20;
Real State::omega = 0.0;
Real State::omega0 = 0.0;
Real State::omega_dot = 0.0;
Real State::driving_rate, State::driving_time;

void State::set_omega(Real o) {
	omega = o;
}

bool State::low_order_variable(int i) {
	return i == pot_index;
}

bool State::smooth_variable(int i) {
	return i == pot_index;
}

Real State::rho_no_floor() const {
	return (*this)[rho_index];
}

Real State::rot_pot(const _3Vec& X) const {
	return rho_no_floor() * rot_phi(X);
}
Real State::rot_phi(const _3Vec& X) const {
	const Real R2 = X[0] * X[0] + X[1] * X[1];
	return -0.5 * R2 * omega * omega;
}

Real State::pot() const {
	return (*this)[pot_index];
}

void State::set_pot(Real a) {
	(*this)[pot_index] = a;
}

State& State::operator=(const State& s) {
	this->Vector<Real, STATE_NF>::operator=(s);
	return *this;
}

State& State::operator=(const Vector<Real, STATE_NF>& s) {
	this->Vector<Real, STATE_NF>::operator=(s);
	return *this;
}

State::State(const State& s) {
	*this = s;
}

Real State::conserved_energy(const _3Vec X) {
	Real egas, egrav, erot;
	egas = et();
	egrav = 0.5 * (pot() - rot_pot(X));
	erot = rot_pot(X);
	return egas + egrav + erot;
}

void State::to_con(const _3Vec& X) {
	(*this)[et_index] += 0.5 * ((*this)[pot_index] + rot_pot(X));
}

void State::from_con(const _3Vec& X) {
	(*this)[et_index] -= 0.5 * ((*this)[pot_index] + rot_pot(X));
}

Real State::phi_eff() const {
	return pot() / rho();
}

Real State::phi(const _3Vec& x) const {
	return (pot() / rho() - rot_pot(x) / (*this).rho_no_floor());
}

Real State::tau() const {
	return (*this)[tau_index];
}

Real State::lz() const {
	return (*this)[lz_index];
}

Real State::set_lz_from_cartesian(const _3Vec& X) {
	set_lz(X[0] * sy() - X[1] * sx());
	return lz();
}

Real State::virial(const _3Vec& x) const {
	Real ek = (vx(x) * vx(x) + vy(x) * vy(x) + vz() * vz()) * 0.5 * rho() + 1.5 * pd();
	Real ep = 0.5 * (pot() - rot_pot(x));
	return 2.0 * ek + ep;
}

void State::set_gamma(Real g) {
	gamma = g;
}

State::~State() {

}

Real State::get_omega() {
	return omega;
}

void State::set_frac(int f, Real a) {
	(*this)[frac_index + f] = a;
}

Real State::frac(int i) const {
	return (*this)[frac_index + i];
}

Vector<Real, STATE_NF> State::x_flux(const _3Vec& X) const {
	Vector<Real, STATE_NF> flux;
	const Real v = vx(X);
	flux = (*this) * v;
	flux[sx_index] += pg(X);
	flux[lz_index] -= X[1] * pg(X);
	flux[et_index] += v * pg(X);
	return flux;
}

Vector<Real, STATE_NF> State::y_flux(const _3Vec& X) const {
	Vector<Real, STATE_NF> flux;
	const Real v = vy(X);
	flux = (*this) * v;
	flux[sy_index] += pg(X);
	flux[lz_index] += X[0] * pg(X);
	flux[et_index] += v * pg(X);
	return flux;
}

Vector<Real, STATE_NF> State::z_flux(const _3Vec& X) const {
	Vector<Real, STATE_NF> flux;
	const Real v = vz();
	flux = (*this) * v;
	flux[sz_index] += pg(X);
	flux[et_index] += v * pg(X);
	return flux;
}

Vector<Real, STATE_NF> State::source(const _3Vec& X, Real t) const {
	Vector<Real, STATE_NF> s = 0.0;
	s[sx_index] += omega * (*this)[sy_index];
	s[sy_index] -= omega * (*this)[sx_index];
	if (omega0 != 0.0) {
		const Real period = 2.0 * M_PI / omega0;
		if (driving_rate != 0.0 && t < driving_time * period) {
			s[lz_index] -= (*this)[lz_index] * driving_rate / period;
		}
	}
	return s;
}

const char* State::field_name(int i) {
	assert(i >= 0);
	assert(i < STATE_NF);
	switch (i) {
	case et_index:
		return "et";
	case sx_index:
		return "sx";
	case sy_index:
		return "sy";
	case lz_index:
		return "lz";
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

void State::enforce_dual_energy_formalism(const _3Vec& X, const State& n1, const State& n2, const State& n3, const State& n4, const State& n5,
		const State& n6) {
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

void State::floor(const _3Vec& X, Real dx) {

	if (et() - ek(X) > 0.1 * et()) {
		(*this)[tau_index] = pow(max(et() - ek(X) - ed(), ei_floor), 1.0 / gamma);
	}
	Real v_x = (*this)[sx_index] / rho();
	Real v_y = (*this)[sy_index] / rho();
	Real v_t, v_r, cos0, sin0, R;
	R = sqrt(X[0] * X[0] + X[1] * X[1]);
	cos0 = X[0] / R;
	sin0 = X[1] / R;
	Real sR = rho() * (v_x * cos0 + v_y * sin0);
	v_t = (*this)[lz_index] / R / rho();
	v_r = sR / rho();
	v_x = v_r * cos0 - v_t * sin0;
	v_y = v_r * sin0 + v_t * cos0;
	Real new_sx = rho() * v_x;
	Real new_sy = rho() * v_y;
	Real w = min(1.0, R / (10.0 * dx));
	(*this)[sx_index] = (1.0 - w) * (*this)[sx_index] + w * new_sx;
	(*this)[sy_index] = (1.0 - w) * (*this)[sy_index] + w * new_sy;
	(*this)[tau_index] = max(pow(ei_floor, 1.0 / gamma), (*this)[tau_index]);
	double total_frac = 0.0;
	for (int i = frac_index; i < frac_index + NFRAC; i++) {
		(*this)[i] /= (*this)[rho_index];
		total_frac += (*this)[i];
	}
	if (total_frac != 0.0) {
		for (int i = frac_index; i < frac_index + NFRAC; i++) {
			(*this)[i] /= total_frac;
		}
	}
	for (int i = frac_index; i < frac_index + NFRAC; i++) {
		(*this)[i] *= (*this)[rho_index];
	}

}

State::State() :
		Vector<Real, STATE_NF>() {
}

State::State(const Vector<Real, STATE_NF>& v) :
		Vector<Real, STATE_NF>(v) {
	return;
}

Real State::rho() const {
	return max(rho_no_floor(), rho_floor);
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
	v_x = (*this)[sx_index] / rho();
	return v_x + X[1] * omega;
}

Real State::vy(const _3Vec& X) const {
	Real v_y;
	v_y = (*this)[sy_index] / rho();
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
	dp_drho = ((8.0 * PhysicalConstants::A) / (3.0 * PhysicalConstants::B)) * x * x / sqrt(x * x + 1.0) + (gamma - 1.0) * ei(X) / rho();
	dp_depsilon = (gamma - 1.0) * rho();
	cs2 = max((pg(X) / (rho() * rho())) * dp_depsilon + dp_drho, 0.0);
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

void State::set_sx(Real a) {
	(*this)[sx_index] = a;
}

void State::set_sy(Real a) {
	(*this)[sy_index] = a;
}

void State::set_lz(Real a) {
	(*this)[lz_index] = a;
}

void State::set_tau(Real a) {
	(*this)[tau_index] = a;
}

void State::set_sz(Real a) {
	(*this)[sz_index] = a;
}

void State::set_et(Real a) {
	(*this)[et_index] = a;
}

Real R(_3Vec x) {
	return sqrt(x[0] * x[0] + x[1] * x[1]);
}

void State::from_prim(const _3Vec& x) {
	double total_frac = 0.0;
	for (int i = frac_index; i < frac_index + NFRAC; i++) {
		total_frac += (*this)[i];
	}
	if (total_frac != 0.0) {
		for (int i = frac_index; i < frac_index + NFRAC; i++) {
			(*this)[i] /= total_frac;
		}
	}
	for (int i = frac_index; i < frac_index + NFRAC; i++) {
		(*this)[i] *= (*this)[rho_index];
	}
	(*this)[sz_index] *= rho();
	(*this)[pot_index] *= rho();
	(*this)[sx_index] -= State::omega * x[1];
	(*this)[sy_index] += State::omega * x[0];
	(*this)[sx_index] *= rho();
	(*this)[sy_index] *= rho();
	(*this)[et_index] += ek(x);
	(*this)[tau_index] = pow((*this)[tau_index], 1.0 / gamma);
	(*this)[lz_index] = (*this)[sy_index] * x[0] - (*this)[sx_index] * x[1];
}

void State::to_prim(const _3Vec& x) {
	double total_frac = 0.0;
	if ((*this)[rho_index] != 0.0) {
		for (int i = frac_index; i < frac_index + NFRAC; i++) {
			(*this)[i] /= (*this)[rho_index];
			total_frac += (*this)[i];
		}
		if (total_frac != 0.0) {
			for (int i = frac_index; i < frac_index + NFRAC; i++) {
				(*this)[i] /= total_frac;
			}
		}
	}
	(*this)[tau_index] = pow((*this)[tau_index], gamma);
	(*this)[et_index] -= ek(x);
	(*this)[sx_index] /= rho();
	(*this)[sy_index] /= rho();
	(*this)[sx_index] += State::omega * x[1];
	(*this)[sy_index] -= State::omega * x[0];
	(*this)[sz_index] /= rho();
	(*this)[pot_index] /= rho();
}
