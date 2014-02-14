#include "../defs.h"
#ifdef ROTATING_DISC

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "rotating_disc.h"
#include "state.h"


const Real State::gamma = EULER_GAMMA;
const int State::d_index = 0;
const int State::sx_index = 1;
const int State::sy_index = 2;
const int State::sz_index = 3;
const int State::et_index = 4;
const int State::tau_index = 5;
const int State::sr_index = 6;
const int State::lz_index = 7;

const Real M_c = M_C;
const Real G = BIG_G;

//3-d torus
//calculates x force of gravity
static double fx(double x, double y, double z) {
	const Real r_cyl = sqrt(x * x + y * y);
	const Real r = sqrt(x * x + y * y + z * z);
	const Real F = -G * M_c / (r * r);
	return F * (x / r);
}
//calculates y force of gravity
static double fy(double x, double y, double z) {
	const Real r_cyl = sqrt(x * x + y * y);
	const Real r = sqrt(x * x + y * y + z * z);
	const Real F = -G * M_c / (r * r);
	return F * (y / r);
}
//calculates z force of gravity
static double fz(double x, double y, double z) {
	const Real r_cyl = sqrt(x * x + y * y);
	const Real r = sqrt(x * x + y * y + z * z);
	const Real F = -G * M_c / (r * r);
	return F * (z / r);
}

// Calculates cylindrical R force of gravity
static double fr(double x, double y, double z) {
	const Real r_cyl = sqrt(x * x + y * y);
	const Real r = sqrt(x * x + y * y + z * z);
	const Real F = -G * M_c / (r * r);
	return F * (r_cyl / r);
}

// This returns cylindrical R.
Real State::R(const _3Vec& X) {
	return sqrt(X[0] * X[0] + X[1] * X[1]);
}

Vector<Real, STATE_NF> State::source(const _3Vec& X) const {
	State s = Vector<Real, STATE_NF>(0.0);

	// This is the radial momentum source term, indep. of gravity.
	// Is it independent of rotation?
	s[sr_index] += (pg(X) + pow(lz() / R(X), 2) / rho()) / R(X);

	// Add half of the Coriolis force for rotating cartesian
	// momentum.
	double Omega = dynamic_cast<RotatingDisc*>(OctNode::get_root())->get_frame_omega();
	s[sx_index] += sy() * Omega;
	s[sy_index] -= sx() * Omega;

	// There won't be any gravity within a certain radius of the center.
	double x_in = dynamic_cast<RotatingDisc*>(OctNode::get_root())->get_x_in();
	double R_inner = x_in * R_OUTER;
	if (R(X) < 0.5 * R_inner)
		return s;

	s[sx_index] += rho() * fx(X[0], X[1], X[2]); //x component of gravity
	s[sy_index] += rho() * fy(X[0], X[1], X[2]); //y component of gravity
	s[sz_index] += rho() * fz(X[0], X[1], X[2]); //z component of gravity

	s[sr_index] += rho() * fr(X[0], X[1], X[2]); // r component of gravity
	return s;
}

const char* State::field_name(int i) {
	static char dstr[3];
	dstr[2] = '\0';
	dstr[0] = 'd';
	assert(i >= 0);
	assert(i < STATE_NF);
	switch (i) {
	case d_index:
		return "d";
	case sz_index:
		return "sz";
	case et_index:
		return "et";
	case sr_index:
		return "sr";
	case lz_index:
		return "lz";
	case sx_index:
		return "sx";
	case sy_index:
		return "sy";
	}
	return "tau";
}

void State::enforce_outflow(const _3Vec& X, const OctFace& f) {
	switch (f) {
	case XU:
		if (vx(X) > 0.0) {
			set_et(et() - 0.5 * sx() * sx() / rho());
			set_sx(0.0);
			(*this)[lz_index] = X[0] * vy(X) * rho();
			(*this)[sr_index] = X[1] * vy(X) * rho() / R(X);
		}
		break;
	case XL:
		if (vx(X) < 0.0) {
			set_et(et() - 0.5 * sx() * sx() / rho());
			set_sx(0.0);
			(*this)[lz_index] = X[0] * vy(X) * rho();
			(*this)[sr_index] = X[1] * vy(X) * rho() / R(X);
		}
		break;
	case YU:
		if (vy(X) > 0.0) {
			set_et(et() - 0.5 * sy() * sy() / rho());
			set_sy(0.0);
			(*this)[lz_index] = -X[1] * vx(X) * rho();
			(*this)[sr_index] = X[0] * vx(X) * rho() / R(X);
		}
		break;
	case YL:
		if (vy(X) < 0.0) {
			set_et(et() - 0.5 * sy() * sy() / rho());
			set_sy(0.0);
			(*this)[lz_index] = -X[1] * vx(X) * rho();
			(*this)[sr_index] = X[0] * vx(X) * rho() / R(X);
		}
		break;
	case ZU:
		if (sz() > 0.0) {
			set_et(et() - 0.5 * sz() * sz() / rho());
			set_sz(0.0);
		}
		break;
	case ZL:
		if (sz() < 0.0) {
			set_et(et() - 0.5 * sz() * sz() / rho());
			set_sz(0.0);
		}
		break;
		(*this)[et_index] += 0.5 * (sr(X) * sr(X) + st(X) * st(X) + sz() * sz()) / rho();
	}
}

void State::floor(const _3Vec& X) {
	const double rho_floor = dynamic_cast<RotatingDisc*>(OctNode::get_root())->get_density_floor();
	const double ei_floor = dynamic_cast<RotatingDisc*>(OctNode::get_root())->get_density_floor();
	(*this)[d_index] = max((*this)[d_index], rho_floor);
	const Real ei0 = et() - ek(X);
	if (ei0 > 0.1 * et()) {
		(*this)[tau_index] = pow(max(ei0, ei_floor), 1.0 / gamma);
	}

	// Floor everything in the center of the grid.

	double x_in = dynamic_cast<RotatingDisc*>(OctNode::get_root())->get_x_in();
	double R_inner = x_in * R_OUTER;
	if (R(X) < 0.5 * R_inner) {
		(*this)[d_index] = rho_floor;
		(*this)[sx_index] = 0.0;
		(*this)[sy_index] = 0.0;
		(*this)[sz_index] = 0.0;
		(*this)[lz_index] = 0.0;
		(*this)[sr_index] = 0.0;
	}

	bool ang_mom_cons = dynamic_cast<RotatingDisc*>(OctNode::get_root())->get_ang_mom_cons();
	// Redefine lz from sx and sy when
	// not conserving angular momentum.
	if (!ang_mom_cons) {
		(*this)[lz_index] = sy() * X[0] - sx() * X[1];
	}

}

void State::reflect_on_z() {
	set_sz(-sz());
}

Real State::poisson_source() const {
	Real a = 4.0 * M_PI * G_CON * rho();
	return a;
}

Real State::refine_value() const {
	return rho();
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

Real State::lz() const {
	return (*this)[lz_index];
}

Real State::sz() const {
	return (*this)[sz_index];
}

Real State::et() const {
	return (*this)[et_index];
}

_3Vec State::V(const _3Vec& X) const {
	_3Vec v;
	v[0] = vx(X);
	v[1] = vy(X);
	v[2] = vz(X);
	return v;
}

Real State::vx(const _3Vec& X) const {
	return (X[0] * sr(X) - X[1] * st(X)) / R(X) / rho();
}

Real State::vy(const _3Vec& X) const {
	return (X[1] * sr(X) + X[0] * st(X)) / R(X) / rho();
}

Real State::vz(const _3Vec&) const {
	return sz() / rho();
}

Real State::ek(const _3Vec& X) const {
	return 0.5 * (sr(X) * sr(X) + st(X) * st(X) + sz() * sz()) / rho();
}

Real State::ei(const _3Vec& X) const {
	const double ei_floor = dynamic_cast<RotatingDisc*>(OctNode::get_root())->get_density_floor();
	return max(et() - ek(X), ei_floor);
}
//this is where pressure is calculated.
// i am changing it to just use a straight polytropic eq of state
Real State::pg(const _3Vec& X) const {

	return KAPPA * pow(rho(), gamma);

	//	const Real ei0 = et() - ek();
	//	if (ei0 < 0.001 * et()) {
	//		return (gamma - 1.0) * pow((*this)[tau_index], gamma);
	//	} else {
	//		return (gamma - 1.0) * ei0;
	//	}
}

Real State::cs(const _3Vec& X) const {
	assert(rho() > 0.0);
	assert(pg(X) >= 0.0);
	return sqrt(gamma * pg(X) / rho());
}

Real State::st(const _3Vec& X) const {
	bool ang_mom_cons = dynamic_cast<RotatingDisc*>(OctNode::get_root())->get_ang_mom_cons();
	double Omega = dynamic_cast<RotatingDisc*>(OctNode::get_root())->get_frame_omega();
	if (ang_mom_cons) {
		return lz() / R(X) - rho() * R(X) * Omega;
	} else {
		return sy() * X[0] / R(X) - sx() * X[1] / R(X) - rho() * R(X) * Omega;
	}
}

Real State::sr(const _3Vec& X) const {
	bool ang_mom_cons = dynamic_cast<RotatingDisc*>(OctNode::get_root())->get_ang_mom_cons();
	if (ang_mom_cons) {
		return (*this)[sr_index];
	} else {
		return sx() * X[0] / R(X) + sy() * X[1] / R(X);
	}
}

Real State::max_abs_x_eigen(const _3Vec& X) const {
	return fabs(vx(X)) + cs(X);
}

Real State::max_abs_y_eigen(const _3Vec& X) const {
	return fabs(vy(X)) + cs(X);
}

Real State::max_abs_z_eigen(const _3Vec& X) const {
	return fabs(vz(X)) + cs(X);
}

Vector<Real, STATE_NF> State::x_flux(const _3Vec& X) const {
	Vector<Real, STATE_NF> flux = 0.0;
	Real v, p, r;
	v = vx(X);
	p = pg(X);
	r = R(X);
	for (int i = 0; i < STATE_NF; i++) {
		flux[i] = (*this)[i] * v;
	}
	flux[sx_index] += p;
	flux[et_index] += v * p;
	flux[sr_index] += X[0] * p / r;
	flux[lz_index] -= X[1] * p;
	return flux;
}

Vector<Real, STATE_NF> State::y_flux(const _3Vec& X) const {
	Vector<Real, STATE_NF> flux = 0.0;
	Real v, p, r;
	v = vy(X);
	p = pg(X);
	r = R(X);
	for (int i = 0; i < STATE_NF; i++) {
		flux[i] = (*this)[i] * v;
	}
	flux[sy_index] += p;
	flux[et_index] += v * p;
	flux[sr_index] += X[1] * p / r;
	flux[lz_index] += X[0] * p;
	return flux;
}

Vector<Real, STATE_NF> State::z_flux(const _3Vec& X) const {
	Vector<Real, STATE_NF> flux = 0.0;
	Real v, p;
	v = vz(X);
	p = pg(X);
	for (int i = 0; i < STATE_NF; i++) {
		flux[i] = (*this)[i] * v;
	}
	flux[sz_index] += p;
	flux[et_index] += v * p;
	return flux;
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

void State::set_lz(Real a) {
	(*this)[lz_index] = a;
}

void State::set_sr(Real a) {
	(*this)[sr_index] = a;
}

void State::set_tau(Real a) {
	(*this)[tau_index] = a;
}

void State::to_prim(const _3Vec& X) {
	double Omega = dynamic_cast<RotatingDisc*>(OctNode::get_root())->get_frame_omega();
	(*this)[et_index] -= ek(X);
	(*this)[sx_index] /= rho();
	(*this)[sy_index] /= rho();
	(*this)[sz_index] /= rho();
	(*this)[sr_index] /= rho();
	(*this)[lz_index] /= rho() * R(X);
	(*this)[lz_index] -= Omega * R(X);
}

void State::to_con(const _3Vec& X) {
	double Omega = dynamic_cast<RotatingDisc*>(OctNode::get_root())->get_frame_omega();
	(*this)[sx_index] *= rho();
	(*this)[sy_index] *= rho();
	(*this)[sz_index] *= rho();
	(*this)[et_index] += ek(X);
	(*this)[sr_index] *= rho();
	(*this)[lz_index] += Omega * R(X);
	(*this)[lz_index] *= rho() * R(X);
}

#endif
