#ifndef EULER0_STATE_H_
#define EULER0_STATE_H_

#include "../defs.h"
#include "../real.h"
#include "state.h"
#include "../vector.h"
#include "../oct_node/oct_face.h"

#ifndef NRHO
#define NRHO 1
#endif

#define G_CON (1.0)

#define STATE_NF (8) //Changed to add angular momentum.

class State: public Vector<Real, STATE_NF> {
public:
	static const Real gamma;
	static const Real ei_floor;
	static const Real rho_floor;
	static const int d_index;
	static const int sx_index;
	static const int sy_index;
	static const int sz_index;
	static const int et_index;
	static const int tau_index;
        static const int lz_index;
        static const int sr_index;
	State();
	State(const Vector<Real, STATE_NF>&);
	virtual Real refine_value() const;
	virtual ~State() {
		return;
	}
	virtual void initialize() {
	}
	static bool low_order_variable(int i) {
		return false;
	}
	static bool smooth_variable(int i) {
		return false;
	}
	virtual Vector<Real, STATE_NF> source(const _3Vec& X) const;
	virtual Real max_abs_x_eigen(const _3Vec& X) const;
	virtual Real max_abs_y_eigen(const _3Vec& X) const;
	virtual Real max_abs_z_eigen(const _3Vec& X) const;
	virtual Vector<Real, STATE_NF> x_flux(const _3Vec& X) const;
	virtual Vector<Real, STATE_NF> y_flux(const _3Vec& X) const;
	virtual Vector<Real, STATE_NF> z_flux(const _3Vec& X) const;
	virtual Real poisson_source() const;
	virtual void to_prim(const _3Vec& X);
	virtual void to_con(const _3Vec& X);
	static const char* field_name(int i);
	virtual void floor(const _3Vec&);
	Real rho() const;
	Real tau() const;
	Real sz() const;
	Real sx() const;
	Real sy() const;
        Real lz() const;
        Real st(const _3Vec& X) const;
        Real sr(const _3Vec& X) const;
	Real vx(const _3Vec& X) const;
	Real vy(const _3Vec& X) const;
	Real vz(const _3Vec& X) const;
	Real ek(const _3Vec& X) const;
        //	_3Vec V() const {
        //		_3Vec v;
        //		v[0] = vx(X);
        //		v[1] = vy(X);
        //		v[2] = vz(X);
        //		return v;
        //	}
        _3Vec V(const _3Vec& X) const;
        static Real R(const _3Vec& X);
	Real et() const;
	Real ei(const _3Vec& X) const;
	Real cs(const _3Vec& X) const;
	Real pg(const _3Vec& X) const;
	Real pot() const;
	Real get_rho(int i) const {
		if (i < 0) {
			return rho();
		} else {
			return (*this)[d_index + i];
		}
	}
	void set_pot(Real);
	void set_rho(Real, int = 0);
	void set_sx(Real);
	void set_sy(Real);
	void set_sz(Real);
	void set_sr(Real);
	void set_lz(Real);
	void set_et(Real);
	void set_tau(Real);
	void reflect_on_z();
	virtual void enforce_outflow(const _3Vec&,const OctFace&);
};

#endif /* EULER_STATE_H_ */
