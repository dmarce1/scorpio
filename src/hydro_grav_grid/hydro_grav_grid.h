/*
 * hydo_grav_grid.h
 *
 *  Created on: Mar 20, 2013
 *      Author: dmarce1
 */

#ifndef HYDO_GRAV_GRID_H_
#define HYDO_GRAV_GRID_H_

#include "../defs.h"
#include "../hydro_grid/hydro_grid.h"
#include "../poisson/poisson.h"
#include "../reconstruct.h"
#include "../indexer3d.h"
#include "../physical_constants.h"

#ifdef HYDRO_GRAV_GRID

class HydroGravGrid: public HydroGrid, public Poisson {
public:
	virtual void init();
	static void run(int argc, char* argv[]);
private:
	static Real poisson_tolerance;
	Reconstruct reconstruct;
	static void to_conserved_energy();
	static void from_conserved_energy();
	static void pot_from_hydro_grid();
	static void set_gravity_source();
	virtual void initialize();
	virtual void compute_dudt(int);
	virtual Real xf(int i) const;
	virtual Real yf(int i) const;
	virtual Real zf(int i) const;
	virtual void set_refine_flags();
	virtual HydroGravGrid* new_octnode() const;
	virtual void deallocate_arrays();
	virtual void allocate_arrays();
	virtual Real get_output_point(int i, int j, int k, int l) const;
	virtual const char* output_field_names(int i) const;
	virtual int nvar_output() const;
	//virtual void max_dt_compute(int);
protected:
    static void pot_to_hydro_grid();
    virtual void create_child(const ChildIndex& c);
	static void set_origin(const _3Vec& o) {
		HydroGrid::set_origin(o);
		MultiGrid::set_origin(o);
	}
	static void solve_poisson();
	virtual void write(FILE* fp) const;
	virtual void read(FILE* fp);
	static void step(Real);
	static void setup_grid_structure();
	static bool check_for_refine();
	static Real hydro_time, poisson_boundary_time, poisson_interior_time;
	Real get_dx() const;
	static void set_poisson_tolerance(Real);
	virtual Real xc(int i) const;
	virtual Real yc(int i) const;
	virtual Real zc(int i) const;
public:

	void send_gforces() {
		HydroGravGrid* sib;
		Vector<Real, 3> dx;
		Real r;
		Array3d<Real, GNX, GNX, GNX> f(true);
		for (int i = 0; i < 6; i++) {
			sib = dynamic_cast<HydroGravGrid*>(get_sibling(i));
			if (sib != NULL) {
				for (Indexer3d i(BW, GNX - BW - 1); !i.end(); i++) {
					f(i) = 0.0;
					for (Indexer3d j(BW, GNX - BW - 1); !j.end(); j++) {
						dx = this->X(i) - sib->X(j);
						r = sqrt(dx.dot(dx));
						f(i) -= (*this)(j).rho() / r;
					}
				}
			}
		}
	}

	_3Vec gforce(int i, int j, int k) const {
		_3Vec g;
		g[0] = gx(i, j, k);
		g[1] = gy(i, j, k);
		g[2] = gz(i, j, k);
		return g;
	}
	Real gx(int i, int j, int k) const {
//		const int o = BW - 1;
//		return -(phi_fx(i + 1 - o, j - o, k - o) - phi_fx(i - o, j - o, k - o)) / get_dx();
		return 0.5 * (get_fx(i - BW + 1, j - BW + 1, k - BW + 1) + get_fx(i - BW + 2, j - BW + 1, k - BW + 1));
	}
	Real gy(int i, int j, int k) const {
//		const int o = BW - 1;
//		return -(phi_fy(i - o, j + 1 - o, k - o) - phi_fy(i - o, j - o, k - o)) / get_dx();
		return 0.5 * (get_fy(i - BW + 1, j - BW + 1, k - BW + 1) + get_fy(i - BW + 1, j - BW + 2, k - BW + 1));
	}
	Real gz(int i, int j, int k) const {
//		const int o = BW - 1;
//		return -(phi_fz(i - o, j - o, k + 1 - o) - phi_fz(i - o, j - o, k - o)) / get_dx();
		return 0.5 * (get_fz(i - BW + 1, j - BW + 1, k - BW + 1) + get_fz(i - BW + 1, j - BW + 1, k - BW + 2));
	}
	HydroGravGrid();
	virtual ~HydroGravGrid();
};

#endif /* HYDO_GRAV_GRID_H_ */

#endif
