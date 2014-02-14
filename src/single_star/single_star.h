/*
 * single_star.h
 *
 *  Created on: Mar 20, 2013
 *      Author: dmarce1
 */

#ifndef SINGLE_STAR_H_
#define SINGLE_STAR_H_

#include "../defs.h"

#ifdef SINGLE_STAR

#include "../hydro_grav_grid/hydro_grav_grid.h"
#include "../hydro_FMM_grid/hydro_FMM_grid.h"

#ifdef HYDRO_GRAV_GRID

#ifdef USE_FMM
class SingleStar: public HydroFMMGrid {
#else
	class SingleStar: public HydroGravGrid {
#endif
private:
	static Vector<Real, 2>* radial_avg_tmp;
	static Vector<Real, 2>* radial_avg;
	static int* radial_bin_cnt;
	static int radial_N;
	static int* radial_bin_cnt_tmp;
	virtual void set_refine_flags();
	Real radius(int i, int j, int k);
public:
	static void read_from_file(const char* str, int* i1, int* i2);
	static void write_to_file(int i1, int i2, const char* idname);
	virtual Real get_output_point(int i, int j, int k, int l) const {
		State U = (*this)(i, j, k);
		_3Vec x = X(i, j, k);
		switch (l) {
		case 0:
			return U.rho();
		case 1:
			return U.frac(1);
		case 2:
			return U.frac(0);
		case 3:
			return U.vx(x);
		case 4:
			return U.vy(x);
		case 5:
			return U.vz();
		case 6:
			return U.temp(x);
		case 7:
			return U.ei(x);
		case 8:
			return D(i,j,k)[State::et_index];
		case 9:
			return U.pg(x)-U.pd();
		case 10:
			return U.pd();

#ifdef USE_FMM
		case 11:
			return get_phi(i, j, k);
		case 12:
			return gx(i, j, k);
		case 13:
			return gy(i, j, k);
		case 14:
			return gz(i, j, k);
#else
			default:
			return get_phi(i - BW + 1, j - BW + 1, k - BW + 1);
#endif
		}
		return 0.0;
	}
	virtual const char* output_field_names(int i) const {
		switch (i) {
		case 0:
			return "rho";
		case 1:
			return "He";
		case 2:
			return "CO";
		case 3:
			return "vx";
		case 4:
			return "vy";
		case 5:
			return "vz";
		case 6:
			return "T";
		case 7:
			return "ei";
		case 8:
			return "de";
		case 9:
			return "pi";
		case 10:
			return "pd";
		case 11:
			return "phi";
		case 12:
			return "gx";
		case 13:
			return "gy";
		case 14:
			return "gz";
		}
		return "";
	}
	virtual int nvar_output() const {
#ifdef USE_FMM
		return 15;
#else
		return 12;
#endif
	}
    static void scf_run();
    static void run(int, char*[]);
	virtual SingleStar* new_octnode() const {
		return new SingleStar;
	}
	virtual void initialize();
	SingleStar();
	virtual ~SingleStar();
};
#endif

#endif /* SINGLE_STAR_H_ */
#endif
