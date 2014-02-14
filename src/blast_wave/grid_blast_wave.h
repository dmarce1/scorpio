/*
 * grid_blast_wave.h
 *
 *  Created on: Mar 12, 2013
 *      Author: dmarce1
 */

#ifndef GRID_BLAST_WAVE_H_
#define GRID_BLAST_WAVE_H_

#include "../hydro_grid/hydro_grid.h"
#ifdef BLAST_WAVE

class GridBlastWave:
public HydroGrid
{
	static Vector<Real, 4>* radial_avg_tmp;
	static Vector<Real, 4>* radial_avg;
	static int* radial_bin_cnt;
	static int radial_N;
	static int* radial_bin_cnt_tmp;
public:

	Array3d<State, GNX, GNX, GNX> analytic;
	static void run(int, char*[]);
	GridBlastWave();
	virtual ~GridBlastWave();
private:
	Real get_dx() const {
		return HydroGrid::get_dx();
	}
	Real xc(int i) const {
		return HydroGrid::xc(i);
	}
	Real yc(int i) const {
		return HydroGrid::yc(i);
	}
	Real zc(int i) const {
		return HydroGrid::zc(i);
	}
	virtual void set_refine_flags();
	virtual GridBlastWave* new_octnode() const;
	void initialize();
	Real get_output_point(int i, int j, int k, int l) const;
	int nvar_output() const;
	const char* output_field_names(int i) const;
};

#endif

#endif /* GRID_BLAST_WAVE_H_ */
