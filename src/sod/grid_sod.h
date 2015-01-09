/*
 * grid_sod.h
 *
 *  Created on: Mar 12, 2013
 *      Author: dmarce1
 */

#ifndef GRID_SOD_H_
#define GRID_SOD_H_

#include "../hydro/hydro.h"

#ifdef SOD

class GridSod:
public Hydro
{
    static Vector<Real, 5>* line_avg_tmp;
    static Vector<Real, 5>* line_avg;
    static int* line_bin_cnt;
    static int line_N;
    static int* line_bin_cnt_tmp;
    static bool oblique;
public:

    Array3d<State, GNX, GNX, GNX> analytic;
    static void run(int, char*[]);
    GridSod();
    virtual ~GridSod();
private:
    Real get_dx() const {
        return Hydro::get_dx();
    }
    Real xc(int i) const {
        return Hydro::xc(i);
    }
    Real yc(int i) const {
        return Hydro::yc(i);
    }
    Real zc(int i) const {
        return Hydro::zc(i);
    }
    Real yzsod(Vector<int, 3> i) const {
        return yzsod(i[0], i[1], i[2]);
    }
    Real xsod(Vector<int, 3> i) const {
        if (oblique) {
            return xsod(i[0], i[1], i[2]);
        } else {
            return X(i)[0];
        }
    }
    Real yzsod(int i, int j, int k) const {
        return sqrt(X(i, j, k).dot(X(i, j, k)) - xsod(i, j, k) * xsod(i, j, k));
    }
    Real xsod(int i, int j, int k) const {
        Real x = 0.0;
        if (oblique) {
            x += X(i, j, k)[0] * (sqrt(3) / 3.0);
            x += X(i, j, k)[1] * (sqrt(3) / 3.0);
            x += X(i, j, k)[2] * (sqrt(3) / 3.0);
        } else {
            x = xc(i);
        }
        return x;

    }
    virtual void flux_physical_bounds(int dir);
    virtual void physical_boundary(int dir);
    virtual void set_refine_flags();
    virtual GridSod* new_octnode() const;
    void initialize();
    Real get_output_point(int i, int j, int k, int l) const;
    int nvar_output() const;
    const char* output_field_names(int i) const;
};

#endif

#endif /* GRID_SOD_H_ */
