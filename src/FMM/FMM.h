/*
 * hydro_FMM_grid.h
 *
 *  Created on: Sep 20, 2013
 *      Author: dmarce1
 */

#ifndef HYDRO_FMM_GRID_H_
#define HYDRO_FMM_GRID_H_

#include "../defs.h"

#ifdef USE_FMM

#define USE_FMM_ANGULAR

#include "../hydro/hydro.h"

#include "expansion.h"

#define FORDER 1
#define FNX (INX+2*FBW)

#define FBW (2*FORDER)

#define FSTAGE (15-4)

typedef struct {
	expansion_nodip_t M;
	_3Vec X;
	bool is_leaf;
} multipole_t;

typedef struct {
	expansion_t phi;
	expansion_m1_t g_lz;
	_3Vec X;
} taylor_t;

typedef struct {
	expansion_t M_dot;
} multipole_dot_t;

typedef struct {
	expansion_t phi_dot;
} taylor_dot_t;

typedef struct {
	Real phi;
	_3Vec g;
} _4force_t;

class FMM: public Hydro {
public:
	Array3d<taylor_t, FNX, FNX, FNX> L;
private:
	Array3d<multipole_t, FNX, FNX, FNX> poles;
	Array3d<multipole_dot_t, FNX, FNX, FNX> poles_dot;
	Array3d<taylor_dot_t, FNX, FNX, FNX> L_dot;
	Array3d<Real, GNX, GNX, GNX> old_pot;
	Array3d<_3Vec, FNX / 2, FNX / 2, FNX / 2> Xp;
	Array3d<Real, GNX, GNX, GNX> dpot;
	Array3d<_4force_t, FNX, FNX, FNX> _4force;
	static MPI_Datatype MPI_multipole_t;
	static MPI_Datatype MPI_multipole_dot_t;
	static MPI_Datatype MPI_taylor_t;
	static MPI_Datatype MPI_taylor_dot_t;
	taylor_t* taylor_buffer;
	taylor_dot_t* taylor_dot_buffer;
	static MPI_Datatype MPI_send_bnd_t[26];
	static MPI_Datatype MPI_recv_bnd_t[26];
	static MPI_Datatype MPI_send_bnd_dot_t[26];
	static MPI_Datatype MPI_recv_bnd_dot_t[26];
	multipole_t* moment_buffer;
	multipole_dot_t* moment_dot_buffer;

	static Real d0_array[INX][INX][INX];
	static Real d1_array[2 * INX + 1][2 * INX + 1][2 * INX + 1][3];
	FMM* neighbors[26];
	typedef void (FMM::*ifunc_t)(int);
	static ifunc_t cs[FSTAGE + 1];
	static ifunc_t cs_dot[FSTAGE + 1];
	static ifunc_t cs_children[5];
	static MPI_Datatype MPI_comm_child_poles_t[8];
	static MPI_Datatype MPI_comm_taylor_t[8];
	static MPI_Datatype MPI_comm_child_poles_dot_t[8];
	static MPI_Datatype MPI_comm_taylor_dot_t[8];
	static MPI_Datatype MPI_comm_child3_t[8];
	static MPI_Datatype MPI_4force_t;
	MPI_Request send_request[26], recv_request[26];
	_4force_t* _4force_buffer;
	Vector<Real, 3> mom_sum;
	bool is_leaf(int i, int j, int k) const;

public:
	static bool solve_on;
	static Vector<Real, 6> momentum_sum();
	static Vector<Real, 4> com_sum();
	Array3d<Real, GNX, GNX, GNX> drho_dt;
	_3Vec g_grav(int i, int j, int k ) const {
		_3Vec g;
		g[0] = gx(i,j,k);
		g[1] = gy(i,j,k);
		g[2] = gz(i,j,k);
		return g;
	}
	_3Vec geff(int i, int j, int k ) const {
		_3Vec g;
		g[0] = gx(i,j,k) + pow(State::get_omega(),2)*xc(i);
		g[1] = gy(i,j,k) + pow(State::get_omega(),2)*yc(j);
		g[2] = gz(i,j,k);
		return g;
	}
	Real gx(int i, int j, int k) const;
#ifdef USE_FMM_ANGULAR
	Real g_lz(int i, int j, int k) const;
#endif
	static _3Vec system_com();
		Real gy(int i, int j, int k) const;
	Real gz(int i, int j, int k) const;
	Real g_energy(int i, int j, int k) const;
	virtual void compute_update(int dir);
	//	virtual void compute_dudt(int dir);
	static void step(Real dt);
    Real set_phi(int i, int j, int k,double);
    Real get_phi(int i, int j, int k) const;
	Real get_dphi_dt(int i, int j, int k) const;
	Real get_drho_dt(int i, int j, int k) const;
	static void FMM_from_children();
	static void FMM_solve();
	static void FMM_solve_dot();
	static void MPI_datatypes_init();
	virtual void allocate_arrays();
	virtual void deallocate_arrays();
	virtual FMM* new_octnode() const;
	virtual void initialize();
	FMM();
	virtual ~FMM();
	void find_neighbors();

	void moments_recv(int);
	void moments_recv_dot(int);

	void moments_recv_wait(int);
	void _4force_recv(int);
	void _4force_recv_wait(int);
	void moments_send_dot(int);
	void moments_send(int);
	void moments_send_wait(int);
	void moments_send_wait_dot(int);
	void _4force_send(int);
	void _4force_send_wait(int);
	void moments_communicate_all(int);
	void moments_communicate_all_dot(int);
	void moments_communicate_wait_all(int);
	void compute_interactions(int);
	void compute_interactions_dot(int);
	void expansion_recv(int);
	void expansion_recv_wait(int);
	void expansion_send(int);
	void expansion_recv_dot(int);
	void expansion_recv_wait_dot(int);
	void expansion_send_dot(int);
	void expansion_send_wait(int);
	void null(int);
	static bool check_for_refine();
	static void store_pot();
	static void update();
	static void account_pot();
	static void to_conserved_energy();
	static void from_conserved_energy();
	static void pot_to_hydro_grid();

	static Real get_phi_at(Real x, Real y, Real z) {
		Real p, tmp;
		FMM* g;
		p = 0.0;
		for (int l = 0; l < get_local_node_cnt(); l++) {
			g = dynamic_cast<FMM*>(get_local_node(l));
			if (MPI_rank() == g->proc()) {
				for (int k = FBW; k < FNX - FBW; k++) {
					if (z >= g->zf(k + BW - FBW) && z < g->zf(k + 1 + BW - FBW)) {
						for (int j = FBW; j < FNX - FBW; j++) {
							if (y >= g->yf(j + BW - FBW) && y < g->yf(j + 1 + BW - FBW)) {
								for (int i = FBW; i < FNX - FBW; i++) {
									if (x >= g->xf(i + BW - FBW) && x < g->xf(i + 1 + BW - FBW)) {
										if (!g->zone_is_refined(i + BW - FBW, j + BW - FBW, k + BW - FBW)) {
											//		p = g->get_phi(i + BW - FBW, j + BW - FBW, k + BW - FBW);
											Real c0 = g->L(i, j, k).phi();
									//		Real c1 = g->L(i, j, k).phi(0);
									//		Real dx = (x - g->xc(i + BW - FBW));
											p = c0;// + c1 * dx;
										}
									}
								}
							}
						}
					}
				}
			}
		}
		tmp = p;
		MPI_Allreduce(&tmp, &p, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );
		return p;
	}

};

#endif /* HYDRO_FMM_GRID_H_ */
#endif
