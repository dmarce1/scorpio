#ifndef GRIDPOISSON_H_
#define GRIDPOISSON_H_

#include "../array3d.h"
#include "../oct_node/oct_node.h"
#include "../oct_node/oct_face.h"
#include "../virtual_process.h"
#include <mpi.h>


#define COMPUTE_STAGE_MAX 43
#define MG_ES_STAGE_MAX 6



class MultiGrid: virtual public OctNode, public VirtualProcess<MultiGrid> {
private:
	static _3Vec origin;
	static ifunc_t cs[COMPUTE_STAGE_MAX + 1];
	static ifunc_t es[MG_ES_STAGE_MAX + 1];
//	static ifunc_t es[MG_ES_STAGE_MAX + 1];
	static MPI_Datatype MPI_send_bnd1_t[6];
	static MPI_Datatype MPI_recv_bnd1_t[6];
	static MPI_Datatype MPI_send_bnd2_t[6];
	static MPI_Datatype MPI_recv_bnd2_t[6];
	static MPI_Datatype MPI_send_bnd3_t[6];
	static MPI_Datatype MPI_recv_bnd3_t[6];
	static MPI_Datatype MPI_comm_child_t[8];
	static MPI_Datatype MPI_recv_amr_child_t[8];
	Array3d<Real, PNX, PNX, PNX> phi0;
	int cx, ax;
	Real dx;
	Vector<int, 3> offset;
	static void MPI_datatypes_init();
	static void redistribute_grids();
	virtual void initialize() = 0;
	virtual void set_refine_flags();
	virtual MultiGrid* new_octnode() const = 0;

	void vdown_init_recv(int);
	void vdown_init_recv_wait(int);
	virtual void vdown_init_compute(int);
	virtual void vdown_init_adjust_send(int);
	void vdown_init_adjust_send_wait(int);
	void vdown_init_adjust_recv(int);
	virtual void vdown_init_adjust_recv_wait(int);
	void vdown_init_store_dphi(int);
	void vdown_init_send(int);
	void vdown_init_send_wait(int);
	void retire_dphi(int);
	void phi_children_recv(int);
	void phi_children_send(int);
	void phi_children_send_wait(int);
	void phi_children_recv_wait(int);
	void vup_init_recv_wait(int);
	void vup_init_recv(int);
	void relax_communicate(int);
	void relax_wait(int);
	void relax_begin_down_loop(int);
	void relax_end_down_loop(int);
	void relax_begin_up_loop(int);
	void relax_end_up_loop(int);
	virtual void relax_compute(int);
	void vup_init_send_wait(int);
	void vup_init_send(int);
	void dphi_amr_boundary_communicate(int);
	void dphi_amr_boundary_wait(int);
	void phi_real_boundary_communicate(int);
	void phi_real_boundary_begin_loop(int);
	void phi_real_boundary_end_loop(int);
	void phi_real_boundary_wait(int);
	void phi_amr_boundary_communicate(int);
	void phi_amr_boundary_wait(int);

	virtual void forces_compute(int);
	void forces_adjust_send(int);
	void forces_adjust_send_wait(int);
	void forces_adjust_recv(int);
	void forces_adjust_recv_wait(int);

	/*	void phi_face_compute(int);
	 void phi_face_adjust_send(int);
	 void phi_face_adjust_send_wait(int);
	 void phi_face_adjust_recv(int);
	 void phi_face_adjust_recv_wait(int);*/

	virtual void residual_error_compute(int);
	//void compute_forces(int);

	void restart_compute();
	virtual void create_child(const ChildIndex& c);
protected:
	Real st0;
	void null(int);
	int amr_cnt;
	MPI_Request amr_request[24];
	int amr_id[24];
	Vector<int, 3> amr_ub[24], amr_lb[24];
	int amr_child_proc[24], amr_face[24];
	MPI_Request send_request[8], recv_request[8];
	Real* mpi_buffer[24];
	Array3d<Real, PNX, PNX, PNX> S;
	Array3d<Real, PNX, PNX, PNX> phi;
	Array3d<Real, PNX, PNX, PNX> fx;
	Array3d<Real, PNX, PNX, PNX> fy;
	Array3d<Real, PNX, PNX, PNX> fz;
	Array3d<Real, PNX, PNX, PNX> dphi1;
	Array3d<Real, PNX, PNX, PNX> dphi;
	static void set_origin(const _3Vec& o) {
		origin = o;
	}
	virtual void write(FILE* fp) const;
	int guard_proc(OctFace f) const;
	virtual void read(FILE* fp);
	virtual void inject_from_parent(ChildIndex);
	virtual void create_multigrid_child(const ChildIndex& c);
	virtual Real get_output_point(int i, int j, int k, int l) const;
	virtual const char* output_field_names(int i) const;
	static Real vcycle();
	static void boundary_communicate();
//	static void compute_faces();
	static bool check_for_refine();
	virtual void deallocate_arrays();
	virtual void allocate_arrays();
	virtual Real xf(int) const;
	virtual int nvar_output() const;
	virtual Real yf(int) const;
	virtual Real zf(int) const;
	Real get_dx() const;
	Real get_dphi(int i, int j, int k) const {
		return dphi(i, j, k);
	}
	void set_dphi(int i, int j, int k, Real d) {
		phi(i, j, k) = d;
	}
	Real get_dphi1(int i, int j, int k) const {
		return dphi1(i, j, k);
	}
	Real get_fx(int i, int j, int k) const;
	Real get_fy(int i, int j, int k) const;
	Real get_fz(int i, int j, int k) const;
	void set_fx(int i, int j, int k, Real a) {
		fx(i, j, k) = a;
	}
	void set_fy(int i, int j, int k, Real a) {
		fy(i, j, k) = a;
	}
	void set_fz(int i, int j, int k, Real a) {
		fz(i, j, k) = a;
	}
	Real get_residual(int i, int j, int k) const;
	Real get_source(int i, int j, int k) const;
	bool poisson_zone_is_refined(int, int, int) const;
public:
	void mult_dx(Real d) {
		if (get_level() == 0) {
			h0 *= d;
			origin *= d;
		}
		dx *= d;
		for (int i = 0; i < 8; i++) {
			if (get_child(i)) {
				dynamic_cast<MultiGrid*>(get_child(i))->MultiGrid::mult_dx(d);
			}
		}
	}
	static Real h0;
	static Real get_phi_at(Real x, Real y, Real z);
	virtual Real xc(int) const;
	virtual Real yc(int) const;
	virtual Real zc(int) const;
	void phi_calc_amr_bounds();
	void set_phi(int, int, int, Real);
	void set_source(int i, int j, int k, Real s);
	Real get_phi(int i, int j, int k) const;
	Real* get_phi_ptr(int i, int j, int k);
	virtual ~MultiGrid();
	virtual void init();
	MultiGrid();
};

#endif /* GRIDPOISSON_H_ */
