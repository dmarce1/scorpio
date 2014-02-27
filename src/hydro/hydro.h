#ifndef GRID_NODE_H_
#define GRID_NODE_H_

#include "../oct_node/oct_node.h"
#include "../reconstruct.h"
#include "../virtual_process.h"
#include <stdlib.h>

#ifdef USE_HYDRO_GRID

#define GRID_CS_SIZE 18
#define GRID_ES_SIZE 14

#ifndef GNX
#define GNX 8
#endif

#ifndef BW
#define BW 2
#endif

class Hydro: public virtual OctNode, public VirtualProcess<Hydro> {
protected:
    static void floor_density();
private:
    static _3Vec origin;
    static bool shadow;
    static ifunc_t cs[GRID_CS_SIZE];
    static ifunc_t cs_children[4];
    static ifunc_t es[GRID_ES_SIZE];
    static bool initialized;
    static MPI_Datatype MPI_interior_t;
    static MPI_Datatype MPI_state_t;
    static MPI_Datatype MPI_guard_send_t[6];
    static MPI_Datatype MPI_guard_recv_t[6];
    static MPI_Datatype MPI_child_t[8];
public:
    static Real _dt, _beta;
private:
    Real time;
    static State FO0;
    static State FO;
    Array3d<Vector<Real, STATE_NF>, GNX, GNX, GNX> F[3];
    Array3d<Vector<Real, STATE_NF>, GNX, GNX, GNX> E;
    bool amr_has[3];
    int amr_cnt[3];
    MPI_Request send_request[8];
    MPI_Request recv_request[8];
    MPI_Request amr_send_request[3][8];
    MPI_Request amr_recv_request[6][8];
    State* mpi_buffer[8];
    State* mpi_amr_buffer[6][8];
    Real dx;
    Vector<int, 3> offset;
    static void mpi_datatypes_initialize();
    virtual void initialize() = 0;
    virtual Hydro* new_octnode() const = 0;
    virtual void physical_boundary(int);
    void inject_from_children_send(int);
    void inject_from_children_recv(int);
    void inject_from_children_send_wait(int);
    void inject_from_children_recv_wait(int);
    void max_dt_bnd_comm(int);
    void max_dt_bnd_recv_wait(int);
    void max_dt_bnd_send_wait(int);
    void flux_bnd_comm(int);
    void flux_compute(int);
    virtual void flux_physical_bounds(int);
    void flux_bnd_recv_wait(int);
    void flux_bnd_send_wait(int);
    void flux_cf_adjust_send(int);
    void flux_cf_adjust_recv(int);
    void flux_cf_adjust_send_wait(int);
    void flux_cf_adjust_recv_wait(int);
    void amr_bnd_send(int);
    void enforce_dual_energy_formalism(int);
    void amr_bnd_send_wait(int);
    void redistribute_send();
    void redistribute_recv();
    void sync(int);
    void error_from_parent_send(int);
    void error_from_parent_recv(int);
    void error_from_parent_send_wait(int);
    void error_from_parent_recv_wait(int);
protected:
    Array3d<State, GNX, GNX, GNX> U;
    Array3d<State, GNX, GNX, GNX> U0;
    Array3d<Vector<Real, STATE_NF>, GNX, GNX, GNX> D;
    virtual void compute_update(int);
    static Reconstruct reconstruct;
    static Real eax;
    Array3d<State, GNX, GNX, GNX> E0;
    static State DFO;
    static State get_flow_off() {
        return FO;
    }
    State get_dfodt() const {
        return DFO;
    }
    static void set_origin(const _3Vec& o) {
        origin = o;
    }
    static _3Vec get_origin() {
        return origin;
    }
    static Real last_dt;
    virtual void write(MPI_File* fh);
    virtual void read(MPI_File* fh);
    State get_flux(int dir, int i, int j, int k) const {
        return F[dir](i, j, k);
    }
    _3Vec Xfx(int, int, int) const;
    _3Vec Xfy(int, int, int) const;
    _3Vec Xfz(int, int, int) const;
    virtual void max_dt_compute(int);
    void add_to_dudt(int, int, int, const Vector<Real, STATE_NF>&);
    Vector<Real, STATE_NF> get_dudt(int, int, int) const;
    static void set_dt(Real);
    static void set_beta(Real);
    static void substep_driver();
    static bool check_for_refine();
    static void inject_from_children();
    static void store();
    static void sum_outflows();
    virtual void compute_dudt(int);
    virtual void compute_flow_off();
    virtual void create_child(const ChildIndex&);
    virtual void inject_from_parent(ChildIndex);
    virtual Real get_output_point(int i, int j, int k, int l) const;
    virtual const char* output_field_names(int) const;
    virtual void set_refine_flags();
    virtual int nvar_output() const;
    virtual void allocate_arrays();
    virtual void deallocate_arrays();
    static void redistribute_grids();
    static void shadow_off();
    static Real max_dt_driver();
    static void boundary_driver();
    bool zone_is_refined(int, int, int) const;
    static void setup_grid_structure(bool = false);
    static Real next_dt(bool* do_output, bool* last_step, int*, Real freq = OUTPUT_TIME_FREQ);
    static void step(Real dt);
    void reduce_dt(Real dt);
public:
    void mult_dx(Real d) {
        if (get_level() == 0) {
            h0 *= d;
            origin *= d;
        }
        dx *= d;
        for (int i = 0; i < 8; i++) {
            if (get_child(i)) {
                dynamic_cast<Hydro*>(get_child(i))->Hydro::mult_dx(d);
            }
        }
    }
    static Real h0;
    _3Vec X(const Vector<int, 3> v) const {
        return X(v[0], v[1], v[2]);
    }
    _3Vec X(int, int, int) const;
    bool using_shadow() const {
        return shadow;
    }
    State get_shadow_error(int i, int j, int k) const {
        if (!shadow) {
            abort();
        }
        return E0(i, j, k);
    }
    Real get_dx() const;
    virtual Real xc(int) const;
    virtual Real yc(int) const;
    virtual Real zc(int) const;
    virtual Real xf(int) const;
    virtual Real yf(int) const;
    virtual Real zf(int) const;
    Vector<State, 4> state_sum() const;
    State& operator()(int, int, int);
    State& operator()(Vector<int, 3>& v) {
        return (*this)(v[0], v[1], v[2]);
    }
    State operator()(int, int, int) const;
    static void run(int, char*[]);
    virtual void init();
    Hydro();
    virtual ~Hydro();
};

#endif

#endif /* GRID_NODE_H_ */
