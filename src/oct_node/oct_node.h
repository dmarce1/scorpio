#ifndef OCOctNode_NODE_H_
#define OCOctNode_NODE_H_
#include <mpi.h>
#include "../assert.h"
#include "../bits.h"
#include "../comm.h"
#include "child_index.h"
#include "grid_output.h"
#include "oct_face.h"
#include "oct_connections.h"

#define MAX_LEVEL 64

class OctNode: public OctConnections {

private:
    static bool initialized;
    static Real grid_time;
    static int maxlev;
    static int id;
    static OctNode** local_node;
    static int local_node_cnt;
    static int node_sums[MAX_LEVEL];
    static int node_counter;
    static int node_total;
    static Real* output_buffer;
    static int max_refine_level;
    int myid;
    int last_processor;
    int processor;
    int next_processor;
    int level;
    Vector<int, 3> location;
    static OctNode* root;
    bool refine[OCT_NCHILD];
    static void reduce_refine_flags();
    int get_leaf_cnt() const;
    void set_proc(int p);
    void set_to_next_proc();
    void use_refine_bits(Bits* bits, int* n);
    void set_refine_bits(Bits* bits, int* n);
    bool use_refine_flags();
    void clear_refine_flags();
    void propagate_refine_flags_up();
    void enforce_proper_nesting();
    void output(grid_output_t* ptr, int, int, Real dtheta = 0.0) const;
protected:
    Vector<int, 3> get_location() const {
        return location;
    }
    static bool check_for_refine();
    static OctNode** get_local_node_ptr();
    static int get_max_level();
    static int get_node_cnt();
    static int get_node_cnt(int lev);
    static Real get_time();
    static void set_time(Real a);
    virtual void allocate_arrays() = 0;
    virtual void deallocate_arrays() = 0;
    virtual int nvar_output() const = 0;
    virtual Real get_output_point(int i, int j, int k, int l) const = 0;
    virtual const char* output_field_names(int) const = 0;
    virtual Real xf(int) const = 0;
    virtual Real yf(int) const = 0;
    virtual Real zf(int) const = 0;
    virtual void initialize() = 0;
    virtual OctNode* new_octnode() const = 0;
    virtual void set_refine_flags() = 0;
    virtual void destroy_child(const ChildIndex&);
    virtual void create_child(const ChildIndex&);
    int get_level() const;
    bool is_real_bound(OctFace f) const;
    int last_proc() const;
    ChildIndex my_child_index() const;
    int next_proc() const;
    virtual void write(MPI_File* fh) = 0;
    virtual void read(MPI_File* fh) = 0;
public:
    bool get_refine_flag(const ChildIndex&) const;
    void set_refine_flag(const ChildIndex&, bool);
    static void initialize_grids();
    static OctNode* get_local_node(int i);
    static int get_local_node_cnt();
    void init();
    bool is_phys_bound(OctFace) const;
    virtual void write_checkpoint(MPI_File*) ;
    virtual void read_checkpoint(MPI_File*);
    static OctNode* get_root();
    static void set_max_level_allowed(int l);
    static int get_max_level_allowed();
    virtual ~OctNode();
    OctNode();
    void allocate_nodes();
    void compute_distribution();
    void find_local_nodes();
    int get_id() const;
    void output(const char* prefix, int counter, int, int, double dtheta = 0.0) const;
    int proc() const;
    bool is_amr_bound(OctFace f) const;
    virtual void expand_grid();

};

#endif
