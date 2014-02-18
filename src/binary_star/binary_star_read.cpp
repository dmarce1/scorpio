#include "binary_star.h"
#include <silo.h>

#ifdef BINARY_STAR
struct node_t {
    double x;
    double y;
    double z;
    double h;
    double rho;
    double CO;
    double He;
    double phi;
    double et, tau;
    double sr, lz, sz;
};

typedef struct {
    DBucdvar* phi;
    DBucdvar* rho;
    DBucdvar* CO;
    DBucdvar* He;
    DBucdvar* tau;
    DBucdvar* sr;
    DBucdvar* lz;
    DBucdvar* sz;
    DBucdvar* et;
} vars_t;

struct frame_t {
    node_t* node_list;
    int node_count;
    double time;
};

void frame_read(frame_t* frame, const char* filename) {
#ifdef USING_MIC
    abort();
    return;
#else

    DBfile* file;
    DBucdmesh* mesh;
    DBzonelist* zones;
    vars_t vars;
    int nels, n, i, ni1, ni2;
    node_t* ptr;

    file = DBOpen(filename, DB_PDB, DB_READ);
    mesh = DBGetUcdmesh(file, "mesh");
    vars.CO = DBGetUcdvar(file, "CO");
    vars.He = DBGetUcdvar(file, "He");
    vars.rho = DBGetUcdvar(file, "rho");
    vars.phi = DBGetUcdvar(file, "phi");
    vars.tau = DBGetUcdvar(file, "tau");
    vars.et = DBGetUcdvar(file, "etot");
    vars.sr = DBGetUcdvar(file, "sR");
    vars.lz = DBGetUcdvar(file, "lz");
    vars.sz = DBGetUcdvar(file, "sz");
    zones = mesh->zones;

    frame->node_count = zones->nzones;
    frame->time = (double) vars.rho->time;
    frame->node_list = new node_t[frame->node_count];
    for (n = 0; n < frame->node_count; n++) {
        ptr = frame->node_list + n;
        ptr->rho = (*reinterpret_cast<double**>(vars.rho->vals))[n];
        ptr->phi = (*reinterpret_cast<double**>(vars.phi->vals))[n];
        ptr->He = (*reinterpret_cast<double**>(vars.He->vals))[n];
        ptr->CO = (*reinterpret_cast<double**>(vars.CO->vals))[n];
        ptr->tau = (*reinterpret_cast<double**>(vars.tau->vals))[n];
        ptr->et = (*reinterpret_cast<double**>(vars.et->vals))[n];
        ptr->sr = (*reinterpret_cast<double**>(vars.sr->vals))[n];
        ptr->lz = (*reinterpret_cast<double**>(vars.lz->vals))[n];
        ptr->sz = (*reinterpret_cast<double**>(vars.sz->vals))[n];
        ni1 = zones->nodelist[8 * n];
        ni2 = zones->nodelist[8 * n + 1];
        ptr->h = reinterpret_cast<double**>(mesh->coords)[0][ni2] - reinterpret_cast<double**>(mesh->coords)[0][ni1];
        ptr->x = reinterpret_cast<double**>(mesh->coords)[0][ni1] + 0.5 * ptr->h;
        ptr->y = reinterpret_cast<double**>(mesh->coords)[1][ni1] + 0.5 * ptr->h;
        ptr->z = reinterpret_cast<double**>(mesh->coords)[2][ni1] + 0.5 * ptr->h;
    }
    DBFreeUcdvar(vars.He);
    DBFreeUcdvar(vars.CO);
    DBFreeUcdvar(vars.rho);
    DBFreeUcdvar(vars.phi);
    DBFreeUcdvar(vars.sr);
    DBFreeUcdvar(vars.sz);
    DBFreeUcdvar(vars.lz);
    DBFreeUcdvar(vars.tau);
    DBFreeUcdvar(vars.et);
    DBFreeUcdmesh(mesh);
    DBClose(file);
#endif
}

void BinaryStar::add_data_point(double x, double y, double z, double h, const State& s, double phi0) {
    // printf("%i\n", get_level());
    //   printf("%e %e\n", get_dx(), h);
#ifndef USE_FMM
    BinaryStar::set_max_level_allowed(get_level() + 1);
    if (fabs(get_dx() - h) < 1.0e-3 * h) {
        int i = int(((x - HydroGrid::xf(0)) / h));
        int j = int(((y - HydroGrid::yf(0)) / h));
        int k = int(((z - HydroGrid::zf(0)) / h));
        U(i, j, k) = s;
        phi(i + 1 - BW, j + 1 - BW, k + 1 - BW) = phi0;
    } else {
        ChildIndex c;
        c.set_x(int(((x - HydroGrid::xf(0)) / get_dx() / double(GNX / 2))));
        c.set_y(int(((y - HydroGrid::yf(0)) / get_dx() / double(GNX / 2))));
        c.set_z(int(((z - HydroGrid::zf(0)) / get_dx() / double(GNX / 2))));
        if (this->get_child(c) == NULL) {
            this->create_child(c);
        }
        dynamic_cast<BinaryStar*>(get_child(c))->add_data_point(x, y, z, h, s, phi0);
    }
#endif
}

void BinaryStar::read_silo(const char* name) {
    if( MPI_rank() != 0 ) {
        return;
    }
    frame_t frame;
    Real hmax;
    _3Vec O, xmin, xmax;
    PhysicalConstants::set_cgs();
    printf("Reading in SILO data...");
    fflush(stdout);
    frame_read(&frame, name);
    xmin = +1.0e+99;
    xmax = -1.0e+99;
    hmax = -1.0e+99;
    for (int i = 0; i < frame.node_count; i++) {
        node_t* n = frame.node_list + i;
        hmax = max(hmax, n->h);
    }
    for (int i = 0; i < frame.node_count; i++) {
        node_t* n = frame.node_list + i;
        xmin[0] = min(xmin[0], n->x - n->h * 0.5);
        xmax[0] = max(xmax[0], n->x + n->h * 0.5);
        xmin[1] = min(xmin[1], n->y - n->h * 0.5);
        xmax[1] = max(xmax[1], n->y + n->h * 0.5);
        xmin[2] = min(xmin[2], n->z - n->h * 0.5);
        xmax[2] = max(xmax[2], n->z + n->h * 0.5);
    }
    printf("%e %e\n", xmin[0], xmax[0]);
    O = xmin + xmax;
    O *= -0.5;
    printf("%e %e %e \n", O[0], O[1], O[2]);
    //  abort();
    double range = (xmax[0] - xmin[0]) / 2.0;
    dynamic_cast<HydroGrid*>(get_root())->HydroGrid::mult_dx(range);
#ifndef USE_FMM
    dynamic_cast<MultiGrid*>(get_root())->MultiGrid::mult_dx(range);
#endif
    set_origin(O);

    for (int i = 0; i < frame.node_count; i++) {
        node_t* n = frame.node_list + i;
        State st;
        st.set_et(n->et);
        st.set_sx(n->sr);
        st.set_sy(n->lz);
        st.set_sz(n->sz);
        st.set_frac(0, n->CO);
        st.set_frac(1, n->He);
        st.set_rho(n->rho);
        st[State::tau_index] = n->tau;
        dynamic_cast<BinaryStar*>(get_root())->add_data_point(n->x, n->y, n->z, n->h, st, n->phi);
        dynamic_cast<BinaryStar*>(get_root())->set_time(frame.time);
    }
    dynamic_cast<BinaryStar*>(get_root())->find_local_nodes();
    printf("done!\n");
    printf("%i grids read\n", dynamic_cast<BinaryStar*>(get_root())->get_node_cnt());
}
#endif
