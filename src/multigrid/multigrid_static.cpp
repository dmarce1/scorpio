#include "../array3d.h"
#include "multigrid.h"
#include <mpi.h>
#include <stdlib.h>

MPI_Datatype MultiGrid::MPI_send_bnd1_t[6];
MPI_Datatype MultiGrid::MPI_recv_bnd1_t[6];
MPI_Datatype MultiGrid::MPI_send_bnd2_t[6];
MPI_Datatype MultiGrid::MultiGrid::MPI_recv_bnd2_t[6];
MPI_Datatype MultiGrid::MPI_send_bnd3_t[6];
MPI_Datatype MultiGrid::MPI_recv_bnd3_t[6];
MPI_Datatype MultiGrid::MPI_comm_child_t[8];
MPI_Datatype MultiGrid::MPI_recv_amr_child_t[8];
Real MultiGrid::h0 = (2.0 * GRID_DIM / Real(PNX - 2));
_3Vec MultiGrid::origin;

MultiGrid::ifunc_t MultiGrid::es[MG_ES_STAGE_MAX + 1] = { &MultiGrid::phi_real_boundary_begin_loop, &MultiGrid::phi_real_boundary_communicate,
        &MultiGrid::phi_real_boundary_wait, &MultiGrid::phi_real_boundary_end_loop, &MultiGrid::phi_amr_boundary_communicate, &MultiGrid::phi_amr_boundary_wait,
        &MultiGrid::null };

MultiGrid::ifunc_t MultiGrid::cs[COMPUTE_STAGE_MAX + 1] = { &MultiGrid::vdown_init_recv, &MultiGrid::vdown_init_recv_wait, &MultiGrid::vdown_init_compute,
        &MultiGrid::vdown_init_adjust_recv, &MultiGrid::vdown_init_adjust_recv_wait, &MultiGrid::vdown_init_store_dphi, &MultiGrid::relax_begin_down_loop,
        &MultiGrid::relax_communicate, &MultiGrid::relax_wait, &MultiGrid::relax_compute, &MultiGrid::relax_end_down_loop, &MultiGrid::vdown_init_send,
        &MultiGrid::vdown_init_send_wait, &MultiGrid::vdown_init_adjust_send, &MultiGrid::vdown_init_adjust_send_wait, &MultiGrid::vup_init_recv,
        &MultiGrid::vup_init_recv_wait, &MultiGrid::relax_begin_up_loop, &MultiGrid::relax_communicate, &MultiGrid::relax_wait, &MultiGrid::relax_compute,
        &MultiGrid::relax_end_up_loop, &MultiGrid::vup_init_send, &MultiGrid::vup_init_send_wait, &MultiGrid::dphi_amr_boundary_communicate,
        &MultiGrid::dphi_amr_boundary_wait, &MultiGrid::retire_dphi, &MultiGrid::phi_children_recv, &MultiGrid::phi_children_recv_wait,
        &MultiGrid::phi_children_send, &MultiGrid::phi_children_send_wait, &MultiGrid::phi_real_boundary_begin_loop, &MultiGrid::phi_real_boundary_communicate,
        &MultiGrid::phi_real_boundary_wait, &MultiGrid::phi_real_boundary_end_loop, &MultiGrid::phi_amr_boundary_communicate, &MultiGrid::phi_amr_boundary_wait,
        &MultiGrid::forces_compute, &MultiGrid::forces_adjust_send, &MultiGrid::forces_adjust_send_wait, &MultiGrid::forces_adjust_recv,
        &MultiGrid::forces_adjust_recv_wait, &MultiGrid::residual_error_compute, &MultiGrid::null };
/*
 MultiGrid::ifunc_t MultiGrid::es[MG_ES_STAGE_MAX + 1] = { &MultiGrid::phi_face_compute, &MultiGrid::phi_face_adjust_send, &MultiGrid::phi_face_adjust_send_wait,
 &MultiGrid::phi_face_adjust_recv, &MultiGrid::phi_face_adjust_recv_wait, &MultiGrid::null };
 */
static void MPI_mask_to_type(MPI_Datatype*, const Array3d<bool, PNX, PNX, PNX>&);

Real MultiGrid::get_phi_at(Real x, Real y, Real z) {
    Real p, tmp;
    MultiGrid* g;
    p = 0.0;
    for (int l = 0; l < get_local_node_cnt(); l++) {
        g = dynamic_cast<MultiGrid*>(get_local_node(l));
        if (MPI_rank() == g->proc()) {
            for (int k = 1; k < PNX - 1; k++) {
                if (z >= g->MultiGrid::zf(k) && z < g->MultiGrid::zf(k + 1)) {
                    for (int j = 1; j < PNX - 1; j++) {
                        if (y >= g->MultiGrid::yf(j) && y < g->MultiGrid::yf(j + 1)) {
                            for (int i = 1; i < PNX - 1; i++) {
                        //        printf("%e %e %e\n", g->MultiGrid::xf(i), x, g->MultiGrid::xf(i + 1));
                               if (x >= g->MultiGrid::xf(i) && x < g->MultiGrid::xf(i + 1)) {
                                    if (!g->poisson_zone_is_refined(i, j, k)) {
                                        p = g->phi(i, j, k);
                                       // printf("%e\n", p);
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

void MultiGrid::redistribute_grids() {
#ifdef HYDRO_GRAV_GRID
    printf("MultiGrid::redistribute_grids should not be called\n");
    abort();
#endif
    get_root()->compute_distribution();
    get_root()->allocate_nodes();
    get_root()->find_local_nodes();
}

Real MultiGrid::vcycle() {
    MPI_datatypes_init();
    Real part, sum;
    MultiGrid** list = new MultiGrid*[get_local_node_cnt()];
    for (int i = 0; i < get_local_node_cnt(); i++) {
        list[i] = dynamic_cast<MultiGrid*>(get_local_node(i));
    }
    run_program(list, get_local_node_cnt(), cs, COMPUTE_STAGE_MAX);
    delete[] list;
    sum = 0.0;
    for (int i = 0; i < OctNode::get_local_node_cnt(); i++) {
        sum += dynamic_cast<MultiGrid*>(OctNode::get_local_node(i))->st0;
    }
    part = sum;
    MPI_Allreduce(&part, &sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );
    return sum;
}

bool MultiGrid::check_for_refine() {
    bool rc;
    rc = OctNode::check_for_refine();
    if (rc) {
        redistribute_grids();
        int count = OctNode::get_local_node_cnt();
        for (int i = 0; i < count; i++) {
            dynamic_cast<MultiGrid*>(OctNode::get_local_node(i))->phi_calc_amr_bounds();
        }
    }
    return rc;
}

void MultiGrid::MPI_datatypes_init() {
    static bool initialized = false;
    if (!initialized) {
        const int lbi = 1;
        const int ubi = PNX - 2;
        const int lbe = 0;
        const int ube = PNX - 1;
        Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_send_bnd1_t + XL, lbi, lbi, lbi, ubi, lbi, ubi, MPI_DOUBLE_PRECISION );
        Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_send_bnd1_t + YL, lbi, ubi, lbi, lbi, lbi, ubi, MPI_DOUBLE_PRECISION );
        Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_send_bnd1_t + ZL, lbi, ubi, lbi, ubi, lbi, lbi, MPI_DOUBLE_PRECISION );
        Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_send_bnd1_t + XU, ubi, ubi, lbi, ubi, lbi, ubi, MPI_DOUBLE_PRECISION );
        Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_send_bnd1_t + YU, lbi, ubi, ubi, ubi, lbi, ubi, MPI_DOUBLE_PRECISION );
        Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_send_bnd1_t + ZU, lbi, ubi, lbi, ubi, ubi, ubi, MPI_DOUBLE_PRECISION );
        Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_recv_bnd1_t + XL, lbe, lbe, lbi, ubi, lbi, ubi, MPI_DOUBLE_PRECISION );
        Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_recv_bnd1_t + YL, lbi, ubi, lbe, lbe, lbi, ubi, MPI_DOUBLE_PRECISION );
        Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_recv_bnd1_t + ZL, lbi, ubi, lbi, ubi, lbe, lbe, MPI_DOUBLE_PRECISION );
        Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_recv_bnd1_t + XU, ube, ube, lbi, ubi, lbi, ubi, MPI_DOUBLE_PRECISION );
        Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_recv_bnd1_t + YU, lbi, ubi, ube, ube, lbi, ubi, MPI_DOUBLE_PRECISION );
        Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_recv_bnd1_t + ZU, lbi, ubi, lbi, ubi, ube, ube, MPI_DOUBLE_PRECISION );
        Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_send_bnd2_t + XL, lbi, lbi, lbi, ubi, lbi, ubi, MPI_DOUBLE_PRECISION );
        Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_send_bnd2_t + YL, lbe, ube, lbi, lbi, lbi, ubi, MPI_DOUBLE_PRECISION );
        Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_send_bnd2_t + ZL, lbe, ube, lbe, ube, lbi, lbi, MPI_DOUBLE_PRECISION );
        Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_send_bnd2_t + XU, ubi, ubi, lbi, ubi, lbi, ubi, MPI_DOUBLE_PRECISION );
        Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_send_bnd2_t + YU, lbe, ube, ubi, ubi, lbi, ubi, MPI_DOUBLE_PRECISION );
        Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_send_bnd2_t + ZU, lbe, ube, lbe, ube, ubi, ubi, MPI_DOUBLE_PRECISION );
        Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_recv_bnd2_t + XL, lbe, lbe, lbi, ubi, lbi, ubi, MPI_DOUBLE_PRECISION );
        Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_recv_bnd2_t + YL, lbe, ube, lbe, lbe, lbi, ubi, MPI_DOUBLE_PRECISION );
        Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_recv_bnd2_t + ZL, lbe, ube, lbe, ube, lbe, lbe, MPI_DOUBLE_PRECISION );
        Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_recv_bnd2_t + XU, ube, ube, lbi, ubi, lbi, ubi, MPI_DOUBLE_PRECISION );
        Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_recv_bnd2_t + YU, lbe, ube, ube, ube, lbi, ubi, MPI_DOUBLE_PRECISION );
        Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_recv_bnd2_t + ZU, lbe, ube, lbe, ube, ube, ube, MPI_DOUBLE_PRECISION );
        Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_send_bnd3_t + XL, lbi, lbi, lbe, ube, lbe, ube, MPI_DOUBLE_PRECISION );
        Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_send_bnd3_t + YL, lbe, ube, lbi, lbi, lbe, ube, MPI_DOUBLE_PRECISION );
        Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_send_bnd3_t + ZL, lbe, ube, lbi, ubi, lbi, lbi, MPI_DOUBLE_PRECISION );
        Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_send_bnd3_t + XU, ubi, ubi, lbe, ube, lbe, ube, MPI_DOUBLE_PRECISION );
        Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_send_bnd3_t + YU, lbe, ube, ubi, ubi, lbe, ube, MPI_DOUBLE_PRECISION );
        Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_send_bnd3_t + ZU, lbe, ube, lbe, ube, ubi, ubi, MPI_DOUBLE_PRECISION );
        Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_recv_bnd3_t + XL, lbe, lbe, lbe, ube, lbe, ube, MPI_DOUBLE_PRECISION );
        Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_recv_bnd3_t + YL, lbe, ube, lbe, lbe, lbe, ube, MPI_DOUBLE_PRECISION );
        Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_recv_bnd3_t + ZL, lbe, ube, lbe, ube, lbe, lbe, MPI_DOUBLE_PRECISION );
        Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_recv_bnd3_t + XU, ube, ube, lbe, ube, lbe, ube, MPI_DOUBLE_PRECISION );
        Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_recv_bnd3_t + YU, lbe, ube, ube, ube, lbe, ube, MPI_DOUBLE_PRECISION );
        Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_recv_bnd3_t + ZU, lbe, ube, lbe, ube, ube, ube, MPI_DOUBLE_PRECISION );

        const int ds = PNX / 2 - 1;
        int ci, xlb, xub, ylb, yub, zlb, zub;

        Array3d<bool, PNX, PNX, PNX> use(true);
        for (int k = 0; k < 2; k++) {
            zlb = lbi + k * (PNX / 2 - 1);
            zub = zlb + ds - 1;
            for (int j = 0; j < 2; j++) {
                ylb = lbi + j * (PNX / 2 - 1);
                yub = ylb + ds - 1;
                for (int i = 0; i < 2; i++) {
                    xlb = lbi + i * (PNX / 2 - 1);
                    xub = xlb + ds - 1;
                    ci = 4 * k + 2 * j + i;
                    for (int k0 = 0; k0 <= PNX - 1; k0++) {
                        for (int j0 = 0; j0 <= PNX - 1; j0++) {
                            for (int i0 = 0; i0 <= PNX - 1; i0++) {
                                use(i0, j0, k0) = false;
                            }
                        }
                    }
                    for (int k0 = zlb; k0 <= zub; k0++) {
                        for (int j0 = ylb; j0 <= yub; j0++) {
                            for (int i0 = xlb; i0 <= xub; i0++) {
                                use(i0, j0, k0) = true;
                            }
                        }
                    }
                    MPI_mask_to_type(MPI_recv_amr_child_t + ci, use);
                    Array3d<Real, PNX, PNX, PNX>::mpi_datatype(MPI_comm_child_t + ci, xlb, xub, ylb, yub, zlb, zub, MPI_DOUBLE_PRECISION );
                }
            }
        }
        initialized = true;
    }
}

static void MPI_mask_to_type(MPI_Datatype* newtype, const Array3d<bool, PNX, PNX, PNX>& use) {
    int indices[PNX * PNX * PNX / 8];
    int blocklens[PNX * PNX * PNX / 8];
    int count, index;
    count = 0;
    blocklens[0] = 0;
    indices[0] = 0;
    int sz = 0;
    for (int k0 = 0; k0 < PNX; k0++) {
        for (int j0 = 0; j0 < PNX; j0++) {
            for (int i0 = 0; i0 < PNX; i0++) {
                if (use(i0, j0, k0) == true) {
                    index = i0 + PNX * (j0 + PNX * k0);
                    if (blocklens[count] + indices[count] == index) {
                        blocklens[count]++;
                        sz++;
                    } else {
                        sz++;
                        count++;
                        blocklens[count] = 1;
                        indices[count] = index;
                    }
                }
            }
        }
    }
    count++;
    MPI_Type_indexed(count, blocklens, indices, MPI_DOUBLE_PRECISION, newtype);
    MPI_Type_commit(newtype);

}
