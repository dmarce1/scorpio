#include <mpi.h>
#include "hydro_grid.h"
#include "../virtual_process.h"
#include "../defs.h"

#ifdef USE_HYDRO_GRID

MPI_Datatype HydroGrid::MPI_interior_t;
MPI_Datatype HydroGrid::MPI_state_t;
MPI_Datatype HydroGrid::MPI_guard_send_t[6];
MPI_Datatype HydroGrid::MPI_guard_recv_t[6];
MPI_Datatype HydroGrid::MPI_child_t[8];
Real HydroGrid::h0 = (2.0 * GRID_DIM / Real(GNX - 2 * BW));
Real HydroGrid::_dt, HydroGrid::_beta, HydroGrid::eax;
bool HydroGrid::initialized = false;
bool HydroGrid::shadow = true;
Real HydroGrid::last_dt = -1.0;
_3Vec HydroGrid::origin;
State HydroGrid::FO0 = Vector<Real, STATE_NF>(0.0);
State HydroGrid::FO = Vector<Real, STATE_NF>(0.0);
State HydroGrid::DFO = Vector<Real, STATE_NF>(0.0);

HydroGrid::ifunc_t HydroGrid::es[GRID_ES_SIZE] = { &HydroGrid::inject_from_children_recv, &HydroGrid::inject_from_children_recv_wait,
        &HydroGrid::inject_from_children_send, &HydroGrid::inject_from_children_send_wait, &HydroGrid::sync, &HydroGrid::flux_bnd_comm,
        &HydroGrid::physical_boundary, &HydroGrid::flux_bnd_recv_wait, &HydroGrid::flux_bnd_send_wait, &HydroGrid::amr_bnd_send, &HydroGrid::amr_bnd_send_wait,
        &HydroGrid::sync, &HydroGrid::enforce_dual_energy_formalism, &HydroGrid::max_dt_compute };

HydroGrid::ifunc_t HydroGrid::cs[GRID_CS_SIZE] = { &HydroGrid::flux_bnd_comm, &HydroGrid::physical_boundary, &HydroGrid::flux_bnd_recv_wait,
        &HydroGrid::flux_bnd_send_wait, &HydroGrid::amr_bnd_send, &HydroGrid::amr_bnd_send_wait, &HydroGrid::flux_compute, &HydroGrid::flux_cf_adjust_recv,
        &HydroGrid::flux_cf_adjust_recv_wait, &HydroGrid::flux_cf_adjust_send, &HydroGrid::flux_cf_adjust_send_wait, &HydroGrid::sync, &HydroGrid::compute_dudt,
        &HydroGrid::error_from_parent_recv, &HydroGrid::error_from_parent_recv_wait, &HydroGrid::error_from_parent_send,
        &HydroGrid::error_from_parent_send_wait, &HydroGrid::compute_update, &HydroGrid::inject_from_children_recv, &HydroGrid::inject_from_children_recv_wait,
        &HydroGrid::inject_from_children_send, &HydroGrid::inject_from_children_send_wait };

HydroGrid::ifunc_t HydroGrid::cs_children[4] = { &HydroGrid::inject_from_children_recv, &HydroGrid::inject_from_children_recv_wait,
        &HydroGrid::inject_from_children_send, &HydroGrid::inject_from_children_send_wait };

bool HydroGrid::check_for_refine() {
    inject_from_children();
    return OctNode::check_for_refine();
}

void HydroGrid::redistribute_grids() {
    HydroGrid** send_list;
    int count, send_count;
    HydroGrid* g;
    get_root()->compute_distribution();
    count = OctNode::get_local_node_cnt();
    send_list = new HydroGrid*[count];
    send_count = 0;
//	printf("B %i MB\n", heap_bytes_used() / 1024 / 1024);
    for (int i = 0; i < count; i++) {
        g = dynamic_cast<HydroGrid*>(OctNode::get_local_node(i));
        g->redistribute_send();
        send_list[i] = g;

    }
    get_root()->allocate_nodes();
    get_root()->find_local_nodes();
    count = OctNode::get_local_node_cnt();
    for (int i = 0; i < count; i++) {
        g = dynamic_cast<HydroGrid*>(OctNode::get_local_node(i));
        g->redistribute_recv();

    }
    count = send_count;
//	printf("D %i MB\n", heap_bytes_used() / 1024 / 1024);
    for (int i = 0; i < count; i++) {
        g = send_list[i];
        if (*(g->send_request) != MPI_REQUEST_NULL ) {
            MPI_Wait(g->send_request, MPI_STATUS_IGNORE );
            delete[] g->mpi_buffer[0];
        }
    }
    count = OctNode::get_local_node_cnt();
    for (int i = 0; i < count; i++) {
        g = dynamic_cast<HydroGrid*>(OctNode::get_local_node(i));
        if (*(g->recv_request) != MPI_REQUEST_NULL ) {
            MPI_Wait(g->recv_request, MPI_STATUS_IGNORE );
        }

    }
    delete[] send_list;
}

void HydroGrid::store() {
    HydroGrid* g;
    FO0 = FO;
    for (int i = 0; i < OctNode::get_local_node_cnt(); i++) {
        g = dynamic_cast<HydroGrid*>(OctNode::get_local_node(i));
        for (int k = BW; k < GNX - BW; k++) {
            for (int j = BW; j < GNX - BW; j++) {
                for (int i = BW; i < GNX - BW; i++) {
                    if (shadow) {
                        g->E0(i, j, k) = Vector<Real, STATE_NF>(0.0);
                    }
                    //          g->U(i, j, k).floor(g->X(i, j, k));
                    g->U0(i, j, k) = g->U(i, j, k);
                }
            }
        }
    }
}

void HydroGrid::shadow_off() {
    shadow = false;
}

Real HydroGrid::next_dt(bool* do_output, bool* last_step, int* ostep_cnt, Real freq) {
    Real dt, tleft, next_output_time;

    *last_step = false;
    dt = max_dt_driver();
//	printf( "%e %e\n", dt, last_dt);
    if (last_dt < 0.0) {
        dt = MAXINITDT * dt;
    } else {
        dt = min(last_dt * MAXDTINC, dt);
    }
    last_dt = dt;
    *do_output = false;
    tleft = TIME_MAX - HydroGrid::get_time();
    next_output_time = Real(*ostep_cnt + 1) * freq;
    if (dt + HydroGrid::get_time() >= next_output_time) {
        dt = next_output_time - HydroGrid::get_time();
        (*ostep_cnt)++;
        *do_output = true;
    } else {
        //	printf( "dt = %e\n", dt);
        //	printf( "next_output_time = %e\n", next_output_time);
        //	printf( "HydroGrid::get_time() = %e\n", HydroGrid::get_time());
        dt = min(dt, (next_output_time - HydroGrid::get_time()) / Real(int((next_output_time - HydroGrid::get_time()) / dt + 1.0)));
        dt *= (1.0 + 1.0e-9);
    }
    if (tleft <= dt) {
        dt = tleft;
        *last_step = true;
        if (next_output_time == TIME_MAX) {
            *do_output = true;
        }
    }
    return dt;
}

void HydroGrid::set_dt(Real a) {
    _dt = a;
}
void HydroGrid::set_beta(Real beta) {
    _beta = beta;
}

void HydroGrid::run(int argc, char* argv[]) {
    Real start_time, dt, tm;
    FILE* fp;
    bool do_output, last_step;
    int step_cnt = 0;
    int ostep_cnt = 0;
    int nnodes;
    Real avg_node_cnt;
    int max_node_cnt, min_node_cnt;

    setup_grid_structure();

    if (OUTPUT_TIME_FREQ <= TIME_MAX) {
        get_root()->output("X", 0.0, GNX, BW);
    }

    avg_node_cnt = 0.0;
    max_node_cnt = min_node_cnt = OctNode::get_node_cnt();

    start_time = MPI_Wtime();
    do {
        dt = next_dt(&do_output, &last_step, &ostep_cnt);
        step(dt);
        if (step_cnt % (GNX - 2 * BW)== 0){
            //	check_for_refine();
            if (check_for_refine()) {
                redistribute_grids();
            }
        }
        step_cnt++;
        if (MPI_rank() == 0) {
            nnodes = OctNode::get_node_cnt();
            avg_node_cnt = (avg_node_cnt * Real(step_cnt - 1) + Real(nnodes)) / Real(step_cnt);
            max_node_cnt = max(max_node_cnt, nnodes);
            min_node_cnt = min(min_node_cnt, nnodes);
            printf("step=%i t=%e dt=%e lmax=%i ngrids=%i avg ngrids=%i", step_cnt, HydroGrid::get_time(), dt, OctNode::get_max_level(), nnodes,
                    (int) avg_node_cnt);
        }
        if (do_output) {
            if (MPI_rank() == 0) {
                printf("*");
            }
            get_root()->output("X", nint(HydroGrid::get_time() / OUTPUT_TIME_FREQ), GNX, BW);
        }
        if (MPI_rank() == 0) {
            printf("\n");
        }
    } while (!last_step);

    if (MPI_rank() == 0) {
        char str[81];
        sprintf(str, "scaling.%i.txt", OctNode::get_max_level_allowed());
        tm = (MPI_Wtime() - start_time);
        fp = fopen(str, "at");
        fprintf(fp, "%i %e %i %i %i\n", MPI_size(), tm, (int) avg_node_cnt, min_node_cnt, max_node_cnt);
        fclose(fp);
    }
}

void HydroGrid::setup_grid_structure(bool one_iter) {
    Real dt;
    HydroGrid::set_time(0.0);
    OctNode::get_root()->find_local_nodes();
    OctNode::initialize_grids();
    dt = max_dt_driver();
    dt *= min(1.0, MAXINITDT);
    for (int l = 0; l < OctNode::get_max_level_allowed()+2; l++) {
        step(dt);
      //  printf( "%e\n", dt );
        HydroGrid::set_time(0.0);
        if (check_for_refine()) {
            redistribute_grids();
        }

        if (MPI_rank() == 0) {
            printf("%i %i\n", OctNode::get_max_level(), OctNode::get_node_cnt());
        }
        OctNode::initialize_grids();
        dt = min(max_dt_driver(), MAXINITDT);
        if (one_iter) {
            break;
        }
    }

}

void HydroGrid::step(Real dt) {
    store();
    HydroGrid::set_dt(dt);
    HydroGrid::set_beta(1.0);
    substep_driver();
    HydroGrid::set_beta(0.25);
    substep_driver();
    HydroGrid::set_beta(2.0 / 3.0);
    substep_driver();
    HydroGrid::set_time(HydroGrid::get_time() + dt);
//    HydroGrid::boundary_driver();
}

void HydroGrid::substep_driver() {
    DFO = Vector<Real, STATE_NF>(0.0);
    HydroGrid** list = new HydroGrid*[get_local_node_cnt()];
    for (int i = 0; i < get_local_node_cnt(); i++) {
        list[i] = dynamic_cast<HydroGrid*>(get_local_node(i));
    }
    run_program(list, get_local_node_cnt(), cs, GRID_CS_SIZE, 3);
    delete[] list;
    Real tmp[STATE_NF];
    for (int i = 0; i < STATE_NF; i++) {
        tmp[i] = DFO[i];
    }
    Real tmp2[STATE_NF];
    MPI_Allreduce(tmp, tmp2, STATE_NF, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );
    for (int i = 0; i < STATE_NF; i++) {
        DFO[i] = tmp2[i];
    }
    FO = (FO + DFO * _dt) * _beta + FO0 * (1.0 - _beta);

}

void HydroGrid::inject_from_children() {
    DFO = Vector<Real, STATE_NF>(0.0);
    HydroGrid** list = new HydroGrid*[get_local_node_cnt()];
    for (int i = 0; i < get_local_node_cnt(); i++) {
        list[i] = dynamic_cast<HydroGrid*>(get_local_node(i));
    }
    run_program(list, get_local_node_cnt(), cs_children, 4, 1);
    delete[] list;

}

Real HydroGrid::max_dt_driver() {
    Real dt_all, dt;
    eax = 1.0e+99;
    HydroGrid** list = new HydroGrid*[get_local_node_cnt()];
    for (int i = 0; i < get_local_node_cnt(); i++) {
        list[i] = dynamic_cast<HydroGrid*>(get_local_node(i));
    }
    run_program(list, get_local_node_cnt(), es, GRID_ES_SIZE, 3);
    delete[] list;
    dt = eax;
    MPI_Allreduce(&dt, &dt_all, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD );
    dt = dt_all;
    return dt;
}

void HydroGrid::boundary_driver() {
    Real dt_all, dt;
    HydroGrid** list = new HydroGrid*[get_local_node_cnt()];
    for (int i = 0; i < get_local_node_cnt(); i++) {
        list[i] = dynamic_cast<HydroGrid*>(get_local_node(i));
    }
    run_program(list, get_local_node_cnt(), es, GRID_ES_SIZE - 1, 3);
    delete[] list;
}

void HydroGrid::mpi_datatypes_initialize() {
    static bool initialized = false;
    if (!initialized) {
        const int lbi0 = BW;
        const int ubi0 = 2 * BW - 1;
        const int lbi1 = GNX - 2 * BW;
        const int ubi1 = GNX - BW - 1;
        const int lbe0 = 0;
        const int ube0 = BW - 1;
        const int lbe1 = GNX - BW;
        const int ube1 = GNX - 1;
        const int hbnd = GNX - BW - 1;
        const int lbnd = BW;
        Vector<int, 3> lb, ub, o;

        MPI_Type_contiguous(sizeof(State), MPI_BYTE, &MPI_state_t);
        MPI_Type_commit(&MPI_state_t);

        for (ChildIndex ci = 0; ci < 8; ci++) {
            lb = BW;
            ub = GNX / 2 - 1;
            o = ci.vec();
            o *= GNX / 2 - BW;
            lb += o;
            ub += o;
            Array3d<Real, GNX, GNX, GNX>::mpi_datatype(&(MPI_child_t[ci]), lb[0], ub[0], lb[1], ub[1], lb[2], ub[2], MPI_state_t);
        }
        Array3d<Real, GNX, GNX, GNX>::mpi_datatype(&MPI_interior_t, BW, GNX - BW - 1, BW, GNX - BW - 1, BW, GNX - BW - 1, MPI_state_t);
        Array3d<Real, GNX, GNX, GNX>::mpi_datatype(MPI_guard_send_t + XL, lbi0, ubi0, lbnd, hbnd, lbnd, hbnd, MPI_state_t);
        Array3d<Real, GNX, GNX, GNX>::mpi_datatype(MPI_guard_send_t + XU, lbi1, ubi1, lbnd, hbnd, lbnd, hbnd, MPI_state_t);
        Array3d<Real, GNX, GNX, GNX>::mpi_datatype(MPI_guard_recv_t + XL, lbe0, ube0, lbnd, hbnd, lbnd, hbnd, MPI_state_t);
        Array3d<Real, GNX, GNX, GNX>::mpi_datatype(MPI_guard_recv_t + XU, lbe1, ube1, lbnd, hbnd, lbnd, hbnd, MPI_state_t);
        Array3d<Real, GNX, GNX, GNX>::mpi_datatype(MPI_guard_send_t + YL, lbnd, hbnd, lbi0, ubi0, lbnd, hbnd, MPI_state_t);
        Array3d<Real, GNX, GNX, GNX>::mpi_datatype(MPI_guard_send_t + YU, lbnd, hbnd, lbi1, ubi1, lbnd, hbnd, MPI_state_t);
        Array3d<Real, GNX, GNX, GNX>::mpi_datatype(MPI_guard_recv_t + YL, lbnd, hbnd, lbe0, ube0, lbnd, hbnd, MPI_state_t);
        Array3d<Real, GNX, GNX, GNX>::mpi_datatype(MPI_guard_recv_t + YU, lbnd, hbnd, lbe1, ube1, lbnd, hbnd, MPI_state_t);
        Array3d<Real, GNX, GNX, GNX>::mpi_datatype(MPI_guard_send_t + ZL, lbnd, hbnd, lbnd, hbnd, lbi0, ubi0, MPI_state_t);
        Array3d<Real, GNX, GNX, GNX>::mpi_datatype(MPI_guard_send_t + ZU, lbnd, hbnd, lbnd, hbnd, lbi1, ubi1, MPI_state_t);
        Array3d<Real, GNX, GNX, GNX>::mpi_datatype(MPI_guard_recv_t + ZL, lbnd, hbnd, lbnd, hbnd, lbe0, ube0, MPI_state_t);
        Array3d<Real, GNX, GNX, GNX>::mpi_datatype(MPI_guard_recv_t + ZU, lbnd, hbnd, lbnd, hbnd, lbe1, ube1, MPI_state_t);
        initialized = true;
    }
}
#endif
