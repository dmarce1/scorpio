#include <mpi.h>
#include "hydro.h"
#include "../virtual_process.h"
#include "../defs.h"

MPI_Datatype Hydro::MPI_interior_t;
MPI_Datatype Hydro::MPI_state_t;
MPI_Datatype Hydro::MPI_guard_send_t[6];
MPI_Datatype Hydro::MPI_guard_recv_t[6];
MPI_Datatype Hydro::MPI_child_t[8];
Real Hydro::h0 = (2.0 * GRID_DIM / Real(GNX - 2 * BW));
Real Hydro::_dt, Hydro::_beta, Hydro::eax;
bool Hydro::initialized = false;
bool Hydro::shadow = false;
Real Hydro::last_dt = -1.0;
_3Vec Hydro::origin;
State Hydro::FO0 = Vector<Real, STATE_NF>(0.0);
State Hydro::FO = Vector<Real, STATE_NF>(0.0);
State Hydro::DFO = Vector<Real, STATE_NF>(0.0);

Hydro::ifunc_t Hydro::es[GRID_ES_SIZE] = { &Hydro::inject_from_children_recv, &Hydro::inject_from_children_recv_wait, &Hydro::inject_from_children_send,
        &Hydro::inject_from_children_send_wait, &Hydro::flux_bnd_comm, &Hydro::physical_boundary, &Hydro::flux_bnd_recv_wait, &Hydro::flux_bnd_send_wait,
        &Hydro::amr_bnd_send, &Hydro::amr_bnd_send_wait, &Hydro::enforce_dual_energy_formalism, &Hydro::max_dt_compute };

Hydro::ifunc_t Hydro::cs[GRID_CS_SIZE] = {
//       &Hydro::flux_bnd_comm,
//      &Hydro::physical_boundary,
//      &Hydro::flux_bnd_recv_wait,
//     &Hydro::flux_bnd_send_wait,
//    &Hydro::amr_bnd_send,
//   &Hydro::amr_bnd_send_wait,
        &Hydro::flux_compute,
        //     &Hydro::flux_cf_adjust_recv,
        //    &Hydro::flux_cf_adjust_recv_wait,
        //   &Hydro::flux_cf_adjust_send,
        //  &Hydro::flux_cf_adjust_send_wait,
        &Hydro::compute_dudt, &Hydro::compute_update,
//       &Hydro::inject_from_children_recv,
//      &Hydro::inject_from_children_recv_wait,
//     &Hydro::inject_from_children_send,
//    &Hydro::inject_from_children_send_wait
        };

Hydro::ifunc_t Hydro::cs_children[4] = { &Hydro::inject_from_children_recv, &Hydro::inject_from_children_recv_wait, &Hydro::inject_from_children_send,
        &Hydro::inject_from_children_send_wait };

bool Hydro::check_for_refine() {
    inject_from_children();
    return OctNode::check_for_refine();
}

void Hydro::redistribute_grids() {
    Hydro** send_list;
    int count, send_count;
    Hydro* g;
    get_root()->compute_distribution();
    count = OctNode::get_local_node_cnt();
    send_list = new Hydro*[count];
    send_count = 0;
//	printf("B %i MB\n", heap_bytes_used() / 1024 / 1024);
    for (int i = 0; i < count; i++) {
        g = dynamic_cast<Hydro*>(OctNode::get_local_node(i));
        g->redistribute_send();
        send_list[i] = g;

    }
    get_root()->allocate_nodes();
    get_root()->find_local_nodes();
    count = OctNode::get_local_node_cnt();
    for (int i = 0; i < count; i++) {
        g = dynamic_cast<Hydro*>(OctNode::get_local_node(i));
        g->redistribute_recv();

    }
    count = send_count;
//	printf("D %i MB\n", heap_bytes_used() / 1024 / 1024);
    for (int i = 0; i < count; i++) {
        g = send_list[i];
        if (*(g->send_request) != MPI_REQUEST_NULL) {
            MPI_Wait(g->send_request, MPI_STATUS_IGNORE);
            delete[] g->mpi_buffer[0];
        }
    }
    count = OctNode::get_local_node_cnt();
    for (int i = 0; i < count; i++) {
        g = dynamic_cast<Hydro*>(OctNode::get_local_node(i));
        if (*(g->recv_request) != MPI_REQUEST_NULL) {
            MPI_Wait(g->recv_request, MPI_STATUS_IGNORE);
        }

    }
    delete[] send_list;
}

void Hydro::store() {
    Hydro* g;
    FO0 = FO;
    for (int i = 0; i < OctNode::get_local_node_cnt(); i++) {
        g = dynamic_cast<Hydro*>(OctNode::get_local_node(i));
//#pragma omp parallel for collapse(2)
        for (int k = BW; k < GNX - BW; k++) {
            for (int j = BW; j < GNX - BW; j++) {
                for (int i = BW; i < GNX - BW; i++) {
                    g->U0(i, j, k) = g->U(i, j, k);
                }
            }
        }
    }
}

void Hydro::shadow_off() {
    shadow = false;
}

Real Hydro::next_dt(bool* do_output, bool* last_step, int* ostep_cnt, Real freq) {
    Real dt, tleft, next_output_time;

    *last_step = false;
    dt = max_dt_driver();
    if (last_dt < 0.0) {
        dt = MAXINITDT * dt;
    } else {
        dt = min(last_dt * MAXDTINC, dt);
    }
    last_dt = dt;
    *do_output = false;
    tleft = TIME_MAX - Hydro::get_time();
    next_output_time = Real(*ostep_cnt + 1) * freq;
    if (dt + Hydro::get_time() >= next_output_time) {
        dt = next_output_time - Hydro::get_time();
        (*ostep_cnt)++;
        *do_output = true;
    } else {
        // printf( "dt = %e\n", dt);
        //	printf( "next_output_time = %e\n", next_output_time);
        // 	printf( "HydroGrid::get_time() = %e\n", get_time());
        dt = min(dt, (next_output_time - Hydro::get_time()) / Real((long long int) ((next_output_time - Hydro::get_time()) / dt + 1.0)));
        dt *= (1.0 + 1.0e-9);
        //    printf( "-------dt = %e\n", dt);
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

void Hydro::set_dt(Real a) {
    _dt = a;
}
void Hydro::set_beta(Real beta) {
    _beta = beta;
}

void Hydro::run(int argc, char* argv[]) {
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
        if (step_cnt % (GNX - 2 * BW) == 0) {
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
            printf("step=%i t=%e dt=%e lmax=%i ngrids=%i avg ngrids=%i", step_cnt, Hydro::get_time(), dt, OctNode::get_max_level(), nnodes, (int) avg_node_cnt);
        }
        if (do_output) {
            if (MPI_rank() == 0) {
                printf("*");
            }
            get_root()->output("X", nint(Hydro::get_time() / OUTPUT_TIME_FREQ), GNX, BW);
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

void Hydro::setup_grid_structure(bool one_iter) {
    Real dt;
    Hydro::set_time(0.0);
    OctNode::get_root()->find_local_nodes();
    OctNode::initialize_grids();
    dt = max_dt_driver();
    dt *= min(1.0, MAXINITDT);
    for (int l = 0; l < OctNode::get_max_level_allowed() + 2; l++) {
        step(dt);
        //  printf( "%e\n", dt );
        Hydro::set_time(0.0);
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

void Hydro::step(Real dt) {
    store();
    Hydro::set_dt(dt);
    Hydro::set_beta(1.0);
    substep_driver();
    Hydro::set_beta(0.25);
    substep_driver();
    Hydro::set_beta(2.0 / 3.0);
    substep_driver();
    Hydro::set_time(Hydro::get_time() + dt);
//    HydroGrid::boundary_driver();
}

void Hydro::substep_driver() {
    DFO = Vector<Real, STATE_NF>(0.0);
    Hydro** list = new Hydro*[get_local_node_cnt()];
    for (int i = 0; i < get_local_node_cnt(); i++) {
        list[i] = dynamic_cast<Hydro*>(get_local_node(i));
    }
    run_program(list, get_local_node_cnt(), cs, GRID_CS_SIZE);
    delete[] list;
    Real tmp[STATE_NF];
    for (int i = 0; i < STATE_NF; i++) {
        tmp[i] = DFO[i];
    }
    Real tmp2[STATE_NF];
    MPI_Allreduce(tmp, tmp2, STATE_NF, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
    for (int i = 0; i < STATE_NF; i++) {
        DFO[i] = tmp2[i];
    }
    DFO[State::et_index] += State::omega * DFO[State::lz_index] + DFO[State::pot_index];
    FO = (FO + DFO * _dt) * _beta + FO0 * (1.0 - _beta);

}


void Hydro::inject_from_children() {
    DFO = Vector<Real, STATE_NF>(0.0);
    Hydro** list = new Hydro*[get_local_node_cnt()];
    for (int i = 0; i < get_local_node_cnt(); i++) {
        list[i] = dynamic_cast<Hydro*>(get_local_node(i));
    }
    run_program(list, get_local_node_cnt(), cs_children, 4);
    delete[] list;

}

Real Hydro::max_dt_driver() {
    Real dt_all, dt;
    eax = 1.0e+99;
    Hydro** list = new Hydro*[get_local_node_cnt()];
    for (int i = 0; i < get_local_node_cnt(); i++) {
        list[i] = dynamic_cast<Hydro*>(get_local_node(i));
    }
    run_program(list, get_local_node_cnt(), es, GRID_ES_SIZE);
    delete[] list;
    dt = eax;
    MPI_Allreduce(&dt, &dt_all, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD);
    dt = dt_all;
    //    printf( "dt=%e\n", dt);
    return dt;
}

void Hydro::boundary_driver() {
    Real dt_all, dt;
    Hydro** list = new Hydro*[get_local_node_cnt()];
    for (int i = 0; i < get_local_node_cnt(); i++) {
        list[i] = dynamic_cast<Hydro*>(get_local_node(i));
    }
    run_program(list, get_local_node_cnt(), es, GRID_ES_SIZE - 2);
    delete[] list;
}

void Hydro::mpi_datatypes_initialize() {
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
