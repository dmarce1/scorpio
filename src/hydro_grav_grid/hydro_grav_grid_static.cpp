#include "hydro_grav_grid.h"
#ifdef HYDRO_GRAV_GRID

Real HydroGravGrid::hydro_time, HydroGravGrid::poisson_boundary_time, HydroGravGrid::poisson_interior_time;
Real HydroGravGrid::poisson_tolerance = 1.0e-6;

void HydroGravGrid::set_poisson_tolerance(Real e) {
    poisson_tolerance = e;
}

void HydroGravGrid::solve_poisson() {
    Real err;
    int iters;
    Real start_time = MPI_Wtime();
    set_gravity_source();
    compute_physical_boundaries();
    set_poisson_tolerance(Poisson::source_sum() / 1.0e+10);
//	printf("%e\n", Poisson::source_sum() / 1.0e+9);
    poisson_boundary_time += MPI_Wtime() - start_time;
    start_time = MPI_Wtime();
    iters = 0;
    do {
        err = vcycle();
        if (MPI_rank() == 0) {
            if (iters > 100) {
                if (iters != 0) {
                    printf("\t       ");
                }
                printf(" iteration=%i residual=%e\n", iters + 1, err);
            }
        }
        iters++;
        if (iters == 1000) {
            if (MPI_rank() == 0) {
                printf("poisson error\n");
            }
            get_root()->output("err", 0, GNX, BW);
            if (get_time() != 0.0) {
                MPI_Barrier(MPI_COMM_WORLD );
                abort();
            }
            break;
        }
    } while (err > poisson_tolerance);
//	compute_faces();
    pot_to_hydro_grid();
    poisson_interior_time += MPI_Wtime() - start_time;
}

bool HydroGravGrid::check_for_refine() {
    bool rc;
    to_conserved_energy();
    rc = OctNode::check_for_refine();
    inject_from_children();
    if (rc) {
        pot_to_hydro_grid();
        HydroGrid::redistribute_grids();
        pot_from_hydro_grid();
        int count = OctNode::get_local_node_cnt();
        for (int i = 0; i < count; i++) {
            dynamic_cast<MultiGrid*>(OctNode::get_local_node(i))->phi_calc_amr_bounds();
        }
        if (MPI_rank() == 0) {
            //		printf("\tresolve");
        }
        if (get_time() != 0.0) {
            solve_poisson();
        }
    }
    from_conserved_energy();
    return rc;

}

void HydroGravGrid::to_conserved_energy() {
    HydroGrid* g;
    for (int n = 0; n < get_local_node_cnt(); n++) {
        g = dynamic_cast<HydroGrid*>(get_local_node(n));
        for (int k = BW - 1; k < GNX - BW + 1; k++) {
            for (int j = BW - 1; j < GNX - BW + 1; j++) {
                for (int i = BW - 1; i < GNX - BW + 1; i++) {
                    (*g)(i, j, k).to_con(g->X(i, j, k));
                }
            }
        }
    }
}

void HydroGravGrid::from_conserved_energy() {
    HydroGrid* g;
    for (int n = 0; n < get_local_node_cnt(); n++) {
        g = dynamic_cast<HydroGrid*>(get_local_node(n));
        for (int k = BW - 1; k < GNX - BW + 1; k++) {
            for (int j = BW - 1; j < GNX - BW + 1; j++) {
                for (int i = BW - 1; i < GNX - BW + 1; i++) {
                    (*g)(i, j, k).from_con(g->X(i, j, k));
                }
            }
        }
    }
}

void HydroGravGrid::pot_to_hydro_grid() {
    const int o = BW - 1;
    Poisson* p;
    HydroGrid* g;
    Real pot;
    for (int n = 0; n < get_local_node_cnt(); n++) {
        p = dynamic_cast<Poisson*>(get_local_node(n));
        g = dynamic_cast<HydroGrid*>(get_local_node(n));
        for (int k = BW - 1; k < GNX - BW + 1; k++) {
            for (int j = BW - 1; j < GNX - BW + 1; j++) {
                for (int i = BW - 1; i < GNX - BW + 1; i++) {
                    pot = p->get_phi(i - o, j - o, k - o);
                    pot *= (*g)(i, j, k).rho();
                    pot += (*g)(i, j, k).rot_pot(g->X(i, j, k));
                    (*g)(i, j, k).set_pot(pot);
                }
            }
        }
    }
}

void HydroGravGrid::pot_from_hydro_grid() {
    const int o = BW - 1;
    Poisson* p;
    HydroGrid* g;
    Real phi;
    for (int n = 0; n < get_local_node_cnt(); n++) {
        p = dynamic_cast<Poisson*>(get_local_node(n));
        g = dynamic_cast<HydroGrid*>(get_local_node(n));
        for (int k = BW; k < GNX - BW; k++) {
            for (int j = BW; j < GNX - BW; j++) {
                for (int i = BW; i < GNX - BW; i++) {
                    phi = (*g)(i, j, k).pot();
                    phi -= (*g)(i, j, k).rot_pot(g->X(i, j, k));
                    phi /= (*g)(i, j, k).rho();
                    p->set_phi(i - o, j - o, k - o, phi);
                }
            }
        }
    }
}

void HydroGravGrid::run(int argc, char* argv[]) {
    Real dt;
    bool do_output, last_step;
    int step_cnt = 0;
    int ostep_cnt = 0;
    int nnodes;

    setup_grid_structure();

    hydro_time = poisson_boundary_time = poisson_interior_time = 0.0;
    if (OUTPUT_TIME_FREQ <= TIME_MAX) {
        get_root()->output("X", 0.0, GNX, BW);
    }
    do {
        dt = next_dt(&do_output, &last_step, &ostep_cnt);
        step(dt);
        if (step_cnt % (GNX - 2 * BW)== 0){
            check_for_refine();
        }
        step_cnt++;
        if (MPI_rank() == 0) {
            nnodes = OctNode::get_node_cnt();
            printf("step=%i t=%e dt=%e lmax=%i ngrids=%i", step_cnt, HydroGrid::get_time(), dt, OctNode::get_max_level(), nnodes);
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
        Vector<State, 4> vec = dynamic_cast<HydroGrid*>(get_root())->state_sum();
        State sum = vec[0];
        _3Vec com = get_center_of_mass();
        if (MPI_rank() == 0) {
            FILE* fp = fopen("sums.dat", "at");
            fprintf(fp, "%e %e %e %e %e %e %e %e\n", get_time(), sum[0], sum[1], sum[2], sum[3], sum[4] + 0.5 * sum[6], sum[State::frac_index + 0],
                    sum[State::frac_index + 1]);
            fclose(fp);
            fp = fopen("com.dat", "at");
            fprintf(fp, "%e %e %e %e\n", get_time(), com[0], com[1], com[2]);
            fclose(fp);
        }
    } while (!last_step);
    if (MPI_rank() == 0) {
        char* str;
        if (asprintf(&str, "time.%i.txt", get_max_level_allowed()) == 0) {
            printf("Unable to create filename\n");
        } else {
            FILE* fp = fopen(str, "at");
            fprintf(fp, "%i %e %e %e %e\n", MPI_size(), hydro_time + poisson_boundary_time + poisson_interior_time, hydro_time, poisson_boundary_time,
                    poisson_interior_time);
            fclose(fp);
            free(str);
        }
    }
}

void HydroGravGrid::set_gravity_source() {
    const int o = BW - 1;
    Poisson* p;
    HydroGrid* g;
    for (int n = 0; n < get_local_node_cnt(); n++) {
        p = dynamic_cast<Poisson*>(get_local_node(n));
        g = dynamic_cast<HydroGrid*>(get_local_node(n));
        for (int k = BW; k < GNX - BW; k++) {
            for (int j = BW; j < GNX - BW; j++) {
                for (int i = BW; i < GNX - BW; i++) {
                    p->set_source(i - o, j - o, k - o, 4.0 * M_PI * PhysicalConstants::G * (*g)(i, j, k).rho());
                }
            }
        }
    }
}

void HydroGravGrid::step(Real dt) {
    Real start_time;
    Real beta[3] = { 1.0, 0.25, 2.0 / 3.0 };
    HydroGrid::set_dt(dt);
    store();
    for (int i = 0; i < 3; i++) {
        if (MPI_rank() == 0) {
            //  printf("\t rk = %i", i + 1);
        }
        HydroGrid::set_beta(beta[i]);
        start_time = MPI_Wtime();
        substep_driver();
        hydro_time += MPI_Wtime() - start_time;
        solve_poisson();
    }
    set_time(get_time() + dt);
}

void HydroGravGrid::setup_grid_structure() {
    Real dt;
    HydroGrid::set_time(0.0);
    OctNode::get_root()->find_local_nodes();
    OctNode::initialize_grids();
    dt = max_dt_driver();
    dt *= min(1.0, MAXINITDT);
    if (MPI_rank() == 0) {
        printf("\t       ");
    }
    solve_poisson();
    for (int l = 0; l <= OctNode::get_max_level_allowed(); l++) {
        //    HydroGravGrid::step(dt);
        HydroGrid::set_time(0.0);
        check_for_refine();
        if (MPI_rank() == 0) {
            printf("%i %i\n", OctNode::get_max_level(), OctNode::get_node_cnt());
        }
        OctNode::initialize_grids();
        dt = min(max_dt_driver(), MAXINITDT);
    }
    if (MPI_rank() == 0) {
        printf("\t       ");
    }
    solve_poisson();

}

#endif
