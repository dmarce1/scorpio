#include "../defs.h"

#ifdef SOD
#include "grid_sod.h"
#include "../indexer3d.h"
#include "exact_sod.h"

sod_init_t sod_init = { 4.0, 1.0, 1.0, 0.1795, State::gamma };

static Real c = 15.0;
static Real etot = 0.0;
static bool diag_mode = false;
static bool initialized1 = false;
bool GridSod::oblique = false;

Vector<Real, 5>* GridSod::line_avg_tmp;
Vector<Real, 5>* GridSod::line_avg;
int* GridSod::line_bin_cnt;
int GridSod::line_N;
int* GridSod::line_bin_cnt_tmp;

void GridSod::set_refine_flags() {

    const Real r0 = pow(1.0e+5 / 0.49, 1.0 / 5.0) * pow(get_time(), 2.0 / 5.0);
    Real r;
    ChildIndex c;
    if (get_level() < 1) {
        for (int i = 0; i < OCT_NCHILD; i++) {
            set_refine_flag(i, true);
        }
    } else if (get_level() < get_max_level_allowed()) {
        for (int i = 0; i < OCT_NCHILD; i++) {
            set_refine_flag(i, false);
        }
        for (int k = BW; k < GNX - BW; k++) {
            c.set_z(2 * k / GNX);
            for (int j = BW; j < GNX - BW; j++) {
                c.set_y(2 * j / GNX);
                for (int i = BW; i < GNX - BW; i++) {
                    c.set_x(2 * i / GNX);
                    if (!get_refine_flag(c)) {
                        Real drho_max = 0.0;
                        Real rho = U(i, j, k).rho();
                        drho_max = max(drho_max, fabs(log(U(i + 1, j, k).rho() / rho)));
                        drho_max = max(drho_max, fabs(log(U(i - 1, j, k).rho() / rho)));
                        drho_max = max(drho_max, fabs(log(U(i, j + 1, k).rho() / rho)));
                        drho_max = max(drho_max, fabs(log(U(i, j - 1, k).rho() / rho)));
                        drho_max = max(drho_max, fabs(log(U(i, j, k + 1).rho() / rho)));
                        drho_max = max(drho_max, fabs(log(U(i, j, k - 1).rho() / rho)));
                        if (drho_max > 0.05) {
                            set_refine_flag(c, true);
                        }
                    }
                }
            }
        }
    }
}

void GridSod::flux_physical_bounds(int dir) {
}

void GridSod::physical_boundary(int dir) {

    // printf( "Roll your own\n");
    for (dir = 0; dir < 3; dir++) {
        Vector<int, 3> lb, ub, k;
        for (int l = 0; l < 2; l++) {
            if (is_phys_bound(2 * dir + l)) {
                lb = BW;
                ub = GNX - BW - 1;
                lb[dir] = 0 + (GNX - BW) * l;
                ub[dir] = BW - 1 + (GNX - BW) * l;
                for (Indexer3d i(lb, ub); !i.end(); i++) {
                    k = i;
                    k[dir] = BW + l * (GNX - 2 * BW - 1);
                    U(i[0], i[1], i[2]) = U(k[0], k[1], k[2]);
                    U(i[0], i[1], i[2]).set_lz_from_cartesian(X(i[0], i[1], i[2]));
                }
            }
        }
    }
    inc_instruction_pointer();
}

GridSod::GridSod() {
}

GridSod* GridSod::new_octnode() const {
    return new GridSod;
}

Real GridSod::get_output_point(int i, int j, int k, int l) const {
    switch (l) {
    case 0:
        return (*this)(i, j, k).rho();
    case 1:
        return (*this)(i, j, k).vx(X(i, j, k));
    case 2:
        return (*this)(i, j, k).vy(X(i, j, k));
    case 3:
        return (*this)(i, j, k).sz() / (*this)(i, j, k).rho();
    case 4:
        return (*this)(i, j, k).ei(X(i, j, k)) * (State::gamma - 1.0);
    default:
        sod_state_t out;
        exact_sod(&out, &sod_init, xsod(i, j, k), get_time());
        switch (l) {
        case 5:
            return out.rho;
        case 6:
            return out.v;
        case 7:
            return 0.0;
        case 8:
            return 0.0;
        case 9:
            return out.p;
        default:
            abort();
            return 0.0;
        }
        break;
    }
}

int GridSod::nvar_output() const {
    return 5;
}

const char* GridSod::output_field_names(int i) const {
    switch (i) {
    case 0:
        return "d_n";
    case 1:
        return "vx_n";
    case 2:
        return "vy_n";
    case 3:
        return "vz_n";
    case 4:
        return "p_n";
    case 5:
        return "d_a";
    case 6:
        return "vx_a";
    case 7:
        return "vy_a";
    case 8:
        return "vz_a";
    case 9:
        return "p_a";
    default:
        assert(false);
        return "error";
    }
}

void GridSod::run(int argc, char* argv[]) {
    //shadow_on();
    Real dxmin = dynamic_cast<GridSod*>(get_root())->get_dx() / Real(1 << get_max_level_allowed());
    _3Vec O = 0.0;
    O[0] = dxmin/4.0;
   // set_origin(O);
    initialized1 = true;
    PhysicalConstants::A = 0.0;
    line_N = int(2.0 * GRID_DIM / dxmin + 1.5);
    line_avg = new Vector<Real, 5> [line_N];
    line_bin_cnt = new int[line_N];
    line_avg_tmp = new Vector<Real, 5> [line_N];
    line_bin_cnt_tmp = new int[line_N];
    for (int i = 0; i < line_N; i++) {
        line_bin_cnt[i] = 0;
        line_avg[i] = Vector<Real, 5>(0.0);
    }
    Hydro::run(argc, argv);
    GridSod* g;
    const int M = (GNX - 2 * BW) * (GNX - 2 * BW) * (GNX - 2 * BW);
    Real ad[M];
    Real ap[M];
    Real av[M];
    Real xpos[M];
    Real l2 = 0.0;
    dxmin = dynamic_cast<GridSod*>(get_root())->get_dx() / Real(1 << get_max_level_allowed());
    FILE* fpc = fopen("center.dat", "wt");
    for (int l = 0; l < get_local_node_cnt(); l++) {
        g = dynamic_cast<GridSod*>(get_local_node(l));
        int j = 0;
        int k;

        j = BW;
        k = BW;
        if (g->yc(j) < g->get_dx() && g->zc(j) < g->get_dx() && g->yc(j) > 0.0 && g->zc(j) > 0.0) {
            for (int i = BW; i < GNX - BW; i++) {
                if (!g->zone_is_refined(i, j, k)) {
                    fprintf(fpc, "%e %e %e %e\n", g->xc(i), (*g)(i, j, k).rho(), (*g)(i, j, k).vx(g->X(i, j, k)), (*g)(i, j, k).vy(g->X(i, j, k)));
                }
            }
        }

        j = 0;
        const Real dv = pow(g->get_dx(), 3);
        for (Indexer3d i(BW, GNX - BW - 1); !i.end(); i++) {
            if (!g->zone_is_refined(i[0], i[1], i[2])) {
                xpos[j] = g->xsod(i);
                sod_state_t state;
                exact_sod(&state, &sod_init, xpos[j], get_time());
                ad[j] = state.rho;
                av[j] = state.v;
                ap[j] = state.p;
                j++;
            }
        }
        j = 0;
        for (Indexer3d i(BW, GNX - BW - 1); !i.end(); i++) {
            int ri;
            if (!oblique) {
                ri = (xpos[j] + GRID_DIM) / dxmin;
            } else {
                ri = (xpos[j] / sqrt(3.0) + GRID_DIM) / dxmin;
            }
            if (!g->zone_is_refined(i[0], i[1], i[2]) && g->yzsod(i) < 0.1) {
                g->line_avg[ri][0] += (*g)(i).rho();
                if (oblique) {
                    g->line_avg[ri][1] += ((*g)(i).vx(g->X(i)) + (*g)(i).vy(g->X(i)) + (*g)(i).vz()) / sqrt(3.0);
                } else {
                    g->line_avg[ri][1] += (*g)(i).vx(g->X(i));

                }
                g->line_avg[ri][2] += (*g)(i).vy(g->X(i));
                g->line_avg[ri][3] += (*g)(i).vz();
                g->line_avg[ri][4] += (*g)(i).pg(g->X(i));
                g->line_bin_cnt[ri]++;
                l2 += fabs((*g)(i).rho() - ad[j]) * pow(g->get_dx(), 3);
                j++;
            }
        }
    }
    fclose(fpc);
    g = dynamic_cast<GridSod*>(get_root());
    for (int i = 0; i < line_N; i++) {
        g->line_avg_tmp[i] = g->line_avg[i];
        g->line_bin_cnt_tmp[i] = g->line_bin_cnt[i];
    }
    MPI_Allreduce(g->line_bin_cnt_tmp, g->line_bin_cnt, g->line_N, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD );
    MPI_Allreduce(g->line_avg_tmp, g->line_avg, g->line_N * 5, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );
    Real tmp = l2;
    MPI_Allreduce(&tmp, &l2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );
    l2 /= 8.0 * GRID_DIM * GRID_DIM * GRID_DIM;
    printf("L2 error = %e\n", l2);
    diag_mode = true;
    get_root()->output("S", 0.0, GNX, BW);
    if (MPI_rank() == 0) {
        Real* xpos = new double[line_N];
        Real* ad = new double[line_N];
        Real* av = new double[line_N];
        Real* ap = new double[line_N];
        for (int i = 0; i < g->line_N; i++) {
            if (oblique) {
                xpos[i] = ((i + 0.5) * dxmin - GRID_DIM) * sqrt(3.0);
            } else {
                xpos[i] = (i + 0.5) * dxmin - GRID_DIM;
            }
            sod_state_t state;
            exact_sod(&state, &sod_init, xpos[i], get_time());
            ad[i] = state.rho;
            av[i] = state.v;
            ap[i] = state.p;
        }
        FILE* fp = fp = fopen("cmp.d.dat", "wt");
        for (int i = 0; i < g->line_N; i++) {
            if (g->line_bin_cnt[i] != 0) {
                g->line_avg[i] /= Real(g->line_bin_cnt[i]);
                fprintf(fp, "%e %e %e\n", -xpos[i], g->line_avg[i][0], ad[i]);
            }
        }
        fclose(fp);
        fp = fp = fopen("cmp.v.dat", "wt");
        for (int i = 0; i < g->line_N; i++) {
            if (g->line_bin_cnt[i] != 0) {
                fprintf(fp, "%e %e %e\n", -xpos[i], g->line_avg[i][1], av[i]);
            }
        }
        fclose(fp);
        fp = fp = fopen("cmp.p.dat", "wt");
        for (int i = 0; i < g->line_N; i++) {
            if (g->line_bin_cnt[i] != 0) {
                fprintf(fp, "%e %e %e\n", -xpos[i], g->line_avg[i][4], ap[i]);
            }
        }
        fclose(fp);
        fp = fp = fopen("cmp.u.dat", "wt");
        for (int i = 0; i < g->line_N; i++) {
            if (g->line_bin_cnt[i] != 0) {
                fprintf(fp, "%e %e %e\n", -xpos[i], g->line_avg[i][4] / (State::gamma - 1.0) / g->line_avg[i][0], ap[i] / (State::gamma - 1.0) / ad[i]);
            }
        }
        fclose(fp);
        fp = fp = fopen("cmp.s.dat", "wt");
        for (int i = 0; i < g->line_N; i++) {
            if (g->line_bin_cnt[i] != 0) {
                fprintf(fp, "%e %e %e\n", -xpos[i], log(pow(g->line_avg[i][4] / (State::gamma - 1.0), (1.0 / State::gamma)) / g->line_avg[i][0]),
                        log(pow(ap[i] / (State::gamma - 1.0), 1.0 / State::gamma) / ad[i]));
            }
        }
        fclose(fp);
        //printf("%e\n", etot);
    }
}

void GridSod::initialize() {
#ifdef SOD
    //  State::set_gamma(7.0 / 5.0);
    State U = Vector<Real, STATE_NF>(0.0);
    Real et;
    for (int k = 0; k < GNX; k++) {
        for (int j = 0; j < GNX; j++) {
            for (int i = 0; i < GNX; i++) {
                if (xsod(i, j, k) < -0.00) {
                    U.set_rho(sod_init.rhol);
                    U.set_et(sod_init.pl / (sod_init.gamma - 1.0));
                } else {
                    U.set_rho(sod_init.rhor);
                    U.set_et(sod_init.pr / (sod_init.gamma - 1.0));
                }
                U.set_tau(pow(U[State::et_index], State::gamma));
                (*this)(i, j, k) = U;
            }
        }
    }
#endif
}

GridSod::~GridSod() {
}

#endif
