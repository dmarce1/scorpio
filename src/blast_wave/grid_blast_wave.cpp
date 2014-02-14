#include "../defs.h"

#ifdef BLAST_WAVE
#include "grid_blast_wave.h"
#include "sedov/sedov.h"
#include "../indexer3d.h"
static Real c = 15.0;
static Real etot = 0.0;
static bool diag_mode = false;
static bool initialized1 = false;

Vector<Real, 4>* GridBlastWave::radial_avg_tmp;
Vector<Real, 4>* GridBlastWave::radial_avg;
int* GridBlastWave::radial_bin_cnt;
int GridBlastWave::radial_N;
int* GridBlastWave::radial_bin_cnt_tmp;

void GridBlastWave::set_refine_flags() {
	if (using_shadow()) {
		printf("set_refine_flags not defined\n");
		abort();
	}
	const Real r0 = pow(1.0e+5 / 0.49, 1.0 / 5.0) * pow(get_time(), 2.0 / 5.0);
	Real r;
	ChildIndex c;
	if (get_level() < max(1, this->get_max_level_allowed() - 2)) {
		for (int i = 0; i < OCT_NCHILD; i++) {
			set_refine_flag(i, true);
		}
	} else if (get_level() < get_max_level_allowed()) {
		for (int k = BW; k < GNX - BW; k++) {
			c.set_z(2 * k / GNX);
			for (int j = BW; j < GNX - BW; j++) {
				c.set_y(2 * j / GNX);
				for (int i = BW; i < GNX - BW; i++) {
					c.set_x(2 * i / GNX);
					if (!get_refine_flag(c)) {
						r = sqrt(X(i, j, k).dot(X(i, j, k)));
						if (fabs(r - r0) < (GNX - 2 * BW + 1) * get_dx() * GRID_CFL_FACTOR) {
							set_refine_flag(c, true);
						}
					}
				}
			}
		}
	}
}

GridBlastWave::GridBlastWave() {
}

GridBlastWave* GridBlastWave::new_octnode() const {
	return new GridBlastWave;
}

Real GridBlastWave::get_output_point(int i, int j, int k, int l) const {
	switch (l) {
	case 0:
		return (*this)(i, j, k).rho();
	case 1:
		return (*this)(i, j, k).sx() / (*this)(i, j, k).rho();
	case 2:
		return (*this)(i, j, k).sy() / (*this)(i, j, k).rho();
	case 3:
		return (*this)(i, j, k).sz() / (*this)(i, j, k).rho();
	case 4:
		return (*this)(i, j, k).ei(X(i,j,k)) / (*this)(i, j, k).rho();
#ifdef USE_RADIATION
	case 5:
		return  (*this)(i, j, k).er();
#endif
	default:
		assert(false);
		return 0.0;
	}
}

int GridBlastWave::nvar_output() const {
#ifdef USE_RADIATION
	return 5 + 1;
#else
	return 5;
#endif
}

const char* GridBlastWave::output_field_names(int i) const {
	switch (i) {
	case 0:
		return "d";
	case 1:
		return "vx";
	case 2:
		return "vy";
	case 3:
		return "vz";
	case 4:
		return "epsilon";
#ifdef USE_RADIATION
	case 5:
		return "erad";
#endif
	default:
		assert(false);
		return "error";
	}
}

void GridBlastWave::run(int argc, char* argv[]) {
	shadow_off();
	Real dxmin = dynamic_cast<GridBlastWave*>(get_root())->get_dx() / Real(1 << get_max_level_allowed());
	initialized1 = true;
	radial_N = int(2.0 * GRID_DIM / dxmin + 1.5);
	radial_avg = new Vector<Real, 4> [radial_N];
	radial_bin_cnt = new int[radial_N];
	radial_avg_tmp = new Vector<Real, 4> [radial_N];
	radial_bin_cnt_tmp = new int[radial_N];
	for (int i = 0; i < radial_N; i++) {
		radial_bin_cnt[i] = 0;
		radial_avg[i] = Vector<Real, 4>(0.0);
	}
#ifdef USE_RADIATION
	HydroRadGrid::run(argc, argv);
#else
	HydroGrid::run(argc, argv);
#endif
	etot = dynamic_cast<GridBlastWave*>(get_root())->state_sum()[0].et();
	printf("Etot = %e\n", etot);
	Real l1, l2, d, tmp;
	GridBlastWave* g;
	const int M = (GNX - 2 * BW) * (GNX - 2 * BW) * (GNX - 2 * BW);
	Real ad[M];
	Real ae[M];
	Real av[M];
	Real rpos[M];
	dxmin = dynamic_cast<GridBlastWave*>(get_root())->get_dx() / Real(1 << get_max_level_allowed());
	int ri;
	for (int l = 0; l < get_local_node_cnt(); l++) {
		g = dynamic_cast<GridBlastWave*>(get_local_node(l));
		g->analytic.allocate();
		int j = 0;
		for (Indexer3d i(BW, GNX - BW - 1); !i.end(); i++) {
			if (!g->zone_is_refined(i[0], i[1], i[2])) {
				rpos[j++] = sqrt(g->X(i).dot(g->X(i)));
			}
		}
		sedov_solution(get_time(), j, rpos, etot, 1.0, ad, ae, av);
		j = 0;
		const Real dv = pow(g->get_dx(), 3);
		for (Indexer3d i(BW, GNX - BW - 1); !i.end(); i++) {
			if (!g->zone_is_refined(i[0], i[1], i[2])) {
				//printf("%e %e %e\n", sqrt(g->X(i).dot(g->X(i))), ad[j], (*g)(i).rho());
				g->analytic(i).set_rho(ad[j]);
				d = ad[j] - (*g)(i).rho();
				l1 += fabs(d) * dv;
				l2 += d * d * dv;
				ri = rpos[j] / dxmin;
				g->radial_avg[ri][0] += (*g)(i).rho();
			//	g->radial_avg[ri][1] += (*g)(i).pg();
			//	g->radial_avg[ri][2] += (*g)(i).speed();
			//	g->radial_avg[ri][3] += (*g)(i).ei() / (*g)(i).rho();
			//	g->radial_bin_cnt[ri]++;
				j++;
			}
		}
	}
	g = dynamic_cast<GridBlastWave*>(get_root());
	for (int i = 0; i < radial_N; i++) {
		g->radial_avg_tmp[i] = g->radial_avg[i];
		g->radial_bin_cnt_tmp[i] = g->radial_bin_cnt[i];
	}
	MPI_Allreduce(g->radial_bin_cnt_tmp, g->radial_bin_cnt, g->radial_N, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce(g->radial_avg_tmp, g->radial_avg, g->radial_N * 4, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );
	tmp = l2;
	MPI_Allreduce(&tmp, &l2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD );
	l2 /= 8.0 * GRID_DIM * GRID_DIM * GRID_DIM;
	printf("L2 error = %e\n", l2);
	diag_mode = true;
	get_root()->output("S", 0.0, GNX, BW);
	if (MPI_rank() == 0) {
		Real* rpos = new double[radial_N];
		Real* ad = new double[radial_N];
		Real* av = new double[radial_N];
		Real* ae = new double[radial_N];
		for (int i = 0; i < g->radial_N; i++) {
			rpos[i] = (i + 0.5) * dxmin;
		}
		sedov_solution(get_time(), radial_N, rpos, etot, 1.0, ad, ae, av);
		FILE* fp = fp = fopen("cmp.d.dat", "wt");
		for (int i = 0; i < g->radial_N; i++) {
			g->radial_avg[i] /= Real(g->radial_bin_cnt[i]);
			fprintf(fp, "%e %e %e\n", (i + 0.5) * dxmin, g->radial_avg[i][0], ad[i]);
		}
		fclose(fp);
		fp = fp = fopen("cmp.v.dat", "wt");
		for (int i = 0; i < g->radial_N; i++) {
			fprintf(fp, "%e %e %e\n", (i + 0.5) * dxmin, g->radial_avg[i][2], av[i]);
		}
		fclose(fp);
		fp = fp = fopen("cmp.p.dat", "wt");
		for (int i = 0; i < g->radial_N; i++) {
			fprintf(fp, "%e %e %e\n", (i + 0.5) * dxmin, g->radial_avg[i][1], ae[i] * ad[i] * (State::gamma - 1.0));
		}
		fclose(fp);
		fp = fp = fopen("cmp.e.dat", "wt");
		for (int i = 0; i < g->radial_N; i++) {
			fprintf(fp, "%e %e %e\n", (i + 0.5) * dxmin, g->radial_avg[i][3], ae[i]);
		}
		fclose(fp);
		//printf("%e\n", etot);
	}
}

void GridBlastWave::initialize() {
#ifdef BLAST_WAVE
	State U = Vector<Real, STATE_NF>(0.0);
	Real r1, this_etot, dv, tmp, v0, r0, r2, dxmin;
	dxmin = dynamic_cast<GridBlastWave*>(get_root())->get_dx() / Real(1 << get_max_level_allowed());
	r0 = 3.5 * dxmin;
	this_etot = 0.0;
	dv = pow(get_dx(), 3);
	for (int k = 0; k < GNX; k++) {
		for (int j = 0; j < GNX; j++) {
			for (int i = 0; i < GNX; i++) {
				r1 = sqrt((xc(i)) * (xc(i)) + (yc(j)) * (yc(j)) + zc(k) * zc(k));
				U.set_rho(1.0);
				v0 = r0 * r0 * r0 * 4.0 / 3.0 * M_PI;
				r2 = max(r0, get_dx());
				if (r1 < r2) {
					U.set_et(1.0e+5 / v0 * 1.0e+5 / 8.908965e+04);
				} else {
#ifdef USE_RADIATION
					U.set_et(1.0);
#endif
				}
				U.set_tau(pow(U.et(), 1.0 / State::gamma));
#ifdef USE_RADIATION
				U.set_er(1.0e-3);
#endif
				(*this)(i, j, k) = U;
				if (!zone_is_refined(i, j, k)) {
					etot += U.et() * dv;
				}
			}
		}
	}
	tmp = etot;
	etot += this_etot;
#endif
}

GridBlastWave::~GridBlastWave() {
}

#endif
