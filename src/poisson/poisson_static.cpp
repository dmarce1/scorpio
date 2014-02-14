#include "poisson.h"

Real** Poisson::i_poles;
Real** Poisson::r_poles;
Real Poisson::com[4];

_3Vec Poisson::get_center_of_mass() {
	_3Vec v;
	for (int i = 0; i < 3; i++) {
		v[i] = com[i + 1];
	}
	return v;
}

void Poisson::compute_com() {
	Real tmp[4];
	for (int i = 0; i < 4; i++) {
		com[i] = 0.0;
	}
	for (int i = 0; i < get_local_node_cnt(); i++) {
		dynamic_cast<Poisson*>(OctNode::get_local_node(i))->accumulate_com();
	}
	for (int i = 0; i < 4; i++) {
		tmp[i] = com[i];
	}
	MPI_Allreduce(tmp, com, 4, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
	for (int i = 1; i < 4; i++) {
		com[i] /= com[0];
	}
}

void Poisson::poles_clear() {
	pole_init();
	int l, m;
	for (l = 0; l <= LMAX; l++) {
		for (m = -l; m <= l; m++) {
			r_poles[l][m] = i_poles[l][m] = 0.0;
		}
	}
}

void Poisson::pole_init() {
	static bool initialized = false;
	if (!initialized) {
		i_poles = new Real*[LMAX + 1];
		r_poles = new Real*[LMAX + 1];
		for (int l = 0; l <= LMAX; l++) {
			i_poles[l] = new Real[2 * l + 1];
			r_poles[l] = new Real[2 * l + 1];
			i_poles[l] += l;
			r_poles[l] += l;
		}
		initialized = true;
	}

}

void Poisson::poles_reduce() {
	int l, m, cnt;
	Real* tmp1, *tmp2;
	tmp1 = new Real[2 * (LMAX + 1) * (LMAX + 1)];
	tmp2 = new Real[2 * (LMAX + 1) * (LMAX + 1)];
	cnt = 0;
	for (l = 0; l <= LMAX; l++) {
		for (m = -l; m <= l; m++) {
			tmp1[cnt++] = r_poles[l][m];
			tmp1[cnt++] = i_poles[l][m];
		}
	}
	assert( cnt == 2 * (LMAX + 1) * (LMAX + 1));
	MPI_Allreduce(tmp1, tmp2, cnt, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
	cnt = 0;
	for (l = 0; l <= LMAX; l++) {
		for (m = -l; m <= l; m++) {
			r_poles[l][m] = tmp2[cnt++];
			i_poles[l][m] = tmp2[cnt++];
		}
	}
	delete[] tmp1;
	delete[] tmp2;
}

void Poisson::compute_physical_boundaries() {
	int count;
	Poisson::poles_clear();
	Poisson::compute_com();
	count = OctNode::get_local_node_cnt();
	for (int i = 0; i < count; i++) {
		dynamic_cast<Poisson*>(OctNode::get_local_node(i))->poles_compute();

	}
	Poisson::poles_reduce();
	for (int i = 0; i < count; i++) {
		dynamic_cast<Poisson*>(OctNode::get_local_node(i))->compute_local_physical_boundaries();
	}
}

void Poisson::run() {
	int iters = 0;
	iters = 0;
	Real boundary_time, solve_time, start_time;
	get_root()->find_local_nodes();
	for (int iters = 0; iters <= OctNode::get_max_level_allowed(); iters++) {
		OctNode::initialize_grids();
		check_for_refine();
		if (MPI_rank() == 0) {
			printf("maxlevel = %i nodes = %i\n", get_max_level(), get_node_cnt());
		}
	}
	OctNode::initialize_grids();

	start_time = MPI_Wtime();
	compute_physical_boundaries();
	boundary_time = MPI_Wtime() - start_time;

	start_time = MPI_Wtime();
	Real err = vcycle();
	for (iters = 0; err > 1.0e-6; iters++) {
		err = vcycle();
		if (MPI_rank() == 0) {
			printf("%i %e\n", iters, err);
		}
	}
	solve_time = MPI_Wtime() - start_time;

	if (MPI_rank() == 0) {
		printf("Boundary time = %e : Solve time = %e\n", boundary_time, solve_time);
	}

	get_root()->output("X", 0, PNX, 1);
}

