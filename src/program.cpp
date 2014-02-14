#include "defs.h"
#include <fenv.h>
#include <stdio.h>
#include "program.h"
#include <stdlib.h>
#include "sod/grid_sod.h"
#include "blast_wave/grid_blast_wave.h"
#include "rotating_disc/rotating_disc.h"
#include "single_star/single_star.h"
#include "binary_star/binary_star.h"
#include "hydro_FMM_grid/hydro_FMM_grid.h"

int Program::run(int argc, char* argv[]) {
	OctNode::set_max_level_allowed(atoi(argv[1]));
#ifdef POISSON_TEST
	PoissonTest::x0 = atof(argv[2]);
	PoissonTest::y0 = atof(argv[3]);
	PoissonTest::z0 = atof(argv[4]);
#endif


#ifndef ROTATING_DISC
	ProblemGrid* root = new ProblemGrid;
	root->init();
	ProblemGrid::run(argc,argv);
	delete root;
#else
	if (MPI_rank() == 0)
		printf("parsing arguments\n");
	int lor = atoi(argv[1]);
	double x_in = atof(argv[2]);
	int kick_mode = atoi(argv[3]);
	double num_orbits = atof(argv[4]);
	printf("argv[5] =");
	printf(argv[5]);
	printf("\n");

	std::string ang_mom_string = argv[5];

	if (MPI_rank() == 0)
		printf("setting max_refine_level...\n");
	OctNode::set_max_level_allowed(lor);
	RotatingDisc* root = new RotatingDisc;
	root->init();
	root->set_x_in(x_in);
	root->set_kick_mode(kick_mode);
	double output_time_freq = root->get_period();
	output_time_freq *= OUTPUT_TIME_FREQ;
	root->set_output_frequency_by_time(output_time_freq);
	if (ang_mom_string == "cart") {
		if (MPI_rank() == 0)
			printf("NOT conserving angular momentum...\n");
		root->set_ang_mom_cons(false);
	} else if (ang_mom_string == "cyl") {
		if (MPI_rank() == 0)
			printf("conserving angular momentum...\n");
		root->set_ang_mom_cons(true);
	} else {
		if (MPI_rank() == 0)
			printf("invalid momentum argument\n");
		abort();
	}
	if (MPI_rank() == 0)

		printf("done setting ang_mom_cons...\n");

	double time_max = root->get_period();
	time_max *= num_orbits;
	root->step_to_time(time_max);
	delete root;
#endif

	return 0;
}

Program::Program() {
#ifndef NDEBUG
	feenableexcept (FE_DIVBYZERO);
	feenableexcept (FE_INVALID);
	feenableexcept (FE_OVERFLOW);
#endif
}

Program::~Program() {
}
