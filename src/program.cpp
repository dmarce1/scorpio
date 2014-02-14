#include "defs.h"
#include <fenv.h>
#include <stdio.h>
#include "program.h"
#include <stdlib.h>
#include "sod/grid_sod.h"
#include "blast_wave/grid_blast_wave.h"
#include "binary_star/binary_star.h"
#include "hydro_FMM_grid/hydro_FMM_grid.h"

int Program::run(int argc, char* argv[]) {
    if( argc > 1 ) {
	OctNode::set_max_level_allowed(atoi(argv[1]));
    } else {
        printf( "Missing maxlevel argument\n");
        abort();
    }


	ProblemGrid* root = new ProblemGrid;
	root->init();
	ProblemGrid::run(argc,argv);
	delete root;
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
