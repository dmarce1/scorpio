#include "defs.h"
#include <fenv.h>
#include <stdio.h>
#include "program.h"
#include <stdlib.h>
#include "sod/grid_sod.h"
#include "blast_wave/grid_blast_wave.h"
#include "binary_star/binary_star.h"
#include "FMM/FMM.h"
#include "parameter_reader.h"

void Program::print_help() {
    if (MPI_rank() == 0) {
        printf("Options\n");
        printf("\t-driving_rate - Angular momentum loss rate at which to drive binary into contact.\n");
        printf("\t-driving_time - Number of periods to drive binary.\n");
        printf("\t-help         - Displays this help file.\n");
        printf("\t-M1=f         - The mass of the accretor star for SCF code.\n");
        printf("\t-M2=f         - The mass of the donor star for SCF code.\n");
        printf("\t-maxlev=i     - The maximum refinement level to i.\n");
        printf("\t-maxgrids=i   - The maximum number of sub-grids to i.\n");
        printf("\t-mingrids=i   - The minimum number of sub-grids to i.\n");
        printf("\t-restart=s    - Invokes the evolution code with restart files with prefix s.\n");
        printf("\t-scf          - Invokes the SCF code\n");
    }
}

int Program::run(int argc, char* argv[]) {
    int maxlev;
    if (read_parameter(argc, argv, "help")) {
        print_help();
        return 0;
    }
    if (read_parameter(argc, argv, "maxlev", &maxlev)) {
        if (MPI_rank() == 0) {
            printf("Maximum refinement level is %i\n", maxlev);
        }
    } else {
        if (MPI_rank() == 0) {
            printf("Missing maxlev argument\n");
            print_help();
        }
        return -1;
    }
    OctNode::set_max_level_allowed(maxlev);

    ProblemGrid* root = new ProblemGrid;
    root->init();
    ProblemGrid::run(argc, argv);
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
