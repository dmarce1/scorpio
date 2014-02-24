#ifndef OPTIONS_____H
#define OPTIONS_____H

//#define BLAST_WAVE
//#define SOD
#define BINARY_STAR

#ifndef ROTATING_DISC
#include "state/state.h"
#else
#include "rotating_disc/state.h"
#endif

#define USE_HYDRO_GRID




#ifdef BINARY_STAR
#define USE_FMM
#define HYDRO_GRAV_GRID
#define NFRAC 2
#include "./binary_star/defs.h"
class BinaryStar;
typedef BinaryStar ProblemGrid;
#endif


#ifdef SINGLE_STAR
#define USE_FMM
#define NFRAC 2
#define HYDRO_GRAV_GRID
//#define USE_FMM
#include "./single_star/defs.h"
class SingleStar;
typedef SingleStar ProblemGrid;
#endif

#ifdef RADIATION_TEST
#define USE_RADIATION
#define EULER_STATE
#include "./radiation_test/defs.h"
class RadiationTest;
typedef RadiationTest ProblemGrid;
#endif

#ifdef POISSON_TEST
#define USE_FMM
#include "./poisson_test/poisson_test_defs.h"
#ifdef USE_FMM
#define USE_HYDRO_GRID
#else
#undef USE_HYDRO_GRID
#endif
class PoissonTest;
typedef PoissonTest ProblemGrid;
#endif

#ifdef BLAST_WAVE
#include "./blast_wave/euler_defs.h"
class GridBlastWave;
typedef GridBlastWave ProblemGrid;
#endif

#ifdef SOD
#include "./sod/euler_defs.h"
class GridSod;
typedef GridSod ProblemGrid;
#endif

#ifdef ROTATING_DISC

#include "rotating_disc/euler_defs.h"
class RotatingDisc;
typedef RotatingDisc ProblemGrid;
#endif

#endif


