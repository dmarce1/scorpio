#ifndef OPTIONS_____H
#define OPTIONS_____H

//#define BLAST_WAVE
//#define SOD
#define BINARY_STAR

#include "./state/state.h"

#ifdef BINARY_STAR
#define USE_FMM
#define HYDRO_GRAV_GRID
#define NFRAC 5
#include "./binary_star/defs.h"
class BinaryStar;
typedef BinaryStar ProblemGrid;
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


#endif


