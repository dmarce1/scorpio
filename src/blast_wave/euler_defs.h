#ifndef OPTIONS_SCF_
#define OPTIONS_SCF_
#include "../real.h"

//#define USE_RADIATION
#define GRID_DIM            (1.28*5.0)
#define MINMOD_THETA        1.3
#define PPM
#define BW 					(3)
#define EULER_GAMMA         (5.0/3.0)
#define GNX 				(8+2*BW)
#define PNX                 (8+2)
#define TIME_MAX           	(1.00e-2)
#define OUTPUT_TIME_FREQ   	(.0025)
#define MAXDTINC            (1.25)
#define GRID_CFL_FACTOR     0.4
#define MAXINITDT           (1.0e-1)

#endif
