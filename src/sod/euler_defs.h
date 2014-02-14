#ifndef OPTIONS_SOD_
#define OPTIONS_SOD_
#include "../real.h"

//#define USE_RADIATION
#define GRID_DIM            (0.15*1.28)
#define MINMOD_THETA        1.0
#define PPM
#define BW 					(3)
#define EULER_GAMMA         (5.0/3.0)
#define GNX 				(8+2*BW)
#define PNX                 (8+2)
#define TIME_MAX           	(0.12)
#define OUTPUT_TIME_FREQ   	(.01)
#define MAXDTINC            (1.25)
#define GRID_CFL_FACTOR     0.4
#define MAXINITDT           (1.0e-1)

#endif
