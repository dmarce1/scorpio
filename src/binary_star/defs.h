#ifndef OPTIONS_SCF_
#define OPTIONS_SCF_
#include <mpi.h>

#include "../real.h"

//#define REFINE_ACC_MORE
//#define USE_FMM
#define DRIVING_TIME        0.0




#ifdef USE_FMM
//#define  RANK_ZERO_HAS_ONE_GRID
#define RANK0_OMP_THREADS 1
#endif

#define DRIVING             0.01
#define GCON 1.0
#define DIAGNOSTIC_FREQ 1
#define DONOR_FILL 0.99
#define INX                 12
#define GRID_DIM            (1.0)
#define MINMOD_THETA        1.3
#define PPM
#define BW 					(3)
#define EULER_GAMMA         (5.0/3.0)
#define GNX 				(INX+2*BW)
#define PNX 				(INX+2)
#define TIME_MAX           	(1.0e+99)
#define OUTPUT_TIME_FREQ   	(1.0/30.0)
#define MAXDTINC            (1.25)
#define GRID_CFL_FACTOR     (1.0/6.0)
#define MAXINITDT           (1.0e-3)
#define ZTWD
//#define NGRID_LIMIT         1000

#endif
