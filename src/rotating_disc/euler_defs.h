#ifndef OPTIONS_SCF_
#define OPTIONS_SCF_
#include "../real.h"

//#define DISCONTINUITY_DETECT
//#define SCF_CODE
//#define MEGABYTES (256)
//#define MAX_POINTS ((MEGABYTES*1024*1024)/1823)
#define MAX_POINTS (1024*1024)
#define DYNAMIC_RANGE 1.0e+2
#ifndef SCF_CODE
//#define READ_FROM_FILE "X.hello.chk"
#endif
//#define BENCH 64
//#defineR0
#define GRAVITY_EULER_STATE
#define USE_POISSON_DRIVE
//#define USE_HYDRO_DRIVE
#define MINMOD_THETA        1.3
//#define Z_REFLECT
//#define POISSON_TEST
#define FRAME_RATE 50
#define PPM
#define BW 					(3)
#define GNX 				(8+2*BW)
#define RK_ORDER          	(3)
#define MAXLEVEL           	(4)
//#define OUTPUT_STEP_FREQ    1
#define MAX_GRIDS           1024
#define MIN_GRIDS            64
//#define OUTPUT_TIME_FREQ   	(2.0*M_PI/1.097069/100.0)
#define DEN_REF
#define REL_REF
#define VIRIAL_TOLERANCE 1.0e-6

#define KL  KR
#define KR 0.130821412355133

#define DYNAMIC_OMEGA

#define MAXDTINC           (1.25)
#define GRID_CFL_FACTOR    0.4
#ifdef SCF_CODE
#define CHKPT_FREQ 10
#else
#define CHKPT_FREQ         (32)
#endif

#define MIRROR_REFINE_Y
//#define MIRROR_REFINE_X


#define MAXINITDT          (1.0e-2)
#define ELL_TOLER          (1.0e-4)
#define ELL_TOLER2          (1.0)


#define KAPPA               (1.0)
#define M_C                 (1.0)
#define BIG_G               (1.0)
#define EULER_GAMMA         (1.333) // = 1+1/n for equilibrium
#define GRID_DIM            (1.5)// goes from  -grid_dim to grid_dim

#define R_OUTER             (1.0)  //Outer edge of torus
#define POLYTROPIC_N        (3.0)   //Determines initial structure
#define A_0                 (1e-2)  //Perturb magnitude

//I would like these, plus omega, to be dynamically calculated at runtime.
#define OUTPUT_TIME_FREQ    (1e-2)  // in orbits
#define TIME_MAX            (40) // in orbits


//#define ANG_MOM_CONS
//#define NO_ROTATE

void init_scf(const char*);
Real scf_phi(Real, Real, Real, Real);
Real scf_rho(Real, Real, Real, Real);
Real scf_et(Real, Real, Real, Real);
void scf_com(Real, Real, Real);


#endif
