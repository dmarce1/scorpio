#ifndef EXACT_SOD__C
#define EXACT_SOD__C


typedef struct {
	double rhol, rhor, pl, pr, gamma;
} sod_init_t;

typedef struct {
	double rho, v, p;
} sod_state_t;

void exact_sod(sod_state_t* out, sod_init_t* in, double x, double t);



#endif
