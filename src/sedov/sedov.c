
void sed_1d_(const __float128* t, const int* N, const __float128* xpos, const __float128* eblast, const __float128* omega,  const __float128*  geom, const __float128* rho0, const __float128* vel0, const __float128* ener0, const __float128* pres0, const __float128* cs0, const __float128* gam0, __float128* den, __float128* ener, __float128* pres, __float128* vel,__float128 *cs);




void sedov_solution( double t, int N, const double* xpos, double E, double rho,  double* dout, double* eout, double* vout ) {
	__float128 *xpos0, *pout, *csout, *dout0, *eout0, *vout0;
	__float128 vel0, gam0, pres0, cs0, omega, geom, ener0, t0, E0, rho0;
	int i;
	t0 = (__float128) t;
	E0 = (__float128) E;
	rho0 = (__float128) rho;
	eout0 = (__float128*) malloc( sizeof( __float128 ) * N );
	dout0 = (__float128*) malloc( sizeof( __float128 ) * N );
	vout0 = (__float128*) malloc( sizeof( __float128 ) * N );
	xpos0 = (__float128*) malloc( sizeof( __float128 ) * N );
	pout = (__float128*) malloc( sizeof( __float128 ) * N );
	csout = (__float128*) malloc( sizeof( __float128 ) * N );
	vel0 = 0.0;
	omega = 0.0;
	geom = 3.0;
	ener0 = 0.0;
	gam0 = 5.0/3.0;
	pres0 = (gam0 - 1.0)*rho0*ener0;
	cs0 = sqrt( gam0 * pres0 / rho0 );
	for( i = 0; i < N; i++ ) {
		xpos0[i] = (__float128) xpos[i];
	}
	sed_1d_( &t0, &N, xpos0, &E0, &omega, &geom, &rho0, &vel0, &ener0, &pres0, &cs0, &gam0, dout0, eout0, pout, vout0, csout );
	free( xpos0 );
	free( pout );
	free( csout );
	for( i = 0; i < N; i++ ) {
		dout[i] = (double) dout0[i];
		eout[i] = (double) eout0[i];
		vout[i] = (double) vout0[i];
	}
	free( dout0 );
	free( vout0 );
	free( eout0 );
}
