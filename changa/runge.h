/*
 * Runge-Kutta integrator originally written for PKDGRAV by Thomas
 * Quinn
 */
void
RungeKutta(void *CTX, 
	   void (*deriv)(void *, double, double *, double*),
	   int nDep,		/* number of dependent variables */
	   double tin,		/* independent variable */
	   double *xin,		/* array of input */
	   double tout,
	   double *xout,		/* array of output */
	   int nSteps);
