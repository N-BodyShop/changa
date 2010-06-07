/*
 * Runge-Kutta integrator originally written for PKDGRAV by Thomas
 * Quinn
 */
#include <math.h>
#include <assert.h>

#define MAXDEP 10
/*
 * Make one R-K step
 */
void
RungeStep(void *CTX,		/* context for extra parameters */
	  /* derivative function takes independent and array of
	     dependent variables, last argument is returned derivates */
	  void (*deriv)(void *, double, double *, double*),
	  int nDep,		/* number of dependent variables */
	  double t,		/* independent variable */
	  double *xin,		/* array of input */
	  double *xout,		/* array of output */
	  double dDelta)	/* step size */
{
    double k1[MAXDEP];
    double k2[MAXDEP];
    double k3[MAXDEP];
    double k4[MAXDEP];
    double xtemp[MAXDEP];
    const double onesixth = 1.0/6.0;
    const double onethird = 1.0/3.0;
    int i;
    
    assert(nDep <= MAXDEP);
    deriv(CTX, t, xin, k1);
    for(i = 0; i < nDep; i++) {
	k1[i] *= dDelta;
	xtemp[i] = xin[i] + 0.5*k1[i];
	}
    deriv(CTX, t + 0.5*dDelta, xtemp, k2);
    for(i = 0; i < nDep; i++) {
	k2[i] *= dDelta;
	xtemp[i] = xin[i] + 0.5*k2[i];
	}
    deriv(CTX, t + 0.5*dDelta, xtemp, k3);
    for(i = 0; i < nDep; i++) {
	k3[i] *= dDelta;
	xtemp[i] = xin[i] + k3[i];
	}
    deriv(CTX, t + dDelta, xtemp, k4);
    for(i = 0; i < nDep; i++) {
	k4[i] *= dDelta;
	xout[i] = xin[i] + onesixth*(k1[i] + k4[i]) + onethird*(k2[i] + k3[i]);
	}
    }

void
RungeKutta(void *CTX, 
	   void (*deriv)(void *, double, double *, double*),
	   int nDep,		/* number of dependent variables */
	   double tin,		/* independent variable */
	   double *xin,		/* array of input */
	   double tout,
	   double *xout,		/* array of output */
	   int nSteps)
{
    double dDelta = (tout - tin)/nSteps;
    int i;
    double xtemp[MAXDEP];
	
    assert(nDep <= MAXDEP);

    for(i = 0; i < nDep; i++) {
	xtemp[i] = xin[i];
	}
    
    for(i = 0; i < nSteps; i++) {
	RungeStep(CTX, deriv, nDep, tin, xtemp, xtemp, dDelta);
	tin += dDelta;
	}
    for(i = 0; i < nDep; i++) {
	xout[i] = xtemp[i];
	}
    }

#if 0
/*
 * Harmonic oscillator test
 */

void csmGrowthFacDeriv(void* csm, double dlnExp, double *dlnDelta, double
		       *dlnDeltadot)
{
    dlnDeltadot[0] = dlnDelta[1];
    dlnDeltadot[1] = -dlnDelta[0];
    }

main()
{
    const int nSteps = 25;
    double dDStart[2];
    double dDEnd[2];
    double dDelta = .8;
    int i;
    void *csm;
    
    dDStart[0] = 0.0;
    dDStart[1] = 1.0;
    for(i = 0; i < nSteps; i++){
	double dlnExp = (i+1)*dDelta;
	
	RungeKutta(csm, (void (*)(void *, double, double *, double*))
	       csmGrowthFacDeriv, 2, 0.0, dDStart, dlnExp, dDEnd,
	       nSteps);
	printf("%g %g %g\n", dlnExp, dDEnd[0], dDEnd[1]);
	}
    }
#endif

