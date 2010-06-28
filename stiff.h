#ifndef STIFF_HINCLUDED
#define STIFF_HINCLUDED

#include "gsl/gsl_odeiv.h"

typedef struct StiffContextStructure {
  double eps;
  gsl_odeiv_step *stepper;
  gsl_odeiv_control *controler;
  gsl_odeiv_evolve *evolver;
  gsl_odeiv_system system;
} STIFF;

/*
 * Integrator Step Headers
 */

STIFF *StiffInit( double eps, int nv, void *Data, 
		  int (*derivs)(double, double [], double [], void *Data),
		  int (*jacobn)(double x, double y[], double dfdx[],
				 double *dfdy, void *Data) );
		   
void StiffFinalize( STIFF *s );
void StiffStep(STIFF *s, double y[], double dydx[], double *xx, double htry, 
		double yscal[], double *hdid, double *hnext );

/* 
 * Root Finder Header
 */

double RootFind(double (*func)(void *Data, double), void *Data, double x1, double x2, double tol);

/*
 * Utils
 */

static double SQR(double a) 
{
    return a*a;
}

static double FMAX(double maxarg1, double maxarg2)
{
    return (maxarg1 > maxarg2 ? maxarg1 : maxarg2);
}

static double FMIN(double minarg1, double minarg2)
{
    return (minarg1 < minarg2 ? minarg1 : minarg2);
}

#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))

#endif

