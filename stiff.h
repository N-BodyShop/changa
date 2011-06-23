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
		  int (*derivs)(double, const double [], double [], void *Data),
		  int (*jacobn)(double x, const double y[], double dfdx[],
				 double *dfdy, void *Data) );
		   
void StiffFinalize( STIFF *s );
void StiffStep(STIFF *s, double y[], double dydx[], double *xx, double htry, 
		double yscal[], double *hdid, double *hnext );

/* 
 * Root Finder Header
 */
#include <gsl/gsl_roots.h>

typedef struct RootFindContextStructure {
  gsl_root_fsolver *solver;
} ROOTFIND;

ROOTFIND *RootFindInit();
void RootFindFinalize(ROOTFIND *r);

double RootFind(ROOTFIND *r, double (*func)(double, void *Data), void *Data,
		double x1, double x2, double tol);

#endif

