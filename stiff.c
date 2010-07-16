#ifndef NOCOOLING
/*
 * Interface to GSL initial value integration package using the Stiff
 * BS routine.
 */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "stiff.h"


STIFF *StiffInit( double eps, int nv, void *Data,
		  int (*derivs)(double, double [], double [], void *Data),
		  int (*jacobn)(double x, double y[], double dfdx[], double *dfdy, void *Data) 
		   ) 
{
  int k,iq;

  STIFF *s;
  const gsl_odeiv_step_type *T = gsl_odeiv_step_bsimp;

  s = (STIFF *) malloc(sizeof(STIFF));
  assert(s!=NULL);

  s->eps = eps;
  s->stepper = gsl_odeiv_step_alloc(T, nv);
  s->controler = gsl_odeiv_control_y_new(0.0, eps);
  s->evolver = gsl_odeiv_evolve_alloc(nv);
  gsl_odeiv_system system = {derivs, jacobn, nv, Data};
  s->system = system;

  return s;
}

void StiffFinalize( STIFF *s ) 
{
  gsl_odeiv_evolve_free (s->evolver);
  gsl_odeiv_control_free (s->controler);
  gsl_odeiv_step_free (s->stepper);
  free(s);
}

void StiffStep(STIFF *s, double y[], double dydx[], double *xx, double htry, 
		double yscal[], double *hdid, double *hnext )
{
  double xend = *xx + htry;
  int status;
  double h;

  h=htry;
  status = gsl_odeiv_evolve_apply(s->evolver, s->controler, s->stepper,
				      &s->system, xx, xend, &h, y);
  *hdid = h;
  *hnext = h;
}

#if 0
/*
 * Jacobn test code
 */ 

int jacobn(double x, /* independent variable */
	   double y[], /* dependent variable */
	   double *dfdy, /* row-ordered matrix returning the Jacobian */
	   double dfdx[], /* derivative of f */
	   void *params)  /* Opaque pointer to non-evolving parameters */
{
  int i;
  const int ndim = 3;

  for (i=0;i<ndim;i++) dfdx[i]=0.0;
  dfdy[0*ndim + 1] = -0.013-1000.0*y[2];
  dfdy[0*ndim + 2]=0.0;
  dfdy[0*ndim + 3] = -1000.0*y[0];
  dfdy[1*ndim + 1]=0.0;
  dfdy[1*ndim + 2] = -2500.0*y[2];
  dfdy[1*ndim + 3] = -2500.0*y[1];
  dfdy[2*ndim + 1] = -0.013-1000.0*y[2];
  dfdy[2*ndim + 2] = -2500.0*y[2];
  dfdy[2*ndim + 3] = -1000.0*y[0]-2500.0*y[1];
}

void derivs(double x, double y[], double dydx[], void *params)
{
  dydx[0] = -0.013*y[0]-1000.0*y[0]*y[2];
  dydx[1] = -2500.0*y[1]*y[2];
  dydx[2] = -0.013*y[0]-1000.0*y[0]*y[2]-2500.0*y[1]*y[2];
}
#endif

#endif

