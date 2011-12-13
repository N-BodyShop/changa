#ifndef NOCOOLING
/*
 * Integrate "stiff" equations similar to chemical equilibrium
 * equations.
 */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "stiff.h"

inline double max(double x, double y) 
{
    if(x > y) return x;
    else return y;
    }

/* implement fortran sign function: return a with the sign of b */
inline double sign(double a, double b) 
{
    double aabs = fabs(a);
    if(b >= 0.0) return aabs;
    else return -aabs;
    }

STIFF *StiffInit( double eps, int nv, void *Data,
		  void (*derivs)(double, const double *, double *, double *,
				void *Data)
		  ) 
{
  STIFF *s;

  s = (STIFF *) malloc(sizeof(STIFF));
  assert(s!=NULL);

  s->nv = nv;
  s->epsmin = eps;
  s->sqreps = 5.0*sqrt(eps);
  s->epscl = 1.0/eps;
  s->epsmax = 10.0;
  s->dtmin = 1e-15;
  s->itermax = 1;
  s->ymin = calloc(nv, sizeof(*(s->ymin)));
  s->y0 = malloc(nv*sizeof(*(s->y0)));
  s->y1 = malloc(nv*sizeof(*(s->y1)));
  s->q = malloc(nv*sizeof(*(s->q)));
  s->d = malloc(nv*sizeof(*(s->d)));
  s->rtau = malloc(nv*sizeof(*(s->rtau)));
  s->ys = malloc(nv*sizeof(*(s->ys)));
  s->qs = malloc(nv*sizeof(*(s->qs)));
  s->rtaus = malloc(nv*sizeof(*(s->rtaus)));
  s->scrarray = malloc(nv*sizeof(*(s->scrarray)));

  s->Data = Data;
  s->derivs = derivs;

  s->epsmax = 10.0;

  return s;
}

void StiffFinalize( STIFF *s ) 
{
    free(s->ymin);
    free(s->y0);
    free(s->y1);
    free(s->q);
    free(s->d);
    free(s->rtau);
    free(s->ys);
    free(s->qs);
    free(s->rtaus);
    free(s->scrarray);
    free(s);
}

void StiffStep(STIFF *s,
	       double y[],	/* dependent variables */
	       double tstart, 	/* start time */
	       double dtg) 	/* time step */
{
/*
C
cd* * *** ** ********** * * * *** **** * ** ** * *
cd
cd chemeq2(dtg, gsub, n, y)
cd
cd original chemeq development:
cd originators: t.r. young nrl 1982
cd vax version: t.r. young nrl code 4040 may 1983
cd workstation: g. patnaik berkeley research jun 1995
cd
cd chemeq2 development: d.r. mott nrl code 6404 may 1999

   Conversion to C: T. Quinn, UW, Dec, 2011
cd
cd Description: Subroutine chemeq2 solves a class of "stiff" ODEs
cd associated with reactive flow problems that cannot be readily
cd solved by the standard classical methods. In contrast to the
cd original chemeq subroutine, this version uses the same
cd quasi-steady-state update for every species regardless of the
cd timescale for that species. An adaptive stepsize is chosen to
cd give accurate results for the fastest changing quantity, and a
cd stability check on the timestep is also available when the
cd corrector is iterated.
cd
cd NOTE: The accuracy-based timestep calculation can be augmented
cd with a stability-based check when at least three corrector
cd iterations are performed. To include this check, "uncomment"
cd the lines that start with "D", or use the compiler flag -d-lines"
cd if available to compile the code including these lines.  If the
cd lines are manually uncommented, the continuation characters
cd must be placed in the correct column. For most problems, the
cd stability check is not needed, and eliminating the calculations
cd and logic associated with the check enhances performance.
cd 
cd The routine assumes that all of the integrated quantites and the
cd time step are positive.
cd
cd argument list definition (name, type, description, input vs. output):
cd dtg  real    the interval of integration or the i
cd              range of the independent variable.
cd              0.0 <= t <= dtg. (global timestep)
cd gsub real    the name of the derivitive function i
cd              evaluator subroutine.
cd n    integer the number of equations to be i
cd              integrated. an error exisis if n is
cd              greater than nd set by the parameter
cd              statement.
cd y(n) real    the initial values at call time i/o
cd              and the final values at return time.
cd
cd Language and limitations: This subroutine is written in standard
cd FORTRAN 77. For high accuracy, this routine should be compiled
cd using whatever "double precision" flag is appropriate for the
cd platform being used (such as "f77 -r8 . . . .)
cd
cd subroutines referenced:
cd
cd gsub;    whose actual name and definition are supplied by the user
cd 	    is called to obtain the derivitive functions.
cd
cd call gsub(y, q, d, t)
cd argument list to gsub;
cd y(n) real current values of the dependent  i
cd           variable.
cd q(n) real calculated formation rates.      o
cd d(n) real calculated loss rates.           o
cd t    real current value of the independent i 
cd           variable. 
cd
*/

    double tn;			/* time within step */
    int i;
    int n = s->nv;
    double *y0 = s->y0;
    double *ymin = s->ymin;
    double *q = s->q;
    double *d = s->d;
    double *rtau = s->rtau;
    double *ys = s->ys;
    double *qs = s->qs;
    double *rtaus = s->rtaus;
    double *scrarray = s->scrarray;
    double *y1 = s->y1;
    double epsmin = s->epsmin;
    double sqreps = s->sqreps;
    double epscl = s->epscl;
    double epsmax = s->epsmax;
    double dtmin = s->dtmin;
    int itermax = s->itermax;
    int gcount = 0;
    int rcount = 0;		/* restart steps */
    double scrtch;
    double ascr;
    double scr1;
    double scr2;
    double dt;			/* timestep used by the integrator */
    double ts;			/* t at start of the chemical timestep */
    double alpha;
    int iter;
    double eps;			/* maximum correction term */
    double rtaub;
    double qt;			/* alpha weighted average of q */
    double pb;
    const double tfd = 1.000008; /* fudge for completion of timestep */
    double rteps;		/* estimate of sqrt(eps) */
    double dto;			/* old timestep; to rescale rtaus */
    
    tn = 0.0;
    for(i = 0; i < n; i++) {
	q[i] = 0.0;
	d[i] = 0.0;
	y0[i] = y[i];
	y[i] = max(y[i], ymin[i]);
	}
    
    s->derivs(tn + tstart, y, q, d, s->Data);
    gcount++;
    
    /*
    C
    c estimate the initial stepsize.
    C
    c strongly increasing functions(q >>> d assumed here) use a step-
    c size estimate proportional to the step needed for the function to
    c reach equilibrium where as functions decreasing or in equilibrium
    c use a stepsize estimate directly proportional to the character-
    c istic stepsize of the function. convergence of the integration
    c scheme is likely since the smallest estimate is chosen for the
    c initial stepsize.
    */
    scrtch  = 1.0e-25;
    for(i = 0; i < n; i++) {
	ascr = fabs(q[i]);
	scr2 = sign(1./y[i],.1*epsmin*ascr - d[i]);
	scr1 = scr2 * d[i];
	scrtch = max(scr1,max(-fabs(ascr-d[i])*scr2,scrtch));
	}
    dt = min(sqreps/scrtch,dtg);
    while(1) {
	/*
	  c the starting values are stored.
	*/
	ts = tn;
	for(i = 0; i < n; i++) {
	    rtau[i] = dt*d[i]/y[i];
	    ys[i] = y[i];
	    qs[i] = q[i];
	    rtaus[i] = rtau[i];
	    }

	/*
	 * find the predictor terms.
	 */
     restart:
	for(i = 0; i < n; i++) {
	    /*
	     * prediction
	     */
	    double rtaui = rtau[i];
	    /*
	    c note that one of two approximations for alpha is chosen:
	    c 1) Pade b for all rtaui (see supporting memo report)
	    c or
	    c 2) Pade a for rtaui<=rswitch,
	    c linear approximation for rtaui > rswitch
	    c (again, see supporting NRL memo report (Mott et al., 2000))
	    c
	    c Option 1): Pade b
	    */
	    alpha = (180.+rtaui*(60.+rtaui*(11.+rtaui)))
		/(360.+ rtaui*(60. + rtaui*(12. + rtaui)));
	    /*
	    c Option 2): Pade a or linear
	    c
	    c if(rtaui.le.rswitch) then
	    c      alpha = (840.+rtaui*(140.+rtaui*(20.+rtaui)))
	    c    &         / (1680. + 40. * rtaui*rtaui)
	    c else
	    c    alpha = 1.-1./rtaui
	    c end if
	    */
	    scrarray[i] = (q[i]-d[i])/(1.0 + alpha*rtaui);
	    }

	iter = 1;
	while(iter <= itermax) {
	    for(i = 0; i < n; i++) {
		/*
		C ym2(i) = ym1(i)
		C ym1(i) = y(i)
		*/
		y[i] = max(ys[i] + dt*scrarray[i], ymin[i]);
		}
	    if(iter == 1) {
		/*
		c the first corrector step advances the time (tentatively) and
		c saves the initial predictor value as y for the timestep
		check later.
		*/
		tn = ts + dt;
		for(i = 0; i < n; i++)
		    y1[i] = y[i];
		}
	    /*
	      evaluate the derivitives for the corrector.
	    */
	    s->derivs(tn + tstart, y, q, d, s->Data);
	    gcount++;
	    eps = 1.0e-10;
	    for(i = 0; i < n; i++) {
		rtaub = .5*(rtaus[i]+dt*d[i]/y[i]);
		/*
		c Same options for calculating alpha as in predictor:
		c
		c Option 1): Pade b
		*/
		alpha = (180.+rtaub*(60.+rtaub*(11.+rtaub)))
		    / (360. + rtaub*(60. + rtaub*(12. + rtaub)));
		/*
		c Option 2): Pade a or linear
		c
		c if(rtaub.le.rswitch)
		c then
		c alpha = (840.+rtaub*(140.+rtaub*(20.+rtaub)))
		c & / (1680. + 40.*rtaub*rtaub)
		c else
		c alpha = 1.- 1./rtaub
		c end if
		*/
		qt = qs[i]*(1. - alpha) + q[i]*alpha;
		pb = rtaub/dt;
		scrarray[i] = (qt - ys[i]*pb) / (1.0 + alpha*rtaub);
		}
	    iter++;
	    }
	/*
	c calculate newf, check for convergence, and limit decreasing
	c functions. the order of the operations in this loop is important.
	*/
	for(i = 0; i < n; i++) {
	    scr2 = max(ys[i] + dt*scrarray[i], 0.0);
	    scr1 = fabs(scr2- y[i]);
	    y[i] = max(scr2, ymin[i]);
	    /*
	    C ym2(i) = ymi(i)
	    C yml(i) = y(i)
	    */
	    if(.25*(ys[i] + y[i]) > ymin[i]) {
		scr1 = scr1/y[i];
		eps = max(.5*(scr1+
			      min(fabs(q[i]-d[i])/(q[i]+d[i]+1.0e-30),scr1)),eps);
		}
	    }
	eps = eps*epscl;
	/* 
	   print out dianostics if stepsize becomes too small.
	*/
	if(dt <= dtmin + 1.0e-16*tn) {
	    fprintf(stderr, "stiffchem: step size too small\n");
	    assert(0);
	    }
	/*
	c check for convergence.
	c
	c The following section is used for the stability check
	C       stab = 0.01
	C if(itermax.ge.3) then
	C       do i=1,n
	C           stab = max(stab, abs(y(i)-yml(i))/
	C       &       (abs(ymi(i)-ym2(i))+1.e-20*y(i)))
	C end do
	C endif
	*/
	if(eps <= epsmax) {
	    /*
	      & .and.stab.le.1.
	    c
	    c Valid step. Return if dtg has been reached.
	    */
	    if(dtg <= tn*tfd) return;
	    }
	else {
	    /*
	      Invalid step; reset tn to ts
	    */
	    tn = ts;
	    }
	/*
	  perform stepsize modifications.
	  estimate sqrt(eps) by newton iteration.
	*/
	rteps = 0.5*(eps + 1.0);
	rteps = 0.5*(rteps + eps/rteps);
	rteps = 0.5*(rteps + eps/rteps);

	dto = dt;
	dt = min(dt*(1.0/rteps+.005), tfd*(dtg - tn));
	/* & ,dto/(stab+.001) */
	/*
	  begin new step if previous step converged.
	*/
	if(eps > epsmax) {
	    /*    & .or. stab. gt. */
	    rcount++;
	    /*
	    c After an unsuccessful step the initial timescales don't
	    c change, but dt does, requiring rtaus to be scaled by the
	    c ratio of the new and old timesteps.
	    */
	    dto = dt/dto;
	    for(i = 1; i < n; i++) {
		rtaus[i] = rtaus[i]*dto;
		}
	    /*
	     * Unsuccessful steps return to line 101 so that the initial
	     * source terms do not get recalculated.
	    */
	    goto restart;
	    }
	/*
	  Successful step; get the source terms for the next step
	  and continue back at line 100
	*/
	gcount++;
	}
    }

#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

/*
 * Interface to GSL Root finder
 */
ROOTFIND *RootFindInit()
{
    ROOTFIND *r;
    const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
    
    r = (ROOTFIND *)malloc(sizeof(ROOTFIND));
    assert(r!=NULL);
    r->solver = gsl_root_fsolver_alloc(T);
    
    return r;
    }

void RootFindFinalize(ROOTFIND *r)
{
    gsl_root_fsolver_free(r->solver);
    free(r);
    }

double RootFind(ROOTFIND *r, double (*func)(double, void *Data), void *Data,
		double x1,  /* lower bracket on the root */
		double x2,  /* upper bracket on the root */
		double tol)  /* Absolute tolerance */
{
    const int itmax = 100;
    const double epsRel = 3.0e-8;
    int iter;
    int status;
    double root;
    
    gsl_function F;
    
    F.function = func;
    F.params = Data;
    
    gsl_root_fsolver_set (r->solver, &F, x1, x2);
    for(iter = 0; iter < itmax; iter++) {
	status = gsl_root_fsolver_iterate(r->solver);
	root = gsl_root_fsolver_root(r->solver);
	x1 = gsl_root_fsolver_x_lower(r->solver);
	x2 = gsl_root_fsolver_x_upper(r->solver);
	status = gsl_root_test_interval (x1, x2, tol, epsRel);
	if(status == GSL_SUCCESS)
	    break;
	}
    assert(iter < itmax);
    return root;
    }

#endif

