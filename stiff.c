#ifndef NOCOOLING
/*
 * Integrate "stiff" equations similar to chemical equilibrium
 * equations.
 *
 * The source code below is taken from Mott, D.R. & Oran, E.S., 2001,
 * "CHEMEQ2: A Solver for the Stiff Ordinary Differential Equations of
 * Chemical Kinetics", Naval Research Laboratory,
 * NRL/MR/6400-01-8553.  The report documentation page states under
 * distribution/availability statement states:
 * "Approved for public release; distribution is unlimited."
 */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "stiff.h"

static inline double max(double x, double y) 
{
    if(x > y) return x;
    else return y;
    }

static inline double min(double a, double b)
{
    if(a < b) return a;
    else return b;
    }

/* implement fortran sign function: return a with the sign of b */
static inline double sign(double a, double b) 
{
    double aabs = fabs(a);
    if(b >= 0.0) return aabs;
    else return -aabs;
    }

/*
 * Set integration parameters and allocate scratch arrays.
 */
STIFF *StiffInit( double eps, int nv, void *Data,
		  void (*derivs)(double, const double *, double *, double *,
				void *Data)
		  ) 
{
  STIFF *s;
  int i;

  s = (STIFF *) malloc(sizeof(STIFF));
  assert(s!=NULL);

  s->nv = nv;
  s->epsmin = eps;
  s->sqreps = 5.0*sqrt(eps);
  s->epscl = 1.0/eps;
  s->epsmax = 10.0;
  s->dtmin = 1e-15;
  s->itermax = 3; /*Increased from 1 to 3 to speed integration by
		    calculating more correctors.  Adjustable parameter*/ 
  s->ymin = malloc(nv*sizeof(*(s->ymin)));
  for(i = 0; i < nv; i++)
      s->ymin[i] = 1e-300;
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

/**
 * Specify minimum values of quantities.
 */
void StiffSetYMin(STIFF *s, const double *ymin)
{
    int i;

    for(i = 0; i < s->nv; i++)
       s->ymin[i] = ymin[i];
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
    /*
     * Local copies of Stiff context
     */
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
    int gcount = 0;		/* count calls to derivs */
    int rcount = 0;		/* count restart steps */
    double scrtch;
    double ascr;
    double scr1;
    double scr2;
    double dt;			/* timestep used by the integrator */
    double ts;			/* t at start of the chemical timestep */
    double alpha;		/* solution parameter used in update */
    int iter;			/* counter for corrector iterations */
    double eps;			/* maximum correction term */
    double rtaub;
    double qt;			/* alpha weighted average of q */
    double pb;
    const double tfd = 1.000008; /* fudge for completion of timestep */
    double rteps;		/* estimate of sqrt(eps) */
    double dto;			/* old timestep; to rescale rtaus */
    double temp;                
    
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
	temp = -fabs(ascr-d[i])*scr2;
	/* If the species is already at the minimum, disregard
	   destruction when calculating step size */
	if (y[i] == ymin[i]) temp = 0.0;
	scrtch = max(scr1,max(temp,scrtch));
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
	    /*	    if(iter == 1) {  Removed from original algorithm
		    so that previous, rather than first, corrector is
		    compared to.  Results in faster integration. */
		/*
		c the first corrector step advances the time (tentatively) and
		c saves the initial predictor value as y1 for the timestep
		check later.
		*/
		tn = ts + dt;
		for(i = 0; i < n; i++)
		    y1[i] = y[i];
		/*		} Close for "if(iter == 1)" above */
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
	c calculate new f, check for convergence, and limit decreasing
	c functions. the order of the operations in this loop is important.
	*/
	for(i = 0; i < n; i++) {
	    scr2 = max(ys[i] + dt*scrarray[i], 0.0);
	    scr1 = fabs(scr2 - y1[i]);
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
	    printf("y[0]: %e, y[1]: %e, y[2]: %e, y[3]: %e, y[4]: %e",y[0],y[1],y[2],y[3],y[4]);
	    printf("q[0]: %e, q[1]: %e, q[2]: %e, q[3]: %e, q[4]: %e",q[0],q[1],q[2],q[3],q[4]);
	    printf("d[0]: %e, d[1]: %e, d[2]: %e, d[3]: %e, d[4]: %e",d[0],d[1],d[2],d[3],d[4]);
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
	    /*    & .or. stab. gt. 1 */
	    rcount++;
	    /*
	    c After an unsuccessful step the initial timescales don't
	    c change, but dt does, requiring rtaus to be scaled by the
	    c ratio of the new and old timesteps.
	    */
	    dto = dt/dto;
	    for(i = 0; i < n; i++) {
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
	s->derivs(tn + tstart, y, q, d, s->Data);
	gcount++;
	}
    }

#ifdef TESTCHEMEQ
/*
 * Test case for the CHEMEQ2 solver.
 */

void csdfe(double t, double *y, double *q, double *d, void *data);

#include <time.h>

int
main(int argc, char** argv) 
{
    /*    
C
C This is the driver program for the seven-species cesium
C mechanism test problem. The code integrates the system
C MXCASE times using differnt values of the chemeq2 variable
C epsmin (set by passing an entry from array EPS through
C CHEMSP before each integration).

C 
C For this example,the external subroutine that calculates the
C source terms is called CSDFE.
C
    */
    char *spsym[7] = {"02-", "CS+", "CS", "CS02", "02", "N2", "NE"} ;
    double eps[15] = {.1, .05, .01, .005, .001, .0005, .0001, .00005,
		      .00001, .000005, .000001, 5.e-7, 1.e-7, 5.e-8, 1.e-8};
/* ,5.e-9, 1.e-9,5.e-10,1.e-10 */
    double tscale;
    
    int mxcase = 9;
    int inlp, ns, na;
    double ti, tf, deltat;
    double yi[7], yf[7];
    double epsil[15];
    int icase;

    /*
Note that the timing routines included maynot work on
all systems. Extra timing options are included as comments.
    */
    double dtime, delta, tarray[2];
    double tnorm;
    delta = 1.;
    /*
    C
    C INITIALIZE CONTROL PARAMETERS.
    C
    C TSCALE is simply a normalization factor for the timing
    C results. It can be used to compare results from differnt
    C machines (by setting it to the time required for that
    C machine to solve a standard problem of some sort) or to
    C simply make the timing results more "friendly."
    */
    tscale = 1.0/1024.;
    
    /*
C INLP allows the user to subdivide the interval over which
C each test is run. For INLP=1,CHEMEQ2 is sent the full
C interval TF-TI (specified below) as the global timestep.
    */

    inlp = 1;
    
    /*
C For this particular test, the electron number density is not
C integrated. The other five reacting species are integrated,
C and the electron density is found through charge conservation.
C This calculation is done within CSDFE. Therefore, NA=5 is
C the number of equations that are integrated, but NS=7 is
C number of species. Species to be integrated must be placed in
C first NA positions within the Y array. CHEMEQ2 only works with
C these first NA entries since NA is passed in the argument list
C below, but all NS values are available to and used by CSDFE.
    */
    ns = 7;
    na = 5;
    /* "TI" - INITIAL TIME, "TF" - FINAL TIME. */
    ti = 0.0;
    tf = 1000.0;
    deltat = (tf - ti)/inlp;
     
     /*
C STORE INITIAL(TI = 0.0) AND FINALT(F = 1000.0)
C 02-
    */
    yi[0] = 5.200e+02;
    yf[0] = 2.59139492061e+04;
    /* CS+ */
    yi[1] = 6.200E+02;
    yf[1] = 7.55718460300e+04;
    /* Cs */
    yi[2] = 1.000E+12;
    yf[2] = 1.53194051722e+03;
    /* CS02 */
    yi[3] = 0.;
    yf[3] = 9.99999923516e+11;
    /* 02 */
    yi[4] = 3.600E+14;
    yf[4] = 3.59000000051e+14;
    /* N2 */
    yi[5] = 1.400E+15;
    yf[5] = 1.40000000000e+15;
    /* NE */
    yi[6] = 1.000E+02;
    yf[6] = 4.96578968239e+04;
    /*
      C LOOP OVER THE TEST CASES.
    */
    for(icase = 0; icase < mxcase; icase++) {
	STIFF *s;
	double cput;
	double t = 0.0;
	double y[7];
	double sum;
	int istep;
	int i;
	
	printf("case %d: eps: %g inlp: %d\n", icase, eps[icase], inlp);
	s = StiffInit(eps[icase], 5, NULL, csdfe);
	cput = clock()/CLOCKS_PER_SEC;
	 /* RESET "Y" TO INITIAL VALUES "YI"'. */
	for(i = 0; i < ns; i++)
	    y[i] = yi[i];
	/*
	  C SET TIMER.
	  T = SECNDS(O.0)
	  delta = dtime(tarray);
	*/
	/*
	 INNER LOOP TO DETERMINE OVERHEAD OR RELATIVE STARTING EFFECIENCY
	 OF ITEGRATION SCHEME BEING TESTED.
	*/
	for(istep = 0; istep < inlp; istep++) {
	    /* CALL INTEGRATOR. */
	    /* chemeq2(deltat, csdfe, na, y); */
	    StiffStep(s, y, t, deltat);
	    t += deltat;
	    }
	/* CALCULATE CPU TIME USED IN THE INTEGRATION PROCESS.
	   delta = dtime(tarray);
	   dsec = tarray[0];
	 */
	/* DSEC = delta */
	cput = clock()/CLOCKS_PER_SEC;
	tnorm = (int)(cput/tscale + .5);
	 /*
	  Calculate final electron density from densities of other
	  charges species
	 */
	y[6] = y[1] - y[0];

	/* CALCULATE RELATIVE ERROR. */
	for(i = 0; i < ns; i++)
	    epsil[i] = fabs(y[i] - yf[i])/min(y[i] , yf[i]);

	sum = 0.0;
	for(i = 0; i < ns; i++)
	    sum = sum + pow(epsil[i],2);

	/*
	  Root-mean-square error is calculated using ns-1 (rather than ns)
	  since N2 is inert.
	*/
	sum = sqrt(sum/(ns-1));

	/*
	  PRINT RESULTS.
	*/
	printf("ti: %g, tf: %g\n", ti, tf);
	for(i = 0; i < ns; i++)
	    printf("%s %g %g %g %g\n", spsym[i], yi[i], yf[i], y[i],
		   epsil[i]);
	printf("sum: %g\n", sum);
	printf("cpu %g tnorm %g\n", cput, tnorm);
	printf("eps %g, cput %g, tnorm %g, sum %g\n", eps[icase], cput,
	       tnorm, sum);

	/* CALL CHEMCT(TF) */
	StiffFinalize(s);
	}
    return 0;
    }

void csdfe(double t, double *y, double *q, double *d, void *data)
{
    /*
cd * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
cd 
cd csdfe(y, q, d, t)
cd 
cd description:
cd derivative function evaluator(gsub) for an atmospheric chemical
cd relaxation test problem involving cesium and cesium ions. format-
cd ion and loss rates are calculated for this set of "stiff ordinary
cd differential equations" that was suggested by by d. edelson of
cd bell laboratories.
cd 
cd argument list definitions:
cd y(i)		r	current values of the functions plus the	i/o
cd 		      	extra data at the end of the array that may be
cd			passed back and forth between "csdfe" and the
cd			main program. locations in y(i) which represent
cd			the functions being advanced should not be
cd			tampered with here.
cd q(i)		r	total formation rates.
cd d(i)		r	total loss rates.
cd t		r	the value of the independent variable.
cd
cd * * * * * * * * * * *
c
c local specifications.
c _ _____ ____ ___
    */
    double ne, n2;
    double o2m, csp, cs, cso2, o2;
    double cr1, cr2, cr3, cr4, cr5, cr6, cr7;
    /*
     * utilize local storage for varibles.
     */
    o2m = y[0];
    csp = y[1];
    cs =y[2];
    cso2 = y[3];
    o2 = y[4];
    n2 = y[5];
    /*
c write(63,*) t
c
c calculate electron density for local use and transmission back to
c the main program via y(7). however in this case this value should
c not be trusted since "chemeq"will not call the gsub" with the
c latest function values after the final step has converged. y(7 )
c will be one iteration behind in this case. y(7) and y(6) are
c examples tho, of how data may be transfered between the gsub" and
c the main program.
    */
    ne = max(csp - o2m, 0.0);
    y[6] = ne;
    /*
c calculate reaction rates.
    */
    cr1 = 5.00e-08*o2m*csp;
    cr2 = 1.00e-12*csp*ne;
    cr3 = 3.24e-03*cs;
    cr4 = 4.00e-01*o2m;
    cr5 = 1.00e-31*o2*cs*(cs + cso2 + n2 + o2);
    cr6 = 1.24e-30*o2*o2*ne;
    cr7 = 1.00e-31*o2*n2*ne;
    /*
if(t.ge.700.) then
c cr4= 0.
c cr6 = 0.
c cr7 = 0.
c end if

c calculate total formation rates (c(i)) and total loss rates (d(i))
c for each species.
c 
c o2m
    */
    q[0] = cr6 + cr7;
    d[0] = cr1 + cr4;
    /*
      Cs+
    */
    q[1] = cr3;
    d[1] = cr1 + cr2;
    /*
      Cs
    */
    q[2] = cr1 + cr2;
    d[2] = cr3 + cr5;
    /*
      cso2
    */
    q[3] = cr5;
    d[3] = 0.0;
    /*
q(4) = q(4) -
d(4) = - l.OOe-31*o2*cs*cso2
    */
    /*
      o2
    */
    q[4] = cr1 + cr4;
    d[4]= cr5 + cr6 + cr7;
    return;
    }

#else
/* TESTCHEMEQ excludes root finder */

/*
 * The following is an implementation of the Brent root finding
 * algorithm.  This is based on code from the GSL which has the
 * following copyright.  The code was modified to avoid the overhead
 * that is somewhat unnecessary for such a simple routine.
 */
/* roots/brent.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Reid Priedhorsky, Brian Gou\
gh
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, US\
A.
 */

const int itmaxRoot = 100;
const double epsRoot = 3.0e-8;

double RootFind(double (*func)(void *Data, double), void *Data, double x1,
                double x2, double tol)
{
    int i;
    double a = x1;
    double fa = (*func)(Data, a);
    double b = x2;
    double fb = (*func)(Data, b);
    double c = x2;
    double fc = fb;
    double d = x2 - x1;
    double e = x2 - x1;

    if ((fa < 0.0 && fb < 0.0) || (fa > 0.0 && fb > 0.0))
    {
        fprintf(stderr, "RootFind: endpoints do not straddle y=0");
        assert(0);
        }

    for(i = 0; i < itmaxRoot; i++) {
        double m;
        double dTol;

        int ac_equal = 0;

        if ((fb < 0 && fc < 0) || (fb > 0 && fc > 0)) {
            ac_equal = 1;
            c = a;
            fc = fa;
            d = b - a;
            e = b - a;
            }

        if (fabs (fc) < fabs (fb)) {
            ac_equal = 1;
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
            }

        dTol = 0.5 * epsRoot * fabs (b) + 0.5*tol;
        m = 0.5 * (c - b);

        if (fb == 0) {
            return b;   /* SUCCESS */
            }
        if (fabs (m) <= dTol) {
            return b; /* SUCCESS */
            }

        if (fabs (e) < dTol || fabs (fa) <= fabs (fb)) {
            d = m;            /* use bisection */
            e = m;
            }
        else {
            double p, q, r;   /* use inverse cubic interpolation */
            double s = fb / fa;

            if (ac_equal) {
                p = 2 * m * s;
                q = 1 - s;
                }
            else {
                q = fa / fc;
                r = fb / fc;
                p = s * (2 * m * q * (q - r) - (b - a) * (r - 1));
                q = (q - 1) * (r - 1) * (s - 1);
                }
            if (p > 0) {
                q = -q;
                }
            else {
                p = -p;
                }

            if (2 * p < min(3 * m * q - fabs (dTol * q), fabs (e * q))) {
                e = d;
                d = p / q;
                }
            else {
                /* interpolation failed, fall back to bisection */
                d = m;
                e = m;
                }
            }
        a = b;
        fa = fb;

        if (fabs (d) > dTol) {
            b += d;
            }
        else {
            b += (m > 0 ? +dTol : -dTol);
            }
        fb = (*func)(Data, b);
        }

    fprintf(stderr, "brent: number of interations exceeded");
    assert(0);
    return 0.0;
    }

#endif /* STIFFTEST */

#endif

