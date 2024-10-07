#ifndef STIFF_HINCLUDED
#define STIFF_HINCLUDED

#ifdef CUDA
#include <cuda_runtime.h>

#define CUDA_DH __device__ __host__
#else
#define CUDA_DH
#endif

#define EPSMAX 10.0
#define DTMIN 1e-15
#define ITERMAX 3.0
#define YMIN0 1e-300

/** @brief Context for stiff integration
*/
typedef struct StiffContextStructure {
    double epsmin;		/* relative accuracy criterion */
    double sqreps;		/* parameter to calculate initial timestep */
    double epscl;		/* 1/epsmin */
    double epsmax;		/* repeat timestep if correction is
				   greater than epsmax*epsmin*y(i) */   
    double dtmin;		/* minimum timestep allowed */
    int itermax;		/* number of corrector iterations */
    int nv; 			/* number of dependent variables */
    double *ymin;		/* minimum value for y */
    double *y0;			/* initial y for global timestep */
    double *y1;			/* predicted value */
    double *q;			/* scratch for creation rate */
    double *d;			/* scratch for destruction rate */
    double *rtau;		/* ratio of timestep to timescale */
    double *ys;			/* initial y for individual timesteps */
    double *qs;			/* initial production rate */
    double *rtaus;		/* initial rtau */
    double *scrarray;
    void *Data; /* Pointer to fixed data used by derivs */
    void (*derivs)(double, const double [], double [], double [], void *Data);
} STIFF;

/*
 * Integrator Step Headers
 */

void StiffInit(STIFF *s,
         double eps, 	/* relative accuracy parameter */
		 int nv,	/* number of dependent variables */
		 void *Data, 	/* pointer to extra data */
		 void (*derivs)(double t, const double yin[],  /* input */
				double yheat[],	/* heating or creation
						   rate */
			       double ycool[],	/* cooling or
						   destruction rate */
				void *Data)
		 );

void StiffStep(STIFF *s, double y[], double tstart, double htry) ;
void StiffSetYMin(STIFF *s, const double *ymin);

#ifndef TESTCHEMEQ
/* 
 * Root Finder Header
 */

CUDA_DH double RootFind(double (*func)(void *Data, double), void *Data,
		double x1, double x2, double tol);

#endif
#endif

