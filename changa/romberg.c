/*
 * A open interval Romberg integrator
 * Originally written for the cosmo module in  PKDGRAV by Tom Quinn
 */
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <assert.h>
#include "romberg.h"

#define MAXLEV 13

/*
 ** Romberg integrator for an open interval.
 */

double
dRombergO(void *CTX,double (*func)(void *, double),double a,double b,
		  double eps)
{
    double tllnew;
    double tll;
    double tlk[MAXLEV+1];
    int n = 1;
    int nsamples = 1;

    tlk[0] = tllnew = (b-a)*(*func)(CTX, 0.5*(b+a));
    if(a == b) return tllnew;
    
    tll = DBL_MAX;

    while((fabs((tllnew-tll)/tllnew) > eps) && (n < MAXLEV)) {
		/*
		 * midpoint rule.
		 */
		double deltax;
		double tlktmp;
		int i;

		nsamples *= 3;
		deltax = (b-a)/nsamples;
		tlktmp = tlk[0];
		tlk[0] = tlk[0]/3.0;
	
		for(i=0;i<nsamples/3;i++) {
			tlk[0] += deltax*(*func)(CTX,a + (3*i + 0.5)*deltax);
			tlk[0] += deltax*(*func)(CTX,a + (3*i + 2.5)*deltax);
			}
    
		/*
		 * Romberg extrapolation.
		 */

		for(i=0;i<n;i++) {
			double tlknew = (pow(9.0, i+1.)*tlk[i] - tlktmp)
				/(pow(9.0, i+1.) - 1.0);
	    
			if(i+1 < n)
			    tlktmp = tlk[i+1];
			tlk[i+1] = tlknew;
			}
		tll = tllnew;
		tllnew = tlk[n];
		n++;
		}

    assert((fabs((tllnew-tll)/(tllnew+eps)) < eps));

    return tllnew;
    }

double
dRombergC(void *CTX,double (*func)(void *, double),double a,double b,
		  double eps)
{
    double tllnew;
    double tll;
    double tlk[MAXLEV+1];
    int n = 1;
    int nsamples = 1;

    tlk[0] = tllnew = (b-a)*0.5*((*func)(CTX, b) + (*func)(CTX, a));
    if(a == b) return tllnew;
    if(tllnew == 0.0) return tllnew;
    
    tll = DBL_MAX;

    while((fabs((tllnew-tll)/tllnew) > eps) && (n < MAXLEV)) {
		/*
		 * midpoint rule.
		 */
		double deltax;
		double tlktmp;
		int i;

		nsamples *= 2;
		deltax = (b-a)/nsamples;
		tlktmp = tlk[0];
		tlk[0] = tlk[0]/2.0;
	
		for(i=0;i<nsamples/2;i++) {
			tlk[0] += deltax*(*func)(CTX,a + (2*i + 1.0)*deltax);
			}
    
		/*
		 * Romberg extrapolation.
		 */

		for(i=0;i<n;i++) {
			double tlknew = (pow(4.0, i+1.)*tlk[i] - tlktmp)
				/(pow(4.0, i+1.) - 1.0);
	    
			tlktmp = tlk[i+1];
			tlk[i+1] = tlknew;
			}
		tll = tllnew;
		tllnew = tlk[n];
		n++;
		}

    assert((fabs((tllnew-tll)/tllnew) < eps));

    return tllnew;
    }
