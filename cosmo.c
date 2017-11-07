/*
 * Cosmology routines originally written for PKDGRAV by Thomas Quinn
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "runge.h"
#include "romberg.h"

#ifdef CRAY_T3D
#include "hyperlib.h"
#endif

#include "cosmo.h"

/*
 * Cosmological module for PKDGRAV.
 * N.B.  This code is being shared with skid and the I.C. generator.
 */

void csmInitialize(CSM *pcsm) 
{
    CSM csm;
    
    csm = (CSM) malloc(sizeof(struct csmContext));
    assert(csm != NULL);
    
    csm->dHubble0 = 0.0;
    csm->dOmega0 = 0.0;
    csm->dLambda = 0.0;
    csm->dOmegaRad = 0.0;
    csm->dOmegab = 0.0;
    csm->dQuintess = 0.0;
    csm->bComove = 0;
    
    *pcsm = csm;
    }

#define EPSCOSMO 1e-7

/*
 ** The cosmological equation of state is entirely determined here.  We
 ** will derive all other quantities from these two functions.
 */

double csmExp2Hub(const CSM csm, double dExp)
{
    double dOmegaCurve = 1.0 - csm->dOmega0 -
	csm->dLambda - csm->dOmegaRad - csm->dQuintess;
    
    assert(dExp > 0.0);
    return csm->dHubble0
	*sqrt(csm->dOmega0*dExp
	      + dOmegaCurve*dExp*dExp
	      + csm->dOmegaRad
	      + csm->dQuintess*dExp*dExp*sqrt(dExp)
	      + csm->dLambda*dExp*dExp*dExp*dExp)/(dExp*dExp);
    }

/*
 ** Return a double dot over a.
 */
double csmExpDot2(CSM csm, double dExp)
{
    return csm->dHubble0*csm->dHubble0
	*(csm->dLambda  - 0.5*csm->dOmega0/(dExp*dExp*dExp)
	  + 0.25*csm->dQuintess/(dExp*sqrt(dExp))
	  - csm->dOmegaRad/(dExp*dExp*dExp*dExp));
    }

double csmTime2Hub(CSM csm,double dTime)
{
	double a = csmTime2Exp(csm,dTime);

	assert(a > 0.0);
	return csmExp2Hub(csm, a);
	}

double csmCosmoTint(const CSM csm, double dY)
{
    double dExp = pow(dY, 2.0/3.0);
    
    assert(dExp > 0.0);
    return 2.0/(3.0*dY*csmExp2Hub(csm, dExp));
    }

double csmExp2Time(CSM csm,double dExp)
{
	double dOmega0 = csm->dOmega0;
	double dHubble0 = csm->dHubble0;
	double a0,A,B,eta;

	if (!csm->bComove) {
		/*
		 ** Invalid call!
		 */
		assert(0);
		}
	if(csm->dLambda == 0.0 && csm->dOmegaRad == 0.0 &&
	   csm->dQuintess == 0.0) {
	    if (dOmega0 == 1.0) {
		    assert(dHubble0 > 0.0);
		    if (dExp == 0.0) return(0.0);
		    return(2.0/(3.0*dHubble0)*pow(dExp,1.5));
		    }
	    else if (dOmega0 > 1.0) {
		    assert(dHubble0 >= 0.0);
		    if (dHubble0 == 0.0) {
			    B = 1.0/sqrt(dOmega0);
			    eta = acos(1.0-dExp); 
			    return(B*(eta-sin(eta)));
			    }
		    if (dExp == 0.0) return(0.0);
		    a0 = 1.0/dHubble0/sqrt(dOmega0-1.0);
		    A = 0.5*dOmega0/(dOmega0-1.0);
		    B = A*a0;
		    eta = acos(1.0-dExp/A);
		    return(B*(eta-sin(eta)));
		    }
	    else if (dOmega0 > 0.0) {
		    assert(dHubble0 > 0.0);
		    if (dExp == 0.0) return(0.0);
		    a0 = 1.0/dHubble0/sqrt(1.0-dOmega0);
		    A = 0.5*dOmega0/(1.0-dOmega0);
		    B = A*a0;
		    eta = acosh(dExp/A+1.0);
		    return(B*(sinh(eta)-eta));
		    }
	    else if (dOmega0 == 0.0) {
		    assert(dHubble0 > 0.0);
		    if (dExp == 0.0) return(0.0);
		    return(dExp/dHubble0);
		    }
	    else {
		    /*
		     ** Bad value.
		     */
		    assert(0);
		    return(0.0);
		}	
	    }
	else {
             /* Set accuracy to 0.01 EPSCOSMO to make Romberg integration
              * more accurate than Newton's method criterion in csmTime2Exp. --JPG
              */
            return dRombergO(csm, (double (*)(const void *, double)) csmCosmoTint,
                             0.0, pow(dExp, 1.5), 0.01*EPSCOSMO);
 	    }
	}

double csmTime2Exp(CSM csm,double dTime)
{
	double dHubble0 = csm->dHubble0;

	if (!csm->bComove) return(1.0);
	else {
	    double dExpOld = 0.0;
	    double dExpNew = dTime*dHubble0;
	    double dDeltaOld = dExpNew;	/* old change in interval */
	    double dUpper = 1.0e38; /* bounds on root */
	    double dLower = 0.0;
	    
	    int it = 0;
	    
	    /*
	     * Root find with Newton's method.
	     */
	    do {
		double dExpNext;
	    	double f = dTime - csmExp2Time(csm, dExpNew);
		double fprime = 1.0/(dExpNew*csmExp2Hub(csm, dExpNew));
		
		if(f*fprime > 0)
		    dLower = dExpNew;
		else 
		    dUpper = dExpNew;
		
		dExpOld = dExpNew;
		dDeltaOld = f/fprime;
		dExpNext = dExpNew + dDeltaOld;
		/* check if bracketed */
		if((dExpNext > dLower) && (dExpNext < dUpper))
		    dExpNew = dExpNext;
		else
		    dExpNew = 0.5*(dUpper + dLower);
		it++;
		assert(it < 40);
		}
	    while (fabs(dExpNew - dExpOld)/dExpNew > EPSCOSMO);
	    return dExpNew;
	    }
	}

double csmComoveDriftInt(const CSM csm, double dIExp)
{
    return -dIExp/(csmExp2Hub(csm, 1.0/dIExp));
    }

/*
 ** Make the substitution y = 1/a to integrate da/(a^2*H(a))
 */
double csmComoveKickInt(const CSM csm, double dIExp)
{
    return -1.0/(csmExp2Hub(csm, 1.0/dIExp));
    }

/*
 ** This function integrates the time dependence of the "drift"-Hamiltonian.
 */
double csmComoveDriftFac(CSM csm,double dTime,double dDelta)
{
	double dOmega0 = csm->dOmega0;
	double dHubble0 = csm->dHubble0;
	double a0,A,B,a1,a2,eta1,eta2;

	if (!csm->bComove) return(dDelta);
	else if(csm->dLambda == 0.0 && csm->dOmegaRad == 0.0 &&
		csm->dQuintess == 0.0) {
		a1 = csmTime2Exp(csm,dTime);
		a2 = csmTime2Exp(csm,dTime+dDelta);
		if (dOmega0 == 1.0) {
			return((2.0/dHubble0)*(1.0/sqrt(a1) - 1.0/sqrt(a2)));
			}
		else if (dOmega0 > 1.0) {
			assert(dHubble0 >= 0.0);
			if (dHubble0 == 0.0) {
				A = 1.0;
				B = 1.0/sqrt(dOmega0);
				}
			else {
				a0 = 1.0/dHubble0/sqrt(dOmega0-1.0);
				A = 0.5*dOmega0/(dOmega0-1.0);
				B = A*a0;
				}
			eta1 = acos(1.0-a1/A);
			eta2 = acos(1.0-a2/A);
			return(B/A/A*(1.0/tan(0.5*eta1) - 1.0/tan(0.5*eta2)));
			}
		else if (dOmega0 > 0.0) {
			assert(dHubble0 > 0.0);
			a0 = 1.0/dHubble0/sqrt(1.0-dOmega0);
			A = 0.5*dOmega0/(1.0-dOmega0);
			B = A*a0;
			eta1 = acosh(a1/A+1.0);
			eta2 = acosh(a2/A+1.0);
			return(B/A/A*(1.0/tanh(0.5*eta1) - 1.0/tanh(0.5*eta2)));
			}
		else if (dOmega0 == 0.0) {
			/*
			 ** YOU figure this one out!
			 */
			assert(0);
			return(0.0);
			}
		else {
			/*
			 ** Bad value?
			 */
			assert(0);
			return(0.0);
			}
		}
	else{
	    return dRombergO(csm,
			     (double (*)(const void *, double)) csmComoveDriftInt,
			     1.0/csmTime2Exp(csm, dTime),
			     1.0/csmTime2Exp(csm, dTime + dDelta), EPSCOSMO);
	    }
	}


/*
 ** This function integrates the time dependence of the "kick"-Hamiltonian.
 */
double csmComoveKickFac(CSM csm,double dTime,double dDelta)
{
	double dOmega0 = csm->dOmega0;
	double dHubble0 = csm->dHubble0;
	double a0,A,B,a1,a2,eta1,eta2;

	if (!csm->bComove) return(dDelta);
	else if(csm->dLambda == 0.0 && csm->dOmegaRad == 0.0
		&& csm->dQuintess == 0.0) {
		a1 = csmTime2Exp(csm,dTime);
		a2 = csmTime2Exp(csm,dTime+dDelta);
		if (dOmega0 == 1.0) {
			return((2.0/dHubble0)*(sqrt(a2) - sqrt(a1)));
			}
		else if (dOmega0 > 1.0) {
			assert(dHubble0 >= 0.0);
			if (dHubble0 == 0.0) {
				A = 1.0;
				B = 1.0/sqrt(dOmega0);
				}
			else {
				a0 = 1.0/dHubble0/sqrt(dOmega0-1.0);
				A = 0.5*dOmega0/(dOmega0-1.0);
				B = A*a0;
				}
			eta1 = acos(1.0-a1/A);
			eta2 = acos(1.0-a2/A);
			return(B/A*(eta2 - eta1));
			}
		else if (dOmega0 > 0.0) {
			assert(dHubble0 > 0.0);
			a0 = 1.0/dHubble0/sqrt(1.0-dOmega0);
			A = 0.5*dOmega0/(1.0-dOmega0);
			B = A*a0;
			eta1 = acosh(a1/A+1.0);
			eta2 = acosh(a2/A+1.0);
			return(B/A*(eta2 - eta1));
			}
		else if (dOmega0 == 0.0) {
			/*
			 ** YOU figure this one out!
			 */
			assert(0);
			return(0.0);
			}
		else {
			/*
			 ** Bad value?
			 */
			assert(0);
			return(0.0);
			}
		}
	else{
	    return dRombergO(csm,
			     (double (*)(const void *, double)) csmComoveKickInt,
			     1.0/csmTime2Exp(csm, dTime),
			     1.0/csmTime2Exp(csm, dTime + dDelta), EPSCOSMO);
	    }
	}

double csmComoveLookbackTime2Exp(CSM csm,double dComoveTime)
{
	if (!csm->bComove) return(1.0);
	else {
	    double dExpOld = 0.0;
	    double dT0 = csmExp2Time(csm, 1.0);
	    double dTime = dT0 - dComoveTime;
	    double dExpNew;
	    int it = 0;
	    
	    if(dTime < EPSCOSMO) dTime = EPSCOSMO;
	    dExpNew = csmTime2Exp(csm, dTime);
	    /*
	     * Root find with Newton's method.
	     */
	    do {
		double dTimeNew = csmExp2Time(csm, dExpNew);
	    	double f = dComoveTime
		    - csmComoveKickFac(csm, dTimeNew, dT0 - dTimeNew);
		double fprime = -1.0/(dExpNew*dExpNew*csmExp2Hub(csm, dExpNew));
		dExpOld = dExpNew;
		dExpNew += f/fprime;
		it++;
		assert(it < 20);
		}
	    while (fabs(dExpNew - dExpOld)/dExpNew > EPSCOSMO);
	    return dExpNew;
	    }
	}

#if 0
/*
 * Code for generating linear growth factors.
 * This is taken from Carroll, Press & Turner Ann. Rev. (1992)
 */

/*
 * The following will not work with general quintessence models.
 * It is replaced with the code that follows.
 */
double csmGrowthFacInt(CSM csm, double dExp)
{
    return pow(dExp*csmExp2Hub(csm, dExp), -3.0);
    }

double csmGrowthFac(CSM csm, double dExp)
{
    return 2.5*csm->dOmega0*csm->dHubble0*csm->dHubble0*csmExp2Hub(csm, dExp)
	*dRombergO(csm, (double (*)(void *, double)) csmGrowthFacInt,
		   0.0, dExp, EPSCOSMO);
    }

/*
 * Time derivative of the growth factor.
 */

double csmGrowthFacDot(CSM csm, double dExp)
{
    double dHubble = csmExp2Hub(csm, dExp);
    
    return 2.5*csm->dOmega0*csm->dHubble0*csm->dHubble0
	*((csmExpDot2(csm, dExp) - dHubble*dHubble)
	  *dRombergO(csm, (double (*)(void *, double)) csmGrowthFacInt,
		     0.0, dExp, EPSCOSMO)
	  + 1.0/(dExp*dExp*dHubble));
    }
#endif

/*
 * delta[1] => deltadot
 */
void csmGrowthFacDeriv(CSM csm, double dlnExp, double *dlnDelta, double
		       *dlnDeltadot)
{
    double dExp = exp(dlnExp);
    double dHubble = csmExp2Hub(csm, dExp);
    
    dlnDeltadot[0] = dlnDelta[1];
    dlnDeltadot[1] = -dlnDelta[1]*dlnDelta[1]
	- dlnDelta[1]*(1.0 + csmExpDot2(csm, dExp)/(dHubble*dHubble))
	+ 1.5*csmExp2Om(csm, dExp);
    }
    
double csmGrowthFac(CSM csm, double dExp)
{
    const double dlnExpStart = -15;
    const int nSteps = 200;
    double dlnExp = log(dExp);
    double dDStart[2];
    double dDEnd[2];
    
    assert(dlnExp > dlnExpStart);
    
    dDStart[0] = dlnExpStart;
    dDStart[1] = 1.0;	/* Growing mode */

    RungeKutta(csm, (void (*)(void *, double, double *, double*))
	       csmGrowthFacDeriv, 2, dlnExpStart, dDStart, dlnExp, dDEnd,
	       nSteps);
    return exp(dDEnd[0]);
    }

double csmGrowthFacDot(CSM csm, double dExp)
{
    const double dlnExpStart = -15;
    const int nSteps = 200;
    double dlnExp = log(dExp);
    double dDStart[2];
    double dDEnd[2];
    
    assert(dlnExp > dlnExpStart);
    
    dDStart[0] = dlnExpStart;
    dDStart[1] = 1.0;	/* Growing mode */

    RungeKutta(csm, (void (*)(void *, double, double *, double*))
	       csmGrowthFacDeriv, 2, dlnExpStart, dDStart, dlnExp, dDEnd,
	       nSteps);
    return dDEnd[1]*csmExp2Hub(csm, dExp)*exp(dDEnd[0]);
    }

/*
 * expansion dependence of Omega_matter
 */

double csmExp2Om(CSM csm, double dExp)
{
    double dHubble = csmExp2Hub(csm, dExp);
    
    return csm->dOmega0*csm->dHubble0*csm->dHubble0
	/(dExp*dExp*dExp*dHubble*dHubble);
    }
