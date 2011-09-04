/*
 * Supernova module originally written for GASOLINE
 * Modified for inclusion in ChaNGa
 *
 * Algorithm was originally implemented in TREESPH but is heavily
 * modified.  Contributors include Neal Katz, Eric Hayashi, Greg Stinson
 */
#include <math.h>
#include "ParallelGravity.h"
#include "feedback.h"
#include "supernova.h"
#include "smooth.h"
#include "physconst.h"

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))
double mod(double a, int b) {return (a-b*floor(a/b));}

/*
 * Supernova module for GASOLINE
 */

void SN::CalcSNIIFeedback(SFEvent *sfEvent,
			    double dTime, /* current time in years */
			    double dDelta, /* length of timestep (years) */
			    FBEffects *fbEffects)
{
    double dMSNTypeII, dNSNTypeII, dMeanMStar;
    /* stellar lifetimes corresponding to beginning and end of 
       current timestep with respect to starbirth time in yrs */
    double dStarLtimeMin = dTime - sfEvent->dTimeForm; 
    double dStarLtimeMax = dStarLtimeMin + dDelta;
    
    /* masses corresponding to these stellar lifetimes in solar masses */
    double dMStarMin = pdva.StarMass(dStarLtimeMax, sfEvent->dMetals); 
    double dMStarMax = pdva.StarMass(dStarLtimeMin, sfEvent->dMetals); 

    double dMtot = imf->CumMass(0.0); /* total mass in stars integrated over IMF */

    /* Quantize feedback */
    if (iNSNIIQuantum && dMStarMin > dMSNIImax && dStarLtimeMin){
	/* blow wind before SN */
	
	fbEffects->dMetals = 0.0; /* zero enrichment from winds */
	fbEffects->dMIron = 0.0;
	fbEffects->dMOxygen = 0.0;
	
	dNSNTypeII = imf->CumNumber(dMSNIImin);
	dNSNTypeII *= sfEvent->dMass/dMtot; /* normalize to star particle mass */
	fbEffects->dMassLoss = dNSNTypeII * dDelta / 1e6; /* 1 M_sol / Myr / SN */
	fbEffects->dEnergy = dNSNTypeII * 3e42 * dDelta/ /* 1e35 erg/s /SN */
	    (MSOLG*fbEffects->dMassLoss); 
	
	} else if (dMStarMin > dMSNIImax  || dMStarMax < dMSNIImin) {
	/* do nothing */
	fbEffects->dMassLoss = 0.0;
	fbEffects->dEnergy = 0.0;
	fbEffects->dMetals = 0.0;
	fbEffects->dMIron = 0.0;
	fbEffects->dMOxygen = 0.0;
	return;
	} else {
	/* supernova */
	double dMStarMinII = max (dMSNIImin, dMStarMin); 
	double dMStarMaxII = min (dMSNIImax, dMStarMax);
	
	CkAssert (dMStarMinII < dMStarMaxII && 
		  dMStarMinII >0.0 && dMStarMaxII > 0.0);
	
	/* cumulative mass of stars with mass greater than dMStarMinII 
	   and dMStarMaxII in solar masses */
	double dCumMMin = imf->CumMass (dMStarMinII);
	double dCumMMax = imf->CumMass (dMStarMaxII);
	
	if(dCumMMax > dCumMMin || dCumMMax < 0) dMSNTypeII = dCumMMin;
	else dMSNTypeII = dCumMMin - dCumMMax; /* mass of stars that go SN II
						  in this timestep normalized to IMF */
	
	dMSNTypeII *= sfEvent->dMass/dMtot; /* REAL MASS of stars that go SNII */
	
	/* cumulative number of stars with mass greater than dMStarMinII and
	   less than dMstarMaxII in solar masses */
	double dCumNMinII = imf->CumNumber (dMStarMinII); 
	double dCumNMaxII = imf->CumNumber (dMStarMaxII);
	
	if(dCumNMaxII > dCumNMinII || dCumNMaxII < 0) dNSNTypeII = dCumNMinII;
	else dNSNTypeII = dCumNMinII - dCumNMaxII;
	dNSNTypeII *= sfEvent->dMass/dMtot; /* number of SNII in star particle */
	
	/* Average Star mass for metalicity calculation.  Metals always distributed
	 * even when energy is quantized, so grab NSNtypeII before it is quantized.
	 */
	if (dNSNTypeII > 0) dMeanMStar = dMSNTypeII/dNSNTypeII;
	else dMeanMStar =0;
	
	/*
	 * Quantized Feedback:  Not recommended until star particles ~< 100 M_sun
	 */
	/* Blow winds before SN (power of winds ~ power of SN) */
	double dMinAge = pdva.Lifetime(dMSNIImax, sfEvent->dMetals); 
	if (dNSNTypeII > 0 && iNSNIIQuantum > 0 && dStarLtimeMin > dMinAge) {
	    /* Make sure only a iNSNIIQuantum number of
	     * SNII go off at a time */
	    double dProb = mod(dNSNTypeII, iNSNIIQuantum)/iNSNIIQuantum;
	    double dRandomNum = (rand()/((double) RAND_MAX));
	    /*	  printf("Random Number = %g\n",dRandomNum);*/
	    if(dRandomNum < dProb) /* SN occurred */ {
		/* Adds missing part to make up quantum */
		dNSNTypeII += (1.-dProb)*iNSNIIQuantum;
		/*printf("NSN: +factor=%g   dNSNTypeII=%g  result=%g fNSN=%g\n",(1.-dProb)*iNSNIIQuantum,dNSNTypeII,dNSNTypeII + (1.-dProb)*iNSNIIQuantum,p->fNSN);*/
		} else {
		dNSNTypeII -= dProb*iNSNIIQuantum;
		/*printf("NSN: -factor=%g   dNSNTypeII=%g  result=%g fNSN=%g\n",dProb*iNSNIIQuantum,dNSNTypeII,dNSNTypeII - dProb*iNSNIIQuantum,p->fNSN);*/
		}
	    if(dNSNTypeII < iNSNIIQuantum) dNSNTypeII = 0;
	    } 
	
	/* decrement mass of star particle by mass of stars that go SN
	   plus mass of SN remnants */
	double dDeltaMSNrem = dNSNTypeII*dMSNrem; /* mass in SN remnants */
	fbEffects->dMassLoss = dMSNTypeII - dDeltaMSNrem;
	
	/* SN specific energy rate to be re-distributed among neighbouring gas
	   particles */
	double dESNTypeII = dNSNTypeII * dESN;
	fbEffects->dEnergy = dESNTypeII/(MSOLG*fbEffects->dMassLoss);   
	
	/* fraction of mass in metals to be re-distributed among neighbouring
	 * gas particles.  Formula 6-8 of  Raiteri, Villata and Navarro, A&A
	 * 315, 105, 1996  are used to calculate SNII yields
	 * Integrate over power law from lowest mass progenitor to highest mass.
	 */
	fbEffects->dMIron = dMFeconst * pow(dMeanMStar, dMFeexp)
	    / (dMEjconst * pow(dMeanMStar, dMEjexp));
	fbEffects->dMOxygen = dMOxconst*pow(dMeanMStar, dMOxexp)
	    / (dMEjconst * pow(dMeanMStar, dMEjexp));
	// Use ratio of Fe to total iron group and O to total non-iron
	// group derived from Asplund et al 2009
	fbEffects->dMetals = 1.06*fbEffects->dMIron + 1.48*fbEffects->dMOxygen;
	}
    }

void SN::CalcSNIaFeedback(SFEvent *sfEvent,
			double dTime, /* current time in years */
			double dDelta, /* length of timestep (years) */
			FBEffects *fbEffects)
{
    /* stellar lifetimes corresponding to beginning and end of 
     * current timestep with respect to starbirth time in yrs */
    double dMinLifetime = dTime - sfEvent->dTimeForm; 
    double dMaxLifetime = dMinLifetime + dDelta;
    
    double dMinMass = pdva.StarMass(dMaxLifetime, sfEvent->dMetals); 
    double dMaxMass = pdva.StarMass(dMinLifetime, sfEvent->dMetals); 
    
    double dTotalMass = imf->CumMass(0.0); /* total mass in stars integrated over IMF */
    
    if (dMinMass < dMBmax/2.) {
	double dMStarMinIa = dMinMass; 
	double dMStarMaxIa = min (dMBmax/2., dMaxMass); 
	
	CkAssert (dMStarMinIa < dMStarMaxIa && dMStarMinIa >0.0 && dMStarMaxIa > 0.0);
	
	/* number of stars that go SNIa */        
	double dNSNTypeIa = NSNIa (dMStarMinIa, dMStarMaxIa); 
	dNSNTypeIa /= dTotalMass;	/* convert to number per solar mass of stars */
	dNSNTypeIa *= sfEvent->dMass; /* convert to absolute number of SNIa */
	
	double dESNTypeIa = dNSNTypeIa * dESN;
	/* decrement mass of star particle by mass of stars that go SN
	   and SN Ia have no remnants.   
	   The MSNrem is the same Chandrasekhar mass that explodes as a SNIa.*/
	fbEffects->dMassLoss = dNSNTypeIa*dMSNrem;
	
	/* SN specific energy rate to be re-distributed among neighbouring gas
	   particles */
	fbEffects->dEnergy = dESNTypeIa/(MSOLG*fbEffects->dMassLoss);   
	
	/* Following Raiteri 1996 who follows Thielemann's 1986 W7 model for
	 * SNIa explosions, the same mass of iron and oxygen is released in
	 * every single explosion.  A comparable amount of Silicon is ejected
	 * to the amount of Oxygen that is ejected.
	 */
	fbEffects->dMIron = dNSNTypeIa*0.63/fbEffects->dMassLoss;
	fbEffects->dMOxygen = dNSNTypeIa*0.13/fbEffects->dMassLoss;
	/* Fraction of mass in metals to be re-distributed among neighbouring
	 * gas particles: assumes fixed amount of metals per
	 * supernovea independent of mass. See Raiteri, Villata and
	 * Navarro, page 108.
	 */
	// Use total metals to Fe and O based on Asplund et al 2009
	fbEffects->dMetals = dNSNTypeIa*(0.63*1.06 + 0.13*1.48)
	    /fbEffects->dMassLoss; 
	} else {
	fbEffects->dMassLoss = 0.0;
	fbEffects->dEnergy = 0.0;    
	fbEffects->dMetals = 0.0; 
	fbEffects->dMIron = 0.0;
	fbEffects->dMOxygen = 0.0;
	}
    return;
    
}
