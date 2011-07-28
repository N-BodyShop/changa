/*
 * Supernova module originally written for GASOLINE
 * Modified for inclusion in ChaNGa
 *
 * Algorithm was originally implemented in TREESPH but is heavily
 * modified.  Contributors include Neal Katz, Eric Hayashi, Greg Stinson
 */
#include <math.h>
#include "ParallelGravity.h"
#include "imf.h"
#include "romberg.h"

// Private function to return result of basic IMF
double MillerScalo::returnimf(double mass)
{
    double dIMF;
    
    if(mass > mmax)
	return 0.0;
    if(mass > m3)
	dIMF = a3*pow(mass, b3);
    else if(mass > m2)
	dIMF = a2*pow(mass, b2);
    else if(mass > m1)
	dIMF = a1*pow(mass, b1);
    else
	dIMF = 0.0;
    return dIMF;
    }

double Kroupa93::returnimf(double mass)
{
    double dIMF;
    
    if(mass > mmax)
	return 0.0;
    if(mass > m3)
	dIMF = a3*pow(mass, b3);
    else if(mass > m2)
	dIMF = a2*pow(mass, b2);
    else if(mass > m1)
	dIMF = a1*pow(mass, b1);
    else
	dIMF = 0.0;
    return dIMF;
    }

// Chabrier imf function
double Chabrier::returnimf(double mass)
{
    double dIMF;

    if(mass > mmax)
	return 0.0;
    if(mass > m2)
	dIMF = a2*pow(mass, b2);
    else if(mass > m1)
	dIMF = a1 * exp(- pow(log10(mass) - log10(m1), 2.0)
			   /(2.0*b1*b1));
    else
	dIMF = 0.0;
   return dIMF;
}

double imfIntM(IMF *theimf, double logMass) 
{
    double mass = pow(10.0, logMass);
    return mass*theimf->returnimf(mass);
    }

/*
 * Cumulative number of stars with mass greater than "mass".
 */
double MillerScalo::CumNumber(double mass)
{
    double dCumN;
    
    if(mass > mmax) return 0;
    if(mass > m3) return a3/b3*(pow(mmax, b3) - pow(mass, b3))/log(10.0);
    else dCumN = a3/b3*(pow(mmax, b3) - pow(m3, b3))/log(10.0); 
    if(mass > m2) {
	dCumN += a2/b2*(pow(m3, b2) - pow(mass, b2))/log(10.0);
	return dCumN;
    }
    else dCumN += a2/b2*(pow(m3, b2) - pow(m2, b2))/log(10.0);

    if(mass > m1) dCumN += a1/b1*(pow(m2, b1) - pow(mass, b1))/log(10.0);
    else dCumN += a1/b1*(pow(m2, b1) - pow(m1, b1))/log(10.0);
    
    return dCumN;
    }

double Kroupa93::CumNumber(double mass)
{
    double dCumN;
    
    if(mass > mmax) return 0;
    if(mass > m3) return a3/b3*(pow(mmax, b3) - pow(mass, b3))/log(10.0);
    else dCumN = a3/b3*(pow(mmax, b3) - pow(m3, b3))/log(10.0); 
    if(mass > m2) {
	dCumN += a2/b2*(pow(m3, b2) - pow(mass, b2))/log(10.0);
	return dCumN;
    }
    else dCumN += a2/b2*(pow(m3, b2) - pow(m2, b2))/log(10.0);

    if(mass > m1) dCumN += a1/b1*(pow(m2, b1) - pow(mass, b1))/log(10.0);
    else dCumN += a1/b1*(pow(m2, b1) - pow(m1, b1))/log(10.0);
    
    return dCumN;
    }

MillerScalo* MillerScalo::clone() const
{
    return new MillerScalo(*this);
    }

Kroupa93* Kroupa93::clone() const
{
    return new Kroupa93(*this);
    }

Chabrier* Chabrier::clone() const
{
    return new Chabrier(*this);
    }

double Chabrier::CumNumber(double mass)
{
    double dCumN;
    
    if(mass > mmax) return 0;
    if(mass > m2) return a2/b2*(pow(mmax, b2) - pow(mass, b2))/log(10.0);
    else dCumN = a2/b2*(pow(mmax, b2) - pow(m2, b2))/log(10.0); 
    /*
     * Introduce the auxilary variable
     * y \equiv (log(m) - log(m1))/(sqrt(2) b1)
     * to evaluate the integral
     */
    {
	double ymin, ymax;
	ymax = (log10(m2) - log10(m1))/(M_SQRT2*b1);
	if(mass > m1)
	    ymin = (log10(mass) - log10(m1))/(M_SQRT2*b1);
	else
	    ymin = 0.0;
	dCumN += a1*b1*sqrt(M_PI)*M_SQRT1_2*(erf(ymax) - erf(ymin));
	}
    return dCumN;
}

/*
 * Cumulative mass of stars with mass greater than "mass".
 */
double MillerScalo::CumMass(double mass)
{
    double dCumM;
    
    if(mass > mmax) return 0;
    if(mass > m3)
	return a3/(b3 + 1)*(pow(mmax, b3 + 1)
			    - pow(mass, b3 + 1))/log(10.0);
    else
	dCumM = a3/(b3 + 1)*(pow(mmax, b3 + 1)
			     - pow(m3, b3 + 1))/log(10.0); 
    if(mass > m2) {
	dCumM += a2/(b2 + 1)*(pow(m3, b2 + 1)
			      - pow(mass, b2 + 1))/log(10.0);
	return dCumM;
	}
    else {
	dCumM += a2/(b2 + 1)*(pow(m3, b2 + 1)
			      - pow(m2, b2 + 1))/log(10.0);
	}
    if(mass > m1)
	dCumM += a1/(b1 + 1)*(pow(m2, b1 + 1)
			      - pow(mass, b1 + 1))/log(10.0);
    else
	dCumM += a1/(b1 + 1)*(pow(m2, b1 + 1)
			      - pow(m1, b1 + 1))/log(10.0);

    return dCumM;
    }

double Kroupa93::CumMass(double mass)
{
    double dCumM;
    
    if(mass > mmax) return 0;
    if(mass > m3)
	return a3/(b3 + 1)*(pow(mmax, b3 + 1)
			    - pow(mass, b3 + 1))/log(10.0);
    else
	dCumM = a3/(b3 + 1)*(pow(mmax, b3 + 1)
			     - pow(m3, b3 + 1))/log(10.0); 
    if(mass > m2) {
	dCumM += a2/(b2 + 1)*(pow(m3, b2 + 1)
			      - pow(mass, b2 + 1))/log(10.0);
	return dCumM;
	}
    else {
	dCumM += a2/(b2 + 1)*(pow(m3, b2 + 1)
			      - pow(m2, b2 + 1))/log(10.0);
	}
    if(mass > m1)
	dCumM += a1/(b1 + 1)*(pow(m2, b1 + 1)
			      - pow(mass, b1 + 1))/log(10.0);
    else
	dCumM += a1/(b1 + 1)*(pow(m2, b1 + 1)
			      - pow(m1, b1 + 1))/log(10.0);

    return dCumM;
    }

double Chabrier::CumMass(double mass)
{
  double dCumM;
    if(mass > m2)
	return a2/(b2 + 1)*(pow(mmax, b2 + 1)
				  - pow(mass, b2 + 1))/log(10.0);
    else
	dCumM = a2/(b2 + 1)*(pow(mmax, b2 + 1)
				   - pow(m2, b2 + 1))/log(10.0); 
    /*
     * Evaluate the integral numerically in log m.
     */
    {
	double xmin;
	const double EPS_IMF = 1e-7;
	if(mass > m1)
	    xmin = log10(mass);
	else
	    xmin = log10(m1);
	dCumM += dRombergO(this, (double (*)(void *, double)) imfIntM, xmin, log10(m2), EPS_IMF);
	}
    return dCumM;
}

#ifdef IMF_TST
#include "testimf.decl.h"

int
main(int argc, char **argv)
{
    int nsamp;
    int i, imax;
    double dlgm;
    double lgm;
    double Ntot;
    double Mtot;
    double dMassT1, dMassT2, dMass, part, sum;
    
    
    MillerScalo msimf;
    Kroupa93 krimf;
    Chabrier chimf;

    sum = 0;
#if 0
    dMassT1 = 3;
    dMassT2 = 8;
    imax = 101;
    for (i = 0; i < imax; i++) {
        dMass = dMassT1 + i*(dMassT2-dMassT1)/(float) imax;
        part = imfSec (MSparam, dMass);
        printf ("%g %g\n", dMass, part);
        sum += part*0.05;
    }
    printf ("sum = %g\n", sum);

    printf ("number SN Ia = %g\n", dNSNIa (MSparam, dMassT1, dMassT2));
    printf ("mass in SN Ia %g\n", dMSNIa (MSparam, dMassT1, dMassT2));
#endif

    CkAssert(argc == 2);
    
    nsamp = atoi(argv[1]);
    dlgm = (2.0 + 1.0)/nsamp;
    
    Ntot = msimf.CumNumber(0.0791);
    Mtot = msimf.CumMass(0.0791);
    printf("# MS Ntot, Mtot = %g %g\n", Ntot, Mtot);
    printf("# MS Fraction of mass in stars that go Type II SN: %g\n",
	   msimf.CumMass(8.0)/Mtot);
    printf("# MS SuperNovea/solar mass of stars formed: %g\n",
	   msimf.CumNumber(8.0)/Mtot);

    Ntot = krimf.CumNumber(0.0791);
    Mtot = krimf.CumMass(0.0791);
    printf("# Kroupa Ntot, Mtot = %g %g\n", Ntot, Mtot);
    printf("# Kroupa Fraction of mass in stars that go Type II SN: %g\n",
	   krimf.CumMass(8.0)/Mtot);
    printf("# Kroupa SuperNovea/solar mass of stars formed: %g\n",
	   krimf.CumNumber(8.0)/Mtot);

    Ntot = chimf.CumNumber(0.0791);
    Mtot = chimf.CumMass(0.0791);
    printf("# Chabrier Ntot, Mtot = %g %g\n", Ntot, Mtot);
    printf("# Chabrier Fraction of mass in stars that go Type II SN: %g\n",
	   chimf.CumMass(8.0)/Mtot);
    printf("# Chabrier SuperNovea/solar mass of stars formed: %g\n",
	   chimf.CumNumber(8.0)/Mtot);
    
    printf("# mass MSIMF MSCumNumber MSCumMass KrIMF KrCumNumber KrCumMass ChIMF ChCumNumber ChCumMass\n");
    for(i = 0; i < nsamp; i++) {
	double mass;
	
	lgm = -1.1 + i*dlgm;

	mass = pow(10.0, lgm);


	printf("%g %g %g %g %g %g %g %g %g %g\n", mass,
	       msimf.returnimf(mass), msimf.CumNumber(mass), msimf.CumMass(mass),
	       krimf.returnimf(mass), krimf.CumNumber(mass), krimf.CumMass(mass),
	       chimf.returnimf(mass), chimf.CumNumber(mass), chimf.CumMass(mass));
	}
    return 0;
    }

#include "testimf.def.h"
#endif
