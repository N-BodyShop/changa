#include <math.h>
#include "ParallelGravity.h"
#include "starlifetime.h"
#include "supernova.h"
#include "romberg.h"

/* Calculation of number and mass of stars that go SN Type Ia in a
   given stellar mass range.  */

#define EPSSNIA 1e-6

/// class to hold information about the secondary
class MSBinary 
{
public:
    double dMass2;  // mass of secondary
    double dGamma_2nd;  // power law index of secondary distribution
    IMF *imf;
    };
    
/// Integrand for secondary mass function: probability of a binary
/// of mass dMassB having a secondary of MSBinary->dMass2
double dMSIMFSecInt(const MSBinary *msb, double dMassB)
{
    // A factor of 1/(dMassB*log(10)) is from the definition of IMF as
    // number per log10(M).  Another factor of 1/dMassB is to
    // normalized the integral over mass ratios.
    return pow(msb->dMass2/dMassB, msb->dGamma_2nd)
	*msb->imf->returnimf(dMassB)/dMassB/dMassB/log(10.0);
    }

/** IMF of secondary in binary system that goes SN Ia
 * The distribution of secondary mass ratios is assumed to be a power
 * law, mu^dGamma as in Matteucci & Greggio (1986)
 */
double dMSIMFSec(const SN *sn, double dMass2)
{
    double dIMFSec;
    double dMass2_2, Msup, Minf;
    const double dGamma_2nd = 2.0;
    const double dNorm_2nd = pow(2.0, 1 + dGamma_2nd)*(1 + dGamma_2nd);  // Normalization
									 // of 2ndary 

    MSBinary msb = {dMass2, dGamma_2nd, sn->imf};

    Msup = dMass2 + sn->dMSNIImin;  // Mass where primary would have gone supernovaII
    // Msup = 16.0;  // Note that the double integral can be tested by
    // uncommenting the above line, setting dFracBinSNIa = 1.0, and
    // comparing NSNIa with the number of stars in the 3.0-16.0 mass interval.
    dMass2_2 = 2.*dMass2;  // Minimum mass of binary
    Minf = (dMass2_2 > sn->dMBmin)?dMass2_2:sn->dMBmin;
    dIMFSec = dRombergO(&msb, (double (*)(const void *, double))dMSIMFSecInt, Minf, Msup, EPSSNIA);
    dIMFSec *= sn->dFracBinSNIa * dNorm_2nd;
    return dIMFSec;
    }

/** calculate number of SN Type Ia a la Raiteri, Villata, Navarro, A&A
  * 315, 105, 1996) Returns number of SN Type Ia that occur during
  * timestep in which dMassT1 and dMassT2 are masses of stars that end
  * their lives at the end and beginning of timestep, respectively
  */
double SN::NSNIa (double dMassT1, double dMassT2) const
{
    CkAssert (dMassT1 < dMassT2);
    // Exclude primaries that go SNII
    if (dMassT1 >= 0.5*dMBmax)
	return 0.0;
    if(dMassT2 > 0.5*dMBmax)
	dMassT2 = 0.5*dMBmax;
    return dRombergO(this, (double (*)(const void *, double))dMSIMFSec, dMassT1, dMassT2, EPSSNIA);
    }

#ifdef SNIA_TST
#include "testsnia.decl.h"
#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))
#include "feedback.h"

int
main(int argc, char **argv)
{
    Chabrier chimf; // Chabrier has power law index of -1.3 which
		    // matches s = -2.35 in Greggio & Renzini, 1983
    // Kroupa93 chimf; // Kroupa is what RVN use.
    // Kroupa01 chimf;
    // MillerScalo chimf;
    SN sn;
    sn.imf = &chimf;
    Padova pdva;
    double z = 0.02;		// metalicity: use value to compare
				// with Greggio & Renzini, 1983

    double nsamp = 100;
    double tfac = log(14.0e9/1e6)/nsamp;  /// equal log interavals from 1Myr
				    /// to 14Gyr
    
    // sn.dFracBinSNIa = 1.00;
    printf("# Total Ia, II supernova: %g %g, SNIa mass range: %g\n",
	   sn.NSNIa(0.8, 8.0)/chimf.CumMass(0.0),
	   (chimf.CumNumber(8.0) - chimf.CumNumber(40.0))/chimf.CumMass(0.0),
	   (chimf.CumNumber(3.0) - chimf.CumNumber(16.0))/chimf.CumMass(0.0));
    {
	// One solar mass formed at t = 0, with metallicity 0.02
	SFEvent sfEvent(1.0, 0.0, 0.02, .005, .005);
	FBEffects fbEffectsII;
	FBEffects fbEffectsIa;
	double t, deltat;
	double dMOxII = 0.0;
	double dMFeII = 0.0;
	double dMOxIa = 0.0;
	double dMFeIa = 0.0;
	
	deltat = 1e6;
	for (int i = 0; i < 10000; i++) {
	    t = i*deltat;
	    sn.CalcSNIIFeedback(&sfEvent, t, deltat, &fbEffectsII);
	    dMOxII += fbEffectsII.dMassLoss*fbEffectsII.dMOxygen;
	    dMFeII += fbEffectsII.dMassLoss*fbEffectsII.dMIron;
	    sn.CalcSNIaFeedback(&sfEvent, t, deltat, &fbEffectsIa);
	    dMOxIa += fbEffectsIa.dMassLoss*fbEffectsIa.dMOxygen;
	    dMFeIa += fbEffectsIa.dMassLoss*fbEffectsIa.dMIron;
	    }
	printf("# Total II Ox: %g, Fe: %g\n", dMOxII, dMFeII);
	printf("# Total Ia Ox: %g, Fe: %g\n", dMOxIa, dMFeIa);
	}
    
    for(int i = 0; i < nsamp; i++) {
	double t = 1e6*exp(i*tfac);
	double deltat = 0.01*t;  /// interval to determine rate
	double dMaxMass = pdva.StarMass(t, z); // stars dying at start
	double dMinMass = pdva.StarMass(t+deltat, z); // stars dying
						      // at end
	double dMStarMinII = max (8.0, dMinMass); // range of SNII
	double dMStarMaxII = min (40.0, dMaxMass);
	double dCumNMinII = chimf.CumNumber(dMStarMinII); 
	double dCumNMaxII = chimf.CumNumber(dMStarMaxII);
	double dMtot = chimf.CumMass(0.0);
	double dNSNII;
	// One solar mass formed at t = 0, with metallicity 0.02
	SFEvent sfEvent(1.0, 0.0, 0.02, .005, .005);
	FBEffects fbEffectsII;
	FBEffects fbEffectsIa;
	
	if(dMaxMass > 8.0 && dMinMass < 40.0) {
	    dNSNII = (dCumNMinII - dCumNMaxII)/dMtot/deltat;
	    }
	else dNSNII = 0.0;
	
	sn.CalcSNIIFeedback(&sfEvent, t, deltat, &fbEffectsII);
	sn.CalcSNIaFeedback(&sfEvent, t, deltat, &fbEffectsIa);
	
	if (dMaxMass > dMinMass)
	    printf ("%g %g %g %g %g %g %g %g %g %g %g %g %g\n", t,
		    sn.NSNIa(dMinMass, dMaxMass)/deltat/dMtot, dNSNII,
		    fbEffectsIa.dEnergy, fbEffectsII.dEnergy,
		    fbEffectsIa.dMassLoss, fbEffectsII.dMassLoss,
		    fbEffectsIa.dMetals, fbEffectsII.dMetals,
		    fbEffectsIa.dMIron, fbEffectsII.dMIron,
		    fbEffectsIa.dMOxygen, fbEffectsII.dMOxygen
		    );
	}
    }
#include "testsnia.def.h"
#endif

