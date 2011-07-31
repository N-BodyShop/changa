#include <math.h>
#include "ParallelGravity.h"
#include "starlifetime.h"
#include "supernova.h"
#include "romberg.h"

/* Calculation of number and mass of stars that go SN Type Ia in a
   given stellar mass range.  */

#define EPSSNIA 1e-7

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
double dMSIMFSecInt(MSBinary *msb, double dMassB) 
{
    return pow(msb->dMass2/dMassB, msb->dGamma_2nd)
	*msb->imf->returnimf(dMassB)/dMassB;
    }

/** IMF of secondary in binary system that goes SN Ia
 * The distribution of secondary mass ratios is assumed to be a power
 * law, mu^dGamma as in Matteucci & Greggio (1986)
 */
double dMSIMFSec(SN *sn, double dMass2)
{
    double dIMFSec;
    double dMass2_2, Msup, Minf;
    const double dGamma_2nd = 2.0;
    const double dNorm_2nd = pow(2.0, 1 + dGamma_2nd)*(1 + dGamma_2nd);  // Normalization
									 // of 2ndary 

    MSBinary msb = {dMass2, dGamma_2nd, sn->imf};

    Msup = dMass2 + 8.0;  // Mass where primary would have gone supernovaII
    dMass2_2 = 2.*dMass2;  // Minimum mass of binary
    Minf = (dMass2_2 > 3.0)?dMass2_2:3.0;
    dIMFSec = dRombergO(&msb, (double (*)(void *, double))dMSIMFSecInt, Minf, Msup, EPSSNIA);
    dIMFSec *= sn->dFracBinSNIa * dNorm_2nd;
    return dIMFSec;
    }

/** calculate number of SN Type Ia a la Raiteri, Villata, Navarro, A&A
  * 315, 105, 1996) Returns number of SN Type Ia that occur during
  * timestep in which dMassT1 and dMassT2 are masses of stars that end
  * their lives at the end and beginning of timestep, respectively
  */
double SN::NSNIa (double dMassT1, double dMassT2)
{
    CkAssert (dMassT1 < dMassT2);
    // Exclude primaries that go SNII
    if (dMassT1 >= 0.5*dMBmax)
	return 0.0;
    if(dMassT2 > 0.5*dMBmax)
	dMassT2 = 0.5*dMBmax;
    return dRombergO(this, (double (*)(void *, double))dMSIMFSec, dMassT1, dMassT2, EPSSNIA);
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
    SN sn;
    sn.imf = &chimf;
    Padova pdva;
    double z = 0.02;		// metalicity: use value to compare
				// with Greggio & Renzini, 1983

    double nsamp = 100;
    double tfac = log(14.0e9/1e6)/nsamp;  /// equal log interavals from 1Myr
				    /// to 14Gyr
    
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

