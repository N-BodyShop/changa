#include <math.h>
#include "ParallelGravity.h"
#include "feedback.h"
#include "starlifetime.h"
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
double dMSIMFSec(Fdbk *fb, double dMass2)
{
    double dIMFSecExp, dIMFSecIntExp, dIMFSec;
    double dMass2_2, Msup, Minf;
    const double dGamma_2nd = 2.0;
    const double dNorm_2nd = pow(2.0, 1 + dGamma_2nd)*(1 + dGamma_2nd);  // Normalization
									 // of 2ndary 

    MSBinary msb = {dMass2, dGamma_2nd, fb->imf};

    Msup = dMass2 + 8.0;  // Mass where primary would have gone supernovaII
    dMass2_2 = 2.*dMass2;  // Minimum mass of binary
    Minf = (dMass2_2 > 3.0)?dMass2_2:3.0;
    dIMFSec = dRombergO(&msb, (double (*)(void *, double))dMSIMFSecInt, Minf, Msup, EPSSNIA);
    dIMFSec *= fb->dFracBinSNIa * dNorm_2nd;
    return dIMFSec;
    }

/** calculate number of SN Type Ia a la Raiteri, Villata, Navarro, A&A
  * 315, 105, 1996) Returns number of SN Type Ia that occur during
  * timestep in which dMassT1 and dMassT2 are masses of stars that end
  * their lives at the end and beginning of timestep, respectively
  */
double Fdbk::NSNIa (double dMassT1, double dMassT2)
{
    CkAssert (dMassT1 < dMassT2);
    return dRombergO(this, (double (*)(void *, double))dMSIMFSec, dMassT1, dMassT2, EPSSNIA);
    }

#if 0
#define STARFORM_HINCLUDED
#define STARFORM
#define GASOLINE
int main ()
{
    SN sn;
    MSPARAM MSparam;
    struct inMSIMFSec mssn;
	//    MSSN mssn;
    
    MSInitialize (&MSparam);
    // crap();
    
    snInitialize (&sn);
	//    snInitConstants (sn);
    mssn.ms = *MSparam;
    mssn.sn = *sn;
    printf ("%g %g\n", mssn.ms.a1, mssn.sn.dFracBinSNIa);
    
	}
#endif

