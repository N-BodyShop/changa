#include <math.h>
#include "ParallelGravity.h"
#include "feedback.h"
#include "starlifetime.h"
#include "romberg.h"

/* Calculation of number and mass of stars that go SN Type Ia in a
   given stellar mass range.  */

#define EPSSNIA 1e-7

double dMSIMFSec(Fdbk *fb, double dMass2)
{
	/* IMF of secondary in binary system that goes SN Ia */
    double dIMFSecExp, dIMFSecIntExp, dIMFSec;
    double dMass2_2, Msup, Minf;
    
    dIMFSecExp = fb->imf->Oneto8Exp() - 3.; /* subtract 1 from exponent in integrand
                                      because MS IMF is per unit log mass,
                                      subtract 2 because of square of ratio of mass
                                      of secondary to mass of binary */
    dIMFSecIntExp = dIMFSecExp + 1.; /* exponent of anti-derivative is +1 */
    Msup = dMass2 + 8;
    dMass2_2 = 2.*dMass2;    
    Minf = (dMass2_2 > 3)?dMass2_2:3;
    dIMFSec = pow (Msup, dIMFSecIntExp) - pow(Minf, dIMFSecIntExp);
    dIMFSec *= fb->dFracBinSNIa * fb->imf->Oneto8PreFactor() * dMass2*dMass2 / dIMFSecIntExp;
    return dIMFSec;
    }

double Fdbk::NSNIa (double dMassT1, double dMassT2)
{
    CkAssert (dMassT1 < dMassT2 && dMassT1 >= dMBmin && dMassT2 <= dMBmax/2.);
    
	/* calculate number of SN Type Ia a la Raiteri, Villata, Navarro, A&A
	   315, 105, 1996) Returns number of SN Type Ia that occur during
	   timestep in which dMassT1 and dMassT2 are masses of stars that end
	   their lives at the end and beginning of timestep, respectively */
    return dRombergO(this, (double (*)(void *, double))dMSIMFSec, dMassT1, dMassT2, EPSSNIA);
	}

/* Not really used any more because mass SNIa ejecta = NSNIa*1.4 */
double MSIMFSecM(Fdbk *fb, double dMass2)
{
	/* Mass times IMF of secondary in binary system that goes SN Ia */
    double dIMFSecExp, dIMFSecIntExp, dIMFSec;
    double dMass2_2, Msup, Minf;
    
    dIMFSecIntExp = fb->imf->Oneto8Exp() - 2.; /* subtract 1 from exponent in integrand
                                      because MS IMF is per unit log mass,
                                      subtract 2 because of square of ratio of mass
                                      of secondary to mass of binary, add 1 to multiply 
                                      IMF by mass */
    /*dIMFSecIntExp = dIMFSecExp + 1.;  exponent of anti-derivative is +1 */
    Msup = dMass2 + 8;
    dMass2_2 = 2.*dMass2;    
    Minf = (dMass2_2 > 3)?dMass2_2:3;
    dIMFSec = pow (Msup, dIMFSecIntExp) - pow(Minf, dIMFSecIntExp);
    dIMFSec *= fb->dFracBinSNIa * fb->imf->Oneto8PreFactor() * dMass2*dMass2*dMass2 / dIMFSecIntExp;
    return dIMFSec;
    }

double Fdbk::MSNIa (double dMassT1, double dMassT2)
{
    CkAssert (dMassT1 < dMassT2 && dMassT1 >= dMBmin/2. && dMassT2 <= dMBmax/2.);
    
	/* calculate mass of stars that go SN Type Ia a la Raiteri, Villata,
	   Navarro, A&A 315, 105, 1996) Returns total mass in stars that go SN
	   Type Ia during timestep in which dMassT1 and dMassT2 are masses of
	   stars that end their lives at the end and beginning of timestep,
	   respectively */
    return dRombergO(this, (double (*)(void *, double)) MSIMFSecM, dMassT1, dMassT2, EPSSNIA);
	}

/* XXX - move dMBmin, dMBmax in dNSNIa and dMSNIa
   XXX - move dFracBinSNIa in dMSIMFSec and dMSIMFSecM
   */



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

