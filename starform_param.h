#ifndef STARFORM_PARAM_HINCLUDED
#define STARFORM_PARAM_HINCLUDED

#include "parameters.h"

class StfmParam {
 public:
    double dDeltaStarForm;	/* timestep in system units */
    int iStarFormRung;		/* rung for star formation */
    double dGmUnit;		/* system mass in grams */
    double dGmPerCcUnit;	/* system density in gm/cc */
    double dSecUnit;		/* system time in seconds */
    double dErgUnit;		/* system energy in ergs */
    double dPhysDenMin;		/* Physical density minimum for star
				   formation (in system units) */
    double dOverDenMin;		/* Overdensity minimum for star formation */
    double dTempMax;		/* Form stars below this temperature
				   EVEN IF the gas is not cooling. */
    double dSoftMin;		/* Jean's length as a fraction of
				   softening at which to form stars*/
    double dCStar;		/* Star formation constant */
    double dStarEff;		/* Fraction of gas mass converted into
				 star mass per timestep. */
    double dInitStarMass;       /* Fixed Initial Star Mass */
    double dMinSpawnStarMass;   /* Minimum Initial Star Mass */
    double dMinGasMass;		/* minimum mass gas before we delete
				   the particle. */
    double dMaxStarMass;	/* maximum mass star particle to form */
    };

inline void operator|(PUP::er &p, StfmParam &param) {
    p|param.dDeltaStarForm;
    p|param.iStarFormRung;
    p|param.dGmUnit;
    p|param.dGmPerCcUnit;
    p|param.dSecUnit;
    p|param.dErgUnit;
    p|param.dPhysDenMin;
    p|param.dOverDenMin;
    p|param.dTempMax;
    p|param.dSoftMin;
    p|param.dCStar;
    p|param.dStarEff;
    p|param.dInitStarMass;
    p|param.dMinSpawnStarMass;
    p|param.dMinGasMass;
    p|param.dMaxStarMass;
    }

void StfmAddParams(StfmParam *stfm, PRM prm);
void StfmCheckParams(StfmParam *stfm, PRM prm, struct parameters &param);

#endif
