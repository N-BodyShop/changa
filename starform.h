#ifndef STARFORM_HINCLUDED
#define STARFORM_HINCLUDED

#include "parameters.h"

class Stfm {
 private:
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
 public:
    void AddParams(PRM prm);
    void CheckParams(PRM prm, struct parameters &param);
    bool isStarFormRung(int aRung) {return aRung <= iStarFormRung;}
    GravityParticle *FormStar(GravityParticle *p,  COOL *Cool, double dTime,
			      double dDelta, double dCosmoFac);
    inline void pup(PUP::er &p);
    };

inline void Stfm::pup(PUP::er &p) {
    p|dDeltaStarForm;
    p|iStarFormRung;
    p|dGmUnit;
    p|dGmPerCcUnit;
    p|dSecUnit;
    p|dErgUnit;
    p|dPhysDenMin;
    p|dOverDenMin;
    p|dTempMax;
    p|dSoftMin;
    p|dCStar;
    p|dStarEff;
    p|dInitStarMass;
    p|dMinSpawnStarMass;
    p|dMinGasMass;
    p|dMaxStarMass;
    }

#endif
