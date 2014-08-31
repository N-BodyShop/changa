/// @file starform.h
/// Declarations for star formation.

#ifndef STARFORM_HINCLUDED
#define STARFORM_HINCLUDED

#include "parameters.h"

/// Parameters and methods to implement star formation.
class Stfm {
 private:
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
    double dMaxStarMass;	/* maximum mass star particle to form */
    int bGasCooling;		/* Can we call cooling for temperature */
 public:
    double dMinGasMass;		/* minimum mass gas before we delete
				   the particle. */
    double dDeltaStarForm;	/* timestep in system units */
    void AddParams(PRM prm);
    void CheckParams(PRM prm, struct parameters &param);
    bool isStarFormRung(int aRung) {return aRung <= iStarFormRung;}
    GravityParticle *FormStar(GravityParticle *p,  COOL *Cool, double dTime,
			      double dDelta, double dCosmoFac, double *T);
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
    p|bGasCooling;
    }

class StarLogEvent
	/* Holds statistics of the star formation event */
{
 public:
    int iOrdStar;
    int iOrdGas;
    double timeForm;
    Vector3D<double> rForm;
    Vector3D<double> vForm;
    double massForm;
    double rhoForm;
    double TForm;
 StarLogEvent() : iOrdGas(-1),	timeForm(0),rForm(0),vForm(0),
	massForm(0),rhoForm(0),TForm(0){}
    StarLogEvent(GravityParticle *p, double dCosmoFac, double TempForm) {
	iOrdGas = p->iOrder;
	// star's iOrder assigned in TreePiece::NewOrder
	timeForm = p->fTimeForm();
	rForm = p->position;
	vForm = p->velocity;
	massForm = p->fMassForm();
	rhoForm = p->fDensity/dCosmoFac;
	TForm = TempForm;
	}
    void pup(PUP::er& p) {
	p | iOrdStar;
	p | iOrdGas;
	p | timeForm;
        p | rForm;
	p | vForm;
	p | massForm;
	p | rhoForm;
	p | TForm;
	}
    };

class StarLog : public PUP::able
{
 public:
    int nOrdered;		/* The number of events that have been
				   globally ordered, incremented by
				   pkdNewOrder() */
    std::string fileName;
    std::vector<StarLogEvent> seTab;		/* The actual table */
 StarLog() : nOrdered(0),fileName("starlog") {}
    void flush();
    PUPable_decl(StarLog);
 StarLog(CkMigrateMessage *m) : PUP::able(m) {}
    void pup(PUP::er& p) {
	PUP::able::pup(p);
	p | nOrdered;
	p | fileName;
	p | seTab;
	}
    };

#endif
