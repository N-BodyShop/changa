#ifndef FEEDBACK_HINCLUDED
#define FEEDBACK_HINCLUDED

#include "parameters.h"
#include "supernova.h"
#include "imf.h"
#include "starlifetime.h"
#include "rand.h"
#include "lymanwerner.h"
#define NFEEDBACKS 5

/**
 * @brief Class to return feedback effects.
 */
class FBEffects {
 public:
   double dEnergy;		/* Energy produced per mass (ergs/gm); note that
				   this might need to be divided up
				   into different forms of energy. */
    double dMassLoss;		/* Mass lost in Solar Masses */
    double dMetals;		/* Fraction of the mass lost in
				   elements heavier than Helium */
    double dMIron;              /* Solar masses of iron ejected */
    double dMOxygen;            /* Solar masses of oxygen ejected */
 FBEffects() : dEnergy(0), dMassLoss(0), dMetals(0), dMIron(0), dMOxygen(0) { }
 FBEffects(double dEnergy, double dMassLoss, double dMetals, double dMIron, double dMOxygen) : 
    dEnergy(dEnergy), dMassLoss(dMassLoss), dMetals(dMetals), dMIron(dMIron), dMOxygen(dMOxygen) { }
     };

/**
 * Class to Characterize Star Formation event.  This is input
 * needed for all feedback effects.
 */
class SFEvent {
 public:
    double dMass;            /* mass in star formation event in solar masses */
    double dTimeForm;      /* time of star formation event in years */
    double dMetals;           /*  metallicity of stars in event */
    double dMFracOxygen;           /*  metallicity of stars in event */
    double dMFracIron;           /*  metallicity of stars in event */
    double dLowNorm;         /* normalization constant for low mass IMF */
#ifdef STOCH24
    double rgdHMStars[24];      /* high mass stars in stochastic IMF */
    SFEvent(double mass, double tform, double mets, double fefrac, double oxfrac, double lownorm, double *hmstars) : 
    dMass(mass), dTimeForm(tform), dMetals(mets), dMFracIron(fefrac), dMFracOxygen(oxfrac), dLowNorm(lownorm) {
        for(int i=0;i<24;i++) rgdHMStars[i]=hmstars[i];
        }
#else
    double rgdHMStars[12];      /* high mass stars in stochastic IMF */
    SFEvent(double mass, double tform, double mets, double fefrac, double oxfrac, double lownorm, double *hmstars) : 
    dMass(mass), dTimeForm(tform), dMetals(mets), dMFracIron(fefrac), dMFracOxygen(oxfrac), dLowNorm(lownorm) {
        for(int i=0;i<12;i++) rgdHMStars[i]=hmstars[i];
        }
#endif

    SFEvent() : dMass(0), dTimeForm(0), dMetals(0), dMFracIron(0), dMFracOxygen(0) { }
    SFEvent(double mass, double tform, double mets, double fefrac, double oxfrac) : 
    dMass(mass), dTimeForm(tform), dMetals(mets), dMFracIron(fefrac), dMFracOxygen(oxfrac) { }
    };

/// @brief Stellar/Supernova feedback parameters and routines.
class Fdbk_Shallow {
 protected:
    double dGmUnit;		/* system mass in grams */
    double dGmPerCcUnit;	/* system density in gm/cc */
    double dErgUnit;		/* system energy in ergs */
    Padova pdva;
 public:
    char achIMF[32];	        /* initial mass function */
    mutable SN sn;
    double dDeltaStarForm;
    double dErgPerGmUnit;	/* system specific energy in ergs/gm */
    double dSecUnit;		/* system time in seconds */
    double dMaxGasMass;		/* Maximum mass of a gas particle */
#ifdef SPLITGAS
    double dInitGasMass;    /* Original mass of a gas particle */
#endif
    int bSNTurnOffCooling;      /* turn off cooling or not */
    int bSmallSNSmooth;	/* smooth SN energy only over blast radius */
    int bShortCoolShutoff;      /* use snowplow time */
    int bAGORAFeedback;         /* Replace stellar feedback with AGORA perscription */
    double dExtraCoolShutoff;      /* multiplicative factor for shutoff time */
    double dRadPreFactor;       /* McKee + Ostriker size constant in system units */
    double dTimePreFactor;      /* McKee + Ostriker time constant in system units */
    int nSmoothFeedback;	/* number of particles to smooth feedback over*/
    double dMaxCoolShutoff;     /* Maximum length of time to shutoff cooling */
    double dEarlyFeedbackFrac;  /* Fraction of SNe II to put in early feedback */
    double dFBInitialMassLoad;  /* Initial Mass Loading in Superbubble feedback*/
    double dMultiPhaseMinTemp;
    double dMultiPhaseMaxTime;
    double dEarlyETotal;  /* Total E in early FB per solar mass of stars */
    };

/// @ brief Derived class to handle the deep copy for the IMF pointer.
/// N.B.  All attributes except imf need to be in Fdbk_Shallow.
class Fdbk : public Fdbk_Shallow {
private:
    Fdbk& operator=(const Fdbk& fb); /* private to prevent use */
    void CalcWindFeedback(SFEvent *sfEvent, double dTime, 
                          double dDelta, FBEffects *fbEffects) const;
    void CalcUVFeedback(SFEvent *sfEvent, double dTime, double dDelta,
                        FBEffects *fbEffects) const;
#ifdef COOLING_MOLECULARH
    double CalcLWFeedback(SFEvent *sfEvent, double dTime, double dDelta, LWDATA *LWData) const;
#endif /*COOLING_MOLECULARH*/

public:
    IMF *imf;
    void AddParams(PRM prm);
    void CheckParams(PRM prm, struct parameters &param);
    void NullFeedback() { imf = new Kroupa01(); } /* Place holder */
    void DoFeedback(GravityParticle *p, double dTime, double dDeltaYr, 
                    FBEffects *fbTotals, LWDATA *LWData, Rand& rndGen) const;
    double NSNIa (double dMassT1, double dMassT2);

    Fdbk() { }
    Fdbk(const Fdbk& fb);
    ~Fdbk() {
	delete imf;
	}
    inline void pup(PUP::er &p);
};

//  "Deep copy" constructer is need because of imf pointer.
inline Fdbk::Fdbk(const Fdbk& fb) : Fdbk_Shallow(fb) {
    imf = fb.imf->clone();
}

inline void Fdbk::pup(PUP::er &p) {
    p | dGmUnit;
    p | dGmPerCcUnit;
    p | dErgUnit;
    p | pdva;
    p(achIMF, 32);
    p | sn;
    p | dDeltaStarForm;
    p | dErgPerGmUnit;
    p | dSecUnit;
    p | dMaxGasMass;
#ifdef SPLITGAS
    p | dInitGasMass;
#endif
    p | bSNTurnOffCooling;
    p | bSmallSNSmooth;
    p | bShortCoolShutoff;
    p | bAGORAFeedback;
    p | dExtraCoolShutoff;
    p | dRadPreFactor;
    p | dTimePreFactor;
    p | nSmoothFeedback;
    p | dMaxCoolShutoff;
    p | dEarlyFeedbackFrac;
    p | dFBInitialMassLoad;
    p | dMultiPhaseMinTemp;
    p | dMultiPhaseMaxTime;
    p | dEarlyETotal;
    p | imf;
    }

enum FBenum{
  FB_SNII=0,
  FB_SNIA,
  FB_WIND,
  FB_UV,
  FB_AGORA
};

#include "smoothparams.h"

/**
 * SmoothParams class for alerting neighboring gas particles when a star particle
 * is about to have an AGORA feedback event
 */

class AGORApreCheckSmoothParams : public SmoothParams
{
    double dTime, dDelta, H, a, gamma, etaCourant, timeToSF, dMsolUnit, dErgPerGmUnit;
    Fdbk fb;
    virtual void fcnSmooth(GravityParticle *p, int nSmooth,
               pqSmoothNode *nList);
    virtual int isSmoothActive(GravityParticle *p) { return p->isStar(); }
    virtual void initSmoothParticle(GravityParticle *p) {}
    virtual void initTreeParticle(GravityParticle *p) {}
    virtual void postTreeParticle(GravityParticle *p) {}
    virtual void initSmoothCache(GravityParticle *p);
    virtual void combSmoothCache(GravityParticle *p1,
                 ExternalSmoothParticle *p2);
    void DistFBMME(GravityParticle *p, int nSmooth, pqSmoothNode *nList);
 public:
    AGORApreCheckSmoothParams() {}
    AGORApreCheckSmoothParams(int _iType, int am, CSM csm, double _dTime, double _dDelta,
                 double _gamma, double _etaCourant, double _timeToSF, Fdbk *feedback,
                 double _dMsolUnit, double _dErgPerGmUnit) :
                    fb (*feedback) {
    iType = _iType;
    activeRung = am;
    bUseBallMax = 0;
    gamma = _gamma;
    etaCourant = _etaCourant;
    timeToSF = _timeToSF;
    dMsolUnit = _dMsolUnit;
    dErgPerGmUnit = _dErgPerGmUnit;
    dTime = _dTime;
    dDelta = _dDelta;
    if(csm->bComove) {
        H = csmTime2Hub(csm,dTime);
        a = csmTime2Exp(csm,dTime);
        }
    else {
        H = 0.0;
        a = 1.0;
        }
    }
    PUPable_decl(AGORApreCheckSmoothParams);
    AGORApreCheckSmoothParams(CkMigrateMessage *m) : SmoothParams(m) {}
    virtual void pup(PUP::er &p) {
        SmoothParams::pup(p);
    p|a;
    p|H;
    p|gamma;
    p|etaCourant;
    p|timeToSF;
    p|dMsolUnit;
    p|dErgPerGmUnit;
    p|fb;
    p|dTime;
    p|dDelta;
    }
    };

/**
 * SmoothParams class for distributing stellar feedback (energy, mass + metals) 
 * to neighboring particles.
 */

class DistStellarFeedbackSmoothParams : public SmoothParams
{
    double dTime, H, a, gamma;
    Fdbk fb;
    virtual void fcnSmooth(GravityParticle *p, int nSmooth,
			   pqSmoothNode *nList);
    virtual int isSmoothActive(GravityParticle *p);
    virtual void initSmoothParticle(GravityParticle *p) {}
    virtual void initTreeParticle(GravityParticle *p);
    virtual void postTreeParticle(GravityParticle *p);
    virtual void initSmoothCache(GravityParticle *p);
    virtual void combSmoothCache(GravityParticle *p1,
				 ExternalSmoothParticle *p2);
    void DistFBMME(GravityParticle *p, int nSmooth, pqSmoothNode *nList);
 public:
    DistStellarFeedbackSmoothParams() {}
    DistStellarFeedbackSmoothParams(int _iType, int am, CSM csm, double _dTime,
				    double _gamma, Fdbk *feedback) : 
				    fb (*feedback) {
	iType = _iType;
	activeRung = am;
	bUseBallMax = 0;
	gamma = _gamma;
	dTime = _dTime;
	if(csm->bComove) {
	    H = csmTime2Hub(csm,dTime);
	    a = csmTime2Exp(csm,dTime);
	    }
	else {
	    H = 0.0;
	    a = 1.0;
	    }
	}
        /*~DistStellarFeedbackSmoothParams() {
	delete fb;
	}*/
    PUPable_decl(DistStellarFeedbackSmoothParams);
    DistStellarFeedbackSmoothParams(CkMigrateMessage *m) : SmoothParams(m) {}
    virtual void pup(PUP::er &p) {
        SmoothParams::pup(p);//Call base class
	p|a;
	p|H;
	p|gamma;
	p|fb;
	p|dTime;
	}
    };

#endif
    
