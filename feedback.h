#ifndef FEEDBACK_HINCLUDED
#define FEEDBACK_HINCLUDED

#include "parameters.h"
#include "imf.h"
#include "starlifetime.h"
#define NFEEDBACKS 4

/*
 * Class to return feedback effects.
 */
class FBEffects {
 public:
   double dEnergy;		/* Energy produced (ergs); note that
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

/*
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
 SFEvent() : dMass(0), dTimeForm(0), dMetals(0), dMFracIron(0), dMFracOxygen(0) { }
 SFEvent(double mass, double tform, double mets, double fefrac, double oxfrac) : 
    dMass(mass), dTimeForm(tform), dMetals(mets), dMFracIron(fefrac), dMFracOxygen(oxfrac) { }
    };

class Fdbk : public PUP::able {
 private:
    Fdbk& operator=(Fdbk& fb) {}
    void CalcSNIIFeedback(SFEvent *sfEvent, double dTime, double dDelta, 
			  FBEffects *fbEffects);
    void CalcSNIaFeedback(SFEvent *sfEvent,double dTime, 
			  double dDelta, FBEffects *fbEffects);
    double NSNIa (double dMassT1, double dMassT2);
    //    double dMSIMFSec(Padova *pdva, double dMass2);
    double MSNIa (double dMassT1, double dMassT2);
    //    double MSIMFSecM(Padova *pdva, double dMass2);
    void CalcWindFeedback(SFEvent *sfEvent, double dTime, 
			  double dDelta, FBEffects *fbEffects);
    void CalcUVFeedback(double dTime, double dDelta, FBEffects *fbEffects);

    std::string achIMF;	        /* initial mass function */
    int iRandomSeed;		/* seed for stochastic quantized feedback */
    double dErgPerGmUnit;	/* system specific energy in ergs/gm */
    double dGmUnit;		/* system mass in grams */
    double dGmPerCcUnit;	/* system density in gm/cc */
    double dErgUnit;		/* system energy in ergs */
    double dInitStarMass; 
    double dMSNrem;		/* mass of SN remnant in M_sun */
    double dMSNIImin;		/* Minimum core collapse SN mass */
    double dMSNIImax;           /* Maximum core collapse SN mass */
    double dMBmin;              /* Minimum Mass of binary that can do SNIa */
    double dMBmax;		/* Maximum mass of binary that an SNIa */
    double dMEjexp;             /* exponent of ejection mass power law */
    double dMEjconst;           /* normalization for ejection mass power law */
    double dMFeexp;             /* exponent of iron productions */
    double dMFeconst;           /* normalization of iron */
    double dMOxexp;             /* exponent of oxygen */
    double dMOxconst;           /* normalization of oxygen */
    double dSNIaMetals;         /* amount of metals produced by SNIa
				   */
    Padova pdva;
 public:
    double dDeltaStarForm;
    double dSecUnit;		/* system time in seconds */
    int iNSNIIQuantum;	/* minimum amount of supernovae */
    double dESN;		/* how much energy comes from a supernova */
    int bSNTurnOffCooling;      /* turn off cooling or not */
    int bSmallSNSmooth;	/* smooth SN energy only over blast radius */
    int bShortCoolShutoff;      /* use snowplow time */
    int dExtraCoolShutoff;      /* multiplicative factor for shutoff time */
    double dRadPreFactor;       /* McKee + Ostriker size constant in system units */
    double dTimePreFactor;      /* McKee + Ostriker time constant in system units */
    int nSmoothFeedback;	/* number of particles to smooth feedback over*/
    double dFracBinSNIa;	/* fraction of binary systems in mass range that 
				   go SNIa  (van den Bergh & McClure, ApJ 425, 205, 1994) */
    IMF *imf;

    void AddParams(PRM prm);
    void CheckParams(PRM prm, struct parameters &param);
    void DoFeedback(GravityParticle *p, double dTime, double dDeltaYr, 
		    FBEffects *fbTotals);
 Fdbk() : achIMF("Kroupa93"),iRandomSeed(1),iNSNIIQuantum(0.),dESN(0.1e51),
	bSmallSNSmooth(1),bShortCoolShutoff(0),dExtraCoolShutoff(0),
	bSNTurnOffCooling(1), nSmoothFeedback(64),
	/* supernova constants */
	dMSNrem(1.4),	/* mass of supernova remnant in solar masses 
			 * Also used for SNIa ejected mass */
	dMSNIImin(8.0),/* Mass above which stars supernova in solar
			  masses */
	dMSNIImax(40.),/* Mass below which stars supernova in solar masses */
	dMBmin(0.8),	/* Minimum mass of binary that can go SNIa */
	
	dMBmax(16.0),	/* Maximum mass of binary that can go SNIa */
	dFracBinSNIa(0.16),/*0.04847),	 fraction of binary systems in
			     appropriate mass range that go SNIa =
			     0.16 (van den Bergh & McClure, ApJ
			     425, 205, 1994) */
	/* normalization constant and exponent in formulae for masses of
	   ejected Fe and O16 as a function of stellar mass taken from
	   Raiteri, Villata and Navarro, A&A 315, 105, 1996 */
	dMEjexp(1.056),dMEjconst(0.7682), dMFeexp(1.864), dMFeconst(2.802e-4),
	dMOxexp(2.721),	dMOxconst(4.586e-4), 
	dSNIaMetals(1.4)  /* Ia's ejecta are entirely metals */
	    { /*pdva = new Padova();*/}

 Fdbk(const Fdbk& fb) : achIMF (fb.achIMF), iRandomSeed(fb.iRandomSeed),
	iNSNIIQuantum (fb.iNSNIIQuantum ), dESN (fb.dESN),
	bSmallSNSmooth (fb.bSmallSNSmooth), 
	bShortCoolShutoff (fb.bShortCoolShutoff),
	dExtraCoolShutoff (fb.dExtraCoolShutoff),
	bSNTurnOffCooling (fb.bSNTurnOffCooling),
	nSmoothFeedback (fb.nSmoothFeedback),
	dMSNrem (fb.dMSNrem),dMSNIImin (fb.dMSNIImin),dMSNIImax (fb.dMSNIImax),
	dMBmin (fb.dMBmin),dMBmax (fb.dMBmax), dFracBinSNIa (fb.dFracBinSNIa),
	dMEjexp (fb.dMEjexp),dMEjconst (fb.dMEjconst), dMFeexp (fb.dMFeexp),
	dMFeconst (fb.dMFeconst), dMOxexp (fb.dMOxexp),dMOxconst (fb.dMOxconst),
	dSNIaMetals (fb.dSNIaMetals), dSecUnit(fb.dSecUnit),
	dGmPerCcUnit(fb.dGmPerCcUnit), dGmUnit(fb.dGmUnit),
	dErgUnit(fb.dErgUnit), dErgPerGmUnit(fb.dErgPerGmUnit),
	dRadPreFactor(fb.dRadPreFactor),dTimePreFactor(fb.dTimePreFactor),
	pdva(fb.pdva) {	imf = fb.imf->clone();	}

    PUPable_decl(Fdbk);
 Fdbk(CkMigrateMessage *m) : PUP::able(m) {}
    ~Fdbk() {
	delete imf;
	}
    inline void pup(PUP::er &p);
    };

inline void Fdbk::pup(PUP::er &p) {
    p | achIMF;
    p | dDeltaStarForm;
    p | iRandomSeed;
    p | dErgPerGmUnit;
    p | dGmUnit;
    p | dGmPerCcUnit;
    p | dErgUnit;
    p | dInitStarMass; 
    p | dMSNrem;
    p | dMSNIImin;
    p | dMSNIImax;
    p | dMBmin;
    p | dMBmax;
    p | dMEjexp;
    p | dMEjconst;
    p | dMFeexp;
    p | dMFeconst;
    p | dMOxexp;
    p | dMOxconst;
    p | dSNIaMetals;
    p | dSecUnit;
    p | iNSNIIQuantum;
    p | dESN;
    p | bSNTurnOffCooling;
    p | bSmallSNSmooth;
    p | bShortCoolShutoff;
    p | dExtraCoolShutoff;
    p | dRadPreFactor;
    p | dTimePreFactor;
    p | nSmoothFeedback;
    p | dFracBinSNIa;
    p | imf;
    }
PUPmarshall(Fdbk);

enum FBenum{
  FB_SNII=0,
  FB_SNIA,
  FB_WIND,
  FB_UV
};

	/* solar mass in grams */
#define MSOLG 1.99e33
        /* seconds per year = 3600*24*365.25 */
#define SEC_YR 3.15576e07;

#endif

    
