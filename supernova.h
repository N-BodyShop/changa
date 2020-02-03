#ifndef SUPERNOVA_HINCLUDED
#define SUPERNOVA_HINCLUDED

#include "starlifetime.h"
#include "imf.h"

class SFEvent;
class FBEffects;

/// Methods for calculating the number and feedback effects of supernova.
class SN 
{
    friend class Fdbk;
    
    double AGORAgasLossPerSN;   /* Amount of gas (in Msun) to be ejected for one supernova */
    double AGORAmetalLossPerSN; /* Amount of metals (in Msun) to be ejected for one supernova */
    double AGORAmetalFracO;     /* Metal fraction of oxygen to be ejected during the event */
    double AGORAmetalFracFe;    /* Metal fraction of iron to be ejected during the event */

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
    Padova pdva;
 public:
    double AGORAsnTime;         /* Time (in yr) at which to set off the single SN */
    double AGORAsnE;            /* Energy (erg) to be released from a single SN */
    double AGORAsnPerMass;      /* Number of supernovae to set off per solar mass of star particle */
    double dESN;		/* how much energy comes from a supernova */
    int bUseStoch;   /* use stochastic IMF */
    double dStochCut;
    int iNSNIIQuantum;	/* minimum amount of supernovae */
    double dFracBinSNIa;	/* fraction of binary systems in mass
				   range that go SNIa  (van den Bergh
				   & McClure, ApJ 425, 205, 1994) */
    IMF *imf;

    SN() {
	/* parameters for AGORA feedback
	   These are specified in paper 4, dataset 2 */
        AGORAsnTime = 5.e6;
        AGORAsnPerMass = 91.;
        AGORAgasLossPerSN = 14.8;
        AGORAmetalLossPerSN = 2.63;
        AGORAmetalFracO = 0.098;
        AGORAmetalFracFe = 0.43;
        AGORAsnE = 1.0e51;

	dMSNrem = 1.4;		/* mass of supernova remnant in solar masses 
                                 * Also used for SNIa ejected mass */
	dMSNIImin = 8.0;	/* Mass above which stars supernova in solar
				   masses */
	dMSNIImax = 40.;	/* Mass below which stars supernova in
				   solar masses */
	dMBmin = 3.0;		/* Minimum mass of binary that can go SNIa */
	dMBmax = 16.0;		/* Maximum mass of binary that can go SNIa */
  /* normalization constant and exponent in formulae for masses of
     ejected Fe and O16 as a function of stellar mass taken from
     Raiteri, Villata and Navarro, A&A 315, 105, 1996 */
	dMEjexp = 1.056;    
	dMEjconst = 0.7682; 
	dMFeexp = 1.864;    
	dMFeconst = 2.802e-4;
	dMOxexp = 2.721;
	dMOxconst = 4.586e-4; 
	dESN = 0.1e51;
    bUseStoch = 1;
    dStochCut= 8.0;
	iNSNIIQuantum = 0;
	dFracBinSNIa = 0.05; 	/* .05 is in line with chemical evolution
				   models of the Milky Way (Francois
				   et al 2004) */
	}
    void CalcAGORAFeedback(SFEvent *sfEvent, double dTime, double dDelta,
                           FBEffects *fbEffects) const;
    void CalcSNIIFeedback(SFEvent *sfEvent, double dTime, double dDelta, 
                          FBEffects *fbEffects) const;
    void CalcSNIaFeedback(SFEvent *sfEvent,double dTime, 
                          double dDelta, FBEffects *fbEffects) const;
    double NSNIa (double dMassT1, double dMassT2) const;
    friend double dMSIMFSec(const SN *sn, double dMass2);
    void pup(PUP::er& p) {
	p|dESN;
    p|bUseStoch;
    p|dStochCut;
	p|dMSNrem;
	p|dMSNIImin;
	p|dMSNIImax;
	p|dMBmin;
	p|dMBmax;
	p|dMEjexp;
	p|dMEjconst;
	p|dMFeexp;
	p|dMFeconst;
	p|dMOxexp;
	p|dMOxconst;
	p|dFracBinSNIa;
	p|iNSNIIQuantum;
	p|pdva;
	// p|*imf;  Think about this
    };
    
};

#endif
