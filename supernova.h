#ifndef SUPERNOVA_HINCLUDED
#define SUPERNOVA_HINCLUDED

#include "starlifetime.h"
#include "imf.h"

class SFEvent;
class FBEffects;

/// Methods for calculating the number and feedback effects of supernova.
class SN 
{
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
    double dESN;		/* how much energy comes from a supernova */
    int bUseStoch;   /* use stochastic IMF */
    int iNSNIIQuantum;	/* minimum amount of supernovae */
    double dFracBinSNIa;	/* fraction of binary systems in mass
				   range that go SNIa  (van den Bergh
				   & McClure, ApJ 425, 205, 1994) */
    IMF *imf;

    SN() {
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
	iNSNIIQuantum = 0;
	dFracBinSNIa = 0.05; 	/* .05 is in line with chemical evolution
				   models of the Milky Way (Francois
				   et al 2004) */
	}
    void CalcSNIIFeedback(SFEvent *sfEvent, double dTime, double dDelta, 
			  FBEffects *fbEffects);
    void CalcSNIaFeedback(SFEvent *sfEvent,double dTime, 
			  double dDelta, FBEffects *fbEffects);
    double NSNIa (double dMassT1, double dMassT2);
    friend double dMSIMFSec(SN *sn, double dMass2);
    void pup(PUP::er& p) {
	p|dESN;
    p|bUseStoch;
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
