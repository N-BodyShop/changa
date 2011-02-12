/*
 * Stellar feedback module originally written for GASOLINE
 * Modified for inclusion in ChaNGa
 *
 * Algorithm was originally implemented in TREESPH but is heavily
 * modified.  Contributers include Neal Katz, Eric Hayashi, Greg Stinson
 */
#include <math.h>
#include "ParallelGravity.h"
#include "feedback.h"
#include "smooth.h"
#include "Sph.h"


///
/// @brief initialize parameters for star formation
///

void Fdbk::AddParams(PRM prm)
{
        achIMF = "Kroupa93";
	prmAddParam(prm,"achIMF", 1, &achIMF, sizeof(char [32]), "achIMF", "<IMF> = Kroupa93");
	iRandomSeed = 1;
	prmAddParam(prm,"iRandomSeed", 1, &iRandomSeed, sizeof(int), "iRand",
		    "<Feedback random Seed> = 1");
	iNSNIIQuantum = 0.;
	prmAddParam(prm,"iNSNIIQuantum", 1, &iNSNIIQuantum, sizeof(int), "snQuant",
		    "<Min # SNII per timestep> = 0.");
	dESN = 0.1e51;
	prmAddParam(prm,"dESN", 2, &dESN, sizeof(double), "snESN",
		    "<Energy of supernova in ergs> = 0.1e51");
	if (dESN > 0.0) bSmallSNSmooth = 1;
        else bSmallSNSmooth = 0;
	prmAddParam(prm,"bSmallSNSmooth", 0, &bSmallSNSmooth, sizeof(int), "bSmallSNSmooth",
		    "<smooth SN ejecta over blast or smoothing radius> = blast radius");
        bShortCoolShutoff = 0;
	prmAddParam(prm,"bShortCoolShutoff", 0, &bShortCoolShutoff,
		    sizeof(int), "bShortCoolShutoff", "<Use snowplow time> = 0");
        dExtraCoolShutoff = 0;
	prmAddParam(prm,"dExtraCoolShutoff", 0, &dExtraCoolShutoff,
		    sizeof(int), "dExtraCoolShutoff", "<Use snowplow time> = 0");
        bSNTurnOffCooling = 1;
	prmAddParam(prm,"bSNTurnOffCooling", 0, &bSNTurnOffCooling,
		    sizeof(int), "bSNTurnOffCooling", "<Do SN turn off cooling> = 1");
	nSmoothFeedback = 64;
	prmAddParam(prm,"nSmoothFeedback",1,&nSmoothFeedback,sizeof(int),"s",
		    "<number of particles to smooth feedback over> = 64");
/* supernova constants */
	
	dMSNrem = 1.4;			/* mass of supernova remnant in solar masses 
						 * Also used for SNIa ejected mass */
    
	dMSNIImin = 8.0;		/* Mass above which stars
					   supernova in solar
					   masses */
	dMSNIImax = 40.;		/* Mass below which stars
					   supernova in solar masses */
	dMBmin = 3.0;			/* Minimum mass of binary that
						   can go SNIa */
	
	dMBmax = 16.0;			/* Maximum mass of binary that
						   can go SNIa */
	dFracBinSNIa = 0.16;/*0.04847;	 fraction of binary systems in
				  appropriate mass range that go SNIa =
				  0.16 (van den Bergh & McClure, ApJ
				  425, 205, 1994) */
	/* normalization constant and exponent in formulae for masses of
	   ejected Fe and O16 as a function of stellar mass taken from
	   Raiteri, Villata and Navarro, A&A 315, 105, 1996 */
	dMEjexp = 1.056;    
	dMEjconst = 0.7682; 
	dMFeexp = 1.864;    
	dMFeconst = 2.802e-4;
	dMOxexp = 2.721;
	dMOxconst = 4.586e-4; 
	dSNIaMetals = 1.4;  /* Ia's ejecta are entirely metals */
	if(achIMF == "MillerScalo") imf = new MillerScalo();
	else if(achIMF == "Chabrier") imf = new Chabrier();
	else imf = new Kroupa93();
	//	pdva = new Padova();
}

void Fdbk::CheckParams(PRM prm, struct parameters &param)
{
#include "physconst.h"

	dSecUnit = param.dSecUnit;
	dGmPerCcUnit = param.dGmPerCcUnit;
	dGmUnit = param.dMsolUnit*MSOLG;
	dErgUnit = GCGS*pow(param.dMsolUnit*MSOLG, 2.0)
	  /(param.dKpcUnit*KPCCM);
	dErgPerGmUnit = GCGS*param.dMsolUnit*MSOLG/(param.dKpcUnit*KPCCM);
	/* eq 3 from McCray + Kafatos (1987) is good both before SN start
	 * exploding and after because Luminosity is essentially constant
	 * Stellar winds:  dM/dt ~ 1e-6 M_sun / yr, v ~ 2000 km/s 
	 * (cf Castor et al (1975) for theory)
	 */
	if (iNSNIIQuantum > 0) {
	  dRadPreFactor = (0.097 / param.dKpcUnit) *
	    pow(param.dGmPerCcUnit/MHYDR,-0.2)*
	    pow(param.dSecUnit/SECONDSPERYEAR / 1e7,0.6);
	  /* Solar metallicity Z from Grevesse & Sauval (1998) SSR 85 161*/
	  /* TOO LONG dTimePreFactor = 4e6*SECONDSPERYEAR/
	    param.dSecUnit * pow(0.018,1.5)*
	    pow(param.dGmPerCcUnit/MHYDR,-0.7);*/
	} else {
	  /* from McKee and Ostriker (1977) ApJ 218 148 */
	  dRadPreFactor = pow(10,1.74)/(param.dKpcUnit*1000.0)*
	    pow(MSOLG*param.dMsolUnit/(param.dMeanMolWeight*MHYDR*pow(KPCCM*param.dKpcUnit,3)),-0.16)*
	    pow(0.0001*GCGS*pow(MSOLG*param.dMsolUnit,2)/(pow(KPCCM*param.dKpcUnit,4)*KBOLTZ),-0.2);
	}
	if (bShortCoolShutoff){        /* end of snowplow */
	  dTimePreFactor = 
	    SECONDSPERYEAR*pow(10,5.92)/(param.dSecUnit)*
	    pow(MSOLG*param.dMsolUnit/(param.dMeanMolWeight*MHYDR*pow(KPCCM*param.dKpcUnit,3)),0.27)*
	    pow(0.0001*GCGS*pow(MSOLG*param.dMsolUnit,2)/(pow(KPCCM*param.dKpcUnit,4)*KBOLTZ),-0.64);
	} else {       /* t_{max}*/
	  dTimePreFactor = dExtraCoolShutoff*SECONDSPERYEAR*pow(10,6.85)/(dSecUnit)*
	    pow(MSOLG*param.dMsolUnit/(param.dMeanMolWeight*MHYDR*pow(KPCCM*param.dKpcUnit,3)),0.32)*
	    pow(0.0001*GCGS*pow(MSOLG*param.dMsolUnit,2)/(pow(KPCCM*param.dKpcUnit,4)*KBOLTZ),-0.70);
	}

}

///
/// form stars main method
///
void Main::StellarFeedback(double dTime, double dDelta) 
{
    int iPhase = 0; // Keeps track of node cache use
    
    if(verbosity)
	CkPrintf("Stellar Feedback ... \n");
    double startTime = CkWallTimer();

    CkReductionMsg *msgFeedback;
    treeProxy.Feedback(*(param.feedback), dTime, dDelta,
			CkCallbackResumeThread((void*&)msgFeedback));
    double *dFeedback = (double *)msgFeedback->getData();
    
    if(verbosity) 
      {
	printf("Feedback totals: mass, energy, metalicity\n");
	for(int i = 0; i < NFEEDBACKS; i++){
	  CkPrintf("feedback %d: %g %g %g\n", i,
		   dFeedback[i*3],
		   dFeedback[i*3 + 1],
		   dFeedback[i*3] != 0.0 ?
		   dFeedback[i*3 + 2]/dFeedback[i*3] : 0.0);
	}
      }
    delete msgFeedback;
    
    if(verbosity)
      CkPrintf("Distribute Stellar Feedback ... ");
    // Need to build tree since we just did addDelParticle.
    // XXX need to check whether a treebuild needs the domain
    // decomposition.  If not, this could be avoided.
    //
    double tolerance = 0.01;	// tolerance for domain decomposition
    sorter.startSorting(dataManagerID, tolerance,
                        CkCallbackResumeThread(), true);
    treeProxy.buildTree(bucketSize, CkCallbackResumeThread());
    DistStellarFeedbackSmoothParams pDSFB(TYPE_STAR, 0, param.csm, dTime, 
					  param.dConstGamma, param.feedback);
    treeProxy.startIterationReSmooth(&pDSFB, CkCallbackResumeThread());
    iPhase++;

    CkAssert(iPhase <= nPhases);
    if(iPhase < nPhases)
	treeProxy.finishNodeCache(nPhases-iPhase, CkCallbackResumeThread());

    CkPrintf("Stellar Feedback Calculated, Wallclock %f secs\n",
	     CkWallTimer() - startTime);
    }

///
/// processor specific method for star formation
/// 
void TreePiece::Feedback(Fdbk fb, double dTime, double dDelta, const CkCallback& cb)
{
    FBEffects fbTotals[NFEEDBACKS];
    double dDeltaYr;
    double dSNIaMassStore;
        
    dTime *= fb.dSecUnit/SEC_YR;
    dDeltaYr = dDelta*fb.dSecUnit/SEC_YR;
    
    /* 0 out accumulator class before starting feedback */
    for(int i = 0; i < NFEEDBACKS; i++) {
	fbTotals[i].dMassLoss = 0.0;
	fbTotals[i].dEnergy = 0.0;
	fbTotals[i].dMetals = 0.0;
	fbTotals[i].dMIron = 0.0;
	fbTotals[i].dMOxygen = 0.0;
	}

    /* loop through particles */
    for(unsigned int i = 1; i <= myNumParticles; ++i) {
	GravityParticle *p = &myParticles[i];
	if(p->isStar() && p->fTimeForm() >= 0.0) 
	  fb.DoFeedback(p, dTime, dDeltaYr, fbTotals);
	else if(p->isGas()){
	  assert(p->u() >= 0.0);
	  assert(p->uPred() >= 0.0);
	  p->fESNrate() = 0.0;	/* reset SN heating rate of gas to zero */
	}
    }

    double dFeedback[NFEEDBACKS*3];
    for(int i = 0; i < NFEEDBACKS; i++) {
	dFeedback[i*3] = fbTotals[i].dMassLoss;
	dFeedback[i*3 + 1] = fbTotals[i].dEnergy;
	dFeedback[i*3 + 2] = fbTotals[i].dMetals;
    }
    contribute(sizeof(dFeedback),dFeedback, CkReduction::sum_double, cb);
}

void Fdbk::DoFeedback(GravityParticle *p, double dTime, double dDeltaYr, 
		      FBEffects *fbTotals)
{
  double dTotMassLoss, dTotMetals, dTotMOxygen, dTotMIron, dDelta;
  dTotMassLoss = dTotMetals = dTotMOxygen = dTotMIron = 0.0;
  dDelta = dDeltaYr*SECONDSPERYEAR / dSecUnit;
  p->fESNrate() = 0.0;
  p->fNSN() = 0.0;
  
    FBEffects fbEffects;
    // Particle properties that will be sent to feedback
    // methods in normal units (M_sun + seconds)
    SFEvent sfEvent(p->fMassForm()*dGmUnit/MSOLG, 
		    p->fTimeForm()*dSecUnit/SECONDSPERYEAR,
		    p->fStarMetals(), p->fStarMFracOxygen(), 
		    p->fStarMFracIron());

    double dSNIaMassStore=0.0;  /* Stores mass loss of Ia so as not to 
				   double count it in wind feedback */
                                    
    /*
     * Call all the effects in order and accumulate them.
     */
    for(int j = 0; j < NFEEDBACKS; j++) {
      double dNSNII = 0;
      switch (j) {
      case FB_SNII:
	CalcSNIIFeedback(&sfEvent, dTime, dDeltaYr, &fbEffects);
	if (dESN > 0) p->fNSN() = fbEffects.dEnergy / dESN;
	break;
      case FB_SNIA:
	CalcSNIaFeedback(&sfEvent, dTime, dDeltaYr, &fbEffects);
	dSNIaMassStore=fbEffects.dMassLoss;
	break;
      case FB_WIND:
	CalcWindFeedback(&sfEvent, dTime, dDeltaYr, &fbEffects);
	fbEffects.dMassLoss -= dSNIaMassStore;
	/*printf("Wind, SNaI Mass Loss: %d   %d\n",fbEffects.dMassLoss,dSNIaMassStore); */
	break;
      case FB_UV:
	CalcUVFeedback(dTime, dDeltaYr, &fbEffects);
	break;
      default:
	assert(0);
      }
      
      // Convert to system units
      fbEffects.dMassLoss *= MSOLG/dGmUnit;
      fbEffects.dEnergy /= dErgPerGmUnit;
      
      dTotMassLoss += fbEffects.dMassLoss;
      p->fESNrate() += fbEffects.dEnergy*fbEffects.dMassLoss;
      dTotMetals += fbEffects.dMetals*fbEffects.dMassLoss;
      dTotMOxygen += fbEffects.dMOxygen*fbEffects.dMassLoss;
      dTotMIron += fbEffects.dMIron*fbEffects.dMassLoss;
      
      fbTotals[j].dMassLoss += fbEffects.dMassLoss;
      fbTotals[j].dEnergy += fbEffects.dEnergy*fbEffects.dMassLoss;
      fbTotals[j].dMetals += fbEffects.dMetals*fbEffects.dMassLoss;
      fbTotals[j].dMIron += fbEffects.dMIron*fbEffects.dMassLoss;
      fbTotals[j].dMOxygen += fbEffects.dMOxygen*fbEffects.dMassLoss;
    }

    /*
     * Modify star particle
     */
    assert(p->mass > dTotMassLoss);
    
    p->mass -= dTotMassLoss;
    p->fMSN() = dTotMassLoss;
    /* The SNMetals and ESNrate used to be specific
       quantities, but we are making them totals as
       they leave the stars so that they are easier
       to divvy up among the gas particles in 
       distSNEnergy in smoothfcn.c.  These quantities
       will be converted back to specific quantities when
       they are parts of gas particles. */
    p->fSNMetals() = dTotMetals;
    p->fMIronOut() = dTotMIron;
    p->fMOxygenOut() = dTotMOxygen;
    p->fESNrate() /= dDelta; /* convert to rate */
}

void Fdbk::CalcWindFeedback(SFEvent *sfEvent, double dTime, /* current time in years */
			    double dDelta, /* length of timestep (years) */
			    FBEffects *fbEffects)
{
    double dMDying;

    /* First determine if dying stars are between 1-8 Msolar
     * stellar lifetimes corresponding to beginning and end of 
     * current timestep with respect to starbirth time in yrs */
    double dMmin=1.0;
    double dMmax=8.0;
    double dMinLifetime = dTime - sfEvent->dTimeForm; 
    double dMaxLifetime = dMinLifetime + dDelta;

    double dMinMass = pdva.StarMass(dMaxLifetime, sfEvent->dMetals); 
    double dMaxMass = pdva.StarMass(dMinLifetime, sfEvent->dMetals); 
    assert(dMaxMass >= dMinMass);

    if (((dMinMass < dMmax) && (dMaxMass > dMmin)) && dMaxMass > dMinMass) {

	/* Mass Fraction returned to ISM taken from Weidemann, 1987, A&A 188 74 
	   then fit to function: MFreturned = 0.86 - exp(-Mass/1.1) */
	double dMassFracReturned=0.86-exp(-((dMaxMass+dMinMass)/2.)/1.1);

	double dMinCumMass = imf->CumMass(dMinMass);
	double dMaxCumMass = imf->CumMass(dMaxMass);
	double dMTot = imf->CumMass(0.0);
	/* Find out mass fraction of dying stars, then multiply by the original
	   mass of the star particle */
	if (dMTot == 0.0){
	  dMDying = 0.0;
	} else { 
	  dMDying = (dMinCumMass - dMaxCumMass)/dMTot;
	}
	dMDying *= sfEvent->dMass;

	/* Figure out feedback effects */
	fbEffects->dMassLoss = dMDying * dMassFracReturned;
	fbEffects->dEnergy = 0.0;    
	/* Use star's metallicity for gas returned */
	fbEffects->dMetals = sfEvent->dMetals; 
	fbEffects->dMIron = sfEvent->dMFracIron; 
	fbEffects->dMOxygen = sfEvent->dMFracOxygen; 

	} else {
	    fbEffects->dMassLoss = 0.0;
	    fbEffects->dEnergy = 0.0;    
	    fbEffects->dMetals = 0.0; 
	    fbEffects->dMIron = 0.0;
	    fbEffects->dMOxygen = 0.0;
	    }
    }

void Fdbk::CalcUVFeedback(double dTime, /* current time in years */
			  double dDelta, /* length of timestep (years) */
			  FBEffects *fbEffects)
{
    fbEffects->dMassLoss = 0.0;
    fbEffects->dEnergy = 0.0;
    fbEffects->dMetals = 0.0;
    fbEffects->dMIron = 0.0;
    fbEffects->dMOxygen = 0.0;
    }

