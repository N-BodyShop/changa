/*
 * Stellar feedback module originally written for GASOLINE
 * Modified for inclusion in ChaNGa
 *
 * Algorithm was originally implemented in TREESPH but is heavily
 * modified.  Contributers include Neal Katz, Eric Hayashi, Greg Stinson
 */
#include <math.h>
#ifdef HAVE_VALUES_H
#include <values.h>
#else
#include <float.h>
#endif

#include "ParallelGravity.h"
#include "feedback.h"
#include "smooth.h"
#include "Sph.h"

#define max(A,B) ((A) > (B) ? (A) : (B))
///
/// @brief initialize parameters for star formation
///

void Fdbk::AddParams(PRM prm)
{
    strcpy(achIMF, "Kroupa93");
    prmAddParam(prm,"achIMF", paramString, &achIMF, 32, "achIMF",
		"<IMF> = Kroupa93");
    sn.iNSNIIQuantum = 0;
    prmAddParam(prm,"iNSNIIQuantum", paramInt, &sn.iNSNIIQuantum, sizeof(int),
		"snQuant", "<Min # SNII per timestep> = 0.");
    sn.dESN = 0.1e51;
    prmAddParam(prm,"dESN", paramDouble, &sn.dESN, sizeof(double), "snESN",
		    "<Energy of supernova in ergs> = 0.1e51");
    bSmallSNSmooth = 1;
    prmAddParam(prm,"bSmallSNSmooth", paramBool, &bSmallSNSmooth, sizeof(int),
		"bSmallSNSmooth",
		"<smooth SN ejecta over blast or smoothing radius> = blast radius");
    bShortCoolShutoff = 0;
    prmAddParam(prm,"bShortCoolShutoff", paramBool, &bShortCoolShutoff,
		    sizeof(int), "bShortCoolShutoff", "<Use snowplow time> = 0");
    dExtraCoolShutoff = 1.0;
    prmAddParam(prm,"dExtraCoolShutoff", paramDouble, &dExtraCoolShutoff,
		sizeof(double), "dExtraCoolShutoff", "<Extend shutoff time> = 1.0");
    dMaxGasMass = FLT_MAX;
    prmAddParam(prm,"dMaxGasMass", paramDouble, &dMaxGasMass,
		sizeof(double), "stMaxGas",
		"<Maximum mass of a gas particle> = FLT_MAX");
    bSNTurnOffCooling = 1;
    prmAddParam(prm,"bSNTurnOffCooling", paramBool, &bSNTurnOffCooling,
		sizeof(int), "bSNTurnOffCooling", "<Do SN turn off cooling> = 1");
    nSmoothFeedback = 64;
    prmAddParam(prm,"nSmoothFeedback", paramInt,&nSmoothFeedback, sizeof(int),
		"s", "<number of particles to smooth feedback over> = 64");
    dMaxCoolShutoff = 3.0e7;
    prmAddParam(prm,"dMaxCoolShutoff", paramDouble, &dMaxCoolShutoff,
		sizeof(double), "fbMaxCoolOff",
		"<Maximum time to shutoff cooling in years> = 3e7");
    bAGORAFeedback = 0;
    prmAddParam(prm, "bAGORAFeedback", paramBool, &bAGORAFeedback, sizeof(int),
            "bAGORAFeedback", "<Replace Type II and Ia supernovae with AGORA perscription> = 0");
}

void Fdbk::CheckParams(PRM prm, struct parameters &param)
{
    if(strcmp(achIMF, "MillerScalo") == 0) imf = new MillerScalo();
    else if(strcmp(achIMF, "Chabrier") == 0) imf = new Chabrier();
    else if(strcmp(achIMF, "Kroupa93") == 0) imf = new Kroupa93();
    else if(strcmp(achIMF, "Kroupa01") == 0) imf = new Kroupa01();
    sn.imf = imf;

#include "physconst.h"
    if(!prmSpecified(prm, "nSmoothFeedback"))
	nSmoothFeedback = param.nSmooth;
    if (sn.dESN > 0.0) bSmallSNSmooth = 1;
    else bSmallSNSmooth = 0;
    param.bDoGas = 1;
    dDeltaStarForm = param.stfm->dDeltaStarForm;
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
    if (sn.iNSNIIQuantum > 0) {
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
    dMaxCoolShutoff *= SECONDSPERYEAR/param.dSecUnit;
    }

///
/// Feedback main method
///
void Main::StellarFeedback(double dTime, double dDelta) 
{
    if(verbosity)
	CkPrintf("Stellar Feedback ... \n");
    double startTime = CkWallTimer();

    CkReductionMsg *msgFeedback;
    treeProxy.Feedback(*(param.feedback), dTime, dDelta,
		       CkCallbackResumeThread((void*&)msgFeedback));
    double *dFeedback = (double *)msgFeedback->getData();
    
    if(verbosity) 
      {
	CkPrintf("Feedback totals: mass, energy, metalicity\n");
	for(int i = 0; i < NFEEDBACKS; i++){
	  CkPrintf("feedback %d: %g %g %g\n", i,
		   dFeedback[i*3],
		   dFeedback[i*3 + 1],
		   dFeedback[i*3] != 0.0 ?
		   dFeedback[i*3 + 2]/dFeedback[i*3] : 0.0);
	}
      }
    delete msgFeedback;
    CkReductionMsg *msgChk;
    treeProxy.massMetalsEnergyCheck(1, CkCallbackResumeThread((void*&)msgChk));
    
    if(verbosity)
      CkPrintf("Distribute Stellar Feedback ... ");
    // Need to build tree since we just did addDelParticle.
    //
    treeProxy.buildTree(bucketSize, CkCallbackResumeThread());
    DistStellarFeedbackSmoothParams pDSFB(TYPE_GAS, 0, param.csm, dTime, 
					  param.dConstGamma, param.feedback);
    double dfBall2OverSoft2 = 4.0*param.dhMinOverSoft*param.dhMinOverSoft;
    treeProxy.startSmooth(&pDSFB, 0, param.feedback->nSmoothFeedback,
			  dfBall2OverSoft2, CkCallbackResumeThread());
    treeProxy.finishNodeCache(CkCallbackResumeThread());

    CkPrintf("Stellar Feedback Calculated, Wallclock %f secs\n",
	     CkWallTimer() - startTime);

    CkReductionMsg *msgChk2;
    treeProxy.massMetalsEnergyCheck(0, CkCallbackResumeThread((void*&)msgChk2));
    double *dTotals = (double *)msgChk->getData();
    double *dTotals2 = (double *)msgChk2->getData();
    int i;
    for(i = 0; i < 5; i++) {
	std::string labels[5] = {"Mass", "Metals", "Oxygen", "Iron", "Energy"};
	if(verbosity > 1)
	    CkPrintf("Total %s: %g\n", labels[i].c_str(), dTotals[i]);

    // Supress energy conservation warning if AGORA feedback is enabled
    if (!param.feedback->bAGORAFeedback) {
        if(fabs(dTotals[i] - dTotals2[i]) > 1e-12*(dTotals[i])) {
            CkError("ERROR: %s not conserved: %.15e != %.15e!\n",
                labels[i].c_str(), dTotals[i], dTotals2[i]);
            }
        }
    }

    delete msgChk;
    delete msgChk2;
    }

///
/// processor specific method for stellar feedback
/// 
void TreePiece::Feedback(Fdbk &fb, double dTime, double dDelta, const CkCallback& cb)
{
    FBEffects fbTotals[NFEEDBACKS];
    double dDeltaYr;
        
    dTime *= fb.dSecUnit/SECONDSPERYEAR ;
    dDeltaYr = max(dDelta,fb.dDeltaStarForm)*fb.dSecUnit/SECONDSPERYEAR ;
    fb.sn.imf = fb.imf;		// point sn imf at our imf
    
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
	  CkAssert(p->u() >= 0.0);
	  CkAssert(p->uPred() >= 0.0);
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

/// @brief Fdbk main method.
/// @param p Star particle doing feedback
/// @param dTime Current time in years.
/// @param dDeltaYr Timestep in years.
/// @param fbTotals pointer to total effects for bookkeeping

void Fdbk::DoFeedback(GravityParticle *p, double dTime, double dDeltaYr, 
		      FBEffects *fbTotals)
{
    double dTotMassLoss, dTotMetals, dTotMOxygen, dTotMIron, dDelta;
    dTotMassLoss = dTotMetals = dTotMOxygen = dTotMIron = 0.0;
    dDelta = dDeltaYr*SECONDSPERYEAR / dSecUnit;
    p->fStarESNrate() = 0.0;
    p->fNSN() = 0.0;
  
    FBEffects fbEffects;
    // Particle properties that will be sent to feedback
    // methods in normal units (M_sun + seconds)
    SFEvent sfEvent(p->fMassForm()*dGmUnit/MSOLG, 
		    p->fTimeForm()*dSecUnit/SECONDSPERYEAR,
		    p->fStarMetals(), p->fStarMFracIron(),
		    p->fStarMFracOxygen());

    double dSNIaMassStore=0.0;  /* Stores mass loss of Ia so as not to 
				   double count it in wind feedback */
                                    
    /*
     * Call all the effects in order and accumulate them.
     */
    for(int j = 0; j < NFEEDBACKS; j++) {
	switch (j) {
	case FB_SNII:
        if (bAGORAFeedback) break;
	    sn.CalcSNIIFeedback(&sfEvent, dTime, dDeltaYr, &fbEffects);
	    if (sn.dESN > 0)
		p->fNSN() = fbEffects.dEnergy*MSOLG*fbEffects.dMassLoss/sn.dESN;
	    break;
	case FB_SNIA:
        if (bAGORAFeedback) break;
	    sn.CalcSNIaFeedback(&sfEvent, dTime, dDeltaYr, &fbEffects);
	    dSNIaMassStore=fbEffects.dMassLoss;
	    break;
	case FB_WIND:
        if (bAGORAFeedback) break;
	    CalcWindFeedback(&sfEvent, dTime, dDeltaYr, &fbEffects);
	    if(dSNIaMassStore < fbEffects.dMassLoss)
		fbEffects.dMassLoss -= dSNIaMassStore;
	    break;
	case FB_UV:
        if (bAGORAFeedback) break;
	    CalcUVFeedback(dTime, dDeltaYr, &fbEffects);
	break;
    case FB_AGORA:
        if (!bAGORAFeedback) break;
        sn.CalcAGORAFeedback(&sfEvent, dTime, dDeltaYr, &fbEffects);
        p->fNSN() = 0.0;
        break;
	default:
	    CkAssert(0);
	    }

	// Convert to system units
	fbEffects.dMassLoss *= MSOLG/dGmUnit;
	fbEffects.dEnergy /= dErgPerGmUnit;
	
	dTotMassLoss += fbEffects.dMassLoss;
	p->fStarESNrate() += fbEffects.dEnergy*fbEffects.dMassLoss;
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
    CkAssert(p->mass > dTotMassLoss);
    
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

   // Convert to a rate, except if we are using AGORA feedback   
   if (!bAGORAFeedback)
        p->fStarESNrate() /= dDelta;
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
    CkAssert(dMaxMass >= dMinMass);

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


/*
 * Methods to distribute stellar feedback
 */
int DistStellarFeedbackSmoothParams::isSmoothActive(GravityParticle *p) 
{
    return p->isStar();
    }

void DistStellarFeedbackSmoothParams::initTreeParticle(GravityParticle *p1)
{
    /* Convert energy and metals to non-specific quantities (not per mass)
     * to make it easier to divvy up SN energy and metals.  
     */
    
    if(TYPETest(p1, TYPE_GAS)){
      p1->fESNrate() *= p1->mass;
      p1->fMetals() *= p1->mass;    
      p1->fMFracOxygen() *= p1->mass;    
      p1->fMFracIron() *= p1->mass;    

     if (fb.bAGORAFeedback) {
        p1->u() *= p1->mass;
        p1->uPred() *= p1->mass;
      }
    }
    
    }

void DistStellarFeedbackSmoothParams::initSmoothCache(GravityParticle *p1)
{
    /*
     * Warning: kludgery.  We need to accumulate mass in the cached
     * particle, but we also need to keep the original mass around.
     * Let's use the curlv field in the cached particle copy to hold the original
     * mass.  Note: original particle curlv's never modified.
     */
    p1->curlv().x = p1->mass;

    /*
     * Zero out accumulated quantities.
     */
    p1->fESNrate() = 0.0;
    p1->fMetals() = 0.0;
    p1->fMFracOxygen() = 0.0;
    p1->fMFracIron() = 0.0;

    if (fb.bAGORAFeedback) {
      p1->u() = 0.0;
      p1->uPred() = 0.0;
    }
    }

void DistStellarFeedbackSmoothParams::combSmoothCache(GravityParticle *p1,
						      ExternalSmoothParticle *p2)
{
    /*
     * See kludgery notice above.
     */
    double fAddedMass = p2->mass - p2->curlv.x;
    
    p1->mass += fAddedMass;
    p1->fESNrate() += p2->fESNrate;
    p1->fMetals() += p2->fMetals;
    p1->fMFracOxygen() += p2->fMFracOxygen;
    p1->fMFracIron() += p2->fMFracIron;
    p1->fTimeCoolIsOffUntil() = max( p1->fTimeCoolIsOffUntil(),
				     p2->fTimeCoolIsOffUntil );
    if (fb.bAGORAFeedback) {
     p1->u() += p2->u;
     p1->uPred() += p2->uPred;
    }
    }

void DistStellarFeedbackSmoothParams::DistFBMME(GravityParticle *p,int nSmooth, pqSmoothNode *nList)
{
    GravityParticle *q;
    double fNorm,ih2,r2,rs,rstot,fNorm_u,fNorm_Pres,fAveDens;
    int i,counter;
    int nHeavy = 0;
    
    if ( p->fMSN() == 0.0 ){return;} /* Is there any feedback mass? */
    CkAssert(TYPETest(p, TYPE_STAR));
    CkAssert(nSmooth > 0);
    ih2 = invH2(p);
    rstot = 0.0;  
    fNorm_u = 0.0;
    fNorm_Pres = 0.0;
    fAveDens = 0.0;
    
    fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
    for (i=0;i<nSmooth;++i) {
	double fDist2 = nList[i].fKey;
	r2 = fDist2*ih2;            
	rs = KERNEL(r2);
	q = nList[i].p;
	if(q->mass > fb.dMaxGasMass) {
	    nHeavy++;
	    continue; /* Skip heavy particles */
	    }
	fNorm_u += q->mass*rs;
        CkAssert(TYPETest(q, TYPE_GAS));
	rs *= fNorm;
	fAveDens += q->mass*rs;
	fNorm_Pres += q->mass*q->uPred()*rs;
	}
    if(fNorm_u == 0.0) {
        CkError("Got %d heavies: no feedback\n", nHeavy);
	}
	    
    CkAssert(fNorm_u > 0.0);  	/* be sure we have at least one neighbor */
    fNorm_Pres *= (gamma-1.0);
    
    fNorm_u = 1./fNorm_u;
    counter=0;
    for (i=0;i<nSmooth;++i) {
	double weight;
	q = nList[i].p;
	if(q->mass > fb.dMaxGasMass) {
	    continue; /* Skip heavy particles */
	    }
	double fDist2 = nList[i].fKey;
	r2 = fDist2*ih2;            
	rs = KERNEL(r2);
	/* Remember: We are dealing with total energy rate and total metal
	 * mass, not energy/gram or metals per gram.  
	 * q->mass is in product to make units work for fNorm_u.
	 */
#ifdef VOLUMEFEEDBACK
	weight = rs*fNorm_u*q->mass/q->fDensity;
#else
	weight = rs*fNorm_u*q->mass;
#endif
	if (p->fNSN() == 0.0) {
	  if (fb.bAGORAFeedback) {
  		q->u() += weight*p->fStarESNrate();
		q->uPred() += weight*p->fStarESNrate();
	  } else q->fESNrate() += weight*p->fStarESNrate();
	}
	q->fMetals() += weight*p->fSNMetals();
	q->fMFracOxygen() += weight*p->fMOxygenOut();
	q->fMFracIron() += weight*p->fMIronOut();
	q->mass += weight*p->fMSN();
	}
    }

void DistStellarFeedbackSmoothParams::fcnSmooth(GravityParticle *p,int nSmooth, pqSmoothNode *nList)
{
    GravityParticle *q;
    double fNorm,ih2,r2,rs,rstot,fNorm_u,fNorm_Pres,fAveDens,f2h2;
    double fBlastRadius,fShutoffTime,fmind;
    double dAge, aFac, dCosmoDenFac;
    int i,counter,imind;
    
    if ( p->fMSN() == 0.0 ){return;}
    
    /* "Simple" ejecta distribution (see function above) */
    DistFBMME(p,nSmooth,nList);
    
    if (p->fNSN() == 0) {return;}
    if ( p->fTimeForm() < 0.0 ) {return;}
    
    /* The following ONLY deals with SNII Energy distribution */
    CkAssert(TYPETest(p, TYPE_STAR));
    ih2 = invH2(p);
    aFac = a;
    dCosmoDenFac = aFac*aFac*aFac;
    rstot = 0.0;  
    fNorm_u = 0.0;
    fNorm_Pres = 0.0;
    fAveDens = 0.0;
    dAge = dTime - p->fTimeForm();
    if (dAge == 0.0) return;
    
    fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
    for (i=0;i<nSmooth;++i) {
	double fDist2 = nList[i].fKey;
	r2 = fDist2*ih2;            
	rs = KERNEL(r2);
	q = nList[i].p;
        CkAssert(TYPETest(q, TYPE_GAS));
	fNorm_u += q->mass*rs;
	rs *= fNorm;
	fAveDens += q->mass*rs;
	fNorm_Pres += q->mass*q->uPred()*rs;
	}
    fNorm_Pres *= (gamma-1.0)/dCosmoDenFac;
    fAveDens /= dCosmoDenFac;
    if (fb.sn.iNSNIIQuantum > 0) {
	/* McCray + Kafatos (1987) ApJ 317 190*/
	fBlastRadius = fb.dRadPreFactor*pow(p->fNSN() / fAveDens, 0.2) * 
	    pow(dAge,0.6)/aFac; /* eq 3 */
	/* TOO LONG    fShutoffTime = dTimePreFactor*pow(p->fMetals, -1.5)*
	   pow(p->fNSN,0.3) / pow(fAveDens,0.7);*/
	}
    else {
	/* from McKee and Ostriker (1977) ApJ 218 148 */
	fBlastRadius = fb.dRadPreFactor*pow(p->fNSN(),0.32)*
	    pow(fAveDens,-0.16)*pow(fNorm_Pres,-0.2)/aFac;
	}
    if (fb.bShortCoolShutoff){
	/* End of snowplow phase */
	fShutoffTime = fb.dTimePreFactor*pow(p->fNSN(),0.31)*
	    pow(fAveDens,0.27)*pow(fNorm_Pres,-0.64);
	} else{        /* McKee + Ostriker 1977 t_{max} */
	fShutoffTime = fb.dTimePreFactor*pow(p->fNSN(),0.32)*
	    pow(fAveDens,0.34)*pow(fNorm_Pres,-0.70);
	}
    /* Shut off cooling for 3 Myr for stellar wind */
    if (p->fNSN() < fb.sn.iNSNIIQuantum)
	fShutoffTime= 3e6 * SECONDSPERYEAR / fb.dSecUnit;
    /* Limit cooling shutoff time */
    if(fShutoffTime > fb.dMaxCoolShutoff)
        fShutoffTime = fb.dMaxCoolShutoff;
    
    fmind = p->fBall*p->fBall;
    imind = 0;
    if ( p->fStarESNrate() > 0.0 ) {
	if(fb.bSmallSNSmooth) {
	    /* Change smoothing radius to blast radius 
	     * so that we only distribute mass, metals, and energy
	     * over that range. 
	     */
	    f2h2 = fBlastRadius*fBlastRadius;
	    ih2 = 4.0/f2h2;
	    }
	
	rstot = 0.0;  
	fNorm_u = 0.0;
	
	for (i=0;i<nSmooth;++i) {
	    double fDist2 = nList[i].fKey;
	    if ( fDist2 < fmind ){imind = i; fmind = fDist2;}
            if ( !fb.bSmallSNSmooth || fDist2 < f2h2 ) {
		r2 = fDist2*ih2;            
		rs = KERNEL(r2);
		q = nList[i].p;
#ifdef VOLUMEFEEDBACK
		fNorm_u += q->mass/q->fDensity*rs;
#else
		fNorm_u += q->mass*rs;
#endif
		assert(TYPETest(q, TYPE_GAS));
		}
	    }
	}
       
    /* If there's no gas particle within blast radius,
       give mass and energy to nearest gas particle. */
    if (fNorm_u ==0.0){
	double fDist2 = nList[imind].fKey;
	r2 = fDist2*ih2;            
	rs = KERNEL(r2);
	/*
	 * N.B. This will be NEGATIVE, but that's OK since it will
	 * cancel out down below.
	 */
#ifdef VOLUMEFEEDBACK
	fNorm_u = nList[imind].p->mass/nList[imind].p->fDensity*rs;
#else
	fNorm_u = nList[imind].p->mass*rs;
#endif
	}
    
    assert(fNorm_u != 0.0);
    fNorm_u = 1./fNorm_u;
    counter=0;
    for (i=0;i<nSmooth;++i) {
	double weight;
	double fDist2 = nList[i].fKey;
	q = nList[i].p;
	if (fb.bSmallSNSmooth) {
	    if ( (fDist2 <= f2h2) || (i == imind) ) {
		if( fb.bSNTurnOffCooling && 
		    (fBlastRadius*fBlastRadius >= fDist2)) {
		    q->fTimeCoolIsOffUntil() = max(q->fTimeCoolIsOffUntil(),
						   dTime + fShutoffTime);
		    }
		
		counter++;  
		r2 = fDist2*ih2;
		rs = KERNEL(r2);
		/* Remember: We are dealing with total energy rate and total metal
		 * mass, not energy/gram or metals per gram.  
		 * q->mass is in product to make units work for fNorm_u.
		 */
#ifdef VOLUMEFEEDBACK
		weight = rs*fNorm_u*q->mass/q->fDensity;
#else
		weight = rs*fNorm_u*q->mass;
#endif
		q->fESNrate() += weight*p->fStarESNrate();
		}
	    } else {
	    double fDist2 = nList[i].fKey;
	    r2 = fDist2*ih2;  
	    rs = KERNEL(r2);
	    /* Remember: We are dealing with total energy rate and total metal
	     * mass, not energy/gram or metals per gram.  
	     * q->mass is in product to make units work for fNorm_u.
	     */
#ifdef VOLUMEFEEDBACK
	    weight = rs*fNorm_u*q->mass/q->fDensity;
#else
	    weight = rs*fNorm_u*q->mass;
#endif
	    q->fESNrate() += weight*p->fESNrate();
	    /*		printf("SNTEST: %d %g %g %g %g\n",q->iOrder,weight,sqrt(q->r[0]*q->r[0]+q->r[1]*q->r[1]+q->r[2]*q->r[2]),q->fESNrate,q->fDensity);*/
	    
	    if ( p->fESNrate() > 0.0 && fb.bSNTurnOffCooling && 
		 (fBlastRadius*fBlastRadius >= fDist2)){
		q->fTimeCoolIsOffUntil() = max(q->fTimeCoolIsOffUntil(),
					       dTime + fShutoffTime);       
		counter++;
		}
	    /*	update mass after everything else so that distribution
		is based entirely upon initial mass of gas particle */
	    } 
	}
    }

void DistStellarFeedbackSmoothParams::postTreeParticle(GravityParticle *p1)
{
    /* Convert energy and metals back to specific quantities (per mass)
       because we are done with our conservative calculations */
    
    if(p1->isGas()){
	p1->fESNrate() /= p1->mass;
	p1->fMetals() /= p1->mass;    
	p1->fMFracIron() /= p1->mass;    
	p1->fMFracOxygen() /= p1->mass;    
	if (fb.bAGORAFeedback) {
		p1->u() /= p1->mass;
		p1->uPred() /= p1->mass;
	}
	}
    
    }

/// @brief total feedback quantities for conservation check
/// 
/// Sums are contributed back to main chare.
/// @param bPreDist Is this before the feedback gets distributed.  In
/// this case the "out" quantities need to be summed.
void
TreePiece::massMetalsEnergyCheck(int bPreDist, const CkCallback& cb)
{
    double dTotals[5];
    for(int j = 0; j < 5; j++)
	dTotals[j] = 0.0;
    
    for(unsigned int i = 1; i <= myNumParticles; ++i) {
	GravityParticle *p = &myParticles[i];
	dTotals[0] += p->mass;
	if(p->isGas()) {
	    dTotals[1] += p->mass*p->fMetals();
	    dTotals[2] += p->mass*p->fMFracOxygen();
	    dTotals[3] += p->mass*p->fMFracIron();
	    dTotals[4] += p->mass*p->fESNrate();
	    }
	if(p->isStar()) {
	    dTotals[1] += p->mass*p->fStarMetals();
	    dTotals[2] += p->mass*p->fStarMFracOxygen();
	    dTotals[3] += p->mass*p->fStarMFracIron();
	    if(bPreDist) { // sum up the quantities that will be distributed
		dTotals[0] += p->fMSN();
		dTotals[1] += p->fSNMetals();
		dTotals[2] += p->fMOxygenOut();
		dTotals[3] += p->fMIronOut();
		dTotals[4] += p->fStarESNrate();
		}
	    }
	}
    contribute(sizeof(dTotals), dTotals, CkReduction::sum_double, cb);
    }

