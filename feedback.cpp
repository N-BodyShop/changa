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
#ifdef SUPERBUBBLE
    bSmallSNSmooth = 0; //Don't use the blastwave smoothing when using superbubble
#else
    bSmallSNSmooth = 1; 
#endif
    prmAddParam(prm,"bSmallSNSmooth", paramBool, &bSmallSNSmooth, sizeof(int),
		"bSmallSNSmooth",
		"<smooth SN ejecta over blast or smoothing radius> = blast radius");
    sn.bUseStoch = 0;
    prmAddParam(prm,"bUseStoch",paramInt, &sn.bUseStoch, sizeof(int),
            "usestoch","<Enable stochastic IMF>");
    sn.dStochCut = 8.0;
    prmAddParam(prm,"dStochCut",paramDouble, &sn.dStochCut, sizeof(double),
            "stochcut","<Cut off mass for stochastic IMF> = 8.0");
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
#ifdef SUPERBUBBLE
    bSNTurnOffCooling = 0; //Don't use cooling shutoffs with superbubble
#else
    bSNTurnOffCooling = 1;
#endif
    prmAddParam(prm,"bSNTurnOffCooling", paramBool, &bSNTurnOffCooling,
		sizeof(int), "bSNTurnOffCooling", "<Do SN turn off cooling> = 1");
#ifdef SUPERBUBBLE
    nSmoothFeedback = 1; //Only use a single neighbour for superbubble feedback
    prmAddParam(prm,"nSmoothFeedback", paramInt,&nSmoothFeedback, sizeof(int),
		"s", "<number of particles to smooth feedback over> = 1");
#else
    nSmoothFeedback = 64;
    prmAddParam(prm,"nSmoothFeedback", paramInt,&nSmoothFeedback, sizeof(int),
		"s", "<number of particles to smooth feedback over> = 64");
#endif
    dMaxCoolShutoff = 3.0e7;
    prmAddParam(prm,"dMaxCoolShutoff", paramDouble, &dMaxCoolShutoff,
		sizeof(double), "fbMaxCoolOff",
		"<Maximum time to shutoff cooling in years> = 3e7");
    dEarlyFeedbackFrac = 0.0;
    prmAddParam(prm,"dEarlyFeedbackFrac", paramDouble, &dEarlyFeedbackFrac,
		sizeof(double), "fbEarlyFrac",
		"<Fraction of SNII energy to put in early feedback> = 0.0");
    dFBInitialMassLoad = 0.0;
    prmAddParam(prm,"dFBInitialMassLoad", paramDouble, &dFBInitialMassLoad,
		sizeof(double), "dFBIML",
		"<Initial Mass Loading for Feedback Ejecta> = 0.0");
	dMultiPhaseMinTemp = 1e5;
	prmAddParam(prm,"dMultiPhaseMinTemp",paramDouble,&dMultiPhaseMinTemp,
				sizeof(double),"multitmin",
				"<Temperature threshold to use multiphase feedback> = 1e6");
    dMultiPhaseMaxTime = 1e8;
    prmAddParam(prm,"dMultiPhaseMaxTime",paramDouble,&dMultiPhaseMaxTime,
                sizeof(double),"evaptmax",
                "<Absolute maximum lifetime for multiphase (yrs)> = 1e8");
#ifdef SPLITGAS
	dInitGasMass = -1;
	prmAddParam(prm,"dInitGasMass", paramDouble, &dInitGasMass,
		    sizeof(double), "stInitGas",
		    "<Initial mass of a gas particle> = -1");
#endif
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
    else CkAbort("Unknown achIMF parameter");
    sn.imf = imf;

#include "physconst.h"
#ifndef SUPERBUBBLE
    /*
     * Make sure that the parameters are sensible for
     * superbubble feedback.
     */
    if(!prmSpecified(prm, "nSmoothFeedback"))
	nSmoothFeedback = param.nSmooth;
    if (sn.dESN > 0.0) bSmallSNSmooth = 1;
    else bSmallSNSmooth = 0;
#endif
    param.bDoGas = 1;
    dDeltaStarForm = param.stfm->dDeltaStarForm;
#if defined(STOCH12) || defined(STOCH24)
    if(!sn.bUseStoch) CkPrintf("Running non-stochastic IMF but stochastic IMF compiled in\n");
#else
    if(sn.bUseStoch) CkAbort("Stochastic IMF requested but not compiled in");
#endif
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
    dMultiPhaseMaxTime *= SECONDSPERYEAR/param.dSecUnit;
    // Normalization of Early Feedback
    /// Total SNe ergs/solar mass of stars
    double dSNETotal = sn.dESN*(imf->CumNumber(sn.dMSNIImin)
                         - imf->CumNumber(sn.dMSNIImax))/imf->CumMass(0.0);
    
    CkPrintf("SNII feedback: %g ergs/solar mass\n", dSNETotal);
    dEarlyETotal = dSNETotal*dEarlyFeedbackFrac;
#ifndef DTADJUST
    CkMustAssert(!bAgoraFeedback, "DTADJUST must be enabled to use AGORA feedback\n");
#endif
    }

/** @brief This routine is called when AGORA feedback is enabled. It checks for any star
 *  particles that will have a feedback event in the next timestep and puts the neighboring
 *  gas particles onto a smaller timestep. Because the amount of energy injected is very
 *  large, particles need to be placed on a small timestep BEFORE the first force calculation 
 *  is done to avoid significant integration errors.
 *
 *  @param dTime The current simulation time (in years)
 *  @param dDelta The size of the next timestep (in years)
 *  @param dTimeToSF The time until the next star formation event (in years)
 */
void Main::AGORAfeedbackPreCheck(double dTime, double dDelta, double dTimeToSF)
{
    AGORApreCheckSmoothParams pAPC(TYPE_GAS, 0, param.csm, dTime, dDelta, 
                         param.dConstGamma, param.dEtaCourant, dTimeToSF, param.feedback,
                         param.dMsolUnit, param.dErgPerGmUnit);
    double dfBall2OverSoft2 = 4.0*param.dhMinOverSoft*param.dhMinOverSoft;
    treeProxy.startSmooth(&pAPC, 0, param.feedback->nSmoothFeedback,
              dfBall2OverSoft2, CkCallbackResumeThread());
    }

void AGORApreCheckSmoothParams::initSmoothCache(GravityParticle *p1)
{
    #ifdef DTADJUST 
    if (p1-> rung >= activeRung) {
            p1->dtNew() = FLT_MAX;
        }
    #endif
    }

void AGORApreCheckSmoothParams::combSmoothCache(GravityParticle *p1,
                              ExternalSmoothParticle *p2)
{
    #ifdef DTADJUST
        if (p2->dtNew < p1->dtNew())
            p1->dtNew() = p2->dtNew;
    #endif
    }

void AGORApreCheckSmoothParams::fcnSmooth(GravityParticle *p, int nSmooth, pqSmoothNode *nList)
{
    GravityParticle *q;
    double dTimeYr, dDeltaYr, dtCourantFac, PoverRho, c_s, uPred, dtNew, pMassMsol, snEerg;
    int i;

    dTimeYr = dTime*fb.dSecUnit/SECONDSPERYEAR ;
    dDeltaYr = max(dDelta, fb.dDeltaStarForm)*fb.dSecUnit/SECONDSPERYEAR ;

    if(p->isStar() && p->fTimeForm() >= 0.0) {
        double dStarLtimeMin = dTimeYr - p->fTimeForm()*fb.dSecUnit/SECONDSPERYEAR;
        double dStarLtimeMax = dStarLtimeMin + dDeltaYr;
        if (dStarLtimeMin < fb.sn.AGORAsnTime && dStarLtimeMax > fb.sn.AGORAsnTime) {
            if (verbosity)
                CkPrintf("AGORA feedback imminent...warning %d neighbors\n", nSmooth);
            for (i=0; i<nSmooth; i++) {
                q = nList[i].p;

                dtCourantFac = etaCourant*a*2.0/1.6;
                pMassMsol = p->mass*dMsolUnit;
                snEerg = fb.sn.AGORAsnE*pMassMsol/fb.sn.AGORAsnPerMass;
                uPred = snEerg/(pMassMsol*MSOLG);
                uPred /= dErgPerGmUnit;
                PoverRho = (gamma-1.0)*uPred;
                c_s = sqrt(gamma*PoverRho);
                dtNew = dtCourantFac*0.5*q->fBall/(2.0*c_s);

#ifdef DTADJUST
                // Reducing the time step size means that the feedback is no longer going to
                // occur in the next step. Reduce dtNew by no larger than a factor of 2 to reduce 
                // the amount of unecessary steps.
                if (q->dt > dtNew) {
                    if (dtNew <= timeToSF/2.0) {
                        q->dtNew() = timeToSF/2.0;
                        }
                    else {
                        q->dtNew() = dtNew;
                        }
                    }
#endif
                }
            }
        }
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
    double tFB = CkWallTimer() - startTime;
    // Overload tAdjust with FB stellar evolution
    timings[PHASE_FEEDBACK].tAdjust += tFB;
    
    if(verbosity)
      CkPrintf("Distribute Stellar Feedback ...\n ");
    // Need to build tree since we just did addDelParticle.
    //
    buildTree(PHASE_FEEDBACK);
    DistStellarFeedbackSmoothParams pDSFB(TYPE_GAS, 0, param.csm, dTime, 
					  param.dConstGamma, param.feedback);
    double dfBall2OverSoft2 = 4.0*param.dhMinOverSoft*param.dhMinOverSoft;
    treeProxy.startSmooth(&pDSFB, 0, param.feedback->nSmoothFeedback,
			  dfBall2OverSoft2, CkCallbackResumeThread());

#ifdef SPLITGAS
    if(param.feedback->dInitGasMass > 0)
    {
        CkReductionMsg *msgCounts;
        treeProxy.SplitGas(param.feedback->dInitGasMass, CkCallbackResumeThread((void*&)msgCounts));
        int *dCounts = (int *)msgCounts->getData();
        if(verbosity)
        CkPrintf("%d Gas Particles Split\n", *dCounts);
        delete msgCounts;
    }
#endif
#ifdef SUPERBUBBLE
    CkPrintf("Promote cold shell to hot... \n");
    double a = csmTime2Exp(param.csm,dTime);
    PromoteToHotGasSmoothParams pPHG(TYPE_GAS, 0, param.dEvapCoeffCode*a, param.dEvapMinTemp,
            param.dErgPerGmUnit, param.dGmPerCcUnit, param.stfm->dDeltaStarForm, dTime);
    treeProxy.startSmooth(&pPHG, 1, param.nSmooth,
              dfBall2OverSoft2, CkCallbackResumeThread());
    ShareWithHotGasSmoothParams pSHG(TYPE_GAS, 0, param.dEvapMinTemp,
            param.dErgPerGmUnit, param.dGmPerCcUnit); 
    treeProxy.startReSmooth(&pSHG, CkCallbackResumeThread());
#endif
    treeProxy.finishNodeCache(CkCallbackResumeThread());

#ifdef SPLITGAS
    addDelParticles();//Don't forget to run an addDelParticles after a split
#endif
    double tDFB = CkWallTimer() - startTime;
    timings[PHASE_FEEDBACK].tuDot += tDFB; // Overload tuDot for feedback.
    CkPrintf("Stellar Feedback Calculated, Wallclock %f secs\n", tFB);

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
void Main::initLWData(){
    dMProxy.initLWData(CkCallbackResumeThread());
    treeProxy.initLWData(CkCallbackResumeThread());
}
void DataManager::initLWData(const CkCallback& cb)
{
    lwInitData(LWData);
    contribute(cb);
}
void TreePiece::initLWData(const CkCallback& cb)
{
    dm = (DataManager*)CkLocalNodeBranch(dataManagerID);
    contribute(cb);
}


///
/// processor specific method for stellar feedback
/// 
void TreePiece::Feedback(const Fdbk &fb, double dTime, double dDelta, const CkCallback& cb)
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
        if(p->isStar()) {
            if(p->fTimeForm() >= 0.0) {
                if(fb.sn.bUseStoch){
                    fb.DoFeedback(p, dTime, dDeltaYr, fbTotals, dm->LWData, rndGen);
                }
                else {
                    fb.DoFeedback(p, dTime, dDeltaYr, fbTotals, NULL, rndGen); /* If not stochastic IMF, LWData table never initialized*/
                }
            }
            else {  // zero out feedback quantities for Sinks
                p->fMSN() = 0.0;
                p->fSNMetals() = 0.0;
                p->fMIronOut() = 0.0;
                p->fMOxygenOut() = 0.0;
                p->fStarESNrate() = 0.0;
                p->fNSN() = 0.0;
            }
        }
	else if(p->isGas()){
	  CkAssert(p->u() >= 0.0);
	  CkAssert(p->uPred() >= 0.0);
	  CkAssert(p->u() < LIGHTSPEED*LIGHTSPEED/fb.dErgPerGmUnit);
	  CkAssert(p->uPred() < LIGHTSPEED*LIGHTSPEED/fb.dErgPerGmUnit);
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
/// @param LWData pointer to the table containing Lyman-Werner data if stochastic IMF sampling is used.
/// @param rndGen Random number generator reference

void Fdbk::DoFeedback(GravityParticle *p, double dTime, double dDeltaYr, 
                      FBEffects *fbTotals,
                      LWDATA *LWData,
                      Rand& rndGen) const
{
    double dTotMassLoss, dTotMetals, dTotMOxygen, dTotMIron, dDelta;
    dTotMassLoss = dTotMetals = dTotMOxygen = dTotMIron = 0.0;
    dDelta = dDeltaYr*SECONDSPERYEAR / dSecUnit;
    p->fStarESNrate() = 0.0;
    p->fNSN() = 0.0;
  
    FBEffects fbEffects;
    // Particle properties that will be sent to feedback
    // methods in normal units (M_sun + seconds)
    
    SFEvent sfEvent;
    if(sn.bUseStoch){
        sfEvent = SFEvent(p->fMassForm()*dGmUnit/MSOLG, 
                p->fTimeForm()*dSecUnit/SECONDSPERYEAR,
                p->fStarMetals(), p->fStarMFracIron(),
                p->fStarMFracOxygen(), p->fLowNorm(), p->rgfHMStars());
    } else {
        sfEvent = SFEvent(p->fMassForm()*dGmUnit/MSOLG, 
                p->fTimeForm()*dSecUnit/SECONDSPERYEAR,
                p->fStarMetals(), p->fStarMFracIron(),
                p->fStarMFracOxygen());
    }

    double dSNIaMassStore=0.0;  /* Stores mass loss of Ia so as not to 
				   double count it in wind feedback */
                                    
    /*
     * Call all the effects in order and accumulate them.
     */
    for(int j = 0; j < NFEEDBACKS; j++) {
	switch (j) {
	case FB_SNII:
	    if (bAGORAFeedback) break;
	    sn.CalcSNIIFeedback(&sfEvent, dTime, dDeltaYr, &fbEffects, rndGen);
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
	    CalcUVFeedback(&sfEvent, dTime, dDeltaYr, &fbEffects);
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

#ifdef COOLING_MOLECULARH /* Calculates LW radiation from a stellar particle of a given age and mass (assumes Kroupa IMF), CC */
    p->dStarLymanWerner() = CalcLWFeedback(&sfEvent, dTime, dDeltaYr, LWData);
#endif /*COOLING_MOLECULARH*/
}

void Fdbk::CalcWindFeedback(SFEvent *sfEvent, double dTime, /* current time in years */
			    double dDelta, /* length of timestep (years) */
			    FBEffects *fbEffects) const
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

    double dMinCumMass;
    double dMaxCumMass;
    /* If using stochastic, the CumMass method doesn't need to be renormalized */
    if(!sn.bUseStoch){
        dMinCumMass = imf->CumMass(dMinMass);
        dMaxCumMass = imf->CumMass(dMaxMass);
    } else {
        dMinCumMass = imf->CumMassStoch(dMinMass, sfEvent->dLowNorm, sfEvent->rgdHMStars, sn.dStochCut);
        dMaxCumMass = imf->CumMassStoch(dMaxMass, sfEvent->dLowNorm, sfEvent->rgdHMStars, sn.dStochCut);
    }
	double dMTot = imf->CumMass(0.0);
	/* Find out mass fraction of dying stars, then multiply by the original
	   mass of the star particle */
    /* If using stochastic, the CumMass method doesn't need to be renormalized */
	if (dMTot == 0.0){
	  dMDying = 0.0;
	} else { 
        if(!sn.bUseStoch){
            dMDying = (dMinCumMass - dMaxCumMass)/dMTot;
            dMDying *= sfEvent->dMass;
        } else {
            dMDying = (dMinCumMass - dMaxCumMass);
        }
    }

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

/// This is currently used for the early feedback scheme: energy is
/// injected before the first SNII go off.
void Fdbk::CalcUVFeedback(SFEvent *sfEvent, /**< Star to process */
                          double dTime, /**< current time in years */
                          double dDelta, /**< length of timestep (years) */
                          FBEffects *fbEffects /**< [out] effects */) const
{
    fbEffects->dMassLoss = 0.0;
    fbEffects->dEnergy = 0.0;
    fbEffects->dMetals = 0.0;
    fbEffects->dMIron = 0.0;
    fbEffects->dMOxygen = 0.0;
    /* stellar lifetimes corresponding to beginning and end of 
       current timestep with respect to starbirth time in yrs */
    double dStarLtimeMin = dTime - sfEvent->dTimeForm; 
    double dStarLtimeMax = dStarLtimeMin + dDelta;
    /* masses corresponding to these stellar lifetimes in solar masses */
    double dMStarMin = pdva.StarMass(dStarLtimeMax, sfEvent->dMetals); 
    double dMStarMax = pdva.StarMass(dStarLtimeMin, sfEvent->dMetals); 
    
    if(dMStarMax < sn.dMSNIImax)
        return;                 // Supernove are going off.
    double dEndEarlyTime = pdva.Lifetime(sn.dMSNIImax, sfEvent->dMetals);
    if(dMStarMin < sn.dMSNIImax)
        dStarLtimeMax = dEndEarlyTime;

    double dEarlyRate = dEarlyETotal*sfEvent->dMass/dEndEarlyTime;
    /// Energy injected during this step (in ergs)
    double dEEarly = dEarlyRate*(dStarLtimeMax - dStarLtimeMin);
    /// no mass is lost: use E = m c^2 and convert to Solar Masses
    if(dEEarly > 0.0) {
        fbEffects->dMassLoss = dEEarly/(MSOLG*LIGHTSPEED*LIGHTSPEED);
        fbEffects->dEnergy = dEEarly/(MSOLG*fbEffects->dMassLoss);
        }
    }

#ifdef COOLING_MOLECULARH
/*----  Lyman-Warner Radiation from young stars*/
double Fdbk::CalcLWFeedback(SFEvent *sfEvent, double dTime, /* current time in years */
			    double dDelta /* length of timestep (years) */ , LWDATA *LWData /*LW Tabular Data*/) const
{
    double dAge1log, dAge2log, dLW1log, dLW2log, dStarAge;
    double dA0old =  70.908586,
      dA1old = -4.0643123;

    dStarAge = dTime - sfEvent->dTimeForm;
    if (dStarAge >= 0 ) { /*Test that it is not a Black Hole */
      if (dStarAge != 0) {
	dAge1log = log10(dStarAge);
      }
      else {
	/*Avoids log of zero error*/
	dAge1log = log10(dDelta);
      }
      if (dAge1log < 7) dAge1log = 7.0;
      dAge2log = log10(dStarAge + dDelta);
      if (dAge2log < 7) dAge2log = 7.0;
      if (!sn.bUseStoch){
          if (dAge1log < 9.0) dLW1log = calcLogSSPLymanWerner(dAge1log,log10(sfEvent->dMass*MSOLG/dGmUnit));
          else dLW1log = dA0old + dA1old*dAge1log + log10(sfEvent->dMass*MSOLG/dGmUnit); /*Close to zero*/
          if (dAge2log < 9.0) dLW2log = calcLogSSPLymanWerner(dAge2log,log10(sfEvent->dMass*MSOLG/dGmUnit));
          else dLW2log = dA0old + dA1old*dAge2log + log10(sfEvent->dMass*MSOLG/dGmUnit); /*Close to zero*/
      }
      else{ // stochastic calculation must be done in Msol and then converted to system units
          if (dAge1log < 9.0){
              double dMax8Log1 = calcLogMax8LymanWerner(dAge1log,log10(sfEvent->dLowNorm*MSOLG/dGmUnit));
              double dStochLog1 = calcLogStochLymanWerner(dAge1log, sfEvent->rgdHMStars, LWData);
              dStochLog1 += log10(MSOLG/dGmUnit); // converting stoch portion to system units
              dLW1log = log10(pow(10, dMax8Log1) + pow(10, dStochLog1));
          }   
          else dLW1log = dA0old + dA1old*dAge1log + log10(sfEvent->dLowNorm*MSOLG/dGmUnit); /*Close to zero*/
          if (dAge2log < 9.0){
              double dMax8Log2 = calcLogMax8LymanWerner(dAge2log,log10(sfEvent->dLowNorm*MSOLG/dGmUnit));
              double dStochLog2 = calcLogStochLymanWerner(dAge2log, sfEvent->rgdHMStars, LWData);
              dStochLog2 += log10(MSOLG/dGmUnit); // converting stoch portion to system units
              dLW2log = log10(pow(10, dMax8Log2) + pow(10, dStochLog2));
          }
          else dLW2log = dA0old + dA1old*dAge2log + log10(sfEvent->dLowNorm*MSOLG/dGmUnit); /*Close to zero*/
      }
      return log10((pow(10,dLW1log) + pow(10,dLW2log))/2);
    }
    else return 0;
}
#endif /*COOLING_MOLECULARH*/

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
#ifdef SUPERBUBBLE
        /*
         * We must scale the energy injection rate by the mass that contains it
         */
      if (p1->massHot() > 0) p1->fESNrate() *= p1->massHot();
      else
#endif
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
#ifdef SUPERBUBBLE
    p1->curlv().y = p1->massHot();
#endif

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
#ifdef SUPERBUBBLE
    double fAddedMassHot = p2->massHot - p2->curlv.y;
    p1->massHot() += fAddedMassHot;
#endif
    p1->fESNrate() += p2->fESNrate;
    p1->fMetals() += p2->fMetals;
    p1->fMFracOxygen() += p2->fMFracOxygen;
    p1->fMFracIron() += p2->fMFracIron;
    p1->fTimeCoolIsOffUntil() = max( p1->fTimeCoolIsOffUntil(),
				     p2->fTimeCoolIsOffUntil );
    p1->dTimeFB() = max( p1->dTimeFB(), p2->dTimeFB );
    if (fb.bAGORAFeedback) {
     p1->u() += p2->u;
     p1->uPred() += p2->uPred;
    }
    }

void DistStellarFeedbackSmoothParams::DistFBMME(GravityParticle *p,int nSmooth, pqSmoothNode *nList)
{
    GravityParticle *q;
    double fNorm,ih2,r2,rs,fNorm_u,fNorm_Pres,fAveDens;
    int i;
    int nHeavy = 0;
    
    if ( p->fMSN() == 0.0 ){return;} /* Is there any feedback mass? */
    if ( p->fTimeForm() < 0 ){return;} /*Check if star particle is a black hole. We store values for BHs in fMSN,FESN,etc that are unrelated so this check is important -- MJT 11/7/13 */
    CkAssert(TYPETest(p, TYPE_STAR));
    CkAssert(nSmooth > 0);
    ih2 = invH2(p);
    fNorm_u = 0.0;
    fNorm_Pres = 0.0;
    fAveDens = 0.0;
    
    fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
    for (i=0;i<nSmooth;++i) {
	double fDist2 = nList[i].fKey;
	r2 = fDist2*ih2;            
#ifdef SUPERBUBBLE
	rs = 1.0; //Use a tophat kernel for superbubble feedback
#else
	rs = KERNEL(r2, nSmooth);
#endif
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
        CkError("WARNING: lonely star skips feedback of mass %g\n", p->fMSN());
        p->mass += p->fMSN(); // mass goes back into star.
        // We could self-enrich to conserve metals, but that could do
        // funny things to the mass metalicity relation.
        // Also prevent any cooling shut-off.
        p->fNSN() = 0.0;
        return;
	}
	    
    CkAssert(fNorm_u > 0.0);  	/* be sure we have at least one neighbor */
    fNorm_Pres *= (gamma-1.0);
    
    fNorm_u = 1./fNorm_u;
    for (i=0;i<nSmooth;++i) {
	double weight;
	q = nList[i].p;
	if(q->mass > fb.dMaxGasMass) {
	    continue; /* Skip heavy particles */
	    }
	double fDist2 = nList[i].fKey;
	r2 = fDist2*ih2;            
#ifdef SUPERBUBBLE
	rs = 1.0; //Use a tophat kernel for superbubble feedback
#else
	rs = KERNEL(r2, nSmooth);
#endif
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
    if(weight > 0) TYPESet(q, TYPE_FEEDBACK);
#ifdef SUPERBUBBLE
    CkAssert(q->uPred() < LIGHTSPEED*LIGHTSPEED/fb.dErgPerGmUnit);
    double Tq = CoolCodeEnergyToTemperature(tp->Cool(), &q->CoolParticle(), q->uPred(),
#ifdef COOLING_GRACKLE
                                            q->fDensity, /* GRACKLE needs density */
#endif
                                            q->fMetals() );
	if(Tq < fb.dMultiPhaseMinTemp && weight > 0 && p->fNSN() > 0.0) { //Only use the multiphase state for cooler particles
		double massHot = q->massHot() + weight*p->fMSN();
		double deltaMassLoad = weight*p->fMSN()*fb.dFBInitialMassLoad;
		if (massHot+deltaMassLoad >= q->mass) { //If deltaMassLoad is bigger than the cold phase, make it all cold.
			deltaMassLoad = q->mass - massHot;
			massHot = q->mass;
			}
		else {
			massHot += deltaMassLoad;
		}
		q->massHot() = massHot;
		CkAssert(q->massHot() >= 0); //Make sure our hot phase has nonzero mass in it
	}
#endif
	}
    }

void DistStellarFeedbackSmoothParams::fcnSmooth(GravityParticle *p,int nSmooth, pqSmoothNode *nList)
{
    GravityParticle *q;
    double fNorm,ih2,r2,rs,fNorm_u,fNorm_Pres,fAveDens,f2h2;
    double fBlastRadius,fShutoffTime,fmind;
    double dAge, aFac, dCosmoDenFac;
    int i,counter,imind;
    
    if ( p->fTimeForm() < 0.0) {return;}  //don't want to do this calculatino for a BH
    if ( p->fMSN() == 0.0 ){return;}
 
    /* "Simple" ejecta distribution (see function above) */
    DistFBMME(p,nSmooth,nList);
    
    if (p->fNSN() == 0) {return;}
    //if ( p->fTimeForm() < 0.0) {return;}
    
    /* The following ONLY deals with SNII Energy distribution */
    CkAssert(TYPETest(p, TYPE_STAR));
    ih2 = invH2(p);
    aFac = a;
    dCosmoDenFac = aFac*aFac*aFac;
    fNorm_u = 0.0;
    fNorm_Pres = 0.0;
    fAveDens = 0.0;
    dAge = dTime - p->fTimeForm();
    if (dAge == 0.0) return;
    
    fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
    for (i=0;i<nSmooth;++i) {
	double fDist2 = nList[i].fKey;
	r2 = fDist2*ih2;            
#ifdef SUPERBUBBLE
	rs = 1.0; //Use a tophat kernel for superbubble feedback
#else
	rs = KERNEL(r2, nSmooth);
#endif
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
	
	fNorm_u = 0.0;
	
	for (i=0;i<nSmooth;++i) {
	    double fDist2 = nList[i].fKey;
	    if ( fDist2 < fmind ){imind = i; fmind = fDist2;}
            if ( !fb.bSmallSNSmooth || fDist2 < f2h2 ) {
		r2 = fDist2*ih2;            
#ifdef SUPERBUBBLE
                rs = 1.0; //Use a tophat kernel for superbubble feedback
#else
		rs = KERNEL(r2, nSmooth);
#endif
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
#ifdef SUPERBUBBLE
	rs = 1.0; //Use a tophat kernel for superbubble feedback
#else
	rs = KERNEL(r2, nSmooth);
#endif
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
		    q->dTimeFB() = dTime;
		    }
		
		counter++;  
		r2 = fDist2*ih2;
#ifdef SUPERBUBBLE
                rs = 1.0; //Use a tophat kernel for superbubble feedback
#else
                rs = KERNEL(r2, nSmooth);
#endif
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
#ifdef SUPERBUBBLE
            rs = 1.0; //Use a tophat kernel for superbubble feedback
#else
	    rs = KERNEL(r2, nSmooth);
#endif
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
	    
	    if ( p->fStarESNrate() > 0.0 && fb.bSNTurnOffCooling && 
		 (fBlastRadius*fBlastRadius >= fDist2)){
		q->fTimeCoolIsOffUntil() = max(q->fTimeCoolIsOffUntil(),
					       dTime + fShutoffTime);       
		q->dTimeFB() = dTime; /* store SN FB time here JMB 2/24/10 */
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
#ifdef SUPERBUBBLE
        if(p1->massHot() > 0) p1->fESNrate() /= p1->massHot();
        else
#endif
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
#ifdef SUPERBUBBLE
	    dTotals[4] += p->massHot()*p->fESNrate();
#endif
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

#ifdef SPLITGAS
void TreePiece::SplitGas(double dInitGasMass, const CkCallback& cb)
{
    int nFormed = 0;
    double norm, uvar, vvar, ux, uy, uz;
    for(unsigned int i = 1; i <= myNumParticles; ++i) {
        GravityParticle *p = &myParticles[i];
        if(p->mass < 1.33*dInitGasMass) continue; //Don't split particles that are too small FOOL
        if(!TYPETest(p, TYPE_GAS)) continue; //Only split heavy gas

        nFormed++;
        norm = 666; // \m/
        while (norm>1.0){ //unit sphere point picking (Marsaglia 1972)
            uvar=2.0*rndGen.dbl()-1.0;  //#random number on [-1,1]
            vvar=2.0*rndGen.dbl()-1.0;
            norm=(uvar*uvar+vvar*vvar);
        }
        norm = sqrt(1.0-norm); //only do one sqrt
        ux=2.0*uvar*norm;
        uy=2.0*vvar*norm;
        uz=1.0-2.0*(uvar*uvar+vvar*vvar);
        p->mass /= 2.0;
#ifdef SUPERBUBBLE
        p->massHot() /= 2.0;
#endif
        GravityParticle daughter = *p;
        extraSPHData extra = *(extraSPHData *)p->extraData;
        daughter.extraData = &extra;
        TYPEReset(&daughter, TYPE_NbrOfACTIVE);
        TYPESet(&daughter, TYPE_GAS);
        daughter.position.x += 0.25*p->fBall*ux;
        daughter.position.y += 0.25*p->fBall*uy;
        daughter.position.z += 0.25*p->fBall*uz;
        daughter.iSplitOrder() = p->iOrder;
        p->position.x -= 0.25*p->fBall*ux;
        p->position.y -= 0.25*p->fBall*uy;
        p->position.z -= 0.25*p->fBall*uz;
        newParticle(&daughter);
    }
    contribute(sizeof(int), &nFormed, CkReduction::sum_int, cb);

}
#endif
