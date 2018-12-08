/*
 * Star forming module originally written for GASOLINE
 * Modified for inclusion in ChaNGa
 *
 * Algorithm was originally implemented in TREESPH but is heavily
 * modified.  Contributers include Neal Katz, Eric Hayashi, Greg Stinson
 */
#include <math.h>
#include <sys/stat.h>
#include "ParallelGravity.h"
#include "starform.h"
#include "smooth.h"
#include "Sph.h"

///
/// @brief initialize parameters for star formation
///

void Stfm::AddParams(PRM prm) 
{
    dOverDenMin = 20.0;
    prmAddParam(prm,"dOverDenMin", paramDouble, &dOverDenMin,
		    sizeof(double), "stODmin",
		    "<Minimum overdensity for forming stars> = 2");
    dPhysDenMin = 0.1;
    prmAddParam(prm,"dPhysDenMin", paramDouble, &dPhysDenMin,
		sizeof(double), "stPDmin",
		"<Minimum physical density for forming stars (atoms/cc)> = .1");
    dStarEff = .3333;
    prmAddParam(prm,"dStarEff", paramDouble, &dStarEff,
		sizeof(double), "stEff",
		"<Fraction of gas converted into stars per timestep> = .3333");
    dInitStarMass = 0.0;
    prmAddParam(prm,"dInitStarMass", paramDouble, &dInitStarMass,
		sizeof(double), "stm0", "<Initial star mass> = 0");
    dMinGasMass = 0.0;
    prmAddParam(prm,"dMinGasMass", paramDouble, &dMinGasMass,
		sizeof(double), "stMinGas",
		"<Minimum mass of a gas particle> = 0.0");
    dMaxStarMass = 0.0;
    prmAddParam(prm,"dMaxStarMass", paramDouble, &dMaxStarMass,
		sizeof(double), "stMaxStarMass",
		"<Maximum amount of star mass a hybrid particle can contain = 0.0");
    dCStar = 0.05;
    prmAddParam(prm,"dCStar", paramDouble, &dCStar,
		sizeof(double), "stCStar",
		"<Star formation coefficient> = 0.1");
    dTempMax = 1.5e4;
    prmAddParam(prm,"dTempMax", paramDouble, &dTempMax,
		sizeof(double), "stTempMax",
		"<Maximum temperature at which star formation occurs> = 1.5e4");
    dSoftMin = 1.0;
    prmAddParam(prm,"dSoftMin", paramDouble, &dSoftMin,
		sizeof(double), "stSoftMin",
		"<Minimum softening for star formation> = 0.0");
    dDeltaStarForm = 1e6;
    prmAddParam(prm,"dDeltaStarForm", paramDouble, &dDeltaStarForm,
		sizeof(double), "dDeltaStarForm",
		"<Minimum SF timestep in years> = 1e6");
    iStarFormRung = 0;
    prmAddParam(prm,"iStarFormRung", paramInt, &iStarFormRung,
		sizeof(int), "iStarFormRung", "<Star Formation Rung> = 0");
#ifdef COOLING_MOLECULARH
    dStarFormEfficiencyH2 = 1.0; /* Set such that SF depends on H2 abundance. Usually 1 to make SF efficiency based off of H2 abundance.  Set to zero for standard, non-H2 SF recipe.*/
    prmAddParam(prm,"dStarFormEfficiencyH2", paramDouble, &dStarFormEfficiencyH2,
		sizeof(double), "dStarFormEfficiencyH2",
		"<Star Formation efficiency as a function of H2> = 0");
#endif
    // Blackhole formation parameters
    bBHForm = 0;
    prmAddParam(prm,"bBHForm",paramBool,&bBHForm,
		sizeof(int),"stBhForm","<enable seed BH formation> = 0");
    dBHFormProb = 0.0;
    prmAddParam(prm,"dBHFormProb",paramDouble,&dBHFormProb, sizeof(double),
		"stBhFormProb", "<seed BH formation probability> = 0.0");
    dInitBHMass = 0.0;
    prmAddParam(prm,"dInitBHMass",paramDouble, &dInitBHMass, sizeof(double),
		"bhm0", "<initial BH mass> = 0.0");
    iRandomSeed = 1;
    prmAddParam(prm,"iRandomSeed", paramInt, &iRandomSeed, sizeof(int),
		"iRand", "<Star formation random Seed> = 1");
    }

/*
 * Verify that the set parameters are OK.
 */
void Stfm::CheckParams(PRM prm, Parameters &param) 
{
    if(param.bStarForm) {
	CkAssert(prmSpecified(prm, "dMsolUnit")
		 && prmSpecified(prm, "dKpcUnit"));
	if(!param.bGasCooling)
	    CkError("Warning: You are not running a cooling EOS with starformation\n");
	param.bDoGas = 1;
	}
    bGasCooling = param.bGasCooling;
    CkAssert((dStarEff > 0.0 && dStarEff < 1.0)
	      || dInitStarMass > 0.0);
    if (dInitStarMass > 0 && dMinGasMass <= 0) {
	/* Only allow 10% underweight star particles */
	dMinGasMass = 0.9*dInitStarMass;
	}
    else
	CkAssert(dMinGasMass > 0.0);
    if (!param.csm->bComove && !prmSpecified(prm, "dOverDenMin")) {
	/*
	 * Overdensity criterion makes no sense for non-comoving simulations
	 */
	dOverDenMin = 0.0;
	}
		
#include "physconst.h"

    dSecUnit = param.dSecUnit;
    dGmPerCcUnit = param.dGmPerCcUnit;
    dGmUnit = param.dMsolUnit*MSOLG;
    dErgUnit = GCGS*pow(param.dMsolUnit*MSOLG, 2.0)
	/(param.dKpcUnit*KPCCM);
    /* convert to system units */
    dPhysDenMin *= MHYDR/dGmPerCcUnit;
    dDeltaStarForm *= SECONDSPERYEAR/param.dSecUnit;

    double testDelta;
    for(testDelta = param.dDelta;
	testDelta >= dDeltaStarForm && dDeltaStarForm > 0.0;
	testDelta *= 0.5 ){
	if ( !(prmSpecified(prm,"iStarFormRung")) )
	    iStarFormRung++;
	}
    if ( testDelta <= param.dDelta ){
	CkPrintf("dDeltaStarForm (set): %g, effectively: %g = %g yrs, iStarFormRung: %i\n",
		 dDeltaStarForm, testDelta,
		 testDelta*param.dSecUnit/SECONDSPERYEAR,
		 iStarFormRung );
	dDeltaStarForm = testDelta;
	}
    else if ( dDeltaStarForm == 0.0 ) {
	CkPrintf(" dDeltaStarForm (set): %g, effectively: 0.0 = 0.0 yrs, iStarFormRung: maxRung\n",
		 dDeltaStarForm );
	
	iStarFormRung = param.iMaxRung;
	}
    else {
	CkPrintf(" dDeltaStarForm (set): %g, effectively:  NO STARS WILL FORM\n",
		 dDeltaStarForm);
	}
    }

void TreePiece::initRand(int iRand, const CkCallback &cb)  
{
   srand(iRand + thisIndex);  // make each piece unique
   contribute(cb);
}

///
/// form stars main method
///
void Main::FormStars(double dTime, double dDelta) 
{
    if(verbosity)
	CkPrintf("Form Stars ... ");
    //
    // Need to build tree since we just did a drift.
    buildTree(PHASE_FEEDBACK);
    double startTime = CkWallTimer();
    DensitySmoothParams pDen(TYPE_GAS, 0);
    double dfBall2OverSoft2 = 4.0*param.dhMinOverSoft*param.dhMinOverSoft;
    treeProxy.startSmooth(&pDen, 1, param.nSmooth, dfBall2OverSoft2,
			  CkCallbackResumeThread());

    CkReductionMsg *msgCounts;
    treeProxy.FormStars(*(param.stfm), dTime, dDelta,
			pow(csmTime2Exp(param.csm, dTime), 3.0), 0,
			CkCallbackResumeThread((void*&)msgCounts));
    int *dCounts = (int *)msgCounts->getData();
    
    int nDelGas = dCounts[1];
    if(verbosity)
	CkPrintf("%d Stars formed, %d gas deleted\n", dCounts[0], dCounts[1]);
    delete msgCounts;
    
    if(nDelGas > 0) {
	if(verbosity)
	    CkPrintf("Distribute Deleted gas\n");
	DistDeletedGasSmoothParams pDGas(TYPE_GAS, 0);
	treeProxy.startReSmooth(&pDGas, CkCallbackResumeThread());

	}

    treeProxy.finishNodeCache(CkCallbackResumeThread());
#ifdef CUDA
    // We didn't do gravity where the registered TreePieces on the
    // DataManager normally get cleared.  Clear them here instead.
    dMProxy.clearRegisteredPieces(CkCallbackResumeThread());
#endif

    addDelParticles();
    double tSF = CkWallTimer() - startTime;
    timings[PHASE_FEEDBACK].tGrav += tSF; // Overload meaning of tGrav
                                          // for SF/Feedback
    CkPrintf("Star Formation Calculated, Wallclock %f secs\n", tSF);
    }

///
/// processor specific method for star formation
/// 

void TreePiece::FormStars(Stfm stfm, double dTime,  double dDelta,
			  double dCosmoFac, double H2FractionForm, const CkCallback& cb) 
{
    int nFormed = 0;
    int nDeleted = 0;
    double dMassFormed = 0.0;
    double TempForm;
    
    // clear indices into starlog table
    iSeTab.clear();
    int myNPartTmp = myNumParticles;  // this could change if newParticle
                                      // is called.
    for(unsigned int i = 1; i <= myNPartTmp; ++i) {
	GravityParticle *p = &myParticles[i];
	if(p->isGas()) {
	    GravityParticle *starp = stfm.FormStar(p, dm->Cool, dTime,
						   dDelta, dCosmoFac, 
						   H2FractionForm, &TempForm);
	    
	    if(starp != NULL) {
		nFormed++;
		dMassFormed += starp->mass;
		newParticle(starp);
                p = &myParticles[i]; // newParticle can change pointers
		CmiLock(dm->lockStarLog);
                // Record current spot in seTab
                iSeTab.push_back(dm->starLog->seTab.size());
		dm->starLog->seTab.push_back(StarLogEvent(starp,dCosmoFac,H2FractionForm,TempForm));
		CmiUnlock(dm->lockStarLog);
		delete (extraStarData *)starp->extraData;
		delete starp;
		if(TYPETest(p, TYPE_DELETED))
		    nDeleted++;
		}
	    }
	}
    
    int counts[2];
    counts[0] = nFormed;
    counts[1] = nDeleted;
    contribute(2*sizeof(int), counts, CkReduction::sum_int, cb);
    }

/*
     taken from TREESPH and modified greatly.
     Uses the following formula for the star formation rate:

              d(ln(rhostar))/dt=cstar*rhogas/tdyn

*/
GravityParticle *Stfm::FormStar(GravityParticle *p,  COOL *Cool, double dTime,
				double dDelta,  // drift timestep
				double dCosmoFac, double H2FractionForm, double *T) 
{
    /*
     * Determine dynamical time.
     */
    double tdyn = 1.0/sqrt(4.0*M_PI*p->fDensity/dCosmoFac);
#ifndef COOLING_NONE
    if(bGasCooling)
#ifdef COOLING_GRACKLE
    	*T = CoolCodeEnergyToTemperature(Cool, &p->CoolParticle(),
                                         p->u(), p->fDensity, p->fMetals());
#else
    	*T = CoolCodeEnergyToTemperature(Cool, &p->CoolParticle(),
				    	p->u(), p->fMetals());
#endif
    else
	*T = 0.0;
#else
    *T = 0.0;
#endif
    /*
     * Determine if this particle satisfies all conditions.
     */
    if(*T > dTempMax)
	return NULL;

    if(p->fDensity < dOverDenMin || p->fDensity/dCosmoFac < dPhysDenMin)
	return NULL;

    double tform = tdyn;
    double dTimeStarForm = (dDelta > dDeltaStarForm ? dDelta : dDeltaStarForm);

#ifdef COOLING_MOLECULARH
    double yH;
    if (p->fMetals() <= 0.1) yH = 1.0 - 4.0*((0.236 + 2.1*(p->fMetals()))/4.0) - p->fMetals();
    else yH = 1.0 - 4.0*((-0.446*(p->fMetals() - 0.1)/0.9 + 0.446)/4.0) - p->fMetals();
    double columnL = sqrt(0.25)*p->fBall;/*correlation length used for H2 shielding, CC*/
    double dMprob;
    if (dStarFormEfficiencyH2 == 0) dMprob  = 1.0 - exp(-dCStar*dTimeStarForm/tform);
    else dMprob = 1.0 - exp(-dCStar*dTimeStarForm/tform*
    			    dStarFormEfficiencyH2*(2.0*(p->CoolParticle().f_H2)/yH));    
    H2FractionForm = 1.0*p->CoolParticle().f_H2/yH;
#else /* COOLING_MOLECULARH */ 
    double dMprob = 1.0 - exp(-dCStar*dTimeStarForm/tform);
    H2FractionForm = 0;   
#endif /* COOLING_MOLECULARH */ 

    /*
     * Decrement mass of particle.
     */
    double dDeltaM;
    if (dInitStarMass > 0.0) 
	dDeltaM = dInitStarMass;
    else 
	dDeltaM = p->mass*dStarEff;

    /* No negative or very tiny masses please! */
    if ( (dDeltaM > p->mass) ) dDeltaM = p->mass;

    if(dMprob*p->mass < dDeltaM*(rand()/((double) RAND_MAX)))
	return NULL;

    /* 
     * Note on number of stars formed:
     * n = log(dMinGasMass/dInitMass)/log(1-dStarEff) = max
     * no. stars 
     * formed per gas particle, e.g. if min gas mass = 10%
     * initial mass,
     * dStarEff = 1/3, max no. stars formed = 6 (round up so
     * gas mass goes below min gas mass)
     */

    GravityParticle *starp = new GravityParticle();
    *starp = StarFromGasParticle(p); /* grab copy before
				       possible deletion */
    /*
     * form star
     */

    starp->fTimeForm() = dTime;
    /* iOrder gets reassigned in NewParticle() */
    starp->iGasOrder() = starp->iOrder;
    /* Seed BH Formation JMB 1/19/09*/
    int newbh = 0;  /* BH tracker */
    if (bBHForm == 1 && starp->fStarMetals() <= 1.0e-6
	&& dBHFormProb > (rand()/((double) RAND_MAX ))) {
	starp->fTimeForm() = -1.0*starp->fTimeForm();
	newbh = 1;
	/* Decrement mass of particle.*/
	if (dInitBHMass > 0) {
	    dDeltaM = dInitBHMass;  /* reassigning dDeltaM to be
				       initBHmass JMB 6/16/09 */
	    /* No negative or very tiny masses please! */
	    if (dDeltaM > p->mass) dDeltaM = p->mass;
	    }
	}


    p->mass -= dDeltaM;
    CkAssert(p->mass >= 0.0);
    starp->mass = dDeltaM;
    starp->fMassForm() = dDeltaM;
    if(p->mass < dMinGasMass) {
	deleteParticle(p);
	}
    /*#ifdef COOLING_MOLECULARH*/ /*Initialize LW radiation for a star particle of that mass and 10^7 (the minimum) years old*/
#ifdef LYMAN_WERNER
    double dAgelog = 7,
      a0 = -84550.812,
      a1 =  54346.066,
      a2 = -13934.144,
      a3 =  1782.1741,
      a4 = -113.68717,
      a5 =  2.8930795;
    starp->dStarLymanWerner() = a0
      + a1*dAgelog
      + a2*dAgelog*dAgelog
      + a3*dAgelog*dAgelog*dAgelog
      + a4*dAgelog*dAgelog*dAgelog*dAgelog
      + a5*dAgelog*dAgelog*dAgelog*dAgelog*dAgelog + log10(dDeltaM);
#endif /*LYMAN_WERNER*/
    /* #endif*/ /*COOLING_MOLECULARH*/

    /* NB: It is important that the star inherit special
       properties of the gas particle such as being a target
       for movies or other tracing
       Thus: Do not remove all the TYPE properties -- just
       the gas specific ones */
    TYPEReset(starp, TYPE_GAS|TYPE_NbrOfACTIVE);
    if(newbh == 0) TYPESet(starp, TYPE_STAR) ; /* if it's a BH make
						   it a SINK  JMB  */
    else TYPESet(starp, TYPE_SINK);
    return starp;
    }

void Main::initStarLog(){
    struct stat statbuf;
    int iSize;
    std::string stLogFile = std::string(param.achOutName) + ".starlog";

    dMProxy.initStarLog(stLogFile,CkCallbackResumeThread());

    if(bIsRestarting) {
	if(!stat(stLogFile.c_str(), &statbuf)) {
	    /* file exists, check number */
	    FILE *fpLog = CmiFopen(stLogFile.c_str(),"r");
	    XDR xdrs;
            if(fpLog == NULL)
                CkAbort("Bad open of starlog file on restart");

	    xdrstdio_create(&xdrs,fpLog,XDR_DECODE);
	    xdr_int(&xdrs, &iSize);
            if(iSize != sizeof(StarLogEvent))
                CkAbort("starlog file format mismatch");
	    xdr_destroy(&xdrs);
	    CmiFclose(fpLog);
	    } else {
	    CkAbort("Simulation restarting with star formation, but starlog file not found");
	    }
	} else {
	FILE *fpLog = CmiFopen(stLogFile.c_str(),"w");
	XDR xdrs;

        if(fpLog == NULL) 
            CkAbort("Can't create starlog file");
	xdrstdio_create(&xdrs,fpLog,XDR_ENCODE);
	iSize = sizeof(StarLogEvent);
	xdr_int(&xdrs, &iSize);
	xdr_destroy(&xdrs);
	CmiFclose(fpLog);
	}

    }

void DataManager::initStarLog(std::string _fileName, const CkCallback &cb) {
    starLog->fileName = _fileName;
    contribute(cb);
    }

/// \brief flush starlog table to disk.
void TreePiece::flushStarLog(const CkCallback& cb) {

    if(verbosity > 3)
	ckout << "TreePiece " << thisIndex << ": Writing output to disk" << endl;
    if (dm == NULL)
        dm = (DataManager*)CkLocalNodeBranch(dataManagerID);
    CmiLock(dm->lockStarLog);
    dm->starLog->flush();
    CmiUnlock(dm->lockStarLog);

    if(thisIndex!=(int)numTreePieces-1) {
	pieces[thisIndex + 1].flushStarLog(cb);
	return;
	}

    cb.send(); // We are done.
    return;
    }

void StarLog::flush(void) {
    if (seTab.size() > 0) {
	FILE* outfile;
	XDR xdrs;
	
	outfile = CmiFopen(fileName.c_str(), "a");
	CkAssert(outfile != NULL);
	xdrstdio_create(&xdrs,outfile,XDR_ENCODE);

	CkAssert(seTab.size() == nOrdered);
    
	for(int iLog = 0; iLog < seTab.size(); iLog++){
	    StarLogEvent *pSfEv = &(seTab[iLog]);
	    xdr_template(&xdrs, &(pSfEv->iOrdStar));
	    xdr_template(&xdrs, &(pSfEv->iOrdGas));
	    xdr_double(&xdrs, &(pSfEv->timeForm));
	    xdr_double(&xdrs, &(pSfEv->rForm[0]));
	    xdr_double(&xdrs, &(pSfEv->rForm[1]));
	    xdr_double(&xdrs, &(pSfEv->rForm[2]));
	    xdr_double(&xdrs, &(pSfEv->vForm[0]));
	    xdr_double(&xdrs, &(pSfEv->vForm[1]));
	    xdr_double(&xdrs, &(pSfEv->vForm[2]));
	    xdr_double(&xdrs, &(pSfEv->massForm));
	    xdr_double(&xdrs, &(pSfEv->rhoForm));
	    xdr_double(&xdrs, &(pSfEv->TForm));
#ifdef COOLING_MOLECULARH 
	    xdr_double(&xdrs, &(pSfEv->H2FracForm));
#endif
	    }
	xdr_destroy(&xdrs);
	int result = CmiFclose(outfile);
	if(result != 0)
            CkAbort("Bad close of starlog");
	seTab.clear();
	nOrdered = 0;
	}
    }
