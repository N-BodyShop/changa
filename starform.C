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
    //bUseStoch=1;
    //prmAddParam(prm, "bUseStoch", paramInt, &bUseStoch,
    //    sizeof(double), "usestoch", "<Enable stochastic IMF = 1>");
    iStarFormRung = 0;
    prmAddParam(prm,"iStarFormRung", paramInt, &iStarFormRung,
		sizeof(int), "iStarFormRung", "<Star Formation Rung> = 0");
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
    if (dInitStarMass > 0) {
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

    bUseStoch = param.feedback->sn.bUseStoch;
    dSecUnit = param.dSecUnit;
    dMsolUnit = param.dMsolUnit;
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
    double startTime = CkWallTimer();
    //
    // Need to build tree since we just did a drift.
#ifdef PUSH_GRAVITY
    treeProxy.buildTree(bucketSize, CkCallbackResumeThread(),true);
#else
    treeProxy.buildTree(bucketSize, CkCallbackResumeThread());
#endif
    double tTB = CkWallTimer() - startTime;
    timings[PHASE_FEEDBACK].tTBuild += tTB;
    DensitySmoothParams pDen(TYPE_GAS, 0);
    double dfBall2OverSoft2 = 4.0*param.dhMinOverSoft*param.dhMinOverSoft;
    treeProxy.startSmooth(&pDen, 1, param.nSmooth, dfBall2OverSoft2,
			  CkCallbackResumeThread());

    CkReductionMsg *msgCounts;
    treeProxy.FormStars(*(param.stfm), dTime, dDelta,
			pow(csmTime2Exp(param.csm, dTime), 3.0),
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
			  double dCosmoFac, const CkCallback& cb) 
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
						   &TempForm);
	    
	    if(starp != NULL) {
		nFormed++;
		dMassFormed += starp->mass;
		newParticle(starp);
                p = &myParticles[i]; // newParticle can change pointers
		CmiLock(dm->lockStarLog);
                // Record current spot in seTab
                iSeTab.push_back(dm->starLog->seTab.size());
		dm->starLog->seTab.push_back(StarLogEvent(starp,dCosmoFac,TempForm));
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
				double dCosmoFac, double *T) 
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
    double dMprob = 1.0 - exp(-dCStar*dTimeStarForm/tform);
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

    /* Stochastically populating from IMF up to target mass 
    *  You end up with an array of individual high mass stars and
    * a normalization constant for the continuous, low mass IMF
    */
    if(bUseStoch){
        /* Setting all high mass stars to default (0) */
        for(int i=0;i<12;i++){
            starp->rgfHMStars(i)=0.0;
            }
        /* Stochastically populate star particle
        * Keep running tally of total mass, but only keep high mass stars. Stars
        * are drawn from the IMF until the running tally exceeds the dDeltaM, then
        * the last star is kept only if the total is closer to dDeltaM with it
        *
        * DrawStar returns stars in units of Msol. HMStars is filled with masses in
        * units of Msol, but depending on context this needs to be converted to
        * system units
        */
        int iArrayLoc = 0;
        double mass_tally = 0.0;
        // Sum of only HMStars, used to find fLowNorm
        double dSumHMStars=0.0;
        while(1){
            //srand?
            /* (Should be very rare) If we form more HMStars than elements in array,
            * wipe everything and start over
            */
            if(iArrayLoc>=12){
                mass_tally=0.0;
                for(int i=0;i<12;i++) starp->rgfHMStars(i)=0;
                iArrayLoc=0;
                dSumHMStars=0.0;
            }
            double num = (rand()/((double) RAND_MAX));
            /* DrawStar returns a mass in Msol, need to convert to system units */
            double new_star = imf->DrawStar(num);
            //CkPrintf("new_star is %f\n",new_star);
            double new_star_unit = new_star/dMsolUnit;
            double test_mass = mass_tally + new_star_unit;
            if(test_mass < dDeltaM){
                mass_tally += new_star_unit;
                if(new_star > 8.0){
                    //CkPrintf("HMStar with mass: %f\n",new_star);
                    starp->rgfHMStars(iArrayLoc)=new_star;
                    dSumHMStars += new_star_unit;
                    iArrayLoc+=1;
                }
            } else if(fabs(dDeltaM - test_mass) < fabs(dDeltaM - mass_tally) ) {
                mass_tally += new_star;
                if(new_star > 8.0){
                    starp->rgfHMStars(iArrayLoc)=new_star;
                    dSumHMStars += new_star_unit;
                }
                break;
            } else break;
        }
        double dTotLowMass=imf->CumMass(0.0)-imf->CumMass(8.0);
        starp->fLowNorm()=(dDeltaM-dSumHMStars)/dTotLowMass;
        //CkPrintf("fLowNorm is %f\n",starp->fLowNorm());
    }

    CkPrintf("Array of HMStars: ");
    for(int i=0;i<12;i++) CkPrintf("%f ", starp->rgfHMStars(i));
    CkPrintf("\n");
    if(!bUseStoch) CkPrintf("bUseStoch = %f\n");
    p->mass -= dDeltaM;
    CkAssert(p->mass >= 0.0);
    starp->mass = dDeltaM;
    starp->fMassForm() = dDeltaM;
    if(p->mass < dMinGasMass) {
	deleteParticle(p);
	}

    /* NB: It is important that the star inherit special
       properties of the gas particle such as being a target
       for movies or other tracing
       Thus: Do not remove all the TYPE properties -- just
       the gas specific ones */
    TYPEReset(starp, TYPE_GAS|TYPE_NbrOfACTIVE);
    TYPESet(starp, TYPE_STAR) ;
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
	    }
	xdr_destroy(&xdrs);
	int result = CmiFclose(outfile);
	if(result != 0)
            CkAbort("Bad close of starlog");
	seTab.clear();
	nOrdered = 0;
	}
    }
