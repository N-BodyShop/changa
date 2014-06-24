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
    CkPrintf("Star Formation Calculated, Wallclock %f secs\n",
	     CkWallTimer() - startTime);
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
    
    for(unsigned int i = 1; i <= myNumParticles; ++i) {
	GravityParticle *p = &myParticles[i];
	if(p->isGas()) {
	    GravityParticle *starp = stfm.FormStar(p, dm->Cool, dTime,
						   dDelta, dCosmoFac, 
						   &TempForm);
	    
	    if(starp != NULL) {
		nFormed++;
		dMassFormed += starp->mass;
		newParticle(starp);
		CmiLock(dm->lockStarLog);
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
    	*T = CoolCodeEnergyToTemperature(Cool, &p->CoolParticle(),
				    	p->u(), p->fMetals());
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
	    FILE *fpLog = fopen(stLogFile.c_str(),"r");
	    XDR xdrs;
	    
	    CkAssert(fpLog != NULL);
	    xdrstdio_create(&xdrs,fpLog,XDR_DECODE);
	    xdr_int(&xdrs, &iSize);
	    CkAssert(iSize == sizeof(StarLogEvent));
	    xdr_destroy(&xdrs);
	    fclose(fpLog);
	    } else {
	    CkAbort("Simulation restarting with star formation, but starlog file not found");
	    }
	} else {
	FILE *fpLog = fopen(stLogFile.c_str(),"w");
	XDR xdrs;

	CkAssert(fpLog != NULL);
	xdrstdio_create(&xdrs,fpLog,XDR_ENCODE);
	iSize = sizeof(StarLogEvent);
	xdr_int(&xdrs, &iSize);
	xdr_destroy(&xdrs);
	fclose(fpLog);
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
	
	outfile = fopen(fileName.c_str(), "a");
	CkAssert(outfile != NULL);
	xdrstdio_create(&xdrs,outfile,XDR_ENCODE);

	CkAssert(seTab.size() == nOrdered);
    
	for(int iLog = 0; iLog < seTab.size(); iLog++){
	    StarLogEvent *pSfEv = &(seTab[iLog]);
	    xdr_int(&xdrs, &(pSfEv->iOrdStar));
	    xdr_int(&xdrs, &(pSfEv->iOrdGas));
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
	int result = fclose(outfile);
	if(result != 0)
	    ckerr << "Bad close of starlog" << endl;
	CkAssert(result == 0);
	seTab.clear();
	nOrdered = 0;
	}
    }
