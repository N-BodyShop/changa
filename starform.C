/*
 * Star forming module originally written for GASOLINE
 * Modified for inclusion in ChaNGa
 *
 * Algorithm was originally implemented in TREESPH but is heavily
 * modified.  Contributers include Neal Katz, Eric Hayashi, Greg Stinson
 */
#include <math.h>
#include "ParallelGravity.h"
#include "starform.h"
#include "smooth.h"

///
/// @brief initialize parameters for star formation
///

void StfmAddParams(StfmParam *stfm, PRM prm) 
{
    stfm->dOverDenMin = 2.0;
    prmAddParam(prm,"dOverDenMin", paramDouble, &stfm->dOverDenMin,
		    sizeof(double), "stODmin",
		    "<Minimum overdensity for forming stars> = 2");
    stfm->dPhysDenMin = 0.1;
    prmAddParam(prm,"dPhysDenMin", paramDouble, &stfm->dPhysDenMin,
		sizeof(double), "stPDmin",
		"<Minimum physical density for forming stars (atoms/cc)> = .1");
    stfm->dStarEff = .3333;
    prmAddParam(prm,"dStarEff", paramDouble, &stfm->dStarEff,
		sizeof(double), "stEff",
		"<Fraction of gas converted into stars per timestep> = .3333");
    stfm->dInitStarMass = 0.0;
    prmAddParam(prm,"dInitStarMass", paramDouble, &stfm->dInitStarMass,
		sizeof(double), "stm0", "<Initial star mass> = 0");
    stfm->dMinGasMass = 0.0;
    prmAddParam(prm,"dMinGasMass", paramDouble, &stfm->dMinGasMass,
		sizeof(double), "stMinGas",
		"<Minimum mass of a gas particle> = 0.0");
    stfm->dMaxStarMass = 0.0;
    prmAddParam(prm,"dMaxStarMass", paramDouble, &stfm->dMaxStarMass,
		sizeof(double), "stMaxStarMass",
		"<Maximum amount of star mass a hybrid particle can contain = 0.0");
    stfm->dCStar = 0.05;
    prmAddParam(prm,"dCStar", paramDouble, &stfm->dCStar,
		sizeof(double), "stCStar",
		"<Star formation coefficient> = 0.1");
    stfm->dTempMax = 1.5e4;
    prmAddParam(prm,"dTempMax", paramDouble, &stfm->dTempMax,
		sizeof(double), "stTempMax",
		"<Maximum temperature at which star formation occurs> = 1.5e4");
    stfm->dSoftMin = 1.0;
    prmAddParam(prm,"dSoftMin", paramDouble, &stfm->dSoftMin,
		sizeof(double), "stSoftMin",
		"<Minimum softening for star formation> = 0.0");
    stfm->dDeltaStarForm = 1e6;
    prmAddParam(prm,"dDeltaStarForm", paramDouble, &stfm->dDeltaStarForm,
		sizeof(double), "dDeltaStarForm",
		"<Minimum SF timestep in years> = 1e6");
    stfm->iStarFormRung = 0;
    prmAddParam(prm,"iStarFormRung", paramInt, &stfm->iStarFormRung,
		sizeof(int), "iStarFormRung", "<Star Formation Rung> = 0");
    }

/*
 * Verify that the set parameters are OK.
 */
void StfmCheckParams(StfmParam *stfm, PRM prm, Parameters &param) 
{
    if(param.bStarForm) {
	CkAssert(prmSpecified(prm, "dMsolUnit")
		 && prmSpecified(prm, "dKpcUnit"));
	if(!param.bGasCooling)
	    CkError("Warning: You are not running a cooling EOS with starformation\n");
	}
    CkAssert((stfm->dStarEff > 0.0 && stfm->dStarEff < 1.0)
	      || stfm->dInitStarMass > 0.0);
    if (stfm->dInitStarMass > 0) {
	/* Only allow 10% underweight star particles */
	stfm->dMinGasMass = 0.9*stfm->dInitStarMass;
	}
    else
	CkAssert(stfm->dMinGasMass > 0.0);
    if (!param.csm->bComove && !prmSpecified(prm, "dOverDenMin")) {
	/*
	 * Overdensity criterion makes no sense for non-comoving simulations
	 */
	stfm->dOverDenMin = 0.0;
	}
		
#include "physconst.h"

    stfm->dSecUnit = param.dSecUnit;
    stfm->dGmPerCcUnit = param.dGmPerCcUnit;
    stfm->dGmUnit = param.dMsolUnit*MSOLG;
    stfm->dErgUnit = GCGS*pow(param.dMsolUnit*MSOLG, 2.0)
	/(param.dKpcUnit*KPCCM);
    /* convert to system units */
    stfm->dPhysDenMin *= MHYDR/stfm->dGmPerCcUnit;
    stfm->dDeltaStarForm *= SECONDSPERYEAR/param.dSecUnit;

    double testDelta;
    for(testDelta = param.dDelta;
	testDelta >= stfm->dDeltaStarForm && stfm->dDeltaStarForm > 0.0;
	testDelta *= 0.5 ){
	if ( !(prmSpecified(prm,"iStarFormRung")) )
	    stfm->iStarFormRung++;
	}
    if ( testDelta <= param.dDelta ){
	CkPrintf("dDeltaStarForm (set): %g, effectively: %g = %g yrs, iStarFormRung: %i",
		 stfm->dDeltaStarForm, testDelta,
		 testDelta*param.dSecUnit/SECONDSPERYEAR,
		 stfm->iStarFormRung );
	stfm->dDeltaStarForm = testDelta;
	}
    else if ( stfm->dDeltaStarForm == 0.0 ) {
	CkPrintf(" dDeltaStarForm (set): %g, effectively: 0.0 = 0.0 yrs, iStarFormRung: maxRung",
		 stfm->dDeltaStarForm );
	
	stfm->iStarFormRung = param.iMaxRung;
	}
    else {
	CkPrintf(" dDeltaStarForm (set): %g, effectively:  NO STARS WILL FORM",
		 stfm->dDeltaStarForm);
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
    // XXX need to check whether a treebuild needs the domain
    // decomposition.  If not, this could be avoided.
    //
    double tolerance = 0.01;	// tolerance for domain decomposition
    sorter.startSorting(dataManagerID, tolerance,
                        CkCallbackResumeThread(), true);
    treeProxy.buildTree(bucketSize, CkCallbackResumeThread());
    DensitySmoothParams pDen(TYPE_GAS, 0);
    treeProxy.startIterationSmooth(&pDen, CkCallbackResumeThread());

    CkReductionMsg *msgCounts;
    treeProxy.FormStars(*(param.stfm), dTime,
			pow(csmTime2Exp(param.csm, dTime), 3.0),
			CkCallbackResumeThread((void*&)msgCounts));
    double *dCounts = (double *)msgCounts->getData();
    
    if(verbosity)
	CkPrintf("%d Stars formed, %d gas deleted\n", dCounts[0], dCounts[1]);
    delete msgCounts;
    
    DistDeletedGasSmoothParams pDGas(TYPE_GAS, 0);
    treeProxy.startIterationSmooth(&pDGas, CkCallbackResumeThread());

    addDelParticles();
    CkPrintf("Star Formation Calculated, Wallclock %f secs\n",
	     CkWallTimer() - startTime);
    }

///
/// processor specific method for star formation
/// 
/*
     taken from TREESPH and modified greatly.
     Uses the following formula for the star formation rate:

              d(ln(rhostar))/dt=cstar*rhogas/tdyn

*/

void TreePiece::FormStars(StfmParam stfm, double dTime, double dCosmoFac,
			  const CkCallback& cb) 
{
    int nFormed = 0;
    int nDeleted = 0;
    double dMassFormed = 0.0;
    
    for(unsigned int i = 1; i <= myNumParticles; ++i) {
	GravityParticle *p = &myParticles[i];
	if(p->isGas()) {
	    double T;
	    double tcool;
	    /*
	     * Determine dynamical time.
	     */
	    double tdyn = 1.0/sqrt(4.0*M_PI*p->fDensity/dCosmoFac);
#ifndef COOLING_NONE
	    T = CoolCodeEnergyToTemperature(dm->Cool, &p->CoolParticle(),
					    p->u(), p->fMetals());
	    double r[3];  // For conversion to C
	    p->position.array_form(r);
	    tcool = p->u()
		/CoolEdotInstantCode(dm->Cool, &p->CoolParticle(), p->u(),
				     p->fDensity, p->fMetals(), r);
#endif
	    /*
	     * Determine if this particle satisfies all conditions.
	     */
	    if(p->fDensity < stfm.dOverDenMin ||
	       p->fDensity/dCosmoFac < stfm.dPhysDenMin)
		continue;
	    
	    double tform;
	    if(tcool < 0.0 || tdyn > tcool || T < stfm.dTempMax)
		tform = tdyn;
	    else
		tform = tcool;

	    double dMprob = 1.0 - exp(-stfm.dCStar*stfm.dDeltaStarForm/tform);

	    /*
	     * Decrement mass of particle.
	     */
	    double dDeltaM;
	    if (stfm.dInitStarMass > 0.0) 
		dDeltaM = stfm.dInitStarMass;
	    else 
		dDeltaM = p->mass*stfm.dStarEff;

	    /* No negative or very tiny masses please! */
	    if ( (dDeltaM > p->mass) ) dDeltaM = p->mass;

	    if(dMprob*p->mass < dDeltaM*(rand()/((double) RAND_MAX)))
		continue;

	    /* 
	     * Note on number of stars formed:
	     * n = log(dMinGasMass/dInitMass)/log(1-dStarEff) = max
	     * no. stars 
	     * formed per gas particle, e.g. if min gas mass = 10%
	     * initial mass,
	     * dStarEff = 1/3, max no. stars formed = 6 (round up so
	     * gas mass goes below min gas mass)
	     */

	    GravityParticle
		starp = StarFromGasParticle(p); /* grab copy before
					       possible deletion */
	    /*
	     * form star
	     */

	    starp.fTimeForm() = dTime;
	    /* iOrder gets reassigned in NewParticle() */
	    starp.iGasOrder() = starp.iOrder;
		
	    p->mass -= dDeltaM;
	    CkAssert(p->mass >= 0.0);
	    starp.mass = dDeltaM;
	    starp.fMassForm() = dDeltaM;
	    if(p->mass < stfm.dMinGasMass) {
		nDeleted++;
		deleteParticle(p);
		}
	    p->mass -= dDeltaM;
	    CkAssert(p->mass >= 0.0);
	    starp.mass = dDeltaM;
	    starp.fMassForm() = dDeltaM;

	    /* NB: It is important that the star inherit special
	       properties of the gas particle such as being a target
	       for movies or other tracing
	       Thus: Do not remove all the TYPE properties -- just
	       the gas specific ones */
	    TYPEReset(&starp, TYPE_GAS|TYPE_NbrOfACTIVE);
	    TYPESet(&starp, TYPE_STAR) ;
    
	    nFormed++;
	    dMassFormed += dDeltaM;
	    newParticle(&starp);
	    delete (extraStarData *)starp.extraData;
	    }
	}
    
    int counts[2];
    counts[0] = nFormed;
    counts[1] = nDeleted;
    contribute(2*sizeof(int), counts, CkReduction::sum_int, cb);
    }

/*
 * Methods to distribute Deleted gas
 */
int DistDeletedGasSmoothParams::isSmoothActive(GravityParticle *p) 
{
    return (TYPETest(p, TYPE_DELETED) && TYPETest(p, iType));
    }

void DistDeletedGasSmoothParams::initSmoothCache(GravityParticle *p) 
{
    if(!TYPETest(p, TYPE_DELETED)) {
	/*
	 * Zero out accumulated quantities.
	 */
	p->mass = 0;
	p->velocity[0] = 0;
	p->velocity[1] = 0;
	p->velocity[2] = 0;
#ifndef COOLING_NONE
	p->u() = 0;
	p->uDot() = 0.0;
#endif
	p->fMetals() = 0.0;
	}
    }

void DistDeletedGasSmoothParams::combSmoothCache(GravityParticle *p1,
					  ExternalSmoothParticle *p2)
{
    /*
     * Distribute u, v, and fMetals for particles returning from cache
     * so that everything is conserved nicely.  
     */
    if(!TYPETest((p1), TYPE_DELETED)) {
	double delta_m = p2->mass;
	double m_new,f1,f2;
	double fTCool; /* time to cool to zero */
	m_new = p1->mass + delta_m;
	if (m_new > 0) {
	    f1 = p1->mass /m_new;
	    f2 = delta_m  /m_new;
	    p1->mass = m_new;
	    p1->velocity = f1*p1->velocity + f2*p2->velocity;            
	    p1->fMetals() = f1*p1->fMetals() + f2*p2->fMetals;
#ifndef COOLING_NONE
	    if(p1->uDot() < 0.0) /* margin of 1% to avoid roundoff
				  * problems */
		fTCool = 1.01*p1->uPred()/p1->uDot(); 
	    p1->u() = f1*p1->u() + f2*p2->u;
	    p1->uPred() = f1*p1->uPred() + f2*p2->uPred;
	    if(p1->uDot() < 0.0)
		p1->uDot() = p1->uPred()/fTCool;
#endif
	    }
	}
    }

void DistDeletedGasSmoothParams::fcnSmooth(GravityParticle *p, int nSmooth,
					   pqSmoothNode *nnList)
{	
    GravityParticle *q;
    double fNorm,ih2,r2,rs,rstot,delta_m,m_new,f1,f2;
    double fTCool; /* time to cool to zero */
    int i;
    CkAssert(TYPETest(p, TYPE_GAS));
    ih2 = invH2(p);
    rstot = 0;        
    for (i=0;i<nSmooth;++i) {
	double fDist2 = nnList[i].fKey;
	q = nnList[i].p;
	if(TYPETest(q, TYPE_DELETED)) continue;
	CkAssert(TYPETest(q, TYPE_GAS));
	r2 = fDist2*ih2;            
        rs = KERNEL(r2);
	rstot += rs;
        }
    if(rstot <= 0.0) {
	if(p->mass == 0.0) /* the particle to be deleted has NOTHING */
	    return;
	/* we have a particle to delete and nowhere to put its mass
	 * => we will keep it around */
	unDeleteParticle(p);
	return;
	}
    CkAssert(rstot > 0.0);
    fNorm = 1./rstot;
    CkAssert(p->mass >= 0.0);
    for (i=0;i<nSmooth;++i) {
	q = nnList[i].p;
	if(TYPETest(q, TYPE_DELETED)) continue;

	double fDist2 = nnList[i].fKey;
	r2 = fDist2*ih2;            
	rs = KERNEL(r2);
	/*
	 * All these quantities are per unit mass.
	 * Exact if only one gas particle being distributed or in serial
	 * Approximate in parallel (small error).
	 */
	delta_m = rs*fNorm*p->mass;
	m_new = q->mass + delta_m;
	/* Cached copies can have zero mass: skip them */
	if (m_new == 0) continue;
	f1 = q->mass /m_new;
	f2 = delta_m  /m_new;
	q->mass = m_new;
	q->velocity = f1*q->velocity + f2*p->velocity;            
	q->fMetals() = f1*q->fMetals() + f2*p->fMetals();
#ifndef COOLING_NONE
	if(q->uDot() < 0.0) /* margin of 1% to avoid roundoff error */
	    fTCool = 1.01*q->uPred()/q->uDot(); 
	q->u() = f1*q->u()+f2*p->u();
	q->uPred() = f1*q->uPred()+f2*p->uPred();
	if(q->uDot() < 0.0) /* make sure we don't shorten cooling time */
	    q->uDot() = q->uPred()/fTCool;
#endif
        }
    }
