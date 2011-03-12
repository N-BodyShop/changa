/*
 * Routines to implement SPH.
 * Main author: James Wadsley, as first implemented in GASOLINE.
 * See Wadsley, J.~W., Stadel, J., Quinn, T.\ 2004.\ Gasoline: a flexible,
 * parallel implementation of TreeSPH.\ New Astronomy 9, 137-158.
 */

#include "ParallelGravity.h"
#include "DataManager.h"
#include "smooth.h"
#include "Sph.h"
#include "physconst.h"

#ifndef MAXPATHLEN
#define MAXPATHLEN PATH_MAX
#endif

void
Main::initSph() 
{
    if(param.bDoGas) {
	ckout << "Calculating densities/divv ...";
	// The following smooths all GAS, and also marks neighbors of
	// actives, and those who have actives as neighbors.
	DenDvDxSmoothParams pDen(TYPE_GAS, 0, param.csm, dTime, 0);
	double startTime = CkWallTimer();
	treeProxy.startIterationSmooth(&pDen, CkCallbackResumeThread());
	iPhase++;
	ckout << " took " << (CkWallTimer() - startTime) << " seconds."
	      << endl;
	if(verbosity)
	    memoryStatsCache();
	double dTuFac = param.dGasConst/(param.dConstGamma-1)
	    /param.dMeanMolWeight;
	double z = 1.0/csmTime2Exp(param.csm, dTime) - 1.0;
	// Update cooling on the datamanager
	dMProxy.CoolingSetTime(z, dTime, CkCallbackResumeThread());
	if(param.bGasCooling)
	    treeProxy.InitEnergy(dTuFac, z, dTime, CkCallbackResumeThread());
	if(verbosity) CkPrintf("Initializing SPH forces\n");
	nActiveSPH = nTotalSPH;
	doSph(0, 0);
	double duDelta[MAXRUNG+1];
	for(int iRung = 0; iRung <= MAXRUNG; iRung++)
	    duDelta[iRung] = 0.5e-7*param.dDelta;
	if(param.bGasCooling)
	    dMProxy.CoolingSetTime(z, dTime, CkCallbackResumeThread());
	treeProxy.updateuDot(0, duDelta, dTime, z, param.bGasCooling, 0,
			     CkCallbackResumeThread());
	}
    }

/*
 * Initialize cooling constants and integration data structures.
 * XXX if this is to work on an SMP configuration, we need to change
 * this so that the integration data structures are per treepiece, not
 * per node.
 */
void Main::initCooling()
{
#ifndef COOLING_NONE
    dMProxy.initCooling(param.dGmPerCcUnit, param.dComovingGmPerCcUnit,
		    param.dErgPerGmUnit, param.dSecUnit, param.dKpcUnit,
		    param.CoolParam, CkCallbackResumeThread());
    
    /* Read in tables from files as necessary */
    int cntTable = 0;
    int nTableRows;
    int nTableColumns;
    char TableFileSuffix[20];
    
    for (;;) {
	CoolTableReadInfo(&param.CoolParam, cntTable, &nTableColumns,
			  TableFileSuffix);
	if (!nTableColumns) break;
    
	cntTable++;
	nTableRows = ReadASCII(TableFileSuffix, nTableColumns, NULL);
	if (nTableRows) {
	    CkAssert(sizeof(double)*nTableRows*nTableColumns <= CL_NMAXBYTETABLE );
	    double *dTableData = (double *)malloc(sizeof(double)*nTableRows*nTableColumns);
	    CkAssert( dTableData != NULL );
	    nTableRows = ReadASCII(TableFileSuffix, nTableColumns, dTableData);
      
	    dMProxy.dmCoolTableRead(dTableData,nTableRows*nTableColumns,
				  CkCallbackResumeThread());
	    free(dTableData);
	    }
	}
    treeProxy.initCoolingData(CkCallbackResumeThread());
#endif
    }

/**
 * Initialized Cooling Read-only data on the DataManager, which
 * doesn't migrate.
 */
void
DataManager::initCooling(double dGmPerCcUnit, double dComovingGmPerCcUnit,
		       double dErgPerGmUnit, double dSecUnit, double dKpcUnit,
		       COOLPARAM inParam, const CkCallback& cb)
{
#ifndef COOLING_NONE
    clInitConstants(Cool, dGmPerCcUnit, dComovingGmPerCcUnit, dErgPerGmUnit,
		    dSecUnit, dKpcUnit, inParam);
    
    CoolInitRatesTable(Cool,inParam);
#endif
    contribute(0, 0, CkReduction::concat, cb);
    }

/**
 * Per thread initialization
 */
void
TreePiece::initCoolingData(const CkCallback& cb)
{
#ifndef COOLING_NONE
    dm = (DataManager*)CkLocalNodeBranch(dataManagerID);
    CoolData = CoolDerivsInit(dm->Cool);
#endif
    contribute(0, 0, CkReduction::concat, cb);
    }

void
DataManager::dmCoolTableRead(double *dTableData, int nData, const CkCallback& cb)
{
#ifndef COOLING_NONE
    CoolTableRead(Cool, nData*sizeof(double), (void *) dTableData);
#endif
    contribute(0, 0, CkReduction::concat, cb);
    }

/*
 * function from PKDGRAV to read an ASCII table
 */
/* Note if dDataOut is NULL it just counts the number of valid input lines */
int Main::ReadASCII(char *extension, int nDataPerLine, double *dDataOut)
{
	FILE *fp;
	int i,ret;
	char achIn[160];
	double *dData;

	if (dDataOut == NULL) 
	    dData = (double *)malloc(sizeof(double)*nDataPerLine);
	else
	    dData = dDataOut;
	
	CkAssert(nDataPerLine > 0 && nDataPerLine <= 10);
	char achFile[MAXPATHLEN];
	sprintf(achFile, "%s.%s", param.achOutName, extension);
	fp = fopen(achFile,"r");
	if (!fp) {
	    if (verbosity)
		CkPrintf("WARNING: Could not open .%s input file:%s\n",
			 extension,achFile);
	    return 0;
	    }

	i = 0;
	while (1) {
	    if (!fgets(achIn,160,fp)) goto Done;
	    switch (nDataPerLine) {
	    case 1:
		ret = sscanf(achIn,"%lf",dData); 
		break;
	    case 2:
		ret = sscanf(achIn,"%lf %lf",dData,dData+1); 
		break;
	    case 3:
		ret = sscanf(achIn,"%lf %lf %lf",dData,dData+1,dData+2); 
		break;
	    case 4:
		ret = sscanf(achIn,"%lf %lf %lf %lf",dData,dData+1,dData+2,dData+3); 
		break;
	    case 5:
		ret = sscanf(achIn,"%lf %lf %lf %lf %lf",dData,dData+1,dData+2,dData+3,dData+4); 
		break;
	    case 6:
		ret = sscanf(achIn,"%lf %lf %lf %lf %lf %lf",dData,dData+1,dData+2,dData+3,dData+4,dData+5); 
		break;
	    case 7:
		ret = sscanf(achIn,"%lf %lf %lf %lf %lf %lf %lf",
			     dData,dData+1,dData+2,dData+3,dData+4,dData+5,dData+6); 
		break;
	    case 8:
		ret = sscanf(achIn,"%lf %lf %lf %lf %lf %lf %lf %lf",
			     dData,dData+1,dData+2,dData+3,dData+4,dData+5,dData+6,dData+7); 
		break;
	    case 9:
		ret = sscanf(achIn,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
			     dData,dData+1,dData+2,dData+3,dData+4,dData+5,dData+6,dData+7,dData+8); 
		break;
	    case 10:
		ret = sscanf(achIn,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
			     dData,dData+1,dData+2,dData+3,dData+4,dData+5,dData+6,dData+7,dData+8,dData+9); 
		break;
	    default:
		ret = EOF;
		CkAssert(0);
		}
	    if (ret != nDataPerLine) goto Done;
	    ++i;
	    if (dDataOut != NULL) dData += nDataPerLine;
	    }
 Done:
	fclose(fp);
	if (dDataOut != NULL && verbosity)
	    printf("Read %i lines from %s\n",i,achFile);
	if (dDataOut == NULL) free(dData);
	return i;
	}

/*
 * Update the cooling functions to the current time.
 * This is on the DataManager to avoid duplication of effort.
 */
void
DataManager::CoolingSetTime(double z, // redshift
			    double dTime, // Time
			    const CkCallback& cb)
{
#ifndef COOLING_NONE
    CoolSetTime( Cool, dTime, z  );
#endif

    contribute(0, 0, CkReduction::concat, cb);
    }

/**
 *  Perform the SPH force calculation.
 *  @param activeRung Timestep rung (and above) on which to perform
 *  SPH
 *  @param bNeedDensity Does the density calculation need to be done?
 *  Defaults to 1
 */
void
Main::doSph(int activeRung, int bNeedDensity) 
{
  if(bNeedDensity) {
    if (param.bFastGas && nActiveSPH < nTotalSPH*param.dFracFastGas) {
	ckout << "Calculating densities/divv on Actives ...";
	// This also marks neighbors of actives
	DenDvDxSmoothParams pDen(TYPE_GAS, activeRung, param.csm, dTime, 1);
	double startTime = CkWallTimer();
	treeProxy.startIterationSmooth(&pDen, CkCallbackResumeThread());
	iPhase++;
	ckout << " took " << (CkWallTimer() - startTime) << " seconds."
	      << endl;

	ckout << "Marking Neighbors ...";
	// This marks particles with actives as neighbors
	MarkSmoothParams pMark(TYPE_GAS, activeRung);
	startTime = CkWallTimer();
	treeProxy.startIterationMarkSmooth(&pMark, CkCallbackResumeThread());
	iPhase++;
	ckout << " took " << (CkWallTimer() - startTime) << " seconds."
	      << endl;
	
	ckout << "Density of Neighbors ...";
	// This does neighbors (but not actives),  It also does no
	// additional marking
	DenDvDxNeighborSmParams pDenN(TYPE_GAS, activeRung, param.csm, dTime);
	startTime = CkWallTimer();
	treeProxy.startIterationSmooth(&pDenN, CkCallbackResumeThread());
	iPhase++;
	ckout << " took " << (CkWallTimer() - startTime) << " seconds."
	      << endl;
	}
    else {
	ckout << "Calculating densities/divv ...";
	// The following smooths all GAS, and also marks neighbors of
	// actives, and those who have actives as neighbors.
	DenDvDxSmoothParams pDen(TYPE_GAS, activeRung, param.csm, dTime, 0);
	double startTime = CkWallTimer();
	treeProxy.startIterationSmooth(&pDen, CkCallbackResumeThread());
	iPhase++;
	ckout << " took " << (CkWallTimer() - startTime) << " seconds."
	      << endl;

	if(verbosity)
	    memoryStatsCache();
	}
      }
    treeProxy.sphViscosityLimiter(param.iViscosityLimiter, activeRung,
			CkCallbackResumeThread());
    if(param.bGasCooling)
	treeProxy.getCoolingGasPressure(param.dConstGamma,
					param.dConstGamma-1,
					CkCallbackResumeThread());
    else
	treeProxy.getAdiabaticGasPressure(param.dConstGamma,
					  param.dConstGamma-1,
					  CkCallbackResumeThread());

    ckout << "Calculating pressure gradients ...";
    PressureSmoothParams pPressure(TYPE_GAS, activeRung, param.csm, dTime,
				   param.dConstAlpha, param.dConstBeta);
    double startTime = CkWallTimer();
    treeProxy.startIterationReSmooth(&pPressure, CkCallbackResumeThread());
    iPhase++;
    ckout << " took " << (CkWallTimer() - startTime) << " seconds."
	  << endl;
    
    treeProxy.ballMax(activeRung, 1.0+param.ddHonHLimit,
		      CkCallbackResumeThread());
    }

/*
 * Initialize energy and ionization state for cooling particles
 */
void TreePiece::InitEnergy(double dTuFac, // T to internal energy
			   double z,	  // redshift
			   double dTime,
			   const CkCallback& cb)
{
#ifndef COOLING_NONE
    COOL *cl;

    dm = (DataManager*)CkLocalNodeBranch(dataManagerID);
    cl = dm->Cool;
#endif

    for(unsigned int i = 1; i <= myNumParticles; ++i) {
	GravityParticle *p = &myParticles[i];
	if (TYPETest(p, TYPE_GAS) && p->rung >= activeRung) {
	    double T,E;
#ifndef COOLING_NONE
	    T = p->u() / dTuFac;
	    CoolInitEnergyAndParticleData(cl, &p->CoolParticle(), &E,
					  p->fDensity, T, p->fMetals() );
	    p->u() = E;
#endif
	    p->uPred() = p->u();
	    }
	}
    contribute(0, 0, CkReduction::concat, cb);
    }

/**
 * Update the cooling rate (uDot)
 */
void TreePiece::updateuDot(int activeRung,
			   double duDelta[MAXRUNG+1], // timesteps
			   double dTime, // current time
			   double z, // current redshift
			   int bCool, // select equation of state
			   int bUpdateState, // update ionization fractions
			   const CkCallback& cb)
{
    double dt; // time in seconds
    
#ifndef COOLING_NONE
    for(unsigned int i = 1; i <= myNumParticles; ++i) {
	GravityParticle *p = &myParticles[i];
	if (TYPETest(p, TYPE_GAS) && p->rung >= activeRung) {
	    dt = CoolCodeTimeToSeconds(dm->Cool, duDelta[p->rung] );
	    double ExternalHeating = p->PdV(); // Will change with star formation
	    if ( bCool ) {
		COOLPARTICLE cp = p->CoolParticle();
		double E = p->u();
		double r[3];  // For conversion to C
		p->position.array_form(r);
#ifdef COOLING_BOLEY
		cp.mrho = pow(p->mass/p->fDensity, 1./3.);
#endif

		CoolIntegrateEnergyCode(dm->Cool, CoolData, &cp, &E,
					ExternalHeating, p->fDensity,
					p->fMetals(), r, dt);
		CkAssert(E > 0.0);
		p->uDot() = (E - p->u())/duDelta[p->rung]; // linear interpolation over interval
		if (bUpdateState) p->CoolParticle() = cp;
		}
	    else { 
		p->uDot() = ExternalHeating;
		}
	    }
	}
#endif
    contribute(0, 0, CkReduction::concat, cb);
    }

/* Set a maximum ball for inverse Nearest Neighbor searching */
void TreePiece::ballMax(int activeRung, double dhFac, const CkCallback& cb)
{
    for(unsigned int i = 1; i <= myNumParticles; ++i) {
	if (TYPETest(&myParticles[i], TYPE_GAS)
	    && myParticles[i].rung >= activeRung) {
	    myParticles[i].fBallMax() = myParticles[i].fBall*dhFac;
	    }
	}
    contribute(0, 0, CkReduction::concat, cb);
    }
    
int DenDvDxSmoothParams::isSmoothActive(GravityParticle *p) 
{
    if(bActiveOnly && p->rung < activeRung)
	return 0;		// not active

    return (TYPETest(p, iType));
    }

// Non-active neighbors of Actives
int DenDvDxNeighborSmParams::isSmoothActive(GravityParticle *p) 
{
    if(p->rung < activeRung && TYPETest(p, iType)
       && TYPETest(p, TYPE_NbrOfACTIVE))
	return 1;

    return 0;
    }

// Only do actives

int MarkSmoothParams::isSmoothActive(GravityParticle *p) 
{
    if(p->rung < activeRung)
	return 0;		// not active

    return (TYPETest(p, iType));
    }

void DenDvDxSmoothParams::initSmoothParticle(GravityParticle *p)
{
    TYPEReset(p, TYPE_NbrOfACTIVE);
    }

void DenDvDxSmoothParams::initTreeParticle(GravityParticle *p)
{
    TYPEReset(p, TYPE_NbrOfACTIVE);
    }

void DenDvDxSmoothParams::initSmoothCache(GravityParticle *p)
{
	}

void DenDvDxSmoothParams::combSmoothCache(GravityParticle *p1,
					  ExternalSmoothParticle *p2)
{
	p1->iType |= p2->iType;
	}

/* Gather only version */
void DenDvDxSmoothParams::fcnSmooth(GravityParticle *p, int nSmooth,
				    pqSmoothNode *nnList)
{
	double ih2,r2,rs,rs1,fDensity,fNorm,fNorm1,vFac;
	double dvxdx, dvxdy, dvxdz, dvydx, dvydy, dvydz, dvzdx, dvzdy, dvzdz;
	double dvx,dvy,dvz,dx,dy,dz,trace;
	GravityParticle *q;
	int i;
	unsigned int qiActive;

	ih2 = invH2(p);
	vFac = 1./(a*a); /* converts v to xdot */
	fNorm = M_1_PI*ih2*sqrt(ih2);
	fNorm1 = fNorm*ih2;	
	fDensity = 0.0;
	dvxdx = 0; dvxdy = 0; dvxdz= 0;
	dvydx = 0; dvydy = 0; dvydz= 0;
	dvzdx = 0; dvzdy = 0; dvzdz= 0;

	qiActive = 0;
	for (i=0;i<nSmooth;++i) {
		double fDist2 = nnList[i].fKey;
		r2 = fDist2*ih2;
		q = nnList[i].p;
		if(q == NULL)
		    CkAbort("NULL neighbor in DenDvDxSmooth");
		if (p->rung >= activeRung)
		    TYPESet(q,TYPE_NbrOfACTIVE); /* important for SPH */
		if(q->rung >= activeRung)
		    qiActive = 1;
		rs = KERNEL(r2);
		fDensity += rs*q->mass;
		rs1 = DKERNEL(r2);
		rs1 *= q->mass;
		dx = nnList[i].dx.x;
		dy = nnList[i].dx.y;
		dz = nnList[i].dx.z;
		dvx = (-p->vPred().x + q->vPred().x)*vFac - dx*H; /* NB: dx = px - qx */
		dvy = (-p->vPred().y + q->vPred().y)*vFac - dy*H;
		dvz = (-p->vPred().z + q->vPred().z)*vFac - dz*H;
		dvxdx += dvx*dx*rs1;
		dvxdy += dvx*dy*rs1;
		dvxdz += dvx*dz*rs1;
		dvydx += dvy*dx*rs1;
		dvydy += dvy*dy*rs1;
		dvydz += dvy*dz*rs1;
		dvzdx += dvz*dx*rs1;
		dvzdy += dvz*dy*rs1;
		dvzdz += dvz*dz*rs1;
		}
	if (qiActive)
	    TYPESet(p,TYPE_NbrOfACTIVE);
		
	p->fDensity = fNorm*fDensity; 
	fNorm1 /= p->fDensity;
	trace = dvxdx+dvydy+dvzdz;
	p->divv() =  fNorm1*trace; /* physical */
	p->curlv().x = fNorm1*(dvzdy - dvydz); 
	p->curlv().y = fNorm1*(dvxdz - dvzdx);
	p->curlv().z = fNorm1*(dvydx - dvxdy);
	}

/* As above, but no marking */
void DenDvDxNeighborSmParams::fcnSmooth(GravityParticle *p, int nSmooth,
				    pqSmoothNode *nnList)
{
	double ih2,r2,rs,rs1,fDensity,fNorm,fNorm1,vFac;
	double dvxdx, dvxdy, dvxdz, dvydx, dvydy, dvydz, dvzdx, dvzdy, dvzdz;
	double dvx,dvy,dvz,dx,dy,dz,trace;
	GravityParticle *q;
	int i;

	ih2 = invH2(p);
	vFac = 1./(a*a); /* converts v to xdot */
	fNorm = M_1_PI*ih2*sqrt(ih2);
	fNorm1 = fNorm*ih2;	
	fDensity = 0.0;
	dvxdx = 0; dvxdy = 0; dvxdz= 0;
	dvydx = 0; dvydy = 0; dvydz= 0;
	dvzdx = 0; dvzdy = 0; dvzdz= 0;

	for (i=0;i<nSmooth;++i) {
		double fDist2 = nnList[i].fKey;
		r2 = fDist2*ih2;
		q = nnList[i].p;
		rs = KERNEL(r2);
		fDensity += rs*q->mass;
		rs1 = DKERNEL(r2);
		rs1 *= q->mass;
		dx = nnList[i].dx.x;
		dy = nnList[i].dx.y;
		dz = nnList[i].dx.z;
		dvx = (-p->vPred().x + q->vPred().x)*vFac - dx*H; /* NB: dx = px - qx */
		dvy = (-p->vPred().y + q->vPred().y)*vFac - dy*H;
		dvz = (-p->vPred().z + q->vPred().z)*vFac - dz*H;
		dvxdx += dvx*dx*rs1;
		dvxdy += dvx*dy*rs1;
		dvxdz += dvx*dz*rs1;
		dvydx += dvy*dx*rs1;
		dvydy += dvy*dy*rs1;
		dvydz += dvy*dz*rs1;
		dvzdx += dvz*dx*rs1;
		dvzdy += dvz*dy*rs1;
		dvzdz += dvz*dz*rs1;
		}
		
	p->fDensity = fNorm*fDensity; 
	fNorm1 /= p->fDensity;
	trace = dvxdx+dvydy+dvzdz;
	p->divv() =  fNorm1*trace; /* physical */
	p->curlv().x = fNorm1*(dvzdy - dvydz); 
	p->curlv().y = fNorm1*(dvxdz - dvzdx);
	p->curlv().z = fNorm1*(dvydx - dvxdy);
	}

void 
TreePiece::sphViscosityLimiter(int bOn, int activeRung, const CkCallback& cb)
{
    int i;
    GravityParticle *p;    

    if (bOn) {
        for(i=1; i<= myNumParticles; ++i) {
	    p = &myParticles[i];
	    /* Only set values for particles with fresh curlv, divv
	       from smooth */
	    if(TYPETest(p, TYPE_GAS) && p->rung >= activeRung) {
		if (p->divv() != 0.0) {         	 
		    p->BalsaraSwitch() = fabs(p->divv())/
			(fabs(p->divv()) + sqrt(p->curlv().lengthSquared()));
		    }
		else { 
		    p->BalsaraSwitch() = 0.0;
		    }
		}
	    }
        }
    else {
        for(i=1; i<= myNumParticles; ++i) {
	    p = &myParticles[i];
	    if(TYPETest(p, TYPE_GAS)) {
		p->BalsaraSwitch() = 1.0;
		}
	    }
        }
    contribute(0, 0, CkReduction::concat, cb);
    }

/* Note: Uses uPred */
void TreePiece::getAdiabaticGasPressure(double gamma, double gammam1,
					const CkCallback &cb)
{
    GravityParticle *p;
    double PoverRho;
    int i;

    for(i=1; i<= myNumParticles; ++i) {
	p = &myParticles[i];
	if (TYPETest(p, TYPE_GAS)) {
	    PoverRho = gammam1*p->uPred();
	    p->PoverRho2() = PoverRho/p->fDensity;
	    p->c() = sqrt(gamma*PoverRho);
	    }
	}
    contribute(0, 0, CkReduction::concat, cb);
    }

/* Note: Uses uPred */
void TreePiece::getCoolingGasPressure(double gamma, double gammam1,
					const CkCallback &cb)
{
    GravityParticle *p;
    double PoverRho;
    int i;
#ifndef COOLING_NONE
    COOL *cl = dm->Cool;

    for(i=1; i<= myNumParticles; ++i) {
	p = &myParticles[i];
	if (TYPETest(p, TYPE_GAS)) {
	    CoolCodePressureOnDensitySoundSpeed(cl, &p->CoolParticle(),
						p->uPred(), p->fDensity(),
						gamma, gammam1, &PoverRho,
						&(p->c()) );
	    p->PoverRho2() = PoverRho/p->fDensity;
	    }
	}
#endif
    contribute(0, 0, CkReduction::concat, cb);
    }

int PressureSmoothParams::isSmoothActive(GravityParticle *p) 
{
    return (TYPETest(p, TYPE_NbrOfACTIVE));
    }

/* Original Particle */
void PressureSmoothParams::initSmoothParticle(GravityParticle *p)
{
	if (p->rung >= activeRung) {
	    p->mumax() = 0.0;
	    p->PdV() = 0.0;
	    }
	}

/* Cached copies of particle */
void PressureSmoothParams::initSmoothCache(GravityParticle *p)
{
	if (p->rung >= activeRung) {
	    p->mumax() = 0.0;
	    p->PdV() = 0.0;
	    p->treeAcceleration = 0.0;
	    }
	}

void PressureSmoothParams::combSmoothCache(GravityParticle *p1,
					  ExternalSmoothParticle *p2)
{
	if (p1->rung >= activeRung) {
	    p1->PdV() += p2->PdV;
	    if (p2->mumax > p1->mumax())
		p1->mumax() = p2->mumax;
	    p1->treeAcceleration += p2->treeAcceleration;
	    }
	}

void PressureSmoothParams::fcnSmooth(GravityParticle *p, int nSmooth,
				    pqSmoothNode *nnList)
{
	GravityParticle *q;
	double ih2,r2,rs1,rq,rp;
	double dx,dy,dz,dvx,dvy,dvz,dvdotdr;
	double pPoverRho2,pPoverRho2f,pMass;
	double qPoverRho2,qPoverRho2f;
	double ph,pc,pDensity,visc,hav,absmu,Accp,Accq;
	double fNorm,fNorm1,aFac,vFac;
	int i;

	pc = p->c();
	pDensity = p->fDensity;
	pMass = p->mass;
	pPoverRho2 = p->PoverRho2();
	pPoverRho2f = pPoverRho2;
	ph = sqrt(0.25*p->fBall*p->fBall);
	ih2 = invH2(p);
	fNorm = 0.5*M_1_PI*ih2/ph;
	fNorm1 = fNorm*ih2;	/* converts to physical u */
	aFac = a;        /* comoving acceleration factor */
	vFac = 1./(a*a); /* converts v to xdot */

	for (i=0;i<nSmooth;++i) {
	    q = nnList[i].p;
	    if ((p->rung < activeRung) && (q->rung < activeRung)) continue;
	    double fDist2 = nnList[i].fKey;
	    r2 = fDist2*ih2;
	    rs1 = DKERNEL(r2);
	    rs1 *= fNorm1;
	    rp = rs1 * pMass;
	    rq = rs1 * q->mass;

	    dx = nnList[i].dx.x;
	    dy = nnList[i].dx.y;
	    dz = nnList[i].dx.z;
	    dvx = p->vPred()[0] - q->vPred()[0];
	    dvy = p->vPred()[1] - q->vPred()[1];
	    dvz = p->vPred()[2] - q->vPred()[2];
	    dvdotdr = vFac*(dvx*dx + dvy*dy + dvz*dz) + fDist2*H;
	    
	    qPoverRho2 = q->PoverRho2();
	    qPoverRho2f = qPoverRho2;

#define PRES_PDV(a,b) (a)
#define PRES_ACC(a,b) (a+b)
#define SWITCHCOMBINE(a,b) (0.5*(a->BalsaraSwitch()+b->BalsaraSwitch()))

	    // Macro to simplify the active/inactive logic
#define SphPressureTermsSymACTIVECODE() \
	    if (dvdotdr>0.0) { \
		PACTIVE( p->PdV() += rq*PRES_PDV(pPoverRho2,qPoverRho2)*dvdotdr; ); \
		QACTIVE( q->PdV() += rp*PRES_PDV(qPoverRho2,pPoverRho2)*dvdotdr; ); \
		PACTIVE( Accp = (PRES_ACC(pPoverRho2f,qPoverRho2f)); ); \
		QACTIVE( Accq = (PRES_ACC(qPoverRho2f,pPoverRho2f)); ); \
		} \
	    else {  \
		hav=0.5*(ph+sqrt(0.25*q->fBall*q->fBall));  /* h mean - using just hp probably ok */  \
		absmu = -hav*dvdotdr*a  \
		    /(fDist2+0.01*hav*hav); /* mu multiply by a to be consistent with physical c */ \
		if (absmu>p->mumax()) p->mumax()=absmu; /* mu terms for gas time step */ \
		if (absmu>q->mumax()) q->mumax()=absmu; \
		/* viscosity term */ \
		visc = SWITCHCOMBINE(p,q)* \
		    (alpha*(pc + q->c()) + beta*2*absmu)	\
		    *absmu/(pDensity + q->fDensity); \
		PACTIVE( p->PdV() += rq*(PRES_PDV(pPoverRho2,q->PoverRho2()) + 0.5*visc)*dvdotdr; ); \
		QACTIVE( q->PdV() += rp*(PRES_PDV(q->PoverRho2(),pPoverRho2) + 0.5*visc)*dvdotdr; ); \
		PACTIVE( Accp = (PRES_ACC(pPoverRho2f,qPoverRho2f) + visc); ); \
		QACTIVE( Accq = (PRES_ACC(qPoverRho2f,pPoverRho2f) + visc); ); \
		} \
	    PACTIVE( Accp *= rq*aFac; );/* aFac - convert to comoving acceleration */ \
	    QACTIVE( Accq *= rp*aFac; ); \
	    PACTIVE( p->treeAcceleration.x -= Accp * dx; ); \
	    PACTIVE( p->treeAcceleration.y -= Accp * dy; ); \
	    PACTIVE( p->treeAcceleration.z -= Accp * dz; ); \
	    QACTIVE( q->treeAcceleration.x += Accq * dx; ); \
	    QACTIVE( q->treeAcceleration.y += Accq * dy; ); \
	    QACTIVE( q->treeAcceleration.z += Accq * dz; );


	    if (p->rung >= activeRung) {
		if (q->rung >= activeRung) {
#define PACTIVE(xxx) xxx
#define QACTIVE(xxx) xxx
		    SphPressureTermsSymACTIVECODE();    
		    }
		else {
#undef QACTIVE
#define QACTIVE(xxx) 
		    SphPressureTermsSymACTIVECODE();    
		    }
		}
	    else if (q->rung >= activeRung) {
#undef PACTIVE
#define PACTIVE(xxx) 
#undef QACTIVE
#define QACTIVE(xxx) xxx
		SphPressureTermsSymACTIVECODE();    
		}
	    }
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
    }

static inline double max(double a, double b) {
    return (a > b ? a : b);
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
    p1->fTimeForm() = max( p1->fTimeForm(),
			   p2->fTimeForm ); /* propagate FB time JMB 2/24/10 */
    }

void DistStellarFeedbackSmoothParams::DistFBMME(GravityParticle *p,int nSmooth, pqSmoothNode *nList)
{
    GravityParticle *q;
    double fNorm,ih2,r2,rs,rstot,fNorm_u,fNorm_Pres,fAveDens;
    int i,counter,imind;
    
    if ( p->fMSN() == 0.0 ){return;} /* Is there any feedback mass? */
    CkAssert(TYPETest(p, TYPE_STAR));
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
	fNorm_u += q->mass*rs;
        CkAssert(TYPETest(q, TYPE_GAS));
	rs *= fNorm;
	fAveDens += q->mass*rs;
	fNorm_Pres += q->mass*q->uPred()*rs;
	}
    fNorm_Pres *= (gamma-1.0);
    
    assert(fNorm_u != 0.0);
    fNorm_u = 1./fNorm_u;
    counter=0;
    for (i=0;i<nSmooth;++i) {
	double weight;
	q = nList[i].p;
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
	if (p->fNSN() == 0.0) q->fESNrate() += weight*p->fESNrate();
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
    if (fb.iNSNIIQuantum > 0) {
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
    if (p->fNSN() < fb.iNSNIIQuantum)
	fShutoffTime= 3e6 * SECONDSPERYEAR / fb.dSecUnit;
    
    fmind = p->fBall*p->fBall;
    imind = 0;
    if ( p->fESNrate() > 0.0 ) {
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
	    if ( fDist2 < f2h2 || !fb.bSmallSNSmooth) {
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
		    q->fTimeForm() = dTime; /* store SN FB time here JMB 2/24/10 */
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
		q->fESNrate() += weight*p->fESNrate();
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
		q->fTimeForm() = dTime;  /* store SN FB time here JMB 2/24/10 */
		counter++;
		}
	    /*	update mass after everything else so that distribution
		is based entirely upon initial mass of gas particle */
	    } 
	}
    /*if(counter>0) printf("%i ",counter);
      if (p->fNSN >0) printf("%i E51:  %g  Dens:  %g  P:  %g  R:  %g shutoff time: %g   StarAge: %g  \n",counter,p->fNSN,fAveDens,fNorm_Pres,fBlastRadius,fShutoffTime,dAge);
      /*if(p->fNSN!= 0.0)printf("E51:  %g  Dens:  %g  P:  %g  R:  %g shutoff time: %g  \n",p->fNSN,fAveDens,fNorm_Pres,fBlastRadius,fShutoffTime);*/
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
        }
    
    }

