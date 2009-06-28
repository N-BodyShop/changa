#include "ParallelGravity.h"
#include "smooth.h"
#include "Sph.h"


void
Main::initSph() 
{
    // Assume Tree is built
    ckout << "Calculating densities ...";
    DensitySmoothParams pDen(TYPE_GAS, 0);
    double startTime = CkWallTimer();
    treeProxy.startIterationSmooth(&pDen, CkCallbackResumeThread());
    ckout << " took " << (CkWallTimer() - startTime) << " seconds."
	  << endl;
    
    if(param.bDoGas) {
	if(verbosity) CkPrintf("Initializing SPH forces\n");
	doSph(0);
	}
    }

void
Main::doSph(int activeRung) 
{
    ckout << "Calculating densities/divv ...";
    DenDvDxSmoothParams pDen(TYPE_GAS, activeRung, param.csm, dTime);
    double startTime = CkWallTimer();
    treeProxy.startIterationSmooth(&pDen, CkCallbackResumeThread());
    ckout << " took " << (CkWallTimer() - startTime) << " seconds."
	  << endl;

    treeProxy.sphViscosityLimiter(param.iViscosityLimiter, activeRung,
			CkCallbackResumeThread());
    treeProxy.getAdiabaticGasPressure(param.dConstGamma, param.dConstGamma-1,
				      CkCallbackResumeThread());

    ckout << "Calculating pressure gradients ...";
    PressureSmoothParams pPressure(TYPE_GAS, activeRung, param.csm, dTime,
				   param.dConstAlpha, param.dConstBeta);
    startTime = CkWallTimer();
    treeProxy.startIterationReSmooth(&pPressure,
				     CkCallbackResumeThread());
    ckout << " took " << (CkWallTimer() - startTime) << " seconds."
	  << endl;
    }

void DenDvDxSmoothParams::initSmoothParticle(GravityParticle *p)
{
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
		double fDist2 = nnList[i].fKey*nnList[i].fKey;
		r2 = fDist2*ih2;
		q = nnList[i].p;
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
	    double fDist2 = nnList[i].fKey*nnList[i].fKey;
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
