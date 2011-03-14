/* Classes to describe smooth functions needed for SPH */
#ifndef __SPH_H
#define __SPH_H

class DenDvDxSmoothParams : public SmoothParams
{
    double a, H; // Cosmological parameters
    int bActiveOnly;
    
    virtual void fcnSmooth(GravityParticle *p, int nSmooth,
			   pqSmoothNode *nList);
    virtual int isSmoothActive(GravityParticle *p);
    virtual void initTreeParticle(GravityParticle *p);
    virtual void postTreeParticle(GravityParticle *p) {} 
    virtual void initSmoothParticle(GravityParticle *p);
    virtual void initSmoothCache(GravityParticle *p);
    virtual void combSmoothCache(GravityParticle *p1,
				 ExternalSmoothParticle *p2);
 public:
    DenDvDxSmoothParams() {}
    DenDvDxSmoothParams(int _iType, int am, CSM csm, double dTime,
			int _bActiveOnly) {
	iType = _iType;
	activeRung = am;
	bActiveOnly = _bActiveOnly;
	if(csm->bComove) {
	    H = csmTime2Hub(csm,dTime);
	    a = csmTime2Exp(csm,dTime);
	    }
	else {
	    H = 0.0;
	    a = 1.0;
	    }
    }
    PUPable_decl(DenDvDxSmoothParams);
    DenDvDxSmoothParams(CkMigrateMessage *m) : SmoothParams(m) {}
    virtual void pup(PUP::er &p) {
        SmoothParams::pup(p);//Call base class
	p|a;
	p|H;
	p|bActiveOnly;
	}
    };

// Like the above class, but only does non active particles marked
// "Neighbor of Active" for the "fast gas" option.
// Also, it doesn't mark any particles.

class DenDvDxNeighborSmParams : public DenDvDxSmoothParams
{
    double a, H; // Cosmological parameters
    int bActiveOnly;
    
    virtual void fcnSmooth(GravityParticle *p, int nSmooth,
			   pqSmoothNode *nList);
    virtual int isSmoothActive(GravityParticle *p);
    virtual void initTreeParticle(GravityParticle *p) {}
    virtual void postTreeParticle(GravityParticle *p) {} 
    virtual void initSmoothParticle(GravityParticle *p) {}
    virtual void initSmoothCache(GravityParticle *p) {}
    virtual void combSmoothCache(GravityParticle *p1,
				 ExternalSmoothParticle *p2) {}
 public:
    DenDvDxNeighborSmParams() {}
    DenDvDxNeighborSmParams(int _iType, int am, CSM csm, double dTime)
	: DenDvDxSmoothParams(_iType, am, csm, dTime, 0) {}
    PUPable_decl(DenDvDxNeighborSmParams);
    DenDvDxNeighborSmParams(CkMigrateMessage *m) : DenDvDxSmoothParams(m) {}
    virtual void pup(PUP::er &p) {
        DenDvDxSmoothParams::pup(p);//Call base class
	}
    };

class MarkSmoothParams : public SmoothParams
{
    virtual void fcnSmooth(GravityParticle *p, int nSmooth,
			   pqSmoothNode *nList) {}
    virtual int isSmoothActive(GravityParticle *p);
    virtual void initTreeParticle(GravityParticle *p) {}
    virtual void postTreeParticle(GravityParticle *p) {} 
    virtual void initSmoothParticle(GravityParticle *p) {}
    virtual void initSmoothCache(GravityParticle *p) {}
    virtual void combSmoothCache(GravityParticle *p1,
				 ExternalSmoothParticle *p2) {}
 public:
    MarkSmoothParams() {}
    MarkSmoothParams(int _iType, int am) {
	iType = _iType;
	activeRung = am;
	}
    PUPable_decl(MarkSmoothParams);
    MarkSmoothParams(CkMigrateMessage *m) : SmoothParams(m) {}
    virtual void pup(PUP::er &p) {
        SmoothParams::pup(p);//Call base class
	}
    };

class PressureSmoothParams : public SmoothParams
{
    double a, H; // Cosmological parameters
    double alpha, beta; // SPH viscosity parameters
    
    virtual void fcnSmooth(GravityParticle *p, int nSmooth,
			   pqSmoothNode *nList);
    virtual int isSmoothActive(GravityParticle *p);
    virtual void initTreeParticle(GravityParticle *p) {}
    virtual void postTreeParticle(GravityParticle *p) {} 
    virtual void initSmoothParticle(GravityParticle *p);
    virtual void initSmoothCache(GravityParticle *p);
    virtual void combSmoothCache(GravityParticle *p1,
				 ExternalSmoothParticle *p2);
 public:
    PressureSmoothParams() {}
    PressureSmoothParams(int _iType, int am, CSM csm, double dTime,
			 double _alpha, double _beta) {
	iType = _iType;
	activeRung = am;
	if(csm->bComove) {
	    H = csmTime2Hub(csm,dTime);
	    a = csmTime2Exp(csm,dTime);
	    }
	else {
	    H = 0.0;
	    a = 1.0;
	    }
	alpha = _alpha;
	beta = _beta;
    }
    PUPable_decl(PressureSmoothParams);
    PressureSmoothParams(CkMigrateMessage *m) : SmoothParams(m) {}
    virtual void pup(PUP::er &p) {
        SmoothParams::pup(p);//Call base class
	p|a;
	p|H;
	p|alpha;
	p|beta;
	}
    };

/*
 * SmoothParams class for distributing deleted gas to neighboring
 * particles.
 */
class DistDeletedGasSmoothParams : public SmoothParams
{
    virtual void fcnSmooth(GravityParticle *p, int nSmooth,
			   pqSmoothNode *nList);
    virtual int isSmoothActive(GravityParticle *p);
    virtual void initSmoothParticle(GravityParticle *p) {};
    virtual void initTreeParticle(GravityParticle *p) {}
    virtual void postTreeParticle(GravityParticle *p) {}
    virtual void initSmoothCache(GravityParticle *p);
    virtual void combSmoothCache(GravityParticle *p1,
				 ExternalSmoothParticle *p2);
 public:
    DistDeletedGasSmoothParams() {}
    DistDeletedGasSmoothParams(int _iType, int am) {
	iType = _iType;
	activeRung = am;
	}
    PUPable_decl(DistDeletedGasSmoothParams);
    DistDeletedGasSmoothParams(CkMigrateMessage *m) : SmoothParams(m) {}
    virtual void pup(PUP::er &p) {
        SmoothParams::pup(p);//Call base class
	}
    };

/*
 * SmoothParams class for distributing stellar feedback (energy, mass + metals) 
 * to neighboring particles.
 */
class DistStellarFeedbackSmoothParams : public SmoothParams
{
    double dTime, H, a, gamma;
    Fdbk *fb;
    virtual void fcnSmooth(GravityParticle *p, int nSmooth,
			   pqSmoothNode *nList);
    virtual int isSmoothActive(GravityParticle *p);
    virtual void initSmoothParticle(GravityParticle *p) {}
    virtual void initTreeParticle(GravityParticle *p);
    virtual void postTreeParticle(GravityParticle *p);
    virtual void initSmoothCache(GravityParticle *p);
    virtual void combSmoothCache(GravityParticle *p1,
				 ExternalSmoothParticle *p2);
    void DistFBMME(GravityParticle *p, int nSmooth, pqSmoothNode *nList);
 public:
    DistStellarFeedbackSmoothParams() {}
    DistStellarFeedbackSmoothParams(int _iType, int am, CSM csm, double _dTime,
				    double _gamma, Fdbk *feedback) : 
    fb (feedback) {
	iType = _iType;
	activeRung = am;
	gamma = _gamma;
	dTime = _dTime;
	if(csm->bComove) {
	    H = csmTime2Hub(csm,dTime);
	    a = csmTime2Exp(csm,dTime);
	    }
	else {
	    H = 0.0;
	    a = 1.0;
	    }
	}
    /*    ~DistStellarFeedbackSmoothParams() {
	delete fb;
	}*/
    PUPable_decl(DistStellarFeedbackSmoothParams);
    DistStellarFeedbackSmoothParams(CkMigrateMessage *m) : SmoothParams(m) {}
    virtual void pup(PUP::er &p) {
        SmoothParams::pup(p);//Call base class
	p|a;
	p|H;
	p|gamma;
	p|fb;
	p|dTime;
	}
    };

#endif
