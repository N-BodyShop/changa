#ifndef __SPH_H
#define __SPH_H

class DenDvDxSmoothParams : public SmoothParams
{
    double a, H; // Cosmological parameters
    
    virtual void fcnSmooth(GravityParticle *p, int nSmooth,
			   pqSmoothNode *nList);
    virtual void initSmoothParticle(GravityParticle *p);
    virtual void initSmoothCache(GravityParticle *p);
    virtual void combSmoothCache(GravityParticle *p1,
				 ExternalSmoothParticle *p2);
 public:
    DenDvDxSmoothParams() {}
    DenDvDxSmoothParams(int _iType, int am, CSM csm, double dTime) {
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
    }
    PUPable_decl(DenDvDxSmoothParams);
    DenDvDxSmoothParams(CkMigrateMessage *m) : SmoothParams(m) {}
    virtual void pup(PUP::er &p) {
        SmoothParams::pup(p);//Call base class
	p|a;
	p|H;
	}
    };

class PressureSmoothParams : public SmoothParams
{
    double a, H; // Cosmological parameters
    double alpha, beta; // SPH viscosity parameters
    
    virtual void fcnSmooth(GravityParticle *p, int nSmooth,
			   pqSmoothNode *nList);
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

#endif
