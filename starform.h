#ifndef STARFORM_HINCLUDED
#define STARFORM_HINCLUDED

#include "parameters.h"
#include "smooth.h"

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
#endif
