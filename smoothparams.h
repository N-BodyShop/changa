#ifndef __SMOOTHPARAMS_H
#define __SMOOTHPARAMS_H
class pqSmoothNode;

// We can make this a base class from which parameters for all smooth
// types can be derived.
class SmoothParams : public PUP::able
{
 public:
    int iType;	// Particle type to smooth over;  "TreeActive"
    int activeRung;
    // Function to apply to smooth particle and neighbors
    virtual void fcnSmooth(GravityParticle *p, int nSmooth, pqSmoothNode *nList) = 0;
    // Particle is doing a neighbor search
    virtual int isSmoothActive(GravityParticle *p) = 0;
    // initialize particles to be smoothed
    virtual void initSmoothParticle(GravityParticle *p) = 0;
    // initialize particles in tree but not smoothed
    virtual void initTreeParticle(GravityParticle *p) = 0;
    // calculation on all tree particles after all walks are done
    virtual void postTreeParticle(GravityParticle *p) = 0;
    // initialize particles as they come into the cache
    virtual void initSmoothCache(GravityParticle *p) = 0;
    // combine cache copy with home particle
    virtual void combSmoothCache(GravityParticle *p1,
				 ExternalSmoothParticle *p2) = 0;
    SmoothParams() {}
    PUPable_abstract(SmoothParams);
    SmoothParams(CkMigrateMessage *m) : PUP::able(m) {}
    virtual void pup(PUP::er &p) {
        PUP::able::pup(p);//Call base class
        p|iType;
        p|activeRung;
	}
    };
#endif
