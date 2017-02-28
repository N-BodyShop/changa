#ifndef __SMOOTHPARAMS_H
#define __SMOOTHPARAMS_H
class pqSmoothNode;

/// @brief A base class from which parameters for all smooth
/// operations can be derived.
class SmoothParams : public PUP::able
{
 public:
    int iType;	///< Particle type to smooth over;  "TreeActive"
    int activeRung;  ///< Currently active rung
    int bUseBallMax;	///< limit fBall growth for bFastGas
    /// Function to apply to smooth particle and neighbors
    virtual void fcnSmooth(GravityParticle *p, int nSmooth, pqSmoothNode *nList) = 0;
    /// Particle is doing a neighbor search
    virtual int isSmoothActive(GravityParticle *p) = 0;
    /// initialize particles to be smoothed
    virtual void initSmoothParticle(GravityParticle *p) = 0;
    /// initialize particles in tree but not smoothed
    virtual void initTreeParticle(GravityParticle *p) = 0;
    /// calculation on all tree particles after all walks are done
    virtual void postTreeParticle(GravityParticle *p) = 0;
    /// @brief initialize particles as they come into the cache
    /// @param p pointer to incoming particle data from remote processor
    virtual void initSmoothCache(GravityParticle *p) = 0;
    /// @brief combine cache copy with home particle
    /// @param p1 pointer to "home" particle data
    /// @param p2 pointer to particle data being flushed back from a
    /// remote treepiece.
    ///
    /// This method enables commutative/associative operations to be
    /// performed on remote data by combining data accumulated by a
    /// remote processor in p2 with the data on the home particle, p1.
    /// Any accumulators used in this function should be initialized
    /// in initSmoothCache() to avoid double counting.
    virtual void combSmoothCache(GravityParticle *p1,
				 ExternalSmoothParticle *p2) = 0;
    // limit ball growth by default
    SmoothParams() { bUseBallMax = 1; }
    PUPable_abstract(SmoothParams);
    SmoothParams(CkMigrateMessage *m) : PUP::able(m) {}
    /// required method for remote entry call.
    virtual void pup(PUP::er &p) {
        PUP::able::pup(p);//Call base class
        p|iType;
        p|activeRung;
	p|bUseBallMax;
	}
    };
#endif
