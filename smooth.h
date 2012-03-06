#ifndef __SMOOTH_H
#define __SMOOTH_H
/*
 * structure for the priority queue.  Also contains information for
 * smooth calculation.
 */
#include <queue>
#include "Compute.h"
#include "State.h"

/// Object for priority queue entry.
class pqSmoothNode
{
 public:
    double fKey;	// distance^2 -> place in priority queue
    Vector3D<double> dx; // displacement of this particle
    GravityParticle *p; // pointer to rest of particle data
    
    inline bool operator<(const pqSmoothNode& n) const {
	return fKey < n.fKey;
	}
    };

	
/// Object to bookkeep a Bucket Smooth Walk.

class NearNeighborState: public State {
public:
    CkVec<pqSmoothNode> *Qs; 
    int nParticlesPending;
    int mynParts; 
    bool started;
    
    NearNeighborState(int nParts, int nSmooth) {
        Qs = new CkVec<pqSmoothNode>[nParts+2];
	mynParts = nParts; 
        }

    void finishBucketSmooth(int iBucket, TreePiece *tp);
    ~NearNeighborState() {
	delete [] Qs; 
        }
};

/// We can make this a base class from which parameters for all smooth
/// types can be derived.
class SmoothParams : public PUP::able
{
 public:
    int iType;	// Particle type to smooth over;  "TreeActive"
    int activeRung;
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
    /// initialize particles as they come into the cache
    virtual void initSmoothCache(GravityParticle *p) = 0;
    /// combine cache copy with home particle
    virtual void combSmoothCache(GravityParticle *p1,
				 ExternalSmoothParticle *p2) = 0;
    SmoothParams() {}
    PUPable_abstract(SmoothParams);
    SmoothParams(CkMigrateMessage *m) : PUP::able(m) {}
    /// required method for remote entry call.
    virtual void pup(PUP::er &p) {
        PUP::able::pup(p);//Call base class
        p|iType;
        p|activeRung;
	}
    };

// XXX so we can access the init function with the Cache unpack
// method.

extern SmoothParams *globalSmoothParams;

/// Class to specify density smooth
class DensitySmoothParams : public SmoothParams
{
    virtual void fcnSmooth(GravityParticle *p, int nSmooth,
			   pqSmoothNode *nList);
    virtual int isSmoothActive(GravityParticle *p);
    virtual void initSmoothParticle(GravityParticle *p);
    virtual void initTreeParticle(GravityParticle *p) {}
    virtual void postTreeParticle(GravityParticle *p) {}
    virtual void initSmoothCache(GravityParticle *p);
    virtual void combSmoothCache(GravityParticle *p1,
				 ExternalSmoothParticle *p2);
 public:
    DensitySmoothParams() {}
    DensitySmoothParams(int _iType, int am) {
	iType = _iType;
	activeRung = am;
	}
    PUPable_decl(DensitySmoothParams);
    DensitySmoothParams(CkMigrateMessage *m) : SmoothParams(m) {}
    virtual void pup(PUP::er &p) {
        SmoothParams::pup(p);//Call base class
	}
    };

/// Super class for Smooth and Resmooth computation.
class SmoothCompute : public Compute
{
 protected:
    TreePiece *tp;

 public:
    SmoothParams *params;

    SmoothCompute(TreePiece *_tp, SmoothParams *_params) : Compute(Smooth){
	params = _params;
	// XXX Assign to global pointer: not thread safe
	globalSmoothParams = params;
        tp = _tp;       // needed in getNewState()
	}
    ~SmoothCompute() { //delete state;
      // delete params;
    }

    virtual 
      void bucketCompare(TreePiece *tp,
			 GravityParticle *p,  // Particle to test
			 GenericTreeNode *node, // bucket
			 GravityParticle *particles, // local particle data
			 Vector3D<double> offset,
			 State *state
			 ) = 0;

    int doWork(GenericTreeNode *node,
	       TreeWalk *tw,
	       State *state,
	       int chunk,
	       int reqID,
	       bool isRoot, 
	       bool &didcomp, int awi);
    void reassoc(void *ce, int ar, Opt *o);    
    void nodeMissedEvent(int reqID, int chunk, State *state, TreePiece *tp);
    


};

/// Class for computation over k nearest neighbors.
class KNearestSmoothCompute : public SmoothCompute 
{
    int nSmooth;
    // limit small smoothing lengths
    int iLowhFix;
    // smoothing to gravitational softening ratio limit
    double dfBall2OverSoft2;
    
public:
    
    KNearestSmoothCompute(TreePiece *_tp, SmoothParams *_params, int nSm,
			  int iLhF, double dfB2OS2)
       : SmoothCompute(_tp, _params){
         nSmooth = nSm;
	 iLowhFix = iLhF;
	 dfBall2OverSoft2 = dfB2OS2;
         }
    ~KNearestSmoothCompute() { //delete state;
	delete params;
    }

    void bucketCompare(TreePiece *tp,
		       GravityParticle *p,  // Particle to test
		       GenericTreeNode *node, // bucket
		       GravityParticle *particles, // local particle data
		       Vector3D<double> offset,
                       State *state
		       ) ;
	    
    void initSmoothPrioQueue(int iBucket, State *state) ;
    int openCriterion(TreePiece *ownerTP, GenericTreeNode *node, int reqID, State *state);
    void startNodeProcessEvent(State *state){ }
    void finishNodeProcessEvent(TreePiece *owner, State *state){ }
    void nodeRecvdEvent(TreePiece *owner, int chunk, State *state, int bucket);
    void recvdParticlesFull(GravityParticle *egp,int num,int chunk,
			int reqID,State *state, TreePiece *tp,
			Tree::NodeKey &remoteBucket);
    void walkDone(State *state) ;

    // this function is used to allocate and initialize a new state object
    // these operations were earlier carried out in the constructor of the
    // class.
    State *getNewState(int d1);
    // Unused.
    State *getNewState(int d1, int d2) {return 0;}
    State *getNewState() {return 0;}
    };

/// Object to bookkeep a Bucket ReSmooth Walk.  This could be merged
/// with the standard smooth if we changed that to using push_heap()
/// and pop_heap()

class ReNearNeighborState: public State {
public:
    CkVec<pqSmoothNode> *Qs;
    int nParticlesPending;
    bool started;
    ReNearNeighborState(int nParts) {
	Qs = new CkVec<pqSmoothNode>[nParts+2];
	}
    void finishBucketSmooth(int iBucket, TreePiece *tp);
    ~ReNearNeighborState() { delete [] Qs; }
};

class ReSmoothCompute : public SmoothCompute 
{
    
public:
    ReSmoothCompute(TreePiece *_tp, SmoothParams *_params) : SmoothCompute(_tp, _params){}

    ~ReSmoothCompute() { //delete state;
	delete params;
    }

    void bucketCompare(TreePiece *tp,
		       GravityParticle *p,  // Particle to test
		       GenericTreeNode *node, // bucket
		       GravityParticle *particles, // local particle data
		       Vector3D<double> offset,
                       State *state
		       ) ;
	    
    int openCriterion(TreePiece *ownerTP, GenericTreeNode *node, int reqID, State *state);
    void startNodeProcessEvent(State *state){ }
    void finishNodeProcessEvent(TreePiece *owner, State *state){ }
    void nodeRecvdEvent(TreePiece *owner, int chunk, State *state, int bucket);
    void recvdParticlesFull(GravityParticle *egp,int num,int chunk,
			int reqID,State *state, TreePiece *tp,
			Tree::NodeKey &remoteBucket);
    void walkDone(State *state) ;

    // this function is used to allocate and initialize a new state object
    // these operations were earlier carried out in the constructor of the
    // class.
    State *getNewState(int d1);
    State *getNewState(int d1, int d2) {return 0;}
    State *getNewState() {return 0;}
    };

/// Computation over "inverse" nearest neighbors.
class MarkSmoothCompute : public SmoothCompute 
{
    
public:
    MarkSmoothCompute(TreePiece *_tp, SmoothParams *_params) : SmoothCompute(_tp, _params){}

    ~MarkSmoothCompute() { //delete state;
	delete params;
    }

    void bucketCompare(TreePiece *tp,
		       GravityParticle *p,  // Particle to test
		       GenericTreeNode *node, // bucket
		       GravityParticle *particles, // local particle data
		       Vector3D<double> offset,
                       State *state
		       ) ;
	    
    int openCriterion(TreePiece *ownerTP, GenericTreeNode *node, int reqID, State *state);
    void startNodeProcessEvent(State *state){ }
    void finishNodeProcessEvent(TreePiece *owner, State *state){ }
    void nodeRecvdEvent(TreePiece *owner, int chunk, State *state, int bucket);
    void recvdParticlesFull(GravityParticle *egp,int num,int chunk,
			int reqID,State *state, TreePiece *tp,
			Tree::NodeKey &remoteBucket);
    void walkDone(State *state) ;

    // this function is used to allocate and initialize a new state object
    // these operations were earlier carried out in the constructor of the
    // class.
    State *getNewState(int d1, int d2) {return 0;}
    State *getNewState(int d1);
    State *getNewState() {return 0;}
    };

/// Object to bookkeep a Bucket MarkSmooth Walk.

class MarkNeighborState: public State {
public:
    int nParticlesPending;
    bool started;
    MarkNeighborState(int nParts) {}
    void finishBucketSmooth(int iBucket, TreePiece *tp);
    ~MarkNeighborState() {}
};

#include "Opt.h"

class SmoothOpt : public Opt{
  public:
  SmoothOpt() : Opt(Local){
    // don't need to open
    // these nodes are your concern
    action_array[0][Internal] = DUMP;
    action_array[0][Bucket] = DUMP;

    action_array[0][Boundary] = DUMP; 
    action_array[0][NonLocal] = DUMP; 
    action_array[0][NonLocalBucket] = DUMP;	
    action_array[0][Cached] = DUMP;	
    action_array[0][CachedBucket] = DUMP;

    action_array[0][Empty] = DUMP;
    action_array[0][CachedEmpty] = DUMP;
    action_array[0][Top] = ERROR;
    action_array[0][Invalid] = ERROR;
    //--------------
    // need to open node
    // local data
    action_array[1][Internal] = KEEP;
    action_array[1][Bucket] = KEEP_LOCAL_BUCKET;
    action_array[1][Boundary] = KEEP;

    // remote data
    action_array[1][NonLocal] = KEEP;
    action_array[1][NonLocalBucket] = KEEP_REMOTE_BUCKET;
    action_array[1][CachedBucket] = KEEP_REMOTE_BUCKET;
    action_array[1][Cached] = KEEP;

    // discard
    action_array[1][Empty] = DUMP;
    action_array[1][CachedEmpty] = DUMP;
    action_array[1][Top] = ERROR;
    action_array[1][Invalid] = ERROR;
  }

};

/* Standard M_4 Kernel */
/* return 1/(h_smooth)^2 for a particle */
inline
double invH2(GravityParticle *p) 
{
    return 4.0/(p->fBall*p->fBall);
    }

inline double KERNEL(double ar2) 
{
    double ak;
    ak = 2.0 - sqrt(ar2);
    if (ar2 < 1.0) ak = (1.0 - 0.75*ak*ar2);
    else ak = 0.25*ak*ak*ak;
    return ak;
    }
inline double DKERNEL(double ar2) 
{
    double adk;
    adk = sqrt(ar2);
    if (ar2 < 1.0) {
	adk = -3 + 2.25*adk;
	}
    else {
	adk = -0.75*(2.0-adk)*(2.0-adk)/adk;
	}
    return adk;
    }

#endif
