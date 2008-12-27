/*
 * structure for the priority queue.  Also contains information for
 * smooth calculation.
 */
#include <queue>
#include "Compute.h"

class pqSmoothNode
{
 public:
    double fKey;	// distance squared -> place in priority queue
    Vector3D<double> dx; // displacement of this particle
    // XXX to be replaced with "ExternalSmoothParticle"
    ExternalGravityParticle *p; // pointer to rest of particle data
    
    inline bool operator<(const pqSmoothNode& n) const {
	return fKey < n.fKey;
	}
    };

	
// Object to bookkeep a Bucket Smooth Walk.

class NearNeighborState: public State {
public:
    std::priority_queue <pqSmoothNode> *Qs;
    int nParticlesPending;
    bool started;
    NearNeighborState(int nParts) {
	Qs = new (std::priority_queue<pqSmoothNode>[nParts+2]);
	}
    void finishBucketSmooth(int iBucket, TreePiece *tp);
    ~NearNeighborState() { delete [] Qs; }
};

class SmoothCompute : public Compute 
{
    int nSmooth;
    void (*fcnSmooth)(GravityParticle *p, int nSmooth, pqSmoothNode *nList);
    void bucketCompare(TreePiece *tp,
		       ExternalGravityParticle *p,  // Particle to test
		       GenericTreeNode *node, // bucket
		       GravityParticle *particles, // local particle data
		       Vector3D<double> offset,
                       State *state
		       ) ;
    TreePiece *tp;
    
public:
 SmoothCompute(TreePiece *_tp, void (*fcn)(GravityParticle *p, int nSmooth,
			   pqSmoothNode *nList),
	       void (*_fcnInit)(ExternalGravityParticle *p),
	       void (*_fcnCombine)(GravityParticle *p1, ExternalGravityParticle *p2),
	       int nSm)
     : Compute(Smooth){
	nSmooth = nSm;
	fcnSmooth = fcn;
        tp = _tp;       // needed in getNewState()
	}
    ~SmoothCompute() { //delete state;
    }
	    
    void initSmoothPrioQueue(int iBucket, State *state) ;

    int openCriterion(TreePiece *ownerTP, GenericTreeNode *node, int reqID, State *state);
    
    int doWork(GenericTreeNode *node,
	       TreeWalk *tw,
	       State *state,
	       int chunk,
	       int reqID,
	       bool isRoot, 
	       bool &didcomp, int awi);
    int startNodeProcessEvent(TreePiece *owner){}
    int finishNodeProcessEvent(TreePiece *owner, State *state){}
    int nodeMissedEvent(int reqID, int chunk, State *state, TreePiece *tp);
    int nodeRecvdEvent(TreePiece *owner, int chunk, State *state, int bucket);
    void recvdParticles(ExternalGravityParticle *egp,int num,int chunk,
			int reqID,State *state, TreePiece *tp,
			Tree::NodeKey &remoteBucket);
    void reassoc(void *ce, int ar, Opt *o);
    void walkDone(State *state) ;

    // this function is used to allocate and initialize a new state object
    // these operations were earlier carried out in the constructor of the
    // class.
    State *getNewState(int d1);

    };

// Object to bookkeep a Bucket ReSmooth Walk.  This could be merged
// with the standard smooth if we changed that to using push_heap()
// and pop_heap()

class ReNearNeighborState: public State {
public:
    std::vector <pqSmoothNode> *Qs;
    int nParticlesPending;
    bool started;
    ReNearNeighborState(int nParts) {
	Qs = new (std::vector<pqSmoothNode>[nParts+2]);
	}
    void finishBucketSmooth(int iBucket, TreePiece *tp);
    ~ReNearNeighborState() { delete [] Qs; }
};

class ReSmoothCompute : public Compute 
{
    void (*fcnSmooth)(GravityParticle *p, int nSmooth, pqSmoothNode *nList);
    void bucketCompare(TreePiece *tp,
		       ExternalGravityParticle *p,  // Particle to test
		       GenericTreeNode *node, // bucket
		       GravityParticle *particles, // local particle data
		       Vector3D<double> offset,
                       State *state
		       ) ;
    TreePiece *tp;
    
public:
 ReSmoothCompute(TreePiece *_tp, void (*fcn)(GravityParticle *p, int nSmooth,
					     pqSmoothNode *nList),
		 void (*_fcnInit)(ExternalGravityParticle *p),
		 void (*_fcnCombine)(GravityParticle *p1, ExternalGravityParticle *p2)
)
     : Compute(Smooth){
	fcnSmooth = fcn;
        tp = _tp;       // needed in getNewState()
	}
    ~ReSmoothCompute() { //delete state;
    }
	    
    int openCriterion(TreePiece *ownerTP, GenericTreeNode *node, int reqID, State *state);
    
    int doWork(GenericTreeNode *node,
	       TreeWalk *tw,
	       State *state,
	       int chunk,
	       int reqID,
	       bool isRoot, 
	       bool &didcomp, int awi);
    int startNodeProcessEvent(TreePiece *owner){}
    int finishNodeProcessEvent(TreePiece *owner, State *state){}
    int nodeMissedEvent(int reqID, int chunk, State *state, TreePiece *tp);
    int nodeRecvdEvent(TreePiece *owner, int chunk, State *state, int bucket);
    void recvdParticles(ExternalGravityParticle *egp,int num,int chunk,
			int reqID,State *state, TreePiece *tp,
			Tree::NodeKey &remoteBucket);
    void reassoc(void *ce, int ar, Opt *o);
    void walkDone(State *state) ;

    // this function is used to allocate and initialize a new state object
    // these operations were earlier carried out in the constructor of the
    // class.
    State *getNewState(int d1);
 
     };
 

void Density(GravityParticle *p,int nSmooth,pqSmoothNode *nnList);
void DensitySym(GravityParticle *p,int nSmooth,pqSmoothNode *nnList);
void initDensity(ExternalGravityParticle *p) ;
void combDensity(GravityParticle *p1, ExternalGravityParticle *p2);

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
