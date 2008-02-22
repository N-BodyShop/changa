/* File defines classes for objects that encapsulate computation
 * */ 

class TreeWalk;
class Opt;
class State;

// this is the computeEntity for PrefetchComputes
// it holds an array of prefetch root bounding boxes
// and the number of elements in this array

struct PrefetchRequestStruct{
  OrientedBox<double> *prefetchReq;
  int numPrefetchReq;
};

/* Computes */
class Compute{
  protected: 
  Opt *opt;
  void *computeEntity;
  int activeRung;
  ComputeType type;

  Compute(ComputeType t) : type(t) {}

  public: 
  void setOpt(Opt *opt);
  // should the dowork method have a state argument?
  // yes, allows listcompute object to keep modifying state
  // which will have within it the checklist, clist and plist
  virtual int doWork(GenericTreeNode *, TreeWalk *tw, State *state, int chunk, int reqID, bool isRoot, bool &didcomp) = 0;
  virtual bool openCriterion(TreePiece *ownerTP, 
                            GenericTreeNode *node, int reqID) = 0;
  // TreeWalk will call this when a request is sent to the CacheManager for
  // a tree node on the Compute's behalf
  // By default, nothing is done.
  // Objects of type PrefetchCompute also call this when they request
  // particles
  // virtual void nodeRequestSent(TreePiece *owner) {};       
  void init(void *cE, int activeRung, Opt *opt);
  //virtual void init(void *cE, int activeRung, Opt *opt) = 0;
  ComputeType getSelfType(){ return type;}
  OptType getOptType();
  int getActiveRung() {return activeRung;}

  // virtual functions to allow for book-keeping
  // these are essentially notifications to the
  // Compute object from the TreeWalk that certain
  // events have taken place - the Compute reacts
  // accordingly.
  virtual int startNodeProcessEvent(TreePiece *owner) = 0;
  virtual int finishNodeProcessEvent(TreePiece *owner) = 0;
  virtual int nodeMissedEvent(TreePiece *owner, int chunk) = 0;


};

class GravityCompute : public Compute{
  // GenericTreeNode *myBucket;
  // int activeRung;

  public:
  GravityCompute() : Compute(Gravity){}
  
  int doWork(GenericTreeNode *, TreeWalk *tw, State *state, int chunk, int reqID, bool isRoot, bool &didcomp);

  bool openCriterion(TreePiece *ownerTP, GenericTreeNode *node, int reqID);
  int computeParticleForces(TreePiece *owner, GenericTreeNode *node, ExternalGravityParticle *part, int reqID);
  int computeNodeForces(TreePiece *owner, GenericTreeNode *nd, int reqID);

  // book keeping on notifications
  int startNodeProcessEvent(TreePiece *owner){
  }
  int finishNodeProcessEvent(TreePiece *owner){
  }
  inline int nodeMissedEvent(TreePiece *owner, int chunk){
    if(getOptType() == Remote){
      owner->addToRemainingChunk(chunk, +1);
    }
  }
  
};

class ListCompute : public Compute{

  public:
  ListCompute() : Compute(List){}
  
  int doWork(GenericTreeNode *, TreeWalk *tw, State *state, int chunk, int reqID, bool isRoot, bool &didcomp);

  bool openCriterion(TreePiece *ownerTP, GenericTreeNode *node, int reqID);

  void addNodeToList(GenericTreeNode *, ListState *);
  void addParticlesToList(ExternalGravityParticle *, int n, ListState *);
  void addParticlesToList(GravityParticle *, int n, ListState *);

  // book keeping on notifications
  int startNodeProcessEvent(TreePiece *owner){
  }
  int finishNodeProcessEvent(TreePiece *owner){
  }
  inline int nodeMissedEvent(TreePiece *owner, int chunk){
    if(getOptType() == Remote){
      owner->addToRemainingChunk(chunk, +1);
    }
  }
  
};

class PrefetchCompute : public Compute{

  public: 
  PrefetchCompute() : Compute(Prefetch) {}
  int doWork(GenericTreeNode *, TreeWalk *tw, State *state, int chunk, int reqID, bool isRoot, bool &didcomp);
  bool openCriterion(TreePiece *ownerTP, GenericTreeNode *node, int reqID);

  // book keeping on notifications
  inline int startNodeProcessEvent(TreePiece *owner){
    return owner->incPrefetchWaiting();
  }
  inline int finishNodeProcessEvent(TreePiece *owner){
    int save = owner->decPrefetchWaiting();
    if(save == 0){
      owner->startRemoteChunk();
    }    
    return save;
  }
  inline int nodeMissedEvent(TreePiece *owner, int chunk){
  }

};


class SPHCompute : public Compute{
  public: 
  SPHCompute() : Compute(SPH){}
};
