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
  PrefetchRequestStruct(OrientedBox<double> *p, int n) : prefetchReq(p), numPrefetchReq(n) {}
};

/* Computes */
class Compute{
  protected: 
  State *state;
  Opt *opt;
  void *computeEntity;
  int activeRung;
  ComputeType type;

  Compute(ComputeType t) : type(t), state(0){}

  public: 
  void setOpt(Opt *opt);
  // should the dowork method have a state argument?
  // yes, allows listcompute object to keep modifying state
  // which will have within it the checklist, clist and plist
  virtual int doWork(GenericTreeNode *, TreeWalk *tw, State *state, int chunk, int reqID, bool isRoot, bool &didcomp, int awi) = 0;
  virtual bool openCriterion(TreePiece *ownerTP, 
                            GenericTreeNode *node, int reqID) = 0;
  // TreeWalk will call this when a request is sent to the CacheManager for
  // a tree node on the Compute's behalf
  // By default, nothing is done.
  // Objects of type PrefetchCompute also call this when they request
  // particles
  // virtual void nodeRequestSent(TreePiece *owner) {};       
  virtual void init(void *cE, int activeRung, Opt *opt);
  virtual void reassoc(void *cE, int aR, Opt *opt){}
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
  virtual void recvdParticles(ExternalGravityParticle *egp,int num,int chunk,int reqID,State *state, TreePiece *tp){}

  virtual ~Compute(){}
  virtual void walkDone(){}
  virtual void deleteComputeEntity(){}
  virtual void setComputeEntity(void *ce){
    computeEntity = ce;
  }

  // Computes are not associated with single buckets anymore
  // This function returns a pointer to the state corresponding to
  // bucket bucketIdx
  // Defaults to returning 0 for each bucket - expected behaviour in 
  // Gravity and Prefetch computes
  virtual State *getResumeState(int bucketIdx);
};

class GravityCompute : public Compute{
  // GenericTreeNode *myBucket;
  // int activeRung;

  public:
  GravityCompute() : Compute(Gravity){}
  
  int doWork(GenericTreeNode *, TreeWalk *tw, State *state, int chunk, int reqID, bool isRoot, bool &didcomp, int awi);

  bool openCriterion(TreePiece *ownerTP, GenericTreeNode *node, int reqID);
  int computeParticleForces(TreePiece *owner, GenericTreeNode *node, ExternalGravityParticle *part, int reqID);
  int computeNodeForces(TreePiece *owner, GenericTreeNode *nd, int reqID);

  // book keeping on notifications
  int startNodeProcessEvent(TreePiece *owner){}
  int finishNodeProcessEvent(TreePiece *owner){}
  int nodeMissedEvent(TreePiece *owner, int chunk);
  
  void recvdParticles(ExternalGravityParticle *egp,int num,int chunk,int reqID,State *state, TreePiece *tp);

  void reassoc(void *cE, int activeRung, Opt *o);

};

class ListCompute : public Compute{

  public:
  ListCompute() : Compute(List){}
  
  int doWork(GenericTreeNode *, TreeWalk *tw, State *state, int chunk, int reqID, bool isRoot, bool &didcomp, int awi);

  bool openCriterion(TreePiece *ownerTP, GenericTreeNode *node, int reqID);

  void addNodeToList(GenericTreeNode *, ListState *);
  void addParticlesToList(ExternalGravityParticle *, int n, ListState *);
  void addParticlesToList(GravityParticle *, int n, ListState *);

  // book keeping on notifications
  int startNodeProcessEvent(TreePiece *owner){}
  int finishNodeProcessEvent(TreePiece *owner){}
  int nodeMissedEvent(TreePiece *owner, int chunk);
};

class PrefetchCompute : public Compute{

  public: 
  PrefetchCompute() : Compute(Prefetch) {
    computeEntity = 0;
  }
  int doWork(GenericTreeNode *, TreeWalk *tw, State *state, int chunk, int reqID, bool isRoot, bool &didcomp, int awi);
  bool openCriterion(TreePiece *ownerTP, GenericTreeNode *node, int reqID);

  // book-keeping on notifications
  int startNodeProcessEvent(TreePiece *owner);
  int finishNodeProcessEvent(TreePiece *owner);
  int nodeMissedEvent(TreePiece *owner, int chunk){}

  //void init(void *cE, int activeRung, Opt *opt);

  void recvdParticles(ExternalGravityParticle *egp,int num,int chunk,int reqID,State *state, TreePiece *tp);
  void deleteComputeEntity(){
  }
  ~PrefetchCompute(){
  }
  /*
  void deleteComputeEntity(){
    delete (PrefetchRequestStruct *)computeEntity;
    CkPrintf("PRS: DELETE\n");
    computeEntity = 0;
  }
  ~PrefetchCompute(){
    if(computeEntity != 0)
      deleteComputeEntity();
  }
  */
  //void walkDone();
};

