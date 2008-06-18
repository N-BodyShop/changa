#ifndef __COMPUTE_H__
#define __COMPUTE_H__

#include "codes.h"
class State;
#include "State.h"

/* File defines classes for objects that encapsulate computation
 * */ 

class TreeWalk;
class Opt;

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
  //State *state;
  Opt *opt;
  void *computeEntity;
  int activeRung;
  ComputeType type;

  Compute(ComputeType t) : type(t) /*state(0)*/{}

  public: 
  void setOpt(Opt *opt);
  // should the dowork method have a state argument?
  // yes, allows listcompute object to keep modifying state
  // which will have within it the checklist, clist and plist
  virtual int doWork(GenericTreeNode *, TreeWalk *tw, State *state, int chunk, int reqID, bool isRoot, bool &didcomp, int awi) = 0;
  virtual bool openCriterion(TreePiece *ownerTP, 
                            GenericTreeNode *node, int reqID, State *state) = 0;
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
  virtual int startNodeProcessEvent(State *state) {}
  virtual int finishNodeProcessEvent(TreePiece *owner, State *state) {}
  //virtual int nodeMissedEvent(TreePiece *owner, int chunk) = 0;
  virtual int nodeMissedEvent(int reqID, int chunk, State *state) {}
  virtual int nodeRecvdEvent(TreePiece *owner, int chunk, State *state, int bucket){}
  virtual void recvdParticles(ExternalGravityParticle *egp,int num,int chunk,int reqID,State *state, TreePiece *tp){}

  virtual ~Compute(){}
  virtual void walkDone(State *state){}
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
  virtual State *getNewState(); 
};

class GravityCompute : public Compute{
  // GenericTreeNode *myBucket;
  // int activeRung;

  public:
  GravityCompute() : Compute(Gravity){}
  
  int doWork(GenericTreeNode *, TreeWalk *tw, State *state, int chunk, int reqID, bool isRoot, bool &didcomp, int awi);

  bool openCriterion(TreePiece *ownerTP, GenericTreeNode *node, int reqID, State *state);
  int computeParticleForces(TreePiece *owner, GenericTreeNode *node, ExternalGravityParticle *part, int reqID);
  int computeNodeForces(TreePiece *owner, GenericTreeNode *nd, int reqID);

  // book keeping on notifications
  int nodeMissedEvent(int reqID, int chunk, State *state);
  int nodeRecvdEvent(TreePiece *owner, int chunk, State *state, int bucket);
  
  void recvdParticles(ExternalGravityParticle *egp,int num,int chunk,int reqID,State *state, TreePiece *tp);

  void reassoc(void *cE, int activeRung, Opt *o);
  State *getNewState();

};

class ListCompute : public Compute{

  public:
  ListCompute() : Compute(List){}
  
  int doWork(GenericTreeNode *, TreeWalk *tw, State *state, int chunk, int reqID, bool isRoot, bool &didcomp, int awi);

  bool openCriterion(TreePiece *ownerTP, GenericTreeNode *node, int reqIDD, State *state);

  void addNodeToList(GenericTreeNode *, ListState *);
  void addParticlesToList(ExternalGravityParticle *, int n, ListState *);
  void addParticlesToList(GravityParticle *, int n, ListState *);

  // book keeping on notifications
  int nodeMissedEvent(int reqID, int chunk, State *state);
};

class PrefetchCompute : public Compute{

  public: 
  PrefetchCompute() : Compute(Prefetch) {
    computeEntity = 0;
  }
  int doWork(GenericTreeNode *, TreeWalk *tw, State *state, int chunk, int reqID, bool isRoot, bool &didcomp, int awi);
  bool openCriterion(TreePiece *ownerTP, GenericTreeNode *node, int reqIDD, State *state);

  // book-keeping on notifications
  int startNodeProcessEvent(State *state);
  int finishNodeProcessEvent(TreePiece *owner, State *state);

  void recvdParticles(ExternalGravityParticle *egp,int num,int chunk,int reqID,State *state, TreePiece *tp);
  
  State *getNewState();
};

#endif
