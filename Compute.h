#ifndef __COMPUTE_H__
#define __COMPUTE_H__

#include "codes.h"
#include "ParallelGravity.h"

class State;
class DoubleWalkState;

/** @file Compute.h
 * Defines classes for objects that encapsulate computation
 */

class TreeWalk;
class Opt;

/// this is the computeEntity for PrefetchComputes
/// it holds an array of prefetch root bounding boxes
/// and the number of elements in this array

struct PrefetchRequestStruct{
  OrientedBox<double> *prefetchReq;
  int numPrefetchReq;
  PrefetchRequestStruct(OrientedBox<double> *p, int n) : prefetchReq(p), numPrefetchReq(n) {}
};

/// @brief Base clase for all tree based computations.
///
/// The Compute object determines what work is to be done at each
/// treenode, as well as what gets done at the beginning and the end
/// of a walk.  The key method is doWork().
class Compute{
  protected:
  Opt *opt;
  void *computeEntity;
  int activeRung;
  ComputeType type;

  Compute(ComputeType t) : type(t) /*state(0)*/{}

  public:
  int nActive;  // accumulate total number of active particles.
  void setOpt(Opt *opt);
  // should the dowork method have a state argument?
  // yes, allows listcompute object to keep modifying state
  // which will have within it the checklist, clist and plist
  /// @brief Work to be done at each node.
  virtual int doWork(GenericTreeNode *, TreeWalk *tw, State *state, int chunk, int reqID, bool isRoot, bool &didcomp, int awi) = 0;
  // should return int, not bool
  virtual int openCriterion(TreePiece *ownerTP,
                            GenericTreeNode *node, int reqID, State *state) = 0;

  virtual void stateReady(State *state, TreePiece *owner, int chunk, int start, int end) {}

  virtual void stateReadyPar(TreePiece *tp, int start, int end,
    CkVec<OffsetNode>& clist, CkVec<RemotePartInfo>& rpilist,
    CkVec<LocalPartInfo>& lpilist){}

  virtual void fillLists(State *state_, TreePiece *tp, int chunk, int start,
    int end, CkVec<OffsetNode>& clistforb, CkVec<RemotePartInfo>& rplistforb,
    CkVec<LocalPartInfo>& lplistforb) {}

  /// @brief Associate computeEntity (target bucket or node),
  /// activeRung and Optimization with this Compute object.
  virtual void init(void *cE, int activeRung, Opt *opt);
  virtual void reassoc(void *cE, int aR, Opt *opt){}
  ComputeType getSelfType(){ return type;}
  OptType getOptType();
  int getActiveRung() {return activeRung;}

  // Default impl is empty. Currently only redefined by ListCompute
  // Initializes state.
  virtual void initState(State *state){}
  /// virtual functions to allow for book-keeping
  /// these are essentially notifications to the
  /// Compute object from the TreeWalk that certain
  /// events have taken place - the Compute reacts
  /// accordingly.
  virtual void startNodeProcessEvent(State *state) {}
  /// Allow book-keeping when finished with a node
  virtual void finishNodeProcessEvent(TreePiece *owner, State *state) {}
  /// Allow book-keeping of a cache miss.
  virtual void nodeMissedEvent(int reqID, int chunk, State *state, TreePiece *tp) {}
  /// Allow book-keeping of a cache receive event.
  virtual void nodeRecvdEvent(TreePiece *owner, int chunk, State *state, int bucket){}
  /// Allow book-keeping of a cache receive event.
  virtual void recvdParticles(ExternalGravityParticle *egp,int num,int chunk,int reqID,State *state, TreePiece *tp, Tree::NodeKey &remoteBucket){}
  /// Allow book-keeping of a cache receive event.
  virtual void recvdParticlesFull(GravityParticle *egp,int num,int chunk,int reqID,State *state, TreePiece *tp, Tree::NodeKey &remoteBucket){}
  virtual ~Compute(){}
  virtual void walkDone(State *state){}
  virtual void setComputeEntity(void *ce){
    computeEntity = ce;
  }

  virtual void *getComputeEntity(){
    return computeEntity;
  }

  virtual State *getNewState(int d1, int d2);
  virtual State *getNewState(int d1);
  virtual State *getNewState();
  virtual void freeState(State *state);
};

#include "SSEdefs.h"

///
/// @brief Class to compute gravity using a "bucket walk".
///
class GravityCompute : public Compute{
#ifdef BENCHMARK_TIME_COMPUTE
  double computeTimePart;
  double computeTimeNode;
#endif

  void updateInterMass(GravityParticle *p, int start, int end, double mass);
  void updateInterMass(GravityParticle *p, int start, int end, GravityParticle *s, Vector3D<cosmoType> &offset);

  public:
  GravityCompute() : Compute(Gravity){
#ifdef BENCHMARK_TIME_COMPUTE
    computeTimePart = 0.0;
    computeTimeNode = 0.0;
#endif
  }
  ~GravityCompute() {
#ifdef BENCHMARK_TIME_COMPUTE
    CkPrintf("Compute time part: %f\n",computeTimePart);
    CkPrintf("Compute time node: %f\n",computeTimeNode);
#endif
  }

  int doWork(GenericTreeNode *, TreeWalk *tw, State *state, int chunk, int reqID, bool isRoot, bool &didcomp, int awi);

  int openCriterion(TreePiece *ownerTP, GenericTreeNode *node, int reqID, State *state);
  int computeParticleForces(TreePiece *owner, GenericTreeNode *node, ExternalGravityParticle *part, int reqID);

  // book keeping on notifications
  void nodeMissedEvent(int reqID, int chunk, State *state, TreePiece *tp);
  void nodeRecvdEvent(TreePiece *owner, int chunk, State *state, int bucket);
  void recvdParticles(ExternalGravityParticle *egp,int num,int chunk,int reqID,State *state, TreePiece *tp, Tree::NodeKey &remoteBucket);

  void reassoc(void *cE, int activeRung, Opt *o);

};

#if INTERLIST_VER > 0
/// @brief Computations for Stadel-style interaction list walk.
///
/// At a given point in the walk, this compares a node against another
/// see if it 1) has to be opened, 2) doesn't need to be opened or 3)
/// undecided, and manipulates the lists in State accordingly.
class ListCompute : public Compute{

  public:
  ListCompute() : Compute(List) {}

  int doWork(GenericTreeNode *, TreeWalk *tw, State *state, int chunk, int reqID, bool isRoot, bool &didcomp, int awi);

  int openCriterion(TreePiece *ownerTP, GenericTreeNode *node, int reqID, State *state);

  // book keeping on notifications
  void nodeMissedEvent(int reqID, int chunk, State *state, TreePiece *tp);
  void nodeRecvdEvent(TreePiece *owner, int chunk, State *state, int bucket);
  void recvdParticles(ExternalGravityParticle *egp,int num,int chunk,int reqID,State *state, TreePiece *tp, Tree::NodeKey &remoteBucket);

  void initState(State *state);

  void stateReady(State *, TreePiece *, int chunk, int start, int end);

  void stateReadyPar(TreePiece *tp, int start, int end,
      CkVec<OffsetNode>& clist, CkVec<RemotePartInfo>& rpilist,
      CkVec<LocalPartInfo>& lpilist);

  void fillLists(State *state_, TreePiece *tp, int chunk, int start,
    int end, CkVec<OffsetNode>& clistforb, CkVec<RemotePartInfo>& rplistforb,
    CkVec<LocalPartInfo>& lplistforb);

//  void printUndlist(DoubleWalkState *state, int level, TreePiece *tp);
//  void printClist(DoubleWalkState *state, int level, TreePiece *tp);
  void reassoc(void *cE, int activeRung, Opt *o);
  State *getNewState(int d1, int d2);
  State *getNewState(int d1);
  State *getNewState();
  void freeState(State *state);
  void freeDoubleWalkState(DoubleWalkState *state);

#ifdef CUDA
#ifdef GPU_LOCAL_TREE_WALK
  void sendLocalTreeWalkTriggerToGpu(State *state, TreePiece *tp, int activeRung, int startBucket, int endBucket);
#endif //GPU_LOCAL_TREE_WALK
  void sendNodeInteractionsToGpu(DoubleWalkState *state, TreePiece *tp);
  void sendPartInteractionsToGpu(DoubleWalkState *state, TreePiece *tp);
#endif

  private:

  void addChildrenToCheckList(GenericTreeNode *node, int reqID, int chunk, int awi, State *s, CheckList &chklist, TreePiece *tp);
  void addNodeToInt(GenericTreeNode *node, int offsetID, DoubleWalkState *s);

  DoubleWalkState *allocDoubleWalkState();

#if defined CHANGA_REFACTOR_PRINT_INTERACTIONS || defined CHANGA_REFACTOR_WALKCHECK_INTERLIST || defined CUDA
  void addRemoteParticlesToInt(ExternalGravityParticle *parts, int n,
			       Vector3D<cosmoType> &offset, DoubleWalkState *s,
			       NodeKey key);
  void addLocalParticlesToInt(GravityParticle *parts, int n,
			      Vector3D<cosmoType> &offset, DoubleWalkState *s,
			      NodeKey key, GenericTreeNode *gtn);
#else
  void addRemoteParticlesToInt(ExternalGravityParticle *parts, int n,
			       Vector3D<cosmoType> &offset, DoubleWalkState *s);
  void addLocalParticlesToInt(GravityParticle *parts, int n,
			      Vector3D<cosmoType> &offset, DoubleWalkState *s);
#endif

#ifdef CUDA
  void getBucketParameters(TreePiece *tp, int bucket, int &bucketStart, int &bucketSize, std::map<NodeKey, int>&lpref);

  public:
  void resetCudaNodeState(DoubleWalkState *state);
  void resetCudaPartState(DoubleWalkState *state);
  void initCudaState(DoubleWalkState *state, int numBuckets, int nodeThreshold, int partThreshold, bool resume);
  
#endif

};
#endif

/// @brief Compute for the remote node prefetch walk.
class PrefetchCompute : public Compute{

  public:
  PrefetchCompute() : Compute(Prefetch) {
    computeEntity = 0;
  }
  virtual int doWork(GenericTreeNode *, TreeWalk *tw, State *state, int chunk, int reqID, bool isRoot, bool &didcomp, int awi);
  int openCriterion(TreePiece *ownerTP, GenericTreeNode *node, int reqIDD, State *state);

  // book-keeping on notifications
  void startNodeProcessEvent(State *state);
  void finishNodeProcessEvent(TreePiece *owner, State *state);

  void recvdParticles(ExternalGravityParticle *egp,int num,int chunk,int reqID,State *state, TreePiece *tp, Tree::NodeKey &remoteBucket);
};

/// when prefetching is disabled from command line, use 
/// DummyPrefetchCompute instead of PrefetchCompute
class DummyPrefetchCompute : public PrefetchCompute {
  public:
    /// Immediately stop the walk.
  int doWork(GenericTreeNode *, TreeWalk *tw, State *state, int chunk, int reqID, bool isRoot, bool &didcomp, int awi){
    return DUMP;
  }
};

/// @brief distingish between the walks that could be running.

enum WalkIndices {
    prefetchAwi = 0,
    interListAwi = 1,
    remoteGravityAwi = 2,
    smoothAwi = 3,
    maxAwi = 4
};
    
/// Object to record a type of active walk. Contains pointers to
/// TreeWalk/Compute/Opt/State (T/C/O/S) combinations
class ActiveWalk {
  public:
  TreeWalk *tw;
  Compute *c;
  Opt *o;
  State *s;
  
  ActiveWalk(TreeWalk *_tw, Compute *_c, Opt *_o, State *state) : 
      tw(_tw), c(_c), o(_o), s(state){}
  ActiveWalk(){}
};

/// @brief Interface for work in LocalTreeTraversal.
class TreeNodeWorker {

  public:
  virtual bool work(GenericTreeNode *node, int level) = 0;
  virtual void doneChildren(GenericTreeNode *node, int level) {}
};

/// @brief Class to build the remote part of the tree.  Fills in
/// Boundary and NonLocal nodes.
class RemoteTreeBuilder : public TreeNodeWorker {
  TreePiece *tp;
  bool requestNonLocalMoments; 

  public:
  RemoteTreeBuilder(TreePiece *owner, bool req) : 
    tp(owner),
    requestNonLocalMoments(req)
  {}

  bool work(GenericTreeNode *node, int level);
  void doneChildren(GenericTreeNode *node, int level);

  private:
  void registerNode(GenericTreeNode *node);

};

/// @brief Class to build the local part of the tree.  Builds Internal nodes.
class LocalTreeBuilder : public TreeNodeWorker {
  TreePiece *tp;

  public:
  LocalTreeBuilder(TreePiece *owner) :
    tp(owner)
  {}

  bool work(GenericTreeNode *node, int level);
  void doneChildren(GenericTreeNode *node, int level);

  private:
  void registerNode(GenericTreeNode *node);
};

/** @brief TreeNodeWorker implementation that just prints out the
 * tree.  This is just for diagnostics.
 */
class LocalTreePrinter : public TreeNodeWorker {
  int index;
  std::ofstream file;
  string description;

  void openFile();

  public:
  LocalTreePrinter(string d, int idx) : 
    index(idx),
    description(d)
  {
    openFile();
  }

  ~LocalTreePrinter(){
    file << "}" << std::endl;
    file.close();
  }

  bool work(GenericTreeNode *node, int level);
  void doneChildren(GenericTreeNode *node, int level);
};

#endif
