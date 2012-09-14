#ifndef __STATE_H__
#define __STATE_H__
#include "ParallelGravity.h"

/// @brief Base class for maintaining the state of a tree walk.
class State {
  public:
    /// Set after our walk is finished, but we are still waiting for
    /// combiner cache flushes to be processed.
    int bWalkDonePending;
    /// The bucket we have started to walk.
    int currentBucket;

    // shifted variable into state. there is an issue of redundancy 
    // here, though. in addition to local state, remote and remote-resume
    // state also have this variable but have no use for it, since only
    // a single copy is required.
    // could have made this the third element in the array below
    /// @brief Keep track of how many buckets are unfinished.  XXX
    /// note the misnomer.
    int myNumParticlesPending;

    // again, redundant variables, since only remote-no-resume
    // walks use this variable to see how many chunks have 
    // been used
    ///
    /// @brief Number of pending chunks.  
    ///
    /// The remote tree walk is divided into chunks for more parallelism.
    /// A chunk is pending wrt a TreePiece until that TreePiece has
    /// finished using it completely.
    int numPendingChunks;

    /// @brief counters to keep track of outstanding remote processor
    // requests tied to each bucket (position 0) and chunk (position 1).
    int *counterArrays[2];
    virtual ~State() {}
};

#if INTERLIST_VER > 0
#if defined CUDA
#include "HostCUDA.h"
#include "DataManager.h"

class DoubleWalkState;

template<typename T>
class GenericList{
  public:
  CkVec<CkVec<T> > lists;
  int totalNumInteractions;

  GenericList() : totalNumInteractions(0) {}

  void reset(){
    // clear all bucket lists:
    for(int i = 0; i < lists.length(); i++){
      lists[i].length() = 0;
    }
    totalNumInteractions = 0;
  }

  void free(){
    for(int i = 0; i < lists.length(); i++){
      lists[i].free();
    }
    lists.free();
    totalNumInteractions = 0;
  }

  void init(int numBuckets, int numper){
    lists.resize(numBuckets);
    for(int i = 0; i < numBuckets; i++){
      lists[i].reserve(numper);
    }
  }

  CudaRequest *serialize(TreePiece *tp);
  void getBucketParameters(TreePiece *tp, 
                           int bucket, 
                           int &bucketStart, int &bucketSize){
                           //std::map<NodeKey, int>&lpref){
	// bucket is listed in this offload
	GenericTreeNode *bucketNode = tp->bucketList[bucket];

	bucketSize = bucketNode->lastParticle - bucketNode->firstParticle + 1;
        bucketStart = bucketNode->bucketArrayIndex;
	CkAssert(bucketStart >= 0);
  }

  void getActiveBucketParameters(TreePiece *tp, 
                           int bucket, 
                           int &bucketStart, int &bucketSize){
                           //std::map<NodeKey, int>&lpref){
	// bucket is listed in this offload
	GenericTreeNode *bucketNode = tp->bucketList[bucket];
        BucketActiveInfo *binfo = &(tp->bucketActiveInfo[bucket]);

	//bucketSize = bucketNode->lastParticle - bucketNode->firstParticle + 1;
        //bucketStart = bucketNode->bucketArrayIndex;
        bucketSize = tp->bucketActiveInfo[bucket].size;
        bucketStart = tp->bucketActiveInfo[bucket].start;
	CkAssert(bucketStart >= 0);
  }

  void push_back(int b, T &ilc, DoubleWalkState *state, TreePiece *tp);
  

};

#endif

///
/// @brief Hold state where both the targets and sources are tree walked.
///
class DoubleWalkState : public State {
  public:
  CheckList *chklists;
  UndecidedLists undlists;
  CkVec<CkVec<OffsetNode> >clists;
  CkVec<CkVec<LocalPartInfo> >lplists;
  CkVec<CkVec<RemotePartInfo> >rplists;
   
  // set once before the first cgr is called for a chunk
  // the idea is to place the chunkRoot (along with replicas)
  // on the remote comp chklist only once per chunk
  //
  // one for each chunk
  bool *placedRoots;
  // to tell a remote-resume state from a remote-no-resume state
  bool resume;

#ifdef CUDA
  int nodeThreshold;
  int partThreshold;

  GenericList<ILCell> nodeLists;
  GenericList<ILPart> particleLists;

#ifdef CUDA_INSTRUMENT_WRS
  double nodeListTime;
  double partListTime;
#endif

  CkVec<CudaMultipoleMoments> *nodes;
  CkVec<CompactPartData> *particles;

  // during 'small' rungs, buckets are marked when
  // they are included for computation in the request's
  // aux. particle array. these markings should be
  // cleared before the assembly of the next request is
  // begun. for this purpose, we keep track of buckets
  // marked during the construction of a request.
  //
  // NB: for large rungs, we don't mark buckets while 
  // compiling requests. for such rungs, since all
  // particles are shipped at the beginning of the iteration,
  // we have them marked at that time. since all particles,
  // are available on the gpu for these rungs, we do not clear 
  // the markings when requests are sent out.
  CkVec<GenericTreeNode *> markedBuckets;

  // TODO : this switch from map to ckvec means that we cannot 
  // use multiple treepieces per processor, since they will all
  // be writing to the nodeArrayIndex field of the CacheManager's nodes.
  // We need a different group that manages GPU memory for this purpose.
  //std::map<NodeKey,int> nodeMap;
  CkVec<GenericTreeNode *> nodeMap;
  std::map<NodeKey,int> partMap;

  bool nodeOffloadReady(){
    return nodeLists.totalNumInteractions >= nodeThreshold;
  }

  bool partOffloadReady(){
    return particleLists.totalNumInteractions >= partThreshold;
  }

#ifdef CUDA_INSTRUMENT_WRS
  void updateNodeThreshold(int t){
    nodeThreshold = t;
  }
  void updatePartThreshold(int t){
    partThreshold = t;
  }
#endif

#endif

  // The lowest nodes reached on paths to each bucket
  // Used to find numBuckets completed when
  // walk returns. Also used to find at which
  // bucket computation should start
  GenericTreeNode *lowestNode;
  int level;

  DoubleWalkState() : chklists(0), lowestNode(0), level(-1)
  {}

#ifdef CUDA_INSTRUMENT_WRS
  void nodeListConstructionTimeStart(){
    nodeListTime = CmiWallTimer();
  }

  double nodeListConstructionTimeStop(){
    return CmiWallTimer()-nodeListTime;
  }

  void partListConstructionTimeStart(){
    partListTime = CmiWallTimer();
  }

  double partListConstructionTimeStop(){
    return CmiWallTimer()-partListTime;
  }

#endif
};
#endif //  INTERLIST_VER 

class NullState : public State {
};

#endif
