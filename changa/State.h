#ifndef __STATE_H__
#define __STATE_H__
#include "ParallelGravity.h"

// less flexible, but probably more efficient
// and allows for sharing of counters between
// states
class State {
  public:
  int *counterArrays[2];
  virtual ~State() {}
};

#if INTERLIST_VER > 0
typedef CkQ<OffsetNode> CheckList;

typedef CkVec<OffsetNode> UndecidedList;
typedef CkVec<UndecidedList> UndecidedLists;

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
                           int &bucketStart, int &bucketSize, 
                           std::map<NodeKey, int>&lpref){
	// bucket b is listed in this offload
	GenericTreeNode *bucketNode = tp->bucketList[bucket];

	bucketSize = bucketNode->lastParticle - bucketNode->firstParticle + 1;
	NodeKey partKey = bucketNode->getKey();
	partKey <<= 1;

	std::map<NodeKey, int>::iterator iter = lpref.find(partKey);
	CkAssert(iter != lpref.end());
	bucketStart = iter->second;
	CkAssert(bucketStart >= 0);
  }

  void push_back(int b, T &ilc, DoubleWalkState *state, TreePiece *tp);
  

};

#endif

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

  CkVec<CudaMultipoleMoments> *nodes;
  CkVec<CompactPartData> *particles;

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
#endif

  // The lowest nodes reached on paths to each bucket
  // Used to find numBuckets completed when
  // walk returns. Also used to find at which
  // bucket computation should start
  GenericTreeNode *lowestNode;
  int level;

  DoubleWalkState() : chklists(0), lowestNode(0), level(-1)
  {}

};


#if defined CUDA
void allocatePinnedHostMemory(void **ptr, int size);

template<typename T>
CudaRequest *GenericList<T>::serialize(TreePiece *tp){
    // for local particles
    std::map<NodeKey, int> &lpref = tp->dm->getLocalPartsOnGpuTable();

    // get count of buckets with interactions first
    int numFilledBuckets = 0;
    int listpos = 0;
    int curbucket = 0;

    double starttime = CmiWallTimer();
    for(int i = 0; i < lists.length(); i++){
      if(lists[i].length() > 0){
        numFilledBuckets++;
      }
    }

    // create flat lists and associated data structures
#ifdef CUDA_USE_CUDAMALLOCHOST
    T *flatlists;
    int *markers, *starts, *sizes;

    allocatePinnedHostMemory((void **)&flatlists, totalNumInteractions*sizeof(T));
    allocatePinnedHostMemory((void **)&markers, (numFilledBuckets+1)*sizeof(int));
    allocatePinnedHostMemory((void **)&starts, (numFilledBuckets)*sizeof(int));
    allocatePinnedHostMemory((void **)&sizes, (numFilledBuckets)*sizeof(int));
#else
    T *flatlists = (T *) malloc(totalNumInteractions*sizeof(T));
    int *markers = (int *) malloc((numFilledBuckets+1)*sizeof(int));
    int *starts = (int *) malloc(numFilledBuckets*sizeof(int));
    int *sizes = (int *) malloc(numFilledBuckets*sizeof(int));
#endif
    int *affectedBuckets = new int[numFilledBuckets];

    // populate flat lists
    int listslen = lists.length();
    for(int i = 0; i < listslen; i++){
      int listilen = lists[i].length();
      if(listilen > 0){
        memcpy(&flatlists[listpos], lists[i].getVec(), listilen*sizeof(T));
        markers[curbucket] = listpos;
        getBucketParameters(tp, i, starts[curbucket], sizes[curbucket], lpref);
        affectedBuckets[curbucket] = i;
        listpos += listilen;
        curbucket++;
      }
    }
    markers[numFilledBuckets] = listpos;
    CkAssert(listpos == totalNumInteractions);

    CudaRequest *request = new CudaRequest;
    request->list = (void *)flatlists;
    request->bucketMarkers = markers;
    request->bucketStarts = starts;
    request->bucketSizes = sizes;
    request->numInteractions = totalNumInteractions;
    request->numBucketsPlusOne = numFilledBuckets+1;
    request->affectedBuckets = affectedBuckets;
    request->tp = (void *)tp;
    request->fperiod = tp->fPeriod.x;

    traceUserBracketEvent(CUDA_SER_LIST, starttime, CmiWallTimer());

    return request;
  }

//#include <typeinfo>
template<typename T>
void GenericList<T>::push_back(int b, T &ilc, DoubleWalkState *state, TreePiece *tp){
    if(lists[b].length() == 0){
        state->counterArrays[0][b]++;
#if COSMO_PRINT_BK > 1
        CkPrintf("[%d] request out bucket %d numAddReq: %d,%d\n", tp->getIndex(), b, tp->sInterListStateRemote->counterArrays[0][b], tp->sInterListStateLocal->counterArrays[0][b]);
#endif
    }
    lists[b].push_back(ilc);
    totalNumInteractions++;
  }
#endif // CUDA
#endif //  INTERLIST_VER 

class NullState : public State {
};

class ListState : public State {
  //public:
  //OffsetNode nodeList;
  //CkVec<LocalPartInfo> localParticleList;
  //CkVec<RemotePartInfo> remoteParticleList;
};

#endif
