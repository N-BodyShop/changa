#ifndef __STATE_H__
#define __STATE_H__
#include "ParallelGravity.h"

// less flexible, but probably more efficient
// and allows for sharing of counters between
// states
class State {
  public:
  int *counterArrays[2];
};

typedef CkQ<OffsetNode> CheckList;

typedef CkVec<OffsetNode> UndecidedList;
typedef CkVec<UndecidedList> UndecidedLists;

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
  bool placedRoots;

#ifdef CUDA
  // arrays of CkVecs - interaction lists
  // fixed size
  CkVec<ILCell> cellLists;
  CkVec<ILPart> partLists;
  // nodes and particles to ship to gpu (null for remote-no-resume and local)
  // pointers to expandable arrays for remote-resume
  CkVec<CudaMultipoleMoments> *nodes;
  CkVec<CompactPartData> *particles;

  // number of nodes/particles per request per TreePiece
  int nodeThreshold;
  int partThreshold;

  // fixed size arrays. elements mark out the boundaries of interaction lists of different buckets
  CkVec<int> nodeInteractionBucketMarkers;
  CkVec<int> partInteractionBucketMarkers;

  // fixed size arrays. mark the beginning of bucket particles in an array on the gpu and
  // size of each bucket
  // nodes
  // which buckets are involved in this offload
  CkVec<int> bucketsNodes;
  CkVec<int> bucketStartMarkersNodes;
  CkVec<int> bucketSizesNodes;

  // parts
  // which buckets are involved in this offload
  CkVec<int> bucketsParts;
  CkVec<int> bucketStartMarkersParts;
  CkVec<int> bucketSizesParts;


  // to tell a remote-resume state from a remote-no-resume state
  bool resume;
#endif

  // The lowest nodes reached on paths to each bucket
  // Used to find numBuckets completed when
  // walk returns. Also used to find at which
  // bucket computation should start
  GenericTreeNode *lowestNode;
  int level;
#ifdef CUDA
  int numInteractions;
#endif

  DoubleWalkState() : chklists(0), lowestNode(0), level(-1)
#ifdef CUDA
  , numInteractions(0)
#endif
  {}

};

class NullState : public State {
};

class ListState : public State {
  //public:
  //OffsetNode nodeList;
  //CkVec<LocalPartInfo> localParticleList;
  //CkVec<RemotePartInfo> remoteParticleList;
};

#endif
