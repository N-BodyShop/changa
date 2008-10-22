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

  // The lowest nodes reached on paths to each bucket
  // Used to find numBuckets completed when
  // walk returns. Also used to find at which 
  // bucket computation should start
  GenericTreeNode *lowestNode;
  int level;

  DoubleWalkState() : chklists(0), lowestNode(0), level(-1) {}
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
