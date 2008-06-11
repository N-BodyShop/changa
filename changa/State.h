#ifndef __STATE_H__
#define __STATE_H__
#include "ParallelGravity.h"

class State {
};

class DoubleWalkState : public State {
  public: 
  CkQ<Tree::NodeKey> chklist;
  CkQ<Tree::NodeKey> clist;
  CkQ<Tree::NodeKey> plist;
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
