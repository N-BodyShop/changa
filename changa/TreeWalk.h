#ifndef __TREEWALK_H__
#define __TREEWALK_H__

#include "codes.h"

class State;
class Compute;
class TreePiece;

/*
class MissRecord{
  public:
  // these fields are saved when the node is missed
  Tree::NodeKey key;            
  State *state;
  int chunk;
  int activeRung;

  WalkType walktype;
  ComputeType computetype;
  OptType opttype;

  MissRecord(Tree::NodeKey & _key, State *_state, int _chunk, int ar, WalkType wt, ComputeType ct, OptType ot) : key(_key), state(_state), chunk(_chunk), activeRung(ar), walktype(wt), computetype(ct), opttype(ot){}
  MissRecord() : state(NULL), chunk(-2), activeRung(-2), walktype(InvalidWalk), computetype(InvalidCompute), opttype(InvalidOpt){}
};
*/

class TreeWalk{
  protected:
  Compute *comp;
  // State *state;
  TreePiece *ownerTP;
  WalkType type;
  // global tree node
  Tree::NodeKey currentGlobalKey;
  // local tree node: do we really need this? - interaction lists?
  GenericTreeNode *currentLocal;

  Tree::NodeKey getChildKey(int i);

  TreeWalk(Compute *_comp, TreePiece *tp, WalkType _type): ownerTP(tp), comp(_comp), type(_type){}
  TreeWalk() : comp(NULL), ownerTP(NULL), type(InvalidWalk){}
  TreeWalk(WalkType t) : comp(NULL), ownerTP(NULL), type(t) {}

  public: 
  // must tell compute the ownerTP so it can perform its openCriterion() test
  TreePiece *getOwnerTP(){return ownerTP;}
  Compute *getCompute(){return comp;}
  virtual void init(Compute *c, TreePiece *owner);
  void reassoc(Compute *c);

  virtual void reset() = 0; 
  virtual void walk(GenericTreeNode *node, State *state, int chunk, int reqID, int activeWalkIndex) = 0;
  // beware of using default implementation - always returns 'false'
  //virtual bool ancestorCheck(GenericTreeNode *node, int reqID) {return false;};
  WalkType getSelfType() {return type;} 
  //virtual void recvdParticles(ExternalGravityParticles *part, int num, int chunk, int reqID, State *state){}

};

class TopDownTreeWalk : public TreeWalk{ 
  private:
#ifndef CHANGA_REFACTOR_WALKCHECK
  void dft(GenericTreeNode *node, State *state, int chunk, int reqID, bool isRoot, int awi);
#else
  void dft(GenericTreeNode *node, State *state, int chunk, int reqID, bool isRoot, int shift, bool doprint);
#endif
  public: 
  TopDownTreeWalk(Compute *_comp, TreePiece *tp):TreeWalk(_comp,tp,TopDown){}
  TopDownTreeWalk() : TreeWalk(TopDown) {}

  //bool ancestorCheck(GenericTreeNode *node, int reqID);
  void walk(GenericTreeNode *node, State *state, int chunk, int reqID, int awi);
  void reset();
};

class BottomUpTreeWalk : public TreeWalk{ 
  public: 
  BottomUpTreeWalk(Compute *_comp, TreePiece *tp):TreeWalk(_comp,tp,BottomUp){}
  BottomUpTreeWalk() : TreeWalk(BottomUp) {}

  void walk(GenericTreeNode *node, State *state, int chunk, int reqID, int awi);
  void reset();
};

/*
class DoubleWalk : public TreeWalk{
  protected:
  // void doubleDft(GenericTreeNode *ancestor, State *state, ... ); 
};
*/
/*
class BucketIteratorWalk : public TreeWalk{
  private:
  void iterate();

  public: 

  BucketIteratorWalk(Compute *_comp, TreePiece *tp):TreeWalk(_comp,tp,BucketIterator){}
  BucketIteratorWalk() : TreeWalk(BucketIterator) {}
  
  void init(Compute *c, TreePiece *owner);

  void reset() = 0; 
  void walk(GenericTreeNode *node, State *state, int chunk, int reqID) = 0;
}
*/

#endif
