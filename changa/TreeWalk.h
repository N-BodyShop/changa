#ifndef __TREEWALK_H__
#define __TREEWALK_H__

#include "codes.h"

class State;
class Compute;
class TreePiece;


class TreeWalk{
  protected:
  Compute *comp;
  TreePiece *ownerTP;
  WalkType type;
    

  Tree::NodeKey getChildKey(int i);

  TreeWalk(Compute *_comp, TreePiece *tp, WalkType _type): ownerTP(tp), comp(_comp), type(_type){}
  TreeWalk() : comp(NULL), ownerTP(NULL), type(InvalidWalk){}
  TreeWalk(WalkType t) : comp(NULL), ownerTP(NULL), type(t){}

  public: 
  // must tell compute the ownerTP so it can perform its openCriterion() test
  TreePiece *getOwnerTP(){return ownerTP;}
  Compute *getCompute(){return comp;}
  virtual void init(Compute *c, TreePiece *owner);
  void reassoc(Compute *c);

  virtual void reset() {}; 
  virtual void walk(GenericTreeNode *node, State *state, int chunk, int reqID, int activeWalkIndex) = 0;
  // when resuming a walk after a missed node is received, use this function
  virtual void resumeWalk(GenericTreeNode *node, State *state, int chunk, int reqID, int activeWalkIndex) {
	  walk(node, state, chunk, reqID, activeWalkIndex);
  }
  // beware of using default implementation - always returns 'false'
  //virtual bool ancestorCheck(GenericTreeNode *node, int reqID) {return false;};
  WalkType getSelfType() {return type;} 

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
};

class BottomUpTreeWalk : public TreeWalk{ 
  public: 
  BottomUpTreeWalk(Compute *_comp, TreePiece *tp):TreeWalk(_comp,tp,BottomUp){}
  BottomUpTreeWalk() : TreeWalk(BottomUp) {}

  void walk(GenericTreeNode *node, State *state, int chunk, int reqID, int awi);
};


#if INTERLIST_VER > 0
// Traverses local tree instead of global one.
// Heads towards a target local bucket until told to DUMP.
// The class is used in the local tree part of the 
// double (interaction list) walk. 
// It is given a new target node and a DoubleWalkState object
// whose checklists and undecidedlists it uses while walking 
// down toward the new target node.
// When it returns, the walk function will have modified 
// the state object to reflect what the next target bucket
// should be.
class LocalTargetWalk : public TreeWalk {

  NodeKey targetKey;
  private:
  void dft(GenericTreeNode *localNode, State *state, int chunk, int reqID, 
                                              bool isRoot, int awi, int level);
  bool processNode(
                   GenericTreeNode *glblNode,
                   State *state, 
                   int chunk, int reqID, 
                   bool isRoot, bool &didcomp, 
                   int awi);
  public:
  LocalTargetWalk(Compute *_comp, TreePiece *tp):TreeWalk(_comp,tp,LocalTarget){}
  LocalTargetWalk() : TreeWalk(LocalTarget) { targetKey = 0;}
  void walk(GenericTreeNode *startAncestor, State *state, int chunk, int reqID, int awi);
  void resumeWalk(GenericTreeNode *node, State *state, int chunk, int reqID, int activeWalkIndex);
  NodeKey getTargetKey() {return targetKey;}
  
};
#endif

#endif
