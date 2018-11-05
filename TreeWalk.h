#ifndef __TREEWALK_H__
#define __TREEWALK_H__

#include "codes.h"

class State;
class Compute;
class TreePiece;

/// @brief Base class for walking trees
class TreeWalk{
  protected:
  Compute *comp;
  TreePiece *ownerTP;
  WalkType type;
#ifdef BENCHMARK_TIME_WALK
  double walkTime, keepTime, finishNodeTime, doWorkTime;
#endif

  TreeWalk(Compute *_comp, TreePiece *tp, WalkType _type): ownerTP(tp), comp(_comp), type(_type){
#ifdef BENCHMARK_TIME_WALK
    walkTime = keepTime = finishNodeTime = doWorkTime = 0.0;
#endif
  }
  TreeWalk() : comp(NULL), ownerTP(NULL), type(InvalidWalk){
#ifdef BENCHMARK_TIME_WALK
    walkTime = keepTime = finishNodeTime = doWorkTime = 0.0;
#endif
  }
  TreeWalk(WalkType t) : comp(NULL), ownerTP(NULL), type(t){
#ifdef BENCHMARK_TIME_WALK
    walkTime = keepTime = finishNodeTime = doWorkTime = 0.0;
#endif
  }
  
  public: 
    virtual ~TreeWalk() {
#ifdef BENCHMARK_TIME_WALK
      CkPrintf("walk,keep,finishNode,doWork time: %lf %lf %lf %lf\n",walkTime,keepTime,finishNodeTime,doWorkTime);
#endif
    }
    
  // must tell compute the ownerTP so it can perform its openCriterion() test
  TreePiece *getOwnerTP(){return ownerTP;}
  Compute *getCompute(){return comp;}
  /// @brief Associate a compute object and a treepiece with this walk.
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

/// @brief Walk a tree starting with the root node.
class TopDownTreeWalk : public TreeWalk{ 
  private:
#ifndef CHANGA_REFACTOR_WALKCHECK
  void dft(GenericTreeNode *node, State *state, int chunk, int reqID, bool isRoot, int awi);
  void bft(GenericTreeNode *node, State *state, int chunk, int reqID, bool isRoot, int awi);
#else
  void dft(GenericTreeNode *node, State *state, int chunk, int reqID, bool isRoot, int shift, bool doprint);
#endif
  public: 
  TopDownTreeWalk(Compute *_comp, TreePiece *tp):TreeWalk(_comp,tp,TopDown){}
  TopDownTreeWalk() : TreeWalk(TopDown) {}

  //bool ancestorCheck(GenericTreeNode *node, int reqID);
  void walk(GenericTreeNode *node, State *state, int chunk, int reqID, int awi);
};

/// Walk a tree starting with a leaf node and working up the tree.
/// This class is used for the k-th nearest neighbor search (Smooth).

class BottomUpTreeWalk : public TreeWalk{ 
  public: 
  BottomUpTreeWalk(Compute *_comp, TreePiece *tp):TreeWalk(_comp,tp,BottomUp){}
  BottomUpTreeWalk() : TreeWalk(BottomUp) {}

  void walk(GenericTreeNode *node, State *state, int chunk, int reqID, int awi);
};


/// Traverses local tree instead of global one.
/// Heads towards a target local bucket until told to DUMP.
/// The class is used in the local tree part of the 
/// double (interaction list) walk. 
/// It is given a new target node and a DoubleWalkState object
/// whose checklists and undecidedlists it uses while walking 
/// down toward the new target node.
/// When it returns, the walk function will have modified 
/// the state object to reflect what the next target bucket
/// should be.
class LocalTargetWalk : public TreeWalk {

  NodeKey targetKey;
  private:
  /// @brief Depth First TreeWalk
  /// @param localNode Node to work on.
  /// @param state DoubleWalkState containing check lists.
  /// @param chunk Chunk of remote tree we are working on.
  /// @param reqID Target bucket index
  /// @param isRoot Is localNode the root.
  /// @param awi Active Walk Index.
  /// @param level Level in the tree.
  void dft(GenericTreeNode *localNode, State *state, int chunk, int reqID, 
                                              bool isRoot, int awi, int level);
  /// @brief process a node from the check list
  /// @brief glblNode Source node to be checked.
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

/// @brief class to walk just the local treepiece.
class LocalTreeTraversal {

  public:
  /// Depth First Treewalk.
  /// @param node Node on which we are working.
  /// @param worker Describes work to be done on each node.
  /// @param level Level of the tree we are on.
  void dft(GenericTreeNode *node, TreeNodeWorker *worker, int level){
    if(worker->work(node,level)){
      for(int i = 0; i < node->numChildren(); i++){
        dft(node->getChildren(i),worker,level+1);
      }
      if(node->numChildren() > 0) worker->doneChildren(node,level);
    }
  }
};

#endif
