#ifndef PE_TREE_MERGER_H
#define PE_TREE_MERGER_H

#include "ParallelGravity.h"

class PETreeMerger : public CBase_PETreeMerger {

  CkVec<GenericTreeNode*> submittedRoots;
  CkVec<TreePiece*> submittedTreePieces;

  NonEmptyTreePieceCounter submittingTreePieceCounter;

  public:
  PETreeMerger() {}
  PETreeMerger(CkMigrateMessage *) {}
  void pup(PUP::er &p) {}

  void mergeNonLocalRequests(GenericTreeNode *root, TreePiece *treePiece);
  void freeTree();

  private:

  GenericTreeNode *mergeWalk(CkVec<GenericTreeNode*> &mergeList, CkVec<TreePiece*> &treePieceList);
  void requestNonLocalMoments(GenericTreeNode *pickedNode, TreePiece *pickedTreePiece);

};

#endif
