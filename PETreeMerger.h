#ifndef PE_TREE_MERGER_H
#define PE_TREE_MERGER_H

#include "ParallelGravity.h"

/// @addtogroup TreeBuild
/// @{

/// @brief Group to coordinate requests for remote moments during the
/// tree build
///
/// This group is used if MERGE_REMOTE_REQUESTS is defined.
class PETreeMerger : public CBase_PETreeMerger {

  CkVec<GenericTreeNode*> submittedRoots;
  CkVec<TreePiece*> submittedTreePieces;

  NonEmptyTreePieceCounter submittingTreePieceCounter;

  public:
  PETreeMerger() {}
  PETreeMerger(CkMigrateMessage *m) : CBase_PETreeMerger(m) {}
  void pup(PUP::er &p) { CBase_PETreeMerger::pup(p); }

  void mergeNonLocalRequests(GenericTreeNode *root, TreePiece *treePiece);
  void freeTree();

  private:

  void mergeWalk(CkVec<GenericTreeNode*> &mergeList, CkVec<TreePiece*> &treePieceList);
  void requestNonLocalMoments(GenericTreeNode *pickedNode, TreePiece *pickedTreePiece);

};
/// @}
#endif
