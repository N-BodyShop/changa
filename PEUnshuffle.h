#ifndef PE_UNSHUFFLE_H
#define PE_UNSHUFFLE_H

#include "ParallelGravity.h"

/// @brief Group to perform unshuffle across an entire PE in order to
/// reduce the number of messages used.
class PEUnshuffle : public CBase_PEUnshuffle 
{
    /// TreePieces on this PE
    CkVec<TreePiece*> vtpLocal;
    /// Count of TreePieces on this PE
    NonEmptyTreePieceCounter cTreePieces;

  public:
    PEUnshuffle() {}
    PEUnshuffle(CkMigrateMessage *m) : CBase_PEUnshuffle(m) {}
    void pup(PUP::er &p) { CBase_PEUnshuffle::pup(p); }

    void sendParticlesDuringDD(TreePiece *treePiece);
};
#endif
