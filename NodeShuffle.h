#ifndef NODE_SHUFFLE_H
#define NODE_SHUFFLE_H
/// @file NodeShuffle.h
/// Bin record of the node-aggregated domain decomposition particle
/// exchange (bNodeShuffle).  Kept separate because the generated
/// ParallelGravity.decl.h names it in entry method signatures and is
/// included from DataManager.h before ParallelGravity.h defines the
/// message and holder classes.
#include "pup.h"

/// @brief One destination TreePiece's share of a NodeShuffleBuf.
///
/// Offsets index the arrays of the NodeShuffleBuf that carries the
/// bin.  On the sending side, srcFirst is the index of the first
/// particle of the bin in the source TreePiece's myParticles.
struct ShuffleBin {
    int destPiece;      ///< destination TreePiece
    int destNode;       ///< node (process) of destPiece
    int iBin;           ///< index of this record in the holder's bin table
    int srcFirst;       ///< source-side: first particle in myParticles
    int srcLoad;        ///< source-side: first entry in myShuffleLoads/Parts
    int srcIndex;       ///< source piece slot in an intra-process holder (NodeShuffleBuf::selfSources)
    int iPart, nPart;   ///< range in particles[]
    int iGas, nGas;     ///< range in pGas[]
    int iStar, nStar;   ///< range in pStar[]
    int iLoad, nLoads;  ///< range in loads[] and parts_per_phase[]
};
PUPbytes(ShuffleBin);

#endif
