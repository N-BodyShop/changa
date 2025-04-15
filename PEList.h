#ifndef PE_LIST_H
#define PE_LIST_H

#include "ParallelGravity.h"

// Need a separate class for each type of request (particle, node, local, remote)
class PEList : public CBase_PEList
{
    /// TreePieces on this PE
    CkVec<TreePiece*> vtpLocal;
    /// Count of TreePieces with particles on this PE
    NonEmptyTreePieceCounter cTreePieces;

    CkVec<ILCell> iList;

    CkVec<int> bucketMarkers;
    int finalBucketMarker;
    CkVec<int> bucketStarts;
    CkVec<int> bucketSizes;

    CkVec<CompactPartData> missedParts;
    CkVec<CudaMultipoleMoments> missedNodes;

    /// Type of request
    int bNode;
    int bRemote;
    int bResume;

    CudaMultipoleMoments *d_localMoments;
    CompactPartData *d_localParts;
    VariablePartData *d_localVars;
    CompactPartData *d_remoteParts;
    CudaMultipoleMoments *d_remoteMoments;
    cudaStream_t stream;

    cudatype fperiod;

  public:
    PEList() {finalBucketMarker = -1;}
    PEList(CkMigrateMessage *m) : CBase_PEList(m) {}
    void pup(PUP::er &p) {}

    void setType(int _bNode, int _bRemote, int _bResume);

    void finishWalk(TreePiece *treePiece, CkCallback &cb);
    void sendList(TreePiece *treePiece, CudaRequest *data);
    void reset();
};

#endif
