#ifndef PE_LIST_H
#define PE_LIST_H

#ifdef CUDA
#include "ParallelGravity.h"

/// @brief PE-level group to coordinate gravity interaction list calculations
/// on the GPU
///
/// There is one instance of this class for each type of walk (particle/node,
/// local, remote, remote-resume). Walks are done at the TreePiece level and
/// the data from each walk is consolidated via the sendList function. Once
/// all of the walks have finished, the finishWalk routine moves the
/// interaction lists to the GPU and starts the computation.
///
/// When GPU_LOCAL_TREE_WALK is enabled, only the remote and remote-resume
/// walks are handled here.
class PEList : public CBase_PEList
{
    /// TreePieces on this PE
    std::vector<TreePiece*> vtpLocal;
    /// Count of TreePieces with particles on this PE
    NonEmptyTreePieceCounter cTreePieces;

    vector<ILCell> iList;

    vector<int> bucketMarkers;
    int finalBucketMarker;
    vector<int> bucketStarts;
    vector<int> bucketSizes;

    vector<CompactPartData> missedParts;
    vector<CudaMultipoleMoments> missedNodes;

    ILCell* iListHost;
    CompactPartData* missedPartsHost;
    CudaMultipoleMoments* missedNodesHost;
    int *bucketMarkersHost;
    int *bucketStartsHost;
    int *bucketSizesHost;

    CudaRequest *request;

    /// Type of request
    int bNode;
    int bRemote;
    int bResume;

    /// Indicate whether remote prefetch data has transferred to the GPU
    int bRemoteReady;
    /// Flags whether the GPU kernel launch was delayed due to the data transfer
    int bKernelDelayed;

    CudaMultipoleMoments *d_localMoments;
    CompactPartData *d_localParts;
    VariablePartData *d_localVars;
    CompactPartData *d_remoteParts;
    CudaMultipoleMoments *d_remoteMoments;
    cudaStream_t stream;
    CkCallback *finishCb;

    cudatype fperiod;

  public:
    PEList(int _bNode, int _bRemote, int _bResume) {
        bNode = _bNode;
        bRemote = _bRemote;
        bResume = _bResume;

	bKernelDelayed = 0;

	finalBucketMarker = -1;
	request = nullptr;
	cudaStreamCreate(&stream);
    }
    PEList(CkMigrateMessage *m) : CBase_PEList(m) {}
    ~PEList() { cudaStreamDestroy(stream); }
    void pup(PUP::er &p) {}

    void finishWalk(TreePiece *treePiece);
    void launchKernel();
    void finishWalkCb();
    void tryLaunchDelayedKernel();
    void sendList(TreePiece *treePiece, CudaRequest *data);
    void reset();
};

#endif
#endif
