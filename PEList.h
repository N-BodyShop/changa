#ifndef PE_LIST_H
#define PE_LIST_H

#ifdef CUDA
#include "ParallelGravity.h"

// To use pinned host memory with a std::vector
template <typename T>
struct PinnedHostAllocator
{
	typedef T value_type;

	T* allocate(size_t n) {
            size_t size = n*sizeof(T);
	    void *ptr = NULL;
	    allocatePinnedHostMemory(&ptr, size);
	    return reinterpret_cast<T*>(ptr);
	}

	void deallocate(T* p, size_t n) {
            freePinnedHostMemory(p);
	}
};

// Need a separate class for each type of request (particle, node, local, remote)
class PEList : public CBase_PEList
{
    /// TreePieces on this PE
    CkVec<TreePiece*> vtpLocal;
    /// Count of TreePieces with particles on this PE
    NonEmptyTreePieceCounter cTreePieces;

    vector<ILCell, PinnedHostAllocator<ILCell>> iList;

    vector<int, PinnedHostAllocator<int>> bucketMarkers;
    int finalBucketMarker;
    vector<int, PinnedHostAllocator<int>> bucketStarts;
    vector<int, PinnedHostAllocator<int>> bucketSizes;

    vector<CompactPartData, PinnedHostAllocator<CompactPartData>> missedParts;
    vector<CudaMultipoleMoments, PinnedHostAllocator<CudaMultipoleMoments>> missedNodes;

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
    PEList() {finalBucketMarker = -1; cudaStreamCreate(&stream); }
    PEList(CkMigrateMessage *m) : CBase_PEList(m) {}
    ~PEList() { cudaStreamDestroy(stream); }
    void pup(PUP::er &p) {}

    void setType(int _bNode, int _bRemote, int _bResume);

    void finishWalk(TreePiece *treePiece, CkCallback &cb);
    void sendList(TreePiece *treePiece, CudaRequest *data);
    void reset();
};

#endif
#endif
