#ifndef PE_LIST_H
#define PE_LIST_H

#ifdef CUDA
#include "ParallelGravity.h"

// To use pinned host memory with a std::vector
template <typename T>
struct PinnedHostAllocator {
    typedef T value_type;

    PinnedHostAllocator() {}

    template <class U>
    PinnedHostAllocator(const PinnedHostAllocator<U>&) {}

    T* allocate(std::size_t n) {
        void* ptr = NULL;
        allocatePinnedHostMemory(&ptr, n * sizeof(T));
        if (!ptr) throw std::bad_alloc();
        return static_cast<T*>(ptr);
    }

    void deallocate(T* p, std::size_t /*n*/) {
        freePinnedHostMemory(p);
    }

    template <typename U>
    struct rebind {
        typedef PinnedHostAllocator<U> other;
    };
};

// Equality operators
template <class T, class U>
bool operator==(const PinnedHostAllocator<T>&, const PinnedHostAllocator<U>&) {
    return true;
}

template <class T, class U>
bool operator!=(const PinnedHostAllocator<T>&, const PinnedHostAllocator<U>&) {
    return false;
}

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

    CudaRequest *request;

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
    CkCallback *finishCb;

    cudatype fperiod;

  public:
    PEList() {finalBucketMarker = -1; cudaStreamCreate(&stream); }
    PEList(CkMigrateMessage *m) : CBase_PEList(m) {}
    ~PEList() { cudaStreamDestroy(stream); }
    void pup(PUP::er &p) {}

    void setType(int _bNode, int _bRemote, int _bResume);

    void finishWalk(TreePiece *treePiece);
    void sendList(TreePiece *treePiece, CudaRequest *data);
    void reset();
};

#endif
#endif
