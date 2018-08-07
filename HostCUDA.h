#ifndef _HOST_CUDA_H_
#define _HOST_CUDA_H_


#include <cuda_runtime.h>
#include "cuda_typedef.h"

#ifdef HAPI_USE_CUDAMALLOCHOST
# ifdef HAPI_MEMPOOL
#  define CUDA_MALLOC(ptr,sz) ptr = hapiPoolMalloc(size)
# else
#  define CUDA_MALLOC(ptr,sz) cudaMallocHost(&(ptr), size)
# endif
#else
# define CUDA_MALLOC(ptr,sz) ptr = malloc(sz)
#endif

#define THREADS_PER_BLOCK 128

#ifdef GPU_LOCAL_TREE_WALK
#define THREADS_PER_WARP 32
#define WARPS_PER_BLOCK (THREADS_PER_BLOCK / THREADS_PER_WARP)
#define WARP_INDEX (threadIdx.x >> 5)
#endif //GPU_LOCAL_TREE_WALK

#ifdef CUDA_2D_TB_KERNEL
#define PARTS_PER_BLOCK 16
#define NODES_PER_BLOCK (THREADS_PER_BLOCK/PARTS_PER_BLOCK)

#define THREADS_PER_BLOCK_PART 128
#define PARTS_PER_BLOCK_PART 16
#define NODES_PER_BLOCK_PART (THREADS_PER_BLOCK_PART/PARTS_PER_BLOCK_PART)
#endif

// FIXME - find appropriate values
#define NUM_INIT_MOMENT_INTERACTIONS_PER_BUCKET 100
#define NUM_INIT_PARTICLE_INTERACTIONS_PER_BUCKET 100

/* defines for Hybrid API buffer indices */
#define LOCAL_MOMENTS        0
#define LOCAL_PARTICLE_CORES  1
#define LOCAL_PARTICLE_VARS      2
#define REMOTE_MOMENTS 3 
#define REMOTE_PARTICLE_CORES 4

#define LOCAL_MOMENTS_IDX        0
#define LOCAL_PARTICLE_CORES_IDX  1
#define LOCAL_PARTICLE_VARS_IDX      2
#define REMOTE_MOMENTS_IDX 0
#define REMOTE_PARTICLE_CORES_IDX 1

#define ILPART 0
#define PART_BUCKET_MARKERS 1
#define PART_BUCKET_START_MARKERS 2
#define PART_BUCKET_SIZES 3
#define ILCELL 0
#define NODE_BUCKET_MARKERS 1
#define NODE_BUCKET_START_MARKERS 2
#define NODE_BUCKET_SIZES 3

#define ILPART_IDX 0
#define PART_BUCKET_MARKERS_IDX 1
#define PART_BUCKET_START_MARKERS_IDX 2
#define PART_BUCKET_SIZES_IDX 3
#define ILCELL_IDX 0
#define NODE_BUCKET_MARKERS_IDX 1
#define NODE_BUCKET_START_MARKERS_IDX 2
#define NODE_BUCKET_SIZES_IDX 3

#define MISSED_MOMENTS 4
#define MISSED_PARTS 4

#define MISSED_MOMENTS_IDX 4
#define MISSED_PARTS_IDX 4

// node moments, particle cores, particle vars
#define DM_TRANSFER_LOCAL_NBUFFERS 3
#define DM_TRANSFER_REMOTE_CHUNK_NBUFFERS 2

// interaction list
// list markers
// bucket starts
// bucket sizes
#define TP_GRAVITY_LOCAL_NBUFFERS 4
#define TP_GRAVITY_LOCAL_NBUFFERS_SMALLPHASE 5

#define TP_NODE_GRAVITY_REMOTE_NBUFFERS 4
#define TP_PART_GRAVITY_REMOTE_NBUFFERS 4

#define TP_NODE_GRAVITY_REMOTE_RESUME_NBUFFERS 5
#define TP_PART_GRAVITY_REMOTE_RESUME_NBUFFERS 5

#define MAX_NBUFFERS 5

// tp_gravity_local uses arrays of particles and nodes already allocated on the gpu
// tp_gravity_remote uses arrays of nodes already on the gpu + particles from an array it supplies
// tp_gravity_remote_resume uses an array each of nodes and particles it supplies
enum kernels {
  DM_TRANSFER_LOCAL=0,
  DM_TRANSFER_REMOTE_CHUNK,
  DM_TRANSFER_BACK,
  DM_TRANSFER_FREE_LOCAL,
  DM_TRANSFER_FREE_REMOTE_CHUNK,
  TP_GRAVITY_LOCAL,
  TP_GRAVITY_REMOTE,
  TP_GRAVITY_REMOTE_RESUME,
  TP_PART_GRAVITY_LOCAL,
  TP_PART_GRAVITY_LOCAL_SMALLPHASE,
  TP_PART_GRAVITY_REMOTE,
  TP_PART_GRAVITY_REMOTE_RESUME,
  EWALD_KERNEL
};


/** @brief Data and parameters for requesting gravity calculations on
 * the GPU. */
typedef struct _CudaRequest{
        /// can either be a ILCell* or an ILPart*
	void *list;
	int *bucketMarkers;     /**< index in the cell or particle
                                 * list for the cells/particles for a
                                 * given bucket */
	int *bucketStarts;      /**< index in the particle array on
                                 * the GPU of the bucket */
	int *bucketSizes;       /**< number of particles in a bucket
                                 * on the GPU  */
	int numInteractions;    /**< Total number of interactions in
                                 * this request  */
	int numBucketsPlusOne;  /**< Number of buckets affected by
                                 *  this request */
        void *tp;               /**< Pointer to TreePiece that made
                                 * this request */

	/// these buckets were finished in this work request
	int *affectedBuckets;
        void *cb;               /**< Callback  */
        void *state;            /**< Pointer to state of walk that
                                 * created this request  */
        cudatype fperiod;       /**< Size of periodic volume  */

        // TODO: remove these later if we don't use COSMO_PRINT_BK.
        /// is this a node or particle computation request?
        bool node;
        /// is this a remote or local computation?
        bool remote;
#ifdef HAPI_INSTRUMENT_WRS
        int tpIndex;
        char phase;
#endif
#ifdef GPU_LOCAL_TREE_WALK
  int firstParticle;
  int lastParticle;
  int rootIdx;
  cosmoType theta;
  cosmoType thetaMono;
  int nReplicas;
  cudatype fperiodY;  // Support periodic boundary condition in more dimensions
  cudatype fperiodZ;  // Support periodic boundary condition in more dimensions
#endif //GPU_LOCAL_TREE_WALK
}CudaRequest;

/** @brief Parameters for the GPU gravity calculations */
typedef struct _ParameterStruct{
  int numInteractions;          /**< Size of the interaction list  */
  int numBucketsPlusOne;        /**< Number of buckets affected (+1
                                 * for fencepost)  */
  cudatype fperiod;             /**< Period of the volume  */
#ifdef GPU_LOCAL_TREE_WALK
  int firstParticle;
  int lastParticle;
  int rootIdx;
  cudatype theta;
  cudatype thetaMono;
  int nReplicas;
  cudatype fperiodY;  // Support periodic boundary condition in more dimensions
  cudatype fperiodZ;  // Support periodic boundary condition in more dimensions
#endif //GPU_LOCAL_TREE_WALK
}ParameterStruct;

#ifdef HAPI_INSTRUMENT_WRS
void DataManagerTransferLocalTree(CudaMultipoleMoments *moments, int nMoments,
                        CompactPartData *compactParts, int nCompactParts,
                        int mype, char phase, void *wrCallback);
void DataManagerTransferRemoteChunk(CudaMultipoleMoments *moments, int nMoments, CompactPartData *compactParts, int nCompactParts, int mype, char phase);
void FreeDataManagerLocalTreeMemory(bool freemom, bool freepart, int pe, char phase);
void FreeDataManagerRemoteChunkMemory(int , void *, bool freemom, bool freepart, int pe, char phase);
void TransferParticleVarsBack(VariablePartData *hostBuffer, size_t size, void *cb, bool, bool, bool, bool, int pe, char phase);
#else
void DataManagerTransferLocalTree(CudaMultipoleMoments *moments, int nMoments,
                        CompactPartData *compactParts, int nCompactParts,
                        int mype, void *wrCallback);
void DataManagerTransferRemoteChunk(CudaMultipoleMoments *moments, int nMoments, CompactPartData *compactParts, int nCompactParts, void *wrCallback);
void FreeDataManagerLocalTreeMemory(bool freemom, bool freepart);
void FreeDataManagerRemoteChunkMemory(int , void *, bool freemom, bool freepart);
/** @brief Transfer forces from the GPU back to the host.
 *  @param hostBuffer Buffer to store results.
 *  @param size hostBuffer size.
 *  @param cb Callback when transfer is done.
 *  @param freemom Boolean: free device buffer with local moment data.
 *  @param freepart Boolean: free device buffer with local particle data.
 *  @param freeRemoteMom Boolean: free device buffer with remote
 *  moment data.
 *  @param freeRemotePart Boolean: free device buffer with remote
 *  particle data.
 */
void TransferParticleVarsBack(VariablePartData *hostBuffer, size_t size, void *cb,
    bool freemom, bool freepart, bool freeRemoteMom, bool freeRemotePart);
#endif

void TreePieceCellListDataTransferLocal(CudaRequest *data);
void TreePieceCellListDataTransferRemote(CudaRequest *data);
void TreePieceCellListDataTransferRemoteResume(CudaRequest *data, CudaMultipoleMoments *missedMoments, int numMissedMoments);


void TreePiecePartListDataTransferLocal(CudaRequest *data);
void TreePiecePartListDataTransferLocalSmallPhase(CudaRequest *data, CompactPartData *parts, int len);
void TreePiecePartListDataTransferRemote(CudaRequest *data);
void TreePiecePartListDataTransferRemoteResume(CudaRequest *data, CompactPartData *missedParticles, int numMissedParticles);

void DummyKernel(void *cb);

#endif
