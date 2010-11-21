#ifndef _HOST_CUDA_H_
#define _HOST_CUDA_H_


#include "cuda_typedef.h"
#include "cuda_runtime.h"
/* Boolean defines */
enum boolean {NO, YES};

#define THREADS_PER_BLOCK 128

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
#define LOCAL_PARTICLE_CORES  0
#define LOCAL_PARTICLE_VARS      1

#define LOCAL_PARTICLE_CORES_IDX  0
#define LOCAL_PARTICLE_VARS_IDX      1

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
#define PART_OFFSETS_IDX 4
#define ILCELL_IDX 0
#define NODE_BUCKET_MARKERS_IDX 1
#define NODE_BUCKET_START_MARKERS_IDX 2
#define NODE_BUCKET_SIZES_IDX 3
#define NODE_OFFSETS_IDX 4

#define MISSED_MOMENTS 4
#define MISSED_PARTS 4

#define MISSED_MOMENTS_IDX 4
#define MISSED_PARTS_IDX 4

// node moments, particle cores, particle vars
#define DM_TRANSFER_LOCAL_NBUFFERS 2
#define DM_TRANSFER_REMOTE_CHUNK_NBUFFERS 0

// interaction list
// list markers
// bucket starts
// bucket sizes
#define TP_GRAVITY_LOCAL_NBUFFERS 5
#define TP_GRAVITY_LOCAL_NBUFFERS_SMALLPHASE 5

#define TP_NODE_GRAVITY_REMOTE_NBUFFERS 5
#define TP_PART_GRAVITY_REMOTE_NBUFFERS 5

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
  TOP_EWALD_KERNEL,
  BOTTOM_EWALD_KERNEL
};


typedef struct _CudaRequest{
        // can be CudaMultipoleMoments or CompactPartData list 
	void *list;
        int *offsets;
	int *bucketMarkers;
	int *bucketStarts;
	int *bucketSizes;

	int numInteractions;
	int numBucketsPlusOne;
        void *tp;

	// these buckets were finished in this work request
	int *affectedBuckets;
        void *cb;
        void *state;
        // call dummy kernel and initiate freeing of remote chunk
        // memory on GPU
        bool callDummy;
        cudatype fperiod;

        // TODO: remove these later. is this a node or particle computation request?
        bool node;
        // is this a remote or local computation 
        bool remote;
#ifdef CUDA_INSTRUMENT_WRS
        int tpIndex;
        char phase;
#endif
}CudaRequest;

typedef struct _ParameterStruct{
  int numInteractions;
  int numBucketsPlusOne;
  int numMissedCores;
  int numEntities;// TODO: can be removed later on
  cudatype fperiod;
}ParameterStruct;

#ifdef CUDA_INSTRUMENT_WRS
void DataManagerTransferLocalTree(CudaMultipoleMoments *moments, int nMoments, CompactPartData *compactParts, int nCompactParts, int mype, char phase);
void DataManagerTransferRemoteChunk(CudaMultipoleMoments *moments, int nMoments, CompactPartData *compactParts, int nCompactParts, int mype, char phase);
void FreeDataManagerLocalTreeMemory(bool freemom, bool freepart, int pe, char phase);
void FreeDataManagerRemoteChunkMemory(int , void *, bool freemom, bool freepart, int pe, char phase);
void TransferParticleVarsBack(VariablePartData *hostBuffer, int size, void *cb, bool, bool, int pe, char phase);
#else
void DataManagerTransferLocalTree(CudaMultipoleMoments *moments, int nMoments, CompactPartData *compactParts, int nCompactParts, int mype);
void DataManagerTransferRemoteChunk(CudaMultipoleMoments *moments, int nMoments, CompactPartData *compactParts, int nCompactParts);
void FreeDataManagerLocalTreeMemory(bool freemom, bool freepart);
void FreeDataManagerRemoteChunkMemory(int , void *, bool freemom, bool freepart);
void TransferParticleVarsBack(VariablePartData *hostBuffer, int size, void *cb, bool, bool);
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
