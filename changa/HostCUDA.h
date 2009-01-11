#ifndef _HOST_CUDA_H_
#define _HOST_CUDA_H_


#include "cuda_typedef.h"
/* Boolean defines */
enum boolean {NO, YES};

/* defines for Hybrid API buffer indices */

#define POST_PREFETCH_MOMENTS        0
#define POST_PREFETCH_PARTICLE_CORES  1
#define PARTICLE_VARS      2

#define ILPART 0
#define ILCELL 0

#define NODE_BUCKET_MARKERS 1
#define NODE_BUCKET_START_MARKERS 2
#define NODE_BUCKET_SIZES 3

#define PART_BUCKET_MARKERS 1
#define PART_BUCKET_START_MARKERS 2
#define PART_BUCKET_SIZES 3

#define MISSED_MOMENTS 4
#define MISSED_PARTS 4

// node moments and particle cores
#define DM_TRANSFER_NBUFFERS 2

// interaction list
// list markers
// bucket starts
// bucket sizes
#define TP_GRAVITY_LOCAL_NBUFFERS 4

#define TP_NODE_GRAVITY_REMOTE_NBUFFERS 4
#define TP_PART_GRAVITY_REMOTE_NBUFFERS 5

#define TP_NODE_GRAVITY_REMOTE_RESUME_NBUFFERS 5
#define TP_PART_GRAVITY_REMOTE_RESUME_NBUFFERS 5

// tp_gravity_local uses arrays of particles and nodes already allocated on the gpu
// tp_gravity_remote uses arrays of nodes already on the gpu + particles from an array it supplies
// tp_gravity_remote_resume uses an array each of nodes and particles it supplies
enum kernels {DM_TRANSFER,
						  TP_GRAVITY_LOCAL,
						  TP_GRAVITY_REMOTE,
						  TP_GRAVITY_REMOTE_RESUME,
						  TP_PART_GRAVITY_LOCAL,
						  TP_PART_GRAVITY_REMOTE,
						  TP_PART_GRAVITY_REMOTE_RESUME};


typedef struct _CudaRequest{
        // can either be a ILCell* or an ILPart*
	void *list;
	int *bucketMarkers;
	int *bucketStarts;
	int *bucketSizes;
	int numInteractions;
	int numBucketsPlusOne;
        void *tp;

	// these buckets were finished in this work request
	int *affectedBuckets;
	// was the last bucket only partially finished?
	bool lastBucketComplete;
        void *cb;
        void *state;
}CudaRequest;

/*
typedef struct _PartListData{
	ILPart *partList;
	int *partListBucketMarkers;
	int *bucketStarts;
	int *bucketSizes;
	int numBucketsPlusOne;
	int numInteractions;

	// these buckets were finished in this work request
	int *affectedBuckets;
	// was the last bucket completely finished?
	bool lastBucketComplete;
        void *cb;
        void *tp; // tp that issued the work request
        void *state;
}PartListData;
*/

// these functions must follow C linkage
// FIXME - also, not defined NVCC
#if defined __cplusplus
extern "C" {
#endif

void DataManagerTransferMoments(CudaMultipoleMoments *moments, int nMoments);
void DataManagerTransferParticles(CompactPartData *compactParts, int nCompactParts);

void TreePieceCellListDataTransferLocal(CudaRequest *data);
void TreePieceCellListDataTransferRemote(CudaRequest *data);
void TreePieceCellListDataTransferRemoteResume(CudaRequest *data, CudaMultipoleMoments *missedMoments, int numMissedMoments);


void TreePiecePartListDataTransferLocal(CudaRequest *data);
// the 'missedParticles' here are actually prefetched particles
void TreePiecePartListDataTransferRemote(CudaRequest *data);
//void TreePiecePartListDataTransferRemote(CudaRequest *data, CompactPartData *missedParticles, int numMissedParticles);
void TreePiecePartListDataTransferRemoteResume(CudaRequest *data, CompactPartData *missedParticles, int numMissedParticles);


#ifdef __cplusplus
}
#endif

#endif
