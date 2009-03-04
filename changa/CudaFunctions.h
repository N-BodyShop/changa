#ifndef _CUDAFUNCTIONS_H_
#define _CUDAFUNCTIONS_H_

#ifdef CUDA
#include "wr.h"
#include "HostCUDA.h"
void TreePieceCellListDataTransferBasic(CudaRequest *data, workRequest *wr);
void TreePiecePartListDataTransferBasic(CudaRequest *data, workRequest *wr);
__global__ void nodeGravityComputation(
		CompactPartData *particleCores,
		VariablePartData *particleVars,
		CudaMultipoleMoments *moments,
		ILCell *ils,
		int numInteractions,
		int *ilmarks,
		int *bucketStarts,
		int *bucketSizes,
		int numBucketsPlusOne, cudatype fperiod, int type,
                int numNodes);

__global__ void particleGravityComputation(
                CompactPartData *targetCores,
                VariablePartData *targetVars,
                CompactPartData *sourceCores,
                ILPart *ils,
		int numInteractions,
                int *ilmarks,
		int *bucketStarts,
		int *bucketSizes,
		int numBucketsPlusOne, cudatype fperiod, int type);

#endif

#endif
