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
		int *ilmarks,
		int *bucketStarts,
		int *bucketSizes,
		cudatype fperiod
                );

__global__ void particleGravityComputation(
                CompactPartData *targetCores,
                VariablePartData *targetVars,
                CompactPartData *sourceCores,
                ILPart *ils,
                int *ilmarks,
		int *bucketStarts,
		int *bucketSizes,
		cudatype fperiod);

#endif

#endif
