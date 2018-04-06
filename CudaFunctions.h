#ifndef _CUDAFUNCTIONS_H_
#define _CUDAFUNCTIONS_H_

#ifdef CUDA
#include "wr.h"
#include "HostCUDA.h"
void TreePieceCellListDataTransferBasic(CudaRequest *data, workRequest *wr);
void TreePiecePartListDataTransferBasic(CudaRequest *data, workRequest *wr);
void TreePieceLocalWalkDataTransferBasic(CudaRequest *data, workRequest *wr);

#ifdef GPU_LOCAL_TREE_WALK
__global__ void gpuLocalTreeWalk(
    CompactPartData *particleCores,
    VariablePartData *particleVars,
    CudaMultipoleMoments* moments,
    int *ilmarks,
    int totalNumOfParticles,
    cudatype theta,
    cudatype thetaMono);
#endif //GPU_LOCAL_TREE_WALK

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

#ifdef CUDA_2D_TB_KERNEL
__global__ void particleGravityComputation(
                CompactPartData *targetCores,
                VariablePartData *targetVars,
                CompactPartData *sourceCores,
                ILCell *ils,
                int *ilmarks,
		int *bucketStarts,
		int *bucketSizes,
		cudatype fperiod);
#else
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

#endif
