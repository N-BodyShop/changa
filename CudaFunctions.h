#ifndef _CUDAFUNCTIONS_H_
#define _CUDAFUNCTIONS_H_

#ifdef CUDA
#include "hapi.h"
#include "HostCUDA.h"
void TreePieceCellListDataTransferBasic(CudaRequest *data, hapiWorkRequest *wr);
void TreePiecePartListDataTransferBasic(CudaRequest *data, hapiWorkRequest *wr);

#ifdef GPU_LOCAL_TREE_WALK
__global__ void gpuLocalTreeWalk(
    CompactPartData *particleCores,
    VariablePartData *particleVars,
    CudaMultipoleMoments* moments,
    int firstParticle,
    int lastParticle,
    int rootIdx,
    cudatype theta,
    cudatype thetaMono,
    int nReplicas,
    cudatype fperiod,
    cudatype fperiodY,
    cudatype fperiodZ);
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
