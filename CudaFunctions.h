#ifndef _CUDAFUNCTIONS_H_
#define _CUDAFUNCTIONS_H_

#ifdef CUDA
#include "hapi.h"
#include "HostCUDA.h"
void DataTransferBasic(CudaRequest *data, CudaDevPtr *ptr, const char* functionTag);
void DataTransferBasicCleanup(CudaDevPtr *ptr, const char* functionTag);

#ifdef GPU_LOCAL_TREE_WALK
__global__ void gpuLocalTreeWalk(
    CudaMultipoleMoments* moments,
    CompactPartData *particleCores,
    VariablePartData *particleVars,
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

__global__ void ZeroVars(VariablePartData *particleVars, int nVars);

__global__ void particleGravityComputation(
                CompactPartData *targetCores,
                VariablePartData *targetVars,
                CompactPartData *sourceCores,
                ILCell *ils,
                int *ilmarks,
		int *bucketStarts,
		int *bucketSizes,
		cudatype fperiod);
#endif

#endif
