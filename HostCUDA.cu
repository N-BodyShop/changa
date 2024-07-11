#ifdef _WIN32
#define NOMINMAX
#endif

#ifdef HAPI_MEMPOOL
#define GPU_MEMPOOL
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
// #include <cutil.h>
#include <assert.h>

#include "CudaFunctions.h"
#include "CUDAMoments.cu"
#include "HostCUDA.h"
#include "EwaldCUDA.h"

#include "hapi.h"
#include "cuda_typedef.h"
#include "cuda/intrinsics/voting.hu"
#include "cuda/intrinsics/shfl.hu"

#ifdef GPU_LOCAL_TREE_WALK
#include "codes.h"
#endif //GPU_LOCAL_TREE_WALK

#ifdef COOLING_MOLECULARH
#include "cooling_metal_H2.h"
#include "stiff.h"
#include "physconst.h"
#endif

#ifdef HAPI_TRACE
#  define HAPI_TRACE_BEGIN()   double trace_start_time = CmiWallTimer()
#  define HAPI_TRACE_END(ID)   traceUserBracketEvent(ID, trace_start_time, CmiWallTimer())
#else
#  define HAPI_TRACE_BEGIN() /* */
#  define HAPI_TRACE_END(ID) /* */
#endif

#define cudaChk(code) cudaErrorDie(code, #code, __FILE__, __LINE__)
inline void cudaErrorDie(cudaError_t retCode, const char* code,
                                              const char* file, int line) {
  if (retCode != cudaSuccess) {
    fprintf(stderr, "Fatal CUDA Error %s at %s:%d.\nReturn value %d from '%s'.",
        cudaGetErrorString(retCode), file, line, retCode, code);
    abort();
  }
}

#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
#include "converse.h"
#endif

__device__ __constant__ EwaldReadOnlyData cachedData[1];
__device__ __constant__ EwtData ewt[NEWH];  


//__constant__ constantData[88];
//
//
#ifdef HAPI_TRACE
extern "C" void traceUserBracketEvent(int e, double beginT, double endT);
extern "C" double CmiWallTimer();
#endif


void allocatePinnedHostMemory(void **ptr, size_t size){
  if(size <= 0){
    *((char **)ptr) = NULL;
#ifdef CUDA_PRINT_ERRORS
    printf("allocatePinnedHostMemory: 0 size!\n");
#endif
    assert(0);
    return;
  }
#ifdef HAPI_MEMPOOL
  hapiMallocHost(ptr, size, true);
#else
  hapiMallocHost(ptr, size, false);
#endif
#ifdef CUDA_PRINT_ERRORS
  printf("allocatePinnedHostMemory: %s size: %zu\n", cudaGetErrorString( cudaGetLastError() ), size);
#endif
}

void freePinnedHostMemory(void *ptr){
  if(ptr == NULL){
#ifdef CUDA_PRINT_ERRORS
    printf("freePinnedHostMemory: NULL ptr!\n");
#endif
    assert(0);
    return;
  }
#ifdef HAPI_MEMPOOL
  hapiFreeHost(ptr, true);
#else
  hapiFreeHost(ptr, false);
#endif
#ifdef CUDA_PRINT_ERRORS
  printf("freePinnedHostMemory: %s\n", cudaGetErrorString( cudaGetLastError() ));
#endif
}

/// @brief Transfer local moments, particle data and acceleration fields to GPU memory
/// @param moments Array of moments
/// @param sMoments Size of moments array
/// @param compactParts Array of particles
/// @param sCompactParts Size of particle array
/// @param varParts Zeroed-out particle acceleration fields
/// @param sVarParts Size of acceleration array
/// @param d_localMoments Uninitalized pointer to moments on GPU
/// @param d_compactParts Uninitalized pointer to particles on GPU
/// @param d_varParts Uninitalized pointer to accelerations on GPU
/// @param stream CUDA stream to handle the memory transfer
/// @param numParticles Total number of particle accelerations to initalize
void DataManagerTransferLocalTree(void *moments, size_t sMoments,
                                  void *compactParts, size_t sCompactParts,
                                  void *varParts, size_t sVarParts,
				  void **d_localMoments, void **d_compactParts, void **d_varParts,
				  cudaStream_t stream, int numParticles,
                                  void *callback) {

#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
  printf("(%d) DM LOCAL TREE moments %zu partcores %zu partvars %zu\n",
           CmiMyPe(),
           sMoments,
           sCompactParts,
           sVarParts
           );
#endif

  HAPI_TRACE_BEGIN();

  cudaChk(cudaMalloc(d_localMoments, sMoments));
  cudaChk(cudaMalloc(d_compactParts, sCompactParts));
  cudaChk(cudaMalloc(d_varParts, sVarParts));

  cudaChk(cudaMemcpyAsync(*d_localMoments, moments, sMoments, cudaMemcpyHostToDevice, stream));
  cudaChk(cudaMemcpyAsync(*d_compactParts, compactParts, sCompactParts, cudaMemcpyHostToDevice, stream));
  cudaChk(cudaMemcpyAsync(*d_varParts, varParts, sVarParts, cudaMemcpyHostToDevice, stream));

  ZeroVars<<<numParticles / THREADS_PER_BLOCK + 1, dim3(THREADS_PER_BLOCK), 0, stream>>>(
      (VariablePartData *) *d_varParts,
      numParticles);
  cudaChk(cudaPeekAtLastError());

  HAPI_TRACE_END(CUDA_XFER_LOCAL);

  hapiAddCallback(stream, callback);
}

/// @brief Transfer remote moments and particle data to GPU memory
/// @param moments Array of remote moments
/// @param sMoments Size of remote moments array
/// @param remoteParts Array of remote particles
/// @param sRemoteParts Size of remote particle array
/// @param d_remoteMoments Uninitalized pointer to remote moments on GPU
/// @param d_remoteParts Uninitalized pointer to remote particles on GPU
/// @param stream CUDA stream to handle the memory transfer
void DataManagerTransferRemoteChunk(void *moments, size_t sMoments,
                                    void *remoteParts, size_t sRemoteParts,
				    void **d_remoteMoments, void **d_remoteParts,
                                    cudaStream_t stream,
                                    void *callback) {

#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
  printf("(%d) DM REMOTE CHUNK moments %zu partcores %zu\n",
        CmiMyPe(),
        sMoments,
        sRemoteParts
        );
#endif

  HAPI_TRACE_BEGIN();

  cudaChk(cudaMalloc(d_remoteMoments, sMoments));
  cudaChk(cudaMalloc(d_remoteParts, sRemoteParts));
  cudaChk(cudaMemcpyAsync(*d_remoteMoments, moments, sMoments, cudaMemcpyHostToDevice, stream));
  cudaChk(cudaMemcpyAsync(*d_remoteParts, remoteParts, sRemoteParts, cudaMemcpyHostToDevice, stream));

  HAPI_TRACE_END(CUDA_XFER_REMOTE);

  hapiAddCallback(stream, callback);
}

/************** Gravity *****************/

/// @brief Initiate a local gravity calculation on the GPU, via an interaction
///        list calculation between nodes, or do a local tree walk
/// @param data CudaRequest object containing parameters for the calculation
void TreePieceCellListDataTransferLocal(CudaRequest *data){
  cudaStream_t stream = data->stream;
  CudaDevPtr devPtr;
  TreePieceDataTransferBasic(data, &devPtr);

#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
  printf("(%d) TRANSFER LOCAL CELL\n", CmiMyPe());
#endif

#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
  printf("TRANSFER LOCAL CELL KERNELSELECT buffers:\nlocal_particles: (0x%x)\nlocal_particle_vars: (0x%x)\nremote_moments: (0x%x)\nil_cell: (0x%x)\n", 
	data->d_localParts,
	data->d_localVars,
	data->d_localMoments,
	devPtr.d_list
      );
#endif

  HAPI_TRACE_BEGIN();
#ifndef CUDA_NO_KERNELS
#ifdef GPU_LOCAL_TREE_WALK
  gpuLocalTreeWalk<<<(data->lastParticle - data->firstParticle + 1)
                     / THREADS_PER_BLOCK + 1, dim3(THREADS_PER_BLOCK), 0, stream>>> (
    data->d_localMoments,
    data->d_localParts,
    data->d_localVars,
    data->firstParticle,
    data->lastParticle,
    data->rootIdx,
    data->theta,
    data->thetaMono,
    data->nReplicas,
    data->fperiod,
    data->fperiodY,
    data->fperiodZ
    );
 #else
    dim3 dimensions = THREADS_PER_BLOCK;
 #ifdef CUDA_2D_TB_KERNEL
    dimensions = dim3(NODES_PER_BLOCK, PARTS_PER_BLOCK);
 #endif
  nodeGravityComputation<<<dim3(data->numBucketsPlusOne-1), dimensions, 0, stream>>> (
    data->d_localParts,
    data->d_localVars,
    data->d_localMoments,
    (ILCell *)devPtr.d_list,
    devPtr.d_bucketMarkers,
    devPtr.d_bucketStarts,
    devPtr.d_bucketSizes,
    data->fperiod
    );
#endif
#endif
  TreePieceDataTransferBasicCleanup(&devPtr);
  cudaChk(cudaPeekAtLastError());
  HAPI_TRACE_END(CUDA_GRAV_LOCAL);

  hapiAddCallback(stream, data->cb);
}

/// @brief Initiate a remote gravity calculation on the GPU between tree nodes
/// @param data CudaRequest object containing parameters for the calculation
void TreePieceCellListDataTransferRemote(CudaRequest *data){
  cudaStream_t stream = data->stream;
  CudaDevPtr devPtr;
  TreePieceDataTransferBasic(data, &devPtr);

#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
  printf("(%d) TRANSFER REMOTE CELL\n", CmiMyPe());
#endif

#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
  printf("TRANSFER REMOTE CELL KERNELSELECT buffers:\nlocal_particles: (0x%x)\nlocal_particle_vars: (0x%x)\nremote_moments: (0x%x)\nil_cell: (0x%x)\n", 
	data->d_localParts,
	data->d_localVars,
	data->d_remoteMoments,
	devPtr.d_list
        );
#endif

  HAPI_TRACE_BEGIN();
#ifndef CUDA_NO_KERNELS
  dim3 dimensions = THREADS_PER_BLOCK;
#ifdef CUDA_2D_TB_KERNEL
  dimensions = dim3(NODES_PER_BLOCK, PARTS_PER_BLOCK);
#endif
  nodeGravityComputation<<<data->numBucketsPlusOne-1, dimensions, 0, stream>>> (
    data->d_localParts,
    data->d_localVars,
    data->d_remoteMoments,
    (ILCell *)devPtr.d_list, 
    devPtr.d_bucketMarkers,
    devPtr.d_bucketStarts,
    devPtr.d_bucketSizes,
    data->fperiod
    ); 
#endif
  TreePieceDataTransferBasicCleanup(&devPtr);
  cudaChk(cudaPeekAtLastError());
  HAPI_TRACE_END(CUDA_GRAV_REMOTE);

  hapiAddCallback(stream, data->cb);
}

/// @brief Initiate a remote resume gravity calculation on the GPU between tree nodes
/// @param data CudaRequest object containing parameters for the calculation
void TreePieceCellListDataTransferRemoteResume(CudaRequest *data){
  cudaStream_t stream = data->stream;
  CudaDevPtr devPtr;
  void *d_missedNodes;
  TreePieceDataTransferBasic(data, &devPtr);

#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
  printf("(%d) TRANSFER REMOTE RESUME CELL\n", CmiMyPe());
#endif

  cudaChk(cudaMalloc(&d_missedNodes, data->sMissed));
  cudaChk(cudaMemcpyAsync(d_missedNodes, data->missedNodes, data->sMissed, cudaMemcpyHostToDevice, stream));

#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
  printf("TRANSFER REMOTE RESUME CELL KERNELSELECT buffers:\nlocal_particles: (0x%x)\nlocal_particle_vars: (0x%x)\nmissed_moments (0x%x)\nil_cell (0x%x)\n", 
        data->d_localParts,
	data->d_localVars,
	d_missedNodes,
	devPtr.d_list
      );
#endif

  HAPI_TRACE_BEGIN();
#ifndef CUDA_NO_KERNELS
    dim3 dimensions = THREADS_PER_BLOCK;
#ifdef CUDA_2D_TB_KERNEL
    dimensions = dim3(NODES_PER_BLOCK, PARTS_PER_BLOCK);
#endif
    nodeGravityComputation<<<data->numBucketsPlusOne-1, dimensions, 0, stream>>> (
      data->d_localParts,
      data->d_localVars,
      (CudaMultipoleMoments *)d_missedNodes,
      (ILCell *)devPtr.d_list,
      devPtr.d_bucketMarkers,
      devPtr.d_bucketStarts,
      devPtr.d_bucketSizes,
      data->fperiod
      );
#endif
  TreePieceDataTransferBasicCleanup(&devPtr);
  cudaChk(cudaFree(d_missedNodes));
  cudaChk(cudaPeekAtLastError());
  HAPI_TRACE_END(CUDA_REMOTE_RESUME);

  hapiAddCallback(stream, data->cb);
}

/// @brief Initiate a small phase local gravity calculation on the GPU between particles
/// @param data CudaRequest object containing parameters for the calculation
void TreePiecePartListDataTransferLocalSmallPhase(CudaRequest *data, CompactPartData *particles, int len){
  cudaStream_t stream = data->stream;
  CudaDevPtr devPtr;
  TreePieceDataTransferBasic(data, &devPtr);

  size_t size = (len) * sizeof(CompactPartData);
  void* bufferHostBuffer;
  void* d_smallParts;

#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
  printf("TreePiecePartListDataTransferLocalSmallPhase KERNELSELECT buffers:\nlocal_particles: (0x%x)\nlocal_particle_vars: (0x%x)\nil_cell: (0x%x)\n",
      data->d_localParts,
      data->d_localVars,
      devPtr.d_list
      );
#endif

#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
  printf("(%d) TRANSFER LOCAL SMALL PHASE  %zu\n",
      CmiMyPe(),
      size
      );
#endif

  HAPI_TRACE_BEGIN();
  allocatePinnedHostMemory(&bufferHostBuffer, size);
#ifdef CUDA_PRINT_ERRORS
  printf("TPPartSmallPhase 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
  memcpy(bufferHostBuffer, particles, size);
  cudaChk(cudaMalloc(&d_smallParts, size));
  cudaChk(cudaMemcpyAsync(d_smallParts, bufferHostBuffer, size, cudaMemcpyHostToDevice, stream));

#ifndef CUDA_NO_KERNELS
#ifdef CUDA_2D_TB_KERNEL
  particleGravityComputation<<<data->numBucketsPlusOne-1, dim3(NODES_PER_BLOCK_PART, PARTS_PER_BLOCK_PART), 0, stream>>> (
    data->d_localParts,
    data->d_localVars,
    (CompactPartData *)d_smallParts,
    (ILCell *)devPtr.d_list,
    devPtr.d_bucketMarkers,
    devPtr.d_bucketStarts,
    devPtr.d_bucketSizes,
    data->fperiod
    );
#else
  particleGravityComputation<<<data->numBucketsPlusOne-1, THREADS_PER_BLOCK, 0, stream>>> (
    data->d_localParts,
    data->d_localVars,
    (CompactPartData *)d_smallParts,
    (ILCell *)devPtr.d_list,
    devPtr.d_bucketMarkers,
    devPtr.d_bucketStarts,
    devPtr.d_bucketSizes,
    data->fperiod
    );
#endif
#endif
  TreePieceDataTransferBasicCleanup(&devPtr);
  cudaChk(cudaPeekAtLastError());
  HAPI_TRACE_END(CUDA_PART_GRAV_LOCAL_SMALL);
  cudaChk(cudaFree(d_smallParts));
  hapiAddCallback(stream, data->cb);
}

/// @brief Initiate a local gravity calculation on the GPU between particles
/// @param data CudaRequest object containing parameters for the calculation
void TreePiecePartListDataTransferLocal(CudaRequest *data){
  cudaStream_t stream = data->stream;
  CudaDevPtr devPtr;
  TreePieceDataTransferBasic(data, &devPtr);

#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
  printf("(%d) TRANSFER LOCAL LARGEPHASE PART\n", CmiMyPe());
#endif

#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
  printf("TreePiecePartListDataTransferLocal buffers:\nlocal_particles: (0x%x)\nlocal_particle_vars: (0x%x)\nil_cell: (0x%x)\n",
        data->d_localParts,
        data->d_localVars,
        devPtr.d_list
        );
#endif

  HAPI_TRACE_BEGIN();
#ifndef CUDA_NO_KERNELS
#ifdef CUDA_2D_TB_KERNEL
  particleGravityComputation<<<data->numBucketsPlusOne-1, dim3(NODES_PER_BLOCK_PART, PARTS_PER_BLOCK_PART), 0, stream>>> (
    data->d_localParts,
    data->d_localVars,
    data->d_localParts,
    (ILCell *)devPtr.d_list,
    devPtr.d_bucketMarkers,
    devPtr.d_bucketStarts,
    devPtr.d_bucketSizes,
    data->fperiod
    );
#else
  particleGravityComputation<<<data->numBucketsPlusOne-1, THREADS_PER_BLOCK, 0, stream>>> (
    data->d_localParts,
    data->d_localVars,
    data->d_localParts,
    (ILPart *)devPtr.d_list,
    devPtr.d_bucketMarkers,
    devPtr.d_bucketStarts,
    devPtr.d_bucketSizes,
    data->fperiod
    );
#endif
#endif
  TreePieceDataTransferBasicCleanup(&devPtr);
  cudaChk(cudaPeekAtLastError());
  HAPI_TRACE_END(CUDA_PART_GRAV_LOCAL);

  hapiAddCallback(stream, data->cb);
}

/// @brief Initiate a remote gravity calculation on the GPU between particles
/// @param data CudaRequest object containing parameters for the calculation
void TreePiecePartListDataTransferRemote(CudaRequest *data){
  cudaStream_t stream = data->stream;
  CudaDevPtr devPtr;
  TreePieceDataTransferBasic(data, &devPtr);

#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
  printf("(%d) TRANSFER REMOTE PART\n", CmiMyPe());
#endif

#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
  printf("TreePiecePartListDataTransferRemote KERNELSELECT buffers:\nlocal_particles: (0x%x)\nlocal_particle_vars: (0x%x)\nil_cell: (0x%x) (0x%x)\n",
        data->d_localParts,
        data->d_localVars,
	data->list,
        devPtr.d_list
        );
#endif

  HAPI_TRACE_BEGIN();
#ifndef CUDA_NO_KERNELS
#ifdef CUDA_2D_TB_KERNEL
  particleGravityComputation<<<data->numBucketsPlusOne-1, dim3(NODES_PER_BLOCK_PART, PARTS_PER_BLOCK_PART), 0, stream>>> (
    data->d_localParts,
    data->d_localVars,
    data->d_remoteParts,
    (ILCell *)devPtr.d_list,
    devPtr.d_bucketMarkers,
    devPtr.d_bucketStarts,
    devPtr.d_bucketSizes,
    data->fperiod
    );
#else
  particleGravityComputation<<<data->numBucketsPlusOne-1, THREADS_PER_BLOCK, 0, stream>>> (
    data->d_localParts,
    data->d_localVars,
    data->d_remoteParts,
    (ILPart *)devPtr.d_list,
    devPtr.d_bucketMarkers,
    devPtr.d_bucketStarts,
    devPtr.d_bucketSizes,
    data->fperiod
    );
#endif
#endif
  TreePieceDataTransferBasicCleanup(&devPtr);
  cudaChk(cudaPeekAtLastError());
  HAPI_TRACE_END(CUDA_PART_GRAV_REMOTE);

  hapiAddCallback(stream, data->cb);
}

/// @brief Initiate a remote gravity calculation on the GPU between particles
/// @param data CudaRequest object containing parameters for the calculation
void TreePiecePartListDataTransferRemoteResume(CudaRequest *data){
  cudaStream_t stream = data->stream;
  CudaDevPtr devPtr;
  void* d_missedParts;
  TreePieceDataTransferBasic(data, &devPtr);

#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
  printf("(%d) TRANSFER REMOTE RESUME PART\n", CmiMyPe());
#endif

  cudaChk(cudaMalloc(&d_missedParts, data->sMissed));
  cudaChk(cudaMemcpyAsync(d_missedParts, data->missedParts, data->sMissed, cudaMemcpyHostToDevice, stream));

#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
  printf("TreePiecePartListDataTransferRemoteResume KERNELSELECT buffers:\nlocal_particles: (0x%x)\nlocal_particle_vars: (0x%x)\nmissed_parts (0x%x)\nil_cell: (0x%x) (0x%x)\n", 
        data->d_localParts,
        data->d_localVars,
        (CompactPartData *)d_missedParts,
	data->list,
        devPtr.d_list
        );
#endif

  HAPI_TRACE_BEGIN();
#ifndef CUDA_NO_KERNELS
#ifdef CUDA_2D_TB_KERNEL
  particleGravityComputation<<<data->numBucketsPlusOne-1, dim3(NODES_PER_BLOCK_PART, PARTS_PER_BLOCK_PART), 0, stream>>> (
    data->d_localParts,
    data->d_localVars,
    (CompactPartData *)d_missedParts,
    (ILCell *)devPtr.d_list,
    devPtr.d_bucketMarkers,
    devPtr.d_bucketStarts,
    devPtr.d_bucketSizes,
    data->fperiod
    );
#else
  particleGravityComputation<<<data->numBucketsPlusOne-1, THREADS_PER_BLOCK, 0, stream>>> (
    data->d_localParts,
    data->d_localVars,
    (CompactPartData *)d_missedParts,
    (ILPart *)devPtr.d_list,
    devPtr.d_bucketMarkers,
    devPtr.d_bucketStarts,
    devPtr.d_bucketSizes,
    data->fperiod
    );
#endif
#endif
  TreePieceDataTransferBasicCleanup(&devPtr);
  cudaChk(cudaFree(d_missedParts));
  cudaChk(cudaPeekAtLastError());
  HAPI_TRACE_END(CUDA_PART_GRAV_REMOTE);

  hapiAddCallback(stream, data->cb);
}

/// @brief Allocate space and copy bucket and interaction list data to
///         device memory
/// @param data CudaRequest object containing parameters for the calculation
/// @param ptr CudaDevPtr object that stores handles to device memory
void TreePieceDataTransferBasic(CudaRequest *data, CudaDevPtr *ptr){
  cudaStream_t stream = data->stream;

  int numBucketsPlusOne = data->numBucketsPlusOne;
  int numBuckets = numBucketsPlusOne-1;
  size_t listSize = (data->numInteractions) * sizeof(ILCell);
  size_t markerSize = (numBucketsPlusOne) * sizeof(int);
  size_t startSize = (numBuckets) * sizeof(int);

  cudaChk(cudaMalloc(&ptr->d_list, listSize));
  cudaChk(cudaMalloc(&ptr->d_bucketMarkers, markerSize));
  cudaChk(cudaMalloc(&ptr->d_bucketStarts, startSize));
  cudaChk(cudaMalloc(&ptr->d_bucketSizes, startSize));
  cudaChk(cudaMemcpyAsync(ptr->d_list, data->list, listSize, cudaMemcpyHostToDevice, stream));
  cudaChk(cudaMemcpyAsync(ptr->d_bucketMarkers, data->bucketMarkers, markerSize, cudaMemcpyHostToDevice, stream));
  cudaChk(cudaMemcpyAsync(ptr->d_bucketStarts, data->bucketStarts, startSize, cudaMemcpyHostToDevice, stream));
  cudaChk(cudaMemcpyAsync(ptr->d_bucketSizes, data->bucketSizes, startSize, cudaMemcpyHostToDevice, stream));

#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
    printf("(%d) TRANSFER BASIC %zu bucket_markers %zu bucket_starts %zu\n",
           CmiMyPe(),
           listSize,
           markerSize,
           startSize
           );
#endif
}

/// @brief Free device memory used for interaction list and bucket data
/// @param ptr CudaDevPtr object that stores handles to device memory
void TreePieceDataTransferBasicCleanup(CudaDevPtr *ptr){
  cudaChk(cudaFree(ptr->d_list));
  cudaChk(cudaFree(ptr->d_bucketMarkers));
  cudaChk(cudaFree(ptr->d_bucketStarts));
  cudaChk(cudaFree(ptr->d_bucketSizes));
}

/** @brief Transfer forces from the GPU back to the host. Also schedules
 *         the freeing of the device buffers used for the force calculation.
 *  @param hostBuffer Buffer to store results.
 *  @param size hostBuffer size.
 *  @param d_varParts Pointer to finalized accelerations on GPU
 *  @param stream CUDA stream to handle the memory transfer
 *  @param cb Callback when transfer is done.
 */
void TransferParticleVarsBack(VariablePartData *hostBuffer, size_t size, void *d_varParts,
                              cudaStream_t stream, void *cb){
  
  HAPI_TRACE_BEGIN();
  cudaChk(cudaMemcpyAsync(hostBuffer, d_varParts, size, cudaMemcpyDeviceToHost, stream));
  HAPI_TRACE_END(CUDA_XFER_BACK);
  hapiAddCallback(stream, cb);
}

/*
void DummyKernel(void *cb){
  hapiWorkRequest* dummy = hapiCreateWorkRequest();

  dummy->setExecParams(1, THREADS_PER_BLOCK);

  dummy->setDeviceToHostCallback(cb);
#ifdef HAPI_TRACE
  dummy->setTraceName("dummyRun");
#endif
  dummy->setRunKernel(run_kernel_DUMMY);
  hapiEnqueue(dummy);

}
*/

/*
 * Kernels
 */

/****
 * GPU Local tree walk (computation integrated)
****/
#ifdef GPU_LOCAL_TREE_WALK
__device__ __forceinline__ void
ldgTreeNode(CUDATreeNode &m, CudaMultipoleMoments *ptr) {
  m.radius        = __ldg(&(ptr->radius));
  m.soft          = __ldg(&(ptr->soft));
  m.totalMass     = __ldg(&(ptr->totalMass));
  m.cm.x          = __ldg(&(ptr->cm.x));
  m.cm.y          = __ldg(&(ptr->cm.y));
  m.cm.z          = __ldg(&(ptr->cm.z));
  m.bucketStart   = __ldg(&(ptr->bucketStart));
  m.bucketSize    = __ldg(&(ptr->bucketSize));
  m.particleCount = __ldg(&(ptr->particleCount));
  m.children[0]   = __ldg(&(ptr->children[0]));
  m.children[1]   = __ldg(&(ptr->children[1]));
  m.type          = __ldg(&(ptr->type));
};

__device__ __forceinline__ void
ldgBucketNode(CUDABucketNode &m, CompactPartData *ptr) {
  m.soft      = __ldg(&(ptr->soft));
  m.totalMass = __ldg(&(ptr->mass));
  m.cm.x      = __ldg(&(ptr->position.x));
  m.cm.y      = __ldg(&(ptr->position.y));
  m.cm.z      = __ldg(&(ptr->position.z));
};

__device__ __forceinline__ void
ldgBucketNode(CUDABucketNode &m, CudaMultipoleMoments *ptr) {
  m.radius    = __ldg(&(ptr->radius));
  m.soft      = __ldg(&(ptr->soft));
  m.totalMass = __ldg(&(ptr->totalMass));
  m.cm.x      = __ldg(&(ptr->cm.x));
  m.cm.y      = __ldg(&(ptr->cm.y));
  m.cm.z      = __ldg(&(ptr->cm.z));

  m.lesser_corner.x   = __ldg(&(ptr->lesser_corner.x));
  m.lesser_corner.y   = __ldg(&(ptr->lesser_corner.y));
  m.lesser_corner.z   = __ldg(&(ptr->lesser_corner.z));

  m.greater_corner.x  = __ldg(&(ptr->greater_corner.x));
  m.greater_corner.y  = __ldg(&(ptr->greater_corner.y));
  m.greater_corner.z  = __ldg(&(ptr->greater_corner.z));
};

__device__ __forceinline__ void
ldgParticle(CompactPartData &m, CompactPartData *ptr) {
  m.mass       = __ldg(&(ptr->mass));
  m.soft       = __ldg(&(ptr->soft));
  m.position.x = __ldg(&(ptr->position.x));
  m.position.y = __ldg(&(ptr->position.y));
  m.position.z = __ldg(&(ptr->position.z));
}

__device__ __forceinline__ void stackInit(int &sp, int* stk, int rootIdx) {
  sp = 0;
  stk[sp] = rootIdx;
}

__device__ __forceinline__ void stackPush(int &sp) {
  ++sp;
}

__device__ __forceinline__ void stackPop(int &sp) {
  --sp;
}

const int stackDepth = 64;

//__launch_bounds__(1024,1)
__global__ void gpuLocalTreeWalk(
  CudaMultipoleMoments *moments,
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
  cudatype fperiodZ) {

  CUDABucketNode  myNode;
  CompactPartData myParticle;

#if __CUDA_ARCH__ >= 700
  // Non-lockstepping code for Volta GPUs
  int sp;
  int stk[stackDepth];
  CUDATreeNode targetNode;

#define SP sp
#define STACK_TOP_INDEX stk[SP]
#define TARGET_NODE targetNode
#else
  // Default lockstepping code
  __shared__ int sp[WARPS_PER_BLOCK];
  __shared__ int stk[WARPS_PER_BLOCK][stackDepth];
  __shared__ CUDATreeNode targetNode[WARPS_PER_BLOCK];

#define SP sp[WARP_INDEX]
#define STACK_TOP_INDEX stk[WARP_INDEX][SP]
#define TARGET_NODE targetNode[WARP_INDEX]
#endif

  CudaVector3D acc = {0,0,0};
  cudatype pot = 0;
  cudatype idt2 = 0;
  CudaVector3D offset = {0,0,0};

  int targetIndex = -1;

  // variables for CUDA_momEvalFmomrcm
  CudaVector3D r;
  cudatype rsq = 0;
  cudatype twoh = 0;

  int flag = 1;
  int critical = stackDepth;
  int cond = 1;

  for(int pidx = blockIdx.x*blockDim.x + threadIdx.x + firstParticle;
      pidx <= lastParticle; pidx += gridDim.x*blockDim.x) {
    
    // initialize the variables belonging to current thread
    int nodePointer = particleCores[pidx].nodeId;
    ldgParticle(myParticle, &particleCores[pidx]);
    ldgBucketNode(myNode, &moments[nodePointer]);

    for(int x = -nReplicas; x <= nReplicas; x++) {
      for(int y = -nReplicas; y <= nReplicas; y++) {
        for(int z = -nReplicas; z <= nReplicas; z++) {
          // generate the offset for the periodic boundary conditions 
          offset.x = x*fperiod;
          offset.y = y*fperiodY;
          offset.z = z*fperiodZ;

          flag = 1;
          critical = stackDepth;
          cond = 1;
#if __CUDA_ARCH__ >= 700
          stackInit(SP, stk, rootIdx);
#else
          stackInit(SP, stk[WARP_INDEX], rootIdx);
#endif
          while(SP >= 0) {
            if (flag == 0 && critical >= SP) {
              flag = 1;
            }

            targetIndex = STACK_TOP_INDEX;
            stackPop(SP);

            if (flag) {
              ldgTreeNode(TARGET_NODE, &moments[targetIndex]);
              // Each warp increases its own TARGET_NODE in the shared memory,
              // so there is no actual data racing here.
              addCudaVector3D(TARGET_NODE.cm, offset, TARGET_NODE.cm);

              int open = CUDA_openCriterionNode(TARGET_NODE, myNode, -1, theta,
                                                thetaMono);
              int action = CUDA_OptAction(open, TARGET_NODE.type);

              critical = SP;
              cond = ((action == KEEP) && (open == CONTAIN || open == INTERSECT));

              if (action == COMPUTE) {
                if (CUDA_openSoftening(TARGET_NODE, myNode)) {
                  r.x = TARGET_NODE.cm.x - myParticle.position.x;
                  r.y = TARGET_NODE.cm.y - myParticle.position.y;
                  r.z = TARGET_NODE.cm.z - myParticle.position.z;

                  rsq = r.x*r.x + r.y*r.y + r.z*r.z;
                  twoh = TARGET_NODE.soft + myParticle.soft;
                  cudatype a, b;
                  if (rsq != 0) {
                    CUDA_SPLINE(rsq, twoh, a, b);
                    idt2 = fmax(idt2, (myParticle.mass + TARGET_NODE.totalMass) * b);

                    pot -= TARGET_NODE.totalMass * a;

                    acc.x += r.x*b*TARGET_NODE.totalMass;
                    acc.y += r.y*b*TARGET_NODE.totalMass;
                    acc.z += r.z*b*TARGET_NODE.totalMass;
                  }
                } else {
                  // compute with the node targetnode
                  r.x = myParticle.position.x - TARGET_NODE.cm.x;
                  r.y = myParticle.position.y - TARGET_NODE.cm.y;
                  r.z = myParticle.position.z - TARGET_NODE.cm.z;

                  rsq = r.x*r.x + r.y*r.y + r.z*r.z;
                  if (rsq != 0) {
                    cudatype dir = rsqrt(rsq);
#if defined (HEXADECAPOLE)
                    CUDA_momEvalFmomrcm(&moments[targetIndex], &r, dir, &acc, &pot);
                    idt2 = fmax(idt2, (myParticle.mass +
                                       moments[targetIndex].totalMass)*dir*dir*dir);
#else
                    cudatype a, b, c, d;
                    twoh = moments[targetIndex].soft + myParticle.soft;
                    CUDA_SPLINEQ(dir, rsq, twoh, a, b, c, d);

                    cudatype qirx = moments[targetIndex].xx*r.x +
                                    moments[targetIndex].xy*r.y +
                                    moments[targetIndex].xz*r.z;
                    cudatype qiry = moments[targetIndex].xy*r.x +
                                    moments[targetIndex].yy*r.y +
                                    moments[targetIndex].yz*r.z;
                    cudatype qirz = moments[targetIndex].xz*r.x +
                                    moments[targetIndex].yz*r.y +
                                    moments[targetIndex].zz*r.z;
                    cudatype qir = 0.5 * (qirx*r.x + qiry*r.y + qirz*r.z);
                    cudatype tr = 0.5 * (moments[targetIndex].xx +
                                         moments[targetIndex].yy +
                                         moments[targetIndex].zz);
                    cudatype qir3 = b*moments[targetIndex].totalMass + d*qir - c*tr;

                    pot -= moments[targetIndex].totalMass * a + c*qir - b*tr;
                    acc.x -= qir3*r.x - c*qirx;
                    acc.y -= qir3*r.y - c*qiry;
                    acc.z -= qir3*r.z - c*qirz;
                    idt2 = fmax(idt2, (myParticle.mass +
                                       moments[targetIndex].totalMass) * b);
#endif //HEXADECAPOLE
                  }
                }
              } else if (action == KEEP_LOCAL_BUCKET) {
                // compute with each particle contained by node targetnode
                int target_firstparticle = TARGET_NODE.bucketStart;
                int target_lastparticle = TARGET_NODE.bucketStart +
                                          TARGET_NODE.bucketSize;
                cudatype a, b;
                for (int i = target_firstparticle; i < target_lastparticle; ++i) {
                  CompactPartData targetParticle = particleCores[i];
                  addCudaVector3D(targetParticle.position, offset, targetParticle.position);

                  r.x = targetParticle.position.x - myParticle.position.x;
                  r.y = targetParticle.position.y - myParticle.position.y;
                  r.z = targetParticle.position.z - myParticle.position.z;

                  rsq = r.x*r.x + r.y*r.y + r.z*r.z;
                  twoh = targetParticle.soft + myParticle.soft;
                  if (rsq != 0) {
                    CUDA_SPLINE(rsq, twoh, a, b);
                    pot -= targetParticle.mass * a;

                    acc.x += r.x*b*targetParticle.mass;
                    acc.y += r.y*b*targetParticle.mass;
                    acc.z += r.z*b*targetParticle.mass;
                    idt2 = fmax(idt2, (myParticle.mass + targetParticle.mass) * b);
                  }
                }
              }

              if (!any(cond)) {
                continue;
              }

              if (!cond) {
                flag = 0;
              } else {
                if (TARGET_NODE.children[1] != -1) {
                  stackPush(SP);
                  STACK_TOP_INDEX = TARGET_NODE.children[1];
                }
                if (TARGET_NODE.children[0] != -1) {
                  stackPush(SP);
                  STACK_TOP_INDEX = TARGET_NODE.children[0];
                }
              }
            }
#if __CUDA_ARCH__ >= 700
            __syncwarp();
#endif
          }
        } // z replicas
      } // y replicas
    } // x replicas

    particleVars[pidx].a.x += acc.x;
    particleVars[pidx].a.y += acc.y;
    particleVars[pidx].a.z += acc.z;
    particleVars[pidx].potential += pot;
    particleVars[pidx].dtGrav = fmax(idt2,  particleVars[pidx].dtGrav);
  }
}
#endif //GPU_LOCAL_TREE_WALK

/**
 * @brief interaction between multipole moments and buckets of particles.
 * @param particleCores Read-only properties of particles.
 * @param particleVars Accumulators of accelerations etc. of particles.
 * @param moments Multipole moments from which to calculate forces.
 * @param ils Cells on the interaction list.  Each Cell has an index into
 *            moments.
 * @param ilmarks Indices into ils for each block.
 * @param bucketStarts Indices into particleCores and particleVars
 *                      for each block
 * @param bucketSizes Size of the bucket for each block
 * @param fPeriod Size of periodic boundary condition.
 */

// 2d thread blocks 
#ifdef CUDA_2D_TB_KERNEL
#define TRANSLATE(x,y) (y*NODES_PER_BLOCK+x)
#ifndef CUDA_2D_FLAT
__device__ __forceinline__ void ldg_moments(CudaMultipoleMoments &m, CudaMultipoleMoments *ptr)
{
  m.radius = __ldg(&(ptr->radius));
  m.soft   = __ldg(&(ptr->soft));
  m.totalMass   = __ldg(&(ptr->totalMass));
  m.cm.x   = __ldg(&(ptr->cm.x));
  m.cm.y   = __ldg(&(ptr->cm.y));
  m.cm.z   = __ldg(&(ptr->cm.z));
#ifdef HEXADECAPOLE
  m.xx   = __ldg(&(ptr->xx));
  m.xy   = __ldg(&(ptr->xy));
  m.xz   = __ldg(&(ptr->xz));
  m.yy   = __ldg(&(ptr->yy));
  m.yz   = __ldg(&(ptr->yz));

  m.xxx   = __ldg(&(ptr->xxx));
  m.xyy   = __ldg(&(ptr->xyy));
  m.xxy   = __ldg(&(ptr->xxy));
  m.yyy   = __ldg(&(ptr->yyy));
  m.xxz   = __ldg(&(ptr->xxz));        
  m.yyz   = __ldg(&(ptr->yyz));        
  m.xyz   = __ldg(&(ptr->xyz));

  m.xxxx   = __ldg(&(ptr->xxxx));
  m.xyyy   = __ldg(&(ptr->xyyy));
  m.xxxy   = __ldg(&(ptr->xxxy));
  m.yyyy   = __ldg(&(ptr->yyyy));
  m.xxxz   = __ldg(&(ptr->xxxz));        
  m.yyyz   = __ldg(&(ptr->yyyz));        
  m.xxyy   = __ldg(&(ptr->xxyy));        
  m.xxyz   = __ldg(&(ptr->xxyz));        
  m.xyyz   = __ldg(&(ptr->xyyz));  
#else
  m.xx   = __ldg(&(ptr->xx));
  m.xy   = __ldg(&(ptr->xy));
  m.xz   = __ldg(&(ptr->xz));
  m.yy   = __ldg(&(ptr->yy));
  m.yz   = __ldg(&(ptr->yz));
  m.zz   = __ldg(&(ptr->zz));
#endif  
}

// we want to limit register usage to be 72 (by observing nvcc output)
// since GK100 has 64K registers, max threads per SM = (64K/72)
// then rounding down to multiple of 128 gives 896 
__launch_bounds__(896,1)
__global__ void nodeGravityComputation(
		CompactPartData *particleCores,
		VariablePartData *particleVars,
		CudaMultipoleMoments* moments,
		ILCell* ils,
		int *ilmarks,
		int *bucketStarts,
		int *bucketSizes,
		cudatype fperiod){
  
  // __shared__ CudaVector3D acc[THREADS_PER_BLOCK];
  // __shared__ cudatype pot[THREADS_PER_BLOCK];
  // __shared__ cudatype idt2[THREADS_PER_BLOCK];
  CudaVector3D acc;
  cudatype pot;
  cudatype idt2;
  __shared__ CudaMultipoleMoments m[NODES_PER_BLOCK];
  __shared__ int offsetID[NODES_PER_BLOCK];
  __shared__ CompactPartData shared_particle_cores[PARTS_PER_BLOCK];

  int
    start = ilmarks[blockIdx.x],
    end = ilmarks[blockIdx.x+1],
    bucketStart = bucketStarts[blockIdx.x];
  int bucketSize = bucketSizes[blockIdx.x];

  /*
  __shared__ int start;
 __shared__ int end;
 __shared__ int bucketStart;
 __shared__ int bucketSize;
 */
  
/*
  if(threadIdx.x == 0 && threadIdx.y == 0){
    start = ilmarks[blockIdx.x];
    end = ilmarks[blockIdx.x+1];
    bucketStart = bucketStarts[blockIdx.x];
    bucketSize = bucketSizes[blockIdx.x];
  }
  __syncthreads();
  */

  char
    tidx = threadIdx.x,
    tidy = threadIdx.y;

  for(int ystart = 0; ystart < bucketSize; ystart += PARTS_PER_BLOCK){
  

    int my_particle_idx = ystart + tidy;
    if(tidx == 0 && my_particle_idx < bucketSize){
      shared_particle_cores[tidy] = particleCores[bucketStart+my_particle_idx];
    }
     
    // __syncthreads(); // wait for leader threads to finish using acc's, pot's of other threads
    // acc[TRANSLATE(tidx,tidy)].x = 0.0;
    // acc[TRANSLATE(tidx,tidy)].y = 0.0;
    // acc[TRANSLATE(tidx,tidy)].z = 0.0;
    // pot[TRANSLATE(tidx,tidy)] = 0.0;
    // idt2[TRANSLATE(tidx,tidy)] = 0.0;
    acc.x = 0, acc.y = 0, acc.z = 0;
    pot = 0;
    idt2 = 0;
    
    
    for(int xstart = start; xstart < end; xstart += NODES_PER_BLOCK){
      int my_cell_idx = xstart + tidx;
      ILCell ilc;

      __syncthreads(); // wait for all threads to finish using 
                       // previous iteration's nodes before reloading
      
      if(tidy == 0 && my_cell_idx < end){
        ilc = ils[my_cell_idx];
        ldg_moments(m[tidx], &moments[ilc.index]);
        // m[tidx] = moments[ilc.index];
        offsetID[tidx] = ilc.offsetID;
      }
      
      __syncthreads(); // wait for nodes to be loaded before using them
      
      if(my_particle_idx < bucketSize && my_cell_idx < end){ // INTERACT
        CudaVector3D r;

        r.x = shared_particle_cores[tidy].position.x -
          ((((offsetID[tidx] >> 22) & 0x7)-3)*fperiod + m[tidx].cm.x);
        r.y = shared_particle_cores[tidy].position.y -
          ((((offsetID[tidx] >> 25) & 0x7)-3)*fperiod + m[tidx].cm.y);
        r.z = shared_particle_cores[tidy].position.z -
          ((((offsetID[tidx] >> 28) & 0x7)-3)*fperiod + m[tidx].cm.z);

        cudatype rsq = r.x*r.x + r.y*r.y + r.z*r.z;

        if(rsq != 0){
          cudatype dir = rsqrt(rsq);

#if defined(HEXADECAPOLE)
          // CUDA_momEvalFmomrcm(&m[tidx], &r, dir, &acc[TRANSLATE(tidx, tidy)], &pot[TRANSLATE(tidx, tidy)]);
          // idt2[TRANSLATE(tidx, tidy)] = fmax(idt2[TRANSLATE(tidx, tidy)],
          //                                (shared_particle_cores[tidy].mass + m[tidx].totalMass)*dir*dir*dir);
          CUDA_momEvalFmomrcm(&m[tidx], &r, dir, &acc, &pot);
          idt2 = fmax(idt2,
                      (shared_particle_cores[tidy].mass + m[tidx].totalMass)*dir*dir*dir);
#else
          cudatype a, b, c, d;
          cudatype
            twoh = m[tidx].soft + shared_particle_cores[tidy].soft;

          // SPLINEQ(dir, rsq, twoh, a, b, c, d);
          // expansion of function below:
          cudatype u,dih;
          if (rsq < twoh*twoh) {
            dih = 2.0/twoh;
            u = dih/dir;
            if (u < 1.0) {
              a = dih*(7.0/5.0 - 2.0/3.0*u*u + 3.0/10.0*u*u*u*u
                  - 1.0/10.0*u*u*u*u*u);
              b = dih*dih*dih*(4.0/3.0 - 6.0/5.0*u*u + 1.0/2.0*u*u*u);
              c = dih*dih*dih*dih*dih*(12.0/5.0 - 3.0/2.0*u);
              d = 3.0/2.0*dih*dih*dih*dih*dih*dih*dir;
            }
            else {
              a = -1.0/15.0*dir + dih*(8.0/5.0 - 4.0/3.0*u*u + u*u*u
                  - 3.0/10.0*u*u*u*u + 1.0/30.0*u*u*u*u*u);
              b = -1.0/15.0*dir*dir*dir + dih*dih*dih*(8.0/3.0 - 3.0*u
                  + 6.0/5.0*u*u - 1.0/6.0*u*u*u);
              c = -1.0/5.0*dir*dir*dir*dir*dir + 3.0*dih*dih*dih*dih*dir
                + dih*dih*dih*dih*dih*(-12.0/5.0 + 1.0/2.0*u);
              d = -dir*dir*dir*dir*dir*dir*dir
                + 3.0*dih*dih*dih*dih*dir*dir*dir
                - 1.0/2.0*dih*dih*dih*dih*dih*dih*dir;
            }
          }
          else {
            a = dir;
            b = a*a*a;
            c = 3.0*b*a*a;
            d = 5.0*c*a*a;
          }

          cudatype
            qirx = m[tidx].xx*r.x + m[tidx].xy*r.y + m[tidx].xz*r.z,
            qiry = m[tidx].xy*r.x + m[tidx].yy*r.y + m[tidx].yz*r.z,
            qirz = m[tidx].xz*r.x + m[tidx].yz*r.y + m[tidx].zz*r.z,
            qir = 0.5*(qirx*r.x + qiry*r.y + qirz*r.z),
            tr = 0.5*(m[tidx].xx + m[tidx].yy + m[tidx].zz),
            qir3 = b*m[tidx].totalMass + d*qir - c*tr;

          pot -= m[tidx].totalMass * a + c*qir - b*tr;
          acc.x -= qir3*r.x - c*qirx;
          acc.y -= qir3*r.y - c*qiry;
          acc.z -= qir3*r.z - c*qirz;
          idt2 = fmax(idt2, (shared_particle_cores[tidy].mass + m[tidx].totalMass)*b);

#endif
        }// end if rsq != 0
      }// end INTERACT
    }// end for each NODE group

    // __syncthreads(); // wait for all threads to finish before results become available

    cudatype sumx, sumy, sumz, poten, idt2max;
    // accumulate forces, potential in global memory data structure
    if (my_particle_idx < bucketSize) {
      sumx = acc.x, sumy = acc.y, sumz = acc.z;
      poten = pot;
      idt2max = idt2;
      for (int offset = NODES_PER_BLOCK/2; offset > 0; offset /= 2) {
        sumx += shfl_down(sumx, offset, NODES_PER_BLOCK);
        sumy += shfl_down(sumy, offset, NODES_PER_BLOCK);
        sumz += shfl_down(sumz, offset, NODES_PER_BLOCK);
        poten += shfl_down(poten, offset, NODES_PER_BLOCK);
        idt2max = fmax(idt2max, shfl_down(idt2max, offset, NODES_PER_BLOCK));
      }
      // if(tidx == 0 && my_particle_idx < bucketSize){
      if (tidx == 0) {
        // sumx = sumy = sumz = 0;
        // poten = 0;
        // idt2max = 0.0;
        // for(int i = 0; i < NODES_PER_BLOCK; i++){
        //   // sumx += acc[TRANSLATE(i,tidy)].x;
        //   // sumy += acc[TRANSLATE(i,tidy)].y;
        //   // sumz += acc[TRANSLATE(i,tidy)].z;
        //   // poten += pot[TRANSLATE(i,tidy)];
        //   idt2max = fmax(idt2[TRANSLATE(i,tidy)], idt2max);
        // }
        particleVars[bucketStart+my_particle_idx].a.x += sumx;
        particleVars[bucketStart+my_particle_idx].a.y += sumy;
        particleVars[bucketStart+my_particle_idx].a.z += sumz;
        particleVars[bucketStart+my_particle_idx].potential += poten;
        particleVars[bucketStart+my_particle_idx].dtGrav = fmax(idt2max,  particleVars[bucketStart+my_particle_idx].dtGrav);
      }
    }

  }// end for each PARTICLE group
}

#else 
__global__ void nodeGravityComputation(
		CompactPartData *particleCores,
		VariablePartData *particleVars,
		CudaMultipoleMoments *moments,
		ILCell *ils,
		int *ilmarks,
		int *bucketStarts,
		int *bucketSizes,
		cudatype fperiod){
  
  __shared__ cudatype accx[THREADS_PER_BLOCK];
  __shared__ cudatype accy[THREADS_PER_BLOCK];
  __shared__ cudatype accz[THREADS_PER_BLOCK];
  __shared__ cudatype pot[THREADS_PER_BLOCK];
  __shared__ cudatype idt2[THREADS_PER_BLOCK];
  //__shared__ cudatype mr[NODES_PER_BLOCK];
  __shared__ cudatype ms[NODES_PER_BLOCK];
  __shared__ cudatype mt[NODES_PER_BLOCK];
  __shared__ cudatype mcmx[NODES_PER_BLOCK];
  __shared__ cudatype mcmy[NODES_PER_BLOCK];
  __shared__ cudatype mcmz[NODES_PER_BLOCK];
  __shared__ cudatype mxx[NODES_PER_BLOCK];
  __shared__ cudatype mxy[NODES_PER_BLOCK];
  __shared__ cudatype mxz[NODES_PER_BLOCK];
  __shared__ cudatype myy[NODES_PER_BLOCK];
  __shared__ cudatype myz[NODES_PER_BLOCK];
  __shared__ cudatype mzz[NODES_PER_BLOCK];
  __shared__ int offsetID[NODES_PER_BLOCK];
  __shared__ CompactPartData shared_particle_cores[PARTS_PER_BLOCK];

  int start = ilmarks[blockIdx.x];
  int end = ilmarks[blockIdx.x+1];
  int bucketStart = bucketStarts[blockIdx.x];
  int bucketSize = bucketSizes[blockIdx.x];

  /*
  __shared__ int start;
 __shared__ int end;
 __shared__ int bucketStart;
 __shared__ int bucketSize;
 */

  int tx, ty;

/*
  if(threadIdx.x == 0 && threadIdx.y == 0){
    start = ilmarks[blockIdx.x];
    end = ilmarks[blockIdx.x+1];
    bucketStart = bucketStarts[blockIdx.x];
    bucketSize = bucketSizes[blockIdx.x];
  }
  __syncthreads();
  */

  int xstart;
  int ystart;
  tx = threadIdx.x;
  ty = threadIdx.y;

  for(ystart = 0; ystart < bucketSize; ystart += PARTS_PER_BLOCK){
  

    int my_particle_idx = ystart + ty;
    if(tx == 0 && my_particle_idx < bucketSize){
      shared_particle_cores[ty] = particleCores[bucketStart+my_particle_idx];
    }
    
    __syncthreads(); // wait for leader threads to finish using acc's, pot's of other threads
    accx[TRANSLATE(tx,ty)] = 0.0;
    accy[TRANSLATE(tx,ty)] = 0.0;
    accz[TRANSLATE(tx,ty)] = 0.0;
    pot[TRANSLATE(tx,ty)] = 0.0;
    idt2[TRANSLATE(tx,ty)] = 0.0;
    
    
    for(xstart = start; xstart < end; xstart += NODES_PER_BLOCK){
      int my_cell_idx = xstart + tx;
      ILCell ilc;

      __syncthreads(); // wait for all threads to finish using 
                       // previous iteration's nodes before reloading
      
      if(ty == 0 && my_cell_idx < end){
        ilc = ils[my_cell_idx];
        //mr[tx] = moments[ilc.index].radius;
        ms[tx] = moments[ilc.index].soft;
        mt[tx] = moments[ilc.index].totalMass;
        mcmx[tx] = moments[ilc.index].cm.x;
        mcmy[tx] = moments[ilc.index].cm.y;
        mcmz[tx] = moments[ilc.index].cm.z;
        mxx[tx] = moments[ilc.index].xx;
        mxy[tx] = moments[ilc.index].xy;
        mxz[tx] = moments[ilc.index].xz;
        myy[tx] = moments[ilc.index].yy;
        myz[tx] = moments[ilc.index].yz;
        mzz[tx] = moments[ilc.index].zz;
        offsetID[tx] = ilc.offsetID;
      }
      
      __syncthreads(); // wait for nodes to be loaded before using them
      
      if(my_particle_idx < bucketSize && my_cell_idx < end){ // INTERACT
        CudaVector3D r;
        cudatype rsq;
        cudatype twoh, a, b, c, d;

        r.x = shared_particle_cores[ty].position.x -
          ((((offsetID[tx] >> 22) & 0x7)-3)*fperiod + mcmx[tx]);
        r.y = shared_particle_cores[ty].position.y -
          ((((offsetID[tx] >> 25) & 0x7)-3)*fperiod + mcmy[tx]);
        r.z = shared_particle_cores[ty].position.z -
          ((((offsetID[tx] >> 28) & 0x7)-3)*fperiod + mcmz[tx]);

        rsq = r.x*r.x + r.y*r.y + r.z*r.z;        
        twoh = ms[tx] + shared_particle_cores[ty].soft;
        if(rsq != 0){
          cudatype dir = 1.0/sqrt(rsq);
          // SPLINEQ(dir, rsq, twoh, a, b, c, d);
          // expansion of function below:
          cudatype u,dih;
          if (rsq < twoh*twoh) {
            dih = 2.0/twoh;
            u = dih/dir;
            if (u < 1.0) {
              a = dih*(7.0/5.0 - 2.0/3.0*u*u + 3.0/10.0*u*u*u*u
                  - 1.0/10.0*u*u*u*u*u);
              b = dih*dih*dih*(4.0/3.0 - 6.0/5.0*u*u + 1.0/2.0*u*u*u);
              c = dih*dih*dih*dih*dih*(12.0/5.0 - 3.0/2.0*u);
              d = 3.0/2.0*dih*dih*dih*dih*dih*dih*dir;
            }
            else {
              a = -1.0/15.0*dir + dih*(8.0/5.0 - 4.0/3.0*u*u + u*u*u
                  - 3.0/10.0*u*u*u*u + 1.0/30.0*u*u*u*u*u);
              b = -1.0/15.0*dir*dir*dir + dih*dih*dih*(8.0/3.0 - 3.0*u
                  + 6.0/5.0*u*u - 1.0/6.0*u*u*u);
              c = -1.0/5.0*dir*dir*dir*dir*dir + 3.0*dih*dih*dih*dih*dir
                + dih*dih*dih*dih*dih*(-12.0/5.0 + 1.0/2.0*u);
              d = -dir*dir*dir*dir*dir*dir*dir
                + 3.0*dih*dih*dih*dih*dir*dir*dir
                - 1.0/2.0*dih*dih*dih*dih*dih*dih*dir;
            }
          }
          else {
            a = dir;
            b = a*a*a;
            c = 3.0*b*a*a;
            d = 5.0*c*a*a;
          }

          cudatype qirx = mxx[tx]*r.x + mxy[tx]*r.y + mxz[tx]*r.z;
          cudatype qiry = mxy[tx]*r.x + myy[tx]*r.y + myz[tx]*r.z;
          cudatype qirz = mxz[tx]*r.x + myz[tx]*r.y + mzz[tx]*r.z;
          cudatype qir = 0.5*(qirx*r.x + qiry*r.y + qirz*r.z);
          cudatype tr = 0.5*(mxx[tx] + myy[tx] + mzz[tx]);
          cudatype qir3 = b*mt[tx] + d*qir - c*tr;

          pot[TRANSLATE(tx, ty)] -= mt[tx] * a + c*qir - b*tr;

          accx[TRANSLATE(tx, ty)] -= qir3*r.x - c*qirx;
          accy[TRANSLATE(tx, ty)] -= qir3*r.y - c*qiry;
          accz[TRANSLATE(tx, ty)] -= qir3*r.z - c*qirz;
          idt2[TRANSLATE(tx, ty)] = fmax(idt2[TRANSLATE(tx, ty)],
                                    (shared_particle_cores[ty].mass + mt[tx]) * b);
        }// end if rsq != 0
      }// end INTERACT
    }// end for each NODE group

    __syncthreads(); // wait for all threads to finish before results become available

    cudatype sumx, sumy, sumz, poten, idt2max;
    sumx = sumy = sumz = poten = idt2max = 0.0;
    // accumulate forces, potential in global memory data structure
    if(tx == 0 && my_particle_idx < bucketSize){
      for(int i = 0; i < NODES_PER_BLOCK; i++){
        sumx += accx[TRANSLATE(i,ty)];
        sumy += accy[TRANSLATE(i,ty)];
        sumz += accz[TRANSLATE(i,ty)];
        poten += pot[TRANSLATE(i,ty)];
        idt2max = fmax(idt2[TRANSLATE(i,ty)], idt2max);
      }
      particleVars[bucketStart+my_particle_idx].a.x += sumx;
      particleVars[bucketStart+my_particle_idx].a.y += sumy;
      particleVars[bucketStart+my_particle_idx].a.z += sumz;
      particleVars[bucketStart+my_particle_idx].potential += poten;
      particleVars[bucketStart+my_particle_idx].dtGrav = fmax(idt2max,  particleVars[bucketStart+my_particle_idx].dtGrav);
    }

  }// end for each PARTICLE group
}
#endif
#else
__global__ void nodeGravityComputation(
		CompactPartData *particleCores,
		VariablePartData *particleVars,
		CudaMultipoleMoments *moments,
		ILCell *ils,
		int *ilmarks,
		int *bucketStarts,
		int *bucketSizes,
		cudatype fperiod){

  // each thread has its own storage for these
  __shared__ CudaVector3D acc[THREADS_PER_BLOCK];
  __shared__ cudatype pot[THREADS_PER_BLOCK];
  __shared__ cudatype idt2[THREADS_PER_BLOCK];
  __shared__ CudaMultipoleMoments m[THREADS_PER_BLOCK];

  __shared__ CompactPartData shared_particle_core;


  // each block is given a bucket to compute
  // each thread in the block computes an interaction of a particle with a node
  // threads must iterate through the interaction lists and sync.
  // then, block leader (first rank in each block) reduces the forces and commits 
  // values to global memory.
  int bucket = blockIdx.x;
  int start = ilmarks[bucket];
  int end = ilmarks[bucket+1];
  int bucketSize = bucketSizes[bucket];
  int bucketStart = bucketStarts[bucket];
  int thread = threadIdx.x;

  CudaVector3D r;
  cudatype rsq;
  cudatype twoh, a, b, c, d;

  for(int particle = 0; particle < bucketSize; particle++){
    if(thread == 0){
      // load shared_particle_core
      shared_particle_core = particleCores[bucketStart+particle];
    }
    __syncthreads();

    acc[thread].x = 0;
    acc[thread].y = 0;
    acc[thread].z = 0;
    pot[thread] = 0;
    idt2[thread] = 0;

    for(int node = start+thread; node < end; node+=THREADS_PER_BLOCK){
      ILCell ilc = ils[node];
      m[thread] = moments[ilc.index];
      int offsetID = ilc.offsetID;

      r.x = shared_particle_core.position.x -
        ((((offsetID >> 22) & 0x7)-3)*fperiod + m[thread].cm.x);
      r.y = shared_particle_core.position.y -
        ((((offsetID >> 25) & 0x7)-3)*fperiod + m[thread].cm.y);
      r.z = shared_particle_core.position.z -
        ((((offsetID >> 28) & 0x7)-3)*fperiod + m[thread].cm.z);

      rsq = r.x*r.x + r.y*r.y + r.z*r.z;        
      twoh = m[thread].soft + shared_particle_core.soft;
      if(rsq != 0){
        cudatype dir = 1.0/sqrt(rsq);
        // SPLINEQ(dir, rsq, twoh, a, b, c, d);
        // expansion of function below:
        cudatype u,dih;
        if (rsq < twoh*twoh) {
          dih = 2.0/twoh;
          u = dih/dir;
          if (u < 1.0) {
            a = dih*(7.0/5.0 - 2.0/3.0*u*u + 3.0/10.0*u*u*u*u
                - 1.0/10.0*u*u*u*u*u);
            b = dih*dih*dih*(4.0/3.0 - 6.0/5.0*u*u + 1.0/2.0*u*u*u);
            c = dih*dih*dih*dih*dih*(12.0/5.0 - 3.0/2.0*u);
            d = 3.0/2.0*dih*dih*dih*dih*dih*dih*dir;
          }
          else {
            a = -1.0/15.0*dir + dih*(8.0/5.0 - 4.0/3.0*u*u + u*u*u
                - 3.0/10.0*u*u*u*u + 1.0/30.0*u*u*u*u*u);
            b = -1.0/15.0*dir*dir*dir + dih*dih*dih*(8.0/3.0 - 3.0*u
                + 6.0/5.0*u*u - 1.0/6.0*u*u*u);
            c = -1.0/5.0*dir*dir*dir*dir*dir + 3.0*dih*dih*dih*dih*dir
              + dih*dih*dih*dih*dih*(-12.0/5.0 + 1.0/2.0*u);
            d = -dir*dir*dir*dir*dir*dir*dir
              + 3.0*dih*dih*dih*dih*dir*dir*dir
              - 1.0/2.0*dih*dih*dih*dih*dih*dih*dir;
          }
        }
        else {
          a = dir;
          b = a*a*a;
          c = 3.0*b*a*a;
          d = 5.0*c*a*a;
        }

        cudatype qirx = m[thread].xx*r.x + m[thread].xy*r.y + m[thread].xz*r.z;
        cudatype qiry = m[thread].xy*r.x + m[thread].yy*r.y + m[thread].yz*r.z;
        cudatype qirz = m[thread].xz*r.x + m[thread].yz*r.y + m[thread].zz*r.z;
        cudatype qir = 0.5*(qirx*r.x + qiry*r.y + qirz*r.z);
        cudatype tr = 0.5*(m[thread].xx + m[thread].yy + m[thread].zz);
        cudatype qir3 = b*m[thread].totalMass + d*qir - c*tr;

        pot[thread] -= m[thread].totalMass * a + c*qir - b*tr;

        acc[thread].x -= qir3*r.x - c*qirx;
        acc[thread].y -= qir3*r.y - c*qiry;
        acc[thread].z -= qir3*r.z - c*qirz;
        idt2[thread] = fmax(idt2[thread], (shared_particle_core.mass + m[thread].totalMass)*b);
      }// end if rsq != 0
    }// for each node in list
    __syncthreads();
    // at this point, the total force on particle is distributed among
    // all active threads;
    // reduce.
    // TODO: make this a parallel reduction

    cudatype sumx, sumy, sumz, poten, idt2max;
    sumx = sumy = sumz = poten = idt2max = 0.0;
    if(thread == 0){
      for(int i = 0; i < THREADS_PER_BLOCK; i++){
        sumx += acc[i].x;
        sumy += acc[i].y;
        sumz += acc[i].z;
        poten += pot[i];
        idt2max = fmax(idt2[i], idt2max);
      }
      particleVars[bucketStart+particle].a.x += sumx;
      particleVars[bucketStart+particle].a.y += sumy;
      particleVars[bucketStart+particle].a.z += sumz;
      particleVars[bucketStart+particle].potential += poten;
      particleVars[bucketStart+particle].dtGrav = fmax(idt2max,  particleVars[bucketStart+particle].dtGrav);
    }
  }// for each particle in bucket

      
}
#endif

/**
 * @brief interaction between source particles and buckets of particles.
 * @param particleCores Read-only properties of target particles.
 * @param particleVars Accumulators of accelerations etc. of target particles.
 * @param sourceCores Properties of source particles.
 * @param ils array of "cells": index into sourceCores and offset for each particle.
 * @param ilmarks Indices into ils for each block.
 * @param bucketStarts Indices into particleCores and particleVars
 *                      for each block
 * @param bucketSizes Size of the bucket for each block
 * @param fPeriod Size of periodic boundary condition.
 */
__device__ __forceinline__ void ldg_cPartData(CompactPartData &m, CompactPartData *ptr)
{
  m.mass         = __ldg(&(ptr->mass));
  m.soft         = __ldg(&(ptr->soft));
  m.position.x   = __ldg(&(ptr->position.x));
  m.position.y   = __ldg(&(ptr->position.y));
  m.position.z   = __ldg(&(ptr->position.z));
}

#ifdef CUDA_2D_TB_KERNEL
#define TRANSLATE_PART(x,y) (y*NODES_PER_BLOCK_PART+x)
//__launch_bounds__(maxThreadsPerBlock, minBlocksPerMultiprocessor)
//maxThreadsPerBlock needs to be a multiple of 128, this can be used 
//to limits the number of  register per thread. More threads less 
//registers  
__launch_bounds__(1408, 1)
__global__ void particleGravityComputation(
		CompactPartData *targetCores,
		VariablePartData *targetVars,
		CompactPartData* sourceCores,
		ILCell *ils,
		int *ilmarks,
		int *bucketStarts,
		int *bucketSizes,
		cudatype fperiod){
  
  //__shared__ CudaVector3D acc[THREADS_PER_BLOCK_PART];
  // __shared__ cudatype pot[THREADS_PER_BLOCK_PART];
  // __shared__ cudatype idt2[THREADS_PER_BLOCK_PART];
  CudaVector3D acc;
  cudatype pot;
  cudatype idt2;
  __shared__ CompactPartData m[NODES_PER_BLOCK_PART];
  __shared__ int offsetID[NODES_PER_BLOCK_PART];
  __shared__ CompactPartData shared_particle_cores[PARTS_PER_BLOCK_PART];

  int start = ilmarks[blockIdx.x];
  int end = ilmarks[blockIdx.x+1];
  int bucketStart = bucketStarts[blockIdx.x];
  int bucketSize = bucketSizes[blockIdx.x];

  /*
  __shared__ int start;
 __shared__ int end;
 __shared__ int bucketStart;
 __shared__ int bucketSize;
 */

  int tx, ty;

/*
  if(threadIdx.x == 0 && threadIdx.y == 0){
    start = ilmarks[blockIdx.x];
    end = ilmarks[blockIdx.x+1];
    bucketStart = bucketStarts[blockIdx.x];
    bucketSize = bucketSizes[blockIdx.x];
  }
  __syncthreads();
  */

  int xstart;
  int ystart;
  tx = threadIdx.x;
  ty = threadIdx.y;

  for(ystart = 0; ystart < bucketSize; ystart += PARTS_PER_BLOCK_PART){
  

    int my_particle_idx = ystart + ty;
    if(tx == 0 && my_particle_idx < bucketSize){
      shared_particle_cores[ty] = targetCores[bucketStart+my_particle_idx];
    }
    
    //__syncthreads(); // wait for leader threads to finish using acc's, pot's of other threads
    //acc[TRANSLATE_PART(tx,ty)].x = 0.0;
    //acc[TRANSLATE_PART(tx,ty)].y = 0.0;
    //acc[TRANSLATE_PART(tx,ty)].z = 0.0;
    //pot[TRANSLATE_PART(tx,ty)] = 0.0;
    //idt2[TRANSLATE_PART(tx,ty)] = 0.0;
    acc.x = 0, acc.y = 0, acc.z = 0;
    pot = 0;
    idt2 = 0;
    
    
    for(xstart = start; xstart < end; xstart += NODES_PER_BLOCK_PART){
      int my_cell_idx = xstart + tx;
      ILCell ilc;

      __syncthreads(); // wait for all threads to finish using 
                       // previous iteration's nodes before reloading
      
      if(ty == 0 && my_cell_idx < end){
        ilc = ils[my_cell_idx];
        ldg_cPartData(m[tx], &sourceCores[ilc.index]);
        //m[tx] = sourceCores[ilc.index];
        offsetID[tx] = ilc.offsetID;
      }
      
      __syncthreads(); // wait for nodes to be loaded before using them
      
      if(my_particle_idx < bucketSize && my_cell_idx < end){ // INTERACT
        CudaVector3D r;
        cudatype rsq;
        cudatype twoh, a, b;

        r.x = (((offsetID[tx] >> 22) & 0x7)-3)*fperiod 
              + m[tx].position.x
              - shared_particle_cores[ty].position.x;
        r.y = (((offsetID[tx] >> 25) & 0x7)-3)*fperiod 
              + m[tx].position.y
              - shared_particle_cores[ty].position.y;
        r.z = (((offsetID[tx] >> 28) & 0x7)-3)*fperiod 
              + m[tx].position.z
              - shared_particle_cores[ty].position.z;

        rsq = r.x*r.x + r.y*r.y + r.z*r.z;        
        twoh = m[tx].soft + shared_particle_cores[ty].soft;
        if(rsq != 0){
          cudatype r1, u, dih, dir;
          r1 = sqrt(rsq);
          if (r1 < (twoh)) {
            dih = 2.0/(twoh);
            u = r1*dih;
            if (u < 1.0) {
              a = dih*(7.0/5.0 - 2.0/3.0*u*u + 3.0/10.0*u*u*u*u
                  - 1.0/10.0*u*u*u*u*u);
              b = dih*dih*dih*(4.0/3.0 - 6.0/5.0*u*u
                  + 1.0/2.0*u*u*u);
            }
            else {
              dir = 1.0/r1;
              a = -1.0/15.0*dir + dih*(8.0/5.0 - 4.0/3.0*u*u +
                  u*u*u - 3.0/10.0*u*u*u*u + 1.0/30.0*u*u*u*u*u);
              b = -1.0/15.0*dir*dir*dir +
                dih*dih*dih*(8.0/3.0 - 3.0*u +
                    6.0/5.0*u*u - 1.0/6.0*u*u*u);
            }
          }
          else {
            a = 1.0/r1;
            b = a*a*a;
          }

          //pot[TRANSLATE_PART(tx, ty)] -= m[tx].mass * a;	

          //acc[TRANSLATE_PART(tx, ty)].x += r.x*b*m[tx].mass;
          //acc[TRANSLATE_PART(tx, ty)].y += r.y*b*m[tx].mass;
          //acc[TRANSLATE_PART(tx, ty)].z += r.z*b*m[tx].mass;
          //idt2[TRANSLATE_PART(tx, ty)] = fmax(idt2[TRANSLATE_PART(tx, ty)],
          //                               (shared_particle_cores[ty].mass + m[tx].mass) * b);

	  pot -= m[tx].mass * a;

          acc.x += r.x*b*m[tx].mass;
          acc.y += r.y*b*m[tx].mass;
          acc.z += r.z*b*m[tx].mass;
          idt2 = fmax(idt2, (shared_particle_cores[ty].mass + m[tx].mass) * b);

        }// end if rsq != 0
      }// end INTERACT
    }// end for each NODE group

    //__syncthreads(); // wait for all threads to finish before results become available

    cudatype sumx, sumy, sumz, poten, idt2max;
    sumx = sumy = sumz = poten = idt2max = 0.0;
    // accumulate forces, potential in global memory data structure
    //if(tx == 0 && my_particle_idx < bucketSize){
    //  for(int i = 0; i < NODES_PER_BLOCK_PART; i++){
    //    sumx += acc[TRANSLATE_PART(i,ty)].x;
    //    sumy += acc[TRANSLATE_PART(i,ty)].y;
    //    sumz += acc[TRANSLATE_PART(i,ty)].z;
    //    poten += pot[TRANSLATE_PART(i,ty)];
    //    idt2max = fmax(idt2[TRANSLATE_PART(i,ty)], idt2max);
    //  }

    // accumulate forces, potential in global memory data structure
    if(my_particle_idx < bucketSize){
      sumx = acc.x;
      sumy = acc.y; 
      sumz = acc.z;
      poten = pot;
      idt2max = idt2;
      for(int offset = NODES_PER_BLOCK/2; offset > 0; offset /= 2){
        sumx += shfl_down(sumx, offset, NODES_PER_BLOCK_PART);
        sumy += shfl_down(sumy, offset, NODES_PER_BLOCK_PART);
        sumz += shfl_down(sumz, offset, NODES_PER_BLOCK_PART);
        poten += shfl_down(poten, offset, NODES_PER_BLOCK_PART);
        idt2max = fmax(idt2max, shfl_down(idt2max, offset, NODES_PER_BLOCK_PART));
      }

      if(tx == 0){
      	targetVars[bucketStart+my_particle_idx].a.x += sumx;
      	targetVars[bucketStart+my_particle_idx].a.y += sumy;
     	targetVars[bucketStart+my_particle_idx].a.z += sumz;
      	targetVars[bucketStart+my_particle_idx].potential += poten;
      	targetVars[bucketStart+my_particle_idx].dtGrav = fmax(idt2max,  targetVars[bucketStart+my_particle_idx].dtGrav);
      }
    }

  }// end for each PARTICLE group
}
#else
__launch_bounds__(896, 1)
__global__ void particleGravityComputation(
                                   CompactPartData *targetCores,
                                   VariablePartData *targetVars,
                                   CompactPartData *sourceCores,
                                   ILPart *ils,
                                   int *ilmarks,
		                   int *bucketStarts,
		                   int *bucketSizes,
		                   cudatype fperiod){

  // each thread has its own storage for these
  __shared__ CudaVector3D acc[THREADS_PER_BLOCK];
  __shared__ cudatype pot[THREADS_PER_BLOCK];
  __shared__ cudatype idt2[THREADS_PER_BLOCK];
  __shared__ CompactPartData source_cores[THREADS_PER_BLOCK];

  __shared__ CompactPartData shared_target_core;


  // each block is given a bucket to compute
  // each thread in the block computes an interaction of a particle with a node
  // threads must iterate through the interaction lists and sync.
  // then, block leader (first rank in each block) reduces the forces and commits 
  // values to global memory.
  int bucket = blockIdx.x;
  int start = ilmarks[bucket];
  int end = ilmarks[bucket+1];
  int bucketSize = bucketSizes[bucket];
  int bucketStart = bucketStarts[bucket];
  int thread = threadIdx.x;

  CudaVector3D r;
  cudatype rsq;
  cudatype twoh, a, b;

  for(int target = 0; target < bucketSize; target++){
    if(thread == 0){
      shared_target_core = targetCores[bucketStart+target];
    }
    __syncthreads();

    acc[thread].x = 0;
    acc[thread].y = 0;
    acc[thread].z = 0;
    pot[thread] = 0;
    idt2[thread] = 0;

    for(int source = start+thread; source < end; source += THREADS_PER_BLOCK){
      ILPart ilp = ils[source]; 
      int oid = ilp.off;
      int num = ilp.num;
      int ilpindex = ilp.index;

      for(int particle = 0; particle < num; particle++){
        source_cores[thread] = sourceCores[ilpindex+particle];

        r.x = (((oid >> 22) & 0x7)-3)*fperiod +
          source_cores[thread].position.x -
          shared_target_core.position.x;

        r.y = (((oid >> 25) & 0x7)-3)*fperiod +
          source_cores[thread].position.y -
          shared_target_core.position.y;

        r.z = (((oid >> 28) & 0x7)-3)*fperiod +
          source_cores[thread].position.z -
          shared_target_core.position.z;

        rsq = r.x*r.x + r.y*r.y + r.z*r.z;
        twoh = source_cores[thread].soft + shared_target_core.soft;
        if(rsq != 0){
          cudatype r1, u,dih,dir;
          r1 = sqrt(rsq);
          if (r1 < (twoh)) {
            dih = 2.0/(twoh);
            u = r1*dih;
            if (u < 1.0) {
              a = dih*(7.0/5.0 - 2.0/3.0*u*u + 3.0/10.0*u*u*u*u
                  - 1.0/10.0*u*u*u*u*u);
              b = dih*dih*dih*(4.0/3.0 - 6.0/5.0*u*u
                  + 1.0/2.0*u*u*u);
            }
            else {
              dir = 1.0/r1;
              a = -1.0/15.0*dir + dih*(8.0/5.0 - 4.0/3.0*u*u +
                  u*u*u - 3.0/10.0*u*u*u*u + 1.0/30.0*u*u*u*u*u);
              b = -1.0/15.0*dir*dir*dir +
                dih*dih*dih*(8.0/3.0 - 3.0*u +
                    6.0/5.0*u*u - 1.0/6.0*u*u*u);
            }
          }
          else {
            a = 1.0/r1;
            b = a*a*a;
          }

          pot[thread] -= source_cores[thread].mass * a;

          acc[thread].x += r.x*b*source_cores[thread].mass;
          acc[thread].y += r.y*b*source_cores[thread].mass;
          acc[thread].z += r.z*b*source_cores[thread].mass;
          idt2[thread] = fmax(idt2[thread], (shared_target_core.mass + source_cores[thread].mass) * b);
        }// if rsq != 0
      }// for each particle in source bucket
    }// for each source bucket 

    __syncthreads();
  
    cudatype sumx, sumy, sumz, poten, idt2max;
    sumx = sumy = sumz = poten = idt2max = 0.0;
    if(thread == 0){
      for(int i = 0; i < THREADS_PER_BLOCK; i++){
        sumx += acc[i].x;
        sumy += acc[i].y;
        sumz += acc[i].z;
        poten += pot[i];
        idt2max = fmax(idt2[i], idt2max);
      }
      targetVars[bucketStart+target].a.x += sumx;
      targetVars[bucketStart+target].a.y += sumy;
      targetVars[bucketStart+target].a.z += sumz;
      targetVars[bucketStart+target].potential += poten;
      targetVars[bucketStart+target].dtGrav = fmax(idt2max,  targetVars[bucketStart+target].dtGrav);
    }

  }// for each target part
}
#endif

__global__ void EwaldKernel(CompactPartData *particleCores, VariablePartData *particleVars, int *markers, int largephase, int First, int Last);

extern unsigned int timerHandle; 

void EwaldHostMemorySetup(EwaldData *h_idata, int nParticles, int nEwhLoop, int largephase) {
  if(largephase)
    allocatePinnedHostMemory((void **)&(h_idata->EwaldMarkers), nParticles*sizeof(int));
  else
    h_idata->EwaldMarkers = NULL;
  allocatePinnedHostMemory((void **)&(h_idata->ewt), nEwhLoop*sizeof(EwtData));
  allocatePinnedHostMemory((void **)&(h_idata->cachedData), sizeof(EwaldReadOnlyData));
}

void EwaldHostMemoryFree(EwaldData *h_idata, int largephase) {
  if(largephase)
    freePinnedHostMemory(h_idata->EwaldMarkers);
  freePinnedHostMemory(h_idata->ewt);
  freePinnedHostMemory(h_idata->cachedData);
}

/** @brief Set up CUDA kernels to perform Ewald sum.
 *  @param d_localParts Local particle data on device
 *  @param d_localVars Local particle accelerations on device
 *  @param h_idata Host data buffers
 *  @param stream CUDA stream to perform GPU operations over
 *  @param cb Callback
 *  @param myIndex Chare index on this node that called this request.
 *  @param largephase Whether to perform large or small phase calculation
 *  
 *  The "top" and "bottom" Ewlad kernels have been combined:
 *    "top" for the real space loop,
 *    "bottom" for the k-space loop.
 *  
 */
void EwaldHost(CompactPartData *d_localParts, VariablePartData *d_localVars,
               EwaldData *h_idata, cudaStream_t stream, void *cb, int myIndex, int largephase)
{
  int n = h_idata->cachedData->n;
  int numBlocks = (int) ceilf((float)n/BLOCK_SIZE);
  int nEwhLoop = h_idata->cachedData->nEwhLoop;
  assert(nEwhLoop <= NEWH);

  size_t size;
  if(largephase) size = n * sizeof(int);
  else size = 0;

  HAPI_TRACE_BEGIN();
  int *d_EwaldMarkers;
  cudaChk(cudaMalloc(&d_EwaldMarkers, size));

  cudaMemcpyAsync(d_EwaldMarkers, h_idata->EwaldMarkers, size, cudaMemcpyHostToDevice, stream);
  cudaMemcpyToSymbolAsync(cachedData, h_idata->cachedData, sizeof(EwaldReadOnlyData), 0, cudaMemcpyHostToDevice, stream);
  cudaMemcpyToSymbolAsync(ewt, h_idata->ewt, nEwhLoop * sizeof(EwtData), 0, cudaMemcpyHostToDevice, stream);

  if (largephase)
      EwaldKernel<<<numBlocks, BLOCK_SIZE, 0, stream>>>(d_localParts, 
		                                         d_localVars,
							 d_EwaldMarkers, 1,
							 h_idata->EwaldRange[0], h_idata->EwaldRange[1]);
  else
      EwaldKernel<<<numBlocks, BLOCK_SIZE, 0, stream>>>(d_localParts, 
                                                         d_localVars,
							 NULL, 0,
							 h_idata->EwaldRange[0], h_idata->EwaldRange[1]);
  HAPI_TRACE_END(CUDA_EWALD);

#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
  printf("[%d] in EwaldHost, enqueued EwaldKernel\n", myIndex);
#endif

  cudaChk(cudaPeekAtLastError());
  hapiAddCallback(stream, cb);
  cudaChk(cudaFree(d_EwaldMarkers));
}

__global__ void EwaldKernel(CompactPartData *particleCores, 
                               VariablePartData *particleVars, 
                               int *markers, int largephase,
                               int First, int Last) {
  /////////////////////////////////////
  ////////////// Ewald TOP ////////////
  /////////////////////////////////////
  int id;
  if(largephase){
    id = blockIdx.x * BLOCK_SIZE + threadIdx.x;
    if(id > Last) return;
    id = markers[id];
  }else{
    id = First + blockIdx.x * BLOCK_SIZE + threadIdx.x;
    if(id > Last) return;
  }

  CompactPartData *p;

  cudatype alphan;
  cudatype fPot, ax, ay, az;
  cudatype x, y, z, r2, dir, dir2, a; 
  cudatype xdif, ydif, zdif; 
  cudatype g0, g1, g2, g3;
  cudatype Q2, Q2mirx, Q2miry, Q2mirz, Q2mir, Qta; 
  int ix, iy, iz, bInHole, bInHolex, bInHolexy;

#ifdef HEXADECAPOLE
  MomcData *mom = &(cachedData->momcRoot);
  MultipoleMomentsData *momQuad = &(cachedData->mm);
  cudatype xx,xxx,xxy,xxz,yy,yyy,yyz,xyy,zz,zzz,xzz,yzz,xy,xyz,xz,yz;
  cudatype g4, g5;
  cudatype Q4mirx,Q4miry,Q4mirz,Q4mir,Q4x,Q4y,Q4z;
  cudatype Q4xx,Q4xy,Q4xz,Q4yy,Q4yz,Q4zz,Q4,Q3x,Q3y,Q3z;
  cudatype Q3mirx,Q3miry,Q3mirz,Q3mir;
  const cudatype onethird = 1.0/3.0;
#else
  MultipoleMomentsData *mom;
  mom = &(cachedData->mm);
#endif

#ifdef HEXADECAPOLE
  Q4xx = 0.5*(mom->xxxx + mom->xxyy + mom->xxzz);
  Q4xy = 0.5*(mom->xxxy + mom->xyyy + mom->xyzz);
  Q4xz = 0.5*(mom->xxxz + mom->xyyz + mom->xzzz);
  Q4yy = 0.5*(mom->xxyy + mom->yyyy + mom->yyzz);
  Q4yz = 0.5*(mom->xxyz + mom->yyyz + mom->yzzz);
  Q4zz = 0.5*(mom->xxzz + mom->yyzz + mom->zzzz);
  Q4 = 0.25*(Q4xx + Q4yy + Q4zz);
  Q3x = 0.5*(mom->xxx + mom->xyy + mom->xzz);
  Q3y = 0.5*(mom->xxy + mom->yyy + mom->yzz);
  Q3z = 0.5*(mom->xxz + mom->yyz + mom->zzz);
#endif

  Q2 = 0.5 * (mom->xx + mom->yy + mom->zz);

  p = &(particleCores[id]);

#ifdef DEBUG
  if (blockIdx.x == 0 && threadIdx.x == 0) {
    printf("Moments\n");
    printf("xx %f, xy %f, xz %f, yy %f, yz %f, zz %f\n", mom->xx, mom->xy, mom->xz,
        mom->yy, mom->yz, mom->zz);
  }
#endif

  ax = 0.0f;
  ay = 0.0f;
  az = 0.0f;

#ifdef HEXADECAPOLE
  xdif = p->position.x - momQuad->cmx; 
  ydif = p->position.y - momQuad->cmy; 
  zdif = p->position.z - momQuad->cmz;
  fPot = momQuad->totalMass*cachedData->k1;
#else
  xdif = p->position.x - mom->cmx; 
  ydif = p->position.y - mom->cmy; 
  zdif = p->position.z - mom->cmz;
  fPot = mom->totalMass*cachedData->k1;
#endif
  for (ix=-(cachedData->nEwReps);ix<=(cachedData->nEwReps);++ix) {  
    bInHolex = (ix >= -cachedData->nReps && ix <= cachedData->nReps);
    x = xdif + ix * cachedData->L;
    for(iy=-(cachedData->nEwReps);iy<=(cachedData->nEwReps);++iy) {
      bInHolexy = (bInHolex && iy >= -cachedData->nReps && iy <= cachedData->nReps);
      y = ydif + iy*cachedData->L;
      for(iz=-(cachedData->nEwReps);iz<=(cachedData->nEwReps);++iz) {
        bInHole = (bInHolexy && iz >= -cachedData->nReps && iz <= cachedData->nReps);
        z = zdif + iz*cachedData->L;
        r2 = x*x + y*y + z*z;
        if (r2 > cachedData->fEwCut2 && !bInHole) continue;
        if (r2 < cachedData->fInner2) {

          /*
           * For small r, series expand about
           * the origin to avoid errors caused
           * by cancellation of large terms.
           * N.B. The following uses expf(), erfcf(), etc. If these ever
           * get changed to use double precision, then fInner2 also needs
           * to be changed. See line 480 in Ewald.cpp.
           */

          alphan = cachedData->ka;
          r2 *= cachedData->alpha2;
          g0 = alphan*((1.0/3.0)*r2 - 1.0);
          alphan *= 2*cachedData->alpha2;
          g1 = alphan*((1.0/5.0)*r2 - (1.0/3.0));
          alphan *= 2*cachedData->alpha2;
          g2 = alphan*((1.0/7.0)*r2 - (1.0/5.0));
          alphan *= 2*cachedData->alpha2;
          g3 = alphan*((1.0/9.0)*r2 - (1.0/7.0));
#ifdef HEXADECAPOLE
	  alphan *= 2*cachedData->alpha2;
	  g4 = alphan*((1.0/11.0)*r2 - (1.0/9.0));
	  alphan *= 2*cachedData->alpha2;
	  g5 = alphan*((1.0/13.0)*r2 - (1.0/11.0));
#endif
        }
        else {
          dir = 1/sqrtf(r2);
          dir2 = dir*dir;
          a = expf(-r2*cachedData->alpha2);
          a *= cachedData->ka*dir2;
          if (bInHole) g0 = -erff(cachedData->alpha/dir);
          else g0 = erfcf(cachedData->alpha/dir);
          g0 *= dir;
          g1 = g0*dir2 + a;
          alphan = 2*cachedData->alpha2;
          g2 = 3*g1*dir2 + alphan*a;
          alphan *= 2*cachedData->alpha2;
          g3 = 5*g2*dir2 + alphan*a;
#ifdef HEXADECAPOLE
	  alphan *= 2*cachedData->alpha2;
	  g4 = 7*g3*dir2 + alphan*a;
	  alphan *= 2*cachedData->alpha2;
	  g5 = 9*g4*dir2 + alphan*a;
#endif
        }
#ifdef HEXADECAPOLE
	xx = 0.5*x*x;
	xxx = onethird*xx*x;
	xxy = xx*y;
	xxz = xx*z;
	yy = 0.5*y*y;
	yyy = onethird*yy*y;
	xyy = yy*x;
	yyz = yy*z;
	zz = 0.5*z*z;
	zzz = onethird*zz*z;
	xzz = zz*x;
	yzz = zz*y;
	xy = x*y;
	xyz = xy*z;
	xz = x*z;
	yz = y*z;
	Q2mirx = mom->xx*x + mom->xy*y + mom->xz*z;
	Q2miry = mom->xy*x + mom->yy*y + mom->yz*z;
	Q2mirz = mom->xz*x + mom->yz*y + mom->zz*z;
	Q3mirx = mom->xxx*xx + mom->xxy*xy + mom->xxz*xz + mom->xyy*yy + mom->xyz*yz + mom->xzz*zz;
	Q3miry = mom->xxy*xx + mom->xyy*xy + mom->xyz*xz + mom->yyy*yy + mom->yyz*yz + mom->yzz*zz;
	Q3mirz = mom->xxz*xx + mom->xyz*xy + mom->xzz*xz + mom->yyz*yy + mom->yzz*yz + mom->zzz*zz;
	Q4mirx = mom->xxxx*xxx + mom->xxxy*xxy + mom->xxxz*xxz + mom->xxyy*xyy + mom->xxyz*xyz +
	  mom->xxzz*xzz + mom->xyyy*yyy + mom->xyyz*yyz + mom->xyzz*yzz + mom->xzzz*zzz;
	Q4miry = mom->xxxy*xxx + mom->xxyy*xxy + mom->xxyz*xxz + mom->xyyy*xyy + mom->xyyz*xyz +
	  mom->xyzz*xzz + mom->yyyy*yyy + mom->yyyz*yyz + mom->yyzz*yzz + mom->yzzz*zzz;
	Q4mirz = mom->xxxz*xxx + mom->xxyz*xxy + mom->xxzz*xxz + mom->xyyz*xyy + mom->xyzz*xyz +
	  mom->xzzz*xzz + mom->yyyz*yyy + mom->yyzz*yyz + mom->yzzz*yzz + mom->zzzz*zzz;
	Q4x = Q4xx*x + Q4xy*y + Q4xz*z;
	Q4y = Q4xy*x + Q4yy*y + Q4yz*z;
	Q4z = Q4xz*x + Q4yz*y + Q4zz*z;
	Q2mir = 0.5*(Q2mirx*x + Q2miry*y + Q2mirz*z) - (Q3x*x + Q3y*y + Q3z*z) + Q4;
	Q3mir = onethird*(Q3mirx*x + Q3miry*y + Q3mirz*z) - 0.5*(Q4x*x + Q4y*y + Q4z*z);
	Q4mir = 0.25*(Q4mirx*x + Q4miry*y + Q4mirz*z);
	Qta = g1*mom->m - g2*Q2 + g3*Q2mir + g4*Q3mir + g5*Q4mir;
	fPot -= g0*mom->m - g1*Q2 + g2*Q2mir + g3*Q3mir + g4*Q4mir;
	ax += g2*(Q2mirx - Q3x) + g3*(Q3mirx - Q4x) + g4*Q4mirx - x*Qta;
	ay += g2*(Q2miry - Q3y) + g3*(Q3miry - Q4y) + g4*Q4miry - y*Qta;
	az += g2*(Q2mirz - Q3z) + g3*(Q3mirz - Q4z) + g4*Q4mirz - z*Qta;
#else
        Q2mirx = mom->xx*x + mom->xy*y + mom->xz*z;
        Q2miry = mom->xy*x + mom->yy*y + mom->yz*z;
        Q2mirz = mom->xz*x + mom->yz*y + mom->zz*z;
        Q2mir = 0.5*(Q2mirx*x + Q2miry*y + Q2mirz*z);
        Qta = g1*mom->totalMass - g2*Q2 + g3*Q2mir;
        fPot -= g0*mom->totalMass - g1*Q2 + g2*Q2mir;

        ax += g2*(Q2mirx) - x*Qta;
        ay += g2*(Q2miry) - y*Qta;
        az += g2*(Q2mirz) - z*Qta;

#endif
      }
    }
  }

  /////////////////////////////////////
  //////////// Ewald Bottom ///////////
  /////////////////////////////////////
  cudatype hdotx, c, s;
  cudatype tempEwt; 
  xdif = ydif = zdif = 0.0; 

  MultipoleMomentsData *momBottom = &(cachedData->mm);

  /*
   ** Scoring for the h-loop (+,*)
   **        Without trig = (10,14)
   **          Trig est.    = 2*(6,11)  same as 1/sqrt scoring.
   **            Total        = (22,36)
   **                                   = 58
   */

  xdif = p->position.x - momBottom->cmx; 
  ydif = p->position.y - momBottom->cmy; 
  zdif = p->position.z - momBottom->cmz; 

  for (int i=0;i<cachedData->nEwhLoop;++i) {
    hdotx = ewt[i].hx * xdif + ewt[i].hy * ydif + ewt[i].hz * zdif;
    c = cosf(hdotx);
    s = sinf(hdotx);    
    fPot += ewt[i].hCfac*c + ewt[i].hSfac*s;    
    tempEwt = ewt[i].hCfac*s - ewt[i].hSfac*c;
    ax += ewt[i].hx * tempEwt;
    ay += ewt[i].hy * tempEwt;
    az += ewt[i].hz * tempEwt;
  }

  particleVars[id].a.x += ax;
  particleVars[id].a.y += ay;
  particleVars[id].a.z += az;
  particleVars[id].potential += fPot;
  
  return;
}

// initialize accelerations and potentials to zero
__global__ void ZeroVars(VariablePartData *particleVars, int nVars) {
    int id;
    id = blockIdx.x * BLOCK_SIZE + threadIdx.x;
    if(id >= nVars) return;

    particleVars[id].a.x = 0.0;
    particleVars[id].a.y = 0.0;
    particleVars[id].a.z = 0.0;
    particleVars[id].potential = 0.0;
    particleVars[id].dtGrav = 0.0;
}

void TreePieceODESolver(CudaSTIFF *d_CudaStiff, double  **y, double tstart, std::vector<double> *dtg, int numParts, cudaStream_t stream) {

    if (numParts == 0) return;

    double *d_y, *d_dtg;
    size_t ySize = numParts * 5 * sizeof(double); // TODO const defined in clIntegrateEnergy
    size_t dtgSize = dtg->size() * sizeof(double);

    cudaChk(cudaMalloc(&d_y, ySize));
    cudaChk(cudaMalloc(&d_dtg, dtgSize));

    double *y_host, *dtg_host, *y_host_out;
    allocatePinnedHostMemory((void **)&y_host, ySize);
    allocatePinnedHostMemory((void **)&y_host_out, ySize);
    allocatePinnedHostMemory((void **)&dtg_host, dtgSize);
    for (int i = 0; i < numParts; ++i) {
        memcpy(y_host + i * 5, y[i], 5 * sizeof(double));
    }
    memcpy(dtg_host, dtg->data(), dtgSize);

    cudaStreamSynchronize(stream);
    cudaChk(cudaMemcpyAsync(d_y, y_host, ySize, cudaMemcpyHostToDevice, stream));
    cudaChk(cudaMemcpyAsync(d_dtg, dtg_host, dtgSize, cudaMemcpyHostToDevice, stream));

    CudaStiffStep<<<numParts / THREADS_PER_BLOCK + 1, dim3(THREADS_PER_BLOCK), 0, stream>>>(d_CudaStiff, d_y, tstart, d_dtg, numParts);
    // Put these in pinned host memory
    cudaChk(cudaMemcpyAsync(y_host_out, d_y, ySize, cudaMemcpyDeviceToHost, stream));
    cudaStreamSynchronize(stream);

    // Copy data back to the original 2D array y
    for (int i = 0; i < numParts; ++i) {
        memcpy(y[i], y_host_out + i * 5, 5 * sizeof(double));
    }

    cudaFree(d_y);
    cudaFree(d_dtg);
    freePinnedHostMemory(y_host);
    freePinnedHostMemory(y_host_out);
    freePinnedHostMemory(dtg_host);
}

// These need to be copied to device memory
__device__ double AP_Gamma_HI_factor[] = { 0.99805271764596307, 0.99877911567687988, 0.99589340865612034,
			  0.99562060764857702, 0.99165170359332663, 0.9900889877822455,
			  0.98483276828954668, 0.97387675312245325, 0.97885673164000397,
			  0.98356305803821331, 0.96655786672182487, 0.9634906824933207,
			  0.95031917373653985, 0.87967606627349137, 0.79917533618355074,
			  0.61276011113763151, 0.16315185162187529, 0.02493663181368239,
			  0.0044013580765645335, 0.00024172553511936628, 1.9576102058649783e-10,
			  0.0, 0.0, 0.0, 0.0, 0.0 };

__device__ double AP_Gamma_HeI_factor[] = { 0.99284882980782224, 0.9946618686265097, 0.98641914356740497,
			   0.98867015777574848, 0.96519214493135597, 0.97188336387980656,
			   0.97529866247535113 , 0.97412477991428936, 0.97904139838765991,
			   0.98368372768570034, 0.96677432842215549, 0.96392622083382651,
			   0.95145730833093178, 0.88213871255482879, 0.80512823597731886,
			   0.62474472739578646, 0.17222786134467002, 0.025861959933038869,
			   0.0045265030237581529, 0.00024724339438128221, 1.3040144591221284e-08,
			   0.0, 0.0, 0.0, 0.0, 0.0};
 
 
__device__ double AP_Gamma_HeII_factor[] = { 0.97990208047216765, 0.98606251822654412,
			0.97657215632444849, 0.97274858503068629, 0.97416108746560681,
			0.97716929017896703, 0.97743607605974214, 0.97555305775319012,
			0.97874250764784809, 0.97849791914637996, 0.95135572977973504,
			0.92948461312852582, 0.89242272355549912, 0.79325512242742746 ,
			0.6683745597121028, 0.51605924897038324, 0.1840253816147828,
			0.035905775349044489, 0.0045537756654992923, 0.00035933897136804514,
			1.2294426136470751e-6, 0.0, 0.0, 0.0, 0.0, 0.0 };

__device__ double AP_Gamma_H2_factor[] = {1.0, 1.0, 1.0,
			  1.0, 1.0, 1.0,
			  1.0, 1.0, 1.0,
			  1.0, 1.0, 1.0,
			  1.0, 1.0, 1.0,
			  1.0, 1.0, 1.0,
			  1.0, 1.0, 1.0,
			  1.0, 1.0, 1.0, 1.0, 1.0}; 

__device__ double cudaClCoolLineH2_HI( double T){ /* H2-H collisionally-induced cooling, Glover & Abel 08, Table 8 */
  double a00 = -16.818342,
    a10 = 37.383713,
    a20 = 58.145166,
    a30 = 48.656103,
    a40 = 20.159831,
    a50 = 3.8479610;
  double a01 = -24.311209,
    a11 = 3.5692468,
    a21 = -11.332860,
    a31 = -27.850082,
    a41 = -21.328264,
    a51 = -4.2519023;
  double a02 = -24.311209,
    a12 = 4.6450521,
    a22 = -3.7209846,
    a32 = 5.9369081,
    a42 = -5.5108047,
    a52 = 1.5538288;
  double xint = 6000, slope = 2.10095, yint = 1.86368e-22;

  if (T <= 100) return pow(10.0, a00 + 
                                 a10*log10(T/1000.0) + 
                                 a20*log10(T/1000.0)*log10(T/1000.0) + 
                                 a30*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0) + 
                                 a40*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0) + 
                                 a50*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0));
  else if (T <= 1000) return pow(10.0, a01 + 
                                 a11*log10(T/1000.0) + 
                                 a21*log10(T/1000.0)*log10(T/1000.0) + 
                                 a31*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0) + 
                                 a41*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0) + 
                                 a51*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0));
  else if (T <= xint) return pow(10.0, a02 + 
                                 a12*log10(T/1000.0) + 
                                 a22*log10(T/1000.0)*log10(T/1000.0) + 
                                 a32*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0) + 
                                 a42*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0) + 
                                 a52*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)) ;
  else return pow(10.0,slope*(log10(T/1000.0) - log10(xint/1000.0)) + log10(yint));
}

__device__ double cudaCLCOOLLINEH2_H2( double T){  /* H2-H2 collisionally-induced cooling, Glover & Abel 08, Table 8*/
  double a0 = -23.962112,
    a1 = 2.09433740,
    a2 = -0.77151436,
    a3 = 0.43693353,
    a4 = -0.14913216,
    a5 = -0.033638326;
 double xint = 6000, slope = 1.34460, yint = 2.19802e-23;
  
 if (T <= xint) return pow(10.0, a0 + 
	     a1*log10(T/1000.0) + 
	     a2*log10(T/1000.0)*log10(T/1000.0) + 
	     a3*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0) + 
	     a4*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0) + 
	     a5*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0));
 else return pow(10.0,slope*(log10(T/1000.0) - log10(xint/1000.0)) + log10(yint));
}

__device__ double cudaCLCOOLLINEH2_HE( double T){ /* H2-He collisionally-induced cooling, Glover & Abel 08, Table 8 */
  double a0 = -23.689237,
    a1 = 2.1892372,
    a2 = -0.81520438,
    a3 = 0.29036281,
    a4 = -0.16596184,
    a5 = 0.19191375;
  double xint = 6000, slope = 1.48703, yint = 4.48145e-23;
  
  if (T <= xint) return pow(10.0, a0 + 
	     a1*log10(T/1000.0) + 
	     a2*log10(T/1000.0)*log10(T/1000.0) + 
	     a3*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0) + 
	     a4*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0) + 
	     a5*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0));
  else return pow(10.0,slope*(log10(T/1000.0) - log10(xint/1000.0)) + log10(yint));
}

__device__ double cudaClCoolLineH2_HII( double T){ /*H2-HII collisionally-induced cooling , Glover & Abel 08, Table 8 */
  double a0 = -21.716699,
    a1 = 1.3865783,
    a2 = -0.37915285,
    a3 = 0.11453688,
    a4 = -0.23214154,
    a5 = 0.058538864;
  double xint = 10000, slope = 0.336011, yint = 1.70474e-21;
  
  if (T <= xint) return pow(10.0, a0 + 
	     a1*log10(T/1000.0) + 
	     a2*log10(T/1000.0)*log10(T/1000.0) + 
	     a3*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0) + 
	     a4*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0) + 
	     a5*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0));
  else return pow(10.0,slope*(log10(T/1000.0) - log10(xint/1000.0)) + log10(yint));
}

__device__ double cudaCLCOOLLINEH2_E( double T){  /*Cooling based on H2-e collisionally-induced dissociation, Glover & Abel 08, Table 8*/
  double a00 = -34.286155,
    a10 = -48.537163,
    a20 = -77.121176,
    a30 = -51.352459,
    a40 = -15.169160,
    a50 = -0.98120322;
  double a01 = -22.190316,
    a11 = 1.5728955,
    a21 = -0.21335100,
    a31 = 0.96149759,
    a41 = -0.91023495,
    a51 = 0.13749749;
 double xint = 10000, slope = 1.07723, yint = 2.28029e-21;

  if (T <= 200) return pow(10.0, a00 + 
                                 a10*log10(T/1000.0) + 
                                 a20*log10(T/1000.0)*log10(T/1000.0) + 
                                 a30*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0) + 
                                 a40*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0) + 
                                 a50*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0));
  else if (T <= xint) return pow(10.0, a01 + 
                                 a11*log10(T/1000.0) + 
                                 a21*log10(T/1000.0)*log10(T/1000.0) + 
                                 a31*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0) + 
                                 a41*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0) + 
                                 a51*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0));
  else return pow(10.0,slope*(log10(T/1000.0) - log10(xint/1000.0)) + log10(yint));
}

__device__ double cudaCLRATEDUSTFORMH2( double z, double clump ) {
  double Rate_dust = 0;
  /* clump = 10.0; Ranges from 2-10 to 30-100, Gendin et al 2008 CC*/ 
  Rate_dust = 3.5e-17*z/ZSOLAR*clump; /*Formation rate coefficient of molecular hydrogen on dust, Gnedin et al 2008, Wolfire 2008, unit of cc per s, divide metallicity by solar metallicity (was 0.0177 in first iterations of code) CC*/  
  return Rate_dust;
 }

/*__device__ double sign(double val) {
    return (val > 0) - (val < 0);
}*/

/*__device__ double sign(double x, double y) {
    return (y >= 0) ? fabs(x) : -fabs(x);
}*/

__device__ double sign(double a, double b)
{
    double aabs = fabs(a);
    if(b >= 0.0) return aabs;
    else return -aabs;
}

__device__ double cudaCLTEMPERATURE( double Y_Total, double E ) {
  return E/(Y_Total*CL_Eerg_gm_degK3_2);
}

__device__ double cudaCLSELFSHIELD (double yH2, double h) {
  double x, column_denH2, omega_H2 = 0.2; 
      if (yH2 <= 0) return 1.0;
      column_denH2 = h*yH2;
      x = column_denH2/5e14;
      return (1 - omega_H2)/(1 + x)/(1 + x) + omega_H2/sqrt(1 + x)*exp(-0.00085*sqrt(1 + x));
}

__device__ double cudaCLDUSTSHIELD (double yHI, double yH2, double z, double h) {
  double column_denHI, column_denH2, sigmad = 2e-21; /*4e-21;*/
      if (yHI < 0) column_denHI = 0;
      else column_denHI = h*yHI;
      if (yH2 < 0) column_denH2 = 0;
      else column_denH2 = h*yH2;
      return exp(-1.0*sigmad*z/ZSOLAR*(column_denHI + 2.0*column_denH2));
}

#define CL_Ccomp0 0.565e-9
#define CL_Tcmb0  2.735
#define CL_Ccomp  (CL_Ccomp0*CL_Tcmb0)

void CudaCoolSetTime( CudaCOOL *cl, double dTime, double z, cudaStream_t stream ) {
    clRatesRedshift<<<1, 1, 0, stream>>>( cl, z, dTime );
}

__device__ void clSetAbundanceTotals(CudaCOOL *cl, double ZMetal, double *pY_H, double *pY_He, double *pY_eMax) {
    double Y_H, Y_He;
    
    if (ZMetal <= 0.1) {
        Y_He = (0.236 + 2.1*ZMetal)/4.0;
	}
    else {
        Y_He = (-0.446*(ZMetal - 0.1)/0.9 + 0.446)/4.0;
	}
    Y_H = 1.0 - Y_He*4.0 - ZMetal; 

    *pY_H = Y_H;
    *pY_He = Y_He;
    *pY_eMax = Y_H+ Y_He*2; /* Ignoring any electrons from metals */

}

__global__ void clRatesRedshift( CudaCOOL *cl, double zIn, double dTimeIn ) {
  int i;
  double xx;
  double zTime;
  CudaUVSPECTRUM *UV,*UV0;
  double Y_H, Y_He, Y_eMax;
  double ten = 10.0,output, expon;

  /* printf("Redshift: %f \n", zIn); */ 
  /* cl->z = 0.0; */ 
  cl->z = zIn; 
  cl->dTime = dTimeIn;
  cl->dComovingGmPerCcUnit = cl->dGmPerCcUnit*pow(1.+zIn,3.);
  cl->dExpand = 1.0/(1.0+zIn);

  cl->R.Cool_Comp = pow((1+zIn)*CL_Ccomp,4.0)*CL_B_gm; 
  cl->R.Tcmb = CL_Tcmb0*(1+zIn);
  clSetAbundanceTotals(cl,0.0,&Y_H,&Y_He,&Y_eMax); /* Hack to estimate Y_H */
  cl->R.Cool_LowTFactor = (cl->bLowTCool ? CL_B_gm*Y_H*Y_H/0.001 : 0 );

  /* Photo-Ionization rates */

  UV = cl->UV;

  if (cl->bUV) {
	  assert( UV != NULL );
	  if (cl->bUVTableUsesTime) {
		  /*
		   ** Table in order of increasing time
		   */
		  zTime = dTimeIn;
		  for ( i=0; i < cl->nUV && zTime >= UV->zTime ; i++,UV++ );
		  }
	  else {
		  /*
		   ** Table in order of high to low redshift 
		   */
		  zTime = zIn;
		  for ( i=0; i < cl->nUV && zTime <= UV->zTime ; i++,UV++ );
		  }
	  }

  if (!cl->bUV || i==0) {
	  cl->R.Rate_Phot_HI = CL_RT_MIN;
	  cl->R.Rate_Phot_HeI = CL_RT_MIN;
	  cl->R.Rate_Phot_HeII = CL_RT_MIN;
	  cl->R.Rate_Phot_H2_cosmo = CL_RT_MIN; 
	  
	  cl->R.Heat_Phot_HI = 0.0;
	  cl->R.Heat_Phot_HeI = 0.0;
	  cl->R.Heat_Phot_HeII = 0.0;
	  cl->R.Heat_Phot_H2 = 0.0; 
	  return;
	  }
  
  UV0=UV-1;
  if (i == cl->nUV ) {
	  cl->R.Rate_Phot_HI = UV0->Rate_Phot_HI;
	  cl->R.Rate_Phot_HeI = UV0->Rate_Phot_HeI;
	  cl->R.Rate_Phot_HeII = UV0->Rate_Phot_HeII;
	  expon = 0.90632725*zTime - 0.16790918*zTime*zTime + 0.010241484*zTime*zTime*zTime - 12.518825;
	  output = pow(ten,expon); /*haardt_madau gal+quasar, z = 0*/ 
	  cl->R.Rate_Phot_H2_cosmo = output;
	  cl->R.Heat_Phot_HI = UV0->Heat_Phot_HI*CL_B_gm;
	  cl->R.Heat_Phot_HeI = UV0->Heat_Phot_HeI*CL_B_gm;
	  cl->R.Heat_Phot_HeII = UV0->Heat_Phot_HeII*CL_B_gm;
	  cl->R.Heat_Phot_H2 = 6.4e-13*CL_B_gm;
	  }
  else {
	  if (cl->bUVTableLinear) { /* use Linear interpolation */	
		  xx = (zTime - UV0->zTime)/(UV->zTime - UV0->zTime);
		  cl->R.Rate_Phot_HI = UV0->Rate_Phot_HI*(1-xx)+UV->Rate_Phot_HI*xx;
		  cl->R.Rate_Phot_HeI = UV0->Rate_Phot_HeI*(1-xx)+UV->Rate_Phot_HeI*xx;
		  cl->R.Rate_Phot_HeII = UV0->Rate_Phot_HeII*(1-xx)+UV->Rate_Phot_HeII*xx;
		  expon = 0.90632725*zTime - 0.16790918*zTime*zTime + 0.010241484*zTime*zTime*zTime - 12.518825;
		  cl->R.Rate_Phot_H2_cosmo = pow(ten,expon); /*haardt_madau gal+quasar, z = 0*/ 
		  cl->R.Heat_Phot_HI = (UV0->Heat_Phot_HI*(1-xx)+UV->Heat_Phot_HI*xx)*CL_B_gm;
		  cl->R.Heat_Phot_HeI = (UV0->Heat_Phot_HeI*(1-xx)+UV->Heat_Phot_HeI*xx)*CL_B_gm;
		  cl->R.Heat_Phot_HeII = (UV0->Heat_Phot_HeII*(1-xx)+UV->Heat_Phot_HeII*xx)*CL_B_gm;
		  cl->R.Heat_Phot_H2 = 6.4e-13*CL_B_gm;
		  }
	  else { /* use Log interpolation with 1+zTime */
		  xx = log((1+zTime)/(1+UV0->zTime))/log((1+UV->zTime)/(1+UV0->zTime));
		  cl->R.Rate_Phot_HI = pow(UV0->Rate_Phot_HI,1-xx)*pow(UV->Rate_Phot_HI,xx);
		  cl->R.Rate_Phot_HeI = pow(UV0->Rate_Phot_HeI,1-xx)*pow(UV->Rate_Phot_HeI,xx);
		  cl->R.Rate_Phot_HeII = pow(UV0->Rate_Phot_HeII,1-xx)*pow(UV->Rate_Phot_HeII,xx);
		  expon = 0.90632725*zTime - 0.16790918*zTime*zTime + 0.010241484*zTime*zTime*zTime - 12.518825;
		  cl->R.Rate_Phot_H2_cosmo = pow(ten,expon);  /*haardt_madau gal+quasar, z = 0*/ 
		  
		  cl->R.Heat_Phot_HI = pow(UV0->Heat_Phot_HI,1-xx)*pow(UV->Heat_Phot_HI,xx)*CL_B_gm;
		  cl->R.Heat_Phot_HeI = pow(UV0->Heat_Phot_HeI,1-xx)*pow(UV->Heat_Phot_HeI,xx)*CL_B_gm;
		  cl->R.Heat_Phot_HeII = pow(UV0->Heat_Phot_HeII,1-xx)*pow(UV->Heat_Phot_HeII,xx)*CL_B_gm;
		  cl->R.Heat_Phot_H2 = 6.4e-13*CL_B_gm;
		  }
	  }
  if (cl->R.Rate_Phot_HI < CL_RT_MIN) cl->R.Rate_Phot_HI = CL_RT_MIN;
  if (cl->R.Rate_Phot_HeI < CL_RT_MIN) cl->R.Rate_Phot_HeI = CL_RT_MIN;
  if (cl->R.Rate_Phot_HeII < CL_RT_MIN) cl->R.Rate_Phot_HeII = CL_RT_MIN;
  if (cl->R.Rate_Phot_H2_cosmo < CL_RT_MIN) cl->R.Rate_Phot_H2_cosmo = CL_RT_MIN; 

  return;
  }

/* Returns Heating - Cooling excluding External Heating, units of ergs s^-1 g^-1 
   Public interface CoolEdotInstantCode */
__device__ double cudaClEdotInstant_Table( CudaCOOL *cl, CudaPERBARYON *Y, CudaRATE *Rate, double rho, 
			    double ZMetal, double *dEdotHeat, double *dEdotCool )
{
  double en_B = rho*CL_B_gm;
  double xTln,wTln0,wTln1;/*,wTln0d,wTln1d;*/
  CudaRATES_T *RT0,*RT1;/*,*RT0d,*RT1d;*/
  int iTln;

  double ne,LowTCool;
  double s_dust, s_self, Rate_Phot_HI;
  s_dust = cudaCLDUSTSHIELD(Y->HI*en_B, Y->H2*en_B, ZMetal, Rate->CorreLength);
  s_self = cudaCLSELFSHIELD(Y->H2*en_B, Rate->CorreLength);
  if (cl->bShieldHI) Rate_Phot_HI = Rate->Phot_HI*s_dust;
  else Rate_Phot_HI = Rate->Phot_HI;

  ne = Y->e*en_B;

  xTln = (Rate->Tln-cl->TlnMin)*cl->rDeltaTln;
  iTln = xTln;
  RT0 = (cl->RT+iTln*TABLEFACTOR);
  RT1 = RT0+TABLEFACTOR; 
  xTln = xTln-iTln;
#ifdef CUBICTABLEINTERP
  RT0d = RT0+1;
  RT1d = RT1+1;
  {
  double x2 = xTln*xTln;
  wTln1 = x2*(3-2*xTln);
  wTln0 = 1-wTln1;
  wTln0d = xTln*(1+xTln*(xTln-2));
  wTln1d = x2*(xTln-1);
/*  wTln1 = xTln;
  wTln0 = 1-xTln;
  wTln0d = 0;
  wTln1d = 0;*/
  }
#else  
  wTln1 = xTln;
  wTln0 = 1-xTln;
#endif

#define DTFRACLOWTCOOL 0.25
  if (Rate->T > cl->R.Tcmb*(1+DTFRACLOWTCOOL))
      LowTCool = TABLEINTERP( Cool_LowT )*cl->R.Cool_LowTFactor*en_B*ZMetal;
  else if (Rate->T < cl->R.Tcmb*(1-DTFRACLOWTCOOL))
      LowTCool = -TABLEINTERP( Cool_LowT )*cl->R.Cool_LowTFactor*en_B*ZMetal;
  else {
      double x = (Rate->T/cl->R.Tcmb-1)*(1./DTFRACLOWTCOOL);
      LowTCool = -TABLEINTERP( Cool_LowT )*cl->R.Cool_LowTFactor*en_B*ZMetal
	  *x*(3-3*fabs(x)+x*x);
      }

  *dEdotCool = 
#ifndef NOCOMPTON
      Y->e * cl->R.Cool_Comp * ( Rate->T - cl->R.Tcmb ) 
#endif
    + 
    ne * (TABLEINTERP( Cool_Brem_1 ) * ( Y->HII + Y->HeII ) +
	  TABLEINTERP( Cool_Brem_2 ) * Y->HeIII +
	  
	  cl->R.Cool_Diel_HeII * Y->HeII * Rate->Diel_HeII +

	  TABLEINTERP( Cool_Line_HI ) * Y->HI +
	  TABLEINTERP( Cool_Line_HeI ) * Y->HeI +
	  TABLEINTERP( Cool_Line_HeII ) * Y->HeII +

	  TABLEINTERP( Cool_Radr_HII ) * Y->HII * Rate->Radr_HII  +
	  TABLEINTERP( Cool_Radr_HeII ) * Y->HeII * Rate->Radr_HeII +
	  TABLEINTERP( Cool_Radr_HeIII ) * Y->HeIII * Rate->Radr_HeIII +

	  cudaCLCOOLLINEH2_E(Rate->T) * Y->H2  * CL_B_gm /*Re-added CL_B_gm 9/16/10 */+
	  cl->R.Cool_Coll_H2 * Y->H2 * Rate->Coll_e_H2 /* s_dust * s_self*/  /*shielding of collisions 11/18/10*/ +
	  cl->R.Cool_Coll_HI * Y->HI * Rate->Coll_HI +
	  cl->R.Cool_Coll_HeI * Y->HeI * Rate->Coll_HeI + 
	  cl->R.Cool_Coll_HeII * Y->HeII * Rate->Coll_HeII )
    +
    cudaClCoolLineH2_HI(Rate->T) /* en_B * Y->H2 * Y->HI * CL_B_gm*/ /*Re-added CL_B_gm 9/16/10:
 cudaClCoolLineH2_HI(Rate->T) [erg cm^3 s^-1 H2atom^-1 HIatom^-1] * en_B [atom/cm^3] * Y->H2 [H2atom atom^-1] * Y->HI [HIatom atom^-1] CL_B_gm [atom g^-1] => [ergs s^-1 g^-1]*/
    +
    cudaCLCOOLLINEH2_H2(Rate->T) /* en_B * Y->H2 * Y->H2 * CL_B_gm */
    +
    cudaCLCOOLLINEH2_HE(Rate->T) /* en_B * Y->H2 * Y->HeI * CL_B_gm */ 
    +
    cudaClCoolLineH2_HII(Rate->T) /* en_B * Y->H2 * Y->HII * CL_B_gm */
    +
    cl->R.Cool_Coll_H2 * Y->H2 * Rate->Coll_HI_H2 * Y->HI *en_B 
    +
    cl->R.Cool_Coll_H2 * Y->H2 * Rate->Coll_H2_H2 * Y->H2 *en_B 
    +
    cl->R.Cool_Coll_H2 * Y->H2 * Rate->Coll_HII_H2 * Y->HII *en_B 
    + 
      LowTCool
#ifndef NOMETALCOOLING
    +
    Rate->Cool_Metal
#endif
    ;

  *dEdotHeat =
#ifndef NOMETALCOOLING
      Rate->Heat_Metal
    +
#endif
    Y->H2 * cl->R.Heat_Phot_H2 * Rate->Phot_H2*s_dust*s_self + /*photon heating and dissociation CC*/
    Y->HI   * cl->R.Heat_Phot_HI * Rate_Phot_HI +
    Y->HeI  * cl->R.Heat_Phot_HeI * Rate->Phot_HeI +
    Y->HeII * cl->R.Heat_Phot_HeII * Rate->Phot_HeII;

  return *dEdotHeat - *dEdotCool;
}

__device__ void cudaClRates_Table( CudaCOOL *cl, CudaRATE *Rate, double T, double rho, double ZMetal, double columnL, double Rate_Phot_H2_stellar) {
  double Tln;
  double xTln,wTln0,wTln1;/*,wTln0d,wTln1d;*/
  CudaRATES_T *RT0,*RT1;/*,*RT0d,*RT1d;*/
  int iTln;

  if (T >= cl->TMax) T=cl->TMax*(1.0 - EPS);   
  if (T < cl->TMin) T=cl->TMin;
  Tln = log(T);

  Rate->T = T;
  Rate->Tln = Tln; 

  xTln = (Tln-cl->TlnMin)*cl->rDeltaTln;
  iTln = xTln;
  RT0 = (cl->RT+iTln*TABLEFACTOR);
  RT1 = RT0+TABLEFACTOR; 
  xTln = xTln-iTln;
#ifdef CUBICTABLEINTERP
  RT0d = RT0+1;
  RT1d = RT1+1;
  {
  double x2 = xTln*xTln;
  wTln1 = x2*(3-2*xTln);
  wTln0 = 1-wTln1;
  wTln0d = xTln*(1+xTln*(xTln-2));
  wTln1d = x2*(xTln-1);
  }
#else  
  wTln1 = xTln;
  wTln0 = 1-xTln;
#endif
  Rate->Coll_HI = TABLEINTERP( Rate_Coll_HI );
  Rate->Coll_HeI = TABLEINTERP( Rate_Coll_HeI );
  Rate->Coll_HeII = TABLEINTERP( Rate_Coll_HeII );
  Rate->Coll_e_H2 = TABLEINTERP( Rate_Coll_e_H2 );
  Rate->Coll_HI_H2 = TABLEINTERP( Rate_Coll_HI_H2 );
  Rate->Coll_H2_H2 = TABLEINTERP( Rate_Coll_H2_H2 );
  Rate->Coll_HII_H2 = TABLEINTERPLIN(Rate_Coll_HII_H2);        
  Rate->Coll_Hm_e = TABLEINTERPLIN(Rate_Coll_Hm_e);          /*gas phase form of H2 */
  Rate->Coll_Hm_HII = TABLEINTERPLIN(Rate_Coll_Hm_HII);           /*gas phase form of H2 */
  Rate->HI_e = TABLEINTERPLIN(Rate_HI_e);          /*gas phase form of H2 */
  Rate->HI_Hm = TABLEINTERPLIN(Rate_HI_Hm);          /*gas phase form of H2 */
  Rate->Radr_HII = TABLEINTERP( Rate_Radr_HII );
  Rate->Radr_HeII = TABLEINTERP( Rate_Radr_HeII );
  Rate->Diel_HeII = TABLEINTERP( Rate_Diel_HeII );
  Rate->Chtr_HeII = TABLEINTERP( Rate_Chtr_HeII );
  Rate->Totr_HeII = Rate->Radr_HeII + Rate->Diel_HeII + Rate->Chtr_HeII;
  Rate->Radr_HeIII = TABLEINTERP( Rate_Radr_HeIII );
  Rate->DustForm_H2 =  cudaCLRATEDUSTFORMH2( ZMetal, cl->dClump );

  Rate->Phot_HI = cl->R.Rate_Phot_HI;
  Rate->Phot_HeI = cl->R.Rate_Phot_HeI;
  Rate->Phot_HeII = cl->R.Rate_Phot_HeII;

  if (Rate_Phot_H2_stellar < 20) {
    Rate_Phot_H2_stellar = 20;
  }
  Rate->Phot_H2 = cl->R.Rate_Phot_H2_cosmo + pow(10,Rate_Phot_H2_stellar - 59.626264)*cl->dMsolUnit*cl->dLymanWernerFrac/cl->dKpcUnit/cl->dKpcUnit/cl->dExpand/cl->dExpand;
  Rate->CorreLength = columnL * cl->dKpcUnit*KPCCM * cl->dExpand; /* From system units to cm*/
  Rate->LymanWernerCode = Rate_Phot_H2_stellar;
  if (cl->bSelfShield) {
      double logen_B;
      logen_B = log10(rho*CL_B_gm);
      if (logen_B > 2.2499) {
	  Rate->Phot_HI = 0;
	  Rate->Phot_HeI = 0;
	  Rate->Phot_HeII = 0;
	  Rate->Phot_H2 = 0;
	  }
      else if (logen_B > -10.25) {
	  double x = (logen_B+10.25)*2.0;
	  int ix;
	  ix = floor(x);
	  x -= ix;
	  Rate->Phot_HI *= (AP_Gamma_HI_factor[ix]*(1-x)+AP_Gamma_HI_factor[ix+1]*x);
	  Rate->Phot_HeI *= (AP_Gamma_HeI_factor[ix]*(1-x)+AP_Gamma_HeI_factor[ix+1]*x);
	  Rate->Phot_HeII *= (AP_Gamma_HeII_factor[ix]*(1-x)+AP_Gamma_HeII_factor[ix+1]*x);
	  Rate->Phot_H2 *= (AP_Gamma_H2_factor[ix]*(1-x)+AP_Gamma_H2_factor[ix+1]*x); /* Untested */
	  }
      }
}

__device__ void cudaClRateMetalTable(CudaCOOL *cl, CudaRATE *Rate, double T, double rho, double Y_H, double ZMetal)
{
  double tempT, tempnH,tempz, nH;
  double Tlog, nHlog; 
  double xTlog, wTlog0, wTlog1, xz, wz0, wz1, xnHlog, wnHlog0, wnHlog1; 
  int    iTlog, iz, inHlog; 
  double Cool000, Cool010, Cool100, Cool110, Cool001, Cool011, Cool101, Cool111; 
  double  Cool00, Cool01, Cool10, Cool11, Cool0, Cool1, Cool;
  
  double  Heat000, Heat010, Heat100, Heat110, Heat001, Heat011, Heat101, Heat111; 
  double  Heat00, Heat01, Heat10, Heat11, Heat0, Heat1, Heat;

  
  if(!cl->bMetal) {
    Rate->Cool_Metal = 0.0; 
    Rate->Heat_Metal = 0.0; 
    return; 
  }
 
  nH = rho*Y_H/M_H;
  
  tempT = T; 
  tempnH = nH;
  tempz = cl->z; 

  if (T >= cl->MetalTMax) tempT = cl->MetalTMax*(1.0-EPS);
  if (T < cl->MetalTMin) tempT=cl->MetalTMin;

  if (nH >= cl->MetalnHMax) tempnH = cl->MetalnHMax*(1.0-EPS);
  if (nH < cl->MetalnHMin) tempnH = cl->MetalnHMin; 
 
  if (cl->z <= cl->MetalzMin) tempz = cl->MetalzMin+EPS;
  /* if redshift is too high or no UV, use the no UV metal cooling table*/
  if (cl->z > cl->MetalzMax || !cl->bUV) tempz = cl->MetalzMax;   

  Tlog = log10(tempT); 
  nHlog = log10(tempnH); 
  
  xz = cl->nzMetalTable -1 - (tempz - cl->MetalzMin)*cl->rDeltaz; 
  iz = xz;   

  xTlog = (Tlog - cl->MetalTlogMin)*cl->rDeltaTlog; 
  assert(xTlog >= 0.0);
  iTlog = xTlog; 

  xnHlog = (nHlog - cl->MetalnHlogMin)*cl->rDeltanHlog; 
  inHlog = xnHlog;
  if (inHlog == cl->nnHMetalTable - 1) inHlog = cl->nnHMetalTable - 2; /*CC; To prevent running over the table.  Should not be used*/
  
  int nnH = cl->nnHMetalTable;
  int nt = cl->nTMetalTable;
  Cool000 = cl->MetalCoolln[(iz * nnH * nt) + (inHlog * nt) + iTlog];
  Cool001 = cl->MetalCoolln[(iz * nnH * nt) + (inHlog * nt) + (iTlog + 1)];
  Cool010 = cl->MetalCoolln[(iz * nnH * nt) + ((inHlog + 1) * nt) + iTlog];
  Cool011 = cl->MetalCoolln[(iz * nnH * nt) + ((inHlog + 1) * nt) + (iTlog + 1)];
  Cool100 = cl->MetalCoolln[((iz + 1) * nnH * nt) + (inHlog * nt) + iTlog];
  Cool101 = cl->MetalCoolln[((iz + 1) * nnH * nt) + (inHlog * nt) + (iTlog + 1)];
  Cool110 = cl->MetalCoolln[((iz + 1) * nnH * nt) + ((inHlog + 1) * nt) + iTlog];
  Cool111 = cl->MetalCoolln[((iz + 1) * nnH * nt) + ((inHlog + 1) * nt) + (iTlog + 1)];
  
  Heat000 = cl->MetalHeatln[(iz * nnH * nt) + (inHlog * nt) + iTlog];
  Heat001 = cl->MetalHeatln[(iz * nnH * nt) + (inHlog * nt) + (iTlog + 1)];
  Heat010 = cl->MetalHeatln[(iz * nnH * nt) + ((inHlog + 1) * nt) + iTlog];
  Heat011 = cl->MetalHeatln[(iz * nnH * nt) + ((inHlog + 1) * nt) + (iTlog + 1)];
  Heat100 = cl->MetalHeatln[((iz + 1) * nnH * nt) + (inHlog * nt) + iTlog];
  Heat101 = cl->MetalHeatln[((iz + 1) * nnH * nt) + (inHlog * nt) + (iTlog + 1)];
  Heat110 = cl->MetalHeatln[((iz + 1) * nnH * nt) + ((inHlog + 1) * nt) + iTlog];
  Heat111 = cl->MetalHeatln[((iz + 1) * nnH * nt) + ((inHlog + 1) * nt) + (iTlog + 1)];

  xz = xz - iz; 
  wz1 = xz; 
  wz0 = 1-xz; 
  
  Cool00 = wz0*Cool000 + wz1*Cool100;
  Cool01 = wz0*Cool001 + wz1*Cool101; 
  Cool10 = wz0*Cool010 + wz1*Cool110; 
  Cool11 = wz0*Cool011 + wz1*Cool111;

  Heat00 = wz0*Heat000 + wz1*Heat100;
  Heat01 = wz0*Heat001 + wz1*Heat101; 
  Heat10 = wz0*Heat010 + wz1*Heat110; 
  Heat11 = wz0*Heat011 + wz1*Heat111;


  xnHlog = xnHlog - inHlog; 
  wnHlog1 = xnHlog; 
  wnHlog0 = 1-xnHlog; 

  Cool0 = wnHlog0*Cool00 + wnHlog1*Cool10; 
  Cool1 = wnHlog0*Cool01 + wnHlog1*Cool11; 

  
  Heat0 = wnHlog0*Heat00 + wnHlog1*Heat10; 
  Heat1 = wnHlog0*Heat01 + wnHlog1*Heat11; 
  
  xTlog = xTlog - iTlog; 
  wTlog1 = xTlog; 
  wTlog0 = 1 - xTlog; 
  Cool = wTlog0*Cool0 + wTlog1*Cool1; 
  Heat = wTlog0*Heat0 + wTlog1*Heat1; 
    /* convert unit to erg/g/sec, time a factor of nH^2/nH, also scale with metalicity */ 
  Rate->Cool_Metal = exp(Cool)*nH*Y_H/M_H * ZMetal/ZSOLAR; 
  Rate->Heat_Metal = exp(Heat)*nH*Y_H/M_H * ZMetal/ZSOLAR;   

}

__device__ void CudaclDerivs(double x, const double *y, double *dGain, double *dLoss,
	     void *Data) {
  CudaclDerivsData *d = (CudaclDerivsData *)Data;
  double T,ne,nHI, nHII, internalheat = 0, externalheat = 0;
  double internalcool = 0;
  double en_B = d->rho*CL_B_gm;
  double s_dust, s_self, nH2, nHminus, Rate_Phot_HI;
  
  d->E = y[0];
  d->Y.HI = y[1];
  d->Y.H2 = y[4];
  d->Y.HII = d->Y_H - d->Y.HI - d->Y.H2*2.0; /*Don't forget about molec H CC*/
  if(d->Y.HII < 0) d->Y.HII = 0;  
 
  d->Y.HeI = y[2];
  d->Y.HeII = y[3];
  d->Y.HeIII = d->Y_He - d->Y.HeI - d->Y.HeII; 
  if(d->Y.HeIII < 0) d->Y.HeIII = 0; 
 
  d->Y.e = d->Y.HII + d->Y.HeII + 2*d->Y.HeIII;
#ifdef Y_EMIN
  if (d->Y.e < Y_EMIN) d->Y.e = Y_EMIN;
#endif

  d->Y.Total = d->Y.e + d->Y_H + d->Y_He + d->ZMetal/MU_METAL - d->Y.H2;  /* H total from cl now -- in future from particle */ 
  T = cudaCLTEMPERATURE( d->Y.Total, d->E );
  CLRATES( d->cl, &d->Rate, T, d->rho, d->ZMetal, d->columnL, d->dLymanWerner);
  s_dust = cudaCLDUSTSHIELD(d->Y.HI*en_B, d->Y.H2*en_B, d->ZMetal, d->Rate.CorreLength); /* dust shielding*/
  s_self = cudaCLSELFSHIELD(d->Y.H2*en_B, d->Rate.CorreLength); /* H2 self shielding*/
  if (d->cl->bShieldHI) Rate_Phot_HI = d->Rate.Phot_HI*s_dust;/* If bShieldHI was set in parameter file, modify the photoionizing flux for HI by the dust shielding*/
  else Rate_Phot_HI = d->Rate.Phot_HI;

  externalheat = d->ExternalHeating;
  if (d->bCool) {
    cudaClRateMetalTable(d->cl, &d->Rate, T, d->rho, d->Y_H, d->ZMetal);
    CLEDOTINSTANT( d->cl, &d->Y, &d->Rate, d->rho, d->ZMetal, &internalheat, &internalcool );
    
    dGain[0] = internalheat;
    dLoss[0] = internalcool;
  }
  if(externalheat > 0.0)
      dGain[0] += externalheat;
  else
      dLoss[0] -= externalheat;

  ne =  en_B*d->Y.e;
  nHI = en_B*d->Y.HI;
  nH2 = en_B*d->Y.H2;
  nHII = en_B*d->Y.HII;
  nHminus = d->Rate.HI_e* d->Y.HI*d->Y.e/(d->Rate.HI_Hm*d->Y.HI + d->Rate.Coll_Hm_HII*d->Y.HI + d->Rate.Coll_Hm_e*d->Y.e);

  dGain[4] = d->Y.HI*nHminus*en_B*d->Rate.HI_Hm + /*gas phase formation of H2, adel 97 */
            (d->Y.HI + 2.0*d->Y.H2)*d->Rate.DustForm_H2*nHI;
  dLoss[4] = d->Y.H2*d->Rate.Phot_H2*s_dust*s_self +
	     d->Y.H2*d->Rate.Coll_e_H2*ne +
	     d->Y.H2*d->Rate.Coll_H2_H2*nH2 +
	     d->Y.H2*d->Rate.Coll_HI_H2*nHI +
             d->Y.H2*d->Rate.Coll_HII_H2*nHII;
  dGain[1] = ne*d->Y.HII*d->Rate.Radr_HII +
             2.0*dLoss[4]; /*Photo or collisionally dissociated H2*/
  dLoss[1] = ne*d->Y.HI*d->Rate.Coll_HI +
             d->Y.HI*Rate_Phot_HI + /*Adding in molec H and shielding, should possibly add shielding for others*/
             2.0*dGain[4]; /*Formation of H2*/
  dGain[2] = ne*d->Y.HeII*d->Rate.Totr_HeII;
  dLoss[2] = ne*d->Y.HeI*d->Rate.Coll_HeI + 
             d->Y.HeI*d->Rate.Phot_HeI;
  dGain[3] = ne*d->Y.HeIII*d->Rate.Radr_HeIII +
             dLoss[2];
  dLoss[3] = ne*d->Y.HeII*d->Rate.Coll_HeII + 
             d->Y.HeII*d->Rate.Phot_HeII +
             dGain[2];
}

__global__ void CudaStiffStep(CudaSTIFF *s0, double *_y, double tstart, double *_dtg, int nVars) {
    int id;
    id = blockIdx.x * BLOCK_SIZE + threadIdx.x;
    if(id >= nVars) return;

    double dtg = _dtg[id]; 

    double tn;			/* time within step */
    int i;

    CudaSTIFF *s = &s0[id];

    /*
     * Local copies of Stiff context
     */
    int n = s->nv;
    double *y = &_y[id * s->nv];

    double *y0 = &s->y0[id * s->nv];
    double *ymin = &s->ymin[id * s->nv];
    double *q = &s->q[id * s->nv];
    double *d = &s->d[id * s->nv];
    double *rtau = &s->rtau[id * s->nv];
    double *ys = &s->ys[id * s->nv];
    double *qs = &s->qs[id * s->nv];
    double *rtaus = &s->rtaus[id * s->nv];
    double *scrarray = &s->scrarray[id * s->nv];
    double *y1 = &s->y1[id * s->nv];

    double epsmin = s->epsmin;
    double sqreps = s->sqreps;
    double epscl = s->epscl;
    double epsmax = s->epsmax;
    double dtmin = s->dtmin;
    int itermax = s->itermax;
    int gcount = 0;		/* count calls to derivs */
    int rcount = 0;		/* count restart steps */
    double scrtch;
    double ascr;
    double scr1;
    double scr2;
    double dt;			/* timestep used by the integrator */
    double ts;			/* t at start of the chemical timestep */
    double alpha;		/* solution parameter used in update */
    int iter;			/* counter for corrector iterations */
    double eps;			/* maximum correction term */
    double rtaub;
    double qt;			/* alpha weighted average of q */
    double pb;
    const double tfd = 1.000008; /* fudge for completion of timestep */
    double rteps;		/* estimate of sqrt(eps) */
    double dto;			/* old timestep; to rescale rtaus */
    double temp;                
    
    tn = 0.0;
    for(i = 0; i < n; i++) {
	q[i] = 0.0;
	d[i] = 0.0;
	y0[i] = y[i];
	y[i] = max(y[i], ymin[i]);
	}
    
    //s->derivs(tn + tstart, y, q, d, s->Data);
    CudaclDerivs(tn + tstart, y, q, d, &((CudaclDerivsData*)s->Data)[id]);
    //printf("q: %g\n", q[0]);
    gcount++;
    
    /*
    C
    c estimate the initial stepsize.
    C
    c strongly increasing functions(q >>> d assumed here) use a step-
    c size estimate proportional to the step needed for the function to
    c reach equilibrium where as functions decreasing or in equilibrium
    c use a stepsize estimate directly proportional to the character-
    c istic stepsize of the function. convergence of the integration
    c scheme is likely since the smallest estimate is chosen for the
    c initial stepsize.
    */
    scrtch  = 1.0e-25;
    for(i = 0; i < n; i++) {
	ascr = fabs(q[i]);
	scr2 = sign(1./y[i],.1*epsmin*ascr - d[i]);
	scr1 = scr2 * d[i];
	temp = -fabs(ascr-d[i])*scr2;
	/* If the species is already at the minimum, disregard
	   destruction when calculating step size */
	if (y[i] == ymin[i]) temp = 0.0;
	scrtch = max(scr1,max(temp,scrtch));
	// this should be a small number, not infinity
	}
    dt = min(sqreps/scrtch,dtg);
    while(1) {
	/*
	  c the starting values are stored.
	*/
	ts = tn;
	for(i = 0; i < n; i++) {
	    rtau[i] = dt*d[i]/y[i];
	    ys[i] = y[i];
	    qs[i] = q[i];
	    rtaus[i] = rtau[i];
	    }

	/*
	 * find the predictor terms.
	 */
     int bLoop = 1;
     while (bLoop) {
     //restart: // doesnt work in CUDA
	for(i = 0; i < n; i++) {
	    /*
	     * prediction
	     */
	    double rtaui = rtau[i];
	    /*
	    c note that one of two approximations for alpha is chosen:
	    c 1) Pade b for all rtaui (see supporting memo report)
	    c or
	    c 2) Pade a for rtaui<=rswitch,
	    c linear approximation for rtaui > rswitch
	    c (again, see supporting NRL memo report (Mott et al., 2000))
	    c
	    c Option 1): Pade b
	    */
	    alpha = (180.+rtaui*(60.+rtaui*(11.+rtaui)))
		/(360.+ rtaui*(60. + rtaui*(12. + rtaui)));
	    /*
	    c Option 2): Pade a or linear
	    c
	    c if(rtaui.le.rswitch) then
	    c      alpha = (840.+rtaui*(140.+rtaui*(20.+rtaui)))
	    c    &         / (1680. + 40. * rtaui*rtaui)
	    c else
	    c    alpha = 1.-1./rtaui
	    c end if
	    */
	    scrarray[i] = (q[i]-d[i])/(1.0 + alpha*rtaui);
	    }

	iter = 1;
	while(iter <= itermax) {
	    for(i = 0; i < n; i++) {
		/*
		C ym2(i) = ym1(i)
		C ym1(i) = y(i)
		*/
		y[i] = max(ys[i] + dt*scrarray[i], ymin[i]);
		}
	    /*	    if(iter == 1) {  Removed from original algorithm
		    so that previous, rather than first, corrector is
		    compared to.  Results in faster integration. */
		/*
		c the first corrector step advances the time (tentatively) and
		c saves the initial predictor value as y1 for the timestep
		check later.
		*/
		tn = ts + dt;
		for(i = 0; i < n; i++)
		    y1[i] = y[i];
		/*		} Close for "if(iter == 1)" above */
	    /*
	      evaluate the derivitives for the corrector.
	    */
	   // s->derivs(tn + tstart, y, q, d, s->Data);
            CudaclDerivs(tn + tstart, y, q, d, &((CudaclDerivsData*)s->Data)[id]);
	    gcount++;
	    eps = 1.0e-10;
	    for(i = 0; i < n; i++) {
		rtaub = .5*(rtaus[i]+dt*d[i]/y[i]);
		/*
		c Same options for calculating alpha as in predictor:
		c
		c Option 1): Pade b
		*/
		alpha = (180.+rtaub*(60.+rtaub*(11.+rtaub)))
		    / (360. + rtaub*(60. + rtaub*(12. + rtaub)));
		/*
		c Option 2): Pade a or linear
		c
		c if(rtaub.le.rswitch)
		c then
		c alpha = (840.+rtaub*(140.+rtaub*(20.+rtaub)))
		c & / (1680. + 40.*rtaub*rtaub)
		c else
		c alpha = 1.- 1./rtaub
		c end if
		*/
		qt = qs[i]*(1. - alpha) + q[i]*alpha;
		pb = rtaub/dt;
		scrarray[i] = (qt - ys[i]*pb) / (1.0 + alpha*rtaub);
		}
	    iter++;
	    }
	/*
	c calculate new f, check for convergence, and limit decreasing
	c functions. the order of the operations in this loop is important.
	*/
	for(i = 0; i < n; i++) {
	    scr2 = max(ys[i] + dt*scrarray[i], 0.0);
	    scr1 = fabs(scr2 - y1[i]);
	    y[i] = max(scr2, ymin[i]);
	    /*
	    C ym2(i) = ymi(i)
	    C yml(i) = y(i)
	    */
	    if(.25*(ys[i] + y[i]) > ymin[i]) {
		scr1 = scr1/y[i];
		eps = max(.5*(scr1+
			      min(fabs(q[i]-d[i])/(q[i]+d[i]+1.0e-30),scr1)),eps);
		}
	    }
	eps = eps*epscl;
	/*
	c check for convergence.
	c
	c The following section is used for the stability check
	C       stab = 0.01
	C if(itermax.ge.3) then
	C       do i=1,n
	C           stab = max(stab, abs(y(i)-yml(i))/
	C       &       (abs(ymi(i)-ym2(i))+1.e-20*y(i)))
	C end do
	C endif
	*/
	if(eps <= epsmax) {
	    /*
	      & .and.stab.le.1.
	    c
	    c Valid step. Return if dtg has been reached.
	    */
	    //printf("%g %g\n", dtg, tn*tfd);
	    //if(1==1){
	    if(dtg <= tn*tfd) {
                for (int i = 0; i < 5; ++i) {
                    _y[id * 5 + i] = y[i];
                }
	        return;
	    }
	    }
	else {
	    /*
	      Invalid step; reset tn to ts
	    */
	    tn = ts;
	    }
	/*
	  perform stepsize modifications.
	  estimate sqrt(eps) by newton iteration.
	*/
	rteps = 0.5*(eps + 1.0);
	rteps = 0.5*(rteps + eps/rteps);
	rteps = 0.5*(rteps + eps/rteps);

	dto = dt;
	dt = min(dt*(1.0/rteps+.005), tfd*(dtg - tn));
	/* & ,dto/(stab+.001) */
	/*
	  begin new step if previous step converged.
	*/
	bLoop = 0;
	if(eps > epsmax) {
	    /*    & .or. stab. gt. 1 */
	    rcount++;
	    /*
	    c After an unsuccessful step the initial timescales don't
	    c change, but dt does, requiring rtaus to be scaled by the
	    c ratio of the new and old timesteps.
	    */
	    dto = dt/dto;
	    for(i = 0; i < n; i++) {
		rtaus[i] = rtaus[i]*dto;
		}
	    /*
	     * Unsuccessful steps return to line 101 so that the initial
	     * source terms do not get recalculated.
	    */
	    //goto restart;
	    bLoop = 1;
	    }
     }
	/*
	  Successful step; get the source terms for the next step
	  and continue back at line 100
	*/
	//s->derivs(tn + tstart, y, q, d, s->Data);
        CudaclDerivs(tn + tstart, y, q, d, &((CudaclDerivsData*)s->Data)[id]);
	gcount++;
	}

    // Copy back the modified data
    for (int i = 0; i < 5; ++i) {
        _y[id * 5 + i] = y[i];
        }
    }
