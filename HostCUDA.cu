#ifdef _WIN32
#define NOMINMAX
#endif

#ifdef HAPI_MEMPOOL
#define GPU_MEMPOOL
#endif
#ifdef CUDA_INSTRUMENT_WRS
#define GPU_INSTRUMENT_WRS
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
// #include <cutil.h>
#include <assert.h>

#include "CudaFunctions.h"
#ifdef HEXADECAPOLE
# include "CUDAMoments.cu"
#endif
#include "HostCUDA.h"
#include "EwaldCUDA.h"

#include "hapi.h"


#ifdef CUDA_TRACE
#  define CUDA_TRACE_BEGIN()   double trace_start_time = CmiWallTimer() 
#  define CUDA_TRACE_END(ID)   traceUserBracketEvent(ID, trace_start_time, CmiWallTimer()) 
#else
#  define CUDA_TRACE_BEGIN() /* */
#  define CUDA_TRACE_END(ID) /* */
#endif


#define CmiMyPe() 0
extern int _Cmi_mype;

__device__ __constant__ EwaldReadOnlyData cachedData[1];
__device__ __constant__ EwtData ewt[NEWH];  


//__constant__ constantData[88];
//
//
#ifdef CUDA_TRACE
extern "C" void traceUserBracketEvent(int e, double beginT, double endT);
extern "C" double CmiWallTimer();
#endif


void allocatePinnedHostMemory(void **ptr, int size){
  if(size <= 0){
    *((char **)ptr) = NULL;
#ifdef CUDA_PRINT_ERRORS
    printf("allocatePinnedHostMemory: 0 size!\n");
#endif
    exit(-1);
    return;
  }
#ifdef HAPI_MEMPOOL
  *ptr = hapiPoolMalloc(size);
#else
  cudaMallocHost(ptr, size);
#endif
#ifdef CUDA_PRINT_ERRORS
  printf("allocatePinnedHostMemory: %s size: %d\n", cudaGetErrorString( cudaGetLastError() ), size);
#endif
}

void freePinnedHostMemory(void *ptr){
  if(ptr == NULL){
#ifdef CUDA_PRINT_ERRORS
    printf("freePinnedHostMemory: NULL ptr!\n");
#endif
    exit(-1);
    return;
  }
  hapiHostFree(ptr);
#ifdef CUDA_PRINT_ERRORS
  printf("freePinnedHostMemory: %s\n", cudaGetErrorString( cudaGetLastError() ));
#endif
}


/******************* Transfers ******************/

void run_DM_TRANSFER_LOCAL(workRequest *wr, cudaStream_t kernel_stream,void** devBuffers) {
#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
  printf("DM_TRANSFER_LOCAL KERNELSELECT\n");
  printf("mom: 0x%x\n", devBuffers[LOCAL_MOMENTS]);
  printf("cores: 0x%x\n", devBuffers[LOCAL_PARTICLE_CORES]);
  printf("vars: 0x%x\n", devBuffers[LOCAL_PARTICLE_VARS]);
#endif
  if( wr->buffers[LOCAL_MOMENTS_IDX].transfer_to_device ){
    hapiHostFree(wr->buffers[LOCAL_MOMENTS_IDX].host_buffer);
  }
#ifdef CUDA_PRINT_ERRORS
  printf("DM_TRANSFER_LOCAL 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
  if( wr->buffers[LOCAL_PARTICLE_CORES_IDX].transfer_to_device ){
	hapiHostFree(wr->buffers[LOCAL_PARTICLE_CORES_IDX].host_buffer);
  }
#ifdef CUDA_PRINT_ERRORS
  printf("DM_TRANSFER_LOCAL 1: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
  if( wr->buffers[LOCAL_PARTICLE_VARS_IDX].transfer_to_device ){
	hapiHostFree(wr->buffers[LOCAL_PARTICLE_VARS_IDX].host_buffer);
  }
#ifdef CUDA_PRINT_ERRORS
  printf("DM_TRANSFER_LOCAL 2: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
}


void run_DM_TRANSFER_REMOTE_CHUNK(workRequest *wr, cudaStream_t kernel_stream,void** devBuffers) {
#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
  printf("DM_TRANSFER_REMOTE_CHUNK, %d KERNELSELECT\n", wr->buffers[REMOTE_MOMENTS_IDX].transfer_to_device);
#endif
  if( wr->buffers[REMOTE_MOMENTS_IDX].transfer_to_device ){
    hapiHostFree(wr->buffers[REMOTE_MOMENTS_IDX].host_buffer);
#ifdef CUDA_PRINT_ERRORS
    printf("DM_TRANSFER_REMOTE_CHUNK 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
  }

  if( wr->buffers[REMOTE_PARTICLE_CORES_IDX].transfer_to_device ){
    hapiHostFree(wr->buffers[REMOTE_PARTICLE_CORES_IDX].host_buffer);
#ifdef CUDA_PRINT_ERRORS
    printf("DM_TRANSFER_REMOTE_CHUNK 1: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
  }
}


void run_DM_TRANSFER_FREE_LOCAL(workRequest *wr, cudaStream_t kernel_stream,void** devBuffers) {
#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
  printf("DM_TRANSFER_FREE_LOCAL KERNELSELECT\n");
  printf("buffers:\nlocal_particles: (0x%x)\nlocal_particle_vars: (0x%x)\nlocal_moments: (0x%x)\n", 
      devBuffers[LOCAL_PARTICLE_CORES],
      devBuffers[LOCAL_PARTICLE_VARS],
      devBuffers[LOCAL_MOMENTS]
      );

#endif
}





/// @brief queue work request to tranfer local moments and particles to GPU
/// @param moments array of moments
/// @param nMoments
/// @param compactParts  Array of particles
/// @param nCompactParts
/// @param mype Only used for debugging
#ifdef CUDA_INSTRUMENT_WRS
void DataManagerTransferLocalTree(CudaMultipoleMoments *moments, int nMoments,
                        CompactPartData *compactParts, int nCompactParts,
                        int mype, char phase, void *wrCallback) {
#else
void DataManagerTransferLocalTree(CudaMultipoleMoments *moments, int nMoments,
                        CompactPartData *compactParts, int nCompactParts,
                        int mype, void *wrCallback) {
#endif

	workRequest* transferKernel = hapiCreateWorkRequest();
	void* bufHostBuffer;
        int size;

        size = (nMoments) * sizeof(CudaMultipoleMoments);

        //double mill = 1e6;
        //printf("(%d) DM local moments: %f mbytes\n", mype, 1.0*size/mill);

#ifdef CUDA_PRINT_ERRORS
        printf("(%d) DMLocal 0000: %s\n", mype, cudaGetErrorString( cudaGetLastError() ));
#endif
        if(size > 0){
#ifdef HAPI_USE_CUDAMALLOCHOST
#ifdef HAPI_MEMPOOL
          bufHostBuffer = hapiPoolMalloc(size);
#else
          cudaMallocHost(&bufHostBuffer, size);
#endif
#else
          bufHostBuffer = malloc(size);
#endif
        }
        else{
          bufHostBuffer = NULL;
        }

#ifdef CUDA_PRINT_ERRORS
        printf("(%d) DMLocal 0: %s hostBuf: 0x%x, size: %d\n", mype, cudaGetErrorString( cudaGetLastError() ), bufHostBuffer, size );
#endif
        memcpy(bufHostBuffer, moments, size);

        transferKernel->addBuffer(bufHostBuffer, size, (size > 0), false,
                                      false, LOCAL_MOMENTS);

        size = (nCompactParts)*sizeof(CompactPartData);

        if(size > 0){
#ifdef HAPI_USE_CUDAMALLOCHOST
#ifdef HAPI_MEMPOOL
          bufHostBuffer = hapiPoolMalloc(size);
#else
          cudaMallocHost(&bufHostBuffer, size);
#endif
#else
          bufHostBuffer = malloc(size);
#endif
#ifdef CUDA_PRINT_ERRORS
          printf("(%d) DMLocal 1: %s\n", mype, cudaGetErrorString( cudaGetLastError() ) );
#endif
          memcpy(bufHostBuffer, compactParts, size);
        }
        else{
          bufHostBuffer = NULL;
        }

        transferKernel->addBuffer(bufHostBuffer, size, (size > 0), false,
                                      false, LOCAL_PARTICLE_CORES);

        VariablePartData *zeroArray;

        size = (nCompactParts)*sizeof(VariablePartData);

        if(size > 0){
#ifdef HAPI_USE_CUDAMALLOCHOST
#ifdef HAPI_MEMPOOL
          bufHostBuffer = hapiPoolMalloc(size);
#else
          cudaMallocHost(&bufHostBuffer, size);
#endif
#else
          bufHostBuffer = malloc(size);
#endif
#ifdef CUDA_PRINT_ERRORS
          printf("(%d) DMLocal 2: %s\n", mype, cudaGetErrorString( cudaGetLastError() ));
#endif
          zeroArray = (VariablePartData *) bufHostBuffer;
          for(int i = 0; i < nCompactParts; i++){
            zeroArray[i].a.x = 0.0;
            zeroArray[i].a.y = 0.0;
            zeroArray[i].a.z = 0.0;
            zeroArray[i].potential = 0.0;
            zeroArray[i].dtGrav = 0.0;
          }
        }
        else{
          bufHostBuffer = NULL;
        }

	transferKernel->addBuffer(bufHostBuffer, size, (size > 0), false, false, LOCAL_PARTICLE_VARS);

	transferKernel->setDeviceToHostCallback(wrCallback);
	transferKernel->setTraceName("xferLocal");
	transferKernel->setRunKernel(run_DM_TRANSFER_LOCAL);
#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
        printf("(%d) DM LOCAL TREE moments %d (%d) partcores %d (%d) partvars %d (%d)\n",
                  CmiMyPe(),
                  transferKernel->buffers[LOCAL_MOMENTS_IDX].size,
                  transferKernel->buffers[LOCAL_MOMENTS_IDX].transfer_to_device,
                  transferKernel->buffers[LOCAL_PARTICLE_CORES_IDX].size,
                  transferKernel->buffers[LOCAL_PARTICLE_CORES_IDX].transfer_to_device,
                  transferKernel->buffers[LOCAL_PARTICLE_VARS_IDX].size,
                  transferKernel->buffers[LOCAL_PARTICLE_VARS_IDX].transfer_to_device
                  );
#endif
#ifdef CUDA_INSTRUMENT_WRS
        transferKernel->chareIndex = mype;
        transferKernel->compType = DM_TRANSFER_LOCAL;
        transferKernel->compPhase = phase;
#endif
	hapiEnqueue(transferKernel);

}

#ifdef CUDA_INSTRUMENT_WRS
  void DataManagerTransferRemoteChunk(CudaMultipoleMoments *moments, int nMoments, CompactPartData *compactParts, int nCompactParts, int mype, char phase, void *wrCallback) {
#else
    void DataManagerTransferRemoteChunk(CudaMultipoleMoments *moments, int nMoments, CompactPartData *compactParts, int nCompactParts, void *wrCallback) {
#endif

  workRequest* transferKernel = hapiCreateWorkRequest();
  void* bufHostBuffer;
  int size;

  size = (nMoments) * sizeof(CudaMultipoleMoments);
  
  //double mill = 1e6;
  //printf("DM remote moments: %f mbytes\n", 1.0*size/mill);
  
  if(size > 0){
    CUDA_MALLOC(bufHostBuffer, size);
#ifdef CUDA_PRINT_ERRORS
    printf("DMRemote 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
    memcpy(bufHostBuffer, moments, size);
  }

  transferKernel->addBuffer(bufHostBuffer, size, (size > 0), false, false,
                                REMOTE_MOMENTS);

  size = (nCompactParts)*sizeof(CompactPartData);  
  //printf("DM remote cores: %f mbytes\n", 1.0*size/mill);

  if(size > 0){
    CUDA_MALLOC(bufHostBuffer, size);
#ifdef CUDA_PRINT_ERRORS
    printf("DMRemote 1: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
    memcpy(bufHostBuffer, compactParts, size);
  }

  transferKernel->addBuffer(bufHostBuffer, size, (size > 0), false, false,
                                REMOTE_PARTICLE_CORES);

#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
  printf("(%d) DM REMOTE CHUNK moments %d (%d) partcores %d (%d)\n",
            CmiMyPe(),
            transferKernel->buffers[REMOTE_MOMENTS_IDX].size,
            transferKernel->buffers[REMOTE_MOMENTS_IDX].transfer_to_device,
            transferKernel->buffers[REMOTE_PARTICLE_CORES_IDX].size,
            transferKernel->buffers[REMOTE_PARTICLE_CORES_IDX].transfer_to_device
            );
#endif

  transferKernel->setDeviceToHostCallback(wrCallback);
  transferKernel->setTraceName("xferRemote");
  transferKernel->setRunKernel(run_DM_TRANSFER_REMOTE_CHUNK);
#ifdef CUDA_INSTRUMENT_WRS
  transferKernel->chareIndex = mype;
  transferKernel->compType = DM_TRANSFER_REMOTE_CHUNK;
  transferKernel->compPhase = phase;
#endif
  hapiEnqueue(transferKernel);

}


/************** Gravity *****************/

void run_TP_GRAVITY_LOCAL(workRequest *wr, cudaStream_t kernel_stream,void** devBuffers) {
  ParameterStruct *ptr = (ParameterStruct *)wr->getUserData();
#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
  printf("TP_GRAVITY_LOCAL KERNELSELECT buffers:\nlocal_particles: (0x%x)\nlocal_particle_vars: (0x%x)\nlocal_moments: (0x%x)\nil_cell: %d (0x%x)\n", 
      devBuffers[LOCAL_PARTICLE_CORES],
      devBuffers[LOCAL_PARTICLE_VARS],
      devBuffers[LOCAL_MOMENTS],
      wr->buffers[ILCELL_IDX].id,
      devBuffers[wr->buffers[ILCELL_IDX].id]
      );
#endif
  if( wr->buffers[ILCELL_IDX].transfer_to_device ){
#ifndef CUDA_NO_KERNELS
    nodeGravityComputation<<<wr->grid_dim, wr->block_dim, wr->shared_mem, kernel_stream>>>
      (
       (CompactPartData *)devBuffers[LOCAL_PARTICLE_CORES],
       (VariablePartData *)devBuffers[LOCAL_PARTICLE_VARS],
       (CudaMultipoleMoments *)devBuffers[LOCAL_MOMENTS],
       (ILCell *)devBuffers[wr->buffers[ILCELL_IDX].id],
       (int *)devBuffers[wr->buffers[NODE_BUCKET_MARKERS_IDX].id],
       (int *)devBuffers[wr->buffers[NODE_BUCKET_START_MARKERS_IDX].id],
       (int *)devBuffers[wr->buffers[NODE_BUCKET_SIZES_IDX].id],
       ptr->fperiod
      );
#endif
    CUDA_TRACE_BEGIN();
    hapiHostFree((ILCell *)wr->buffers[ILCELL_IDX].host_buffer);
    hapiHostFree((int *)wr->buffers[NODE_BUCKET_MARKERS_IDX].host_buffer);
    hapiHostFree((int *)wr->buffers[NODE_BUCKET_START_MARKERS_IDX].host_buffer);
    hapiHostFree((int *)wr->buffers[NODE_BUCKET_SIZES_IDX].host_buffer);

#ifdef CUDA_PRINT_ERRORS
    printf("TP_GRAVITY_LOCAL 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
    CUDA_TRACE_END(CUDA_LOCAL_NODE_KERNEL);

  }
}

void run_TP_PART_GRAVITY_LOCAL_SMALLPHASE(workRequest *wr, cudaStream_t kernel_stream,void** devBuffers) {
  ParameterStruct *ptr = (ParameterStruct *)wr->getUserData();
#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
  printf("TP_PART_GRAVITY_LOCAL_SMALLPHASE KERNELSELECT buffers:\nlocal_particles: (0x%x)\nlocal_particle_vars: (0x%x)\nil_cell: %d (0x%x)\n", 
      devBuffers[LOCAL_PARTICLE_CORES],
      devBuffers[LOCAL_PARTICLE_VARS],
      wr->buffers[ILPART_IDX].id,
      devBuffers[wr->buffers[ILPART_IDX].id]
      );
#endif
  if( wr->buffers[ILPART_IDX].transfer_to_device ){
#ifndef CUDA_NO_KERNELS
#ifdef CUDA_2D_TB_KERNEL
    particleGravityComputation<<<wr->grid_dim, wr->block_dim, wr->shared_mem, kernel_stream>>>
      (
       (CompactPartData *)devBuffers[LOCAL_PARTICLE_CORES],
       (VariablePartData *)devBuffers[LOCAL_PARTICLE_VARS],
       (CompactPartData *)devBuffers[wr->buffers[MISSED_PARTS_IDX].id],
       (ILCell *)devBuffers[wr->buffers[ILPART_IDX].id],
       (int *)devBuffers[wr->buffers[PART_BUCKET_MARKERS_IDX].id],
       (int *)devBuffers[wr->buffers[PART_BUCKET_START_MARKERS_IDX].id],
       (int *)devBuffers[wr->buffers[PART_BUCKET_SIZES_IDX].id],
       ptr->fperiod 
      );
#else
    particleGravityComputation<<<wr->grid_dim, wr->block_dim, wr->shared_mem, kernel_stream>>>
      (
       (CompactPartData *)devBuffers[LOCAL_PARTICLE_CORES],
       (VariablePartData *)devBuffers[LOCAL_PARTICLE_VARS],
       (CompactPartData *)devBuffers[wr->buffers[MISSED_PARTS_IDX].id],
       (ILPart *)devBuffers[wr->buffers[ILPART_IDX].id],
       (int *)devBuffers[wr->buffers[PART_BUCKET_MARKERS_IDX].id],
       (int *)devBuffers[wr->buffers[PART_BUCKET_START_MARKERS_IDX].id],
       (int *)devBuffers[wr->buffers[PART_BUCKET_SIZES_IDX].id],
       ptr->fperiod 
      );
#endif
#endif


    CUDA_TRACE_BEGIN();
    hapiHostFree((CompactPartData *)wr->buffers[MISSED_PARTS_IDX].host_buffer);
    hapiHostFree(wr->buffers[ILPART_IDX].host_buffer);
    hapiHostFree((int *)wr->buffers[PART_BUCKET_MARKERS_IDX].host_buffer);
    hapiHostFree((int *)wr->buffers[PART_BUCKET_START_MARKERS_IDX].host_buffer);
    hapiHostFree((int *)wr->buffers[PART_BUCKET_SIZES_IDX].host_buffer);

#ifdef CUDA_PRINT_ERRORS
    printf("TP_PART_GRAVITY_LOCAL_SMALLPHASE 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
    CUDA_TRACE_END(CUDA_LOCAL_PART_KERNEL);
  }
}

void run_TP_PART_GRAVITY_LOCAL(workRequest *wr, cudaStream_t kernel_stream,void** devBuffers) {
  ParameterStruct *ptr = (ParameterStruct *)wr->getUserData();
#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
  printf("TP_PART_GRAVITY_LOCAL KERNELSELECT buffers:\nlocal_particles: (0x%x)\nlocal_particle_vars: (0x%x)\nil_cell: %d (0x%x)\n", 
      devBuffers[LOCAL_PARTICLE_CORES],
      devBuffers[LOCAL_PARTICLE_VARS],
      wr->buffers[ILPART_IDX].id,
      devBuffers[wr->buffers[ILPART_IDX].id]
      );
#endif
  if( wr->buffers[ILPART_IDX].transfer_to_device ){
#ifndef CUDA_NO_KERNELS
#ifdef CUDA_2D_TB_KERNEL
    particleGravityComputation<<<wr->grid_dim, wr->block_dim, wr->shared_mem, kernel_stream>>>
      (
       (CompactPartData *)devBuffers[LOCAL_PARTICLE_CORES],
       (VariablePartData *)devBuffers[LOCAL_PARTICLE_VARS],
       (CompactPartData *)devBuffers[LOCAL_PARTICLE_CORES],
       (ILCell *)devBuffers[wr->buffers[ILPART_IDX].id],
       (int *)devBuffers[wr->buffers[PART_BUCKET_MARKERS_IDX].id],
       (int *)devBuffers[wr->buffers[PART_BUCKET_START_MARKERS_IDX].id],
       (int *)devBuffers[wr->buffers[PART_BUCKET_SIZES_IDX].id],
       ptr->fperiod
      );
#else
    particleGravityComputation<<<wr->grid_dim, wr->block_dim, wr->shared_mem, kernel_stream>>>
      (
       (CompactPartData *)devBuffers[LOCAL_PARTICLE_CORES],
       (VariablePartData *)devBuffers[LOCAL_PARTICLE_VARS],
       (CompactPartData *)devBuffers[LOCAL_PARTICLE_CORES],
       (ILPart *)devBuffers[wr->buffers[ILPART_IDX].id],
       (int *)devBuffers[wr->buffers[PART_BUCKET_MARKERS_IDX].id],
       (int *)devBuffers[wr->buffers[PART_BUCKET_START_MARKERS_IDX].id],
       (int *)devBuffers[wr->buffers[PART_BUCKET_SIZES_IDX].id],
       ptr->fperiod
      );
#endif
#endif

    CUDA_TRACE_BEGIN();
    hapiHostFree(wr->buffers[ILPART_IDX].host_buffer);
    hapiHostFree((int *)wr->buffers[PART_BUCKET_MARKERS_IDX].host_buffer);
    hapiHostFree((int *)wr->buffers[PART_BUCKET_START_MARKERS_IDX].host_buffer);
    hapiHostFree((int *)wr->buffers[PART_BUCKET_SIZES_IDX].host_buffer);

#ifdef CUDA_PRINT_ERRORS
    printf("TP_PART_GRAVITY_LOCAL 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
    CUDA_TRACE_END(CUDA_LOCAL_PART_KERNEL);
  }
}

void run_TP_GRAVITY_REMOTE(workRequest *wr, cudaStream_t kernel_stream,void** devBuffers) {
  ParameterStruct *ptr = (ParameterStruct *)wr->getUserData();
#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
  printf("TP_GRAVITY_REMOTE KERNELSELECT buffers:\nlocal_particles: (0x%x)\nlocal_particle_vars: (0x%x)\nremote_moments: (0x%x)\nil_cell: %d (0x%x)\n", 
      devBuffers[LOCAL_PARTICLE_CORES],
      devBuffers[LOCAL_PARTICLE_VARS],
      devBuffers[REMOTE_MOMENTS],
      wr->buffers[ILCELL_IDX].id,
      devBuffers[wr->buffers[ILCELL_IDX].id]
      );
#endif
  if( wr->buffers[ILCELL_IDX].transfer_to_device ){
#ifndef CUDA_NO_KERNELS
    nodeGravityComputation<<<wr->grid_dim, wr->block_dim, wr->shared_mem, kernel_stream>>>
      (
       (CompactPartData *)devBuffers[LOCAL_PARTICLE_CORES],
       (VariablePartData *)devBuffers[LOCAL_PARTICLE_VARS],
       (CudaMultipoleMoments *)devBuffers[REMOTE_MOMENTS],
       (ILCell *)devBuffers[wr->buffers[ILCELL_IDX].id],
       (int *)devBuffers[wr->buffers[NODE_BUCKET_MARKERS_IDX].id],
       (int *)devBuffers[wr->buffers[NODE_BUCKET_START_MARKERS_IDX].id],
       (int *)devBuffers[wr->buffers[NODE_BUCKET_SIZES_IDX].id],
       ptr->fperiod
      );
#endif


    CUDA_TRACE_BEGIN();
    hapiHostFree((ILCell *)wr->buffers[ILCELL_IDX].host_buffer);
    hapiHostFree((int *)wr->buffers[NODE_BUCKET_MARKERS_IDX].host_buffer);
    hapiHostFree((int *)wr->buffers[NODE_BUCKET_START_MARKERS_IDX].host_buffer);
    hapiHostFree((int *)wr->buffers[NODE_BUCKET_SIZES_IDX].host_buffer);
#ifdef CUDA_PRINT_ERRORS
    printf("TP_GRAVITY_REMOTE 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
    CUDA_TRACE_END(CUDA_REMOTE_NODE_KERNEL);
  }
}

void run_TP_PART_GRAVITY_REMOTE(workRequest *wr, cudaStream_t kernel_stream,void** devBuffers) {
  ParameterStruct *ptr = (ParameterStruct *)wr->getUserData();
#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
  printf("TP_PART_GRAVITY_REMOTE KERNELSELECT buffers:\nlocal_particles: (0x%x)\nlocal_particle_vars: (0x%x)\nil_cell: %d (0x%x)\n", 
      devBuffers[LOCAL_PARTICLE_CORES],
      devBuffers[LOCAL_PARTICLE_VARS],
      wr->buffers[ILPART_IDX].id,
      devBuffers[wr->buffers[ILPART_IDX].id]
      );
#endif
  if( wr->buffers[ILPART_IDX].transfer_to_device ){
#ifndef CUDA_NO_KERNELS
#ifdef CUDA_2D_TB_KERNEL
    particleGravityComputation<<<wr->grid_dim, wr->block_dim, wr->shared_mem, kernel_stream>>>
      (
       (CompactPartData *)devBuffers[LOCAL_PARTICLE_CORES],
       (VariablePartData *)devBuffers[LOCAL_PARTICLE_VARS],
       (CompactPartData *)devBuffers[REMOTE_PARTICLE_CORES],
       (ILCell *)devBuffers[wr->buffers[ILPART_IDX].id],
       (int *)devBuffers[wr->buffers[PART_BUCKET_MARKERS_IDX].id],
       (int *)devBuffers[wr->buffers[PART_BUCKET_START_MARKERS_IDX].id],
       (int *)devBuffers[wr->buffers[PART_BUCKET_SIZES_IDX].id],
       ptr->fperiod
      );
#else
    particleGravityComputation<<<wr->grid_dim, wr->block_dim, wr->shared_mem, kernel_stream>>>
      (
       (CompactPartData *)devBuffers[LOCAL_PARTICLE_CORES],
       (VariablePartData *)devBuffers[LOCAL_PARTICLE_VARS],
       (CompactPartData *)devBuffers[REMOTE_PARTICLE_CORES],
       (ILPart *)devBuffers[wr->buffers[ILPART_IDX].id],
       (int *)devBuffers[wr->buffers[PART_BUCKET_MARKERS_IDX].id],
       (int *)devBuffers[wr->buffers[PART_BUCKET_START_MARKERS_IDX].id],
       (int *)devBuffers[wr->buffers[PART_BUCKET_SIZES_IDX].id],
       ptr->fperiod
      );
#endif
#endif


    CUDA_TRACE_BEGIN();
    hapiHostFree(wr->buffers[ILPART_IDX].host_buffer);
    hapiHostFree((int *)wr->buffers[PART_BUCKET_MARKERS_IDX].host_buffer);
    hapiHostFree((int *)wr->buffers[PART_BUCKET_START_MARKERS_IDX].host_buffer);
    hapiHostFree((int *)wr->buffers[PART_BUCKET_SIZES_IDX].host_buffer);

#ifdef CUDA_PRINT_ERRORS
    printf("TP_PART_GRAVITY_REMOTE 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
    CUDA_TRACE_END(CUDA_REMOTE_PART_KERNEL);
  }
}

void run_TP_GRAVITY_REMOTE_RESUME(workRequest *wr, cudaStream_t kernel_stream,void** devBuffers) {
  ParameterStruct *ptr = (ParameterStruct *)wr->getUserData();
#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
  printf("TP_GRAVITY_REMOTE_RESUME KERNELSELECT buffers:\nlocal_particles: (0x%x)\nlocal_particle_vars: (0x%x)\nmissed_moments: %d (0x%x)\nil_cell: %d (0x%x)\n", 
      devBuffers[LOCAL_PARTICLE_CORES],
      devBuffers[LOCAL_PARTICLE_VARS],
      wr->buffers[MISSED_MOMENTS_IDX].id,
      devBuffers[wr->buffers[MISSED_MOMENTS_IDX].id],
      wr->buffers[ILCELL_IDX].id,
      devBuffers[wr->buffers[ILCELL_IDX].id]
      );
#endif
  if( wr->buffers[ILCELL_IDX].transfer_to_device ){
#ifndef CUDA_NO_KERNELS
    nodeGravityComputation<<<wr->grid_dim, wr->block_dim, wr->shared_mem, kernel_stream>>>
      (
       (CompactPartData *)devBuffers[LOCAL_PARTICLE_CORES],
       (VariablePartData *)devBuffers[LOCAL_PARTICLE_VARS],
       (CudaMultipoleMoments *)devBuffers[wr->buffers[MISSED_MOMENTS_IDX].id],
       (ILCell *)devBuffers[wr->buffers[ILCELL_IDX].id],
       (int *)devBuffers[wr->buffers[NODE_BUCKET_MARKERS_IDX].id],
       (int *)devBuffers[wr->buffers[NODE_BUCKET_START_MARKERS_IDX].id],
       (int *)devBuffers[wr->buffers[NODE_BUCKET_SIZES_IDX].id],
       ptr->fperiod
      );
#endif


    CUDA_TRACE_BEGIN();
    hapiHostFree((CudaMultipoleMoments *)wr->buffers[MISSED_MOMENTS_IDX].host_buffer);
    hapiHostFree((ILCell *)wr->buffers[ILCELL_IDX].host_buffer);
    hapiHostFree((int *)wr->buffers[NODE_BUCKET_MARKERS_IDX].host_buffer);
    hapiHostFree((int *)wr->buffers[NODE_BUCKET_START_MARKERS_IDX].host_buffer);
    hapiHostFree((int *)wr->buffers[NODE_BUCKET_SIZES_IDX].host_buffer);
#ifdef CUDA_PRINT_ERRORS
    printf("TP_GRAVITY_REMOTE_RESUME 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
    CUDA_TRACE_END(CUDA_REMOTE_RESUME_NODE_KERNEL);
  }
}

void run_TP_PART_GRAVITY_REMOTE_RESUME(workRequest *wr, cudaStream_t kernel_stream,void** devBuffers) {
  ParameterStruct *ptr = (ParameterStruct *)wr->getUserData();
#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
  printf("TP_PART_GRAVITY_REMOTE_RESUME KERNELSELECT buffers:\nlocal_particles: (0x%x)\nlocal_particle_vars: (0x%x)\nmissed_parts: %d (0x%x)\nil_cell: %d (0x%x)\n", 
      devBuffers[LOCAL_PARTICLE_CORES],
      devBuffers[LOCAL_PARTICLE_VARS],
      wr->buffers[MISSED_PARTS_IDX].id,
      devBuffers[wr->buffers[MISSED_PARTS_IDX].id],
      wr->buffers[ILPART_IDX].id,
      devBuffers[wr->buffers[ILPART_IDX].id]
      );
#endif
  if( wr->buffers[ILPART_IDX].transfer_to_device ){
#ifndef CUDA_NO_KERNELS
#ifdef CUDA_2D_TB_KERNEL
    particleGravityComputation<<<wr->grid_dim, wr->block_dim, wr->shared_mem, kernel_stream>>>
      (
       (CompactPartData *)devBuffers[LOCAL_PARTICLE_CORES],
       (VariablePartData *)devBuffers[LOCAL_PARTICLE_VARS],
       (CompactPartData *)devBuffers[wr->buffers[MISSED_PARTS_IDX].id],
       (ILCell *)devBuffers[wr->buffers[ILPART_IDX].id],
       (int *)devBuffers[wr->buffers[PART_BUCKET_MARKERS_IDX].id],
       (int *)devBuffers[wr->buffers[PART_BUCKET_START_MARKERS_IDX].id],
       (int *)devBuffers[wr->buffers[PART_BUCKET_SIZES_IDX].id],
       ptr->fperiod
      );
#else
    particleGravityComputation<<<wr->grid_dim, wr->block_dim, wr->shared_mem, kernel_stream>>>
      (
       (CompactPartData *)devBuffers[LOCAL_PARTICLE_CORES],
       (VariablePartData *)devBuffers[LOCAL_PARTICLE_VARS],
       (CompactPartData *)devBuffers[wr->buffers[MISSED_PARTS_IDX].id],
       (ILPart *)devBuffers[wr->buffers[ILPART_IDX].id],
       (int *)devBuffers[wr->buffers[PART_BUCKET_MARKERS_IDX].id],
       (int *)devBuffers[wr->buffers[PART_BUCKET_START_MARKERS_IDX].id],
       (int *)devBuffers[wr->buffers[PART_BUCKET_SIZES_IDX].id],
       ptr->fperiod
      );
#endif
#endif


    CUDA_TRACE_BEGIN();
    hapiHostFree((CompactPartData *)wr->buffers[MISSED_PARTS_IDX].host_buffer);
    hapiHostFree(wr->buffers[ILPART_IDX].host_buffer);
    hapiHostFree((int *)wr->buffers[PART_BUCKET_MARKERS_IDX].host_buffer);
    hapiHostFree((int *)wr->buffers[PART_BUCKET_START_MARKERS_IDX].host_buffer);
    hapiHostFree((int *)wr->buffers[PART_BUCKET_SIZES_IDX].host_buffer);
#ifdef CUDA_PRINT_ERRORS
    printf("TP_PART_GRAVITY_REMOTE_RESUME 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
    CUDA_TRACE_END(CUDA_REMOTE_RESUME_PART_KERNEL);
  }
}



void TreePieceCellListDataTransferLocal(CudaRequest *data){
	int numBlocks = data->numBucketsPlusOne-1;

	workRequest* gravityKernel = hapiCreateWorkRequest();

#ifdef CUDA_2D_TB_KERNEL
	gravityKernel->setExecParams(numBlocks, dim3(NODES_PER_BLOCK, PARTS_PER_BLOCK));
#else
	gravityKernel->setexecParams(numBlocks, THREADS_PER_BLOCK);
#endif

	TreePieceCellListDataTransferBasic(data, gravityKernel);
#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
        printf("(%d) TRANSFER LOCAL CELL\n", CmiMyPe());
#endif

	gravityKernel->setDeviceToHostCallback(data->cb);
	gravityKernel->setTraceName("gravityLocal");
	gravityKernel->setRunKernel(run_TP_GRAVITY_LOCAL);
#ifdef CUDA_INSTRUMENT_WRS
        gravityKernel->chareIndex = data->tpIndex;
        gravityKernel->compType = TP_GRAVITY_LOCAL;
        gravityKernel->compPhase = data->phase;
#endif
	hapiEnqueue(gravityKernel);
}

void TreePieceCellListDataTransferRemote(CudaRequest *data){
	int numBlocks = data->numBucketsPlusOne-1;

	workRequest* gravityKernel = hapiCreateWorkRequest();

#ifdef CUDA_2D_TB_KERNEL
	gravityKernel->setExecParams(numBlocks, dim3(NODES_PER_BLOCK, PARTS_PER_BLOCK));
#else
	gravityKernel->setExecParams(numBlocks, THREADS_PER_BLOCK);
#endif

	TreePieceCellListDataTransferBasic(data, gravityKernel);
#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
        printf("(%d) TRANSFER REMOTE CELL\n", CmiMyPe());
#endif

	gravityKernel->setDeviceToHostCallback(data->cb);
	gravityKernel->setTraceName("gravityRemote");
	gravityKernel->setRunKernel(run_TP_GRAVITY_REMOTE);
#ifdef CUDA_INSTRUMENT_WRS
        gravityKernel->chareIndex = data->tpIndex;
        gravityKernel->compType = TP_GRAVITY_REMOTE;
        gravityKernel->compPhase = data->phase;
#endif
	hapiEnqueue(gravityKernel);
}

void TreePieceCellListDataTransferRemoteResume(CudaRequest *data, CudaMultipoleMoments *missedMoments, int numMissedMoments){
  int numBlocks = data->numBucketsPlusOne-1;
  int size;

  workRequest* gravityKernel = hapiCreateWorkRequest();
  void* bufferHostBuffer;
  bool transfer;

#ifdef CUDA_2D_TB_KERNEL
  gravityKernel->setExecParams(numBlocks, dim3(NODES_PER_BLOCK, PARTS_PER_BLOCK));
#else
  gravityKernel->setExecParams(numBlocks, THREADS_PER_BLOCK);
#endif

  TreePieceCellListDataTransferBasic(data, gravityKernel);
  transfer = gravityKernel->buffers[ILCELL_IDX].transfer_to_device;

  size = (numMissedMoments) * sizeof(CudaMultipoleMoments);

  if(transfer){
    CUDA_MALLOC(bufferHostBuffer, size);
#ifdef CUDA_PRINT_ERRORS
    printf("TPRR 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
    memcpy(bufferHostBuffer, missedMoments, size);
  }
#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
  printf("(%d) TRANSFER REMOTE RESUME CELL %d (%d)\n", CmiMyPe(),
        size, transfer);
#endif

  gravityKernel->addBuffer(bufferHostBuffer, size, transfer, false, transfer);

  ParameterStruct *ptr = (ParameterStruct *)gravityKernel->getUserData();
  ptr->numEntities = numMissedMoments;

  gravityKernel->setDeviceToHostCallback(data->cb);
  gravityKernel->setTraceName("remoteResume");
  gravityKernel->setRunKernel(run_TP_GRAVITY_REMOTE_RESUME);
#ifdef CUDA_INSTRUMENT_WRS
        gravityKernel->chareIndex = data->tpIndex;
        gravityKernel->compType = TP_GRAVITY_REMOTE_RESUME;
        gravityKernel->compPhase = data->phase;
#endif
  hapiEnqueue(gravityKernel);
}

void TreePieceCellListDataTransferBasic(CudaRequest *data, workRequest *gravityKernel){
	int numBucketsPlusOne = data->numBucketsPlusOne;
        int numBuckets = numBucketsPlusOne-1;
        int size = (data->numInteractions) * sizeof(ILCell);
        bool transfer = size > 0;

	gravityKernel->addBuffer(data->list, size, transfer, false, transfer);

        size = (numBucketsPlusOne) * sizeof(int);
	gravityKernel->addBuffer(data->bucketMarkers, (transfer ? size : 0),
	                             transfer, false, transfer);

	size = (numBuckets) * sizeof(int);
	gravityKernel->addBuffer(data->bucketStarts, size, transfer, false, transfer);

	gravityKernel->addBuffer(data->bucketSizes, size, transfer, false, transfer);

        ParameterStruct *ptr = (ParameterStruct *)malloc(sizeof(ParameterStruct));
        ptr->numInteractions = data->numInteractions;
        ptr->numBucketsPlusOne = numBucketsPlusOne;
        ptr->fperiod = data->fperiod;

        gravityKernel->setUserData(ptr, sizeof(*ptr));

        free(ptr);

#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
        printf("(%d) TRANSFER BASIC cells %d (%d) bucket_markers %d (%d) bucket_starts %d (%d) bucket_sizes %d (%d)\n",
            CmiMyPe(),
            gravityKernel->buffers[ILCELL_IDX].size,
            gravityKernel->buffers[ILCELL_IDX].transfer_to_device,
            gravityKernel->buffers[NODE_BUCKET_MARKERS_IDX].size,
            gravityKernel->buffers[NODE_BUCKET_MARKERS_IDX].transfer_to_device,
            gravityKernel->buffers[NODE_BUCKET_START_MARKERS_IDX].size,
            gravityKernel->buffers[NODE_BUCKET_START_MARKERS_IDX].transfer_to_device,
            gravityKernel->buffers[NODE_BUCKET_SIZES_IDX].size,
            gravityKernel->buffers[NODE_BUCKET_SIZES_IDX].transfer_to_device
            );
#endif
}

void TreePiecePartListDataTransferLocalSmallPhase(CudaRequest *data, CompactPartData *particles, int len){
	int numBlocks = data->numBucketsPlusOne-1;

	workRequest* gravityKernel = hapiCreateWorkRequest();
	void* bufferHostBuffer;
        int size;
        ParameterStruct *ptr;
        bool transfer;

#ifdef CUDA_2D_TB_KERNEL
	gravityKernel->setExecParams(numBlocks, dim3(NODES_PER_BLOCK_PART, PARTS_PER_BLOCK_PART));
#else
	gravityKernel->setExecParams(numBlocks, THREADS_PER_BLOCK);
#endif

	TreePiecePartListDataTransferBasic(data, gravityKernel);
        transfer = gravityKernel->buffers[ILPART_IDX].transfer_to_device;

        size = (len) * sizeof(CompactPartData);

#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
        printf("(%d) TRANSFER LOCAL SMALL PHASE  %d (%d)\n",
            CmiMyPe(),
            size,
            transfer
            );
#endif

        if(transfer){
          CUDA_MALLOC(bufferHostBuffer, size);
#ifdef CUDA_PRINT_ERRORS
          printf("TPPartSmallPhase 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
          memcpy(bufferHostBuffer, particles, size);
        }
        gravityKernel->addBuffer(bufferHostBuffer, size, transfer, false, transfer);

        ptr = (ParameterStruct *)gravityKernel->getUserData();
        ptr->numMissedCores = len;

	gravityKernel->setDeviceToHostCallback(data->cb);
	gravityKernel->setTraceName("partGravityLocal");
	gravityKernel->setRunKernel(run_TP_PART_GRAVITY_LOCAL_SMALLPHASE);
#ifdef CUDA_INSTRUMENT_WRS
        gravityKernel->chareIndex = data->tpIndex;
        gravityKernel->compType = TP_PART_GRAVITY_LOCAL_SMALLPHASE;
        gravityKernel->compPhase = data->phase;
#endif
	hapiEnqueue(gravityKernel);
}

void TreePiecePartListDataTransferLocal(CudaRequest *data){
	int numBlocks = data->numBucketsPlusOne-1;

	workRequest* gravityKernel = hapiCreateWorkRequest();

#ifdef CUDA_2D_TB_KERNEL
	gravityKernel->setExecParams(numBlocks, dim3(NODES_PER_BLOCK_PART, PARTS_PER_BLOCK_PART));
#else
	gravityKernel->setExecParams(numBlocks, THREADS_PER_BLOCK);
#endif

	TreePiecePartListDataTransferBasic(data, gravityKernel);

	gravityKernel->setDeviceToHostCallback(data->cb);
	gravityKernel->setTraceName("partGravityLocal");
	gravityKernel->setRunKernel(run_TP_PART_GRAVITY_LOCAL);
#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
        printf("(%d) TRANSFER LOCAL LARGEPHASE PART\n", CmiMyPe());
#endif
#ifdef CUDA_INSTRUMENT_WRS
        gravityKernel->chareIndex = data->tpIndex;
        gravityKernel->compType = TP_PART_GRAVITY_LOCAL;
        gravityKernel->compPhase = data->phase;
#endif
	hapiEnqueue(gravityKernel);
}

void TreePiecePartListDataTransferRemote(CudaRequest *data){
	int numBlocks = data->numBucketsPlusOne-1;

	workRequest* gravityKernel = hapiCreateWorkRequest();

#ifdef CUDA_2D_TB_KERNEL
	gravityKernel->setExecParams(numBlocks, dim3(NODES_PER_BLOCK_PART, PARTS_PER_BLOCK_PART));
#else
	gravityKernel->setExecParams(numBlocks, THREADS_PER_BLOCK);
#endif

	TreePiecePartListDataTransferBasic(data, gravityKernel);

	gravityKernel->setDeviceToHostCallback(data->cb);
	gravityKernel->setTraceName("partGravityRemote");
	gravityKernel->setRunKernel(run_TP_PART_GRAVITY_REMOTE);
#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
        printf("(%d) TRANSFER REMOTE PART\n", CmiMyPe());
#endif
#ifdef CUDA_INSTRUMENT_WRS
        gravityKernel->chareIndex = data->tpIndex;
        gravityKernel->compType = TP_PART_GRAVITY_REMOTE;
        gravityKernel->compPhase = data->phase;
#endif
	hapiEnqueue(gravityKernel);
}

void TreePiecePartListDataTransferRemoteResume(CudaRequest *data, CompactPartData *missedParts, int numMissedParts){
	int numBlocks = data->numBucketsPlusOne-1;
        int size;
        ParameterStruct *ptr;
        bool transfer;

	workRequest* gravityKernel = hapiCreateWorkRequest();
	void* bufferHostBuffer;

#ifdef CUDA_2D_TB_KERNEL
	gravityKernel->setExecParams(numBlocks, dim3(NODES_PER_BLOCK_PART, PARTS_PER_BLOCK_PART));
#else
	gravityKernel->setExecParams(numBlocks, THREADS_PER_BLOCK);
#endif

	TreePiecePartListDataTransferBasic(data, gravityKernel);
        transfer = gravityKernel->buffers[ILPART_IDX].transfer_to_device;

        size = (numMissedParts) * sizeof(CompactPartData);

#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
        printf("(%d) TRANSFER REMOTE RESUME PART %d (%d)\n",
            CmiMyPe(),
            size,
            transfer
            );
#endif

        if(transfer){
          CUDA_MALLOC(bufferHostBuffer, size);
#ifdef CUDA_PRINT_ERRORS
          printf("TPPartRR 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
          memcpy(bufferHostBuffer, missedParts, size);
        }
        gravityKernel->addBuffer(bufferHostBuffer, size, transfer, false, transfer);

        ptr = (ParameterStruct *)gravityKernel->getUserData();
        ptr->numMissedCores = numMissedParts;

	gravityKernel->setDeviceToHostCallback(data->cb);
	gravityKernel->setTraceName("partGravityRemote");
	gravityKernel->setRunKernel(run_TP_PART_GRAVITY_REMOTE_RESUME);
#ifdef CUDA_INSTRUMENT_WRS
        gravityKernel->chareIndex = data->tpIndex;
        gravityKernel->compType = TP_PART_GRAVITY_REMOTE_RESUME;
        gravityKernel->compPhase = data->phase;
#endif
	hapiEnqueue(gravityKernel);
}

void TreePiecePartListDataTransferBasic(CudaRequest *data, workRequest *gravityKernel){
	int numInteractions = data->numInteractions;
	int numBucketsPlusOne = data->numBucketsPlusOne;
        int numBuckets = numBucketsPlusOne-1;
        int size;
        bool transfer;

	size = (numInteractions) * sizeof(ILCell);
        transfer = size > 0;
	gravityKernel->addBuffer(data->list, size, transfer, false, transfer);

	size = (numBucketsPlusOne) * sizeof(int);
	gravityKernel->addBuffer(data->bucketMarkers, (transfer ? size : 0), transfer, false, transfer);

	size = (numBuckets) * sizeof(int);
	gravityKernel->addBuffer(data->bucketStarts, size, transfer, false, transfer);

	gravityKernel->addBuffer(data->bucketSizes, size, transfer, false, transfer);
        
        ParameterStruct *ptr = (ParameterStruct *)malloc(sizeof(ParameterStruct));
        ptr->numInteractions = data->numInteractions;
        ptr->numBucketsPlusOne = numBucketsPlusOne;
        ptr->fperiod = data->fperiod;
                                       
        gravityKernel->setUserData(ptr, sizeof(*ptr));

        free(ptr);

#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
        printf("(%d) TRANSFER BASIC PART parts %d (%d) bucket_markers %d (%d) bucket_starts %d (%d) bucket_sizes %d (%d)\n",
            CmiMyPe(),
            gravityKernel->buffers[ILPART_IDX].size,
            gravityKernel->buffers[ILPART_IDX].transfer_to_device,
            gravityKernel->buffers[PART_BUCKET_MARKERS_IDX].size,
            gravityKernel->buffers[PART_BUCKET_MARKERS_IDX].transfer_to_device,
            gravityKernel->buffers[PART_BUCKET_START_MARKERS_IDX].size,
            gravityKernel->buffers[PART_BUCKET_START_MARKERS_IDX].transfer_to_device,
            gravityKernel->buffers[PART_BUCKET_SIZES_IDX].size,
            gravityKernel->buffers[PART_BUCKET_SIZES_IDX].transfer_to_device
            );
#endif

}

/* kernels:

TP_GRAVITY_LOCAL,
TP_GRAVITY_REMOTE,
TP_GRAVITY_REMOTE_RESUME,
TP_PART_GRAVITY_LOCAL,
TP_PART_GRAVITY_REMOTE,
TP_PART_GRAVITY_REMOTE_RESUME

 */

#ifdef CUDA_INSTRUMENT_WRS
void FreeDataManagerLocalTreeMemory(bool freemom, bool freepart, int index, char phase){
#else
void FreeDataManagerLocalTreeMemory(bool freemom, bool freepart){
#endif
  workRequest* gravityKernel = hapiCreateWorkRequest();

  gravityKernel->addBuffer(NULL, 0, false, false, freemom, LOCAL_MOMENTS);

  gravityKernel->addBuffer(NULL, 0, false, false, freepart, LOCAL_PARTICLE_CORES);

  gravityKernel->setTraceName("freeLocal");
  gravityKernel->setRunKernel(run_DM_TRANSFER_FREE_LOCAL);
  //printf("DM TRANSFER FREE LOCAL\n");
#ifdef CUDA_INSTRUMENT_WRS
  gravityKernel->chareIndex = index;
  gravityKernel->compType = DM_TRANSFER_FREE_LOCAL;
  gravityKernel->compPhase = phase;
#endif
  hapiEnqueue(gravityKernel);

}

// this function begins the transfer of the next
// pending chunk on dm, once we are assured that
// the previous chunk's allocated memory has been
// freed
void initiateNextChunkTransfer(void *dm);

void run_DM_TRANSFER_FREE_REMOTE_CHUNK(workRequest *wr, cudaStream_t kernel_stream,void** devBuffers) {
#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
  printf("DM_TRANSFER_FREE_REMOTE_CHUNK KERNELSELECT\n");

#endif
  initiateNextChunkTransfer(wr->getUserData());
}


#ifdef CUDA_INSTRUMENT_WRS
void FreeDataManagerRemoteChunkMemory(int chunk, void *dm, bool freemom, bool freepart, int index, char phase){
#else
void FreeDataManagerRemoteChunkMemory(int chunk, void *dm, bool freemom, bool freepart){
#endif
  workRequest* gravityKernel = hapiCreateWorkRequest();

  gravityKernel->addBuffer(NULL, 0, false, false, freemom, REMOTE_MOMENTS);

  gravityKernel->addBuffer(NULL, 0, false, false, freepart, REMOTE_PARTICLE_CORES);

  gravityKernel->setTraceName("freeRemote");
  gravityKernel->setRunKernel(run_DM_TRANSFER_FREE_REMOTE_CHUNK);

  // save a pointer to the data manager so that
  // the next chunk's transfer can be initiated once
  // the memory for this chunk has been freed
  gravityKernel->user_data = static_cast<char*>(dm);
  //printf("DM TRANSFER FREE REMOTE CHUNK\n");
#ifdef CUDA_INSTRUMENT_WRS
  gravityKernel->chareIndex = index;
  gravityKernel->compType = DM_TRANSFER_FREE_REMOTE_CHUNK;
  gravityKernel->compPhase = phase;
#endif
  hapiEnqueue(gravityKernel);

}


void run_DM_TRANSFER_BACK(workRequest *wr, cudaStream_t kernel_stream,void** devBuffers) {
#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
  printf("DM_TRANSFER_BACK: 0x%x KERNELSELECT\n", devBuffers[LOCAL_PARTICLE_VARS]);
#endif
}

/* Schedule the transfer of the accelerations back from the GPU to the host.
 * This also schedules the freeing of the device buffers used for the
 * force calculation.
 */
#ifdef CUDA_INSTRUMENT_WRS
void TransferParticleVarsBack(VariablePartData *hostBuffer, int size, void *cb,
     bool freemom, bool freepart, bool freeRemoteMom, bool, freeRemotePart,
     int index, char phase){
#else
void TransferParticleVarsBack(VariablePartData *hostBuffer, int size, void *cb,
     bool freemom, bool freepart, bool freeRemoteMom, bool freeRemotePart){
#endif
  workRequest* gravityKernel = hapiCreateWorkRequest();

#ifdef CUDA_PRINT_TRANSFER_BACK_PARTICLES
  printf("Enqueue kernel to transfer particles back from: 0x%x\n", 0/*devBuffers[LOCAL_PARTICLE_VARS]*/);
#endif

  /* The buffers are: 1) forces on local particles, 2) Local multipole moments
   * 3) Local particle data, 4) Remote multipole moments and 5) remote
   * particle data.  Buffers 2-5 are here to get their device buffer
   * deallocated.
   */

  /* schedule buffers for transfer to the GPU */
  gravityKernel->addBuffer(hostBuffer, size, false, (size > 0), (size > 0), LOCAL_PARTICLE_VARS);

  gravityKernel->addBuffer(NULL, 0, false, false, freemom, LOCAL_MOMENTS);

  gravityKernel->addBuffer(NULL, 0, false, false, freepart, LOCAL_PARTICLE_CORES);

  gravityKernel->addBuffer(NULL, 0, false, false, freeRemoteMom, REMOTE_MOMENTS);

  gravityKernel->addBuffer(NULL, 0, false, false, freeRemotePart, REMOTE_PARTICLE_CORES);

  gravityKernel->setDeviceToHostCallback(cb);
  gravityKernel->setTraceName("transferBack");
  gravityKernel->setRunKernel(run_DM_TRANSFER_BACK);
#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
  printf("(%d) DM TRANSFER BACK\n", CmiMyPe());
#endif
#ifdef CUDA_INSTRUMENT_WRS
  gravityKernel->chareIndex = index;
  gravityKernel->compType = DM_TRANSFER_BACK;
  gravityKernel->compPhase = phase;
#endif
  hapiEnqueue(gravityKernel);
}

/*
void DummyKernel(void *cb){
  workRequest* dummy = hapiCreateWorkRequest();

  dummy->setExecParams(1, THREADS_PER_BLOCK);

  dummy->setDeviceToHostCallback(cb);
  dummy->setTraceName("dummyRun");
  dummy->setRunKernel(run_kernel_DUMMY);
  hapiEnqueue(dummy);

}
*/

/*
 * Kernels
 */

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
        sumx += __shfl_down(sumx, offset, NODES_PER_BLOCK);
        sumy += __shfl_down(sumy, offset, NODES_PER_BLOCK);
        sumz += __shfl_down(sumz, offset, NODES_PER_BLOCK);      
        poten += __shfl_down(poten, offset, NODES_PER_BLOCK);
        idt2max = fmax(idt2max, __shfl_down(idt2max, offset, NODES_PER_BLOCK));
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
        sumx += __shfl_down(sumx, offset, NODES_PER_BLOCK_PART);
        sumy += __shfl_down(sumy, offset, NODES_PER_BLOCK_PART);
        sumz += __shfl_down(sumz, offset, NODES_PER_BLOCK_PART);
        poten += __shfl_down(poten, offset, NODES_PER_BLOCK_PART);
        idt2max = fmax(idt2max, __shfl_down(idt2max, offset, NODES_PER_BLOCK_PART));
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

void run_EWALD_KERNEL_Large(workRequest *wr, cudaStream_t kernel_stream,void** devBuffers) {
  int *Ewaldptr = (int*)wr->getUserData();
  cudaMemcpyToSymbol(cachedData, wr->buffers[EWALD_READ_ONLY_DATA_IDX].host_buffer, sizeof(EwaldReadOnlyData), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(ewt, wr->buffers[EWALD_TABLE_IDX].host_buffer, NEWH * sizeof(EwtData), 0, cudaMemcpyHostToDevice);
  EwaldKernel<<<wr->grid_dim, wr->block_dim, wr->shared_mem, kernel_stream>>>
    ((CompactPartData *)devBuffers[LOCAL_PARTICLE_CORES],
     (VariablePartData *)devBuffers[LOCAL_PARTICLE_VARS],
     (int *)devBuffers[wr->buffers[PARTICLE_TABLE_IDX].id], 1,
     Ewaldptr[0], Ewaldptr[1]);
#ifdef CUDA_PRINT_ERRORS
  printf("EWALD_KERNEL: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
}

void run_EWALD_KERNEL_Small(workRequest *wr, cudaStream_t kernel_stream,void** devBuffers) {
  int *Ewaldptr = (int*)wr->getUserData();
  cudaMemcpyToSymbol(cachedData, wr->buffers[EWALD_READ_ONLY_DATA_IDX].host_buffer, sizeof(EwaldReadOnlyData), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol(ewt, wr->buffers[EWALD_TABLE_IDX].host_buffer, NEWH * sizeof(EwtData), 0, cudaMemcpyHostToDevice);  
  EwaldKernel<<<wr->grid_dim, wr->block_dim, wr->shared_mem, kernel_stream>>>
    ((CompactPartData *)devBuffers[LOCAL_PARTICLE_CORES],
     (VariablePartData *)devBuffers[LOCAL_PARTICLE_VARS],
     (int *)NULL, 0,
     Ewaldptr[0], Ewaldptr[1]);
#ifdef CUDA_PRINT_ERRORS
  printf("EWALD_KERNEL: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
}

extern unsigned int timerHandle; 

void EwaldHostMemorySetup(EwaldData *h_idata, int nParticles, int nEwhLoop, int largephase) {
#ifdef HAPI_MEMPOOL
  if(largephase)
    h_idata->EwaldMarkers = (int *) hapiPoolMalloc((nParticles)*sizeof(int));
  else
    h_idata->EwaldMarkers = NULL;
  h_idata->ewt = (EwtData *) hapiPoolMalloc((nEwhLoop)*sizeof(EwtData));
  h_idata->cachedData = (EwaldReadOnlyData *) hapiPoolMalloc(sizeof(EwaldReadOnlyData));
#else
  if(largephase)
    cudaMallocHost((void**)&(h_idata->EwaldMarkers), (nParticles)*sizeof(int));
  else
    h_idata->EwaldMarkers = NULL;
  cudaMallocHost((void**)&(h_idata->ewt), nEwhLoop * sizeof(EwtData));
  cudaMallocHost((void**)&(h_idata->cachedData), sizeof(EwaldReadOnlyData));
#endif
}

void EwaldHostMemoryFree(EwaldData *h_idata, int largephase) {
#ifdef HAPI_MEMPOOL
  if(largephase)
    hapiPoolFree(h_idata->EwaldMarkers); 
  hapiPoolFree(h_idata->ewt); 
  hapiPoolFree(h_idata->cachedData); 
#else
  if(largephase)
    cudaFreeHost(h_idata->EwaldMarkers); 
  cudaFreeHost(h_idata->ewt); 
  cudaFreeHost(h_idata->cachedData); 
#endif
}

/** @brief Set up CUDA kernels to perform Ewald sum.
 *  @param h_idata Host data buffers
 *  @param cb Callback
 *  @param myIndex Chare index on this node that called this request.
 *  
 *  The "top" and "bottom" Ewlad kernels have been combined:
 *    "top" for the real space loop,
 *    "bottom" for the k-space loop.
 *  
 */
#ifdef CUDA_INSTRUMENT_WRS
void EwaldHost(EwaldData *h_idata, void *cb, int myIndex, char phase, int largephase)
#else
void EwaldHost(EwaldData *h_idata, void *cb, int myIndex, int largephase)
#endif
{

  int n = h_idata->cachedData->n;
  int numBlocks = (int) ceilf((float)n/BLOCK_SIZE);
  int nEwhLoop = h_idata->cachedData->nEwhLoop;
  assert(nEwhLoop <= NEWH);

  workRequest* EwaldKernel = hapiCreateWorkRequest();

  EwaldKernel->setExecParams(numBlocks, BLOCK_SIZE);

  //Range of particles on GPU
  EwaldKernel->setUserData(h_idata->EwaldRange, sizeof(h_idata->EwaldRange));

  /* schedule buffers for transfer to the GPU */
  bool transferToDevice, freeBuffer;
  size_t size;
  if(largephase) transferToDevice = true;
  else transferToDevice = false;
  if(largephase) freeBuffer = true;
  else freeBuffer = false;
  if(largephase) size = n * sizeof(int);
  else size = 0;
  EwaldKernel->addBuffer(h_idata->EwaldMarkers, size, transferToDevice,
                             false, freeBuffer);

  EwaldKernel->addBuffer(h_idata->cachedData, sizeof(EwaldReadOnlyData),
                                false, false, false,
                                NUM_GRAVITY_BUFS + EWALD_READ_ONLY_DATA);

  EwaldKernel->addBuffer(h_idata->ewt, nEwhLoop * sizeof(EwtData), false, false,
                         false, NUM_GRAVITY_BUFS + EWALD_TABLE);

  /* See NUM_BUFFERS define in
   * charm/src/arch/cuda/hybridAPI/cuda-hybrid-api.cu
   * BufferIDs larger than this are assigned by the runtime.
   */
  assert(NUM_GRAVITY_BUFS + EWALD_TABLE < 256);

  EwaldKernel->setDeviceToHostCallback(cb);
  if(largephase){
    EwaldKernel->setRunKernel(run_EWALD_KERNEL_Large);
    EwaldKernel->setTraceName("EwaldLarge");
  }else{
    EwaldKernel->setRunKernel(run_EWALD_KERNEL_Small);
    EwaldKernel->setTraceName("EwaldSmall");
  }
#ifdef CUDA_INSTRUMENT_WRS
  EwaldKernel->chareIndex = myIndex;
  EwaldKernel->compType = EWALD_KERNEL;
  EwaldKernel->compPhase = phase;
#endif
  hapiEnqueue(EwaldKernel);
#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
  printf("[%d] in EwaldHost, enqueued EwaldKernel\n", myIndex);
#endif
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
  cudatype g0, g1, g2, g3, g4, g5;
  cudatype Q2, Q2mirx, Q2miry, Q2mirz, Q2mir, Qta; 
  int ix, iy, iz, bInHole, bInHolex, bInHolexy;

#ifdef HEXADECAPOLE
  MomcData *mom = &(cachedData->momcRoot);
  MultipoleMomentsData *momQuad = &(cachedData->mm);
  cudatype xx,xxx,xxy,xxz,yy,yyy,yyz,xyy,zz,zzz,xzz,yzz,xy,xyz,xz,yz;
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
           * to be changed. See line 480 in Ewald.C.
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
	  alphan *= 2*cachedData->alpha2;
	  g4 = alphan*((1.0/11.0)*r2 - (1.0/9.0));
	  alphan *= 2*cachedData->alpha2;
	  g5 = alphan*((1.0/13.0)*r2 - (1.0/11.0));
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
	  alphan *= 2*cachedData->alpha2;
	  g4 = 7*g3*dir2 + alphan*a;
	  alphan *= 2*cachedData->alpha2;
	  g5 = 9*g4*dir2 + alphan*a;
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
