#ifdef _WIN32
#define NOMINMAX
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
// #include <cutil.h>

#include "CudaFunctions.h"
#include "HostCUDA.h"
#include "wr.h"

#define BLOCK_SIZE 64

extern workRequestQueue *wrQueue;
extern void **devBuffers;
extern cudaStream_t kernel_stream;

//__constant__ constantData[88];
//

void allocatePinnedHostMemory(void **ptr, int size){
  cudaMallocHost(ptr, size);
#ifdef CUDA_PRINT_ERRORS
  printf("allocatePinnedHostMemory: %s size: %d\n", cudaGetErrorString( cudaGetLastError() ), size);
#endif
}

void freePinnedHostMemory(void *ptr){
  cudaFreeHost(ptr);
#ifdef CUDA_PRINT_ERRORS
  printf("freePinnedHostMemory: %s\n", cudaGetErrorString( cudaGetLastError() ));
#endif
}

void DataManagerTransferLocalTree(CudaMultipoleMoments *moments, int nMoments, CompactPartData *compactParts, int nCompactParts) {

	workRequest transferKernel;
	dataInfo *buf;
        int size;

        // operation will not invoke a kernel
	transferKernel.dimGrid = dim3(0);
	transferKernel.dimBlock = dim3(0);
	transferKernel.smemSize = 0;

	transferKernel.nBuffers = DM_TRANSFER_LOCAL_NBUFFERS;

	/* schedule two buffers for transfer to the GPU */
	transferKernel.bufferInfo = (dataInfo *) malloc(transferKernel.nBuffers * sizeof(dataInfo));

	buf = &(transferKernel.bufferInfo[LOCAL_MOMENTS_IDX]);
	buf->bufferID = LOCAL_MOMENTS;
	buf->transferToDevice = YES;
	buf->transferFromDevice = NO;
	buf->freeBuffer = NO;
        size = (nMoments) * sizeof(CudaMultipoleMoments);
	buf->size = size;

#ifdef CUDA_PRINT_ERRORS
        printf("DMLocal 0000: %s\n", cudaGetErrorString( cudaGetLastError() ));
#endif
#ifdef CUDA_USE_CUDAMALLOCHOST
        cudaMallocHost(&buf->hostBuffer, size);
#else
        buf->hostBuffer = malloc(size);
#endif

#ifdef CUDA_PRINT_ERRORS
        printf("DMLocal 0: %s hostBuf: 0x%x, size: %d\n", cudaGetErrorString( cudaGetLastError() ), buf->hostBuffer, size );
#endif
        memcpy(buf->hostBuffer, moments, size);

	buf = &(transferKernel.bufferInfo[LOCAL_PARTICLE_CORES_IDX]);
	buf->bufferID = LOCAL_PARTICLE_CORES;
	buf->transferToDevice = YES;
	buf->transferFromDevice = NO;
	buf->freeBuffer = NO;
        size = (nCompactParts)*sizeof(CompactPartData);
        buf->size = size;

#ifdef CUDA_USE_CUDAMALLOCHOST
        cudaMallocHost(&buf->hostBuffer, size);
#else
        buf->hostBuffer = malloc(size);
#endif
#ifdef CUDA_PRINT_ERRORS
        printf("DMLocal 1: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
        memcpy(buf->hostBuffer, compactParts, size);

        VariablePartData *zeroArray;

	buf = &(transferKernel.bufferInfo[LOCAL_PARTICLE_VARS_IDX]);
	buf->bufferID = LOCAL_PARTICLE_VARS;
	buf->transferToDevice = YES;
	buf->transferFromDevice = NO;
	buf->freeBuffer = NO;
        size = (nCompactParts)*sizeof(VariablePartData);
	buf->size = size;

#ifdef CUDA_USE_CUDAMALLOCHOST
        cudaMallocHost(&buf->hostBuffer, size);
#else
        buf->hostBuffer = malloc(size);
#endif
#ifdef CUDA_PRINT_ERRORS
        printf("DMLocal 2: %s\n", cudaGetErrorString( cudaGetLastError() ));
#endif
        zeroArray = (VariablePartData *) buf->hostBuffer;
        for(int i = 0; i < nCompactParts; i++){
          zeroArray[i].a.x = 0.0;
          zeroArray[i].a.y = 0.0;
          zeroArray[i].a.z = 0.0;
          zeroArray[i].potential = 0.0;
        }

	transferKernel.callbackFn = 0;
	transferKernel.id = DM_TRANSFER_LOCAL;
	enqueue(wrQueue, &transferKernel);

}

void DataManagerTransferRemoteChunk(CudaMultipoleMoments *moments, int nMoments, CompactPartData *compactParts, int nCompactParts) {

  workRequest transferKernel;
  dataInfo *buf;
  int size;

  // operation will not invoke a kernel
  transferKernel.dimGrid = dim3(0);
  transferKernel.dimBlock = dim3(0);
  transferKernel.smemSize = 0;

  transferKernel.nBuffers = DM_TRANSFER_REMOTE_CHUNK_NBUFFERS;

  transferKernel.bufferInfo = (dataInfo *) malloc(transferKernel.nBuffers * sizeof(dataInfo));

  buf = &(transferKernel.bufferInfo[REMOTE_MOMENTS_IDX]);
  buf->bufferID = REMOTE_MOMENTS;
  buf->transferToDevice = YES;
  buf->transferFromDevice = NO;
  buf->freeBuffer = NO;
  size = (nMoments) * sizeof(CudaMultipoleMoments);
  buf->size = size;

  if(size > 0){
#ifdef CUDA_USE_CUDAMALLOCHOST
    cudaMallocHost(&buf->hostBuffer, size);
#else
    buf->hostBuffer = malloc(size);
#endif
#ifdef CUDA_PRINT_ERRORS
    printf("DMRemote 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
    memcpy(buf->hostBuffer, moments, size);
  }

  buf = &(transferKernel.bufferInfo[REMOTE_PARTICLE_CORES_IDX]);
  buf->bufferID = REMOTE_PARTICLE_CORES;
  buf->transferToDevice = YES;
  buf->transferFromDevice = NO;
  buf->freeBuffer = NO;
  size = (nCompactParts)*sizeof(CompactPartData);  
  buf->size = size;

  if(size > 0){
#ifdef CUDA_USE_CUDAMALLOCHOST
    cudaMallocHost(&buf->hostBuffer, size);
#else
    buf->hostBuffer = malloc(size);
#endif
#ifdef CUDA_PRINT_ERRORS
    printf("DMRemote 1: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
    memcpy(buf->hostBuffer, compactParts, size);
  }

  transferKernel.callbackFn = 0;
  transferKernel.id = DM_TRANSFER_REMOTE_CHUNK;
  enqueue(wrQueue, &transferKernel);

}


void TreePieceCellListDataTransferLocal(CudaRequest *data){
	int numBlocks = data->numBucketsPlusOne-1;

	workRequest gravityKernel;
	//dataInfo *buffer, *partCoreBuffer;

	gravityKernel.dimGrid = dim3(numBlocks);
	gravityKernel.dimBlock = dim3(BLOCK_SIZE);
	gravityKernel.smemSize = 0;

	gravityKernel.nBuffers = TP_GRAVITY_LOCAL_NBUFFERS;

	/* schedule buffers for transfer to the GPU */
	gravityKernel.bufferInfo = (dataInfo *) malloc(gravityKernel.nBuffers * sizeof(dataInfo));

	TreePieceCellListDataTransferBasic(data, &gravityKernel);

	gravityKernel.callbackFn = data->cb;
	gravityKernel.id = TP_GRAVITY_LOCAL;
	enqueue(wrQueue, &gravityKernel);
}

void TreePieceCellListDataTransferRemote(CudaRequest *data){
	int numBlocks = data->numBucketsPlusOne-1;

	workRequest gravityKernel;
	//dataInfo *buffer, *partCoreBuffer;

	gravityKernel.dimGrid = dim3(numBlocks);
	gravityKernel.dimBlock = dim3(BLOCK_SIZE);
	gravityKernel.smemSize = 0;

	gravityKernel.nBuffers = TP_NODE_GRAVITY_REMOTE_NBUFFERS;

	/* schedule buffers for transfer to the GPU */
	gravityKernel.bufferInfo = (dataInfo *) malloc(gravityKernel.nBuffers * sizeof(dataInfo));

	TreePieceCellListDataTransferBasic(data, &gravityKernel);

	gravityKernel.callbackFn = data->cb;
	gravityKernel.id = TP_GRAVITY_REMOTE;
	enqueue(wrQueue, &gravityKernel);
}

void TreePieceCellListDataTransferRemoteResume(CudaRequest *data, CudaMultipoleMoments *missedMoments, int numMissedMoments){
  int numBlocks = data->numBucketsPlusOne-1;
  int size;

  workRequest gravityKernel;
  dataInfo *buffer;

  gravityKernel.dimGrid = dim3(numBlocks);
  gravityKernel.dimBlock = dim3(BLOCK_SIZE);
  gravityKernel.smemSize = 0;

  gravityKernel.nBuffers = TP_NODE_GRAVITY_REMOTE_RESUME_NBUFFERS;

  /* schedule buffers for transfer to the GPU */
  gravityKernel.bufferInfo = (dataInfo *) malloc(gravityKernel.nBuffers * sizeof(dataInfo));

  TreePieceCellListDataTransferBasic(data, &gravityKernel);

  buffer = &(gravityKernel.bufferInfo[MISSED_MOMENTS_IDX]);
  buffer->bufferID = -1;
  buffer->transferToDevice = YES;
  buffer->transferFromDevice = NO;
  buffer->freeBuffer = YES;
  size = (numMissedMoments) * sizeof(CudaMultipoleMoments);
  buffer->size = size;

  if(size > 0){
#ifdef CUDA_USE_CUDAMALLOCHOST
    cudaMallocHost(&buffer->hostBuffer, size);
#else
    buffer->hostBuffer = malloc(size);
#endif
#ifdef CUDA_PRINT_ERRORS
    printf("TPRR 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
    memcpy(buffer->hostBuffer, missedMoments, size);
  }

  ParameterStruct *ptr = (ParameterStruct *)gravityKernel.userData;
  ptr->numEntities = numMissedMoments;

  gravityKernel.callbackFn = data->cb;
  gravityKernel.id = TP_GRAVITY_REMOTE_RESUME;
  enqueue(wrQueue, &gravityKernel);
}

void TreePieceCellListDataTransferBasic(CudaRequest *data, workRequest *gravityKernel){
	dataInfo *buffer;
	int numBucketsPlusOne = data->numBucketsPlusOne;
        int numBuckets = numBucketsPlusOne-1;
        int size;

	buffer = &(gravityKernel->bufferInfo[ILCELL_IDX]);
	buffer->bufferID = -1;
	buffer->transferToDevice = YES;
	buffer->transferFromDevice = NO;
	buffer->freeBuffer = YES;
        size = (data->numInteractions) * sizeof(ILCell);
	buffer->size = size;
        buffer->hostBuffer = data->list;

	buffer = &(gravityKernel->bufferInfo[NODE_BUCKET_MARKERS_IDX]);
	buffer->bufferID = -1;
	buffer->transferToDevice = YES;
	buffer->transferFromDevice = NO;
	buffer->freeBuffer = YES;
        size = (numBucketsPlusOne) * sizeof(int);
        buffer->size = size;
        buffer->hostBuffer = data->bucketMarkers;

	buffer = &(gravityKernel->bufferInfo[NODE_BUCKET_START_MARKERS_IDX]);
	buffer->bufferID = -1;
	buffer->transferToDevice = YES;
	buffer->transferFromDevice = NO;
	buffer->freeBuffer = YES;
	size = (numBuckets) * sizeof(int);
        buffer->size = size;
        buffer->hostBuffer = data->bucketStarts;

	buffer = &(gravityKernel->bufferInfo[NODE_BUCKET_SIZES_IDX]);
	buffer->bufferID = -1;
	buffer->transferToDevice = YES;
	buffer->transferFromDevice = NO;
	buffer->freeBuffer = YES;
        buffer->size = size;
        buffer->hostBuffer = data->bucketSizes;

        ParameterStruct *ptr = (ParameterStruct *)malloc(sizeof(ParameterStruct));
        ptr->numInteractions = data->numInteractions;
        ptr->numBucketsPlusOne = numBucketsPlusOne;
        ptr->fperiod = data->fperiod;
        gravityKernel->userData = ptr;
}

void TreePiecePartListDataTransferLocal(CudaRequest *data){
	int numBlocks = data->numBucketsPlusOne-1;

	workRequest gravityKernel;
	//dataInfo *buffer, *partCoreBuffer;

	gravityKernel.dimGrid = dim3(numBlocks);
	gravityKernel.dimBlock = dim3(BLOCK_SIZE);
	gravityKernel.smemSize = 0;

	gravityKernel.nBuffers = TP_GRAVITY_LOCAL_NBUFFERS;

	/* schedule buffers for transfer to the GPU */
	gravityKernel.bufferInfo = (dataInfo *) malloc(gravityKernel.nBuffers * sizeof(dataInfo));

	TreePiecePartListDataTransferBasic(data, &gravityKernel);

	gravityKernel.callbackFn = data->cb;
	gravityKernel.id = TP_PART_GRAVITY_LOCAL;
	enqueue(wrQueue, &gravityKernel);
}

void TreePiecePartListDataTransferRemote(CudaRequest *data){
	int numBlocks = data->numBucketsPlusOne-1;

	workRequest gravityKernel;

	gravityKernel.dimGrid = dim3(numBlocks);
	gravityKernel.dimBlock = dim3(BLOCK_SIZE);
	gravityKernel.smemSize = 0;

	gravityKernel.nBuffers = TP_PART_GRAVITY_REMOTE_NBUFFERS;

	gravityKernel.bufferInfo = (dataInfo *) malloc(gravityKernel.nBuffers * sizeof(dataInfo));

	TreePiecePartListDataTransferBasic(data, &gravityKernel);

	gravityKernel.callbackFn = data->cb;
	gravityKernel.id = TP_PART_GRAVITY_REMOTE;
	enqueue(wrQueue, &gravityKernel);
}

void TreePiecePartListDataTransferRemoteResume(CudaRequest *data, CompactPartData *missedParts, int numMissedParts){
	int numBlocks = data->numBucketsPlusOne-1;
        int size;
        ParameterStruct *ptr;

	workRequest gravityKernel;
        dataInfo *buffer;

	gravityKernel.dimGrid = dim3(numBlocks);
	gravityKernel.dimBlock = dim3(BLOCK_SIZE);
	gravityKernel.smemSize = 0;

	gravityKernel.nBuffers = TP_PART_GRAVITY_REMOTE_RESUME_NBUFFERS;

	/* schedule buffers for transfer to the GPU */
	gravityKernel.bufferInfo = (dataInfo *) malloc(gravityKernel.nBuffers * sizeof(dataInfo));

	TreePiecePartListDataTransferBasic(data, &gravityKernel);

	buffer = &(gravityKernel.bufferInfo[MISSED_PARTS_IDX]);
	buffer->bufferID = -1;
	buffer->transferToDevice = YES;
	buffer->transferFromDevice = NO;
	buffer->freeBuffer = YES;
        size = (numMissedParts) * sizeof(CompactPartData);
        buffer->size = size;

        if(size > 0){
#ifdef CUDA_USE_CUDAMALLOCHOST
          cudaMallocHost(&buffer->hostBuffer, size);
#else
          buffer->hostBuffer = malloc(size);
#endif
#ifdef CUDA_PRINT_ERRORS
          printf("TPPartRR 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
          memcpy(buffer->hostBuffer, missedParts, size);
        }

        ptr = (ParameterStruct *)gravityKernel.userData;
        ptr->numMissedCores = numMissedParts;

	gravityKernel.callbackFn = data->cb;
	gravityKernel.id = TP_PART_GRAVITY_REMOTE_RESUME;
	enqueue(wrQueue, &gravityKernel);
}

void TreePiecePartListDataTransferBasic(CudaRequest *data, workRequest *gravityKernel){
	dataInfo *buffer;

	int numInteractions = data->numInteractions;
	int numBucketsPlusOne = data->numBucketsPlusOne;
        int numBuckets = numBucketsPlusOne-1;
        int size;

	buffer = &(gravityKernel->bufferInfo[ILPART_IDX]);
	buffer->bufferID = -1;
	buffer->transferToDevice = YES;
	buffer->transferFromDevice = NO;
	buffer->freeBuffer = YES;
	size = (numInteractions) * sizeof(ILPart);
        buffer->size = size;
        buffer->hostBuffer = data->list;

	buffer = &(gravityKernel->bufferInfo[PART_BUCKET_MARKERS_IDX]);
	buffer->bufferID = -1;
	buffer->transferToDevice = YES;
	buffer->transferFromDevice = NO;
	buffer->freeBuffer = YES;
	size = (numBucketsPlusOne) * sizeof(int);
        buffer->size = size;
        buffer->hostBuffer = data->bucketMarkers;

	buffer = &(gravityKernel->bufferInfo[PART_BUCKET_START_MARKERS_IDX]);
	buffer->bufferID = -1;
	buffer->transferToDevice = YES;
	buffer->transferFromDevice = NO;
	buffer->freeBuffer = YES;
	buffer->hostBuffer = data->bucketStarts;
	size = (numBuckets) * sizeof(int);
        buffer->size = size;
        buffer->hostBuffer = data->bucketStarts;

	buffer = &(gravityKernel->bufferInfo[PART_BUCKET_SIZES_IDX]);
	buffer->bufferID = -1;
	buffer->transferToDevice = YES;
	buffer->transferFromDevice = NO;
	buffer->freeBuffer = YES;
        buffer->size = size;
        buffer->hostBuffer = data->bucketSizes;
        
        ParameterStruct *ptr = (ParameterStruct *)malloc(sizeof(ParameterStruct));
        ptr->numInteractions = data->numInteractions;
        ptr->numBucketsPlusOne = numBucketsPlusOne;
        ptr->fperiod = data->fperiod;
                                       
        gravityKernel->userData = ptr;
}

/*
__global__ void EwaldTopKernel(GravityParticleData *particleTable) {
*/

/* kernels:

TP_GRAVITY_LOCAL,
TP_GRAVITY_REMOTE,
TP_GRAVITY_REMOTE_RESUME,
TP_PART_GRAVITY_LOCAL,
TP_PART_GRAVITY_REMOTE,
TP_PART_GRAVITY_REMOTE_RESUME

 */

void FreeDataManagerLocalTreeMemory(){
  workRequest gravityKernel;
  dataInfo *buffer;

  gravityKernel.dimGrid = dim3(0);
  gravityKernel.dimBlock = dim3(0);
  gravityKernel.smemSize = 0;

  gravityKernel.nBuffers = DM_TRANSFER_LOCAL_NBUFFERS-1;

  /* schedule buffers for transfer to the GPU */
  gravityKernel.bufferInfo = (dataInfo *) malloc((DM_TRANSFER_LOCAL_NBUFFERS-1) * sizeof(dataInfo));

  buffer = &(gravityKernel.bufferInfo[LOCAL_MOMENTS_IDX]);
  buffer->bufferID = LOCAL_MOMENTS;
  buffer->transferToDevice = NO ;
  buffer->transferFromDevice = NO;
  buffer->freeBuffer = YES;
  buffer->hostBuffer = 0;
  buffer->size = 0;

  buffer = &(gravityKernel.bufferInfo[LOCAL_PARTICLE_CORES_IDX]);
  buffer->bufferID = LOCAL_PARTICLE_CORES;
  buffer->transferToDevice = NO ;
  buffer->transferFromDevice = NO;
  buffer->freeBuffer = YES;
  buffer->hostBuffer = 0;
  buffer->size = 0;

  gravityKernel.callbackFn = 0;
  gravityKernel.id = DM_TRANSFER_FREE_LOCAL;
  enqueue(wrQueue, &gravityKernel);

}

void FreeDataManagerRemoteChunkMemory(int chunk){
  workRequest gravityKernel;
  dataInfo *buffer;

  gravityKernel.dimGrid = dim3(0);
  gravityKernel.dimBlock = dim3(0);
  gravityKernel.smemSize = 0;

  gravityKernel.nBuffers = DM_TRANSFER_REMOTE_CHUNK_NBUFFERS;

  /* schedule buffers for transfer to the GPU */
  gravityKernel.bufferInfo = (dataInfo *) malloc((DM_TRANSFER_REMOTE_CHUNK_NBUFFERS) * sizeof(dataInfo));

  buffer = &(gravityKernel.bufferInfo[REMOTE_MOMENTS_IDX]);
  buffer->bufferID = REMOTE_MOMENTS;
  buffer->transferToDevice = NO ;
  buffer->transferFromDevice = NO;
  buffer->freeBuffer = YES;
  buffer->hostBuffer = 0;
  buffer->size = 0;

  buffer = &(gravityKernel.bufferInfo[REMOTE_PARTICLE_CORES_IDX]);
  buffer->bufferID = REMOTE_PARTICLE_CORES;
  buffer->transferToDevice = NO ;
  buffer->transferFromDevice = NO;
  buffer->freeBuffer = YES;
  buffer->hostBuffer = 0;
  buffer->size = 0;

  gravityKernel.callbackFn = 0;
  gravityKernel.id = DM_TRANSFER_FREE_REMOTE_CHUNK;
  enqueue(wrQueue, &gravityKernel);

}

void TransferParticleVarsBack(VariablePartData *hostBuffer, int size, void *cb){
  workRequest gravityKernel;
  dataInfo *buffer;

#ifdef CUDA_PRINT_TRANSFER_BACK_PARTICLES
  printf("Enqueue kernel to transfer particles back from: 0x%x\n", devBuffers[LOCAL_PARTICLE_VARS]);
#endif
  gravityKernel.dimGrid = dim3(0);
  gravityKernel.dimBlock = dim3(0);
  gravityKernel.smemSize = 0;

  gravityKernel.nBuffers = 1;

  /* schedule buffers for transfer to the GPU */
  gravityKernel.bufferInfo = (dataInfo *) malloc(1 * sizeof(dataInfo));

  buffer = &(gravityKernel.bufferInfo[0]);
  buffer->bufferID = LOCAL_PARTICLE_VARS;
  buffer->transferToDevice = NO ;
  buffer->transferFromDevice = YES;
  buffer->freeBuffer = YES;
  buffer->hostBuffer = hostBuffer;
  buffer->size = size;

  gravityKernel.callbackFn = cb;
  gravityKernel.id = DM_TRANSFER_BACK;
  enqueue(wrQueue, &gravityKernel);
}

/*
void DummyKernel(void *cb){
  workRequest dummy;
  //dataInfo *buffer;

  dummy.dimGrid = dim3(1);
  dummy.dimBlock = dim3(BLOCK_SIZE);
  dummy.smemSize = 0;
  dummy.nBuffers = 0;
  dummy.bufferInfo = 0; //(dataInfo *) malloc(1 * sizeof(dataInfo));

  dummy.callbackFn = cb;
  dummy.id = DUMMY;
  enqueue(wrQueue, &dummy);

}
*/

// kernel selector function
void kernelSelect(workRequest *wr) {

  ParameterStruct *ptr;
  int size;
  switch (wr->id) {

    case DM_TRANSFER_LOCAL:
#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
      printf("DM_TRANSFER_LOCAL\n");
      printf("mom: 0x%x\n", devBuffers[LOCAL_MOMENTS]);
      printf("cores: 0x%x\n", devBuffers[LOCAL_PARTICLE_CORES]);
      printf("vars: 0x%x\n", devBuffers[LOCAL_PARTICLE_VARS]);
#endif
#ifdef CUDA_USE_CUDAMALLOCHOST
      cudaFreeHost(wr->bufferInfo[LOCAL_MOMENTS_IDX].hostBuffer);
#ifdef CUDA_PRINT_ERRORS
      printf("DM_TRANSFER_LOCAL 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
      cudaFreeHost(wr->bufferInfo[LOCAL_PARTICLE_CORES_IDX].hostBuffer);
#ifdef CUDA_PRINT_ERRORS
      printf("DM_TRANSFER_LOCAL 1: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
      cudaFreeHost(wr->bufferInfo[LOCAL_PARTICLE_VARS_IDX].hostBuffer);
#ifdef CUDA_PRINT_ERRORS
      printf("DM_TRANSFER_LOCAL 2: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
#else
      free(wr->bufferInfo[LOCAL_MOMENTS_IDX].hostBuffer);
      free(wr->bufferInfo[LOCAL_PARTICLE_CORES_IDX].hostBuffer);
      free(wr->bufferInfo[LOCAL_PARTICLE_VARS_IDX].hostBuffer);
#endif
      break;

    case DM_TRANSFER_REMOTE_CHUNK:
#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
      printf("DM_TRANSFER_REMOTE_CHUNK\n");
#endif
      size = wr->bufferInfo[REMOTE_MOMENTS_IDX].size;
      if(size > 0){
#ifdef CUDA_USE_CUDAMALLOCHOST
        cudaFreeHost(wr->bufferInfo[REMOTE_MOMENTS_IDX].hostBuffer);
#ifdef CUDA_PRINT_ERRORS
        printf("DM_TRANSFER_REMOTE_CHUNK 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
#else
        free(wr->bufferInfo[REMOTE_MOMENTS_IDX].hostBuffer);
#endif
      }

      size = wr->bufferInfo[REMOTE_PARTICLE_CORES_IDX].size;
      if(size > 0){
#ifdef CUDA_USE_CUDAMALLOCHOST
        cudaFreeHost(wr->bufferInfo[REMOTE_PARTICLE_CORES_IDX].hostBuffer);
#ifdef CUDA_PRINT_ERRORS
        printf("DM_TRANSFER_REMOTE_CHUNK 1: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
#else
        free(wr->bufferInfo[REMOTE_PARTICLE_CORES_IDX].hostBuffer);
#endif
      }
      break;

    case DM_TRANSFER_FREE_LOCAL:
#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
      printf("DM_TRANSFER_FREE_LOCAL\n");
      printf("buffers:\nlocal_particles: (0x%x)\nlocal_particle_vars: (0x%x)\nlocal_moments: (0x%x)\nil_cell: %d (0x%x)\n", 
          devBuffers[LOCAL_PARTICLE_CORES],
          devBuffers[LOCAL_PARTICLE_VARS],
          devBuffers[LOCAL_MOMENTS],
          wr->bufferInfo[ILCELL_IDX].bufferID,
          devBuffers[wr->bufferInfo[ILCELL_IDX].bufferID]
          );

#endif
      break;
    case DM_TRANSFER_FREE_REMOTE_CHUNK:
#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
      printf("DM_TRANSFER_FREE_REMOTE_CHUNK\n");

#endif
      break;
    case DM_TRANSFER_BACK:
#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
      printf("DM_TRANSFER_BACK: 0x%x\n", devBuffers[LOCAL_PARTICLE_VARS]);
#endif
      break;

/*
    case DUMMY: 
#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
      printf("DUMMY\n");
#endif
      break;
*/

    case TP_GRAVITY_LOCAL:
      ptr = (ParameterStruct *)wr->userData;
#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
      printf("TP_GRAVITY_LOCAL buffers:\nlocal_particles: (0x%x)\nlocal_particle_vars: (0x%x)\nlocal_moments: (0x%x)\nil_cell: %d (0x%x)\n", 
          devBuffers[LOCAL_PARTICLE_CORES],
          devBuffers[LOCAL_PARTICLE_VARS],
          devBuffers[LOCAL_MOMENTS],
          wr->bufferInfo[ILCELL_IDX].bufferID,
          devBuffers[wr->bufferInfo[ILCELL_IDX].bufferID]
          );
#endif
#ifndef CUDA_NO_KERNELS
      nodeGravityComputation<<<wr->dimGrid, wr->dimBlock, wr->smemSize, kernel_stream>>>
        (
         (CompactPartData *)devBuffers[LOCAL_PARTICLE_CORES],
         (VariablePartData *)devBuffers[LOCAL_PARTICLE_VARS],
         (CudaMultipoleMoments *)devBuffers[LOCAL_MOMENTS],
         (ILCell *)devBuffers[wr->bufferInfo[ILCELL_IDX].bufferID],
         ptr->numInteractions,
         (int *)devBuffers[wr->bufferInfo[NODE_BUCKET_MARKERS_IDX].bufferID],
         (int *)devBuffers[wr->bufferInfo[NODE_BUCKET_START_MARKERS_IDX].bufferID],
         (int *)devBuffers[wr->bufferInfo[NODE_BUCKET_SIZES_IDX].bufferID],
         ptr->numBucketsPlusOne,
         ptr->fperiod,
         0,
         ptr->numEntities
        );
#endif

#ifdef CUDA_USE_CUDAMALLOCHOST
      cudaFreeHost((ILCell *)wr->bufferInfo[ILCELL_IDX].hostBuffer);
      cudaFreeHost((int *)wr->bufferInfo[NODE_BUCKET_MARKERS_IDX].hostBuffer);
      cudaFreeHost((int *)wr->bufferInfo[NODE_BUCKET_START_MARKERS_IDX].hostBuffer);
      cudaFreeHost((int *)wr->bufferInfo[NODE_BUCKET_SIZES_IDX].hostBuffer);
#ifdef CUDA_PRINT_ERRORS
      printf("TP_GRAVITY_LOCAL 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
#else
      free((ILCell *)wr->bufferInfo[ILCELL_IDX].hostBuffer);
      free((int *)wr->bufferInfo[NODE_BUCKET_MARKERS_IDX].hostBuffer);
      free((int *)wr->bufferInfo[NODE_BUCKET_START_MARKERS_IDX].hostBuffer);
      free((int *)wr->bufferInfo[NODE_BUCKET_SIZES_IDX].hostBuffer);
#endif
      free((ParameterStruct *)ptr);
      break;

    case TP_PART_GRAVITY_LOCAL:
      ptr = (ParameterStruct *)wr->userData;
#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
      printf("TP_PART_GRAVITY_LOCAL buffers:\nlocal_particles: (0x%x)\nlocal_particle_vars: (0x%x)\nil_cell: %d (0x%x)\n", 
          devBuffers[LOCAL_PARTICLE_CORES],
          devBuffers[LOCAL_PARTICLE_VARS],
          wr->bufferInfo[ILPART_IDX].bufferID,
          devBuffers[wr->bufferInfo[ILPART_IDX].bufferID]
          );
#endif
#ifndef CUDA_NO_KERNELS
      particleGravityComputation<<<wr->dimGrid, wr->dimBlock, wr->smemSize, kernel_stream>>>
        (
         (CompactPartData *)devBuffers[LOCAL_PARTICLE_CORES],
         (VariablePartData *)devBuffers[LOCAL_PARTICLE_VARS],
         (CompactPartData *)devBuffers[LOCAL_PARTICLE_CORES],
         (ILPart *)devBuffers[wr->bufferInfo[ILPART_IDX].bufferID],
         ptr->numInteractions,
         (int *)devBuffers[wr->bufferInfo[PART_BUCKET_MARKERS_IDX].bufferID],
         (int *)devBuffers[wr->bufferInfo[PART_BUCKET_START_MARKERS_IDX].bufferID],
         (int *)devBuffers[wr->bufferInfo[PART_BUCKET_SIZES_IDX].bufferID],
         ptr->numBucketsPlusOne,
         ptr->fperiod,0 
        );
#endif

#ifdef CUDA_USE_CUDAMALLOCHOST
      cudaFreeHost((ILPart *)wr->bufferInfo[ILPART_IDX].hostBuffer);
      cudaFreeHost((int *)wr->bufferInfo[PART_BUCKET_MARKERS_IDX].hostBuffer);
      cudaFreeHost((int *)wr->bufferInfo[PART_BUCKET_START_MARKERS_IDX].hostBuffer);
      cudaFreeHost((int *)wr->bufferInfo[PART_BUCKET_SIZES_IDX].hostBuffer);
#ifdef CUDA_PRINT_ERRORS
      printf("TP_PART_GRAVITY_LOCAL 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
#else
      free((ILPart *)wr->bufferInfo[ILPART_IDX].hostBuffer);
      free((int *)wr->bufferInfo[PART_BUCKET_MARKERS_IDX].hostBuffer);
      free((int *)wr->bufferInfo[PART_BUCKET_START_MARKERS_IDX].hostBuffer);
      free((int *)wr->bufferInfo[PART_BUCKET_SIZES_IDX].hostBuffer);
#endif
      free((ParameterStruct *)ptr);
      break;

    case TP_GRAVITY_REMOTE:
      ptr = (ParameterStruct *)wr->userData;
#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
      printf("TP_GRAVITY_REMOTE buffers:\nlocal_particles: (0x%x)\nlocal_particle_vars: (0x%x)\nremote_moments: (0x%x)\nil_cell: %d (0x%x)\n", 
          devBuffers[LOCAL_PARTICLE_CORES],
          devBuffers[LOCAL_PARTICLE_VARS],
          devBuffers[REMOTE_MOMENTS],
          wr->bufferInfo[ILCELL_IDX].bufferID,
          devBuffers[wr->bufferInfo[ILCELL_IDX].bufferID]
          );
#endif
#ifndef CUDA_NO_KERNELS
      nodeGravityComputation<<<wr->dimGrid, wr->dimBlock, wr->smemSize, kernel_stream>>>
        (
         (CompactPartData *)devBuffers[LOCAL_PARTICLE_CORES],
         (VariablePartData *)devBuffers[LOCAL_PARTICLE_VARS],
         (CudaMultipoleMoments *)devBuffers[REMOTE_MOMENTS],
         (ILCell *)devBuffers[wr->bufferInfo[ILCELL_IDX].bufferID],
         ptr->numInteractions,
         (int *)devBuffers[wr->bufferInfo[NODE_BUCKET_MARKERS_IDX].bufferID],
         (int *)devBuffers[wr->bufferInfo[NODE_BUCKET_START_MARKERS_IDX].bufferID],
         (int *)devBuffers[wr->bufferInfo[NODE_BUCKET_SIZES_IDX].bufferID],
         ptr->numBucketsPlusOne,
         ptr->fperiod,
         1,
         ptr->numEntities
        );
#endif

#ifdef CUDA_USE_CUDAMALLOCHOST
      cudaFreeHost((ILCell *)wr->bufferInfo[ILCELL_IDX].hostBuffer);
      cudaFreeHost((int *)wr->bufferInfo[NODE_BUCKET_MARKERS_IDX].hostBuffer);
      cudaFreeHost((int *)wr->bufferInfo[NODE_BUCKET_START_MARKERS_IDX].hostBuffer);
      cudaFreeHost((int *)wr->bufferInfo[NODE_BUCKET_SIZES_IDX].hostBuffer);
#ifdef CUDA_PRINT_ERRORS
      printf("TP_GRAVITY_REMOTE 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
#else
      free((ILCell *)wr->bufferInfo[ILCELL_IDX].hostBuffer);
      free((int *)wr->bufferInfo[NODE_BUCKET_MARKERS_IDX].hostBuffer);
      free((int *)wr->bufferInfo[NODE_BUCKET_START_MARKERS_IDX].hostBuffer);
      free((int *)wr->bufferInfo[NODE_BUCKET_SIZES_IDX].hostBuffer);
#endif
      free((ParameterStruct *)ptr);
      break;

    case TP_PART_GRAVITY_REMOTE:
      ptr = (ParameterStruct *)wr->userData;
#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
      printf("TP_PART_GRAVITY_REMOTE buffers:\nlocal_particles: (0x%x)\nlocal_particle_vars: (0x%x)\nil_cell: %d (0x%x)\n", 
          devBuffers[LOCAL_PARTICLE_CORES],
          devBuffers[LOCAL_PARTICLE_VARS],
          wr->bufferInfo[ILPART_IDX].bufferID,
          devBuffers[wr->bufferInfo[ILPART_IDX].bufferID]
          );
#endif
#ifndef CUDA_NO_KERNELS
      particleGravityComputation<<<wr->dimGrid, wr->dimBlock, wr->smemSize, kernel_stream>>>
        (
         (CompactPartData *)devBuffers[LOCAL_PARTICLE_CORES],
         (VariablePartData *)devBuffers[LOCAL_PARTICLE_VARS],
         (CompactPartData *)devBuffers[REMOTE_PARTICLE_CORES],
         (ILPart *)devBuffers[wr->bufferInfo[ILPART_IDX].bufferID],
         ptr->numInteractions,
         (int *)devBuffers[wr->bufferInfo[PART_BUCKET_MARKERS_IDX].bufferID],
         (int *)devBuffers[wr->bufferInfo[PART_BUCKET_START_MARKERS_IDX].bufferID],
         (int *)devBuffers[wr->bufferInfo[PART_BUCKET_SIZES_IDX].bufferID],
         ptr->numBucketsPlusOne,
         ptr->fperiod,1
        );
#endif

#ifdef CUDA_USE_CUDAMALLOCHOST
      cudaFreeHost((ILPart *)wr->bufferInfo[ILPART_IDX].hostBuffer);
      cudaFreeHost((int *)wr->bufferInfo[PART_BUCKET_MARKERS_IDX].hostBuffer);
      cudaFreeHost((int *)wr->bufferInfo[PART_BUCKET_START_MARKERS_IDX].hostBuffer);
      cudaFreeHost((int *)wr->bufferInfo[PART_BUCKET_SIZES_IDX].hostBuffer);
#ifdef CUDA_PRINT_ERRORS
      printf("TP_PART_GRAVITY_REMOTE 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
#else
      free((ILPart *)wr->bufferInfo[ILPART_IDX].hostBuffer);
      free((int *)wr->bufferInfo[PART_BUCKET_MARKERS_IDX].hostBuffer);
      free((int *)wr->bufferInfo[PART_BUCKET_START_MARKERS_IDX].hostBuffer);
      free((int *)wr->bufferInfo[PART_BUCKET_SIZES_IDX].hostBuffer);
#endif
      free((ParameterStruct *)ptr);
      break;

    case TP_GRAVITY_REMOTE_RESUME:
      ptr = (ParameterStruct *)wr->userData;
#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
      printf("TP_GRAVITY_REMOTE_RESUME buffers:\nlocal_particles: (0x%x)\nlocal_particle_vars: (0x%x)\nmissed_moments: %d (0x%x)\nil_cell: %d (0x%x)\n", 
          devBuffers[LOCAL_PARTICLE_CORES],
          devBuffers[LOCAL_PARTICLE_VARS],
          wr->bufferInfo[MISSED_MOMENTS_IDX].bufferID,
          devBuffers[wr->bufferInfo[MISSED_MOMENTS_IDX].bufferID],
          wr->bufferInfo[ILCELL_IDX].bufferID,
          devBuffers[wr->bufferInfo[ILCELL_IDX].bufferID]
          );
#endif
#ifndef CUDA_NO_KERNELS
      nodeGravityComputation<<<wr->dimGrid, wr->dimBlock, wr->smemSize, kernel_stream>>>
        (
         (CompactPartData *)devBuffers[LOCAL_PARTICLE_CORES],
         (VariablePartData *)devBuffers[LOCAL_PARTICLE_VARS],
         (CudaMultipoleMoments *)devBuffers[wr->bufferInfo[MISSED_MOMENTS_IDX].bufferID],
         (ILCell *)devBuffers[wr->bufferInfo[ILCELL_IDX].bufferID],
         ptr->numInteractions,
         (int *)devBuffers[wr->bufferInfo[NODE_BUCKET_MARKERS_IDX].bufferID],
         (int *)devBuffers[wr->bufferInfo[NODE_BUCKET_START_MARKERS_IDX].bufferID],
         (int *)devBuffers[wr->bufferInfo[NODE_BUCKET_SIZES_IDX].bufferID],
         ptr->numBucketsPlusOne,
         ptr->fperiod,
         2,
         ptr->numEntities
        );
#endif

#ifdef CUDA_USE_CUDAMALLOCHOST
      cudaFreeHost((CudaMultipoleMoments *)wr->bufferInfo[MISSED_MOMENTS_IDX].hostBuffer);
      cudaFreeHost((ILCell *)wr->bufferInfo[ILCELL_IDX].hostBuffer);
      cudaFreeHost((int *)wr->bufferInfo[NODE_BUCKET_MARKERS_IDX].hostBuffer);
      cudaFreeHost((int *)wr->bufferInfo[NODE_BUCKET_START_MARKERS_IDX].hostBuffer);
      cudaFreeHost((int *)wr->bufferInfo[NODE_BUCKET_SIZES_IDX].hostBuffer);
#ifdef CUDA_PRINT_ERRORS
      printf("TP_GRAVITY_REMOTE_RESUME 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
#else
      free((CudaMultipoleMoments *)wr->bufferInfo[MISSED_MOMENTS_IDX].hostBuffer);
      free((ILCell *)wr->bufferInfo[ILCELL_IDX].hostBuffer);
      free((int *)wr->bufferInfo[NODE_BUCKET_MARKERS_IDX].hostBuffer);
      free((int *)wr->bufferInfo[NODE_BUCKET_START_MARKERS_IDX].hostBuffer);
      free((int *)wr->bufferInfo[NODE_BUCKET_SIZES_IDX].hostBuffer);
#endif
      free((ParameterStruct *)ptr);
      break;

    case TP_PART_GRAVITY_REMOTE_RESUME:
      ptr = (ParameterStruct *)wr->userData;
#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
      printf("TP_PART_GRAVITY_REMOTE_RESUME buffers:\nlocal_particles: (0x%x)\nlocal_particle_vars: (0x%x)\nmissed_parts: %d (0x%x)\nil_cell: %d (0x%x)\n", 
          devBuffers[LOCAL_PARTICLE_CORES],
          devBuffers[LOCAL_PARTICLE_VARS],
          wr->bufferInfo[MISSED_PARTS_IDX].bufferID,
          devBuffers[wr->bufferInfo[MISSED_PARTS_IDX].bufferID],
          wr->bufferInfo[ILPART_IDX].bufferID,
          devBuffers[wr->bufferInfo[ILPART_IDX].bufferID]
          );
#endif
#ifndef CUDA_NO_KERNELS
      particleGravityComputation<<<wr->dimGrid, wr->dimBlock, wr->smemSize, kernel_stream>>>
        (
         (CompactPartData *)devBuffers[LOCAL_PARTICLE_CORES],
         (VariablePartData *)devBuffers[LOCAL_PARTICLE_VARS],
         (CompactPartData *)devBuffers[wr->bufferInfo[MISSED_PARTS_IDX].bufferID],
         (ILPart *)devBuffers[wr->bufferInfo[ILPART_IDX].bufferID],
         ptr->numMissedCores,
         (int *)devBuffers[wr->bufferInfo[PART_BUCKET_MARKERS_IDX].bufferID],
         (int *)devBuffers[wr->bufferInfo[PART_BUCKET_START_MARKERS_IDX].bufferID],
         (int *)devBuffers[wr->bufferInfo[PART_BUCKET_SIZES_IDX].bufferID],
         ptr->numBucketsPlusOne,
         ptr->fperiod,2
        );
#endif

#ifdef CUDA_USE_CUDAMALLOCHOST
      cudaFreeHost((CompactPartData *)wr->bufferInfo[MISSED_PARTS_IDX].hostBuffer);
      cudaFreeHost((ILPart *)wr->bufferInfo[ILPART_IDX].hostBuffer);
      cudaFreeHost((int *)wr->bufferInfo[PART_BUCKET_MARKERS_IDX].hostBuffer);
      cudaFreeHost((int *)wr->bufferInfo[PART_BUCKET_START_MARKERS_IDX].hostBuffer);
      cudaFreeHost((int *)wr->bufferInfo[PART_BUCKET_SIZES_IDX].hostBuffer);
#ifdef CUDA_PRINT_ERRORS
      printf("TP_PART_GRAVITY_REMOTE_RESUME 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
#else
      free((CompactPartData *)wr->bufferInfo[MISSED_PARTS_IDX].hostBuffer);
      free((ILPart *)wr->bufferInfo[ILPART_IDX].hostBuffer);
      free((int *)wr->bufferInfo[PART_BUCKET_MARKERS_IDX].hostBuffer);
      free((int *)wr->bufferInfo[PART_BUCKET_START_MARKERS_IDX].hostBuffer);
      free((int *)wr->bufferInfo[PART_BUCKET_SIZES_IDX].hostBuffer);
#endif
      free((ParameterStruct *)ptr);
      break;

    default:
      printf("error: id %d not valid\n", wr->id);
      break;

  }
}

/*
 * Kernels
 */

#define GROUP(t)  ((t)/MAX_THREADS_PER_GROUP)
#define GROUP_INDEX(t) ((t)%MAX_THREADS_PER_GROUP)

__global__ void nodeGravityComputation(
		CompactPartData *particleCores,
		VariablePartData *particleVars,
		CudaMultipoleMoments *moments,
		ILCell *ils,
		int numInteractions,
		int *ilmarks,
		int *bucketStarts,
		int *bucketSizes,
		int numBucketsPlusOne, cudatype fperiod, int type, int numNodes){

  // each thread has its own storage for these
  __shared__ CudaVector3D acc[THREADS_PER_BLOCK];
  __shared__ cudatype pot[THREADS_PER_BLOCK];
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

#if defined  __DEVICE_EMULATION__ && defined CUDA_EMU_KERNEL_NODE_PRINTS
  int tp, id;
  if(thread == 0)
    printf("numMissedCores: %d, bucketSize: %d, start: %d, end: %d\n", numInteractions, bucketSize, start, end);
#endif


#ifdef __DEVICE_EMULATION__
  if(blockIdx.x == 0){
    //printf("t: %d, blen: %d, llen: %d\n", threadIdx.x, blen, llen);
    //printf("group: %d, particle: %d, ngroups: %d, groupSize: %d\n", group, particle, ngroups, groupSize);
  }
#endif

  for(int particle = 0; particle < bucketSize; particle++){
    if(thread == 0){
      // load shared_particle_core
      shared_particle_core = particleCores[bucketStart+particle];
    }
    __syncthreads();

#if defined __DEVICE_EMULATION__ && defined CUDA_EMU_KERNEL_NODE_PRINTS
    id = shared_particle_core.id;
    tp = shared_particle_core.tp;
    if(thread == 0){
      printf("targetVars: 0x%x targetVars[%d,%d].a.x = %f\n",
                             particleVars, tp,id,
                             particleVars[bucketStart+particle].a.x);
    }
#endif

    acc[thread].x = 0;
    acc[thread].y = 0;
    acc[thread].z = 0;
    pot[thread] = 0;

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
      }// end if rsq != 0
#if defined  __DEVICE_EMULATION__ && defined CUDA_EMU_KERNEL_NODE_PRINTS
      if(1){
        CudaMultipoleMoments mm = m[thread];
        printf("NODE target particle (%d,%d) type %d with moments[%d], acc: (%f,%f,%f)\n", tp, id, 
            type,
            ilc.index,
            acc[thread].x, 
            acc[thread].y,
            acc[thread].z);
        printf("{m=%f,s=%f,r=%f}\n", mm.totalMass, mm.soft, mm.radius);
        printf("{xx=%f,xy=%f,xz=%f}\n", mm.xx, mm.xy, mm.xz);
        printf("{yy=%f,yz=%f,zz=%f}\n", mm.yy, mm.yz, mm.zz);
      }
#endif

    }// for each node in list
    __syncthreads();
    // at this point, the total force on particle is distributed among
    // all active threads;
    // reduce.
    // TODO: make this a parallel reduction

    cudatype sumx, sumy, sumz, poten;
    sumx = sumy = sumz = poten = 0.0;
    if(thread == 0){
      for(int i = 0; i < THREADS_PER_BLOCK; i++){
        sumx += acc[i].x;
        sumy += acc[i].y;
        sumz += acc[i].z;
        poten += pot[i];
      }
#if defined  __DEVICE_EMULATION__ && defined CUDA_EMU_KERNEL_NODE_PRINTS
      printf("aggregate acc[%d,%d] = (%f,%f,%f)\n", tp, id, sumx, sumy, sumz);
      /*
      printf("before (%d,%d) type: NODE %d acc: (%f,%f,%f)\n", tp, id, type,
                                                    particleVars[bucketStart+particle].a.x,
                                                    particleVars[bucketStart+particle].a.y,
                                                    particleVars[bucketStart+particle].a.z
                                                    );
      */
#endif
      particleVars[bucketStart+particle].a.x += sumx;
      particleVars[bucketStart+particle].a.y += sumy;
      particleVars[bucketStart+particle].a.z += sumz;
      particleVars[bucketStart+particle].potential += poten;
      
#if defined  __DEVICE_EMULATION__ && defined CUDA_EMU_KERNEL_NODE_PRINTS
      printf("after (%d,%d) type: NODE %d acc: (%f,%f,%f)\n", tp, id, type,
                                                    particleVars[bucketStart+particle].a.x,
                                                    particleVars[bucketStart+particle].a.y,
                                                    particleVars[bucketStart+particle].a.z
                                                    );
#endif
    }


    /*
    cudatype sum = 0.0;
    if(thread == 0){
      for(int i = 0; i < THREADS_PER_BLOCK; i++){
        sum += acc[i].x;
      }
      particleVars[bucketStart+particle].a.x += sum;

      sum = 0;
      for(int i = 0; i < THREADS_PER_BLOCK; i++){
        sum += acc[i].y;
      }
      particleVars[bucketStart+particle].a.y += sum;

      sum = 0;
      for(int i = 0; i < THREADS_PER_BLOCK; i++){
        sum += acc[i].z;
      }
      particleVars[bucketStart+particle].a.z += sum;

      sum = 0;
      for(int i = 0; i < THREADS_PER_BLOCK; i++){
        sum += pot[i];
      }
      particleVars[bucketStart+particle].potential += sum;
    }
    */

  }// for each particle in bucket

      
}

__global__ void particleGravityComputation(
                                   CompactPartData *targetCores,
                                   VariablePartData *targetVars,
                                   CompactPartData *sourceCores,
                                   ILPart *ils,
		                   int numInteractions,
                                   int *ilmarks,
		                   int *bucketStarts,
		                   int *bucketSizes,
		                   int numBucketsPlusOne, cudatype fperiod, int type){

  // each thread has its own storage for these
  __shared__ CudaVector3D acc[THREADS_PER_BLOCK];
  __shared__ cudatype pot[THREADS_PER_BLOCK];
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

#if defined  __DEVICE_EMULATION__ && defined CUDA_EMU_KERNEL_PART_PRINTS
  int tp, id;
#endif

  CudaVector3D r;
  cudatype rsq;
  cudatype twoh, a, b;

#if defined  __DEVICE_EMULATION__ && defined CUDA_EMU_KERNEL_PART_PRINTS
  if(thread == 0)
    printf("numMissedCores: %d, bucketSize: %d, start: %d, end: %d\n", numInteractions, bucketSize, start, end);
#endif

  for(int target = 0; target < bucketSize; target++){
    if(thread == 0){
      shared_target_core = targetCores[bucketStart+target];
    }
    __syncthreads();
#if defined __DEVICE_EMULATION__ && defined CUDA_EMU_KERNEL_PART_PRINTS
    id = shared_target_core.id;
    tp = shared_target_core.tp;
    if(thread == 0){
      printf("targetVars: 0x%x targetVars[%d,%d].a.x = %f\n",
                             targetVars, tp,id,
                             targetVars[bucketStart+target].a.x);
    }
#endif

    acc[thread].x = 0;
    acc[thread].y = 0;
    acc[thread].z = 0;
    pot[thread] = 0;

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
          //SPLINE(rsq, twoh, a, b);
          //SPLINE(r2, twoh, a, b);
          //expanded below:
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
        }// if rsq != 0
#if 0 && defined  __DEVICE_EMULATION__ && defined CUDA_EMU_KERNEL_PART_PRINTS
        if(tp == 1 && type == 2){
          //printf("target particle (%d,%d) type %d with sourceCores[%d+%d]\n", tp, id, type, ilpindex,particle);
          printf("target particle (%d,%d) type %d acc: (%f,%f,%f)\n", tp, id, type, acc[thread].x, 
                                                                                    acc[thread].y,
                                                                                    acc[thread].z);
        }
#endif
      }// for each particle in source bucket
    }// for each source bucket 

    __syncthreads();
  
    cudatype sumx, sumy, sumz, poten;
    sumx = sumy = sumz = poten = 0.0;
    if(thread == 0){
      for(int i = 0; i < THREADS_PER_BLOCK; i++){
        sumx += acc[i].x;
        sumy += acc[i].y;
        sumz += acc[i].z;
        poten += pot[i];
      }
#if defined  __DEVICE_EMULATION__ && defined CUDA_EMU_KERNEL_PART_PRINTS
      //printf("before acc[%d,%d] = (%f,%f,%f)\n", tp, id, sumx, sumy, sumz);
      printf("before (%d,%d) type: %d acc: (%f,%f,%f)\n", tp, id, type,
                                                    targetVars[bucketStart+target].a.x,
                                                    targetVars[bucketStart+target].a.y,
                                                    targetVars[bucketStart+target].a.z
                                                    );
#endif
      targetVars[bucketStart+target].a.x += sumx;
      targetVars[bucketStart+target].a.y += sumy;
      targetVars[bucketStart+target].a.z += sumz;
      targetVars[bucketStart+target].potential += poten;
      
#if defined  __DEVICE_EMULATION__ && defined CUDA_EMU_KERNEL_PART_PRINTS
      printf("after (%d,%d) type: %d acc: (%f,%f,%f)\n", tp, id, type,
                                                    targetVars[bucketStart+target].a.x,
                                                    targetVars[bucketStart+target].a.y,
                                                    targetVars[bucketStart+target].a.z
                                                    );
#endif
    }

    /*
    cudatype sum = 0.0;
    if(thread == 0){
      for(int i = 0; i < THREADS_PER_BLOCK; i++){
#if defined  __DEVICE_EMULATION__ && defined CUDA_EMU_KERNEL_PART_PRINTS
        printf("acc[%d][%d].x = %f\n", bucket, i, acc[i].x);
#endif
        sum += acc[i].x;
      }
      targetVars[bucketStart+target].a.x += sum;

      sum = 0;
      for(int i = 0; i < THREADS_PER_BLOCK; i++){
#if defined  __DEVICE_EMULATION__ && defined CUDA_EMU_KERNEL_PART_PRINTS
        printf("acc[%d][%d].y = %f\n", bucket, i, acc[i].y);
#endif
        sum += acc[i].y;
      }
      targetVars[bucketStart+target].a.y += sum;

      sum = 0;
      for(int i = 0; i < THREADS_PER_BLOCK; i++){
#if defined  __DEVICE_EMULATION__ && defined CUDA_EMU_KERNEL_PART_PRINTS
        printf("acc[%d][%d].z = %f\n", bucket, i, acc[i].z);
#endif
        sum += acc[i].z;
      }
      targetVars[bucketStart+target].a.z += sum;

      sum = 0;
      for(int i = 0; i < THREADS_PER_BLOCK; i++){
        sum += pot[i];
      }
      targetVars[bucketStart+target].potential += sum;
    }
    */

  }// for each target part
}

