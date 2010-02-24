#ifdef _WIN32
#define NOMINMAX
#endif

#ifdef CUDA_MEMPOOL
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

#include "CudaFunctions.h"
#include "HostCUDA.h"
#include "EwaldCUDA.h"

#include "wr.h"

#define BLOCK_SIZE_GRAV 64

//extern workRequestQueue *wrQueue;
//extern void **devBuffers;
//extern cudaStream_t kernel_stream;

__device__ __constant__ EwaldReadOnlyData cachedData[1];
__device__ __constant__ EwtData ewt[sizeof(EwtData) * NEWH];  


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
#ifdef CUDA_MEMPOOL
  *ptr = hapi_poolMalloc(size);
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
#ifdef CUDA_MEMPOOL
  hapi_poolFree(ptr);
#else
  delayedFree(ptr);
#endif
#ifdef CUDA_PRINT_ERRORS
  printf("freePinnedHostMemory: %s\n", cudaGetErrorString( cudaGetLastError() ));
#endif
}

#ifdef CUDA_INSTRUMENT_WRS
void DataManagerTransferLocalTree(CudaMultipoleMoments *moments, int nMoments, CompactPartData *compactParts, int nCompactParts, int mype, char phase) {
#else
void DataManagerTransferLocalTree(CudaMultipoleMoments *moments, int nMoments, CompactPartData *compactParts, int nCompactParts, int mype) {
#endif

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
	buf->freeBuffer = NO;
        size = (nMoments) * sizeof(CudaMultipoleMoments);
	buf->size = size;
	buf->transferToDevice = size > 0 ? YES : NO;
	buf->transferFromDevice = NO;

        //double mill = 1e6;
        //printf("(%d) DM local moments: %f mbytes\n", mype, 1.0*size/mill);

#ifdef CUDA_PRINT_ERRORS
        printf("(%d) DMLocal 0000: %s\n", mype, cudaGetErrorString( cudaGetLastError() ));
#endif
        if(size > 0){
#ifdef CUDA_USE_CUDAMALLOCHOST
#ifdef CUDA_MEMPOOL
          buf->hostBuffer = hapi_poolMalloc(size);
#else
          cudaMallocHost(&buf->hostBuffer, size);
#endif
#else
          buf->hostBuffer = malloc(size);
#endif
        }

#ifdef CUDA_PRINT_ERRORS
        printf("(%d) DMLocal 0: %s hostBuf: 0x%x, size: %d\n", mype, cudaGetErrorString( cudaGetLastError() ), buf->hostBuffer, size );
#endif
        memcpy(buf->hostBuffer, moments, size);

	buf = &(transferKernel.bufferInfo[LOCAL_PARTICLE_CORES_IDX]);
	buf->bufferID = LOCAL_PARTICLE_CORES;
	buf->freeBuffer = NO;
        size = (nCompactParts)*sizeof(CompactPartData);
        buf->size = size;
	buf->transferToDevice = size > 0 ? YES : NO;
	buf->transferFromDevice = NO;

        if(size > 0){
#ifdef CUDA_USE_CUDAMALLOCHOST
#ifdef CUDA_MEMPOOL
          buf->hostBuffer = hapi_poolMalloc(size);
#else
          cudaMallocHost(&buf->hostBuffer, size);
#endif
#else
          buf->hostBuffer = malloc(size);
#endif
#ifdef CUDA_PRINT_ERRORS
          printf("(%d) DMLocal 1: %s\n", mype, cudaGetErrorString( cudaGetLastError() ) );
#endif
          memcpy(buf->hostBuffer, compactParts, size);
        }

        VariablePartData *zeroArray;

	buf = &(transferKernel.bufferInfo[LOCAL_PARTICLE_VARS_IDX]);
	buf->bufferID = LOCAL_PARTICLE_VARS;
	buf->freeBuffer = NO;
        size = (nCompactParts)*sizeof(VariablePartData);
	buf->size = size;
	buf->transferToDevice = size > 0 ? YES : NO;
	buf->transferFromDevice = NO;

        if(size > 0){
#ifdef CUDA_USE_CUDAMALLOCHOST
#ifdef CUDA_MEMPOOL
          buf->hostBuffer = hapi_poolMalloc(size);
#else
          cudaMallocHost(&buf->hostBuffer, size);
#endif
#else
          buf->hostBuffer = malloc(size);
#endif
#ifdef CUDA_PRINT_ERRORS
          printf("(%d) DMLocal 2: %s\n", mype, cudaGetErrorString( cudaGetLastError() ));
#endif
          zeroArray = (VariablePartData *) buf->hostBuffer;
          for(int i = 0; i < nCompactParts; i++){
            zeroArray[i].a.x = 0.0;
            zeroArray[i].a.y = 0.0;
            zeroArray[i].a.z = 0.0;
            zeroArray[i].potential = 0.0;
          }
        }

	transferKernel.callbackFn = 0;
	transferKernel.id = DM_TRANSFER_LOCAL;
#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
        printf("DM LOCAL TREE\n");
#endif
#ifdef CUDA_INSTRUMENT_WRS
        transferKernel.chareIndex = mype;
        transferKernel.compType = transferKernel.id;
        transferKernel.compPhase = phase; 
#endif
	enqueue(wrQueue, &transferKernel);

}

#ifdef CUDA_INSTRUMENT_WRS
void DataManagerTransferRemoteChunk(CudaMultipoleMoments *moments, int nMoments, CompactPartData *compactParts, int nCompactParts, int mype, char phase) {
#else
void DataManagerTransferRemoteChunk(CudaMultipoleMoments *moments, int nMoments, CompactPartData *compactParts, int nCompactParts) {
#endif

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
  buf->transferFromDevice = NO;
  buf->freeBuffer = NO;
  size = (nMoments) * sizeof(CudaMultipoleMoments);
  buf->size = size;
  buf->transferToDevice = (size > 0 ? YES : NO);
  
  //double mill = 1e6;
  //printf("DM remote moments: %f mbytes\n", 1.0*size/mill);
  
  if(size > 0){
#ifdef CUDA_USE_CUDAMALLOCHOST
#ifdef CUDA_MEMPOOL
    buf->hostBuffer = hapi_poolMalloc(size);
#else
    cudaMallocHost(&buf->hostBuffer, size);
#endif
#else
    buf->hostBuffer = malloc(size);
#endif
#ifdef CUDA_PRINT_ERRORS
    printf("DMRemote 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
    memcpy(buf->hostBuffer, moments, size);
  }
#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
  printf("DM REMOTE CHUNK node size: %d\n", size);
#endif

  buf = &(transferKernel.bufferInfo[REMOTE_PARTICLE_CORES_IDX]);
  buf->bufferID = REMOTE_PARTICLE_CORES;
  buf->transferFromDevice = NO;
  buf->freeBuffer = NO;
  size = (nCompactParts)*sizeof(CompactPartData);  
  buf->size = size;
  buf->transferToDevice = (size > 0 ? YES : NO);
  //printf("DM remote cores: %f mbytes\n", 1.0*size/mill);

  if(size > 0){
#ifdef CUDA_USE_CUDAMALLOCHOST
#ifdef CUDA_MEMPOOL
    buf->hostBuffer = hapi_poolMalloc(size);
#else
    cudaMallocHost(&buf->hostBuffer, size);
#endif
#else
    buf->hostBuffer = malloc(size);
#endif
#ifdef CUDA_PRINT_ERRORS
    printf("DMRemote 1: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
    memcpy(buf->hostBuffer, compactParts, size);
  }
#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
  printf("DM REMOTE CHUNK particle size: %d\n", size);
#endif

  transferKernel.callbackFn = 0;
  transferKernel.id = DM_TRANSFER_REMOTE_CHUNK;
#ifdef CUDA_INSTRUMENT_WRS
  transferKernel.chareIndex = mype;
  transferKernel.compType = transferKernel.id;
  transferKernel.compPhase = phase; 
#endif
  enqueue(wrQueue, &transferKernel);

}


void TreePieceCellListDataTransferLocal(CudaRequest *data){
	int numBlocks = data->numBucketsPlusOne-1;

	workRequest gravityKernel;
	//dataInfo *buffer, *partCoreBuffer;
        //ParameterStruct *pmtr;

	gravityKernel.dimGrid = dim3(numBlocks);
	gravityKernel.dimBlock = dim3(BLOCK_SIZE_GRAV);
	gravityKernel.smemSize = 0;

	gravityKernel.nBuffers = TP_GRAVITY_LOCAL_NBUFFERS;

	/* schedule buffers for transfer to the GPU */
	gravityKernel.bufferInfo = (dataInfo *) malloc(gravityKernel.nBuffers * sizeof(dataInfo));

	TreePieceCellListDataTransferBasic(data, &gravityKernel);

	gravityKernel.callbackFn = data->cb;
	gravityKernel.id = TP_GRAVITY_LOCAL;
#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
        printf("TRANSFER LOCAL CELL\n");
#endif
#ifdef CUDA_INSTRUMENT_WRS
        gravityKernel.chareIndex = data->tpIndex;
        gravityKernel.compType = gravityKernel.id;
        gravityKernel.compPhase = data->phase; 
#endif
	enqueue(wrQueue, &gravityKernel);
}

void TreePieceCellListDataTransferRemote(CudaRequest *data){
	int numBlocks = data->numBucketsPlusOne-1;

	workRequest gravityKernel;
	//dataInfo *buffer, *partCoreBuffer;

	gravityKernel.dimGrid = dim3(numBlocks);
	gravityKernel.dimBlock = dim3(BLOCK_SIZE_GRAV);
	gravityKernel.smemSize = 0;

	gravityKernel.nBuffers = TP_NODE_GRAVITY_REMOTE_NBUFFERS;

	/* schedule buffers for transfer to the GPU */
	gravityKernel.bufferInfo = (dataInfo *) malloc(gravityKernel.nBuffers * sizeof(dataInfo));

	TreePieceCellListDataTransferBasic(data, &gravityKernel);

	gravityKernel.callbackFn = data->cb;
	gravityKernel.id = TP_GRAVITY_REMOTE;
#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
        printf("TRANSFER REMOTE CELL\n");
#endif
#ifdef CUDA_INSTRUMENT_WRS
        gravityKernel.chareIndex = data->tpIndex;
        gravityKernel.compType = gravityKernel.id;
        gravityKernel.compPhase = data->phase; 
#endif
	enqueue(wrQueue, &gravityKernel);
}

void TreePieceCellListDataTransferRemoteResume(CudaRequest *data, CudaMultipoleMoments *missedMoments, int numMissedMoments){
  int numBlocks = data->numBucketsPlusOne-1;
  int size;

  workRequest gravityKernel;
  dataInfo *buffer;
  bool transfer;

  gravityKernel.dimGrid = dim3(numBlocks);
  gravityKernel.dimBlock = dim3(BLOCK_SIZE_GRAV);
  gravityKernel.smemSize = 0;

  gravityKernel.nBuffers = TP_NODE_GRAVITY_REMOTE_RESUME_NBUFFERS;

  /* schedule buffers for transfer to the GPU */
  gravityKernel.bufferInfo = (dataInfo *) malloc(gravityKernel.nBuffers * sizeof(dataInfo));

  TreePieceCellListDataTransferBasic(data, &gravityKernel);
  transfer = gravityKernel.bufferInfo[ILCELL_IDX].transferToDevice == YES;

  buffer = &(gravityKernel.bufferInfo[MISSED_MOMENTS_IDX]);
  buffer->bufferID = -1;
  size = (numMissedMoments) * sizeof(CudaMultipoleMoments);
  buffer->size = size;
  buffer->transferToDevice = transfer ? YES : NO;
  buffer->freeBuffer = transfer ? YES : NO;
  buffer->transferFromDevice = NO;

  if(transfer){
#ifdef CUDA_USE_CUDAMALLOCHOST
#ifdef CUDA_MEMPOOL
    buffer->hostBuffer = hapi_poolMalloc(size);
#else
    cudaMallocHost(&buffer->hostBuffer, size);
#endif
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
#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
        printf("TRANSFER REMOTE RESUME CELL\n");
#endif
#ifdef CUDA_INSTRUMENT_WRS
        gravityKernel.chareIndex = data->tpIndex;
        gravityKernel.compType = gravityKernel.id;
        gravityKernel.compPhase = data->phase; 
#endif
  enqueue(wrQueue, &gravityKernel);
}

void TreePieceCellListDataTransferBasic(CudaRequest *data, workRequest *gravityKernel){
	dataInfo *buffer;
	int numBucketsPlusOne = data->numBucketsPlusOne;
        int numBuckets = numBucketsPlusOne-1;
        int size;
        bool transfer;

	buffer = &(gravityKernel->bufferInfo[ILCELL_IDX]);
	buffer->bufferID = -1;
        size = (data->numInteractions) * sizeof(ILCell);
	buffer->size = size;

        transfer = size > 0;

	buffer->transferToDevice = transfer ? YES : NO;
	buffer->freeBuffer = transfer ? YES : NO;
	buffer->transferFromDevice = NO;
        buffer->hostBuffer = data->list;

	buffer = &(gravityKernel->bufferInfo[NODE_BUCKET_MARKERS_IDX]);
	buffer->bufferID = -1;
        size = (numBucketsPlusOne) * sizeof(int);
        buffer->size = transfer ? size : 0;
	buffer->transferToDevice = transfer ? YES : NO;
	buffer->freeBuffer = transfer ? YES : NO;
	buffer->transferFromDevice = NO;
        buffer->hostBuffer = data->bucketMarkers;

	buffer = &(gravityKernel->bufferInfo[NODE_BUCKET_START_MARKERS_IDX]);
	buffer->bufferID = -1;
	size = (numBuckets) * sizeof(int);
        buffer->size = size;
	buffer->transferToDevice = transfer ? YES : NO;
	buffer->freeBuffer = transfer ? YES : NO;
	buffer->transferFromDevice = NO;
        buffer->hostBuffer = data->bucketStarts;

	buffer = &(gravityKernel->bufferInfo[NODE_BUCKET_SIZES_IDX]);
	buffer->bufferID = -1;
        buffer->size = size;
	buffer->transferToDevice = transfer ? YES : NO;
	buffer->freeBuffer = transfer ? YES : NO;
	buffer->transferFromDevice = NO;
        buffer->hostBuffer = data->bucketSizes;

        ParameterStruct *ptr = (ParameterStruct *)malloc(sizeof(ParameterStruct));
        ptr->numInteractions = data->numInteractions;
        ptr->numBucketsPlusOne = numBucketsPlusOne;
        ptr->fperiod = data->fperiod;

        gravityKernel->userData = ptr;
}

void TreePiecePartListDataTransferLocalSmallPhase(CudaRequest *data, CompactPartData *particles, int len){
	int numBlocks = data->numBucketsPlusOne-1;

	workRequest gravityKernel;
        dataInfo *buffer;
        int size;
        ParameterStruct *ptr;
        bool transfer;

	gravityKernel.dimGrid = dim3(numBlocks);
	gravityKernel.dimBlock = dim3(BLOCK_SIZE_GRAV);
	gravityKernel.smemSize = 0;

	gravityKernel.nBuffers = TP_GRAVITY_LOCAL_NBUFFERS_SMALLPHASE;

	/* schedule buffers for transfer to the GPU */
	gravityKernel.bufferInfo = (dataInfo *) malloc(gravityKernel.nBuffers * sizeof(dataInfo));

	TreePiecePartListDataTransferBasic(data, &gravityKernel);
        transfer = gravityKernel.bufferInfo[ILPART_IDX].transferToDevice == YES;

	buffer = &(gravityKernel.bufferInfo[MISSED_PARTS_IDX]);
	buffer->bufferID = -1;
        size = (len) * sizeof(CompactPartData);
        buffer->size = size;
	buffer->transferToDevice = transfer ? YES : NO;
	buffer->freeBuffer = transfer ? YES : NO;
	buffer->transferFromDevice = NO;

        if(transfer){
#ifdef CUDA_USE_CUDAMALLOCHOST
#ifdef CUDA_MEMPOOL
          buffer->hostBuffer = hapi_poolMalloc(size);
#else
          cudaMallocHost(&buffer->hostBuffer, size);
#endif
#else
          buffer->hostBuffer = malloc(size);
#endif
#ifdef CUDA_PRINT_ERRORS
          printf("TPPartSmallPhase 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
          memcpy(buffer->hostBuffer, particles, size);
        }

        ptr = (ParameterStruct *)gravityKernel.userData;
        ptr->numMissedCores = len;

	gravityKernel.callbackFn = data->cb;
	gravityKernel.id = TP_PART_GRAVITY_LOCAL_SMALLPHASE;
#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
        printf("TRANSFER LOCAL SMALLPHASE PART\n");
#endif
#ifdef CUDA_INSTRUMENT_WRS
        gravityKernel.chareIndex = data->tpIndex;
        gravityKernel.compType = gravityKernel.id;
        gravityKernel.compPhase = data->phase; 
#endif
	enqueue(wrQueue, &gravityKernel);
}

void TreePiecePartListDataTransferLocal(CudaRequest *data){
	int numBlocks = data->numBucketsPlusOne-1;

	workRequest gravityKernel;
	//dataInfo *buffer, *partCoreBuffer;

	gravityKernel.dimGrid = dim3(numBlocks);
	gravityKernel.dimBlock = dim3(BLOCK_SIZE_GRAV);
	gravityKernel.smemSize = 0;

	gravityKernel.nBuffers = TP_GRAVITY_LOCAL_NBUFFERS;

	/* schedule buffers for transfer to the GPU */
	gravityKernel.bufferInfo = (dataInfo *) malloc(gravityKernel.nBuffers * sizeof(dataInfo));

	TreePiecePartListDataTransferBasic(data, &gravityKernel);

	gravityKernel.callbackFn = data->cb;
	gravityKernel.id = TP_PART_GRAVITY_LOCAL;
#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
        printf("TRANSFER LOCAL LARGEPHASE PART\n");
#endif
#ifdef CUDA_INSTRUMENT_WRS
        gravityKernel.chareIndex = data->tpIndex;
        gravityKernel.compType = gravityKernel.id;
        gravityKernel.compPhase = data->phase; 
#endif
	enqueue(wrQueue, &gravityKernel);
}

void TreePiecePartListDataTransferRemote(CudaRequest *data){
	int numBlocks = data->numBucketsPlusOne-1;

	workRequest gravityKernel;

	gravityKernel.dimGrid = dim3(numBlocks);
	gravityKernel.dimBlock = dim3(BLOCK_SIZE_GRAV);
	gravityKernel.smemSize = 0;

	gravityKernel.nBuffers = TP_PART_GRAVITY_REMOTE_NBUFFERS;

	gravityKernel.bufferInfo = (dataInfo *) malloc(gravityKernel.nBuffers * sizeof(dataInfo));

	TreePiecePartListDataTransferBasic(data, &gravityKernel);

	gravityKernel.callbackFn = data->cb;
	gravityKernel.id = TP_PART_GRAVITY_REMOTE;
#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
        printf("TRANSFER REMOTE PART\n");
#endif
#ifdef CUDA_INSTRUMENT_WRS
        gravityKernel.chareIndex = data->tpIndex;
        gravityKernel.compType = gravityKernel.id;
        gravityKernel.compPhase = data->phase; 
#endif
	enqueue(wrQueue, &gravityKernel);
}

void TreePiecePartListDataTransferRemoteResume(CudaRequest *data, CompactPartData *missedParts, int numMissedParts){
	int numBlocks = data->numBucketsPlusOne-1;
        int size;
        ParameterStruct *ptr;
        bool transfer;

	workRequest gravityKernel;
        dataInfo *buffer;

	gravityKernel.dimGrid = dim3(numBlocks);
	gravityKernel.dimBlock = dim3(BLOCK_SIZE_GRAV);
	gravityKernel.smemSize = 0;

	gravityKernel.nBuffers = TP_PART_GRAVITY_REMOTE_RESUME_NBUFFERS;

	/* schedule buffers for transfer to the GPU */
	gravityKernel.bufferInfo = (dataInfo *) malloc(gravityKernel.nBuffers * sizeof(dataInfo));

	TreePiecePartListDataTransferBasic(data, &gravityKernel);
        transfer = gravityKernel.bufferInfo[ILPART_IDX].transferToDevice == YES;

	buffer = &(gravityKernel.bufferInfo[MISSED_PARTS_IDX]);
	buffer->bufferID = -1;
        size = (numMissedParts) * sizeof(CompactPartData);
        buffer->size = size;
	buffer->transferToDevice = transfer ? YES : NO;
	buffer->freeBuffer = transfer ? YES : NO;
	buffer->transferFromDevice = NO;

        if(transfer){
#ifdef CUDA_USE_CUDAMALLOCHOST
#ifdef CUDA_MEMPOOL
          buffer->hostBuffer = hapi_poolMalloc(size);
#else
          cudaMallocHost(&buffer->hostBuffer, size);
#endif
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
#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
        printf("TRANSFER REMOTE RESUME PART\n");
#endif
#ifdef CUDA_INSTRUMENT_WRS
        gravityKernel.chareIndex = data->tpIndex;
        gravityKernel.compType = gravityKernel.id;
        gravityKernel.compPhase = data->phase; 
#endif
	enqueue(wrQueue, &gravityKernel);
}

void TreePiecePartListDataTransferBasic(CudaRequest *data, workRequest *gravityKernel){
	dataInfo *buffer;

	int numInteractions = data->numInteractions;
	int numBucketsPlusOne = data->numBucketsPlusOne;
        int numBuckets = numBucketsPlusOne-1;
        int size;
        bool transfer;

	buffer = &(gravityKernel->bufferInfo[ILPART_IDX]);
	buffer->bufferID = -1;
	size = (numInteractions) * sizeof(ILPart);
        buffer->size = size;

        transfer = size > 0;

	buffer->transferToDevice = transfer ? YES : NO;
	buffer->freeBuffer = transfer ? YES : NO;
	buffer->transferFromDevice = NO;

        buffer->hostBuffer = data->list;

	buffer = &(gravityKernel->bufferInfo[PART_BUCKET_MARKERS_IDX]);
	buffer->bufferID = -1;
	size = (numBucketsPlusOne) * sizeof(int);
        buffer->size = transfer ? size : 0;
	buffer->transferToDevice = transfer ? YES : NO;
	buffer->freeBuffer = transfer ? YES : NO;
	buffer->transferFromDevice = NO;
        buffer->hostBuffer = data->bucketMarkers;

	buffer = &(gravityKernel->bufferInfo[PART_BUCKET_START_MARKERS_IDX]);
	buffer->bufferID = -1;
	buffer->hostBuffer = data->bucketStarts;
	size = (numBuckets) * sizeof(int);
        buffer->size = size;
	buffer->transferToDevice = transfer ? YES : NO;
	buffer->freeBuffer = transfer ? YES : NO;
	buffer->transferFromDevice = NO;
        buffer->hostBuffer = data->bucketStarts;

	buffer = &(gravityKernel->bufferInfo[PART_BUCKET_SIZES_IDX]);
	buffer->bufferID = -1;
        buffer->size = size;
	buffer->transferToDevice = transfer ? YES : NO;
	buffer->freeBuffer = transfer ? YES : NO;
	buffer->transferFromDevice = NO;
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

#ifdef CUDA_INSTRUMENT_WRS
void FreeDataManagerLocalTreeMemory(bool freemom, bool freepart, int index, char phase){
#else
void FreeDataManagerLocalTreeMemory(bool freemom, bool freepart){
#endif
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
  buffer->freeBuffer = freemom ? YES : NO;
  buffer->hostBuffer = 0;
  buffer->size = 0;

  buffer = &(gravityKernel.bufferInfo[LOCAL_PARTICLE_CORES_IDX]);
  buffer->bufferID = LOCAL_PARTICLE_CORES;
  buffer->transferToDevice = NO ;
  buffer->transferFromDevice = NO;
  buffer->freeBuffer = freepart ? YES : NO;
  buffer->hostBuffer = 0;
  buffer->size = 0;

  gravityKernel.callbackFn = 0;
  gravityKernel.id = DM_TRANSFER_FREE_LOCAL;
  //printf("DM TRANSFER FREE LOCAL\n");
#ifdef CUDA_INSTRUMENT_WRS
  gravityKernel.chareIndex = index;
  gravityKernel.compType = gravityKernel.id;
  gravityKernel.compPhase = phase; 
#endif
  enqueue(wrQueue, &gravityKernel);

}

// this function begins the transfer of the next
// pending chunk on dm, once we are assured that
// the previous chunk's allocated memory has been
// freed
void initiateNextChunkTransfer(void *dm);

#ifdef CUDA_INSTRUMENT_WRS
void FreeDataManagerRemoteChunkMemory(int chunk, void *dm, bool freemom, bool freepart, int index, char phase){
#else
void FreeDataManagerRemoteChunkMemory(int chunk, void *dm, bool freemom, bool freepart){
#endif
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
  buffer->freeBuffer = (freemom ? YES : NO);
  buffer->size = 0;

  buffer = &(gravityKernel.bufferInfo[REMOTE_PARTICLE_CORES_IDX]);
  buffer->bufferID = REMOTE_PARTICLE_CORES;
  buffer->transferToDevice = NO ;
  buffer->transferFromDevice = NO;
  buffer->freeBuffer = (freepart ? YES : NO);
  buffer->size = 0;

  gravityKernel.callbackFn = 0;
  gravityKernel.id = DM_TRANSFER_FREE_REMOTE_CHUNK;

  // save a pointer to the data manager so that
  // the next chunk's transfer can be initiated once
  // the memory for this chunk has been freed
  gravityKernel.userData = dm;
  //printf("DM TRANSFER FREE REMOTE CHUNK\n");
#ifdef CUDA_INSTRUMENT_WRS
  gravityKernel.chareIndex = index;
  gravityKernel.compType = gravityKernel.id;
  gravityKernel.compPhase = phase; 
#endif
  enqueue(wrQueue, &gravityKernel);

}

#ifdef CUDA_INSTRUMENT_WRS
void TransferParticleVarsBack(VariablePartData *hostBuffer, int size, void *cb, int index, char phase){
#else
void TransferParticleVarsBack(VariablePartData *hostBuffer, int size, void *cb){
#endif
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
  buffer->hostBuffer = hostBuffer;
  buffer->size = size;
  buffer->transferFromDevice = size > 0 ? YES : NO;
  buffer->freeBuffer = size > 0 ? YES : NO;

  gravityKernel.callbackFn = cb;
  gravityKernel.id = DM_TRANSFER_BACK;
#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
  printf("DM TRANSFER BACK\n");
#endif
#ifdef CUDA_INSTRUMENT_WRS
  gravityKernel.chareIndex = index;
  gravityKernel.compType = gravityKernel.id;
  gravityKernel.compPhase = phase; 
#endif
  enqueue(wrQueue, &gravityKernel);
}

/*
void DummyKernel(void *cb){
  workRequest dummy;
  //dataInfo *buffer;

  dummy.dimGrid = dim3(1);
  dummy.dimBlock = dim3(BLOCK_SIZE_GRAV);
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

#ifdef CUDA_TRACE
  double st;
#endif
  ParameterStruct *ptr;
  switch (wr->id) {

    case DM_TRANSFER_LOCAL:
#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
      printf("DM_TRANSFER_LOCAL KERNELSELECT\n");
      printf("mom: 0x%x\n", devBuffers[LOCAL_MOMENTS]);
      printf("cores: 0x%x\n", devBuffers[LOCAL_PARTICLE_CORES]);
      printf("vars: 0x%x\n", devBuffers[LOCAL_PARTICLE_VARS]);
#endif
#ifdef CUDA_USE_CUDAMALLOCHOST
      if(wr->bufferInfo[LOCAL_MOMENTS_IDX].transferToDevice == YES){
#ifdef CUDA_MEMPOOL
        hapi_poolFree(wr->bufferInfo[LOCAL_MOMENTS_IDX].hostBuffer);
#else
        delayedFree(wr->bufferInfo[LOCAL_MOMENTS_IDX].hostBuffer);
#endif
      }
#ifdef CUDA_PRINT_ERRORS
      printf("DM_TRANSFER_LOCAL 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
      if(wr->bufferInfo[LOCAL_PARTICLE_CORES_IDX].transferToDevice == YES){
#ifdef CUDA_MEMPOOL
        hapi_poolFree(wr->bufferInfo[LOCAL_PARTICLE_CORES_IDX].hostBuffer);
#else
        delayedFree(wr->bufferInfo[LOCAL_PARTICLE_CORES_IDX].hostBuffer);
#endif
      }
#ifdef CUDA_PRINT_ERRORS
      printf("DM_TRANSFER_LOCAL 1: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
      if(wr->bufferInfo[LOCAL_PARTICLE_VARS_IDX].transferToDevice == YES){
#ifdef CUDA_MEMPOOL
        hapi_poolFree(wr->bufferInfo[LOCAL_PARTICLE_VARS_IDX].hostBuffer);
#else
        delayedFree(wr->bufferInfo[LOCAL_PARTICLE_VARS_IDX].hostBuffer);
#endif
      }
#ifdef CUDA_PRINT_ERRORS
      printf("DM_TRANSFER_LOCAL 2: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
#else
      if(wr->bufferInfo[LOCAL_MOMENTS_IDX].transferToDevice == YES){
        free(wr->bufferInfo[LOCAL_MOMENTS_IDX].hostBuffer);
      }
      if(wr->bufferInfo[LOCAL_PARTICLE_CORES_IDX].transferToDevice == YES){
        free(wr->bufferInfo[LOCAL_PARTICLE_CORES_IDX].hostBuffer);
      }
      if(wr->bufferInfo[LOCAL_PARTICLE_VARS_IDX].transferToDevice == YES){
        free(wr->bufferInfo[LOCAL_PARTICLE_VARS_IDX].hostBuffer);
      }
#endif
      break;

    case DM_TRANSFER_REMOTE_CHUNK:
#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
      printf("DM_TRANSFER_REMOTE_CHUNK, %d KERNELSELECT\n", wr->bufferInfo[REMOTE_MOMENTS_IDX].transferToDevice == YES);
#endif
      if(wr->bufferInfo[REMOTE_MOMENTS_IDX].transferToDevice == YES){
#ifdef CUDA_USE_CUDAMALLOCHOST
#ifdef CUDA_MEMPOOL
        hapi_poolFree(wr->bufferInfo[REMOTE_MOMENTS_IDX].hostBuffer);
#else
        delayedFree(wr->bufferInfo[REMOTE_MOMENTS_IDX].hostBuffer);
#endif
#ifdef CUDA_PRINT_ERRORS
        printf("DM_TRANSFER_REMOTE_CHUNK 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
#else
        free(wr->bufferInfo[REMOTE_MOMENTS_IDX].hostBuffer);
#endif
      }

      if(wr->bufferInfo[REMOTE_PARTICLE_CORES_IDX].transferToDevice == YES){
#ifdef CUDA_USE_CUDAMALLOCHOST
#ifdef CUDA_MEMPOOL
        hapi_poolFree(wr->bufferInfo[REMOTE_PARTICLE_CORES_IDX].hostBuffer);
#else
        delayedFree(wr->bufferInfo[REMOTE_PARTICLE_CORES_IDX].hostBuffer);
#endif
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
      printf("DM_TRANSFER_FREE_LOCAL KERNELSELECT\n");
      printf("buffers:\nlocal_particles: (0x%x)\nlocal_particle_vars: (0x%x)\nlocal_moments: (0x%x)\n", 
          devBuffers[LOCAL_PARTICLE_CORES],
          devBuffers[LOCAL_PARTICLE_VARS],
          devBuffers[LOCAL_MOMENTS]
          );

#endif
      break;

    case DM_TRANSFER_FREE_REMOTE_CHUNK:
#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
      printf("DM_TRANSFER_FREE_REMOTE_CHUNK KERNELSELECT\n");

#endif
      initiateNextChunkTransfer(wr->userData);
      break;

    case DM_TRANSFER_BACK:
#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
      printf("DM_TRANSFER_BACK: 0x%x KERNELSELECT\n", devBuffers[LOCAL_PARTICLE_VARS]);
#endif
      break;

    case TP_GRAVITY_LOCAL:
      ptr = (ParameterStruct *)wr->userData;
#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
      printf("TP_GRAVITY_LOCAL KERNELSELECT buffers:\nlocal_particles: (0x%x)\nlocal_particle_vars: (0x%x)\nlocal_moments: (0x%x)\nil_cell: %d (0x%x)\n", 
          devBuffers[LOCAL_PARTICLE_CORES],
          devBuffers[LOCAL_PARTICLE_VARS],
          devBuffers[LOCAL_MOMENTS],
          wr->bufferInfo[ILCELL_IDX].bufferID,
          devBuffers[wr->bufferInfo[ILCELL_IDX].bufferID]
          );
#endif
      if(wr->bufferInfo[ILCELL_IDX].transferToDevice == YES){
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

#ifdef CUDA_TRACE
    st = CmiWallTimer();
#endif
#ifdef CUDA_USE_CUDAMALLOCHOST
#ifdef CUDA_MEMPOOL
        hapi_poolFree((ILCell *)wr->bufferInfo[ILCELL_IDX].hostBuffer);
        hapi_poolFree((int *)wr->bufferInfo[NODE_BUCKET_MARKERS_IDX].hostBuffer);
        hapi_poolFree((int *)wr->bufferInfo[NODE_BUCKET_START_MARKERS_IDX].hostBuffer);
        hapi_poolFree((int *)wr->bufferInfo[NODE_BUCKET_SIZES_IDX].hostBuffer);
#else
        delayedFree((ILCell *)wr->bufferInfo[ILCELL_IDX].hostBuffer);
        delayedFree((int *)wr->bufferInfo[NODE_BUCKET_MARKERS_IDX].hostBuffer);
        delayedFree((int *)wr->bufferInfo[NODE_BUCKET_START_MARKERS_IDX].hostBuffer);
        delayedFree((int *)wr->bufferInfo[NODE_BUCKET_SIZES_IDX].hostBuffer);
#endif
#ifdef CUDA_PRINT_ERRORS
        printf("TP_GRAVITY_LOCAL 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
#else
        free((ILCell *)wr->bufferInfo[ILCELL_IDX].hostBuffer);
        free((int *)wr->bufferInfo[NODE_BUCKET_MARKERS_IDX].hostBuffer);
        free((int *)wr->bufferInfo[NODE_BUCKET_START_MARKERS_IDX].hostBuffer);
        free((int *)wr->bufferInfo[NODE_BUCKET_SIZES_IDX].hostBuffer);
#endif
#ifdef CUDA_TRACE
      traceUserBracketEvent(CUDA_LOCAL_NODE_KERNEL, st, CmiWallTimer()); 
#endif
      }
      free((ParameterStruct *)ptr);

      break;

    case TP_PART_GRAVITY_LOCAL_SMALLPHASE:

      ptr = (ParameterStruct *)wr->userData;
#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
      printf("TP_PART_GRAVITY_LOCAL_SMALLPHASE KERNELSELECT buffers:\nlocal_particles: (0x%x)\nlocal_particle_vars: (0x%x)\nil_cell: %d (0x%x)\n", 
          devBuffers[LOCAL_PARTICLE_CORES],
          devBuffers[LOCAL_PARTICLE_VARS],
          wr->bufferInfo[ILPART_IDX].bufferID,
          devBuffers[wr->bufferInfo[ILPART_IDX].bufferID]
          );
#endif
      if(wr->bufferInfo[ILPART_IDX].transferToDevice == YES){
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
           ptr->fperiod,3 
          );
#endif

#ifdef CUDA_TRACE
    st = CmiWallTimer();
#endif
#ifdef CUDA_USE_CUDAMALLOCHOST
#ifdef CUDA_MEMPOOL
        hapi_poolFree((CompactPartData *)wr->bufferInfo[MISSED_PARTS_IDX].hostBuffer);
        hapi_poolFree((ILPart *)wr->bufferInfo[ILPART_IDX].hostBuffer);
        hapi_poolFree((int *)wr->bufferInfo[PART_BUCKET_MARKERS_IDX].hostBuffer);
        hapi_poolFree((int *)wr->bufferInfo[PART_BUCKET_START_MARKERS_IDX].hostBuffer);
        hapi_poolFree((int *)wr->bufferInfo[PART_BUCKET_SIZES_IDX].hostBuffer);
#else
        delayedFree((CompactPartData *)wr->bufferInfo[MISSED_PARTS_IDX].hostBuffer);
        delayedFree((ILPart *)wr->bufferInfo[ILPART_IDX].hostBuffer);
        delayedFree((int *)wr->bufferInfo[PART_BUCKET_MARKERS_IDX].hostBuffer);
        delayedFree((int *)wr->bufferInfo[PART_BUCKET_START_MARKERS_IDX].hostBuffer);
        delayedFree((int *)wr->bufferInfo[PART_BUCKET_SIZES_IDX].hostBuffer);
#endif
#ifdef CUDA_PRINT_ERRORS
        printf("TP_PART_GRAVITY_LOCAL_SMALLPHASE 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
#else
        free((CompactPartData *)wr->bufferInfo[MISSED_PARTS_IDX].hostBuffer);
        free((ILPart *)wr->bufferInfo[ILPART_IDX].hostBuffer);
        free((int *)wr->bufferInfo[PART_BUCKET_MARKERS_IDX].hostBuffer);
        free((int *)wr->bufferInfo[PART_BUCKET_START_MARKERS_IDX].hostBuffer);
        free((int *)wr->bufferInfo[PART_BUCKET_SIZES_IDX].hostBuffer);
#endif
#ifdef CUDA_TRACE
      traceUserBracketEvent(CUDA_LOCAL_PART_KERNEL, st, CmiWallTimer()); 
#endif
      }
      free((ParameterStruct *)ptr);

      break;


    case TP_PART_GRAVITY_LOCAL:

      ptr = (ParameterStruct *)wr->userData;
#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
      printf("TP_PART_GRAVITY_LOCAL KERNELSELECT buffers:\nlocal_particles: (0x%x)\nlocal_particle_vars: (0x%x)\nil_cell: %d (0x%x)\n", 
          devBuffers[LOCAL_PARTICLE_CORES],
          devBuffers[LOCAL_PARTICLE_VARS],
          wr->bufferInfo[ILPART_IDX].bufferID,
          devBuffers[wr->bufferInfo[ILPART_IDX].bufferID]
          );
#endif
      if(wr->bufferInfo[ILPART_IDX].transferToDevice == YES){
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

#ifdef CUDA_TRACE
    st = CmiWallTimer();
#endif
#ifdef CUDA_USE_CUDAMALLOCHOST
#ifdef CUDA_MEMPOOL
        hapi_poolFree((ILPart *)wr->bufferInfo[ILPART_IDX].hostBuffer);
        hapi_poolFree((int *)wr->bufferInfo[PART_BUCKET_MARKERS_IDX].hostBuffer);
        hapi_poolFree((int *)wr->bufferInfo[PART_BUCKET_START_MARKERS_IDX].hostBuffer);
        hapi_poolFree((int *)wr->bufferInfo[PART_BUCKET_SIZES_IDX].hostBuffer);
#else
        delayedFree((ILPart *)wr->bufferInfo[ILPART_IDX].hostBuffer);
        delayedFree((int *)wr->bufferInfo[PART_BUCKET_MARKERS_IDX].hostBuffer);
        delayedFree((int *)wr->bufferInfo[PART_BUCKET_START_MARKERS_IDX].hostBuffer);
        delayedFree((int *)wr->bufferInfo[PART_BUCKET_SIZES_IDX].hostBuffer);
#endif
#ifdef CUDA_PRINT_ERRORS
        printf("TP_PART_GRAVITY_LOCAL 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
#else
        free((ILPart *)wr->bufferInfo[ILPART_IDX].hostBuffer);
        free((int *)wr->bufferInfo[PART_BUCKET_MARKERS_IDX].hostBuffer);
        free((int *)wr->bufferInfo[PART_BUCKET_START_MARKERS_IDX].hostBuffer);
        free((int *)wr->bufferInfo[PART_BUCKET_SIZES_IDX].hostBuffer);
#endif
#ifdef CUDA_TRACE
      traceUserBracketEvent(CUDA_LOCAL_PART_KERNEL, st, CmiWallTimer()); 
#endif
      }
      free((ParameterStruct *)ptr);

      break;

    case TP_GRAVITY_REMOTE:

      ptr = (ParameterStruct *)wr->userData;
#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
      printf("TP_GRAVITY_REMOTE KERNELSELECT buffers:\nlocal_particles: (0x%x)\nlocal_particle_vars: (0x%x)\nremote_moments: (0x%x)\nil_cell: %d (0x%x)\n", 
          devBuffers[LOCAL_PARTICLE_CORES],
          devBuffers[LOCAL_PARTICLE_VARS],
          devBuffers[REMOTE_MOMENTS],
          wr->bufferInfo[ILCELL_IDX].bufferID,
          devBuffers[wr->bufferInfo[ILCELL_IDX].bufferID]
          );
#endif
      if(wr->bufferInfo[ILCELL_IDX].transferToDevice == YES){
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

#ifdef CUDA_TRACE
    st = CmiWallTimer();
#endif
#ifdef CUDA_USE_CUDAMALLOCHOST
#ifdef CUDA_MEMPOOL
        hapi_poolFree((ILCell *)wr->bufferInfo[ILCELL_IDX].hostBuffer);
        hapi_poolFree((int *)wr->bufferInfo[NODE_BUCKET_MARKERS_IDX].hostBuffer);
        hapi_poolFree((int *)wr->bufferInfo[NODE_BUCKET_START_MARKERS_IDX].hostBuffer);
        hapi_poolFree((int *)wr->bufferInfo[NODE_BUCKET_SIZES_IDX].hostBuffer);
#else
        delayedFree((ILCell *)wr->bufferInfo[ILCELL_IDX].hostBuffer);
        delayedFree((int *)wr->bufferInfo[NODE_BUCKET_MARKERS_IDX].hostBuffer);
        delayedFree((int *)wr->bufferInfo[NODE_BUCKET_START_MARKERS_IDX].hostBuffer);
        delayedFree((int *)wr->bufferInfo[NODE_BUCKET_SIZES_IDX].hostBuffer);
#endif
#ifdef CUDA_PRINT_ERRORS
        printf("TP_GRAVITY_REMOTE 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
#else
        free((ILCell *)wr->bufferInfo[ILCELL_IDX].hostBuffer);
        free((int *)wr->bufferInfo[NODE_BUCKET_MARKERS_IDX].hostBuffer);
        free((int *)wr->bufferInfo[NODE_BUCKET_START_MARKERS_IDX].hostBuffer);
        free((int *)wr->bufferInfo[NODE_BUCKET_SIZES_IDX].hostBuffer);
#endif
#ifdef CUDA_TRACE
      traceUserBracketEvent(CUDA_REMOTE_NODE_KERNEL, st, CmiWallTimer()); 
#endif
      }
      free((ParameterStruct *)ptr);

      break;

    case TP_PART_GRAVITY_REMOTE:

      ptr = (ParameterStruct *)wr->userData;
#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
      printf("TP_PART_GRAVITY_REMOTE KERNELSELECT buffers:\nlocal_particles: (0x%x)\nlocal_particle_vars: (0x%x)\nil_cell: %d (0x%x)\n", 
          devBuffers[LOCAL_PARTICLE_CORES],
          devBuffers[LOCAL_PARTICLE_VARS],
          wr->bufferInfo[ILPART_IDX].bufferID,
          devBuffers[wr->bufferInfo[ILPART_IDX].bufferID]
          );
#endif
      if(wr->bufferInfo[ILPART_IDX].transferToDevice == YES){
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

#ifdef CUDA_TRACE
    st = CmiWallTimer();
#endif
#ifdef CUDA_USE_CUDAMALLOCHOST
#ifdef CUDA_MEMPOOL
        hapi_poolFree((ILPart *)wr->bufferInfo[ILPART_IDX].hostBuffer);
        hapi_poolFree((int *)wr->bufferInfo[PART_BUCKET_MARKERS_IDX].hostBuffer);
        hapi_poolFree((int *)wr->bufferInfo[PART_BUCKET_START_MARKERS_IDX].hostBuffer);
        hapi_poolFree((int *)wr->bufferInfo[PART_BUCKET_SIZES_IDX].hostBuffer);
#else
        delayedFree((ILPart *)wr->bufferInfo[ILPART_IDX].hostBuffer);
        delayedFree((int *)wr->bufferInfo[PART_BUCKET_MARKERS_IDX].hostBuffer);
        delayedFree((int *)wr->bufferInfo[PART_BUCKET_START_MARKERS_IDX].hostBuffer);
        delayedFree((int *)wr->bufferInfo[PART_BUCKET_SIZES_IDX].hostBuffer);
#endif
#ifdef CUDA_PRINT_ERRORS
        printf("TP_PART_GRAVITY_REMOTE 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
#else
        free((ILPart *)wr->bufferInfo[ILPART_IDX].hostBuffer);
        free((int *)wr->bufferInfo[PART_BUCKET_MARKERS_IDX].hostBuffer);
        free((int *)wr->bufferInfo[PART_BUCKET_START_MARKERS_IDX].hostBuffer);
        free((int *)wr->bufferInfo[PART_BUCKET_SIZES_IDX].hostBuffer);
#endif
#ifdef CUDA_TRACE
      traceUserBracketEvent(CUDA_REMOTE_PART_KERNEL, st, CmiWallTimer()); 
#endif
      }
      free((ParameterStruct *)ptr);

      break;

    case TP_GRAVITY_REMOTE_RESUME:

      ptr = (ParameterStruct *)wr->userData;
#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
      printf("TP_GRAVITY_REMOTE_RESUME KERNELSELECT buffers:\nlocal_particles: (0x%x)\nlocal_particle_vars: (0x%x)\nmissed_moments: %d (0x%x)\nil_cell: %d (0x%x)\n", 
          devBuffers[LOCAL_PARTICLE_CORES],
          devBuffers[LOCAL_PARTICLE_VARS],
          wr->bufferInfo[MISSED_MOMENTS_IDX].bufferID,
          devBuffers[wr->bufferInfo[MISSED_MOMENTS_IDX].bufferID],
          wr->bufferInfo[ILCELL_IDX].bufferID,
          devBuffers[wr->bufferInfo[ILCELL_IDX].bufferID]
          );
#endif
      if(wr->bufferInfo[ILCELL_IDX].transferToDevice == YES){
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

#ifdef CUDA_TRACE
    st = CmiWallTimer();
#endif
#ifdef CUDA_USE_CUDAMALLOCHOST
#ifdef CUDA_MEMPOOL
        hapi_poolFree((CudaMultipoleMoments *)wr->bufferInfo[MISSED_MOMENTS_IDX].hostBuffer);
        hapi_poolFree((ILCell *)wr->bufferInfo[ILCELL_IDX].hostBuffer);
        hapi_poolFree((int *)wr->bufferInfo[NODE_BUCKET_MARKERS_IDX].hostBuffer);
        hapi_poolFree((int *)wr->bufferInfo[NODE_BUCKET_START_MARKERS_IDX].hostBuffer);
        hapi_poolFree((int *)wr->bufferInfo[NODE_BUCKET_SIZES_IDX].hostBuffer);
#else
        delayedFree((CudaMultipoleMoments *)wr->bufferInfo[MISSED_MOMENTS_IDX].hostBuffer);
        delayedFree((ILCell *)wr->bufferInfo[ILCELL_IDX].hostBuffer);
        delayedFree((int *)wr->bufferInfo[NODE_BUCKET_MARKERS_IDX].hostBuffer);
        delayedFree((int *)wr->bufferInfo[NODE_BUCKET_START_MARKERS_IDX].hostBuffer);
        delayedFree((int *)wr->bufferInfo[NODE_BUCKET_SIZES_IDX].hostBuffer);
#endif
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
#ifdef CUDA_TRACE
      traceUserBracketEvent(CUDA_REMOTE_RESUME_NODE_KERNEL, st, CmiWallTimer()); 
#endif
      }
      free((ParameterStruct *)ptr);

      break;

    case TP_PART_GRAVITY_REMOTE_RESUME:

      ptr = (ParameterStruct *)wr->userData;
#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
      printf("TP_PART_GRAVITY_REMOTE_RESUME KERNELSELECT buffers:\nlocal_particles: (0x%x)\nlocal_particle_vars: (0x%x)\nmissed_parts: %d (0x%x)\nil_cell: %d (0x%x)\n", 
          devBuffers[LOCAL_PARTICLE_CORES],
          devBuffers[LOCAL_PARTICLE_VARS],
          wr->bufferInfo[MISSED_PARTS_IDX].bufferID,
          devBuffers[wr->bufferInfo[MISSED_PARTS_IDX].bufferID],
          wr->bufferInfo[ILPART_IDX].bufferID,
          devBuffers[wr->bufferInfo[ILPART_IDX].bufferID]
          );
#endif
      if(wr->bufferInfo[ILPART_IDX].transferToDevice == YES){
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

#ifdef CUDA_TRACE
    st = CmiWallTimer();
#endif
#ifdef CUDA_USE_CUDAMALLOCHOST
#ifdef CUDA_MEMPOOL
        hapi_poolFree((CompactPartData *)wr->bufferInfo[MISSED_PARTS_IDX].hostBuffer);
        hapi_poolFree((ILPart *)wr->bufferInfo[ILPART_IDX].hostBuffer);
        hapi_poolFree((int *)wr->bufferInfo[PART_BUCKET_MARKERS_IDX].hostBuffer);
        hapi_poolFree((int *)wr->bufferInfo[PART_BUCKET_START_MARKERS_IDX].hostBuffer);
        hapi_poolFree((int *)wr->bufferInfo[PART_BUCKET_SIZES_IDX].hostBuffer);
#else
        delayedFree((CompactPartData *)wr->bufferInfo[MISSED_PARTS_IDX].hostBuffer);
        delayedFree((ILPart *)wr->bufferInfo[ILPART_IDX].hostBuffer);
        delayedFree((int *)wr->bufferInfo[PART_BUCKET_MARKERS_IDX].hostBuffer);
        delayedFree((int *)wr->bufferInfo[PART_BUCKET_START_MARKERS_IDX].hostBuffer);
        delayedFree((int *)wr->bufferInfo[PART_BUCKET_SIZES_IDX].hostBuffer);
#endif
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
#ifdef CUDA_TRACE
      traceUserBracketEvent(CUDA_REMOTE_RESUME_PART_KERNEL, st, CmiWallTimer()); 
#endif
      }
      free((ParameterStruct *)ptr);

      break;

    case TOP_EWALD_KERNEL:
      cudaMemcpyToSymbol(cachedData, wr->bufferInfo[EWALD_READ_ONLY_DATA].hostBuffer, sizeof(EwaldReadOnlyData), 0, cudaMemcpyHostToDevice);
      EwaldTopKernel<<<wr->dimGrid, wr->dimBlock, wr->smemSize, kernel_stream>>>
        ((GravityParticleData*)devBuffers[wr->bufferInfo[PARTICLE_TABLE].bufferID]);
#ifdef CUDA_PRINT_ERRORS
      printf("TOP_EWALD_KERNEL: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
      

      break; 
    case BOTTOM_EWALD_KERNEL:
      cudaMemcpyToSymbol(ewt, wr->bufferInfo[EWALD_TABLE].hostBuffer, NEWH * sizeof(EwtData), 0, cudaMemcpyHostToDevice);
      EwaldBottomKernel<<<wr->dimGrid, wr->dimBlock, 
        wr->smemSize, kernel_stream>>>
          ((GravityParticleData
            *)devBuffers[wr->bufferInfo[PARTICLE_TABLE].bufferID]);
#ifdef CUDA_PRINT_ERRORS
      printf("BOTTOM_EWALD_KERNEL : %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
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

__global__ void EwaldTopKernel(GravityParticleData *particleTable);
__global__ void EwaldBottomKernel(GravityParticleData *particleTable);

extern unsigned int timerHandle; 

void EwaldHostMemorySetup(EwaldData *h_idata, int nParticles, int nEwhLoop,
    void *cb) {
/*
  pinnedMemReq reqs; 

  int nBuffers = 3; 
  size_t *sizes = (size_t *)malloc(nBuffers * sizeof(size_t)); 
  void ***hostPtrs = (void ***) malloc(nBuffers * sizeof(void **)); 
  hostPtrs[0] = (void **) &(h_idata->p); 
  hostPtrs[1] = (void **) &(h_idata->ewt); 
  hostPtrs[2] = (void **) &(h_idata->cachedData); 
  sizes[0] = (nParticles+1) * sizeof(GravityParticleData);
  sizes[1] = nEwhLoop * sizeof(EwtData);
  sizes[2] = sizeof(EwaldReadOnlyData);


  reqs.sizes = sizes; 
  reqs.hostPtrs = hostPtrs; 
  reqs.callbackFn = cb;
  reqs.nBuffers = nBuffers; 

  pinnedMallocHost(&reqs); 
*/
#ifdef CUDA_MEMPOOL
  h_idata->p = (GravityParticleData *) hapi_poolMalloc((nParticles+1)*sizeof(GravityParticleData));
  h_idata->ewt = (EwtData *) hapi_poolMalloc((nEwhLoop)*sizeof(EwtData));
  h_idata->cachedData = (EwaldReadOnlyData *) hapi_poolMalloc(sizeof(EwaldReadOnlyData));
#else
  cudaMallocHost((void**)&(h_idata->p), (nParticles+1) * sizeof(GravityParticleData));
  cudaMallocHost((void**)&(h_idata->ewt), nEwhLoop * sizeof(EwtData));
  cudaMallocHost((void**)&(h_idata->cachedData), sizeof(EwaldReadOnlyData));
#endif
}

void EwaldHostMemoryFree(EwaldData *h_idata) {
#ifdef CUDA_MEMPOOL
  hapi_poolFree(h_idata->p); 
  hapi_poolFree(h_idata->ewt); 
  hapi_poolFree(h_idata->cachedData); 
#else
  cudaFreeHost(h_idata->p); 
  cudaFreeHost(h_idata->ewt); 
  cudaFreeHost(h_idata->cachedData); 
#endif
}

#ifdef CUDA_INSTRUMENT_WRS
void EwaldHost(EwaldData *h_idata, void *cb, int myIndex, char phase) {
#else
void EwaldHost(EwaldData *h_idata, void *cb, int myIndex) {
#endif

  int n = h_idata->cachedData->n;
  int numBlocks = (int) ceilf((float)n/BLOCK_SIZE);
  int nEwhLoop = h_idata->cachedData->nEwhLoop;

  workRequest topKernel;
  dataInfo *particleTableInfo, *cachedDataInfo, *ewtInfo;  

  topKernel.dimGrid = dim3(numBlocks); 
  topKernel.dimBlock = dim3(BLOCK_SIZE); 
  topKernel.smemSize = 0; 

  topKernel.nBuffers = 2; 

  /* schedule two buffers for transfer to the GPU */ 
  topKernel.bufferInfo = 
    (dataInfo *) malloc(topKernel.nBuffers * sizeof(dataInfo));

  particleTableInfo = &(topKernel.bufferInfo[0]);
  particleTableInfo->bufferID = NUM_GRAVITY_BUFS + BUFFERS_PER_CHARE * myIndex + PARTICLE_TABLE; 
  particleTableInfo->transferToDevice = YES; 
  particleTableInfo->transferFromDevice = NO; 
  particleTableInfo->freeBuffer = NO; 
  particleTableInfo->hostBuffer = h_idata->p; 
  particleTableInfo->size = (n+1) * sizeof(GravityParticleData); 

  cachedDataInfo = &(topKernel.bufferInfo[1]); 
  cachedDataInfo->bufferID = NUM_GRAVITY_BUFS + BUFFERS_PER_CHARE * myIndex + EWALD_READ_ONLY_DATA;
  cachedDataInfo->transferToDevice = NO; 
  cachedDataInfo->transferFromDevice = NO; 
  cachedDataInfo->freeBuffer = NO; 
  cachedDataInfo->hostBuffer = h_idata->cachedData; 
  cachedDataInfo->size = sizeof(EwaldReadOnlyData); 

  topKernel.callbackFn = NULL;
  topKernel.id = TOP_EWALD_KERNEL; 
#ifdef CUDA_INSTRUMENT_WRS
  topKernel.chareIndex = myIndex;
  topKernel.compType = topKernel.id;
  topKernel.compPhase = phase; 
#endif
  enqueue(wrQueue, &topKernel); 
#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
  printf("[%d] in EwaldHost, enqueued TopKernel\n", myIndex);
#endif

  workRequest bottomKernel; 

  bottomKernel.dimGrid = dim3(numBlocks); 
  bottomKernel.dimBlock = dim3(BLOCK_SIZE); 
  bottomKernel.smemSize = 0; 

  bottomKernel.nBuffers = 3; 

  bottomKernel.bufferInfo = 
    (dataInfo *) malloc(bottomKernel.nBuffers * sizeof(dataInfo));

  particleTableInfo = &(bottomKernel.bufferInfo[0]);
  particleTableInfo->bufferID = NUM_GRAVITY_BUFS + BUFFERS_PER_CHARE * myIndex + PARTICLE_TABLE; 
  particleTableInfo->transferToDevice = NO; 
  particleTableInfo->transferFromDevice = YES; 
  particleTableInfo->freeBuffer = YES; 
  particleTableInfo->hostBuffer = h_idata->p; 
  particleTableInfo->size = (n+1) * sizeof(GravityParticleData); 

  cachedDataInfo = &(bottomKernel.bufferInfo[1]); 
  cachedDataInfo->bufferID = NUM_GRAVITY_BUFS + BUFFERS_PER_CHARE * myIndex + EWALD_READ_ONLY_DATA; 
  cachedDataInfo->transferToDevice = NO; 
  cachedDataInfo->transferFromDevice = NO; 
  cachedDataInfo->freeBuffer = NO; 
  cachedDataInfo->hostBuffer = h_idata->cachedData; 
  cachedDataInfo->size = sizeof(EwaldReadOnlyData); 

  ewtInfo = &(bottomKernel.bufferInfo[2]); 
  ewtInfo->bufferID = NUM_GRAVITY_BUFS + BUFFERS_PER_CHARE * myIndex + EWALD_TABLE; 
  ewtInfo->transferToDevice = NO; 
  ewtInfo->transferFromDevice = NO; 
  ewtInfo->freeBuffer = NO; 
  ewtInfo->hostBuffer = h_idata->ewt; 
  ewtInfo->size = nEwhLoop * sizeof(EwtData); 

  bottomKernel.callbackFn = cb;
  bottomKernel.id = BOTTOM_EWALD_KERNEL; 
#ifdef CUDA_INSTRUMENT_WRS
  bottomKernel.chareIndex = myIndex;
  bottomKernel.compType = bottomKernel.id;
  bottomKernel.compPhase = phase; 
#endif
  enqueue(wrQueue, &bottomKernel); 
#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
  printf("[%d] in EwaldHost, enqueued BotKernel\n", myIndex);
#endif
}

__global__ void EwaldTopKernel(GravityParticleData *particleTable) {

  GravityParticleData *p;  
  MultipoleMomentsData *mom; 

  float alphan;
  float fPot, ax, ay, az;
  float x, y, z, r2, dir, dir2, a; 
  float xdif, ydif, zdif; 
  float g0, g1, g2, g3; 
  float Q2, Q2mirx, Q2miry, Q2mirz, Q2mir, Qta; 
  int id, ix, iy, iz, bInHole, bInHolex, bInHolexy;

  mom = &(cachedData->mm); 
  Q2 = 0.5 * (mom->xx + mom->yy + mom->zz); 


  id = blockIdx.x * BLOCK_SIZE + threadIdx.x + 1;
  if (id > cachedData->n) {
    return;
  }
  p = &(particleTable[id]);


#ifdef DEBUG
  if (blockIdx.x == 0 && threadIdx.x == 0) {
    printf("Moments\n");
    printf("xx %f, xy %f, xz %f, yy %f, yz %f, zz %f\n", mom->xx, mom->xy, mom->xz,
        mom->yy, mom->yz, mom->zz);
  }
#endif

  // if (p->rung < activeRung) return; // only for multistepping

  ax = 0.0f;
  ay = 0.0f;
  az = 0.0f;

  xdif = p->position_x - mom->cmx; 
  ydif = p->position_y - mom->cmy; 
  zdif = p->position_z - mom->cmz;
  fPot = mom->totalMass*cachedData->k1;
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
        }

        Q2mirx = mom->xx*x + mom->xy*y + mom->xz*z;
        Q2miry = mom->xy*x + mom->yy*y + mom->yz*z;
        Q2mirz = mom->xz*x + mom->yz*y + mom->zz*z;
        Q2mir = 0.5*(Q2mirx*x + Q2miry*y + Q2mirz*z);
        Qta = g1*mom->totalMass - g2*Q2 + g3*Q2mir;
        fPot -= g0*mom->totalMass - g1*Q2 + g2*Q2mir;

        ax += g2*(Q2mirx) - x*Qta;
        ay += g2*(Q2miry) - y*Qta;
        az += g2*(Q2mirz) - z*Qta;
      }
    }
  }

  p->potential += fPot;
  p->acceleration_x += ax;
  p->acceleration_y += ay;
  p->acceleration_z += az;

}

__global__ void EwaldBottomKernel(GravityParticleData *particleTable) {

  int i, id;
  float hdotx, c, s, fPot; 
  float ax, ay, az; 
  float tempEwt; 
  float xdif, ydif, zdif; 
  GravityParticleData *p; 
  MultipoleMomentsData *mom; 

  mom = &(cachedData->mm); 

  id = blockIdx.x * BLOCK_SIZE + threadIdx.x + 1;
  if (id > cachedData->n) {
    return;
  }
  p = &(particleTable[id]);

  /*
   ** Scoring for the h-loop (+,*)
   **        Without trig = (10,14)
   **          Trig est.    = 2*(6,11)  same as 1/sqrt scoring.
   **            Total        = (22,36)
   **                                   = 58
   */

  fPot=0.0f;
  ax=0.0f;
  ay=0.0f;
  az=0.0f;

  xdif = p->position_x - mom->cmx; 
  ydif = p->position_y - mom->cmy; 
  zdif = p->position_z - mom->cmz; 

  for (i=0;i<cachedData->nEwhLoop;++i) {
    hdotx = ewt[i].hx * xdif + ewt[i].hy * ydif + ewt[i].hz * zdif;
    c = cosf(hdotx);
    s = sinf(hdotx);    
    fPot += ewt[i].hCfac*c + ewt[i].hSfac*s;    
    tempEwt = ewt[i].hCfac*s - ewt[i].hSfac*c;
    ax += ewt[i].hx * tempEwt;
    ay += ewt[i].hy * tempEwt;
    az += ewt[i].hz * tempEwt;
  }

  p->potential += fPot;
  p->acceleration_x += ax;
  p->acceleration_y += ay;
  p->acceleration_z += az;

}
