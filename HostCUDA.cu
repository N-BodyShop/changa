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
#include <assert.h>

#include "CudaFunctions.h"
#ifdef HEXADECAPOLE
# include "CUDAMoments.cu"
#endif
#include "HostCUDA.h"
#include "EwaldCUDA.h"

#include "wr.h"

#define CmiMyPe() _Cmi_mype
extern int _Cmi_mype;

//extern workRequestQueue *wrQueue;
//extern void **devBuffers;
//extern cudaStream_t kernel_stream;

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
void DataManagerTransferLocalTree(CudaMultipoleMoments *moments, int nMoments, CompactPartData *compactParts, int nCompactParts, int mype, char phase, void *wrCallback) {
#else
  void DataManagerTransferLocalTree(CudaMultipoleMoments *moments, int nMoments, CompactPartData *compactParts, int nCompactParts, int mype, void *wrCallback) {
#endif

	workRequest transferKernel;
	dataInfo *buf;
        int size;

        // operation will not invoke a kernel
	transferKernel.dimGrid = dim3(0);
	transferKernel.dimBlock = dim3(0);
	transferKernel.smemSize = 0;

	transferKernel.nBuffers = DM_TRANSFER_LOCAL_NBUFFERS;

	/* schedule buffers for transfer to the GPU */
	transferKernel.bufferInfo = (dataInfo *) malloc(transferKernel.nBuffers * sizeof(dataInfo));

	buf = &(transferKernel.bufferInfo[LOCAL_MOMENTS_IDX]);
	buf->bufferID = LOCAL_MOMENTS;
	buf->freeBuffer = false;
        size = (nMoments) * sizeof(CudaMultipoleMoments);
	buf->size = size;
	buf->transferToDevice = size > 0;
	buf->transferFromDevice = false;

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
        else{
          buf->hostBuffer = NULL;
        }

#ifdef CUDA_PRINT_ERRORS
        printf("(%d) DMLocal 0: %s hostBuf: 0x%x, size: %d\n", mype, cudaGetErrorString( cudaGetLastError() ), buf->hostBuffer, size );
#endif
        memcpy(buf->hostBuffer, moments, size);

	buf = &(transferKernel.bufferInfo[LOCAL_PARTICLE_CORES_IDX]);
	buf->bufferID = LOCAL_PARTICLE_CORES;
	buf->freeBuffer = false;
        size = (nCompactParts)*sizeof(CompactPartData);
        buf->size = size;
	buf->transferToDevice = size > 0;
	buf->transferFromDevice = false;

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
        else{
          buf->hostBuffer = NULL;
        }

        VariablePartData *zeroArray;

	buf = &(transferKernel.bufferInfo[LOCAL_PARTICLE_VARS_IDX]);
	buf->bufferID = LOCAL_PARTICLE_VARS;
	buf->freeBuffer = false;
        size = (nCompactParts)*sizeof(VariablePartData);
	buf->size = size;
	buf->transferToDevice = size > 0;
	buf->transferFromDevice = false;

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
            zeroArray[i].dtGrav = 0.0;
          }
        }
        else{
          buf->hostBuffer = NULL;
        }

	transferKernel.callbackFn = wrCallback;
	transferKernel.id = DM_TRANSFER_LOCAL;
#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
        printf("(%d) DM LOCAL TREE moments %d (%d) partcores %d (%d) partvars %d (%d)\n",
                  CmiMyPe(),
                  transferKernel.bufferInfo[LOCAL_MOMENTS_IDX].size, 
                  transferKernel.bufferInfo[LOCAL_MOMENTS_IDX].transferToDevice,
                  transferKernel.bufferInfo[LOCAL_PARTICLE_CORES_IDX].size, 
                  transferKernel.bufferInfo[LOCAL_PARTICLE_CORES_IDX].transferToDevice,
                  transferKernel.bufferInfo[LOCAL_PARTICLE_VARS_IDX].size, 
                  transferKernel.bufferInfo[LOCAL_PARTICLE_VARS_IDX].transferToDevice
                  );
#endif
#ifdef CUDA_INSTRUMENT_WRS
        transferKernel.chareIndex = mype;
        transferKernel.compType = transferKernel.id;
        transferKernel.compPhase = phase; 
#endif
	enqueue(&transferKernel);

}

#ifdef CUDA_INSTRUMENT_WRS
  void DataManagerTransferRemoteChunk(CudaMultipoleMoments *moments, int nMoments, CompactPartData *compactParts, int nCompactParts, int mype, char phase, void *wrCallback) {
#else
    void DataManagerTransferRemoteChunk(CudaMultipoleMoments *moments, int nMoments, CompactPartData *compactParts, int nCompactParts, void *wrCallback) {
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
  buf->transferFromDevice = false;
  buf->freeBuffer = false;
  size = (nMoments) * sizeof(CudaMultipoleMoments);
  buf->size = size;
  buf->transferToDevice = (size > 0);
  
  //double mill = 1e6;
  //printf("DM remote moments: %f mbytes\n", 1.0*size/mill);
  
  if(size > 0){
    CUDA_MALLOC(buf->hostBuffer, size);
#ifdef CUDA_PRINT_ERRORS
    printf("DMRemote 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
    memcpy(buf->hostBuffer, moments, size);
  }

  buf = &(transferKernel.bufferInfo[REMOTE_PARTICLE_CORES_IDX]);
  buf->bufferID = REMOTE_PARTICLE_CORES;
  buf->transferFromDevice = false;
  buf->freeBuffer = false;
  size = (nCompactParts)*sizeof(CompactPartData);  
  buf->size = size;
  buf->transferToDevice = (size > 0);
  //printf("DM remote cores: %f mbytes\n", 1.0*size/mill);

  if(size > 0){
    CUDA_MALLOC(buf->hostBuffer, size);
#ifdef CUDA_PRINT_ERRORS
    printf("DMRemote 1: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
    memcpy(buf->hostBuffer, compactParts, size);
  }

#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
  printf("(%d) DM REMOTE CHUNK moments %d (%d) partcores %d (%d)\n",
            CmiMyPe(),
            transferKernel.bufferInfo[REMOTE_MOMENTS_IDX].size, 
            transferKernel.bufferInfo[REMOTE_MOMENTS_IDX].transferToDevice,
            transferKernel.bufferInfo[REMOTE_PARTICLE_CORES_IDX].size, 
            transferKernel.bufferInfo[REMOTE_PARTICLE_CORES_IDX].transferToDevice
            );
#endif

  transferKernel.callbackFn = wrCallback;
  transferKernel.id = DM_TRANSFER_REMOTE_CHUNK;
#ifdef CUDA_INSTRUMENT_WRS
  transferKernel.chareIndex = mype;
  transferKernel.compType = transferKernel.id;
  transferKernel.compPhase = phase; 
#endif
  enqueue(&transferKernel);

}


void TreePieceCellListDataTransferLocal(CudaRequest *data){
	int numBlocks = data->numBucketsPlusOne-1;

	workRequest gravityKernel;
	//dataInfo *buffer, *partCoreBuffer;
        //ParameterStruct *pmtr;

	gravityKernel.dimGrid = dim3(numBlocks);
#ifdef CUDA_2D_TB_KERNEL
        gravityKernel.dimBlock = dim3(NODES_PER_BLOCK,PARTS_PER_BLOCK);
#else
	gravityKernel.dimBlock = dim3(THREADS_PER_BLOCK);
#endif
	gravityKernel.smemSize = 0;

	gravityKernel.nBuffers = TP_GRAVITY_LOCAL_NBUFFERS;

	/* schedule buffers for transfer to the GPU */
	gravityKernel.bufferInfo = (dataInfo *) malloc(gravityKernel.nBuffers * sizeof(dataInfo));

	TreePieceCellListDataTransferBasic(data, &gravityKernel);
#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
        printf("(%d) TRANSFER LOCAL CELL\n", CmiMyPe());
#endif

	gravityKernel.callbackFn = data->cb;
	gravityKernel.id = TP_GRAVITY_LOCAL;
#ifdef CUDA_INSTRUMENT_WRS
        gravityKernel.chareIndex = data->tpIndex;
        gravityKernel.compType = gravityKernel.id;
        gravityKernel.compPhase = data->phase; 
#endif
	enqueue(&gravityKernel);
}

void TreePieceCellListDataTransferRemote(CudaRequest *data){
	int numBlocks = data->numBucketsPlusOne-1;

	workRequest gravityKernel;
	//dataInfo *buffer, *partCoreBuffer;

	gravityKernel.dimGrid = dim3(numBlocks);
#ifdef CUDA_2D_TB_KERNEL
        gravityKernel.dimBlock = dim3(NODES_PER_BLOCK,PARTS_PER_BLOCK);
#else
	gravityKernel.dimBlock = dim3(THREADS_PER_BLOCK);
#endif
	gravityKernel.smemSize = 0;

	gravityKernel.nBuffers = TP_NODE_GRAVITY_REMOTE_NBUFFERS;

	/* schedule buffers for transfer to the GPU */
	gravityKernel.bufferInfo = (dataInfo *) malloc(gravityKernel.nBuffers * sizeof(dataInfo));

	TreePieceCellListDataTransferBasic(data, &gravityKernel);
#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
        printf("(%d) TRANSFER REMOTE CELL\n", CmiMyPe());
#endif

	gravityKernel.callbackFn = data->cb;
	gravityKernel.id = TP_GRAVITY_REMOTE;
#ifdef CUDA_INSTRUMENT_WRS
        gravityKernel.chareIndex = data->tpIndex;
        gravityKernel.compType = gravityKernel.id;
        gravityKernel.compPhase = data->phase; 
#endif
	enqueue(&gravityKernel);
}

void TreePieceCellListDataTransferRemoteResume(CudaRequest *data, CudaMultipoleMoments *missedMoments, int numMissedMoments){
  int numBlocks = data->numBucketsPlusOne-1;
  int size;

  workRequest gravityKernel;
  dataInfo *buffer;
  bool transfer;

  gravityKernel.dimGrid = dim3(numBlocks);
#ifdef CUDA_2D_TB_KERNEL
        gravityKernel.dimBlock = dim3(NODES_PER_BLOCK,PARTS_PER_BLOCK);
#else
  gravityKernel.dimBlock = dim3(THREADS_PER_BLOCK);
#endif
  gravityKernel.smemSize = 0;

  gravityKernel.nBuffers = TP_NODE_GRAVITY_REMOTE_RESUME_NBUFFERS;

  /* schedule buffers for transfer to the GPU */
  gravityKernel.bufferInfo = (dataInfo *) malloc(gravityKernel.nBuffers * sizeof(dataInfo));

  TreePieceCellListDataTransferBasic(data, &gravityKernel);
  transfer = gravityKernel.bufferInfo[ILCELL_IDX].transferToDevice;

  buffer = &(gravityKernel.bufferInfo[MISSED_MOMENTS_IDX]);
  buffer->bufferID = -1;
  size = (numMissedMoments) * sizeof(CudaMultipoleMoments);
  buffer->size = size;
  buffer->transferToDevice = transfer;
  buffer->freeBuffer = transfer;
  buffer->transferFromDevice = false;

  if(transfer){
    CUDA_MALLOC(buffer->hostBuffer, size);
#ifdef CUDA_PRINT_ERRORS
    printf("TPRR 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
    memcpy(buffer->hostBuffer, missedMoments, size);
  }
#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
  printf("(%d) TRANSFER REMOTE RESUME CELL %d (%d)\n", CmiMyPe(),
        buffer->size, buffer->transferToDevice);
#endif

  ParameterStruct *ptr = (ParameterStruct *)gravityKernel.userData;
  ptr->numEntities = numMissedMoments;

  gravityKernel.callbackFn = data->cb;
  gravityKernel.id = TP_GRAVITY_REMOTE_RESUME;
#ifdef CUDA_INSTRUMENT_WRS
        gravityKernel.chareIndex = data->tpIndex;
        gravityKernel.compType = gravityKernel.id;
        gravityKernel.compPhase = data->phase; 
#endif
  enqueue(&gravityKernel);
}

void TreePieceCellListDataTransferBasic(CudaRequest *data, workRequest *gravityKernel){
	dataInfo *buffer;
	int numBucketsPlusOne = data->numBucketsPlusOne;
        int numBuckets = numBucketsPlusOne-1;
        int size = (data->numInteractions) * sizeof(ILCell);
        bool transfer = size > 0;

	buffer = &(gravityKernel->bufferInfo[ILCELL_IDX]);
	buffer->bufferID = -1;
	buffer->size = size;
	buffer->transferToDevice = transfer;
	buffer->freeBuffer = transfer;
	buffer->transferFromDevice = false;
        buffer->hostBuffer = data->list;

	buffer = &(gravityKernel->bufferInfo[NODE_BUCKET_MARKERS_IDX]);
	buffer->bufferID = -1;
        size = (numBucketsPlusOne) * sizeof(int);
        buffer->size = transfer ? size : 0;
	buffer->transferToDevice = transfer;
	buffer->freeBuffer = transfer;
	buffer->transferFromDevice = false;
        buffer->hostBuffer = data->bucketMarkers;

	buffer = &(gravityKernel->bufferInfo[NODE_BUCKET_START_MARKERS_IDX]);
	buffer->bufferID = -1;
	size = (numBuckets) * sizeof(int);
        buffer->size = size;
	buffer->transferToDevice = transfer;
	buffer->freeBuffer = transfer;
	buffer->transferFromDevice = false;
        buffer->hostBuffer = data->bucketStarts;

	buffer = &(gravityKernel->bufferInfo[NODE_BUCKET_SIZES_IDX]);
	buffer->bufferID = -1;
        buffer->size = size;
	buffer->transferToDevice = transfer;
	buffer->freeBuffer = transfer;
	buffer->transferFromDevice = false;
        buffer->hostBuffer = data->bucketSizes;

        ParameterStruct *ptr = (ParameterStruct *)malloc(sizeof(ParameterStruct));
        ptr->numInteractions = data->numInteractions;
        ptr->numBucketsPlusOne = numBucketsPlusOne;
        ptr->fperiod = data->fperiod;

        gravityKernel->userData = ptr;

#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
        printf("(%d) TRANSFER BASIC cells %d (%d) bucket_markers %d (%d) bucket_starts %d (%d) bucket_sizes %d (%d)\n",
            CmiMyPe(),
            gravityKernel->bufferInfo[ILCELL_IDX].size, 
            gravityKernel->bufferInfo[ILCELL_IDX].transferToDevice,
            gravityKernel->bufferInfo[NODE_BUCKET_MARKERS_IDX].size, 
            gravityKernel->bufferInfo[NODE_BUCKET_MARKERS_IDX].transferToDevice,
            gravityKernel->bufferInfo[NODE_BUCKET_START_MARKERS_IDX].size, 
            gravityKernel->bufferInfo[NODE_BUCKET_START_MARKERS_IDX].transferToDevice,
            gravityKernel->bufferInfo[NODE_BUCKET_SIZES_IDX].size, 
            gravityKernel->bufferInfo[NODE_BUCKET_SIZES_IDX].transferToDevice
            );
#endif
}

void TreePiecePartListDataTransferLocalSmallPhase(CudaRequest *data, CompactPartData *particles, int len){
	int numBlocks = data->numBucketsPlusOne-1;

	workRequest gravityKernel;
        dataInfo *buffer;
        int size;
        ParameterStruct *ptr;
        bool transfer;

	gravityKernel.dimGrid = dim3(numBlocks);
#ifdef CUDA_2D_TB_KERNEL
	gravityKernel.dimBlock = dim3(NODES_PER_BLOCK_PART,PARTS_PER_BLOCK_PART);
#else
	gravityKernel.dimBlock = dim3(THREADS_PER_BLOCK);
#endif
	gravityKernel.smemSize = 0;

	gravityKernel.nBuffers = TP_GRAVITY_LOCAL_NBUFFERS_SMALLPHASE;

	/* schedule buffers for transfer to the GPU */
	gravityKernel.bufferInfo = (dataInfo *) malloc(gravityKernel.nBuffers * sizeof(dataInfo));

	TreePiecePartListDataTransferBasic(data, &gravityKernel);
        transfer = gravityKernel.bufferInfo[ILPART_IDX].transferToDevice;

	buffer = &(gravityKernel.bufferInfo[MISSED_PARTS_IDX]);
	buffer->bufferID = -1;
        size = (len) * sizeof(CompactPartData);
        buffer->size = size;
	buffer->transferToDevice = transfer;
	buffer->freeBuffer = transfer;
	buffer->transferFromDevice = false;

#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
        printf("(%d) TRANSFER LOCAL SMALL PHASE  %d (%d)\n",
            CmiMyPe(),
            buffer->size, 
            buffer->transferToDevice
            );
#endif

        if(transfer){
          CUDA_MALLOC(buffer->hostBuffer, size);
#ifdef CUDA_PRINT_ERRORS
          printf("TPPartSmallPhase 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
          memcpy(buffer->hostBuffer, particles, size);
        }

        ptr = (ParameterStruct *)gravityKernel.userData;
        ptr->numMissedCores = len;

	gravityKernel.callbackFn = data->cb;
	gravityKernel.id = TP_PART_GRAVITY_LOCAL_SMALLPHASE;
#ifdef CUDA_INSTRUMENT_WRS
        gravityKernel.chareIndex = data->tpIndex;
        gravityKernel.compType = gravityKernel.id;
        gravityKernel.compPhase = data->phase; 
#endif
	enqueue(&gravityKernel);
}

void TreePiecePartListDataTransferLocal(CudaRequest *data){
	int numBlocks = data->numBucketsPlusOne-1;

	workRequest gravityKernel;
	//dataInfo *buffer, *partCoreBuffer;

	gravityKernel.dimGrid = dim3(numBlocks);
#ifdef CUDA_2D_TB_KERNEL
	gravityKernel.dimBlock = dim3(NODES_PER_BLOCK_PART,PARTS_PER_BLOCK_PART);
#else
	gravityKernel.dimBlock = dim3(THREADS_PER_BLOCK);
#endif
	gravityKernel.smemSize = 0;

	gravityKernel.nBuffers = TP_GRAVITY_LOCAL_NBUFFERS;

	/* schedule buffers for transfer to the GPU */
	gravityKernel.bufferInfo = (dataInfo *) malloc(gravityKernel.nBuffers * sizeof(dataInfo));

	TreePiecePartListDataTransferBasic(data, &gravityKernel);

	gravityKernel.callbackFn = data->cb;
	gravityKernel.id = TP_PART_GRAVITY_LOCAL;
#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
        printf("(%d) TRANSFER LOCAL LARGEPHASE PART\n", CmiMyPe());
#endif
#ifdef CUDA_INSTRUMENT_WRS
        gravityKernel.chareIndex = data->tpIndex;
        gravityKernel.compType = gravityKernel.id;
        gravityKernel.compPhase = data->phase; 
#endif
	enqueue(&gravityKernel);
}

void TreePiecePartListDataTransferRemote(CudaRequest *data){
	int numBlocks = data->numBucketsPlusOne-1;

	workRequest gravityKernel;

	gravityKernel.dimGrid = dim3(numBlocks);
#ifdef CUDA_2D_TB_KERNEL
	gravityKernel.dimBlock = dim3(NODES_PER_BLOCK_PART,PARTS_PER_BLOCK_PART);
#else
	gravityKernel.dimBlock = dim3(THREADS_PER_BLOCK);
#endif
	gravityKernel.smemSize = 0;

	gravityKernel.nBuffers = TP_PART_GRAVITY_REMOTE_NBUFFERS;

	gravityKernel.bufferInfo = (dataInfo *) malloc(gravityKernel.nBuffers * sizeof(dataInfo));

	TreePiecePartListDataTransferBasic(data, &gravityKernel);

	gravityKernel.callbackFn = data->cb;
	gravityKernel.id = TP_PART_GRAVITY_REMOTE;
#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
        printf("(%d) TRANSFER REMOTE PART\n", CmiMyPe());
#endif
#ifdef CUDA_INSTRUMENT_WRS
        gravityKernel.chareIndex = data->tpIndex;
        gravityKernel.compType = gravityKernel.id;
        gravityKernel.compPhase = data->phase; 
#endif
	enqueue(&gravityKernel);
}

void TreePiecePartListDataTransferRemoteResume(CudaRequest *data, CompactPartData *missedParts, int numMissedParts){
	int numBlocks = data->numBucketsPlusOne-1;
        int size;
        ParameterStruct *ptr;
        bool transfer;

	workRequest gravityKernel;
        dataInfo *buffer;

	gravityKernel.dimGrid = dim3(numBlocks);
#ifdef CUDA_2D_TB_KERNEL
	gravityKernel.dimBlock = dim3(NODES_PER_BLOCK_PART,PARTS_PER_BLOCK_PART);
#else
	gravityKernel.dimBlock = dim3(THREADS_PER_BLOCK);
#endif
	gravityKernel.smemSize = 0;

	gravityKernel.nBuffers = TP_PART_GRAVITY_REMOTE_RESUME_NBUFFERS;

	/* schedule buffers for transfer to the GPU */
	gravityKernel.bufferInfo = (dataInfo *) malloc(gravityKernel.nBuffers * sizeof(dataInfo));

	TreePiecePartListDataTransferBasic(data, &gravityKernel);
        transfer = gravityKernel.bufferInfo[ILPART_IDX].transferToDevice;

	buffer = &(gravityKernel.bufferInfo[MISSED_PARTS_IDX]);
	buffer->bufferID = -1;
        size = (numMissedParts) * sizeof(CompactPartData);
        buffer->size = size;
	buffer->transferToDevice = transfer;
	buffer->freeBuffer = transfer;
	buffer->transferFromDevice = false;

#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
        printf("(%d) TRANSFER REMOTE RESUME PART %d (%d)\n",
            CmiMyPe(),
            buffer->size, 
            buffer->transferToDevice
            );
#endif

        if(transfer){
          CUDA_MALLOC(buffer->hostBuffer, size);
#ifdef CUDA_PRINT_ERRORS
          printf("TPPartRR 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
          memcpy(buffer->hostBuffer, missedParts, size);
        }

        ptr = (ParameterStruct *)gravityKernel.userData;
        ptr->numMissedCores = numMissedParts;

	gravityKernel.callbackFn = data->cb;
	gravityKernel.id = TP_PART_GRAVITY_REMOTE_RESUME;
#ifdef CUDA_INSTRUMENT_WRS
        gravityKernel.chareIndex = data->tpIndex;
        gravityKernel.compType = gravityKernel.id;
        gravityKernel.compPhase = data->phase; 
#endif
	enqueue(&gravityKernel);
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
	size = (numInteractions) * sizeof(ILCell);
        buffer->size = size;
        transfer = size > 0;
	buffer->transferToDevice = transfer;
	buffer->freeBuffer = transfer;
	buffer->transferFromDevice = false;
        buffer->hostBuffer = data->list;

	buffer = &(gravityKernel->bufferInfo[PART_BUCKET_MARKERS_IDX]);
	buffer->bufferID = -1;
	size = (numBucketsPlusOne) * sizeof(int);
        buffer->size = transfer ? size : 0;
	buffer->transferToDevice = transfer;
	buffer->freeBuffer = transfer;
	buffer->transferFromDevice = false;
        buffer->hostBuffer = data->bucketMarkers;

	buffer = &(gravityKernel->bufferInfo[PART_BUCKET_START_MARKERS_IDX]);
	buffer->bufferID = -1;
	buffer->hostBuffer = data->bucketStarts;
	size = (numBuckets) * sizeof(int);
        buffer->size = size;
	buffer->transferToDevice = transfer;
	buffer->freeBuffer = transfer;
	buffer->transferFromDevice = false;
        buffer->hostBuffer = data->bucketStarts;

	buffer = &(gravityKernel->bufferInfo[PART_BUCKET_SIZES_IDX]);
	buffer->bufferID = -1;
        buffer->size = size;
	buffer->transferToDevice = transfer;
	buffer->freeBuffer = transfer;
	buffer->transferFromDevice = false;
        buffer->hostBuffer = data->bucketSizes;
        
        ParameterStruct *ptr = (ParameterStruct *)malloc(sizeof(ParameterStruct));
        ptr->numInteractions = data->numInteractions;
        ptr->numBucketsPlusOne = numBucketsPlusOne;
        ptr->fperiod = data->fperiod;
                                       
        gravityKernel->userData = ptr;

#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
        printf("(%d) TRANSFER BASIC PART parts %d (%d) bucket_markers %d (%d) bucket_starts %d (%d) bucket_sizes %d (%d)\n",
            CmiMyPe(),
            gravityKernel->bufferInfo[ILPART_IDX].size, 
            gravityKernel->bufferInfo[ILPART_IDX].transferToDevice,
            gravityKernel->bufferInfo[PART_BUCKET_MARKERS_IDX].size, 
            gravityKernel->bufferInfo[PART_BUCKET_MARKERS_IDX].transferToDevice,
            gravityKernel->bufferInfo[PART_BUCKET_START_MARKERS_IDX].size, 
            gravityKernel->bufferInfo[PART_BUCKET_START_MARKERS_IDX].transferToDevice,
            gravityKernel->bufferInfo[PART_BUCKET_SIZES_IDX].size, 
            gravityKernel->bufferInfo[PART_BUCKET_SIZES_IDX].transferToDevice
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
  buffer->transferToDevice = false ;
  buffer->transferFromDevice = false;
  buffer->freeBuffer = freemom;
  buffer->hostBuffer = 0;
  buffer->size = 0;

  buffer = &(gravityKernel.bufferInfo[LOCAL_PARTICLE_CORES_IDX]);
  buffer->bufferID = LOCAL_PARTICLE_CORES;
  buffer->transferToDevice = false ;
  buffer->transferFromDevice = false;
  buffer->freeBuffer = freepart;
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
  enqueue(&gravityKernel);

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
  buffer->transferToDevice = false;
  buffer->transferFromDevice = false;
  buffer->freeBuffer = freemom;
  buffer->size = 0;

  buffer = &(gravityKernel.bufferInfo[REMOTE_PARTICLE_CORES_IDX]);
  buffer->bufferID = REMOTE_PARTICLE_CORES;
  buffer->transferToDevice = false;
  buffer->transferFromDevice = false;
  buffer->freeBuffer = freepart;
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
  enqueue(&gravityKernel);

}

#ifdef CUDA_INSTRUMENT_WRS
void TransferParticleVarsBack(VariablePartData *hostBuffer, int size, void *cb, bool freemom, bool freepart, int index, char phase){
#else
void TransferParticleVarsBack(VariablePartData *hostBuffer, int size, void *cb, bool freemom, bool freepart){
#endif
  workRequest gravityKernel;
  dataInfo *buffer;

#ifdef CUDA_PRINT_TRANSFER_BACK_PARTICLES
  printf("Enqueue kernel to transfer particles back from: 0x%x\n", devBuffers[LOCAL_PARTICLE_VARS]);
#endif
  gravityKernel.dimGrid = dim3(0);
  gravityKernel.dimBlock = dim3(0);
  gravityKernel.smemSize = 0;

  gravityKernel.nBuffers = 3;

  /* schedule buffers for transfer to the GPU */
  gravityKernel.bufferInfo = (dataInfo *) malloc(3 * sizeof(dataInfo));

  buffer = &(gravityKernel.bufferInfo[LOCAL_PARTICLE_VARS_IDX]);
  buffer->bufferID = LOCAL_PARTICLE_VARS;
  buffer->transferToDevice = false;
  buffer->hostBuffer = hostBuffer;
  buffer->size = size;
  buffer->transferFromDevice = size > 0;
  buffer->freeBuffer = size > 0;

  buffer = &(gravityKernel.bufferInfo[LOCAL_MOMENTS_IDX]);
  buffer->bufferID = LOCAL_MOMENTS;
  buffer->transferToDevice = false;
  buffer->transferFromDevice = false;
  buffer->freeBuffer = freemom;
  buffer->hostBuffer = NULL;
  buffer->size = 0;

  buffer = &(gravityKernel.bufferInfo[LOCAL_PARTICLE_CORES_IDX]);
  buffer->bufferID = LOCAL_PARTICLE_CORES;
  buffer->transferToDevice = false ;
  buffer->transferFromDevice = false;
  buffer->freeBuffer = freepart;
  buffer->hostBuffer = NULL;
  buffer->size = 0;


  gravityKernel.callbackFn = cb;
  gravityKernel.id = DM_TRANSFER_BACK;
#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
  printf("(%d) DM TRANSFER BACK\n", CmiMyPe());
#endif
#ifdef CUDA_INSTRUMENT_WRS
  gravityKernel.chareIndex = index;
  gravityKernel.compType = gravityKernel.id;
  gravityKernel.compPhase = phase; 
#endif
  enqueue(&gravityKernel);
}

/*
void DummyKernel(void *cb){
  workRequest dummy;
  //dataInfo *buffer;

  dummy.dimGrid = dim3(1);
  dummy.dimBlock = dim3(THREADS_PER_BLOCK);
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
  cudaStream_t kernel_stream = getKernelStream(); 
  void** devBuffers = getdevBuffers(); 
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
      if( wr->bufferInfo[LOCAL_MOMENTS_IDX].transferToDevice ){
#ifdef CUDA_MEMPOOL
        hapi_poolFree(wr->bufferInfo[LOCAL_MOMENTS_IDX].hostBuffer);
#else
        delayedFree(wr->bufferInfo[LOCAL_MOMENTS_IDX].hostBuffer);
#endif
      }
#ifdef CUDA_PRINT_ERRORS
      printf("DM_TRANSFER_LOCAL 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
      if( wr->bufferInfo[LOCAL_PARTICLE_CORES_IDX].transferToDevice ){
#ifdef CUDA_MEMPOOL
        hapi_poolFree(wr->bufferInfo[LOCAL_PARTICLE_CORES_IDX].hostBuffer);
#else
        delayedFree(wr->bufferInfo[LOCAL_PARTICLE_CORES_IDX].hostBuffer);
#endif
      }
#ifdef CUDA_PRINT_ERRORS
      printf("DM_TRANSFER_LOCAL 1: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
      if( wr->bufferInfo[LOCAL_PARTICLE_VARS_IDX].transferToDevice ){
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
      if( wr->bufferInfo[LOCAL_MOMENTS_IDX].transferToDevice ){
        free(wr->bufferInfo[LOCAL_MOMENTS_IDX].hostBuffer);
      }
      if( wr->bufferInfo[LOCAL_PARTICLE_CORES_IDX].transferToDevice ){
        free(wr->bufferInfo[LOCAL_PARTICLE_CORES_IDX].hostBuffer);
      }
      if( wr->bufferInfo[LOCAL_PARTICLE_VARS_IDX].transferToDevice ){
        free(wr->bufferInfo[LOCAL_PARTICLE_VARS_IDX].hostBuffer);
      }
#endif
      break;

    case DM_TRANSFER_REMOTE_CHUNK:
#ifdef CUDA_NOTIFY_DATA_TRANSFER_DONE
      printf("DM_TRANSFER_REMOTE_CHUNK, %d KERNELSELECT\n", wr->bufferInfo[REMOTE_MOMENTS_IDX].transferToDevice);
#endif
      if( wr->bufferInfo[REMOTE_MOMENTS_IDX].transferToDevice ){
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

      if( wr->bufferInfo[REMOTE_PARTICLE_CORES_IDX].transferToDevice ){
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
      if( wr->bufferInfo[ILCELL_IDX].transferToDevice ){
#ifndef CUDA_NO_KERNELS
        nodeGravityComputation<<<wr->dimGrid, wr->dimBlock, wr->smemSize, kernel_stream>>>
          (
           (CompactPartData *)devBuffers[LOCAL_PARTICLE_CORES],
           (VariablePartData *)devBuffers[LOCAL_PARTICLE_VARS],
           (CudaMultipoleMoments *)devBuffers[LOCAL_MOMENTS],
           (ILCell *)devBuffers[wr->bufferInfo[ILCELL_IDX].bufferID],
           (int *)devBuffers[wr->bufferInfo[NODE_BUCKET_MARKERS_IDX].bufferID],
           (int *)devBuffers[wr->bufferInfo[NODE_BUCKET_START_MARKERS_IDX].bufferID],
           (int *)devBuffers[wr->bufferInfo[NODE_BUCKET_SIZES_IDX].bufferID],
           ptr->fperiod
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
      if( wr->bufferInfo[ILPART_IDX].transferToDevice ){
#ifndef CUDA_NO_KERNELS
#ifdef CUDA_2D_TB_KERNEL
        particleGravityComputation<<<wr->dimGrid, wr->dimBlock, wr->smemSize, kernel_stream>>>
          (
           (CompactPartData *)devBuffers[LOCAL_PARTICLE_CORES],
           (VariablePartData *)devBuffers[LOCAL_PARTICLE_VARS],
           (CompactPartData *)devBuffers[wr->bufferInfo[MISSED_PARTS_IDX].bufferID],
           (ILCell *)devBuffers[wr->bufferInfo[ILPART_IDX].bufferID],
           (int *)devBuffers[wr->bufferInfo[PART_BUCKET_MARKERS_IDX].bufferID],
           (int *)devBuffers[wr->bufferInfo[PART_BUCKET_START_MARKERS_IDX].bufferID],
           (int *)devBuffers[wr->bufferInfo[PART_BUCKET_SIZES_IDX].bufferID],
           ptr->fperiod 
          );
#else
        particleGravityComputation<<<wr->dimGrid, wr->dimBlock, wr->smemSize, kernel_stream>>>
          (
           (CompactPartData *)devBuffers[LOCAL_PARTICLE_CORES],
           (VariablePartData *)devBuffers[LOCAL_PARTICLE_VARS],
           (CompactPartData *)devBuffers[wr->bufferInfo[MISSED_PARTS_IDX].bufferID],
           (ILPart *)devBuffers[wr->bufferInfo[ILPART_IDX].bufferID],
           (int *)devBuffers[wr->bufferInfo[PART_BUCKET_MARKERS_IDX].bufferID],
           (int *)devBuffers[wr->bufferInfo[PART_BUCKET_START_MARKERS_IDX].bufferID],
           (int *)devBuffers[wr->bufferInfo[PART_BUCKET_SIZES_IDX].bufferID],
           ptr->fperiod 
          );
#endif
#endif

#ifdef CUDA_TRACE
    st = CmiWallTimer();
#endif
#ifdef CUDA_USE_CUDAMALLOCHOST
#ifdef CUDA_MEMPOOL
        hapi_poolFree((CompactPartData *)wr->bufferInfo[MISSED_PARTS_IDX].hostBuffer);
        hapi_poolFree(wr->bufferInfo[ILPART_IDX].hostBuffer);
        hapi_poolFree((int *)wr->bufferInfo[PART_BUCKET_MARKERS_IDX].hostBuffer);
        hapi_poolFree((int *)wr->bufferInfo[PART_BUCKET_START_MARKERS_IDX].hostBuffer);
        hapi_poolFree((int *)wr->bufferInfo[PART_BUCKET_SIZES_IDX].hostBuffer);
#else
        delayedFree((CompactPartData *)wr->bufferInfo[MISSED_PARTS_IDX].hostBuffer);
        delayedFree(wr->bufferInfo[ILPART_IDX].hostBuffer);
        delayedFree((int *)wr->bufferInfo[PART_BUCKET_MARKERS_IDX].hostBuffer);
        delayedFree((int *)wr->bufferInfo[PART_BUCKET_START_MARKERS_IDX].hostBuffer);
        delayedFree((int *)wr->bufferInfo[PART_BUCKET_SIZES_IDX].hostBuffer);
#endif
#ifdef CUDA_PRINT_ERRORS
        printf("TP_PART_GRAVITY_LOCAL_SMALLPHASE 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
#else
        free((CompactPartData *)wr->bufferInfo[MISSED_PARTS_IDX].hostBuffer);
        free(wr->bufferInfo[ILPART_IDX].hostBuffer);
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
      if( wr->bufferInfo[ILPART_IDX].transferToDevice ){
#ifndef CUDA_NO_KERNELS
#ifdef CUDA_2D_TB_KERNEL
        particleGravityComputation<<<wr->dimGrid, wr->dimBlock, wr->smemSize, kernel_stream>>>
          (
           (CompactPartData *)devBuffers[LOCAL_PARTICLE_CORES],
           (VariablePartData *)devBuffers[LOCAL_PARTICLE_VARS],
           (CompactPartData *)devBuffers[LOCAL_PARTICLE_CORES],
           (ILCell *)devBuffers[wr->bufferInfo[ILPART_IDX].bufferID],
           (int *)devBuffers[wr->bufferInfo[PART_BUCKET_MARKERS_IDX].bufferID],
           (int *)devBuffers[wr->bufferInfo[PART_BUCKET_START_MARKERS_IDX].bufferID],
           (int *)devBuffers[wr->bufferInfo[PART_BUCKET_SIZES_IDX].bufferID],
           ptr->fperiod
          );
#else
        particleGravityComputation<<<wr->dimGrid, wr->dimBlock, wr->smemSize, kernel_stream>>>
          (
           (CompactPartData *)devBuffers[LOCAL_PARTICLE_CORES],
           (VariablePartData *)devBuffers[LOCAL_PARTICLE_VARS],
           (CompactPartData *)devBuffers[LOCAL_PARTICLE_CORES],
           (ILPart *)devBuffers[wr->bufferInfo[ILPART_IDX].bufferID],
           (int *)devBuffers[wr->bufferInfo[PART_BUCKET_MARKERS_IDX].bufferID],
           (int *)devBuffers[wr->bufferInfo[PART_BUCKET_START_MARKERS_IDX].bufferID],
           (int *)devBuffers[wr->bufferInfo[PART_BUCKET_SIZES_IDX].bufferID],
           ptr->fperiod
          );
#endif
#endif

#ifdef CUDA_TRACE
    st = CmiWallTimer();
#endif
#ifdef CUDA_USE_CUDAMALLOCHOST
#ifdef CUDA_MEMPOOL
        hapi_poolFree(wr->bufferInfo[ILPART_IDX].hostBuffer);
        hapi_poolFree((int *)wr->bufferInfo[PART_BUCKET_MARKERS_IDX].hostBuffer);
        hapi_poolFree((int *)wr->bufferInfo[PART_BUCKET_START_MARKERS_IDX].hostBuffer);
        hapi_poolFree((int *)wr->bufferInfo[PART_BUCKET_SIZES_IDX].hostBuffer);
#else
        delayedFree(wr->bufferInfo[ILPART_IDX].hostBuffer);
        delayedFree((int *)wr->bufferInfo[PART_BUCKET_MARKERS_IDX].hostBuffer);
        delayedFree((int *)wr->bufferInfo[PART_BUCKET_START_MARKERS_IDX].hostBuffer);
        delayedFree((int *)wr->bufferInfo[PART_BUCKET_SIZES_IDX].hostBuffer);
#endif
#ifdef CUDA_PRINT_ERRORS
        printf("TP_PART_GRAVITY_LOCAL 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
#else
        free(wr->bufferInfo[ILPART_IDX].hostBuffer);
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
      if( wr->bufferInfo[ILCELL_IDX].transferToDevice ){
#ifndef CUDA_NO_KERNELS
        nodeGravityComputation<<<wr->dimGrid, wr->dimBlock, wr->smemSize, kernel_stream>>>
          (
           (CompactPartData *)devBuffers[LOCAL_PARTICLE_CORES],
           (VariablePartData *)devBuffers[LOCAL_PARTICLE_VARS],
           (CudaMultipoleMoments *)devBuffers[REMOTE_MOMENTS],
           (ILCell *)devBuffers[wr->bufferInfo[ILCELL_IDX].bufferID],
           (int *)devBuffers[wr->bufferInfo[NODE_BUCKET_MARKERS_IDX].bufferID],
           (int *)devBuffers[wr->bufferInfo[NODE_BUCKET_START_MARKERS_IDX].bufferID],
           (int *)devBuffers[wr->bufferInfo[NODE_BUCKET_SIZES_IDX].bufferID],
           ptr->fperiod
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
      if( wr->bufferInfo[ILPART_IDX].transferToDevice ){
#ifndef CUDA_NO_KERNELS
#ifdef CUDA_2D_TB_KERNEL
        particleGravityComputation<<<wr->dimGrid, wr->dimBlock, wr->smemSize, kernel_stream>>>
          (
           (CompactPartData *)devBuffers[LOCAL_PARTICLE_CORES],
           (VariablePartData *)devBuffers[LOCAL_PARTICLE_VARS],
           (CompactPartData *)devBuffers[REMOTE_PARTICLE_CORES],
           (ILCell *)devBuffers[wr->bufferInfo[ILPART_IDX].bufferID],
           (int *)devBuffers[wr->bufferInfo[PART_BUCKET_MARKERS_IDX].bufferID],
           (int *)devBuffers[wr->bufferInfo[PART_BUCKET_START_MARKERS_IDX].bufferID],
           (int *)devBuffers[wr->bufferInfo[PART_BUCKET_SIZES_IDX].bufferID],
           ptr->fperiod
          );
#else
        particleGravityComputation<<<wr->dimGrid, wr->dimBlock, wr->smemSize, kernel_stream>>>
          (
           (CompactPartData *)devBuffers[LOCAL_PARTICLE_CORES],
           (VariablePartData *)devBuffers[LOCAL_PARTICLE_VARS],
           (CompactPartData *)devBuffers[REMOTE_PARTICLE_CORES],
           (ILPart *)devBuffers[wr->bufferInfo[ILPART_IDX].bufferID],
           (int *)devBuffers[wr->bufferInfo[PART_BUCKET_MARKERS_IDX].bufferID],
           (int *)devBuffers[wr->bufferInfo[PART_BUCKET_START_MARKERS_IDX].bufferID],
           (int *)devBuffers[wr->bufferInfo[PART_BUCKET_SIZES_IDX].bufferID],
           ptr->fperiod
          );
#endif
#endif

#ifdef CUDA_TRACE
    st = CmiWallTimer();
#endif
#ifdef CUDA_USE_CUDAMALLOCHOST
#ifdef CUDA_MEMPOOL
        hapi_poolFree(wr->bufferInfo[ILPART_IDX].hostBuffer);
        hapi_poolFree((int *)wr->bufferInfo[PART_BUCKET_MARKERS_IDX].hostBuffer);
        hapi_poolFree((int *)wr->bufferInfo[PART_BUCKET_START_MARKERS_IDX].hostBuffer);
        hapi_poolFree((int *)wr->bufferInfo[PART_BUCKET_SIZES_IDX].hostBuffer);
#else
        delayedFree(wr->bufferInfo[ILPART_IDX].hostBuffer);
        delayedFree((int *)wr->bufferInfo[PART_BUCKET_MARKERS_IDX].hostBuffer);
        delayedFree((int *)wr->bufferInfo[PART_BUCKET_START_MARKERS_IDX].hostBuffer);
        delayedFree((int *)wr->bufferInfo[PART_BUCKET_SIZES_IDX].hostBuffer);
#endif
#ifdef CUDA_PRINT_ERRORS
        printf("TP_PART_GRAVITY_REMOTE 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
#else
        free(wr->bufferInfo[ILPART_IDX].hostBuffer);
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
      if( wr->bufferInfo[ILCELL_IDX].transferToDevice ){
#ifndef CUDA_NO_KERNELS
        nodeGravityComputation<<<wr->dimGrid, wr->dimBlock, wr->smemSize, kernel_stream>>>
          (
           (CompactPartData *)devBuffers[LOCAL_PARTICLE_CORES],
           (VariablePartData *)devBuffers[LOCAL_PARTICLE_VARS],
           (CudaMultipoleMoments *)devBuffers[wr->bufferInfo[MISSED_MOMENTS_IDX].bufferID],
           (ILCell *)devBuffers[wr->bufferInfo[ILCELL_IDX].bufferID],
           (int *)devBuffers[wr->bufferInfo[NODE_BUCKET_MARKERS_IDX].bufferID],
           (int *)devBuffers[wr->bufferInfo[NODE_BUCKET_START_MARKERS_IDX].bufferID],
           (int *)devBuffers[wr->bufferInfo[NODE_BUCKET_SIZES_IDX].bufferID],
           ptr->fperiod
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
      if( wr->bufferInfo[ILPART_IDX].transferToDevice ){
#ifndef CUDA_NO_KERNELS
#ifdef CUDA_2D_TB_KERNEL
        particleGravityComputation<<<wr->dimGrid, wr->dimBlock, wr->smemSize, kernel_stream>>>
          (
           (CompactPartData *)devBuffers[LOCAL_PARTICLE_CORES],
           (VariablePartData *)devBuffers[LOCAL_PARTICLE_VARS],
           (CompactPartData *)devBuffers[wr->bufferInfo[MISSED_PARTS_IDX].bufferID],
           (ILCell *)devBuffers[wr->bufferInfo[ILPART_IDX].bufferID],
           (int *)devBuffers[wr->bufferInfo[PART_BUCKET_MARKERS_IDX].bufferID],
           (int *)devBuffers[wr->bufferInfo[PART_BUCKET_START_MARKERS_IDX].bufferID],
           (int *)devBuffers[wr->bufferInfo[PART_BUCKET_SIZES_IDX].bufferID],
           ptr->fperiod
          );
#else
        particleGravityComputation<<<wr->dimGrid, wr->dimBlock, wr->smemSize, kernel_stream>>>
          (
           (CompactPartData *)devBuffers[LOCAL_PARTICLE_CORES],
           (VariablePartData *)devBuffers[LOCAL_PARTICLE_VARS],
           (CompactPartData *)devBuffers[wr->bufferInfo[MISSED_PARTS_IDX].bufferID],
           (ILPart *)devBuffers[wr->bufferInfo[ILPART_IDX].bufferID],
           (int *)devBuffers[wr->bufferInfo[PART_BUCKET_MARKERS_IDX].bufferID],
           (int *)devBuffers[wr->bufferInfo[PART_BUCKET_START_MARKERS_IDX].bufferID],
           (int *)devBuffers[wr->bufferInfo[PART_BUCKET_SIZES_IDX].bufferID],
           ptr->fperiod
          );
#endif
#endif

#ifdef CUDA_TRACE
    st = CmiWallTimer();
#endif
#ifdef CUDA_USE_CUDAMALLOCHOST
#ifdef CUDA_MEMPOOL
        hapi_poolFree((CompactPartData *)wr->bufferInfo[MISSED_PARTS_IDX].hostBuffer);
        hapi_poolFree(wr->bufferInfo[ILPART_IDX].hostBuffer);
        hapi_poolFree((int *)wr->bufferInfo[PART_BUCKET_MARKERS_IDX].hostBuffer);
        hapi_poolFree((int *)wr->bufferInfo[PART_BUCKET_START_MARKERS_IDX].hostBuffer);
        hapi_poolFree((int *)wr->bufferInfo[PART_BUCKET_SIZES_IDX].hostBuffer);
#else
        delayedFree((CompactPartData *)wr->bufferInfo[MISSED_PARTS_IDX].hostBuffer);
        delayedFree(wr->bufferInfo[ILPART_IDX].hostBuffer);
        delayedFree((int *)wr->bufferInfo[PART_BUCKET_MARKERS_IDX].hostBuffer);
        delayedFree((int *)wr->bufferInfo[PART_BUCKET_START_MARKERS_IDX].hostBuffer);
        delayedFree((int *)wr->bufferInfo[PART_BUCKET_SIZES_IDX].hostBuffer);
#endif
#ifdef CUDA_PRINT_ERRORS
        printf("TP_PART_GRAVITY_REMOTE_RESUME 0: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
#else
        free((CompactPartData *)wr->bufferInfo[MISSED_PARTS_IDX].hostBuffer);
        free(wr->bufferInfo[ILPART_IDX].hostBuffer);
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
        ((GravityParticleData*)devBuffers[wr->bufferInfo[PARTICLE_TABLE].bufferID],
         wr->bufferInfo[PARTICLE_TABLE].size/sizeof(GravityParticleData));
#ifdef CUDA_PRINT_ERRORS
      printf("TOP_EWALD_KERNEL: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
      

      break; 
    case BOTTOM_EWALD_KERNEL:
      cudaMemcpyToSymbol(ewt, wr->bufferInfo[EWALD_TABLE].hostBuffer, NEWH * sizeof(EwtData), 0, cudaMemcpyHostToDevice);
      EwaldBottomKernel<<<wr->dimGrid, wr->dimBlock, 
        wr->smemSize, kernel_stream>>>
          ((GravityParticleData *)devBuffers[wr->bufferInfo[PARTICLE_TABLE].bufferID],
            wr->bufferInfo[PARTICLE_TABLE].size/sizeof(GravityParticleData));
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
__global__ void nodeGravityComputation(
		CompactPartData *particleCores,
		VariablePartData *particleVars,
		CudaMultipoleMoments *moments,
		ILCell *ils,
		int *ilmarks,
		int *bucketStarts,
		int *bucketSizes,
		cudatype fperiod){
  
  __shared__ CudaVector3D acc[THREADS_PER_BLOCK];
  __shared__ cudatype pot[THREADS_PER_BLOCK];
  __shared__ cudatype idt2[THREADS_PER_BLOCK];
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
    
    __syncthreads(); // wait for leader threads to finish using acc's, pot's of other threads
    acc[TRANSLATE(tidx,tidy)].x = 0.0;
    acc[TRANSLATE(tidx,tidy)].y = 0.0;
    acc[TRANSLATE(tidx,tidy)].z = 0.0;
    pot[TRANSLATE(tidx,tidy)] = 0.0;
    idt2[TRANSLATE(tidx,tidy)] = 0.0;
    
    
    for(int xstart = start; xstart < end; xstart += NODES_PER_BLOCK){
      int my_cell_idx = xstart + tidx;
      ILCell ilc;

      __syncthreads(); // wait for all threads to finish using 
                       // previous iteration's nodes before reloading
      
      if(tidy == 0 && my_cell_idx < end){
        ilc = ils[my_cell_idx];
        m[tidx] = moments[ilc.index];
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
          CUDA_momEvalFmomrcm(&m[tidx], &r, dir, &acc[TRANSLATE(tidx, tidy)], &pot[TRANSLATE(tidx, tidy)]);
          idt2[TRANSLATE(tidx, tidy)] = fmax(idt2[TRANSLATE(tidx, tidy)],
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

          pot[TRANSLATE(tidx, tidy)] -= m[tidx].totalMass * a + c*qir - b*tr;

          acc[TRANSLATE(tidx, tidy)].x -= qir3*r.x - c*qirx;
          acc[TRANSLATE(tidx, tidy)].y -= qir3*r.y - c*qiry;
          acc[TRANSLATE(tidx, tidy)].z -= qir3*r.z - c*qirz;
          idt2[TRANSLATE(tidx, tidy)] = fmax(idt2[TRANSLATE(tidx, tidy)],
                                        (shared_particle_cores[tidy].mass + m[tidx].totalMass)*b);
#endif
        }// end if rsq != 0
      }// end INTERACT
    }// end for each NODE group

    __syncthreads(); // wait for all threads to finish before results become available

    cudatype sumx, sumy, sumz, poten, idt2max;
    sumx = sumy = sumz = poten = idt2max = 0.0;
    // accumulate forces, potential in global memory data structure
    if(tidx == 0 && my_particle_idx < bucketSize){
      for(int i = 0; i < NODES_PER_BLOCK; i++){
        sumx += acc[TRANSLATE(i,tidy)].x;
        sumy += acc[TRANSLATE(i,tidy)].y;
        sumz += acc[TRANSLATE(i,tidy)].z;
        poten += pot[TRANSLATE(i,tidy)];
        idt2max = fmax(idt2[TRANSLATE(i,tidy)], idt2max);
      }
      particleVars[bucketStart+my_particle_idx].a.x += sumx;
      particleVars[bucketStart+my_particle_idx].a.y += sumy;
      particleVars[bucketStart+my_particle_idx].a.z += sumz;
      particleVars[bucketStart+my_particle_idx].potential += poten;
      particleVars[bucketStart+my_particle_idx].dtGrav = fmax(idt2max,  particleVars[bucketStart+my_particle_idx].dtGrav);
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
      particleVars[bucketStart+particle].dtGrav = fmax(idt2max,  particleVars[bucketStart+particle].dtGrav);
      
#if defined  __DEVICE_EMULATION__ && defined CUDA_EMU_KERNEL_NODE_PRINTS
      printf("after (%d,%d) type: NODE %d acc: (%f,%f,%f)\n", tp, id, type,
                                                    particleVars[bucketStart+particle].a.x,
                                                    particleVars[bucketStart+particle].a.y,
                                                    particleVars[bucketStart+particle].a.z
                                                    );
#endif
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
 
#ifdef CUDA_2D_TB_KERNEL
#define TRANSLATE_PART(x,y) (y*NODES_PER_BLOCK_PART+x)
__global__ void particleGravityComputation(
		CompactPartData *targetCores,
		VariablePartData *targetVars,
		CompactPartData *sourceCores,
		ILCell *ils,
		int *ilmarks,
		int *bucketStarts,
		int *bucketSizes,
		cudatype fperiod){
  
  __shared__ CudaVector3D acc[THREADS_PER_BLOCK_PART];
  __shared__ cudatype pot[THREADS_PER_BLOCK_PART];
  __shared__ cudatype idt2[THREADS_PER_BLOCK_PART];
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
    
    __syncthreads(); // wait for leader threads to finish using acc's, pot's of other threads
    acc[TRANSLATE_PART(tx,ty)].x = 0.0;
    acc[TRANSLATE_PART(tx,ty)].y = 0.0;
    acc[TRANSLATE_PART(tx,ty)].z = 0.0;
    pot[TRANSLATE_PART(tx,ty)] = 0.0;
    idt2[TRANSLATE_PART(tx,ty)] = 0.0;
    
    
    for(xstart = start; xstart < end; xstart += NODES_PER_BLOCK_PART){
      int my_cell_idx = xstart + tx;
      ILCell ilc;

      __syncthreads(); // wait for all threads to finish using 
                       // previous iteration's nodes before reloading
      
      if(ty == 0 && my_cell_idx < end){
        ilc = ils[my_cell_idx];
        m[tx] = sourceCores[ilc.index];
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

          pot[TRANSLATE_PART(tx, ty)] -= m[tx].mass * a;

          acc[TRANSLATE_PART(tx, ty)].x += r.x*b*m[tx].mass;
          acc[TRANSLATE_PART(tx, ty)].y += r.y*b*m[tx].mass;
          acc[TRANSLATE_PART(tx, ty)].z += r.z*b*m[tx].mass;
          idt2[TRANSLATE_PART(tx, ty)] = fmax(idt2[TRANSLATE_PART(tx, ty)],
                                         (shared_particle_cores[ty].mass + m[tx].mass) * b);

        }// end if rsq != 0
      }// end INTERACT
    }// end for each NODE group

    __syncthreads(); // wait for all threads to finish before results become available

    cudatype sumx, sumy, sumz, poten, idt2max;
    sumx = sumy = sumz = poten = idt2max = 0.0;
    // accumulate forces, potential in global memory data structure
    if(tx == 0 && my_particle_idx < bucketSize){
      for(int i = 0; i < NODES_PER_BLOCK_PART; i++){
        sumx += acc[TRANSLATE_PART(i,ty)].x;
        sumy += acc[TRANSLATE_PART(i,ty)].y;
        sumz += acc[TRANSLATE_PART(i,ty)].z;
        poten += pot[TRANSLATE_PART(i,ty)];
        idt2max = fmax(idt2[TRANSLATE_PART(i,ty)], idt2max);
      }
      targetVars[bucketStart+my_particle_idx].a.x += sumx;
      targetVars[bucketStart+my_particle_idx].a.y += sumy;
      targetVars[bucketStart+my_particle_idx].a.z += sumz;
      targetVars[bucketStart+my_particle_idx].potential += poten;
      targetVars[bucketStart+my_particle_idx].dtGrav = fmax(idt2max,  targetVars[bucketStart+my_particle_idx].dtGrav);
    }

  }// end for each PARTICLE group
}
#else
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
      targetVars[bucketStart+target].dtGrav = fmax(idt2max,  targetVars[bucketStart+target].dtGrav);
      
#if defined  __DEVICE_EMULATION__ && defined CUDA_EMU_KERNEL_PART_PRINTS
      printf("after (%d,%d) type: %d acc: (%f,%f,%f)\n", tp, id, type,
                                                    targetVars[bucketStart+target].a.x,
                                                    targetVars[bucketStart+target].a.y,
                                                    targetVars[bucketStart+target].a.z
                                                    );
#endif
    }

  }// for each target part
}
#endif

__global__ void EwaldTopKernel(GravityParticleData *particleTable, int nPart);
__global__ void EwaldBottomKernel(GravityParticleData *particleTable, int nPart);

extern unsigned int timerHandle; 

void EwaldHostMemorySetup(EwaldData *h_idata, int nParticles, int nEwhLoop) {
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
  assert(nEwhLoop <= NEWH);

  workRequest topKernel;
  dataInfo *particleTableInfo, *cachedDataInfo, *ewtInfo;  

  topKernel.dimGrid = dim3(numBlocks); 
  topKernel.dimBlock = dim3(BLOCK_SIZE); 
  topKernel.smemSize = 0; 

  topKernel.nBuffers = 2; 

  /* schedule two buffers for transfer to the GPU */ 
  topKernel.bufferInfo = 
    (dataInfo *) malloc(topKernel.nBuffers * sizeof(dataInfo));

  particleTableInfo = &(topKernel.bufferInfo[PARTICLE_TABLE]);
  particleTableInfo->bufferID = NUM_GRAVITY_BUFS + PARTICLE_TABLE; 
  particleTableInfo->transferToDevice = true; 
  particleTableInfo->transferFromDevice = false; 
  particleTableInfo->freeBuffer = false; 
  particleTableInfo->hostBuffer = h_idata->p; 
  particleTableInfo->size = n * sizeof(GravityParticleData); 

  cachedDataInfo = &(topKernel.bufferInfo[EWALD_READ_ONLY_DATA]); 
  cachedDataInfo->bufferID = NUM_GRAVITY_BUFS + EWALD_READ_ONLY_DATA;
  cachedDataInfo->transferToDevice = false; 
  cachedDataInfo->transferFromDevice = false; 
  cachedDataInfo->freeBuffer = false; 
  cachedDataInfo->hostBuffer = h_idata->cachedData; 
  cachedDataInfo->size = sizeof(EwaldReadOnlyData); 

  topKernel.callbackFn = NULL;
  topKernel.id = TOP_EWALD_KERNEL; 
#ifdef CUDA_INSTRUMENT_WRS
  topKernel.chareIndex = myIndex;
  topKernel.compType = topKernel.id;
  topKernel.compPhase = phase; 
#endif
  enqueue(&topKernel); 
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

  particleTableInfo = &(bottomKernel.bufferInfo[PARTICLE_TABLE]);
  particleTableInfo->bufferID = NUM_GRAVITY_BUFS + PARTICLE_TABLE; 
  particleTableInfo->transferToDevice = false; 
  particleTableInfo->transferFromDevice = true; 
  particleTableInfo->freeBuffer = true; 
  particleTableInfo->hostBuffer = h_idata->p; 
  particleTableInfo->size = n * sizeof(GravityParticleData); 

  cachedDataInfo = &(bottomKernel.bufferInfo[EWALD_READ_ONLY_DATA]); 
  cachedDataInfo->bufferID = NUM_GRAVITY_BUFS + EWALD_READ_ONLY_DATA; 
  cachedDataInfo->transferToDevice = false; 
  cachedDataInfo->transferFromDevice = false; 
  cachedDataInfo->freeBuffer = false; 
  cachedDataInfo->hostBuffer = h_idata->cachedData; 
  cachedDataInfo->size = sizeof(EwaldReadOnlyData); 

  ewtInfo = &(bottomKernel.bufferInfo[EWALD_TABLE]); 
  ewtInfo->bufferID = NUM_GRAVITY_BUFS + EWALD_TABLE; 
  ewtInfo->transferToDevice = false; 
  ewtInfo->transferFromDevice = false; 
  ewtInfo->freeBuffer = false; 
  ewtInfo->hostBuffer = h_idata->ewt; 
  ewtInfo->size = nEwhLoop * sizeof(EwtData); 

  bottomKernel.callbackFn = cb;
  bottomKernel.id = BOTTOM_EWALD_KERNEL; 
#ifdef CUDA_INSTRUMENT_WRS
  bottomKernel.chareIndex = myIndex;
  bottomKernel.compType = bottomKernel.id;
  bottomKernel.compPhase = phase; 
#endif
  enqueue(&bottomKernel); 
#ifdef CUDA_VERBOSE_KERNEL_ENQUEUE
  printf("[%d] in EwaldHost, enqueued BotKernel\n", myIndex);
#endif
}

__global__ void EwaldTopKernel(GravityParticleData *particleTable,
           int nPart) {

  GravityParticleData *p;

  cudatype alphan;
  cudatype fPot, ax, ay, az;
  cudatype x, y, z, r2, dir, dir2, a; 
  cudatype xdif, ydif, zdif; 
  cudatype g0, g1, g2, g3, g4, g5;
  cudatype Q2, Q2mirx, Q2miry, Q2mirz, Q2mir, Qta; 
  int id, ix, iy, iz, bInHole, bInHolex, bInHolexy;

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

  id = blockIdx.x * BLOCK_SIZE + threadIdx.x;
  if (id >= nPart) {
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

  ax = 0.0f;
  ay = 0.0f;
  az = 0.0f;

#ifdef HEXADECAPOLE
  xdif = p->position_x - momQuad->cmx; 
  ydif = p->position_y - momQuad->cmy; 
  zdif = p->position_z - momQuad->cmz;
  fPot = momQuad->totalMass*cachedData->k1;
#else
  xdif = p->position_x - mom->cmx; 
  ydif = p->position_y - mom->cmy; 
  zdif = p->position_z - mom->cmz;
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

  p->potential += fPot;
  p->acceleration_x += ax;
  p->acceleration_y += ay;
  p->acceleration_z += az;

}

__global__ void EwaldBottomKernel(GravityParticleData *particleTable, int nPart) {

  int i, id;
  cudatype hdotx, c, s, fPot; 
  cudatype ax, ay, az; 
  cudatype tempEwt; 
  cudatype xdif, ydif, zdif; 
  GravityParticleData *p; 
  MultipoleMomentsData *mom; 

  mom = &(cachedData->mm); 

  id = blockIdx.x * BLOCK_SIZE + threadIdx.x;
  if (id >= nPart) {
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
