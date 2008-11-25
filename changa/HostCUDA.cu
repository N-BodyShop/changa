#ifdef _WIN32
#define NOMINMAX
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "wr.h"
// #include <cutil.h>

#include "cuda_runtime.h"

#define BLOCK_SIZE 128

extern workRequestQueue *wrQueue;
extern void **devBuffers;
extern cudaStream_t kernel_stream;

__constant__ constantData[88];

void DataManagerTransfer(CudaMultipoleMoments *moments, int nMoments, CompactPartData *compactParts, int nCompactParts, int bufferID) {

	int numBlocks = (int) ceilf((float)n/BLOCK_SIZE);

	//workRequest *transferKernel = (workRequest*) malloc(sizeof(workRequest));
	workRequest transferKernel;
	dataInfo *momentBuffer, *partCoreBuffer;

	transferKernel.dimGrid = dim3(numBlocks);
	transferKernel.dimBlock = dim3(BLOCK_SIZE);
	transferKernel.smemSize = 0;

	transferKernel.nBuffers = DM_TRANSFER_NBUFFERS;

	/* schedule two buffers for transfer to the GPU */
	transferKernel.bufferInfo = (dataInfo *) malloc(transferKernel.nBuffers * sizeof(dataInfo));

	momentBuffer = &(transferKernel.bufferInfo[POST_PREFETCH_MOMENTS]);
	momentBuffer->bufferID = bufferID+POST_PREFETCH_MOMENTS;
	momentBuffer->transferToDevice = YES;
	momentBuffer->transferFromDevice = NO;
	momentBuffer->freeBuffer = NO;
	momentBuffer->hostBuffer = moments;
	momentBuffer->size = (nMoments) * sizeof(CudaMultipoleMoments);

	partCoreBuffer = &(topKernel.bufferInfo[POST_PREFETCH_PARTICLE_CORES]);
	partCoreBuffer->bufferID = bufferID+POST_PREFETCH_PARTICLE_CORES;
	partCoreBuffer->transferToDevice = YES;
	partCoreBuffer->transferFromDevice = NO;
	partCoreBuffer->freeBuffer = NO;
	partCoreBuffer->hostBuffer = compactParts;
	partCoreBuffer->size = (nCompactParts)*sizeof(CompactPartData);

	transferKernel.callbackFn = 0;
	transferKernel.id = DM_TRANSFER;
	transferKernel.executing = 0;
	enqueue(wrQueue, &transferKernel);

}

void TreePieceCellListDataTransferLocal(CellListData *data){
	int numBlocks = (int) ceilf((float)n/BLOCK_SIZE);

	workRequest gravityKernel;
	dataInfo *buffer, *partCoreBuffer;

	// FIXME - size properly
	gravityKernel.dimGrid = dim3(numBlocks);
	gravityKernel.dimBlock = dim3(BLOCK_SIZE);
	gravityKernel.smemSize = 0;

	gravityKernel.nBuffers = TP_GRAVITY_LOCAL_NBUFFERS;

	/* schedule buffers for transfer to the GPU */
	gravityKernel.bufferInfo = (dataInfo *) malloc(gravityKernel.nBuffers * sizeof(dataInfo));

	TreePieceCellListDataTransferBasic(data, &gravityKernel);

	gravityKernel.callbackFn = 0;
	gravityKernel.id = TP_GRAVITY_LOCAL;
	gravityKernel.executing = 0;
	enqueue(wrQueue, &gravityKernel);
}

void TreePieceCellListDataTransferRemote(CellListData *data){
	int numBlocks = (int) ceilf((float)n/BLOCK_SIZE);

	workRequest gravityKernel;
	dataInfo *buffer, *partCoreBuffer;

	// FIXME - size properly
	gravityKernel.dimGrid = dim3(numBlocks);
	gravityKernel.dimBlock = dim3(BLOCK_SIZE);
	gravityKernel.smemSize = 0;

	gravityKernel.nBuffers = TP_NODE_GRAVITY_REMOTE_NBUFFERS;

	/* schedule buffers for transfer to the GPU */
	gravityKernel.bufferInfo = (dataInfo *) malloc(gravityKernel.nBuffers * sizeof(dataInfo));

	TreePieceCellListDataTransferBasic(data, &gravityKernel);

	gravityKernel.callbackFn = 0;
	gravityKernel.id = TP_GRAVITY_REMOTE;
	gravityKernel.executing = 0;
	enqueue(wrQueue, &gravityKernel);
}

void TreePieceCellListDataTransferRemoteResume(CellListData *data, CudaMultipoleMoments *missedMoments, int numMissedMoments){
	int numBlocks = (int) ceilf((float)n/BLOCK_SIZE);

	//workRequest *gravityKernel = (workRequest*) malloc(sizeof(workRequest));
	workRequest gravityKernel;
	dataInfo *buffer, *partCoreBuffer;

	// FIXME - size properly
	gravityKernel.dimGrid = dim3(numBlocks);
	gravityKernel.dimBlock = dim3(BLOCK_SIZE);
	gravityKernel.smemSize = 0;

	gravityKernel.nBuffers = TP_NODE_GRAVITY_REMOTE_RESUME_NBUFFERS;

	/* schedule buffers for transfer to the GPU */
	gravityKernel.bufferInfo = (dataInfo *) malloc(gravityKernel.nBuffers * sizeof(dataInfo));

	TreePieceCellListDataTransferBasic(data, &gravityKernel);

	buffer = &(gravityKernel.bufferInfo[MISSED_MOMENTS]);
	// FIXME - needs to be unique?
	buffer->bufferID = bufferID+MISSED_MOMENTS;
	buffer->transferToDevice = YES;
	buffer->transferFromDevice = NO;
	buffer->freeBuffer = YES;
	buffer->hostBuffer = missedMoments;
	buffer->size = (numMissedMoments) * sizeof(int);

	gravityKernel.callbackFn = 0;
	gravityKernel.id = TP_GRAVITY_REMOTE_RESUME;
	gravityKernel.executing = 0;
	enqueue(wrQueue, &gravityKernel);
}

void TreePieceCellListDataTransferBasic(CellListData *data, workRequest *gravityKernel){
	dataInfo *buffer;
	int numInteractions = data->numInteractions;
	int numBucketsPlusOne = data->numBucketsPlusOne;
	ILCell *cellList = data->cellList;
	int *cellListBucketMarkers = data->partListBucketMarkers;
	int *bucketStarts = data->bucketStarts;
	int *bucketSizes = data->bucketSizes;

	buffer = &(gravityKernel->bufferInfo[ILCELL]);
	// FIXME - needs to be unique?
	buffer->bufferID = bufferID+ILCELL;
	buffer->transferToDevice = YES;
	buffer->transferFromDevice = NO;
	buffer->freeBuffer = YES;
	buffer->hostBuffer = cellList;
	buffer->size = (numInteractions) * sizeof(ILCell);

	buffer = &(gravityKernel->bufferInfo[NODE_INTLIST_BUCKET_MARKERS]);
	// FIXME - needs to be unique?
	buffer->bufferID = bufferID+NODE_INTLIST_BUCKET_MARKERS;
	buffer->transferToDevice = YES;
	buffer->transferFromDevice = NO;
	buffer->freeBuffer = YES;
	buffer->hostBuffer = cellListBucketMarkers;
	buffer->size = (numBucketsPlusOne) * sizeof(int);

	buffer = &(gravityKernel->bufferInfo[NODE_INTLIST_BUCKET_STARTS]);
	// FIXME - needs to be unique?
	buffer->bufferID = bufferID+NODE_INTLIST_BUCKET_STARTS;
	buffer->transferToDevice = YES;
	buffer->transferFromDevice = NO;
	buffer->freeBuffer = YES;
	buffer->hostBuffer = bucketStarts;
	buffer->size = (numBucketsPlusOne-1) * sizeof(int);

	buffer = &(gravityKernel->bufferInfo[NODE_INTLIST_BUCKET_SIZES]);
	// FIXME - needs to be unique?
	buffer->bufferID = bufferID+NODE_INTLIST_BUCKET_SIZES;
	buffer->transferToDevice = YES;
	buffer->transferFromDevice = NO;
	buffer->freeBuffer = YES;
	buffer->hostBuffer = bucketSizes;
	buffer->size = (numBucketsPlusOne-1) * sizeof(int);
}

void TreePiecePartListDataTransferLocal(PartListData *data){
	int numBlocks = (int) ceilf((float)n/BLOCK_SIZE);

	workRequest gravityKernel;
	dataInfo *buffer, *partCoreBuffer;

	// FIXME - size properly
	gravityKernel.dimGrid = dim3(numBlocks);
	gravityKernel.dimBlock = dim3(BLOCK_SIZE);
	gravityKernel.smemSize = 0;

	gravityKernel.nBuffers = TP_GRAVITY_LOCAL_NBUFFERS;

	/* schedule buffers for transfer to the GPU */
	gravityKernel.bufferInfo = (dataInfo *) malloc(gravityKernel.nBuffers * sizeof(dataInfo));

	TreePiecePartListDataTransferBasic(data, gravityKernel);

	gravityKernel.callbackFn = 0;
	gravityKernel.id = TP_PART_GRAVITY_LOCAL;
	gravityKernel.executing = 0;
	enqueue(wrQueue, &gravityKernel);
}

void TreePiecePartListDataTransferRemote(PartListData *data, CompactPartData *missedParts, int numMissedParts){
	int numBlocks = (int) ceilf((float)n/BLOCK_SIZE);

	//workRequest *gravityKernel = (workRequest*) malloc(sizeof(workRequest));
	workRequest gravityKernel;
	dataInfo *buffer, *partCoreBuffer;

	// FIXME - size properly
	gravityKernel.dimGrid = dim3(numBlocks);
	gravityKernel.dimBlock = dim3(BLOCK_SIZE);
	gravityKernel.smemSize = 0;

	gravityKernel.nBuffers = TP_PART_GRAVITY_REMOTE_NBUFFERS;

	/* schedule buffers for transfer to the GPU */
	gravityKernel.bufferInfo = (dataInfo *) malloc(gravityKernel.nBuffers * sizeof(dataInfo));

	TreePiecePartListDataTransferBasic(data, gravityKernel);

	buffer = &(gravityKernel.bufferInfo[MISSED_PARTS]);
	// FIXME - needs to be unique?
	buffer->bufferID = bufferID+MISSED_PARTS;
	buffer->transferToDevice = YES;
	buffer->transferFromDevice = NO;
	buffer->freeBuffer = YES;
	buffer->hostBuffer = missedParts;
	buffer->size = (numMissedParts) * sizeof(int);

	gravityKernel.callbackFn = 0;
	gravityKernel.id = TP_PART_GRAVITY_REMOTE;
	gravityKernel.executing = 0;
	enqueue(wrQueue, &gravityKernel);
}

void TreePiecePartListDataTransferRemoteResume(PartListData *data, CompactPartData *missedParts, int numMissedParts){
	int numBlocks = (int) ceilf((float)n/BLOCK_SIZE);

	//workRequest *gravityKernel = (workRequest*) malloc(sizeof(workRequest));
	workRequest gravityKernel;
	dataInfo *buffer, *partCoreBuffer;

	// FIXME - size properly
	gravityKernel.dimGrid = dim3(numBlocks);
	gravityKernel.dimBlock = dim3(BLOCK_SIZE);
	gravityKernel.smemSize = 0;

	gravityKernel.nBuffers = TP_PART_GRAVITY_REMOTE_RESUME_NBUFFERS;

	/* schedule buffers for transfer to the GPU */
	gravityKernel.bufferInfo = (dataInfo *) malloc(gravityKernel.nBuffers * sizeof(dataInfo));

	TreePiecePartListDataTransferBasic(data, gravityKernel);

	buffer = &(gravityKernel.bufferInfo[MISSED_PARTS]);
	// FIXME - needs to be unique?
	buffer->bufferID = bufferID+MISSED_PARTS;
	buffer->transferToDevice = YES;
	buffer->transferFromDevice = NO;
	buffer->freeBuffer = YES;
	buffer->hostBuffer = missedParts;
	buffer->size = (numMissedParts) * sizeof(int);

	gravityKernel.callbackFn = 0;
	gravityKernel.id = TP_PART_GRAVITY_REMOTE_RESUME;
	gravityKernel.executing = 0;
	enqueue(wrQueue, &gravityKernel);
}

void TreePiecePartListDataTransferBasic(PartListData *data, workRequest *gravityKernel){
	dataInfo *buffer;

	int numInteractions = data->numInteractions;
	int numBucketsPlusOne = data->numBucketsPlusOne;
	ILPart *partList = data->partList;
	int *partListBucketMarkers = data->partListBucketMarkers;
	int *bucketStarts = data->bucketStarts;
	int *bucketSizes = data->bucketSizes;

	buffer = &(gravityKernel->bufferInfo[ILPART]);
	// FIXME - needs to be unique?
	buffer->bufferID = bufferID+ILPART;
	buffer->transferToDevice = YES;
	buffer->transferFromDevice = NO;
	buffer->freeBuffer = YES;
	buffer->hostBuffer = partList;
	buffer->size = (numInteractions) * sizeof(ILPart);

	buffer = &(gravityKernel->bufferInfo[PART_INTLIST_BUCKET_MARKERS]);
	// FIXME - needs to be unique?
	buffer->bufferID = bufferID+PART_INTLIST_BUCKET_MARKERS;
	buffer->transferToDevice = YES;
	buffer->transferFromDevice = NO;
	buffer->freeBuffer = YES;
	buffer->hostBuffer = partListBucketMarkers;
	buffer->size = (numBucketsPlusOne) * sizeof(int);

	buffer = &(gravityKernel->bufferInfo[PART_INTLIST_BUCKET_STARTS]);
	// FIXME - needs to be unique?
	buffer->bufferID = bufferID+PART_INTLIST_BUCKET_STARTS;
	buffer->transferToDevice = YES;
	buffer->transferFromDevice = NO;
	buffer->freeBuffer = YES;
	buffer->hostBuffer = bucketStarts;
	buffer->size = (numBucketsPlusOne-1) * sizeof(int);

	buffer = &(gravityKernel->bufferInfo[PART_INTLIST_BUCKET_SIZES]);
	// FIXME - needs to be unique?
	buffer->bufferID = bufferID+PART_INTLIST_BUCKET_SIZES;
	buffer->transferToDevice = YES;
	buffer->transferFromDevice = NO;
	buffer->freeBuffer = YES;
	buffer->hostBuffer = bucketSizes;
	buffer->size = (numBucketsPlusOne-1) * sizeof(int);
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

extern void DeleteHostMoments(CudaMultipoleMoments *array);
extern void DeleteHostParticles(CompactPartData *array);
static boolean dmTransferDone = NO;

void kernelSelect(workRequest *wr) {

  switch (wr->id) {
  case DM_TRANSFER:
	  dmTransferDone = YES;
	  DeleteHostMoments(wr->bufferInfo[POST_PREFETCH_MOMENTS].hostBuffer);
	  DeleteHostParticles(wr->bufferInfo[POST_PREFETCH_PARTICLE_CORES].hostBuffer);
	  break;

  case TP_GRAVITY_LOCAL:
	  // FIXME - fix arguments
	  if(dmTransferDone){
		  //GravityKernel<<<wr->dimGrid, wr->dimBlock, wr->smemSize, kernel_stream>>>
		  //((GravityParticleData *)devBuffers[wr->bufferInfo[PARTICLE_TABLE].bufferID],
			//(EwtData *)devBuffers[wr->bufferInfo[EWALD_TABLE].bufferID]);
	  }
	  else{
		  // fix buffer transfer flags and re-enqueue
		  NoTransferEnqueueNodeBasic(wr);
	  }
    break;

  case TP_PART_GRAVITY_LOCAL:
  	  // FIXME - fix arguments
  	  if(dmTransferDone){
  		  //GravityKernel<<<wr->dimGrid, wr->dimBlock, wr->smemSize, kernel_stream>>>
  		  //((GravityParticleData *)devBuffers[wr->bufferInfo[PARTICLE_TABLE].bufferID],
  			//(EwtData *)devBuffers[wr->bufferInfo[EWALD_TABLE].bufferID]);
  	  }
  	  else{
  		  // fix buffer transfer flags and re-enqueue
  		  NoTransferEnqueuePartBasic(wr);
  	  }
      break;

  case TP_GRAVITY_REMOTE:
	  if(dmTransferDone){
		  // invoke kernel
	  }
	  else{
		  // fix buffer transfer flags and re-enqueue
		  NoTransferEnqueueNodeBasic(wr);
	  }
	  break;

  case TP_PART_GRAVITY_REMOTE:
	  if(dmTransferDone){
		  // invoke kernel
	  }
	  else{
		  // fix buffer transfer flags and re-enqueue
		  NoTransferEnqueuePartBasic(wr);

	  }
	  break;

	  // FIXME - remote resume for parts and nodes
	  // after that, set callbacks
	  // finally, do book-keeping, remaining mindful of the lastBucketComplete flag in the state

  default:
    printf("error: id %d not valid\n", wr->id);
    break;
  }
}

void NoTransferEnqueueNodeBasic(workRequest *gravityKernel){
	dataInfo *buffer;

	// all buffer transfer parameters set to NO because there already exists a workRequest that will do all
	// of this
	workRequest *newGravityKernel = (workRequest*) malloc(sizeof(workRequest));
	*newGravityKernel = *gravityKernel;

	buffer = &(newGravityKernel->bufferInfo[ILCELL]);
	buffer->transferToDevice = NO;
	buffer->transferFromDevice = NO;
	buffer->freeBuffer = NO;

	buffer = &(newGravityKernel->bufferInfo[NODE_INTLIST_BUCKET_MARKERS]);
	buffer->transferToDevice = NO;
	buffer->transferFromDevice = NO;
	buffer->freeBuffer = NO;

	buffer = &(newGravityKernel->bufferInfo[NODE_INTLIST_BUCKET_STARTS]);
	buffer->transferToDevice = NO;
	buffer->transferFromDevice = NO;
	buffer->freeBuffer = NO;

	buffer = &(newGravityKernel->bufferInfo[NODE_INTLIST_BUCKET_SIZES]);
	buffer->transferToDevice = NO;
	buffer->transferFromDevice = NO;
	buffer->freeBuffer = NO;

	enqueue(wrQueue, newGravityKernel);
}

void NoTransferEnqueuePartBasic(workRequest *gravityKernel){
	dataInfo *buffer;

	// all buffer transfer parameters set to NO because there already exists a workRequest that will do all
	// of this
	workRequest *newGravityKernel = (workRequest*) malloc(sizeof(workRequest));
	*newGravityKernel = *gravityKernel;

	buffer = &(newGravityKernel->bufferInfo[ILPART]);
	buffer->transferToDevice = NO;
	buffer->transferFromDevice = NO;
	buffer->freeBuffer = NO;

	buffer = &(newGravityKernel->bufferInfo[PART_INTLIST_BUCKET_MARKERS]);
	buffer->transferToDevice = NO;
	buffer->transferFromDevice = NO;
	buffer->freeBuffer = NO;

	buffer = &(newGravityKernel->bufferInfo[PART_INTLIST_BUCKET_STARTS]);
	buffer->transferToDevice = NO;
	buffer->transferFromDevice = NO;
	buffer->freeBuffer = NO;

	buffer = &(newGravityKernel->bufferInfo[PART_INTLIST_BUCKET_SIZES]);
	buffer->transferToDevice = NO;
	buffer->transferFromDevice = NO;
	buffer->freeBuffer = NO;

	enqueue(wrQueue, gravityKernel);
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
		int numBucketsPlusOne){
  // FIXME - register usage, etc.

  CudaVector3D priv_acc[THREADS_PER_BLOCK];
  cudatype priv_pot[THREADS_PER_BLOCK];
  CudaMultipoleMoments shared_moments[THREADS_PER_BLOCK];
  CompactPartData shared_particle_cores[quotceil(THREADS_PER_BLOCK, MAX_THREADS_PER_GROUP)];


  // each block is given a bucket to compute
  // each thread in the block computes an interaction of a particle with a node
  int bucket = blockIdx.x;
  int start = ilmarks[bucket];
  int end = ilmarks[bucket+1];
  int bucketSize = bucketSizes[bucket];
  int bucketStart = bucketStarts[bucket];

  // length of cell interaction list for this bucket
  int llen = ilmarks[bucket+1] - ilmarks[bucket];

  CudaVector3D r;
  cudatype rsq;
  cudatype twoh, a, b, c, d;
  CudaMultipoleMoments m;

  int node, particle, group;
  group = GROUP(thread);
  particle = group;

  int ngroups = quotceil(THREADS_PER_BLOCK,MAX_THREADS_PER_GROUP);
  int groupSize = (group < ngroups-1 ? MAX_THREADS_PER_GROUP : THREADS_PER_BLOCK - MAX_THREADS_PER_GROUP*(ngroups-1));

#ifdef __DEVICE_EMULATION__
  if(blockIdx.x == 0){
    printf("t: %d, blen: %d, llen: %d\n", threadIdx.x, blen, llen);
    printf("group: %d, particle: %d, ngroups: %d, groupSize: %d\n", group, particle, ngroups, groupSize);
  }
#endif

  while(particle < bucketSize){
	  node = GROUP_INDEX(thread);

#ifdef __DEVICE_EMULATION__
    if(blockIdx.x == 0){
      printf("shared_particle_cores[%d] = particleCores[%d]\n", group, bmarks[start+bucket]+particle);
    }
#endif

    if(node == 0){
    	shared_particle_cores[group] = particleCores[bucketStart+particle];
    }

    priv_acc[thread].x = 0;
    priv_acc[thread].y = 0;
    priv_acc[thread].z = 0;
    priv_pot[thread] = 0;

    while(node < llen){
#ifdef __DEVICE_EMULATION__
      if(blockIdx.x == 0){
        printf("t: %d, particle: %d, node: %d\n", threadIdx.x, particle, node);
      }
#endif
      // FIXME - ought to keep ilmarks[blockIdx.x] in a register
      // all threads in a group access different mem locations
      // however, threads accessing the same node in different
      // groups access the same memory location
      // therefore, there ought to be 8 accesses for every 4 groups
      // FIXME - can this be reduced somehow?
#ifdef __DEVICE_EMULATION__
      if(blockIdx.x == 0){
        printf("shared_moments[%d] = moments[%d]\n", threadIdx.x, ils[ilmarks[bucket]+node].index);
      }
#endif

      shared_moments[threadIdx.x] =
                          moments[ils[ilmarks[bucket]+node].index];
      // above analysis applies here as well. again, ilmarks[blockIdx.x]
      // should be in a register.
#ifdef __DEVICE_EMULATION__
      if(blockIdx.x == 0){
        printf("int offsetID = ils[%d].offsetID\n", ilmarks[blockIdx.x]+node);
      }
#endif
      int offsetID = ils[ilmarks[bucket]+node].offsetID;
      m = shared_moments[thread];

      r.x = shared_particle_cores[group].position.x -
                            ((((offsetID >> 22) & 0x7)-3)*fperiod + m.cm.x);
      r.y = shared_particle_cores[group].position.y -
                            ((((offsetID >> 25) & 0x7)-3)*fperiod + m.cm.y);
      r.z = shared_particle_cores[group].position.z -
                            ((((offsetID >> 28) & 0x7)-3)*fperiod + m.cm.z);

      // FIXME - special functions here?
      rsq = r.x*r.x + r.y*r.y + r.z*r.z;
      twoh = m.soft + shared_particle_cores[group].soft;
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
       cudatype qirx = m.xx*r.x + m.xy*r.y + m.xz*r.z;
       cudatype qiry = m.xy*r.x + m.yy*r.y + m.yz*r.z;
       cudatype qirz = m.xz*r.x + m.yz*r.y + m.zz*r.z;
       cudatype qir = 0.5*(qirx*r.x + qiry*r.y + qirz*r.z);
       cudatype tr = 0.5*(m.xx + m.yy + m.zz);
       cudatype qir3 = b*m.totalMass + d*qir - c*tr;


       priv_pot[thread] -= m.totalMass * a + c*qir - b*tr;

       priv_acc[thread].x -= qir3*r.x - c*qirx;
       priv_acc[thread].y -= qir3*r.y - c*qiry;
       priv_acc[thread].z -= qir3*r.z - c*qirz;
      }
      node += groupSize;
    }
    // before moving on to the next particle,
    // we must add up the accelerations due to different threads in the
    // same group and commit them to global memory.
    // since there are at most only MAX_THREADS_PER_GROUP-1 additions
    // to perform, we do them serially.
    //__syncthreads();
    cudatype sum = 0.0;
    for(int i = MAX_THREADS_PER_GROUP*group; i < MAX_THREADS_PER_GROUP*group + groupSize; i++){
      sum += priv_acc[i].x;
    }
    sum += particleVars[bmarks[start+bucket]+particle].a.x;
    particleVars[bmarks[start+bucket]+particle].a.x = sum;

    sum = 0.0;
    for(int i = MAX_THREADS_PER_GROUP*group; i < MAX_THREADS_PER_GROUP*group + groupSize; i++){
      sum += priv_acc[i].y;
    }
    sum += particleVars[bmarks[start+bucket]+particle].a.y;
    particleVars[bmarks[start+bucket]+particle].a.y = sum;

    sum = 0.0;
    for(int i = MAX_THREADS_PER_GROUP*group; i < MAX_THREADS_PER_GROUP*group + groupSize; i++){
      sum += priv_acc[i].z;
    }
    sum += particleVars[bmarks[start+bucket]+particle].a.z;
    particleVars[bmarks[start+bucket]+particle].a.z = sum;

    sum = 0.0;
    for(int i = MAX_THREADS_PER_GROUP*group; i < MAX_THREADS_PER_GROUP*group + groupSize; i++){
      sum += priv_pot[i];
    }
    sum += particleVars[bmarks[start+bucket]+particle].potential;
    particleVars[bmarks[start+bucket]+particle].potential = sum;

    particle += ngroups;
  }
}

__global__ void particleGravityComputation(
                                   CompactPartData *particleCores,
                                   VariablePartData *particleVars,
				   int *bmarks,
                                   ILPart *ils,
                                   int *ilmarks,
                                   int start, int end,
                                   int x, cudatype fperiod){

  CudaVector3D priv_acc[THREADS_PER_BLOCK];
  cudatype priv_pot[THREADS_PER_BLOCK];
  CompactPartData shared_source_cores[THREADS_PER_BLOCK];
  CompactPartData shared_target_cores[quotceil(THREADS_PER_BLOCK, MAX_THREADS_PER_GROUP)];

  // size of bucket in particles
  // add start to blockIdx.x here because bmarks is defined
  // over entire set of buckets
  int blen = bmarks[start+blockIdx.x+1] - bmarks[start+blockIdx.x];

  // length of cell interaction list for this bucket
  // no need to do the same here because ilmarks is defined
  // only for the current chunk of buckets
  int llen = ilmarks[blockIdx.x+1] - ilmarks[blockIdx.x];
#ifdef __DEVICE_EMULATION__
#endif

  CudaVector3D r;
  cudatype rsq;
  cudatype twoh, a, b;

  int source_particle, target_particle;
  int group;
  target_particle = threadIdx.x/MAX_THREADS_PER_GROUP;
  group = threadIdx.x/MAX_THREADS_PER_GROUP;

  int ngroups = quotceil(THREADS_PER_BLOCK, MAX_THREADS_PER_GROUP);
  int groupSize = (target_particle < ngroups-1 ? MAX_THREADS_PER_GROUP : THREADS_PER_BLOCK - MAX_THREADS_PER_GROUP*(ngroups-1));

  while(target_particle < blen){
    // All threads in group, and in warp, access same locations,
    // so no use asking only one thread to do the load - all
    // requests ought to be merged, without any SIMD divergence.
    shared_target_cores[group] = particleCores[bmarks[start+blockIdx.x]+target_particle];
    // FIXME - no synchronization needed here?
    priv_acc[threadIdx.x].x = 0;
    priv_acc[threadIdx.x].y = 0;
    priv_acc[threadIdx.x].z = 0;

    priv_pot[threadIdx.x] = 0;

    source_particle = threadIdx.x%MAX_THREADS_PER_GROUP;
    while(source_particle < llen){
#ifdef __DEVICE_EMULATION__
#endif
      // again, accesses to ilmarks[blockIdx.x] by every thread
      // in block should be coalesced
      // FIXME - ought to keep ilmarks[blockIdx.x] in a register
      // in addition, due to different 'node' values in threads,
      // all threads in a group access different mem locations
      // however, threads accessing the same node in different
      // groups access the same memory location
      // therefore, there ought to be 8 accesses for every 4 groups
      // FIXME - can this be reduced somehow?
      shared_source_cores[threadIdx.x] =
                          particleCores[ils[ilmarks[blockIdx.x]+source_particle].index];
      // above analysis applies here as well. again, ilmarks[blockIdx.x]
      // should be in a register.
      int oid = ils[ilmarks[blockIdx.x]+source_particle].off;

      r.x = (((oid >> 22) & 0x7)-3)*fperiod +
            shared_source_cores[threadIdx.x].position.x -
            shared_target_cores[group].position.x;

      r.y = (((oid >> 25) & 0x7)-3)*fperiod +
            shared_source_cores[threadIdx.x].position.y -
            shared_target_cores[group].position.y;

      r.z = (((oid >> 28) & 0x7)-3)*fperiod +
            shared_source_cores[threadIdx.x].position.z -
            shared_target_cores[group].position.z;

      // FIXME - special functions here?
      rsq = r.x*r.x + r.y*r.y + r.z*r.z;
      twoh = shared_source_cores[threadIdx.x].soft + shared_target_cores[group].soft;
      if(rsq != 0){

       //SPLINE(rsq, twoh, a, b);
       //SPLINE(r2, twoh, a, b);
       //expanded below:
	cudatype r1, u,dih,dir;
        // FIXME -  can we use rsq and twoh*twoh below instead?
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

       priv_pot[threadIdx.x] -= shared_source_cores[threadIdx.x].mass * a;

       priv_acc[threadIdx.x].x += r.x*b*shared_source_cores[threadIdx.x].mass;
       priv_acc[threadIdx.x].y += r.y*b*shared_source_cores[threadIdx.x].mass;
       priv_acc[threadIdx.x].z += r.z*b*shared_source_cores[threadIdx.x].mass;
      }
      source_particle += groupSize;
    }
    // before moving on to the next particle,
    // we must add up the accelerations due to different threads in the
    // same group and commit them to global memory.
    // since there are at most only MAX_THREADS_PER_GROUP-1 additions
    // to perform, we do them serially.
    cudatype sum = 0.0;
    for(int i = MAX_THREADS_PER_GROUP*group; i < MAX_THREADS_PER_GROUP*group + groupSize; i++){
      sum += priv_acc[i].x;
    }
    sum += particleVars[bmarks[start+blockIdx.x]+target_particle].a.x;
    particleVars[bmarks[start+blockIdx.x]+target_particle].a.x = sum;

    sum = 0.0;
    for(int i = MAX_THREADS_PER_GROUP*group; i < MAX_THREADS_PER_GROUP*group + groupSize; i++){
      sum += priv_acc[i].y;
    }
    sum += particleVars[bmarks[start+blockIdx.x]+target_particle].a.y;
    particleVars[bmarks[start+blockIdx.x]+target_particle].a.y = sum;

    sum = 0.0;
    for(int i = MAX_THREADS_PER_GROUP*group; i < MAX_THREADS_PER_GROUP*group + groupSize; i++){
      sum += priv_acc[i].z;
    }
    sum += particleVars[bmarks[start+blockIdx.x]+target_particle].a.z;
    particleVars[bmarks[start+blockIdx.x]+target_particle].a.z = sum;

    sum = 0.0;
    for(int i = MAX_THREADS_PER_GROUP*group; i < MAX_THREADS_PER_GROUP*group + groupSize; i++){
      sum += priv_pot[i];
    }
    sum += particleVars[bmarks[start+blockIdx.x]+target_particle].potential;
    particleVars[bmarks[start+blockIdx.x]+target_particle].potential = sum;

    target_particle +=  ngroups;
  }
}

