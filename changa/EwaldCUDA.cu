

#ifdef _WIN32
#define NOMINMAX
#endif

#include "EwaldCUDA.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "wr.h"
// #include <cutil.h>

#include "cuda_runtime.h"

#define NEWH 80
#define BLOCK_SIZE 128

__device__ __constant__ EwaldReadOnlyData cachedData[1];
__device__ __constant__ EwtData ewt[sizeof(EwtData) * NEWH];  

__global__ void EwaldTopKernel(GravityParticleData *particleTable/*, EwaldReadOnlyData *cachedData*/);
__global__ void EwaldBottomKernel(GravityParticleData *particleTable, EwtData *ewt/*, EwaldReadOnlyData *cachedData*/);

void EwaldHostMemorySetup(EwaldData *h_idata, int nParticles, int nEwhLoop) {
  
  cudaMallocHost((void **) &(h_idata->p), 
		 (nParticles+1) * sizeof(GravityParticleData));
  cudaMallocHost((void **) &(h_idata->ewt),
		 nEwhLoop * sizeof(EwtData));
  cudaMallocHost((void **) &(h_idata->cachedData), 
		 sizeof(EwaldReadOnlyData)); 
}

void EwaldHostMemoryFree(EwaldData *h_idata) {
  cudaFreeHost(h_idata->p); 
  cudaFreeHost(h_idata->ewt); 
  cudaFreeHost(h_idata->cachedData); 
}

void EwaldHost(EwaldData *h_idata, void *cb, int myIndex) {

  int n = h_idata->cachedData->n;
  int numBlocks = (int) ceilf((float)n/BLOCK_SIZE);
  int nEwhLoop = h_idata->cachedData->nEwhLoop;

  workRequest topKernel;
  dataInfo *particleTableInfo, *cachedDataInfo, *ewtInfo;  
  // topKernel = (workRequest*) malloc(sizeof(workRequest));

  topKernel.dimGrid = dim3(numBlocks); 
  topKernel.dimBlock = dim3(BLOCK_SIZE); 
  topKernel.smemSize = 0; 

  topKernel.nBuffers = 2; 

  /* schedule two buffers for transfer to the GPU */ 
  topKernel.bufferInfo = 
    (dataInfo *) malloc(topKernel.nBuffers * sizeof(dataInfo));
  
  particleTableInfo = &(topKernel.bufferInfo[0]);
  particleTableInfo->bufferID = BUFFERS_PER_CHARE * myIndex + PARTICLE_TABLE; 
  particleTableInfo->transferToDevice = YES; 
  particleTableInfo->transferFromDevice = NO; 
  particleTableInfo->freeBuffer = NO; 
  particleTableInfo->hostBuffer = h_idata->p; 
  particleTableInfo->size = (n+1) * sizeof(GravityParticleData); 

  cachedDataInfo = &(topKernel.bufferInfo[1]); 
  cachedDataInfo->bufferID = BUFFERS_PER_CHARE * myIndex + EWALD_READ_ONLY_DATA;
  cachedDataInfo->transferToDevice = NO; 
  cachedDataInfo->transferFromDevice = NO; 
  cachedDataInfo->freeBuffer = NO; 
  cachedDataInfo->hostBuffer = h_idata->cachedData; 
  cachedDataInfo->size = sizeof(EwaldReadOnlyData); 

  topKernel.callbackFn = NULL;
  topKernel.id = TOP_KERNEL; 
  enqueue(wrQueue, &topKernel); 

  workRequest bottomKernel; 
  // bottomKernel = (workRequest*) malloc(sizeof(workRequest)); 

  bottomKernel.dimGrid = dim3(numBlocks); 
  bottomKernel.dimBlock = dim3(BLOCK_SIZE); 
  bottomKernel.smemSize = 0; 

  bottomKernel.nBuffers = 3; 

  bottomKernel.bufferInfo = 
    (dataInfo *) malloc(bottomKernel.nBuffers * sizeof(dataInfo));
  
  particleTableInfo = &(bottomKernel.bufferInfo[0]);
  particleTableInfo->bufferID = BUFFERS_PER_CHARE * myIndex + PARTICLE_TABLE; 
  particleTableInfo->transferToDevice = NO; 
  particleTableInfo->transferFromDevice = YES; 
  particleTableInfo->freeBuffer = YES; 
  particleTableInfo->hostBuffer = h_idata->p; 
  particleTableInfo->size = (n+1) * sizeof(GravityParticleData); 
  
  cachedDataInfo = &(bottomKernel.bufferInfo[1]); 
  cachedDataInfo->bufferID = BUFFERS_PER_CHARE * myIndex + EWALD_READ_ONLY_DATA; 
  cachedDataInfo->transferToDevice = NO; 
  cachedDataInfo->transferFromDevice = NO; 
  cachedDataInfo->freeBuffer = NO; 
  cachedDataInfo->hostBuffer = h_idata->cachedData; 
  cachedDataInfo->size = sizeof(EwaldReadOnlyData); 
  
  ewtInfo = &(bottomKernel.bufferInfo[2]); 
  ewtInfo->bufferID = BUFFERS_PER_CHARE * myIndex + EWALD_TABLE; 
  ewtInfo->transferToDevice = NO; 
  ewtInfo->transferFromDevice = NO; 
  ewtInfo->freeBuffer = NO; 
  ewtInfo->hostBuffer = h_idata->ewt; 
  ewtInfo->size = nEwhLoop * sizeof(EwtData); 

  bottomKernel.callbackFn = cb;
  bottomKernel.id = BOTTOM_KERNEL; 

  enqueue(wrQueue, &bottomKernel); 
}

__global__ void EwaldTopKernel(GravityParticleData *particleTable/*, EwaldReadOnlyData *cachedData*/) {

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

__global__ void EwaldBottomKernel(GravityParticleData *particleTable/*, EwtData *ewt, __constant__ EwaldReadOnlyData *cachedData*/) {

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
  ** 	Without trig = (10,14)
  **	    Trig est.	 = 2*(6,11)  same as 1/sqrt scoring.
  **		Total        = (22,36)
  **					 = 58
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

void kernelSelect(workRequest *wr) {

  switch (wr->id) {
  case TOP_KERNEL: 

    cudaMemcpyToSymbol(cachedData, wr->bufferInfo[EWALD_READ_ONLY_DATA].hostBuffer, sizeof(EwaldReadOnlyData)); 
    EwaldTopKernel<<<wr->dimGrid, wr->dimBlock, wr->smemSize, kernel_stream>>>
      ((GravityParticleData *)devBuffers[wr->bufferInfo[PARTICLE_TABLE].bufferID]
/*, 
										   (EwaldReadOnlyData *)devBuffers[wr->bufferInfo[EWALD_READ_ONLY_DATA].bufferID]*/);     
    break; 
  case BOTTOM_KERNEL:
    cudaMemcpyToSymbol(ewt, wr->bufferInfo[EWALD_TABLE].hostBuffer, NEWH * sizeof(EwtData)); 
    EwaldBottomKernel<<<wr->dimGrid, wr->dimBlock, 
      wr->smemSize, kernel_stream>>>
      ((GravityParticleData *)devBuffers[wr->bufferInfo[PARTICLE_TABLE].bufferID]/*, 
       (EwtData *)devBuffers[wr->bufferInfo[EWALD_TABLE].bufferID], 
								    (EwaldReadOnlyData *)devBuffers[wr->bufferInfo[EWALD_READ_ONLY_DATA].bufferID]*/);
    break; 
  default:
    printf("error: id %d not valid\n", wr->id); 
    break; 
  }
}


