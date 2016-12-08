#ifndef _EWALD_CUDA_H_
#define _EWALD_CUDA_H_ 

#include "HostCUDA.h"

/* defines for Hybrid API buffer indices */ 

#define EWALD_READ_ONLY_DATA  0
#define EWALD_TABLE           1
/* This has to be larger than the first two because there will be one
 * of these for each Chare on the node. */
#define PARTICLE_TABLE        2 

/* Defines for place in the bufferInfo array */
#define PARTICLE_TABLE_IDX    0 
#define EWALD_READ_ONLY_DATA_IDX  1
#define EWALD_TABLE_IDX           2

#define BUFFERS_PER_CHARE     3

#define NEWH 80
#define BLOCK_SIZE 128

/* See "defines for Hybrid API buffer indices" in HostCUDA.h for this
 * number. */
#define NUM_GRAVITY_BUFS 5

/** @brief Data for the Ewald h loop in the CUDA kernel
 */
typedef struct {
  cudatype hx, hy, hz; 
  cudatype hCfac, hSfac; 
} EwtData;

/** @brief CUDA version of complete MultipoleMoments for Ewald
 */
typedef struct {
#ifndef HEXADECAPOLE
  cudatype xx, xy, xz, yy, yz, zz; 
#endif
  cudatype totalMass; 
  cudatype cmx, cmy, cmz; 

} MultipoleMomentsData; 

/** @brief CUDA version of MOMC for Ewald
 */
typedef struct {
  cudatype m;
  cudatype xx,yy,xy,xz,yz;
  cudatype xxx,xyy,xxy,yyy,xxz,yyz,xyz;
  cudatype xxxx,xyyy,xxxy,yyyy,xxxz,yyyz,xxyy,xxyz,xyyz;
  cudatype zz;
  cudatype xzz,yzz,zzz;
  cudatype xxzz,xyzz,xzzz,yyzz,yzzz,zzzz;
} MomcData;

/** @brief Parameters and data for Ewald in the CUDA kernel
 */
typedef struct {
  MultipoleMomentsData mm; 
  MomcData momcRoot;
  
  int n, nReps, nEwReps, nEwhLoop;
  cudatype L, fEwCut, alpha, alpha2, k1, ka, fEwCut2, fInner2;

} EwaldReadOnlyData; 

/** @brief Particle data for the CUDA Ewald kernels
 */

typedef struct {
  int *ewaldmarkers;
  EwtData *ewt; 
  EwaldReadOnlyData *cachedData;  
} EwaldData; 

void EwaldHostMemorySetup(EwaldData *h_idata, int NumActiveParticles, int nEwhLoop); 
void EwaldHostMemoryFree(EwaldData *h_idata); 
#ifdef CUDA_INSTRUMENT_WRS
void EwaldHost(EwaldData *h_idata, void *cb, int myIndex, char phase); 
#else
void EwaldHost(EwaldData *h_idata, void *cb, int myIndex); 
#endif

__global__ void EwaldTopKernel(CompactPartData *particleCores, VariablePartData *particleVars, int *markers, int nPart);
__global__ void EwaldBottomKernel(CompactPartData *particleCores, VariablePartData *particleVars, int *markers, int nPart);

#endif

