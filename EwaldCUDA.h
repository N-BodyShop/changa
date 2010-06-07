#ifndef _EWALD_CUDA_H_
#define _EWALD_CUDA_H_ 

#include "HostCUDA.h"

/* defines for Hybrid API buffer indices */ 

#define PARTICLE_TABLE        0 
#define EWALD_READ_ONLY_DATA  1
#define EWALD_TABLE           2

#define BUFFERS_PER_CHARE     3

#define NEWH 80
#define BLOCK_SIZE 128
#define NUM_GRAVITY_BUFS 10

typedef struct {
  float hx, hy, hz; 
  float hCfac, hSfac; 
} EwtData;

typedef struct {
  float xx, xy, xz, yy, yz, zz; 
  float totalMass; 
  float cmx, cmy, cmz; 

} MultipoleMomentsData; 

typedef struct {
  MultipoleMomentsData mm; 
  
  int n, nReps, nEwReps, nEwhLoop;
  float L, fEwCut, alpha, alpha2, k1, ka, fEwCut2, fInner2;

} EwaldReadOnlyData; 

typedef struct {
  float position_x,position_y,position_z;
  float acceleration_x, acceleration_y, acceleration_z; 
  float potential;
} GravityParticleData;

typedef struct {
  GravityParticleData *p;
  EwtData *ewt; 
  EwaldReadOnlyData *cachedData;
} EwaldData; 

void EwaldHostMemorySetup(EwaldData *h_idata, int nParticles, int nEwhLoop, void *cb); 
void EwaldHostMemoryFree(EwaldData *h_idata); 
#ifdef CUDA_INSTRUMENT_WRS
void EwaldHost(EwaldData *h_idata, void *cb, int myIndex, char phase); 
#else
void EwaldHost(EwaldData *h_idata, void *cb, int myIndex); 
#endif

__global__ void EwaldTopKernel(GravityParticleData *particleTable);
__global__ void EwaldBottomKernel(GravityParticleData *particleTable);

#endif

