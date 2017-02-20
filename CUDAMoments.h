#ifndef CUDAMoments_h
#define CUDAMoments_h 1

#include "cuda_typedef.h"

/** CUDA version of momEvalFmomrcm. */
#ifdef CUDA_UNIT_TEST
__global__ void
#else
__device__ inline void __attribute__(( always_inline ))
#endif
CUDA_momEvalFmomrcm(const CudaMultipoleMoments* _m,
                    const CudaVector3D* _r,
                    cudatype dir,
                    CudaVector3D* out,
                    cudatype* pot);

#endif  /* CUDAMoments_h */
