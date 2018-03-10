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

#ifdef CAMBRIDGE

__device__ inline bool __attribute__(( always_inline ))
CUDA_intersect(CUDABucketNode &b, CudaSphere &s);

__device__ inline bool __attribute__(( always_inline ))
CUDA_intersect(CudaSphere &s1, CudaSphere &s2);

__device__ inline bool __attribute__(( always_inline ))
CUDA_contains(const CudaSphere &s, const CudaVector3D &v);

__device__ inline bool __attribute__(( always_inline ))
CUDA_contained(const CUDATreeNode &b, const CudaSphere &s);

__device__ inline int __attribute__(( always_inline ))
CUDA_openSoftening(CUDATreeNode &node, CUDABucketNode &myNode);

__device__ inline CudaVector3D __attribute__(( always_inline ))
CUDA_openCriterionNode(CUDATreeNode &node,
                    CUDABucketNode &myNode,
                    int localIndex,
                    cudatype theta,
                    cudatype thetaMono);

__device__ inline void __attribute__(( always_inline ))
CUDA_SPLINEQ(cudatype invr, cudatype r2, cudatype twoh, cudatype& a,
       cudatype& b,cudatype& c,cudatype& d);

__device__ inline void __attribute__(( always_inline ))
CUDA_SPLINE(cudatype r2, cudatype twoh, cudatype &a, cudatype &b); 

__device__ inline int __attribute__(( always_inline ))
CUDA_OptAction(int fakeOpen, int nodetype);

#endif

#endif  /* CUDAMoments_h */
