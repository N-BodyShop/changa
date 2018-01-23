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
cuda_intersect(CUDATreeNode &b, CudaSphere &s);

__device__ inline bool __attribute__(( always_inline ))
cuda_intersect(CudaSphere &s1, CudaSphere &s2);

__device__ inline bool __attribute__(( always_inline ))
cuda_contains(const CudaSphere &s, const CudaVector3D &v);

__device__ inline bool __attribute__(( always_inline ))
cuda_contained(const CUDATreeNode &b, const CudaSphere &s);

__device__ inline int __attribute__(( always_inline ))
cuda_openSoftening(CUDATreeNode &node, CUDATreeNode &myNode);

__device__ inline int __attribute__(( always_inline ))
cuda_encodeOffset(int reqID, int x, int y, int z);

__device__ inline int __attribute__(( always_inline ))
cuda_reEncodeOffset(int reqID, int offsetID);

__device__ inline CudaVector3D __attribute__(( always_inline ))
cuda_decodeOffset(int reqID, CudaVector3D fPeriod);

__device__ inline CudaVector3D __attribute__(( always_inline ))
cuda_openCriterionNode(CUDATreeNode &node,
                    CUDATreeNode &myNode,
                    int localIndex,
                    cudatype theta,
                    cudatype thetaMono);

__device__ inline void __attribute__(( always_inline ))
cuda_SPLINEQ(cudatype invr, cudatype r2, cudatype twoh, cudatype& a,
       cudatype& b,cudatype& c,cudatype& d);

__device__ inline void __attribute__(( always_inline ))
cuda_SPLINE(cudatype r2, cudatype twoh, cudatype &a, cudatype &b); 

__device__ inline int __attribute__(( always_inline ))
cuda_OptAction(int fakeOpen, int nodetype);

#endif

#endif  /* CUDAMoments_h */
