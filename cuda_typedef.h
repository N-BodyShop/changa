#ifndef CUDA_TYPEDEF_H_
#define CUDA_TYPEDEF_H_

/** @file cuda_typedef.h
 *
 * Definitions of types for the CUDA port.
 */

#include "cosmoType.h"

/** @brief floating point type on the GPU */
typedef float cudatype;
/** @brief floating point type on the host */
typedef cosmoType hosttype;

// set these to appropriate values (in millions)
// local
#define NODE_INTERACTIONS_PER_REQUEST_L 1.0
#define PART_INTERACTIONS_PER_REQUEST_L 0.1
// remote, no-resume
#define NODE_INTERACTIONS_PER_REQUEST_RNR 1.0
#define PART_INTERACTIONS_PER_REQUEST_RNR 0.1
// remote, resume
#define NODE_INTERACTIONS_PER_REQUEST_RR 1.0
#define PART_INTERACTIONS_PER_REQUEST_RR 0.1

#ifdef CUDA_STATS
#define CUDA_SER_TREE 9900
#define CUDA_SER_LIST 9901

#define CUDA_LOCAL_NODE_KERNEL 9910
#define CUDA_REMOTE_NODE_KERNEL 9911
#define CUDA_REMOTE_RESUME_NODE_KERNEL 9912
#define CUDA_LOCAL_PART_KERNEL 9913
#define CUDA_REMOTE_PART_KERNEL 9914
#define CUDA_REMOTE_RESUME_PART_KERNEL 9915

#endif

// TODO: Fix small phase code
#define TP_LARGE_PHASE_THRESHOLD_DEFAULT 0.0
#define AVG_SOURCE_PARTICLES_PER_ACTIVE 10 

/** @brief 3D vector of cudatype.
 */
typedef struct CudaVector3D{
  cudatype x,y,z;
#if __cplusplus && !defined __CUDACC__
  inline CudaVector3D& operator=(Vector3D<hosttype> &a){
    x = a.x;
    y = a.y;
    z = a.z;
    return *this;
  }
  
  inline Vector3D<hosttype> operator+(Vector3D<hosttype> &v){
    return Vector3D<hosttype>(x + v.x, y + v.y, z + v.z);
  }

  CudaVector3D(Vector3D<hosttype> &o){
    x = o.x;
    y = o.y;
    z = o.z;
  }
  
  CudaVector3D(){}
#endif
}CudaVector3D;

#ifdef GPU_LOCAL_TREE_WALK
typedef struct CudaSphere {
  /// The origin of this sphere
  CudaVector3D origin;
  /// The radius of this sphere
  cudatype radius;
}CudaSphere;

enum CudaNodeType {
    CudaInvalid = 1,
    CudaBucket = 2,
    CudaInternal = 3,
    CudaBoundary = 4,
    CudaNonLocal = 5,
    CudaEmpty = 6,
    CudaTop = 7,
    CudaNonLocalBucket = 8,
    CudaCached = 9,
    CudaCachedBucket = 10,
    CudaCachedEmpty = 11
};
#endif //GPU_LOCAL_TREE_WALK

/** @brief Version of MultipoleMoments using cudatype
 */
typedef struct CudaMultipoleMoments{
  cudatype radius;
  cudatype soft;
  cudatype totalMass;
  CudaVector3D cm;

#ifdef GPU_LOCAL_TREE_WALK
  // We need tree node's spatial and tree structural information to do GPU
  // tree walk. The spatial info is used for force computation, and the
  // structural data is needed in tree traversal.
  CudaVector3D lesser_corner;
  CudaVector3D greater_corner;
  int bucketStart;
  int bucketSize;
  int particleCount;
  int nodeArrayIndex;
  int children[2];
  int type;
#endif //GPU_LOCAL_TREE_WALK

#ifdef HEXADECAPOLE
  cudatype xx, xy, xz, yy, yz;
  cudatype xxx,xyy,xxy,yyy,xxz,yyz,xyz;
  cudatype xxxx,xyyy,xxxy,yyyy,xxxz,yyyz,xxyy,xxyz,xyyz;
#else
  cudatype xx, xy, xz, yy, yz, zz;
#endif

#if __cplusplus && !defined __CUDACC__
  CudaMultipoleMoments(){}
  CudaMultipoleMoments(MultipoleMoments &mom){
    *this = mom;
  }
  inline CudaMultipoleMoments& operator=(MultipoleMoments &m){
    radius = m.radius;
    soft = m.soft;
    totalMass = m.totalMass;

    cm = m.cm;
#if ! defined(HEXADECAPOLE)
    xx = m.xx;
    xy = m.xy;
    xz = m.xz;
    yy = m.yy;
    yz = m.yz;
    zz = m.zz;
#else
    xx = m.mom.xx;
    yy = m.mom.yy;
    xy = m.mom.xy;
    xz = m.mom.xz;
    yz = m.mom.yz;
    xxx = m.mom.xxx;
    xyy = m.mom.xyy;
    xxy = m.mom.xxy;
    yyy = m.mom.yyy;
    xxz = m.mom.xxz;
    yyz = m.mom.yyz;
    xyz = m.mom.xyz;
    xxxx = m.mom.xxxx;
    xyyy = m.mom.xyyy;
    xxxy = m.mom.xxxy;
    yyyy = m.mom.yyyy;
    xxxz = m.mom.xxxz;
    yyyz = m.mom.yyyz;
    xxyy = m.mom.xxyy;
    xxyz = m.mom.xxyz;
    xyyz = m.mom.xyyz;
#endif

    return *this;
  }
#endif
}CudaMultipoleMoments;

#ifdef GPU_LOCAL_TREE_WALK
struct CUDATreeNode
{
  cudatype radius;
  cudatype soft;
  cudatype totalMass;
  CudaVector3D cm;

  int bucketStart;
  int bucketSize;
  int particleCount;
  int children[2];
  int type;
};

struct CUDABucketNode
{
  cudatype radius;
  cudatype soft;
  cudatype totalMass;
  CudaVector3D cm;

  CudaVector3D lesser_corner;
  CudaVector3D greater_corner;
};
#endif //GPU_LOCAL_TREE_WALK

/** @brief Bucket of particles on the interaction list for the GPU.
 */
typedef struct ILPart{
  /** Index of particles on GPU. */
  int index;
  /** Encodes the replica for periodic boundaries */
  int off;
  /** Number of particles in the bucket. */
  int num;

#if __cplusplus && !defined __CUDACC__
  ILPart() {}
  //ILPart() : index(-1), numParticles(-1) {}
  ILPart(int i, int o, int n) : index(i), off(o), num(n) {}
#endif
}ILPart;

/** @brief Cell on the interaction list for the GPU
 */
typedef struct ILCell{
  /** Index of this cell's moments on the GPU. */
  int index;
  /** Encodes the replica for periodic boundaries */
  int offsetID;
#if __cplusplus && !defined __CUDACC__
  ILCell() {}
  //ILCell() :index(-1), offsetID(-1) {}
  ILCell(int ind, int off) : index(ind), offsetID(off) {}
#endif
}ILCell;

/**
 *  @brief Particle data needed on the GPU to calculate gravity.
 */
typedef struct CompactPartData{
  cudatype mass;
  cudatype soft;
  CudaVector3D position;

#ifdef GPU_LOCAL_TREE_WALK
  int nodeId;
#endif

#if __cplusplus && !defined __CUDACC__
  CompactPartData(){}
  CompactPartData(ExternalGravityParticle &egp){
    *this = egp;
  }
  CompactPartData(cudatype m, cudatype s, Vector3D<hosttype> &rr) : mass(m), soft(s), position(rr){}

  inline CompactPartData& operator=(ExternalGravityParticle &gp){
    mass = gp.mass;
    soft = gp.soft;
    position = gp.position;
    return *this;
  }
#endif
}CompactPartData;

/**
 *  @brief Particle data that gets calculated by the GPU.
 */
typedef struct VariablePartData{
  CudaVector3D a;
  cudatype potential;
  cudatype dtGrav;

#ifdef GPU_LOCAL_TREE_WALK
  int numOfNodesTraversed;
  int numOfParticlesTraversed;
#endif //GPU_LOCAL_TREE_WALK
}VariablePartData;


#endif /* CUDA_TYPEDEF_H_*/
