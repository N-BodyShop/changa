#ifndef CUDA_TYPEDEF_H_
#define CUDA_TYPEDEF_H_

/** @file cuda_typedef.h
 *
 * Definitions of types for the CUDA port.
 */

#include "cosmoType.h"

#ifdef HEXADECAPOLE
/** @brief floating point type on the GPU */
typedef double cudatype;
#else
/** @brief floating point type on the GPU */
typedef float cudatype;
#endif
/** @brief floating point type on the host */
typedef cosmoType hosttype;

// set these to appropriate values (in millions)
// local
#define NODE_INTERACTIONS_PER_REQUEST_L 0.1
#define PART_INTERACTIONS_PER_REQUEST_L 0.1
// remote, no-resume
#define NODE_INTERACTIONS_PER_REQUEST_RNR 0.1
#define PART_INTERACTIONS_PER_REQUEST_RNR 0.1
// remote, resume
#define NODE_INTERACTIONS_PER_REQUEST_RR 0.1
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

#define TP_LARGE_PHASE_THRESHOLD_DEFAULT 0.3
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

/** @brief Version of MultipoleMoments using cudatype
 */
typedef struct CudaMultipoleMoments{
  cudatype radius;
  cudatype soft;
  cudatype totalMass;
  CudaVector3D cm;

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

typedef struct CompactPartData{
  cudatype mass;
  cudatype soft;
  CudaVector3D position;
#if defined CUDA_EMU_KERNEL_NODE_PRINTS || defined CUDA_EMU_KERNEL_PART_PRINTS
  int tp, id;
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

typedef struct VariablePartData{
  CudaVector3D a;
  cudatype potential;
}VariablePartData;

typedef struct PartData{
  CompactPartData core;
  CudaVector3D a;
  cudatype potential;
  cudatype dtGrav;

#if __cplusplus && !defined __CUDACC__
  PartData(){}
  PartData(GravityParticle &gp){
    *this = gp;
  }
  PartData(CompactPartData &cpd, Vector3D<hosttype> &_a, cudatype p, cudatype dtg) : core(cpd), a(_a), potential(p), dtGrav(dtg) {}

  inline PartData& operator=(GravityParticle &gp){
    core = gp;
    a.x = 0.0;
    a.y = 0.0;
    a.z = 0.0;
    potential = 0.0;
    dtGrav = 0.0;
    return *this;
  }
#endif
}PartData;

#if 0
#ifdef __cplusplus
// work request data structures

template <class T>
struct CudaGroupRequest{
	T *intlist;
	/*
	CkVec<int> bucketMarkers;
	CkVec<int> bucketStarts;
	CkVec<int> buckets;
	CkVec<int> bucketSizes;
	*/

	int *bucketMarkers;
	int *bucketStarts;
	int *buckets;
	int *bucketSizes;

	int numInteractions; // number of interactions in intlist
	int numBucketsPlusOne; // number of buckets involved in the work request
	bool lastIsPartial; // is the last bucket partially computed?

	TreePiece *tp;
	State *state;

	virtual void cleanUp(){
		delete [] intlist;
	}

	CudaGroupRequest(int tpBuckets){
		bucketMarkers = new int[tpBuckets+1];
		bucketStarts = new int [tpBuckets];
		buckets = new int [tpBuckets];
		bucketSizes  = new int[tpBuckets];
	}
};

template <class S, class T, int size>
struct CudaGroupMissedRequest : public CudaGroupRequest <T> {
	S *missed;
	int numMissed;

	void cleanUp(){
		CudaGroupRequest::cleanUp();
		delete [] missed;
	}

	CudaGroupMissedRequest(int b){
		CudaGroupRequest(b);
		missed = new S[size];
	}
};

#endif // __cplusplus
#endif

#endif /* CUDA_TYPEDEF_H_*/
