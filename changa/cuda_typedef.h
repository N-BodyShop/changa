#ifndef CUDA_TYPEDEF_H_
#define CUDA_TYPEDEF_H_

typedef float cudatype;
typedef double hosttype;

// FIXME - set these to appropriate values
// remote, no-resume/local
#define NODE_INTERACTIONS_PER_REQUEST_RNR 3
#define PART_INTERACTIONS_PER_REQUEST_RNR 3
// remote, resume
#define NODE_INTERACTIONS_PER_REQUEST_RR 3
#define PART_INTERACTIONS_PER_REQUEST_RR 3

#define MAX_PARTICLES_PER_BUCKET 10
#define THREADS_PER_BLOCK 128

typedef struct CudaVector3D{
  cudatype x,y,z;
#if __cplusplus && !defined __CUDACC__
  inline CudaVector3D& operator=(Vector3D<hosttype> &a){
    x = a.x;
    y = a.y;
    z = a.z;
    return *this;
  }
  CudaVector3D(Vector3D<hosttype> &o){
    x = o.x;
    y = o.y;
    z = o.z;
  }
  CudaVector3D(){}
#endif
}CudaVector3D;

typedef struct CudaMultipoleMoments{
  cudatype radius;
  cudatype soft;
  cudatype totalMass;
  CudaVector3D cm;
  cudatype xx, xy, xz, yy, yz, zz;
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
    xx = m.xx;
    xy = m.xy;
    xz = m.xz;
    yy = m.yy;
    yz = m.yz;
    zz = m.zz;

    return *this;
  }
#endif
}CudaMultipoleMoments;

typedef struct ILPart{
  int index;
  int off;
  int num;

#if __cplusplus && !defined __CUDACC__
  ILPart() {}
  //ILPart() : index(-1), numParticles(-1) {}
  ILPart(int i, int o, int n) : index(i), off(o), num(n) {}
#endif
}ILPart;

typedef struct ILCell{
  int index;
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
