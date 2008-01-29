#ifndef CELL_TYPEDEF_H_
#define CELL_TYPEDEF_H_

#include "spert.h"

#ifdef __cplusplus
#include "Vector3D.h"
#endif

typedef float cellSPEtype;

typedef struct CellVector3D {
  cellSPEtype x, y, z;
#ifdef __cplusplus
  inline CellVector3D& operator=(Vector3D<double> &o) {
    x = o.x;
    y = o.y;
    z = o.z;
    return *this;
  }
  inline CellVector3D& operator+=(Vector3D<double> &o) {
    x += o.x;
    y += o.y;
    z += o.z;
    return *this;
  }
#endif
} CellVector3D;

typedef struct CellVector4D {
  cellSPEtype x, y, z;
  cellSPEtype padding;
#ifdef __cplusplus
  inline CellVector4D& operator=(Vector3D<double> &o) {
    x = o.x;
    y = o.y;
    z = o.z;
    return *this;
  }
  inline CellVector4D& operator+=(Vector3D<double> &o) {
    x += o.x;
    y += o.y;
    z += o.z;
    return *this;
  }
#endif
} CellVector4D;

typedef struct CellMultipoleMoments {
  /// A physical size for this multipole expansion, calculated by an external function using some other information
  cellSPEtype radius;
  cellSPEtype soft;        /* Effective softening */

  /// The center of mass (zeroth order multipole)
  CellVector3D cm;
  /// The total mass represented by this expansion
  cellSPEtype totalMass;
  //Tensor for higher order moments goes here
  cellSPEtype xx, xy, xz, yy, yz, zz;
#ifdef __cplusplus
  inline CellMultipoleMoments& operator=(MultipoleMoments &o) {
    radius = o.radius;
    soft = o.soft;
    totalMass = o.totalMass;
    cm = o.cm;
    xx = o.xx;
    xy = o.xy;
    xz = o.xz;
    yy = o.yy;
    yz = o.yz;
    zz = o.zz;
    return *this;
  }
#endif
} CellMultipoleMoments;

typedef struct CellExternalGravityParticle {
  CellVector4D position;
  cellSPEtype mass;
  cellSPEtype soft;
#ifdef __cplusplus
  inline CellExternalGravityParticle& operator=(ExternalGravityParticle &o) {
    mass = o.mass;
    soft = o.soft;
    position = o.position;
    return *this;
  }
#endif
} CellExternalGravityParticle;

typedef struct CellGravityParticle {
  CellExternalGravityParticle core;
  cellSPEtype potential;
  cellSPEtype dtGrav;
  CellVector4D treeAcceleration;
#ifdef __cplusplus
  inline CellGravityParticle& operator=(GravityParticle &o) {
    core = o;
    treeAcceleration.x = 0.0;
    treeAcceleration.y = 0.0;
    treeAcceleration.z = 0.0;
    potential = 0.0;
    dtGrav = 0.0;
    return *this;
  }
#endif
} CellGravityParticle;

typedef struct CellContainer {
  int numInt;
  int numExt;
  void *data;
} CellContainer;

typedef struct CellEWT {
  cellSPEtype hx,hy,hz;
  cellSPEtype hCfac,hSfac;
#ifdef __cplusplus
  inline CellEWT& operator=(EWT &o) {
    hx = o.hx;
    hy = o.hy;
    hz = o.hz;
    hCfac = o.hCfac;
    hSfac = o.hSfac;
    return *this;
  }
#endif
} CellEWT;

typedef struct CellEwaldContainer {
  CellMultipoleMoments rootMoments;
  cellSPEtype fEwCut;
  cellSPEtype fPeriod;
  int numPart;
  int nEwhLoop;
  int nReps;
  int pad1,pad2,pad3;
} CellEwaldContainer;

#ifdef __cplusplus
inline CellVector3D& operator+(CellVector3D &p, Vector3D<double> &o) {
  p.x += o.x;
  p.y += o.y;
  p.z += o.z;
  return p;
}

inline Vector3D<double>& operator+(Vector3D<double> &p, CellVector3D &o) {
  p.x += o.x;
  p.y += o.y;
  p.z += o.z;
  return p;
}

inline CellVector4D& operator+(CellVector4D &p, Vector3D<double> &o) {
  p.x += o.x;
  p.y += o.y;
  p.z += o.z;
  return p;
}

inline Vector3D<double>& operator+(Vector3D<double> &p, CellVector4D &o) {
  p.x += o.x;
  p.y += o.y;
  p.z += o.z;
  return p;
}

class CellRequest {
public:
  CellGravityParticle *activeData;
  void *roData;
  GravityParticle **particles;
  int numActiveData;
  TreePiece *tp;

 CellRequest(CellGravityParticle *ad, int num, void *ro, GravityParticle **part, TreePiece *t) : activeData(ad), roData(ro), numActiveData(num), particles(part), tp(t) {}
};

class CellGroupRequest {
public:
  TreePiece *tp;
  int bucket;
  GravityParticle **particles;
  
  CellGroupRequest(TreePiece *_tp, int b, GravityParticle **p) : tp(_tp), bucket(b), particles(p) {}
};

class CellEwaldRequest {
 public:
  cellSPEtype *woData;
  CellEwaldContainer *roData;
  TreePiece *tp;
  GravityParticle **particles;
  int numActiveData;
  int firstBucket, lastBucket;

  CellEwaldRequest(cellSPEtype *output, int num, CellEwaldContainer *input, GravityParticle **part, TreePiece *t, int first, int last) : woData(output), numActiveData(num), roData(input), particles(part), tp(t), firstBucket(first), lastBucket(last) {}
};

class CellComputation {
 public:
  CProxyElement_TreePiece owner;
  dummyMsg *msg;

  CellComputation() {}
  CellComputation(CProxyElement_TreePiece tp, dummyMsg *m) : owner(tp), msg(m) {}
};
#endif

#endif /*CELL_TYPEDEF_H_*/
