#ifndef CELL_TYPEDEF_H_
#define CELL_TYPEDEF_H_

typedef float cellSPEtype;

typedef struct CellVector3D {
  cellSPEtype x, y, z;
} CellVector3D;

typedef struct CellMultipoleMoments {
    /// A physical size for this multipole expansion, calculated by an external function using some other information
    cellSPEtype radius;
    cellSPEtype soft;        /* Effective softening */

    /// The total mass represented by this expansion
    cellSPEtype totalMass;
    /// The center of mass (zeroth order multipole)
    CellVector3D cm;
    //Tensor for higher order moments goes here
    cellSPEtype xx, xy, xz, yy, yz, zz;
} CellMultipoleMoments;

typedef struct CellExternalGravityParticle {
  cellSPEtype mass;
  cellSPEtype soft;
  CellVector3D position;
} CellExternalGravityParticle;

typedef struct CellGravityParticle {
  ExternalGravityParticle core;
  CellVector3D treeAcceleration;
  cellSPEtype potential;
  cellSPEtype dtGrav;
} CellGravityParticle;

#ifdef __cplusplus
inline CellVector3D& operator+(CellVector3D &p, Vector3d<double> &o) {
  p.x += o.x;
  p.y += o.y;
  p.z += o.z;
  return p;
}

inline CellVector3D& operator=(CellVector3D &p, Vector3D<double> &o) {
  p.x = o.x;
  p.y = o.y;
  p.z = o.z;
  return p;
}

inline Vector3D<double>& operator+(Vector3D<double> &p, CellVector3D &o) {
  p.x += o.x;
  p.y += o.y;
  p.z += o.z;
  return p;
}

inline CellExternalGravityParticle& operator=(CellExternalGravityParticle &p, ExternalGravityParticle &o) {
  p.mass = o.mass;
  p.soft = o.soft;
  p.position = o.position;
  return p;
}

inline CellGravityParticle& operator=(CellGravityParticle &p, GravityParticle &o) {
  p.core = o;
  p.treeAcceleration.x = 0.0;
  p.treeAcceleration.y = 0.0;
  p.treeAcceleration.z = 0.0;
  p.potential = 0.0;
  p.dtGrav = 0.0;
  return p;
}

inline CellMultipoleMoments& operator=(CellMultipoleMoments &p, MultipoleMoments &o) {
  p.radius = o.radius;
  p.soft = o.soft;
  p.totalMass = o.totalMass;
  p.cm = o.cm;
  p.xx = o.xx;
  p.xy = o.xy;
  p.xz = o.xz;
  p.yy = o.yy;
  p.yz = o.yz;
  p.zz = o.zz;
  return p;
}

class CellRequest {
  CellGravityParticle *activeData;
  void *roData;
  GravityParticle **particles;
  int numActiveData;

  CellRequest(CellGrabityParticle *ad, int num, void *ro) : activeData(ad), roData(ro), numActiveData(num) {}
};

class CellGroupRequest {
  TreePiece *tp;
  int bucket;
  GravityParticle **particles;
  
  CellGroupRequest(TreePiece *_tp, int b, GravityParticle **p) : tp(_tp), bucket(b), particle(p) {}
};

#endif

#endif /*CELL_TYPEDEF_H_*/
