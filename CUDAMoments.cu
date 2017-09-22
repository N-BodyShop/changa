#ifdef CUDA_UNIT_TEST
# include "CUDAMoments.h"
#else
# include "cuda_typedef.h"
#endif

#ifdef CAMBRIDGE
#include "codes.h"
#endif

#ifdef CUDA_UNIT_TEST
__global__ void
#else
__device__ inline void __attribute__(( always_inline ))
#endif
CUDA_momEvalFmomrcm(const CudaMultipoleMoments* _m,
                    const CudaVector3D* _r,
                    cudatype dir,
                    CudaVector3D* out,
                    cudatype* pot)
{
  /* Added in these values for standalone function. */
  CudaMultipoleMoments m = *_m;
  CudaVector3D r = *_r;

  /* The following code is adapted from from `momEvalFmomrcm` in
     "moments.c"; the changes make the code work within the inputs
     available here, and (hopefullly) make the code a little more
     readable. */
  const cudatype onethird = 1.0 / 3.0;

  /* -> Build the reciprocal-of-radius and scaling-factor values. */
  cudatype
    /* in `momEvalFmomrcm`, `u` is a parameter, and the value passed a
       MultipoleMoments::radius instance (see point(s) of call at
       `nodeBucketForce` in "gravity.h").  `momEvalFmomrcm` also
       multiplies the parameter by `dir` prior to use. */
    u = dir * m.radius;

  /* -> Build the "g" terms, whose purpose is probably apparent to those
     who actually understand the math...  */
  cudatype
    g0 = dir,
    g2 = 3 * dir * u * u,
    g3 = 5 * g2 * u,
    g4 = 7 * g3 * u;


  /* -> "Calculate the trace-free distance terms." */
  cudatype
    x = r.x * dir,
    y = r.y * dir,
    z = r.z * dir,
    xx = 0.5f * x * x,
    xy = x * y,
    xz = x * z,
    yy = 0.5f * y * y,
    yz = y * z,
    zz = 0.5f * z * z,
    xxx = x * (onethird*xx - zz),
    xxz = z * (xx - onethird * zz),
    yyy = y * (onethird*yy - zz),
    yyz = z*(yy - onethird*zz);

  /* replace intermediates used above with their "final" values... */
  xx -= zz;
  yy -= zz;

  /* ...and finish with the trace-free terms. */
  cudatype
    xxy = y * xx,
    xyy = x * yy,
    xyz = xy * z;

  /* -> "Now calculate the interaction up to Hexadecapole order." */
  cudatype
    tx = g4 * ( m.xxxx*xxx + m.xyyy*yyy + m.xxxy*xxy +
                m.xxxz*xxz + m.xxyy*xyy + m.xxyz*xyz +
                m.xyyz*yyz ),
    ty = g4 * ( m.xyyy*xyy + m.xxxy*xxx + m.yyyy*yyy +
                m.yyyz*yyz + m.xxyy*xxy + m.xxyz*xxz +
                m.xyyz*xyz ),
    tz = g4 * (- m.xxxx*xxz - (m.xyyy + m.xxxy)*xyz
               - m.yyyy*yyz + m.xxxz*xxx + m.yyyz*yyy
               - m.xxyy*(xxz + yyz) + m.xxyz*xxy + m.xyyz*xyy);

  g4 = 0.25*(tx*x + ty*y + tz*z);

  /* Note that these variables have already been initialized; we're re-using them. */
  xxx = g3 * (m.xxx*xx + m.xyy*yy + m.xxy*xy + m.xxz*xz + m.xyz*yz);
  xxy = g3 * (m.xyy*xy + m.xxy*xx + m.yyy*yy + m.yyz*yz + m.xyz*xz);
  xxz = g3 * (-(m.xxx + m.xyy)*xz - (m.xxy + m.yyy)*yz + m.xxz*xx + m.yyz*yy + m.xyz*xy);

  g3 = onethird * (xxx*x + xxy*y + xxz*z);

  xx = g2*(m.xx*x + m.xy*y + m.xz*z);
  xy = g2*(m.yy*y + m.xy*x + m.yz*z);
  xz = g2*(-(m.xx + m.yy)*z + m.xz*x + m.yz*y);

  g2 = 0.5f*(xx*x + xy*y + xz*z);
  g0 *= m.totalMass;

  /* store the calculated potential  */
  *pot += -(g0 + g2 + g3 + g4);

  g0 += 5*g2 + 7*g3 + 9*g4;
  /* and the calculated acceleration. */
  out->x += dir*(xx + xxx + tx - x*g0);
  out->y += dir*(xy + xxy + ty - y*g0);
  out->z += dir*(xz + xxz + tz - z*g0);
}

#ifdef CAMBRIDGE

#define OPENING_GEOMETRY_FACTOR (2.0 / sqrt(3.0))

__device__ inline bool __attribute__(( always_inline ))
cuda_intersect(CudaMultipoleMoments &b, CudaSphere &s) {
  cudatype dsq = 0.0;
  cudatype rsq = s.radius * s.radius;
  cudatype delta;
  if((delta = b.lesser_corner.x - s.origin.x) > 0)
      dsq += delta * delta;
    else if((delta = s.origin.x - b.greater_corner.x) > 0)
      dsq += delta * delta;
    if(rsq < dsq)
      return false;
    if((delta = b.lesser_corner.y - s.origin.y) > 0)
      dsq += delta * delta;
    else if((delta = s.origin.y - b.greater_corner.y) > 0)
      dsq += delta * delta;
    if(rsq < dsq)
      return false;
    if((delta = b.lesser_corner.z - s.origin.z) > 0)
      dsq += delta * delta;
    else if((delta = s.origin.z - b.greater_corner.z) > 0)
      dsq += delta * delta;
    return (dsq <= s.radius * s.radius);
}

__device__ inline bool __attribute__(( always_inline ))
cuda_intersect(CudaSphere &s1, CudaSphere &s2) {
  CudaVector3D diff;
  cudatype dist;
  minusCudaVector3D(s1.origin, s2.origin, diff);
  dist = sqrt(diff.x*diff.x + diff.y*diff.y + diff.z*diff.z);
  return (dist <= s1.radius + s2.radius);
}


__device__ inline bool __attribute__(( always_inline ))
cuda_contains(const CudaSphere &s, const CudaVector3D &v) {
  CudaVector3D diff;
  cudatype dist;
  minusCudaVector3D(s.origin, v, diff);
  dist = sqrt(diff.x*diff.x + diff.y*diff.y + diff.z*diff.z);
  return (dist <= s.radius);
}

__device__ inline bool __attribute__(( always_inline ))
cuda_contained(const CudaMultipoleMoments &b, const CudaSphere &s) {
  CudaVector3D vec;
  cudatype delta1;
  cudatype delta2;

  delta1 = b.lesser_corner.x - s.origin.x;
  delta2 = s.origin.x - b.greater_corner.x;
  if(delta1 > 0)
    vec.x = b.greater_corner.x;
  else if(delta2 > 0)
    vec.x = b.lesser_corner.x;
  else{
    delta1 = -delta1;
    delta2 = -delta2;
    if(delta1>=delta2)
      vec.x = b.lesser_corner.x;
    else
      vec.x = b.greater_corner.x;
  }
  
  delta1 = b.lesser_corner.y - s.origin.y;
  delta2 = s.origin.y - b.greater_corner.y;
  if(delta1 > 0)
    vec.y = b.greater_corner.y;
  else if(delta2 > 0)
    vec.y = b.lesser_corner.y;
  else{
    delta1 = -delta1;
    delta2 = -delta2;
    if(delta1>=delta2)
      vec.y = b.lesser_corner.y;
    else
      vec.y = b.greater_corner.y;
  }

  delta1 = b.lesser_corner.z - s.origin.z;
  delta2 = s.origin.z - b.greater_corner.z;
  if(delta1 > 0)
    vec.z = b.greater_corner.z;
  else if(delta2 > 0)
    vec.z = b.lesser_corner.z;
  else{
    delta1 = -delta1;
    delta2 = -delta2;
    if(delta1>=delta2)
      vec.z = b.lesser_corner.z;
    else
      vec.z = b.greater_corner.z;
  }

  if(cuda_contains(s,vec))
    return true;
  else
    return false;

}

__device__ inline int __attribute__(( always_inline ))
cuda_openSoftening(CudaMultipoleMoments &node, CudaMultipoleMoments &myNode,
      CudaVector3D &offset) {
  CudaSphere s;
  addCudaVector3D(node.cm, offset, s.origin);
  s.radius = 2.0*node.soft;

  CudaSphere myS;
  assignCudaVector3D(node.cm, myS.origin);
  myS.radius = 2.0*node.soft;

  if(cuda_intersect(myS, s))
    return true;
  return cuda_intersect(myNode, s);
}

__device__ inline int __attribute__(( always_inline ))
cuda_encodeOffset(int reqID, int x, int y, int z) {
  // Replica in each direction is mapped to 0-7 range
  int offsetcode = (x + 3) | ((y+3) << 3) | ((z+3) << 6);

  // Assume we only have 32 bits to work with (minus sign bit)
  // 4 million buckets is our limitation
  return reqID | (offsetcode << 22);
}

__device__ inline int __attribute__(( always_inline ))
cuda_reEncodeOffset(int reqID, int offsetID)
{
    const int offsetmask = 0x1ff << 22;

    return reqID | (offsetmask & offsetID);
}

__device__ inline CudaVector3D __attribute__(( always_inline ))
cuda_decodeOffset(int reqID, cudatype fPeriod) {
  int offsetcode = reqID >> 22;
  int x = (offsetcode & 0x7) - 3;
  int y = ((offsetcode >> 3) & 0x7) - 3;
  int z = ((offsetcode >> 6) & 0x7) - 3;

  CudaVector3D offset;
  offset.x = x*fPeriod;
  offset.y = y*fPeriod;
  offset.z = z*fPeriod;

  return offset;
}

__device__ inline int __attribute__(( always_inline ))
cuda_openCriterionNode(CudaMultipoleMoments &node,
                    CudaMultipoleMoments &myNode,
                    CudaVector3D &offset,
                    int localIndex,
                    cudatype theta,
                    cudatype thetaMono) {
  const int nMinParticleNode = 6;
//  if(node.bucketSize <= nMinParticleNode) {
  if(node.particleCount <= nMinParticleNode) {
    return 1;
  }

  // Note that some of this could be pre-calculated into an "opening radius"
  cudatype radius = OPENING_GEOMETRY_FACTOR * node.radius / theta;

  if(radius < node.radius)
      radius = node.radius;

  CudaSphere s;
  addCudaVector3D(node.cm, offset, s.origin);
  s.radius = radius;

  if(myNode.type==cuda_Bucket || myNode.type==cuda_CachedBucket || myNode.type==cuda_NonLocalBucket){
    if(cuda_intersect(myNode, s))
        return 1;
    else
#ifdef HEXADECAPOLE
    {
      // Well separated, now check softening
      if(!cuda_openSoftening(node, myNode, offset)) {
        return 0;   // passed both tests: will be a Hex interaction
      } else {      // Open as monopole?
        radius = OPENING_GEOMETRY_FACTOR*node.radius/thetaMono;
        CudaSphere sM;
        addCudaVector3D(node.cm, offset, sM.origin);
        s.radius = radius;
        if(cuda_intersect(myNode, sM))
          return 1;
        else
           return 0;
      }
    }
#else
      return 0;
#endif
    } else {
      if(cuda_intersect(myNode, s)){
        if(cuda_contained(myNode,s))
          return 1;
        else
          return -1;
      }
      else
#ifdef HEXADECAPOLE
      {
      // Well separated, now check softening
      if(!cuda_openSoftening(node, myNode, offset)) {
        return 0;   // passed both tests: will be a Hex interaction
      }
      else {      // Open as monopole?
        radius = OPENING_GEOMETRY_FACTOR*node.radius/thetaMono;
        CudaSphere sM;
        addCudaVector3D(node.cm, offset, sM.origin);
        s.radius = radius;
        if(cuda_intersect(myNode, sM))
          return 1;
        else
          return 0;
        }
      }
#else
      return 0;
#endif
    }
}

__device__ inline void __attribute__(( always_inline ))
cuda_SPLINEQ(cudatype invr, cudatype r2, cudatype twoh, cudatype& a,
       cudatype& b,cudatype& c,cudatype& d) {
  cudatype u,dih,dir=(invr);
  if ((r2) < (twoh)*(twoh)) {
    dih = cudatype(2.0)/(twoh);
    u = dih/dir;
    if (u < cudatype(1.0)) {
      a = dih*(cudatype(7.0)/cudatype(5.0) 
         - cudatype(2.0)/cudatype(3.0)*u*u 
         + cudatype(3.0)/cudatype(10.0)*u*u*u*u
         - cudatype(1.0)/cudatype(10.0)*u*u*u*u*u);
      b = dih*dih*dih*(cudatype(4.0)/cudatype(3.0) 
           - cudatype(6.0)/cudatype(5.0)*u*u 
           + cudatype(1.0)/cudatype(2.0)*u*u*u);
      c = dih*dih*dih*dih*dih*(cudatype(12.0)/cudatype(5.0) 
             - cudatype(3.0)/cudatype(2.0)*u);
      d = cudatype(3.0)/cudatype(2.0)*dih*dih*dih*dih*dih*dih*dir;
    }
    else {
      a = cudatype(-1.0)/cudatype(15.0)*dir 
  + dih*(cudatype(8.0)/cudatype(5.0) 
         - cudatype(4.0)/cudatype(3.0)*u*u + u*u*u
         - cudatype(3.0)/cudatype(10.0)*u*u*u*u 
         + cudatype(1.0)/cudatype(30.0)*u*u*u*u*u);
      b = cudatype(-1.0)/cudatype(15.0)*dir*dir*dir 
  + dih*dih*dih*(cudatype(8.0)/cudatype(3.0) - cudatype(3.0)*u 
           + cudatype(6.0)/cudatype(5.0)*u*u 
           - cudatype(1.0)/cudatype(6.0)*u*u*u);
      c = cudatype(-1.0)/cudatype(5.0)*dir*dir*dir*dir*dir 
  + cudatype(3.0)*dih*dih*dih*dih*dir
  + dih*dih*dih*dih*dih*(cudatype(-12.0)/cudatype(5.0) 
             + cudatype(1.0)/cudatype(2.0)*u);
      d = -dir*dir*dir*dir*dir*dir*dir
  + cudatype(3.0)*dih*dih*dih*dih*dir*dir*dir
  - cudatype(1.0)/cudatype(2.0)*dih*dih*dih*dih*dih*dih*dir;
    }
  }
  else {
    a = dir;
    b = a*a*a;
    c = cudatype(3.0)*b*a*a;
    d = cudatype(5.0)*c*a*a;
  }
}

__device__ inline void __attribute__(( always_inline ))
cuda_SPLINE(cudatype r2, cudatype twoh, cudatype &a, cudatype &b) {
  cudatype r, u,dih,dir;
  r = sqrt(r2);

  if (r < (twoh)) {
    dih = (2.0)/(twoh);
    u = r*dih;
    if (u < (1.0)) {
      a = dih*((7.0)/(5.0) 
         - (2.0)/(3.0)*u*u 
         + (3.0)/(10.0)*u*u*u*u
         - (1.0)/(10.0)*u*u*u*u*u);
      b = dih*dih*dih*((4.0)/(3.0) 
           - (6.0)/(5.0)*u*u 
           + (1.0)/(2.0)*u*u*u);
    }
    else {
      dir = (1.0)/r;
      a = (-1.0)/(15.0)*dir 
  + dih*((8.0)/(5.0) 
         - (4.0)/(3.0)*u*u + u*u*u
         - (3.0)/(10.0)*u*u*u*u 
         + (1.0)/(30.0)*u*u*u*u*u);
      b = (-1.0)/(15.0)*dir*dir*dir 
  + dih*dih*dih*((8.0)/(3.0) - (3.0)*u 
           + (6.0)/(5.0)*u*u 
           - (1.0)/(6.0)*u*u*u);
    }
  }
  else {
    a = (1.0)/r;
    b = a*a*a;
  }
}

__device__ inline int __attribute__(( always_inline ))
cuda_OptAction(int fakeOpen, int nodetype) {
  if (fakeOpen == 0) {
    if (nodetype == cuda_Internal || nodetype == cuda_Bucket || nodetype == cuda_Boundary || nodetype == cuda_NonLocalBucket) {
      return COMPUTE;
    } else if (nodetype == cuda_NonLocal || nodetype == cuda_Cached || nodetype == cuda_CachedBucket || nodetype == cuda_Empty || nodetype == cuda_CachedEmpty) {
      return DUMP;
    } else if (nodetype == cuda_Top || nodetype == cuda_Invalid) {
      return ERROR;
    } else {
      printf("ERROR in cuda_OptAction\n");
      return -1;
    }
  } else {
    if (nodetype == cuda_Internal || nodetype == cuda_Boundary) {
      return KEEP;
    } else if (nodetype == cuda_Bucket) {
      return KEEP_LOCAL_BUCKET;
    } else if (nodetype == cuda_NonLocal || nodetype == cuda_NonLocalBucket || nodetype == cuda_CachedBucket || nodetype == cuda_Cached || nodetype == cuda_Empty ||
              nodetype == cuda_CachedEmpty) {
      return DUMP;
    } else if (nodetype == cuda_Top || nodetype == cuda_Invalid) {
      return ERROR;
    } else {
      printf("ERROR in cuda_OptAction\n");
      return -1;
    }
  }
}

#endif


