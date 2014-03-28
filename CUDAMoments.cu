#ifdef CUDA_UNIT_TEST
# include "CUDAMoments.h"
#else
# include "cuda_typedef.h"
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
  const cudatype onethird = 1.0f / 3.0f;

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
