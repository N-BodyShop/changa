#ifndef __GRAVITY_H__
#define __GRAVITY_H__

#include "TreeNode.h"
#include "GenericTreeNode.h"
#include "Space.h"
#include "SSEdefs.h"
#include <vector>

extern cosmoType theta;
extern cosmoType thetaMono;

/// A queued node contribution for a staged bucket. Moments are copied by value
/// so the queue does not depend on the tree node's lifetime.
struct BucketNodeTerm {
  MultipoleMoments m;
  Vector3D<cosmoType> cm;
};

/// A source particle queued for the same batched kernel as nodes.
/// base is offset + source position, precomputed in the scalar order.
struct BucketPartTerm {
  Vector3D<cosmoType> base;
  cosmoType mass;
  cosmoType soft;
};

/// Queued terms per flush. Scales with the compile-time SIMD width.
/// The fork has no explicit batching. CHANGA_TERM_BATCH_MULT may be
/// overridden at build time to tune the batch size. The default of 8 gives
/// 16 terms per flush for SSE double and 32 for SSE float/AVX double.
#ifndef CHANGA_TERM_BATCH_MULT
#define CHANGA_TERM_BATCH_MULT 8
#endif
enum { TERM_BATCH_MAX = CHANGA_TERM_BATCH_MULT * SSE_VECTOR_WIDTH };

/// FORCE_SOA_FIELDS is the number of per-particle cosmoType arrays in
/// ActiveSoA (px through az), so one bind() buffer holds
/// FORCE_SOA_FIELDS * nPad values.
enum { FORCE_SOA_FIELDS = 10 };

/// Padded SoA of one bucket's active particles, with FORCE_SOA_FIELDS arrays
/// of nPad lanes each. Lanes from n to nPad are zeroed so whole SIMD vectors
/// can be loaded; the kernels' masks keep those lanes from contributing, and
/// scatter() never writes them back. Owning storage lives in buf; the
/// no-bucket fallbacks bind the pointers to stack memory via bind().
/// @var n     Number of active particles.
/// @var nPad  n rounded up to a multiple of SSE_VECTOR_WIDTH.
/// @var buf   Owned storage (empty when bound to caller memory).
/// @var px..az  Pointers into buf (or bound memory): position, softening,
///            mass, and dtGrav inputs; potential and acceleration outputs.
struct ActiveSoA {
  int n;
  int nPad;
  std::vector<cosmoType> buf;
  cosmoType *px, *py, *pz, *soft, *mass, *dtGrav, *potential, *ax, *ay, *az;

  /// Default constructor, initializes n and nPad to zero.
  ActiveSoA() : n(0), nPad(0) {}

  /// Rounds count up to a whole number of SIMD vectors.
  /// @param n     Particle count.
  /// @return      n rounded up to a multiple of SSE_VECTOR_WIDTH.
  static int padded(int n) {
    return ((n + SSE_VECTOR_WIDTH - 1) / SSE_VECTOR_WIDTH) * SSE_VECTOR_WIDTH;
  }

  /// Points the fields at base, which must hold FORCE_SOA_FIELDS *
  /// padded(count) values.
  void bind(cosmoType *base, int count) {
    n = count;
    nPad = padded(count);
    px = base;
    py = px + nPad;
    pz = py + nPad;
    soft = pz + nPad;
    mass = soft + nPad;
    dtGrav = mass + nPad;
    potential = dtGrav + nPad;
    ax = potential + nPad;
    ay = ax + nPad;
    az = ay + nPad;
  }

  /// Owns storage for count particles.
  void allocate(int count) {
    buf.assign(FORCE_SOA_FIELDS * padded(count), 0.0);
    bind(buf.data(), count);
  }

  /// Clears buffer and resets n and nPad to zero after release
  void release() {
    std::vector<cosmoType>().swap(buf);
    n = nPad = 0;
  }

  /// Copies the AoS particles act[0..n) into the SoA in one pass over the
  /// padded lanes, zeroing lanes from n to nPad.
  void gather(GravityParticle **act) {
    for (int i = 0; i < nPad; ++i) {
      const GravityParticle *p = i < n ? act[i] : nullptr;
      px[i] = p ? p->position.x : 0.0;
      py[i] = p ? p->position.y : 0.0;
      pz[i] = p ? p->position.z : 0.0;
      soft[i] = p ? p->soft : 0.0;
      mass[i] = p ? p->mass : 0.0;
      dtGrav[i] = p ? p->dtGrav : 0.0;
      potential[i] = p ? p->potential : 0.0;
      ax[i] = p ? p->treeAcceleration.x : 0.0;
      ay[i] = p ? p->treeAcceleration.y : 0.0;
      az[i] = p ? p->treeAcceleration.z : 0.0;
    }
  }

  /// Writes acceleration, potential and dtGrav back to the AoS.
  void scatter(GravityParticle **act) const {
    for (int i = 0; i < n; ++i) {
      GravityParticle *p = act[i];
      p->treeAcceleration.x = ax[i];
      p->treeAcceleration.y = ay[i];
      p->treeAcceleration.z = az[i];
      p->potential = potential[i];
      p->dtGrav = dtGrav[i];
    }
  }
};

/// Per-bucket staging of active particles. One is owned by each TreePiece
/// bucket. It is gathered once when the bucket is staged and scattered back
/// when it is suspended or finished (finishBucket).
/// Node and particle terms share one walk-ordered queue: order[k] >= 0
/// indexes nodes, order[k] < 0 indexes parts as ~order[k]. flushTerms
/// applies them in that order, so each particle's accumulation matches the
/// scalar walk.
/// @var staged  True while soa holds gathered particles for this bucket.
/// @var done    True once unstageBucket has run for this walk; the bucket is
///              then never restaged and the AoS path is used. Cleared at the
///              start of each walk (TreePiece.cpp).
/// @var activeRung  Lowest rung that counts as active for this bucket.
/// @var soa     Padded SoA of the active particles.
/// @var nodes   Queued node terms (copied moments, walk order).
/// @var parts   Queued particle terms (offset + position precomputed).
/// @var order   Walk-order index into nodes (>= 0) or parts (~value < 0).
/// @var act     Active particles in gather order; scatter target.
struct BucketSoA {
  bool staged;
  bool done;
  int activeRung;
  ActiveSoA soa;
  std::vector<BucketNodeTerm> nodes;
  std::vector<BucketPartTerm> parts;
  std::vector<int> order;
  std::vector<GravityParticle*> act;
  BucketSoA() : staged(false), done(false), activeRung(0) {}
};

/*
** see (A1) and (A2) of TREESPH: A UNIFICATION OF SPH WITH THE 
** HIERARCHICAL TREE METHOD by Lars Hernquist and Neal Katz.
** APJ Supplemant Series 70:416-446, 1989
** 
** Higher derivative terms c and d for use with quadrupole spline
** softening (Joachim Stadel, Dec. 94).
*/
// Lane-wise m ? t : f. A pure copy, so every selected value is bit-exact.
// In-line where masks used for regime selection in SPLINEQ function.
inline SSEcosmoType sseSelect(SSEMask m, SSEcosmoType t, SSEcosmoType f) {
  stdx::where(m, f) = t;
  return f;
}
// m ? x : +0.0 per lane (zeroes masked-off lanes, which may hold inf/NaN).
inline SSEcosmoType sseAnd(SSEMask m, SSEcosmoType x) {
  return sseSelect(m, x, SSEcosmoType(COSMO_CONST(0.0)));
}
// Same lane semantics as Charm's _mm_max_pd(a, b): a > b ? a : b.
inline SSEcosmoType sseMax(SSEcosmoType a, SSEcosmoType b) {
  return sseSelect(a > b, a, b);
}

inline
void SPLINEQ(SSEcosmoType invr, SSEcosmoType r2, SSEcosmoType twoh,
	     SSEcosmoType& a, SSEcosmoType& b, 
	     SSEcosmoType& c, SSEcosmoType& d)
{
  SSEcosmoType dir = invr;
  SSEMask select0 = r2 < twoh * twoh;
  // Common case: no lane inside the softening radius. Same arithmetic as the
  // a0..d0 regime below, so results are bit-identical.
  if (stdx::none_of(select0)) {
    a = dir;
    b = dir*dir*dir;
    c = COSMO_CONST(3.0)*b*dir*dir;
    d = COSMO_CONST(5.0)*c*dir*dir;
    return;
  }
  // Every regime is evaluated on all lanes and then selected per lane.
  // Lanes outside a regime may hold inf/NaN; the selects discard them.
  SSEcosmoType dih = COSMO_CONST(2.0)/twoh;
  SSEcosmoType u = dih/dir;
  SSEMask select1 = u < COSMO_CONST(1.0);

  // Operation order matches the fork's SSE SPLINEQ (direct powers, left to
  // right) so results are bit-identical to the fork.
  SSEcosmoType a0 = dir;
  SSEcosmoType b0 = a0*a0*a0;
  SSEcosmoType c0 = COSMO_CONST(3.0)*b0*a0*a0;
  SSEcosmoType d0 = COSMO_CONST(5.0)*c0*a0*a0;

  SSEcosmoType a1 = dih*(COSMO_CONST(7.0)/COSMO_CONST(5.0)
			 - COSMO_CONST(2.0)/COSMO_CONST(3.0)*u*u
			 + COSMO_CONST(3.0)/COSMO_CONST(10.0)*u*u*u*u
			 - COSMO_CONST(1.0)/COSMO_CONST(10.0)*u*u*u*u*u);
  SSEcosmoType b1 = dih*dih*dih*(COSMO_CONST(4.0)/COSMO_CONST(3.0)
				 - COSMO_CONST(6.0)/COSMO_CONST(5.0)*u*u
				 + COSMO_CONST(1.0)/COSMO_CONST(2.0)*u*u*u);
  SSEcosmoType c1 = dih*dih*dih*dih*dih*(COSMO_CONST(12.0)/COSMO_CONST(5.0)
					 - COSMO_CONST(3.0)/COSMO_CONST(2.0)*u);
  SSEcosmoType d1 = COSMO_CONST(3.0)/COSMO_CONST(2.0)*dih*dih*dih*dih*dih
    *dih*dir;

  SSEcosmoType a2 = COSMO_CONST(-1.0)/COSMO_CONST(15.0)*dir
    + dih*(COSMO_CONST(8.0)/COSMO_CONST(5.0)
	   - COSMO_CONST(4.0)/COSMO_CONST(3.0)*u*u + u*u*u
	   - COSMO_CONST(3.0)/COSMO_CONST(10.0)*u*u*u*u
	   + COSMO_CONST(1.0)/COSMO_CONST(30.0)*u*u*u*u*u);
  SSEcosmoType b2 = COSMO_CONST(-1.0)/COSMO_CONST(15.0)*dir*dir*dir
    + dih*dih*dih*(COSMO_CONST(8.0)/COSMO_CONST(3.0) - COSMO_CONST(3.0)*u
		   + COSMO_CONST(6.0)/COSMO_CONST(5.0)*u*u
		   - COSMO_CONST(1.0)/COSMO_CONST(6.0)*u*u*u);
  SSEcosmoType c2 = COSMO_CONST(-1.0)/COSMO_CONST(5.0)*dir*dir*dir*dir*dir
    + COSMO_CONST(3.0)*dih*dih*dih*dih*dir
    + dih*dih*dih*dih*dih*(COSMO_CONST(-12.0)/COSMO_CONST(5.0)
			   + COSMO_CONST(1.0)/COSMO_CONST(2.0)*u);
  SSEcosmoType d2 = -dir*dir*dir*dir*dir*dir*dir
    + COSMO_CONST(3.0)*dih*dih*dih*dih*dir*dir*dir
    - COSMO_CONST(1.0)/COSMO_CONST(2.0)*dih*dih*dih*dih*dih*dih*dir;

  a = sseSelect(select0, sseSelect(select1, a1, a2), a0);
  b = sseSelect(select0, sseSelect(select1, b1, b2), b0);
  c = sseSelect(select0, sseSelect(select1, c1, c2), c0);
  d = sseSelect(select0, sseSelect(select1, d1, d2), d0);
}

/// Spline softening kernel, single SIMD register. Computes the potential
/// factor a and force factor b for each lane. With u = 2r / twoh, the
/// regimes are r >= twoh (plain 1/r), 1 <= u < 2 (outer spline), and
/// u < 1 (inner spline). Every regime is evaluated and selected per lane.
/// @param r2    Squared distance per lane.
/// @param twoh  Softening length (2h) per lane.
/// @param a     [out] Potential factor per lane.
/// @param b     [out] Force factor per lane (a*a*a in the plain regime).
inline
void SPLINE(SSEcosmoType r2, SSEcosmoType twoh,
	    SSEcosmoType &a, SSEcosmoType &b)
{
  SSEcosmoType r = sqrt(r2);
  SSEMask select0 = r < twoh;
  if (stdx::none_of(select0)) {
    a = COSMO_CONST(1.0)/r;
    b = a*a*a;
    return;
  }
  // Every regime is evaluated on all lanes and then selected per lane.
  // Lanes outside a regime may hold inf/NaN; the selects discard them.
  SSEcosmoType dih = COSMO_CONST(2.0)/twoh;
  SSEcosmoType u = r * dih;
  SSEcosmoType dir = COSMO_CONST(1.0)/r;
  SSEMask select1 = u < COSMO_CONST(1.0);

  // Operation order matches the fork's SSE SPLINE (direct powers, left to
  // right) so results are bit-identical to the fork.
  SSEcosmoType a0 = dir;
  SSEcosmoType b0 = a0 * a0 * a0;

  SSEcosmoType a1 = dih*(COSMO_CONST(7.0)/COSMO_CONST(5.0)
			 - COSMO_CONST(2.0)/COSMO_CONST(3.0)*u*u
			 + COSMO_CONST(3.0)/COSMO_CONST(10.0)*u*u*u*u
			 - COSMO_CONST(1.0)/COSMO_CONST(10.0)*u*u*u*u*u);
  SSEcosmoType b1 = dih*dih*dih*(COSMO_CONST(4.0)/COSMO_CONST(3.0)
				 - COSMO_CONST(6.0)/COSMO_CONST(5.0)*u*u
				 + COSMO_CONST(1.0)/COSMO_CONST(2.0)*u*u*u);

  SSEcosmoType a2 = COSMO_CONST(-1.0)/COSMO_CONST(15.0)*dir
    + dih*(COSMO_CONST(8.0)/COSMO_CONST(5.0)
	   - COSMO_CONST(4.0)/COSMO_CONST(3.0)*u*u + u*u*u
	   - COSMO_CONST(3.0)/COSMO_CONST(10.0)*u*u*u*u
	   + COSMO_CONST(1.0)/COSMO_CONST(30.0)*u*u*u*u*u);
  SSEcosmoType b2 = COSMO_CONST(-1.0)/COSMO_CONST(15.0)*dir*dir*dir
    + dih*dih*dih*(COSMO_CONST(8.0)/COSMO_CONST(3.0) - COSMO_CONST(3.0)*u
		   + COSMO_CONST(6.0)/COSMO_CONST(5.0)*u*u
		   - COSMO_CONST(1.0)/COSMO_CONST(6.0)*u*u*u);

  a = sseSelect(select0, sseSelect(select1, a1, a2), a0);
  b = sseSelect(select0, sseSelect(select1, b1, b2), b0);
}

//
// Return true if the soften nodes overlap, or if the source node's
// softening overlaps the bounding box; i.e. the forces involve softening
// @param node Node to test
// @param myNode
// @param offset Periodic offset applied to node
//
inline
int openSoftening(Tree::GenericTreeNode *node, Tree::GenericTreeNode *myNode,
		  Vector3D<cosmoType> offset)
{
  Sphere<cosmoType> s(node->moments.cm + offset, 2.0*node->moments.soft);
  Sphere<cosmoType> myS(myNode->moments.cm, 2.0*myNode->moments.soft);
  if(Space::intersect(myS, s))
      return true;
  return Space::intersect(myNode->boundingBox, s);
}

#ifdef CMK_VERSION_BLUEGENE
static int forProgress = 0;
#endif

// Masked SIMD kernel over padded SoA lanes. Terms are applied in queue
// order for every register, so each lane accumulates exactly as the scalar
// walk does. Node terms follow nodeBucketForce; particle terms follow
// partBucketForce.
/// @param nodes  Queued node terms, indexed by order[k] >= 0 (may be null
///               when nt has no nodes).
/// @param parts  Queued particle terms, indexed by ~order[k] < 0.
/// @param order  Walk-order list of term indices, length nt.
/// @param nt     Number of queued terms to apply.
/// @param s      Padded SoA whose force fields are accumulated in place.
/// @param nPad   Lane count of s (a multiple of SSE_VECTOR_WIDTH).
inline void termBatchSoA(const BucketNodeTerm *nodes,
			 const BucketPartTerm *parts,
			 const int *order, int nt,
			 const ActiveSoA &s, int nPad) {
  Vector3D<SSEcosmoType> r;
  SSEcosmoType rsq;
  SSEcosmoType twoh;
  SSEcosmoType a,b,c,d;

  for (int i=0; i<nPad; i+=SSE_VECTOR_WIDTH) {
#ifdef CMK_VERSION_BLUEGENE
    if (++forProgress > 200) {
      forProgress = 0;
#ifdef COSMO_EVENTS
      traceUserEvents(networkProgressUE);
#endif
      CmiNetworkProgress();
    }
#endif
    Vector3D<SSEcosmoType> 
      packedPos(SSELoad(SSEcosmoType, s.px, i, ),
		SSELoad(SSEcosmoType, s.py, i, ),
		SSELoad(SSEcosmoType, s.pz, i, )); 
    Vector3D<SSEcosmoType> 
      packedAcc(SSELoad(SSEcosmoType, s.ax, i, ),
		SSELoad(SSEcosmoType, s.ay, i, ),
		SSELoad(SSEcosmoType, s.az, i, )); 
    SSEcosmoType SSELoad(packedPotential, s.potential, i, ); 
    SSEcosmoType SSELoad(packedMass, s.mass, i, );  
    SSEcosmoType SSELoad(packedDtGrav, s.dtGrav, i, );
    SSEcosmoType SSELoad(packedSoft, s.soft, i, );
    for (int k=0; k<nt; ++k) {
    if (order[k] < 0) {
      const BucketPartTerm &p = parts[~order[k]];
      Vector3D<SSEcosmoType> baseV(SSEcosmoType(p.base.x),
				   SSEcosmoType(p.base.y),
				   SSEcosmoType(p.base.z));
      // Operand order matches the scalar path: (offset + source) - packedPos.
      // The fork's SSE path computes -packedPos + offset + part
      r = baseV - packedPos;
      rsq = r.lengthSquared();
      twoh = SSEcosmoType(p.soft) + packedSoft;
      SSEMask select = rsq > COSMO_CONST(0.0);
      // Matches partBucketForce: no active lane, so skip the spline and update.
      if (stdx::none_of(select)) continue;
      SPLINE(rsq, twoh, a, b);
      a = sseAnd(select, a);
      b = sseAnd(select, b);
      SSEcosmoType idt2 = (packedMass + SSEcosmoType(p.mass)) * b;
      idt2 = sseMax(idt2, packedDtGrav);
      packedAcc += r * (b * SSEcosmoType(p.mass));
      packedPotential -= SSEcosmoType(p.mass) * a;
      packedDtGrav = idt2;
      continue;
    }
    const MultipoleMoments &m = nodes[order[k]].m;
    const Vector3D<cosmoType> &cm = nodes[order[k]].cm;
    r = packedPos - cm;
    rsq = r.lengthSquared();
    SSEMask select = rsq > COSMO_CONST(0.0);
    // Self and padded lanes only: nothing to accumulate.
    if (stdx::none_of(select)) continue;
    /* Masked-off lanes: 1/sqrt(0) is inf for self and padded lanes. */
    SSEcosmoType dir = sseAnd(select, COSMO_CONST(1.0)/sqrt(rsq));
#ifdef HEXADECAPOLE
    SSEcosmoType magai;
    momEvalFmomrcmVec(&m.mom, m.getRadius(), dir, r,
		      &packedPotential, &packedAcc, &magai);
    SSEcosmoType idt2 = (packedMass + m.totalMass)*dir*dir*dir;
#else
    twoh = CONVERT_TO_COSMO_TYPE m.soft + packedSoft;
    SPLINEQ(dir, rsq, twoh, a, b, c, d);
    a = sseAnd(select, a);
    b = sseAnd(select, b);
    c = sseAnd(select, c);
    d = sseAnd(select, d);
    Vector3D<SSEcosmoType> rowX(SSEcosmoType(CONVERT_TO_COSMO_TYPE m.xx),
				SSEcosmoType(CONVERT_TO_COSMO_TYPE m.xy),
				SSEcosmoType(CONVERT_TO_COSMO_TYPE m.xz));
    Vector3D<SSEcosmoType> rowY(SSEcosmoType(CONVERT_TO_COSMO_TYPE m.xy),
				SSEcosmoType(CONVERT_TO_COSMO_TYPE m.yy),
				SSEcosmoType(CONVERT_TO_COSMO_TYPE m.yz));
    Vector3D<SSEcosmoType> rowZ(SSEcosmoType(CONVERT_TO_COSMO_TYPE m.xz),
				SSEcosmoType(CONVERT_TO_COSMO_TYPE m.yz),
				SSEcosmoType(CONVERT_TO_COSMO_TYPE m.zz));
    Vector3D<SSEcosmoType> qi(vdot3(rowX, r), vdot3(rowY, r), vdot3(rowZ, r));
    SSEcosmoType qir = COSMO_CONST(0.5)*vdot3(qi, r);
    SSEcosmoType tr = COSMO_CONST(0.5)*(CONVERT_TO_COSMO_TYPE m.xx 
					+ CONVERT_TO_COSMO_TYPE m.yy 
					+ CONVERT_TO_COSMO_TYPE m.zz);
    SSEcosmoType qir3 = b*CONVERT_TO_COSMO_TYPE m.totalMass + d*qir - c*tr;
    packedPotential -= CONVERT_TO_COSMO_TYPE m.totalMass *a 
      + c * qir - b * tr;      
    packedAcc -= r*qir3 - qi*c;
    SSEcosmoType idt2 = (packedMass + CONVERT_TO_COSMO_TYPE m.totalMass)*b;
#endif
    idt2 = sseMax(idt2, packedDtGrav);
    packedDtGrav = idt2;
    }
    SSEStore(packedAcc.x, s.ax, i, );
    SSEStore(packedAcc.y, s.ay, i, );
    SSEStore(packedAcc.z, s.az, i, );
    SSEStore(packedPotential, s.potential, i, );
    SSEStore(packedDtGrav, s.dtGrav, i, ); 
  }
}

// Masked SIMD partBucketForce kernel over padded SoA lanes.
// Runs one particle source through termBatchSoA as a single-term batch.
inline void partKernelSoA(ExternalGravityParticle *part, const ActiveSoA &s,
			  int nPad, Vector3D<cosmoType> offset) {
  BucketPartTerm t = {offset + part->position, part->mass, part->soft};
  int idx = -1;
  termBatchSoA(nullptr, &t, &idx, 1, s, nPad);
}

inline int partBucketForce(ExternalGravityParticle *part, 
			   Tree::GenericTreeNode *req, 
			   GravityParticle **activeParticles, 
			   Vector3D<cosmoType> offset, int nActiveParts) {
  ActiveSoA s;
  s.bind((cosmoType*)alloca(FORCE_SOA_FIELDS * ActiveSoA::padded(nActiveParts)
			    * sizeof(cosmoType)), nActiveParts);
  s.gather(activeParticles);
  partKernelSoA(part, s, s.nPad, offset);
  s.scatter(activeParticles);
  return nActiveParts;
}

inline int partBucketForce(ExternalGravityParticle *part, 
			   Tree::GenericTreeNode *req, 
			   GravityParticle *particles, 
			   Vector3D<cosmoType> offset, int activeRung) {
  int nActiveParts = 0; 
  GravityParticle **activeParticles = 
    (GravityParticle**)alloca((req->lastParticle - req->firstParticle 
			       + 1 + FORCE_INPUT_LIST_PAD) 
			      * sizeof(GravityParticle*));
  for (int j=req->firstParticle; j <= req->lastParticle; ++j) {
    if (particles[j].rung >= activeRung) 
      activeParticles[nActiveParts++] = &particles[j];
  }
  
  int ret = partBucketForce(part, req, activeParticles, offset, nActiveParts); 
  return ret; 
}

//
// Calculated forces on active particles in a bucket due to the
// multipole of a TreeNode.  Return number of multipoles evaluated.
//

// Runs one node's moments through termBatchSoA as a single-term batch. Used by
// the unstaged nodeBucketForce path; staged buckets queue terms instead.
inline void nodeKernelSoA(const MultipoleMoments &m, Vector3D<cosmoType> cm,
			  const ActiveSoA &s, int nPad) {
  BucketNodeTerm t = {m, cm};
  int idx = 0;
  termBatchSoA(&t, nullptr, &idx, 1, s, nPad);
}

/// Applies all queued node and particle terms to the bucket's SoA in walk
/// order.
inline void flushTerms(BucketSoA &bs) {
  if (bs.order.empty())
    return;
  termBatchSoA(bs.nodes.data(), bs.parts.data(), bs.order.data(),
	       (int)bs.order.size(), bs.soa, bs.soa.nPad);
  bs.nodes.clear();
  bs.parts.clear();
  bs.order.clear();
}

//
// Calculated forces on active particles in a bucket due to the
// multipole of a TreeNode.  Return number of multipoles evaluated.
//
inline
int nodeBucketForce(Tree::GenericTreeNode *node, 
		    Tree::GenericTreeNode *req,  
		    GravityParticle *particles, 
		    Vector3D<cosmoType> offset,    
		    int activeRung)
{
  const MultipoleMoments &m = node->moments;
  Vector3D<cosmoType> cm(m.cm + offset);
  int nActiveParts = 0;

  GravityParticle **activeParticles = 
    (GravityParticle**)alloca((req->lastParticle - req->firstParticle 
		       + 1 + FORCE_INPUT_LIST_PAD) 
		      * sizeof(GravityParticle*));
  for (int j=req->firstParticle; j <= req->lastParticle; ++j) {
    if (particles[j].rung >= activeRung) 
      activeParticles[nActiveParts++] = &particles[j];
  }

#ifdef HEXADECAPOLE
  if(openSoftening(node, req, offset)) {
    ExternalGravityParticle tmpPart;
    tmpPart.mass = m.totalMass;
    tmpPart.soft = m.soft;
    tmpPart.position = m.cm;
    int ret = partBucketForce(&tmpPart, req, activeParticles, 
			   offset, nActiveParts);
    return ret; 
    }
#endif
  ActiveSoA s;
  s.bind((cosmoType*)alloca(FORCE_SOA_FIELDS * ActiveSoA::padded(nActiveParts)
			    * sizeof(cosmoType)), nActiveParts);
  s.gather(activeParticles);
  nodeKernelSoA(m, cm, s, s.nPad);
  s.scatter(activeParticles);
  return nActiveParts;
}

/// @brief Gravity opening criterion for a bucket walk.
/// @param node Source node to be tested
/// @param bucketNode Target bucket
/// @param offset Offset of periodic replica
/// @return True if the node's opening radius intersects the
/// boundingBox of the bucket, i.e. the node needs to be opened.

inline bool
openCriterionBucket(Tree::GenericTreeNode *node,
                   Tree::GenericTreeNode *bucketNode,
                   Vector3D<cosmoType> offset // Offset of node
                   ) {

#if COSMO_STATS > 0
  node->used = true;
#endif
  // Always open node if this many particles or fewer.
  const int nMinParticleNode = 6;
  if(node->particleCount <= nMinParticleNode) {
      return true;
      }

  // Note that some of this could be pre-calculated into an "opening radius"
  cosmoType radius = TreeStuff::opening_geometry_factor
    * node->moments.getRadius() / theta;
  if(radius < node->moments.getRadius())
      radius = node->moments.getRadius();

  Sphere<cosmoType> s(node->moments.cm + offset, radius);
  
#ifdef HEXADECAPOLE
  if(!Space::intersect(bucketNode->boundingBox, s)) {
      // Well separated, now check softening
      if(!openSoftening(node, bucketNode, offset)) {
	  return false; // passed both tests: will be a Hex interaction
      }
      else {        // Open as monopole?
        radius = TreeStuff::opening_geometry_factor
	  * node->moments.getRadius() / thetaMono;
      Sphere<cosmoType> sM(node->moments.cm + offset, radius);
      return Space::intersect(bucketNode->boundingBox, sM);
      }
      }
  return true;
#else
  return Space::intersect(bucketNode->boundingBox, s);
#endif
}

/// @brief Gravity opening criterion for "double walk".
/// @param node Source node to be tested
/// @param myNode Target node
/// @param offset Offset for periodic replica
/// @return -1 if there is an intersection
///    0 if no intersection
///    +1 if completely contained
///
/// This is the criterion described in Stadel(2001): if no
/// intersection, then the node is valid for all particles in myNode.
/// If completely contained then it needs to be opened for all
/// particles in myNode. If a partial intersection, then node needs
/// to be checked by the children of myNode. If myNode is already a
/// bucket, then only returns 0 or 1 are possible.

inline int openCriterionNode(Tree::GenericTreeNode *node,
                    Tree::GenericTreeNode *myNode,
                    Vector3D<cosmoType> offset
                    ) {

#if COSMO_STATS > 0
  node->used = true;
#endif
  // Always open node if this many particles or fewer.
  const int nMinParticleNode = 6;
  if(node->particleCount <= nMinParticleNode) {
      return 1;
      }

  // Note that some of this could be pre-calculated into an "opening radius"
  cosmoType radius = TreeStuff::opening_geometry_factor
    * node->moments.getRadius() / theta;
  if(radius < node->moments.getRadius())
      radius = node->moments.getRadius();

  Sphere<cosmoType> s(node->moments.cm + offset, radius);

  if(myNode->getType()==Tree::Bucket
     || myNode->getType()==Tree::CachedBucket
     || myNode->getType()==Tree::NonLocalBucket){
    if(Space::intersect(myNode->boundingBox, s))
        return 1;
    else
#ifdef HEXADECAPOLE
        {
        // Well separated, now check softening
        if(!openSoftening(node, myNode, offset)) {
            return 0;   // passed both tests: will be a Hex interaction
            }
        else {      // Open as monopole?
          radius = TreeStuff::opening_geometry_factor
	    * node->moments.getRadius() / thetaMono;
            Sphere<cosmoType> sM(node->moments.cm + offset, radius);
            if(Space::intersect(myNode->boundingBox, sM))
                return 1;
            else
                return 0;
            }
        }
#else
    return 0;
#endif
    }
    else{
        if(Space::intersect(myNode->boundingBox, s)){
            if(Space::contained(myNode->boundingBox,s))
                return 1;
            else
                return -1;
        }
        else
#ifdef HEXADECAPOLE
            {
            // Well separated, now check softening
            if(!openSoftening(node, myNode, offset)) {
                return 0;   // passed both tests: will be a Hex interaction
                }
            else {      // Open as monopole?
                radius = TreeStuff::opening_geometry_factor
		  * node->moments.getRadius() / thetaMono;
                Sphere<cosmoType> sM(node->moments.cm + offset, radius);
                if(Space::intersect(myNode->boundingBox, sM))
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

// Bucket-level staging. stageBucket gathers the bucket's active particles
// into bs once; the bucket overloads below run on that SoA until
// unstageBucket scatters the results back to the AoS particles and marks the
// bucket done (see BucketSoA::done). Used when the walk for a bucket ends.
inline void stageBucket(BucketSoA &bs, Tree::GenericTreeNode *req,
			GravityParticle *particles, int activeRung) {
  bs.activeRung = activeRung;
  bs.act.clear();
  bs.nodes.clear();
  bs.parts.clear();
  bs.order.clear();
  bs.nodes.reserve(TERM_BATCH_MAX);
  bs.parts.reserve(TERM_BATCH_MAX);
  bs.order.reserve(TERM_BATCH_MAX);
  for (int j = req->firstParticle; j <= req->lastParticle; ++j) {
    if (particles[j].rung >= activeRung)
      bs.act.push_back(&particles[j]);
  }
  bs.soa.allocate((int)bs.act.size());
  bs.soa.gather(bs.act.data());
  bs.staged = true;
}

// Flushes queued terms and writes the SoA back to the AoS, releasing the
// SoA. Leaves the bucket restageable (unlike unstageBucket) so it can be
// suspended around Ewald, which adds to the AoS directly.
inline void suspendBucket(BucketSoA &bs) {
  if (bs.staged) {
    flushTerms(bs);
    bs.soa.scatter(bs.act.data());
  }
  bs.staged = false;
  bs.soa.release();
  std::vector<GravityParticle*>().swap(bs.act);
}

inline void unstageBucket(BucketSoA &bs) {
  suspendBucket(bs);
  bs.done = true;
}

/// Bucket overload of partBucketForce. When bs is staged the source term is
/// queued (flushed at TERM_BATCH_MAX); otherwise it runs the scalar path.
/// @param part       Source particle.
/// @param req        Bucket node the walk is evaluating.
/// @param particles  AoS particles of the TreePiece (unstaged path).
/// @param offset     Periodic offset applied to the source.
/// @param activeRung Lowest active rung (unstaged path).
/// @param bs         Bucket staging state.
/// @return           Number of active particles in the bucket.
inline int partBucketForce(ExternalGravityParticle *part,
			   Tree::GenericTreeNode *req,
			   GravityParticle *particles,
			   Vector3D<cosmoType> offset, int activeRung,
			   BucketSoA &bs) {
  // Staged: queue the source term in walk order.
  // Unstaged: scalar-equivalent path.
  if (bs.staged) {
    bs.order.push_back(~(int)bs.parts.size());
    bs.parts.push_back(BucketPartTerm{offset + part->position,
				      part->mass, part->soft});
    if ((int)bs.order.size() >= TERM_BATCH_MAX)
      flushTerms(bs);
    return bs.soa.n;
  }
  return partBucketForce(part, req, particles, offset, activeRung);
}

/// Bucket overload of nodeBucketForce. When bs is staged the node moments are
/// queued (or a particle term if softening opens the node); otherwise it runs
/// the scalar path.
/// @param node       Source tree node.
/// @param req        Bucket node the walk is evaluating.
/// @param particles  AoS particles of the TreePiece (unstaged path).
/// @param offset     Periodic offset applied to the node.
/// @param activeRung Lowest active rung (unstaged path).
/// @param bs         Bucket staging state.
/// @return           Number of active particles in the bucket.
inline int nodeBucketForce(Tree::GenericTreeNode *node,
			   Tree::GenericTreeNode *req,
			   GravityParticle *particles,
			   Vector3D<cosmoType> offset, int activeRung,
			   BucketSoA &bs) {
  // Staged: queue the node (or a particle term if softening opens it), as in
  // scalar nodeBucketForce. The HEXADECAPOLE branch matches the scalar path.
  if (bs.staged) {
    const MultipoleMoments &m = node->moments;
    Vector3D<cosmoType> cm(m.cm + offset);
#ifdef HEXADECAPOLE
    if(openSoftening(node, req, offset)) {
      bs.order.push_back(~(int)bs.parts.size());
      bs.parts.push_back(BucketPartTerm{offset + m.cm,
					m.totalMass, m.soft});
      if ((int)bs.order.size() >= TERM_BATCH_MAX)
	flushTerms(bs);
      return bs.soa.n;
    }
#endif
    bs.order.push_back((int)bs.nodes.size());
    bs.nodes.push_back(BucketNodeTerm{m, cm});
    if ((int)bs.order.size() >= TERM_BATCH_MAX)
      flushTerms(bs);
    return bs.soa.n;
  }
  return nodeBucketForce(node, req, particles, offset, activeRung);
}

#endif
