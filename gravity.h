#ifndef __GRAVITY_H__
#define __GRAVITY_H__

#include "TreeNode.h"
#include "GenericTreeNode.h"
#include "Space.h"
#include "SSEdefs.h"

extern double theta;
extern double thetaMono;

/*
** see (A1) and (A2) of TREESPH: A UNIFICATION OF SPH WITH THE 
** HIERARCHICAL TREE METHOD by Lars Hernquist and Neal Katz.
** APJ Supplemant Series 70:416-446, 1989
** 
** Higher derivative terms c and d for use with quadrupole spline
** softening (Joachim Stadel, Dec. 94).
*/
#ifndef __SSE2__
inline
void SPLINEQ(cosmoType invr,cosmoType r2,cosmoType twoh,cosmoType& a,
	     cosmoType& b,cosmoType& c,cosmoType& d)
{
  cosmoType u,dih,dir=(invr);
  if ((r2) < (twoh)*(twoh)) {
    dih = COSMO_CONST(2.0)/(twoh);
    u = dih/dir;
    if (u < COSMO_CONST(1.0)) {
      a = dih*(COSMO_CONST(7.0)/COSMO_CONST(5.0) 
	       - COSMO_CONST(2.0)/COSMO_CONST(3.0)*u*u 
	       + COSMO_CONST(3.0)/COSMO_CONST(10.0)*u*u*u*u
	       - COSMO_CONST(1.0)/COSMO_CONST(10.0)*u*u*u*u*u);
      b = dih*dih*dih*(COSMO_CONST(4.0)/COSMO_CONST(3.0) 
		       - COSMO_CONST(6.0)/COSMO_CONST(5.0)*u*u 
		       + COSMO_CONST(1.0)/COSMO_CONST(2.0)*u*u*u);
      c = dih*dih*dih*dih*dih*(COSMO_CONST(12.0)/COSMO_CONST(5.0) 
			       - COSMO_CONST(3.0)/COSMO_CONST(2.0)*u);
      d = COSMO_CONST(3.0)/COSMO_CONST(2.0)*dih*dih*dih*dih*dih*dih*dir;
    }
    else {
      a = COSMO_CONST(-1.0)/COSMO_CONST(15.0)*dir 
	+ dih*(COSMO_CONST(8.0)/COSMO_CONST(5.0) 
	       - COSMO_CONST(4.0)/COSMO_CONST(3.0)*u*u + u*u*u
	       - COSMO_CONST(3.0)/COSMO_CONST(10.0)*u*u*u*u 
	       + COSMO_CONST(1.0)/COSMO_CONST(30.0)*u*u*u*u*u);
      b = COSMO_CONST(-1.0)/COSMO_CONST(15.0)*dir*dir*dir 
	+ dih*dih*dih*(COSMO_CONST(8.0)/COSMO_CONST(3.0) - COSMO_CONST(3.0)*u 
		       + COSMO_CONST(6.0)/COSMO_CONST(5.0)*u*u 
		       - COSMO_CONST(1.0)/COSMO_CONST(6.0)*u*u*u);
      c = COSMO_CONST(-1.0)/COSMO_CONST(5.0)*dir*dir*dir*dir*dir 
	+ COSMO_CONST(3.0)*dih*dih*dih*dih*dir
	+ dih*dih*dih*dih*dih*(COSMO_CONST(-12.0)/COSMO_CONST(5.0) 
			       + COSMO_CONST(1.0)/COSMO_CONST(2.0)*u);
      d = -dir*dir*dir*dir*dir*dir*dir
	+ COSMO_CONST(3.0)*dih*dih*dih*dih*dir*dir*dir
	- COSMO_CONST(1.0)/COSMO_CONST(2.0)*dih*dih*dih*dih*dih*dih*dir;
    }
  }
  else {
    a = dir;
    b = a*a*a;
    c = COSMO_CONST(3.0)*b*a*a;
    d = COSMO_CONST(5.0)*c*a*a;
  }
}
#else
inline
void SPLINEQ(SSEcosmoType invr, SSEcosmoType r2, SSEcosmoType twoh,
	     SSEcosmoType& a, SSEcosmoType& b, 
	     SSEcosmoType& c, SSEcosmoType& d)
{
  SSEcosmoType dir;
  dir = invr; 
  SSEcosmoType select0 = r2 < twoh * twoh;
  int compare0 = movemask(select0);
  // make common case fast, at the cost of some redundant code
  //   in the infrequent case
  if (!compare0) {
    a = dir;
    b = a*a*a;
    c = COSMO_CONST(3.0)*b*a*a;
    d = COSMO_CONST(5.0)*c*a*a;
  }
  else {
    SSEcosmoType u, dih; 
    SSEcosmoType select1; 
    SSEcosmoType a0,b0,c0,d0,a1,b1,c1,d1,a2,b2,c2,d2;
    int compare1; 

    dih = COSMO_CONST(2.0)/twoh;
    u = dih/dir;

    if ((~compare0) & cosmoMask) {
      a0 = dir;
      b0 = a0*a0*a0;
      c0 = COSMO_CONST(3.0)*b0*a0*a0;
      d0 = COSMO_CONST(5.0)*c0*a0*a0;
    }

    select1 = u < COSMO_CONST(1.0);
    compare1 = movemask(select1);
    if (compare1) {
      a1 = dih*(COSMO_CONST(7.0)/COSMO_CONST(5.0) 
	       - COSMO_CONST(2.0)/COSMO_CONST(3.0)*u*u 
	       + COSMO_CONST(3.0)/COSMO_CONST(10.0)*u*u*u*u
	       - COSMO_CONST(1.0)/COSMO_CONST(10.0)*u*u*u*u*u);
      b1 = dih*dih*dih*(COSMO_CONST(4.0)/COSMO_CONST(3.0) 
		       - COSMO_CONST(6.0)/COSMO_CONST(5.0)*u*u 
		       + COSMO_CONST(1.0)/COSMO_CONST(2.0)*u*u*u);
      c1 = dih*dih*dih*dih*dih*(COSMO_CONST(12.0)/COSMO_CONST(5.0) 
			       - COSMO_CONST(3.0)/COSMO_CONST(2.0)*u);
      d1 = COSMO_CONST(3.0)/COSMO_CONST(2.0)*dih*dih*dih*dih*dih*dih*dir;
    }  
    if ((~compare1) & cosmoMask) {
      a2 = COSMO_CONST(-1.0)/COSMO_CONST(15.0)*dir 
	+ dih*(COSMO_CONST(8.0)/COSMO_CONST(5.0) 
	       - COSMO_CONST(4.0)/COSMO_CONST(3.0)*u*u + u*u*u
	       - COSMO_CONST(3.0)/COSMO_CONST(10.0)*u*u*u*u 
	       + COSMO_CONST(1.0)/COSMO_CONST(30.0)*u*u*u*u*u);
      b2 = COSMO_CONST(-1.0)/COSMO_CONST(15.0)*dir*dir*dir 
	+ dih*dih*dih*(COSMO_CONST(8.0)/COSMO_CONST(3.0) - COSMO_CONST(3.0)*u 
		       + COSMO_CONST(6.0)/COSMO_CONST(5.0)*u*u 
		       - COSMO_CONST(1.0)/COSMO_CONST(6.0)*u*u*u);
      c2 = COSMO_CONST(-1.0)/COSMO_CONST(5.0)*dir*dir*dir*dir*dir 
	+ COSMO_CONST(3.0)*dih*dih*dih*dih*dir
	+ dih*dih*dih*dih*dih*(COSMO_CONST(-12.0)/COSMO_CONST(5.0) 
			       + COSMO_CONST(1.0)/COSMO_CONST(2.0)*u);
      d2 = -dir*dir*dir*dir*dir*dir*dir
	+ COSMO_CONST(3.0)*dih*dih*dih*dih*dir*dir*dir
	- COSMO_CONST(1.0)/COSMO_CONST(2.0)*dih*dih*dih*dih*dih*dih*dir;
    }

    a = andnot(select0, a0) 
      | (select0 & ((select1 & a1) | andnot(select1, a2))); 
    b = andnot(select0, b0) 
      | (select0 & ((select1 & b1) | andnot(select1, b2))); 
    c = andnot(select0, c0) 
      | (select0 & ((select1 & c1) | andnot(select1, c2))); 
    d = andnot(select0, d0) 
      | (select0 & ((select1 & d1) | andnot(select1, d2))); 
  }
}
#endif

#ifndef __SSE2__
inline
void SPLINE(cosmoType r2, cosmoType twoh, cosmoType &a, cosmoType &b)
{
  cosmoType r, u,dih,dir;
  r = sqrt(r2);

  if (r < (twoh)) {
    dih = COSMO_CONST(2.0)/(twoh);
    u = r*dih;
    if (u < COSMO_CONST(1.0)) {
      a = dih*(COSMO_CONST(7.0)/COSMO_CONST(5.0) 
	       - COSMO_CONST(2.0)/COSMO_CONST(3.0)*u*u 
	       + COSMO_CONST(3.0)/COSMO_CONST(10.0)*u*u*u*u
	       - COSMO_CONST(1.0)/COSMO_CONST(10.0)*u*u*u*u*u);
      b = dih*dih*dih*(COSMO_CONST(4.0)/COSMO_CONST(3.0) 
		       - COSMO_CONST(6.0)/COSMO_CONST(5.0)*u*u 
		       + COSMO_CONST(1.0)/COSMO_CONST(2.0)*u*u*u);
    }
    else {
      dir = COSMO_CONST(1.0)/r;
      a = COSMO_CONST(-1.0)/COSMO_CONST(15.0)*dir 
	+ dih*(COSMO_CONST(8.0)/COSMO_CONST(5.0) 
	       - COSMO_CONST(4.0)/COSMO_CONST(3.0)*u*u + u*u*u
	       - COSMO_CONST(3.0)/COSMO_CONST(10.0)*u*u*u*u 
	       + COSMO_CONST(1.0)/COSMO_CONST(30.0)*u*u*u*u*u);
      b = COSMO_CONST(-1.0)/COSMO_CONST(15.0)*dir*dir*dir 
	+ dih*dih*dih*(COSMO_CONST(8.0)/COSMO_CONST(3.0) - COSMO_CONST(3.0)*u 
		       + COSMO_CONST(6.0)/COSMO_CONST(5.0)*u*u 
		       - COSMO_CONST(1.0)/COSMO_CONST(6.0)*u*u*u);
    }
  }
  else {
    a = COSMO_CONST(1.0)/r;
    b = a*a*a;
  }
}
#else 
inline
void SPLINE(SSEcosmoType r2, SSEcosmoType twoh, 
	    SSEcosmoType &a, SSEcosmoType &b)
{
  SSEcosmoType r;
  r = sqrt(r2);
  SSEcosmoType select0 = r < twoh;
  int compare0 = movemask(select0);
  // make common case fast, at the cost of some redundant code
  //   in the infrequent case
  if (!compare0) {
    a = COSMO_CONST(1.0)/r; 
    b = a * a * a; 
  }
  else {
    SSEcosmoType u, dih, dir; 
    SSEcosmoType a0, b0, a1, b1, a2, b2;
    SSEcosmoType select1; 
    int compare1; 

    dih = COSMO_CONST(2.0)/twoh;
    u = r * dih;
    
    if ((~compare0) & cosmoMask) {
      a0 = COSMO_CONST(1.0)/r; 
      b0 = a0 * a0 * a0; 
    }

    select1 = u < COSMO_CONST(1.0);
    compare1 = movemask(select1);
    if (compare1) {
      a1 = dih*(COSMO_CONST(7.0)/COSMO_CONST(5.0) 
		- COSMO_CONST(2.0)/COSMO_CONST(3.0)*u*u 
		+ COSMO_CONST(3.0)/COSMO_CONST(10.0)*u*u*u*u
		- COSMO_CONST(1.0)/COSMO_CONST(10.0)*u*u*u*u*u);
      b1 = dih*dih*dih*(COSMO_CONST(4.0)/COSMO_CONST(3.0) 
			- COSMO_CONST(6.0)/COSMO_CONST(5.0)*u*u 
			+ COSMO_CONST(1.0)/COSMO_CONST(2.0)*u*u*u);     
    }
    if ((~compare1) & cosmoMask) {
      dir = COSMO_CONST(1.0)/r;
      a2 = COSMO_CONST(-1.0)/COSMO_CONST(15.0)*dir 
	+ dih*(COSMO_CONST(8.0)/COSMO_CONST(5.0) 
	       - COSMO_CONST(4.0)/COSMO_CONST(3.0)*u*u + u*u*u
	       - COSMO_CONST(3.0)/COSMO_CONST(10.0)*u*u*u*u 
	       + COSMO_CONST(1.0)/COSMO_CONST(30.0)*u*u*u*u*u);
      b2 = COSMO_CONST(-1.0)/COSMO_CONST(15.0)*dir*dir*dir 
	+ dih*dih*dih*(COSMO_CONST(8.0)/COSMO_CONST(3.0) - COSMO_CONST(3.0)*u 
		       + COSMO_CONST(6.0)/COSMO_CONST(5.0)*u*u 
		       - COSMO_CONST(1.0)/COSMO_CONST(6.0)*u*u*u);
    }

    a = andnot(select0, a0) 
      | (select0 & ((select1 & a1) | andnot(select1, a2))); 
    b = andnot(select0, b0) 
      | (select0 & ((select1 & b1) | andnot(select1, b2))); 
  }	
}
#endif

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

#ifndef __SSE2__
inline int partBucketForce(ExternalGravityParticle *part, 
			   Tree::GenericTreeNode *req, 
			   GravityParticle *particles, 
			   Vector3D<cosmoType> offset, int activeRung) {
  int computed = 0;
  Vector3D<cosmoType> r;
  cosmoType rsq;
  cosmoType twoh, a, b;

  for(int j = req->firstParticle; j <= req->lastParticle; ++j) {
    if (particles[j].rung >= activeRung) {
#ifdef CMK_VERSION_BLUEGENE
      if (++forProgress > 200) {
        forProgress = 0;
#ifdef COSMO_EVENTS
        traceUserEvents(networkProgressUE);
#endif
        CmiNetworkProgress();
      }
#endif
      computed++;
      r = offset + part->position - particles[j].position;
      rsq = r.lengthSquared();
      twoh = part->soft + particles[j].soft;
      if(rsq != 0) {
        SPLINE(rsq, twoh, a, b);
	cosmoType idt2 = (particles[j].mass + part->mass)*b; // (timescale)^-2
	// of interaction
        particles[j].treeAcceleration += r * (b * part->mass);
        particles[j].potential -= part->mass * a;
	if(idt2 > particles[j].dtGrav)
	  particles[j].dtGrav = idt2;
      }
    }
  }
  return computed;
}
#else
inline int partBucketForce(ExternalGravityParticle *part, 
			   Tree::GenericTreeNode *req, 
			   GravityParticle **activeParticles, 
			   Vector3D<cosmoType> offset, int nActiveParts) {
  Vector3D<SSEcosmoType> r;
  SSEcosmoType rsq; 
  SSEcosmoType twoh, a, b; 

  for (int i=0; i<nActiveParts; i+=SSE_VECTOR_WIDTH) {
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
      packedPos(SSELoad(SSEcosmoType, activeParticles, i, ->position.x),
		SSELoad(SSEcosmoType, activeParticles, i, ->position.y),
		SSELoad(SSEcosmoType, activeParticles, i, ->position.z)); 
    SSELoad(SSEcosmoType packedSoft, activeParticles, i, ->soft); 

    r = -packedPos + offset + part->position; 
    rsq = r.lengthSquared();
    twoh = part->soft + packedSoft; 
    SSEcosmoType select = rsq > COSMO_CONST(0.0); 
    int compare = movemask(select); 
    if(compare) {
      SPLINE(rsq, twoh, a, b);
      if ((~compare) & cosmoMask) {
	a = select & a; 
	b = select & b; 
      }
      SSEcosmoType SSELoad(packedMass, activeParticles, i, ->mass);  
      SSEcosmoType SSELoad(packedDtGrav, activeParticles, i, ->dtGrav);
      Vector3D<SSEcosmoType> 
	packedAcc(SSELoad(SSEcosmoType, 
			  activeParticles, i, ->treeAcceleration.x),
		  SSELoad(SSEcosmoType, 
			  activeParticles, i, ->treeAcceleration.y),
		  SSELoad(SSEcosmoType, 
			  activeParticles, i, ->treeAcceleration.z)); 
      SSEcosmoType SSELoad(packedPotential, activeParticles, i, ->potential); 
      SSEcosmoType idt2 = (packedMass + part->mass) * b;       
      idt2 = max(idt2, packedDtGrav); 
      packedAcc += r * (b * part->mass);
      packedPotential -= part->mass * a; 
      SSEStore(packedAcc.x, activeParticles, i, ->treeAcceleration.x);
      SSEStore(packedAcc.y, activeParticles, i, ->treeAcceleration.y);
      SSEStore(packedAcc.z, activeParticles, i, ->treeAcceleration.z);
      SSEStore(packedPotential, activeParticles, i, ->potential);
      SSEStore(idt2, activeParticles, i, ->dtGrav); 
    }
  }
  return nActiveParts;
}

inline int partBucketForce(ExternalGravityParticle *part, 
			   Tree::GenericTreeNode *req, 
			   GravityParticle *particles, 
			   Vector3D<cosmoType> offset, int activeRung) {
  int nActiveParts = 0; 
  GravityParticle dummyPart = particles[0];
  dummyPart.soft = 0.0;

  GravityParticle **activeParticles = 
    (GravityParticle**)alloca((req->lastParticle - req->firstParticle 
			       + 1 + FORCE_INPUT_LIST_PAD) 
			      * sizeof(GravityParticle*));
  for (int j=req->firstParticle; j <= req->lastParticle; ++j) {
    if (particles[j].rung >= activeRung) 
      activeParticles[nActiveParts++] = &particles[j];
  }
  
  activeParticles[nActiveParts] = &dummyPart; 
#ifdef COSMO_FLOAT
  activeParticles[nActiveParts+1] = &dummyPart; 
  activeParticles[nActiveParts+2] = &dummyPart; 
#endif

  return partBucketForce(part, req, activeParticles, offset, nActiveParts); 
}

#endif


//
// Calculated forces on active particles in a bucket due to the
// multipole of a TreeNode.  Return number of multipoles evaluated.
//
#if  !defined(__SSE2__) 
inline
int nodeBucketForce(Tree::GenericTreeNode *node, 
		    Tree::GenericTreeNode *req,  
		    GravityParticle *particles, 
		    Vector3D<cosmoType> offset,    
		    int activeRung)

{
  int computed = 0;
  Vector3D<cosmoType> r;
  cosmoType rsq;
  cosmoType twoh, a, b, c, d;
  MultipoleMoments &m = node->moments;

  Vector3D<cosmoType> cm(m.cm + offset);

#ifdef HEXADECAPOLE
  if(openSoftening(node, req, offset)) {
    ExternalGravityParticle tmpPart;
    tmpPart.mass = m.totalMass;
    tmpPart.soft = m.soft;
    tmpPart.position = m.cm;
    return partBucketForce(&tmpPart, req, particles, offset, activeRung);
  }
#endif
  for(int j = req->firstParticle; j <= req->lastParticle; ++j) {
    if (particles[j].rung >= activeRung) {
      particles[j].interMass += m.totalMass;
#ifdef CMK_VERSION_BLUEGENE
      if (++forProgress > 200) {
        forProgress = 0;
#ifdef COSMO_EVENTS
        traceUserEvents(networkProgressUE);
#endif
        CmiNetworkProgress();
      }
#endif
      computed++;
      r = Vector3D<cosmoType>(particles[j].position) - cm;
      rsq = r.lengthSquared();
      cosmoType dir = COSMO_CONST(1.0)/sqrt(rsq);
#ifdef HEXADECAPOLE
      momEvalMomr(&m.mom, dir, -r.x, -r.y, -r.z, &particles[j].potential,
		  &particles[j].treeAcceleration.x,
		  &particles[j].treeAcceleration.y,
		  &particles[j].treeAcceleration.z);
      cosmoType idt2 = (particles[j].mass + m.totalMass)*dir*dir*dir;
#else
      twoh = CONVERT_TO_COSMO_TYPE(m.soft + particles[j].soft);
      SPLINEQ(dir, rsq, twoh, a, b, c, d);
      cosmoType qirx = CONVERT_TO_COSMO_TYPE m.xx*r.x 
	+ CONVERT_TO_COSMO_TYPE m.xy*r.y + CONVERT_TO_COSMO_TYPE m.xz*r.z;
      cosmoType qiry = CONVERT_TO_COSMO_TYPE m.xy*r.x 
	+ CONVERT_TO_COSMO_TYPE m.yy*r.y + CONVERT_TO_COSMO_TYPE m.yz*r.z;
      cosmoType qirz = CONVERT_TO_COSMO_TYPE m.xz*r.x 
	+ CONVERT_TO_COSMO_TYPE m.yz*r.y + CONVERT_TO_COSMO_TYPE m.zz*r.z;
      cosmoType qir = COSMO_CONST(0.5)*(qirx*r.x + qiry*r.y + qirz*r.z);
      cosmoType tr = COSMO_CONST(0.5)*(CONVERT_TO_COSMO_TYPE m.xx 
				       + CONVERT_TO_COSMO_TYPE m.yy 
				       + CONVERT_TO_COSMO_TYPE m.zz);
      cosmoType qir3 = b*CONVERT_TO_COSMO_TYPE m.totalMass + d*qir - c*tr;
      particles[j].potential -= CONVERT_TO_COSMO_TYPE m.totalMass * a 
	+ c*qir - b*tr;
      particles[j].treeAcceleration.x -= qir3*r.x - c*qirx;
      particles[j].treeAcceleration.y -= qir3*r.y - c*qiry;
      particles[j].treeAcceleration.z -= qir3*r.z - c*qirz;
      cosmoType idt2 = (CONVERT_TO_COSMO_TYPE particles[j].mass 
			+ CONVERT_TO_COSMO_TYPE m.totalMass)*b;
#endif
      if(idt2 > particles[j].dtGrav)
        particles[j].dtGrav = idt2;
    }
  }
  return computed;
}

#elif defined(__SSE2__) 
inline
int nodeBucketForce(Tree::GenericTreeNode *node, 
		    Tree::GenericTreeNode *req,  
		    GravityParticle *particles, 
		    Vector3D<cosmoType> offset,    
		    int activeRung)

{
  Vector3D<SSEcosmoType> r; 
  SSEcosmoType rsq;
  SSEcosmoType twoh;
  SSEcosmoType a,b,c,d;
  MultipoleMoments m = node->moments;
  Vector3D<cosmoType> cm(m.cm + offset);
  int nActiveParts = 0;
  GravityParticle dummyPart = particles[0];
  dummyPart.soft = 0.0;

  GravityParticle **activeParticles = 
    (GravityParticle**)alloca((req->lastParticle - req->firstParticle 
			       + 1 + FORCE_INPUT_LIST_PAD) 
			      * sizeof(GravityParticle*));
  for (int j=req->firstParticle; j <= req->lastParticle; ++j) {
    if (particles[j].rung >= activeRung) 
      activeParticles[nActiveParts++] = &particles[j];
  }

  activeParticles[nActiveParts] = &dummyPart; 
#ifdef COSMO_FLOAT
  activeParticles[nActiveParts+1] = &dummyPart; 
  activeParticles[nActiveParts+2] = &dummyPart; 
#endif

#ifdef HEXADECAPOLE
  if(openSoftening(node, req, offset)) {
    ExternalGravityParticle tmpPart;
    tmpPart.mass = m.totalMass;
    tmpPart.soft = m.soft;
    tmpPart.position = m.cm;
    return partBucketForce(&tmpPart, req, activeParticles, 
			   offset, nActiveParts);
  }
#endif
  for (int i=0; i<nActiveParts; i+=SSE_VECTOR_WIDTH) {
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
      packedPos(SSELoad(SSEcosmoType, activeParticles, i, ->position.x),
		SSELoad(SSEcosmoType, activeParticles, i, ->position.y),
		SSELoad(SSEcosmoType, activeParticles, i, ->position.z)); 
    r = packedPos - cm; 
    rsq = r.lengthSquared();
    SSEcosmoType dir = COSMO_CONST(1.0)/sqrt(rsq);
    Vector3D<SSEcosmoType> 
      packedAcc(SSELoad(SSEcosmoType, 
			activeParticles, i, ->treeAcceleration.x),
		SSELoad(SSEcosmoType, 
			activeParticles, i, ->treeAcceleration.y),
		SSELoad(SSEcosmoType, 
			activeParticles, i, ->treeAcceleration.z)); 
    SSEcosmoType SSELoad(packedPotential, activeParticles, i, ->potential); 
    SSEcosmoType SSELoad(packedMass, activeParticles, i, ->mass);  
    SSEcosmoType SSELoad(packedDtGrav, activeParticles, i, ->dtGrav);
#ifdef HEXADECAPOLE
    momEvalMomr(&m.mom, dir, -r.x, -r.y, -r.z, &packedPotential,
		&packedAcc.x,
		&packedAcc.y,
		&packedAcc.z);
    SSEcosmoType idt2 = (packedMass + m.totalMass)*dir*dir*dir;
#else
    SSELoad(SSEcosmoType packedSoft, activeParticles, i, ->soft); 
    twoh = CONVERT_TO_COSMO_TYPE m.soft + packedSoft;
    SPLINEQ(dir, rsq, twoh, a, b, c, d);
    SSEcosmoType qirx = CONVERT_TO_COSMO_TYPE m.xx*r.x 
      + CONVERT_TO_COSMO_TYPE m.xy*r.y + CONVERT_TO_COSMO_TYPE m.xz*r.z;
    SSEcosmoType qiry = CONVERT_TO_COSMO_TYPE m.xy*r.x 
      + CONVERT_TO_COSMO_TYPE m.yy*r.y + CONVERT_TO_COSMO_TYPE m.yz*r.z;
    SSEcosmoType qirz = CONVERT_TO_COSMO_TYPE m.xz*r.x 
      + CONVERT_TO_COSMO_TYPE m.yz*r.y + CONVERT_TO_COSMO_TYPE m.zz*r.z;
    SSEcosmoType qir = COSMO_CONST(0.5)*(qirx*r.x + qiry*r.y + qirz*r.z);
    SSEcosmoType tr = COSMO_CONST(0.5)*(CONVERT_TO_COSMO_TYPE m.xx 
					+ CONVERT_TO_COSMO_TYPE m.yy 
					+ CONVERT_TO_COSMO_TYPE m.zz);
    SSEcosmoType qir3 = b*CONVERT_TO_COSMO_TYPE m.totalMass + d*qir - c*tr;
    packedPotential -= CONVERT_TO_COSMO_TYPE m.totalMass *a 
      + c * qir - b * tr;      
    packedAcc.x -= qir3*r.x - c*qirx;
    packedAcc.y -= qir3*r.y - c*qiry;
    packedAcc.z -= qir3*r.z - c*qirz;
    SSEcosmoType idt2 = (packedMass + CONVERT_TO_COSMO_TYPE m.totalMass)*b;
#endif
    SSEStore(packedAcc.x, activeParticles, i, ->treeAcceleration.x);
    SSEStore(packedAcc.y, activeParticles, i, ->treeAcceleration.y);
    SSEStore(packedAcc.z, activeParticles, i, ->treeAcceleration.z);
    SSEStore(packedPotential, activeParticles, i, ->potential);
    idt2 = max(idt2, packedDtGrav);   
    SSEStore(idt2, activeParticles, i, ->dtGrav); 
  }
  return nActiveParts;
}
#endif

/// @brief Gravity opening criterion for a bucket walk.
/// @param node Source node to be tested
/// @param bucketNode Target bucket
/// @param offset Offset of periodic replica
/// @param localIndex Index of requesting TreePiece
/// @return True if the node's opening radius intersects the
/// boundingBox of the bucket, i.e. the node needs to be opened.

inline bool
openCriterionBucket(Tree::GenericTreeNode *node,
                   Tree::GenericTreeNode *bucketNode,
                   Vector3D<double> offset, // Offset of node
                   int localIndex // requesting TreePiece
                   ) {
  // mark the node as used by the requesting TreePiece
  node->markUsedBy(localIndex);

#if COSMO_STATS > 0
  node->used = true;
#endif
  // Always open node if this many particles or fewer.
  const int nMinParticleNode = 6;
  if(node->particleCount < nMinParticleNode) {
      return true;
      }

  // Note that some of this could be pre-calculated into an "opening radius"
  double radius = TreeStuff::opening_geometry_factor * node->moments.radius / theta;
  if(radius < node->moments.radius)
      radius = node->moments.radius;

  Sphere<double> s(node->moments.cm + offset, radius);
  
#ifdef HEXADECAPOLE
  if(!Space::intersect(bucketNode->boundingBox, s)) {
      // Well separated, now check softening
      if(!openSoftening(node, bucketNode, offset)) {
	  return false; // passed both tests: will be a Hex interaction
      }
      else {        // Open as monopole?
        radius = TreeStuff::opening_geometry_factor*node->moments.radius/thetaMono;
      Sphere<double> sM(node->moments.cm + offset, radius);
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
/// @param localIndex Index of requesting TreePiece
/// @return -1 if there is an intersection
///    0 if no intersection
///    +1 if completely contained
///
/// This is the criterion described in Stadel(2001): if no
/// intersection, then the node is valid for all particles in myNode.
/// If completely contained then it needs to be opened for all
/// particles in myNode.  If a partial intersection, then node needs
/// to be checked by the children of myNode.  If myNode is already a
/// bucket, then only returns 0 or 1 are possible.

inline int openCriterionNode(Tree::GenericTreeNode *node,
                    Tree::GenericTreeNode *myNode,
                    Vector3D<double> offset,
                    int localIndex // requesting TreePiece
                    ) {
  // mark the node as used by this TreePiece
  node->markUsedBy(localIndex);

#if COSMO_STATS > 0
  node->used = true;
#endif
  // Always open node if this many particles or fewer.
  const int nMinParticleNode = 6;
  if(node->particleCount < nMinParticleNode) {
      return 1;
      }

  // Note that some of this could be pre-calculated into an "opening radius"
  double radius = TreeStuff::opening_geometry_factor * node->moments.radius / theta;
  if(radius < node->moments.radius)
      radius = node->moments.radius;

  Sphere<double> s(node->moments.cm + offset, radius);

  if(myNode->getType()==Tree::Bucket || myNode->getType()==Tree::CachedBucket || myNode->getType()==Tree::NonLocalBucket){
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
          radius = TreeStuff::opening_geometry_factor*node->moments.radius/thetaMono;
            Sphere<double> sM(node->moments.cm + offset, radius);
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
                radius = TreeStuff::opening_geometry_factor*node->moments.radius/thetaMono;
                Sphere<double> sM(node->moments.cm + offset, radius);
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

#endif
