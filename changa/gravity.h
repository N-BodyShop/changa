#ifndef __GRAVITY_H__
#define __GRAVITY_H__

#include "TreeNode.h"
#include "GenericTreeNode.h"
#include "Space.h"

extern double theta;
extern double thetaMono;

#if  !defined(__SSE2__) && defined(COSMO_FLOAT)
#define CONVERT_TO_COSMO_TYPE (float)  
#elif !defined(__SSE2__) && !defined(COSMO_FLOAT)
#define CONVERT_TO_COSMO_TYPE 
#elif defined(__SSE2__) && defined(COSMO_FLOAT)
#include "SSE-Float.h"
typedef SSEFloat SSEcosmoType; 
#elif defined(__SSE2__) && !defined(COSMO_FLOAT)
#include "SSE-Double.h"
typedef SSEDouble SSEcosmoType; 
#endif

#ifdef COSMO_FLOAT
typedef float cosmoType;
#define COSMO_CONST(val) val##f
#else
typedef double cosmoType; 
#define COSMO_CONST(val) val
#endif

/*
** see (A1) and (A2) of TREESPH: A UNIFICATION OF SPH WITH THE 
** HIERARCHICAL TREE METHOD by Lars Hernquist and Neal Katz.
** APJ Supplemant Series 70:416-446, 1989
** 
** Higher derivative terms c and d for use with quadrupole spline
** softening (Joachim Stadel, Dec. 94).
*/
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

/*
inline void SPLINEM(double invr,double r2,double twoh,double& a,double& b)
{
  double u,dih,dir=(invr);
  if ((r2) < (twoh)*(twoh)) {
    dih = 2.0/(twoh);
    u = dih/dir;
    if (u < 1.0) {
      a = dih*(7.0/5.0 - 2.0/3.0*u*u + 3.0/10.0*u*u*u*u
	       - 1.0/10.0*u*u*u*u*u);
      b = dih*dih*dih*(4.0/3.0 - 6.0/5.0*u*u + 1.0/2.0*u*u*u);
    }
    else {
      a = -1.0/15.0*dir + dih*(8.0/5.0 - 4.0/3.0*u*u + u*u*u
			       - 3.0/10.0*u*u*u*u + 1.0/30.0*u*u*u*u*u);
      b = -1.0/15.0*dir*dir*dir 
	+ dih*dih*dih*(8.0/3.0 - 3.0*u + 6.0/5.0*u*u - 1.0/6.0*u*u*u);
    }
  }
  else {
    a = dir;
    b = a*a*a;
  }
}
*/

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

//
// Return true if the soften nodes overlap, i.e. the forces involve softening
// @param node
// @param myNode
// @param offset Periodic offset applied to node
//
inline
int openSoftening(Tree::GenericTreeNode *node, Tree::GenericTreeNode *myNode,
		  Vector3D<cosmoType> offset)
{
  Sphere<cosmoType> s(node->moments.cm + offset, 2.0*node->moments.soft);
  Sphere<cosmoType> myS(myNode->moments.cm, 2.0*myNode->moments.soft);
  return Space::intersect(myS, s);
}

#ifdef CMK_VERSION_BLUEGENE
static int forProgress = 0;
#endif

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
  MultipoleMoments m = node->moments;

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
#ifdef HEXADECAPOLE

      cosmoType dir = COSMO_CONST(1.0)/sqrt(rsq);
      momEvalMomr(&m.mom, dir, -r.x, -r.y, -r.z, &particles[j].potential,
		  &particles[j].treeAcceleration.x,
		  &particles[j].treeAcceleration.y,
		  &particles[j].treeAcceleration.z);

      double idt2 = (particles[j].mass + m.totalMass)*dir*dir*dir;
#else
      twoh = CONVERT_TO_COSMO_TYPE(m.soft + particles[j].soft);
      cosmoType dir = COSMO_CONST(1.0)/sqrt(rsq);
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

#elif defined(__SSE2__) && !defined(HEXADECAPOLE) && !defined(COSMO_FLOAT)

inline
int nodeBucketForce(Tree::GenericTreeNode *node,
                    Tree::GenericTreeNode *req,
                    GravityParticle *particles,
                    Vector3D<double> offset,
                    int activeRung)

{
  int computed = 0;
  SSEDouble rx,ry,rz,select2;
  SSEDouble rsq,twoh,a,b,c,d,a1,a2,b1,b2,c1,c2,d1,d2;
  MultipoleMoments m = node->moments;
  SSEDouble mtotalMass(m.totalMass);
  GravityParticle dummyPart = particles[0];
  Vector3D<double> cm(m.cm + offset);
  SSEDouble tr(0.5*(m.xx + m.yy + m.zz));

  GravityParticle **input = (GravityParticle**)alloca((req->lastParticle - req->firstParticle + 2) * sizeof(GravityParticle*));

  int activeParts = 0;
  for (int j=req->firstParticle; j <= req->lastParticle; ++j) {
    if (particles[j].rung >= activeRung) input[activeParts++] = &particles[j];
  }
  input[activeParts] = &dummyPart;

  for (int i=0; i<activeParts; i+=2) {
    computed+=2;

    SSEDouble accX(0.0),accY(0.0),accZ(0.0),potential(0.0);
      
    SSEDouble posX(input[i]->position.x, input[i+1]->position.x);
    SSEDouble posY(input[i]->position.y, input[i+1]->position.y);
    SSEDouble posZ(input[i]->position.z, input[i+1]->position.z);
    SSEDouble cmx(cm.x); SSEDouble cmy(cm.y); SSEDouble cmz(cm.z); 

    rx = posX - cmx;
    ry = posY - cmy;
    rz = posZ - cmz;
     
    rsq= rx*rx+(ry*ry+rz*rz);
    SSEDouble k(m.soft);SSEDouble t(input[i]->soft,input[i+1]->soft);
    twoh = k+ t ;
   
    SSEDouble dir =  (1.0) / sqrt(rsq);
    
    SSEDouble dir2 = dir*dir;
    SSEDouble dir4 = dir2*dir2;
    SSEDouble select = rsq < (twoh * twoh);
    int compare =movemask(select);

    if (compare) {
      SSEDouble dih = 2.0/twoh;
      SSEDouble u = dih /dir;
      select2 = u < 1.0;
      int compare2 = movemask(select2);
      SSEDouble u2 = u *u;
      SSEDouble u4 = u2 * u2;
      SSEDouble dih2 =dih *dih;
      SSEDouble dih4 =dih2 *dih2;
      if (compare2) {
        a = dih * ( (7.0/5.0) - ((2.0/3.0)* u2) + ((3.0/10.0)*u4) - ((1.0/10.0)* u4 *u));
        b = dih2*dih*(4.0/3.0-(6.0/5.0*u2) + (1.0/2.0*u2*u));
        c = dih4 * dih *(12.0/5.0-(3.0/2.0*u));
        d = 3.0/2.0*dih4*dih2*dir;
      }
       
 
      if ((~compare2) & 0x3) {
	a2 = -1.0/15.0*dir + dih *(8.0/5.0-(4.0/3.0*u2)+u2*u - 3.0/10.0*u4 + 1.0/30.0*u4*u); 
	b2 = -1.0/15.0*dir2*dir + dih2*dih*(8.0/3.0 - 3.0 * u + 6.0/5.0*u2 -1.0/6.0*u2*u);
	c2 = -1.0/5.0*dir4*dir + 3.0 * dih4 * dir + dih4*dih*(-12.0/5.0 + 1.0/2.0*u);
	d2 = -1.0*dir4*dir2*dir + 3.0 *dih4 * dir2 * dir - 1.0/2.0*dih4*dih2*dir;
      }

      a1 = ((select2 & a) | andnot(select2, a2));
      b1 = ((select2 & b) | andnot(select2, b2));
      c1 = ((select2 & c) | andnot(select2, c2));
      d1 = ((select2 & d) | andnot(select2, d2));
    }
    if ((~compare) & 0x3) {
      a2 = dir;
      b2 = a2* a2*a2;
      c2 = 3.0*b2*a2*a2;
      d2 = 5.0*c2*a2*a2;
    }
      
    a = ((select & a1) | andnot(select, a2));
    b = ((select & b1) | andnot(select, b2));
    c = ((select & c1) | andnot(select, c2));
    d = ((select & d1) | andnot(select, d2));

    SSEDouble mxx(m.xx);SSEDouble mxy(m.xy);SSEDouble myy(m.yy); SSEDouble myz(m.yz);SSEDouble mzz(m.zz);SSEDouble mxz(m.xz);   
    SSEDouble qirx = mxx*rx + mxy*ry + mxz*rz;
    SSEDouble qiry = mxy*rx + myy*ry + myz*rz;
    SSEDouble qirz = mxz*rx + myz*ry + mzz*rz;
    SSEDouble qir = 0.5 * (qirx*rx+qiry*ry+qirz*rz);
    SSEDouble qir3 = m.totalMass * b  + d*qir - c *tr;

    SSEDouble temppotential(input[i]->potential,input[i+1]->potential);
    potential = temppotential - (m.totalMass*a + c*qir - b*tr);

    SSEDouble tempaccX(input[i]->treeAcceleration.x,input[i+1]->treeAcceleration.x);
    accX = tempaccX -( (qir3 *rx) -(c *qirx));
    SSEDouble tempaccY(input[i]->treeAcceleration.y,input[i+1]->treeAcceleration.y);
    accY = tempaccY -( (qir3 *ry) -(c *qiry));
    SSEDouble tempaccZ(input[i]->treeAcceleration.z,input[i+1]->treeAcceleration.z);
    accZ = tempaccZ - (qir3 *rz -c *qirz);

    SSEDouble tempmass(input[i]->mass,input[i+1]->mass);
    SSEDouble idt2 = (m.totalMass + tempmass) * b;
    SSEDouble dtGrav(input[i]->dtGrav, input[i+1]->dtGrav);
    idt2 = max(idt2,dtGrav);
    storel(&input[i]->potential, potential);
    storeh(&input[i+1]->potential, potential);

    storel(&input[i]->treeAcceleration.x, accX);
    storeh(&input[i+1]->treeAcceleration.x, accX);
    storel(&input[i]->treeAcceleration.y, accY);
    storeh(&input[i+1]->treeAcceleration.y, accY);
    storel(&input[i]->treeAcceleration.z, accZ);
    storeh(&input[i+1]->treeAcceleration.z, accZ);
    storel(&input[i]->dtGrav, idt2);
    storeh(&input[i+1]->dtGrav, idt2);
  }

  return computed;
}  
     
#elif  defined(__SSE2__) && defined(COSMO_FLOAT)  && !defined(HEXADECAPOLE)
inline
int nodeBucketForce(Tree::GenericTreeNode *node,
                    Tree::GenericTreeNode *req,
                    GravityParticle *particles,
                    Vector3D<float> offset,

                    int activeRung)

{
  int computed = 0;
  SSEFloat rx,ry,rz,select2;
  SSEFloat rsq,twoh,a,b,c,d,a1,a2,b1,b2,c1,c2,d1,d2;
  MultipoleMoments m = node->moments;
  SSEFloat mtotalMass(m.totalMass);
  GravityParticle dummyPart = particles[0];
  Vector3D<float> cm(m.cm + offset);
  SSEFloat tr(0.5*(m.xx + m.yy + m.zz));
  GravityParticle **input = (GravityParticle**)alloca((req->lastParticle - req->firstParticle + 5) * sizeof(GravityParticle*));
  int activeParts = 0;
  for (int j=req->firstParticle; j <= req->lastParticle; ++j) {

    if (particles[j].rung >= activeRung) input[activeParts++] = &particles[j];
  }
  input[activeParts] = &dummyPart;

  int rem = 3-(activeParts % 4);  
  for(int cnt =0;cnt<rem;cnt++)
    input[++activeParts] = &dummyPart;
  for (int i=0; i<activeParts; i+=4) {
    computed+=4;
    SSEFloat accX(0.0),accY(0.0),accZ(0.0),potential(0.0);
    SSEFloat posX(input[i]->position.x, input[i+1]->position.x,input[i+2]->position.x, input[i+3]->position.x);
    SSEFloat posY(input[i]->position.y, input[i+1]->position.y,input[i+2]->position.y, input[i+3]->position.y);
    SSEFloat posZ(input[i]->position.z, input[i+1]->position.z,input[i+2]->position.z, input[i+3]->position.z);
    SSEFloat cmx(cm.x); SSEFloat cmy(cm.y); SSEFloat cmz(cm.z);
    rx = posX - cmx;
    ry = posY - cmy;
    rz = posZ - cmz;
    rsq= rx*rx+(ry*ry+rz*rz);
    SSEFloat k(m.soft);	SSEFloat t(input[i]->soft,input[i+1]->soft,input[i+2]->soft,input[i+3]->soft);
    twoh = k+ t ; 
    SSEFloat dir =  (1.0) / sqrt(rsq);
    SSEFloat dir2 = dir*dir;
    SSEFloat dir4 = dir2*dir2;
    SSEFloat select = rsq < (twoh * twoh);
    int compare =movemask(select);
    if (compare) {
      SSEFloat dih = 2.0/twoh;
      SSEFloat u = dih /dir;
      select2 = u < 1.0;
      int compare2 = movemask(select2);
      SSEFloat u2 = u *u;
      SSEFloat u4 = u2 * u2;
      SSEFloat dih2 =dih *dih;
      SSEFloat dih4 =dih2 *dih2;

      if (compare2) {
        a = dih * ( (7.0/5.0) - ((2.0/3.0)* u2) + ((3.0/10.0)*u4) - ((1.0/10.0)* u4 *u));
        b = dih2*dih*(4.0/3.0-(6.0/5.0*u2) + (1.0/2.0*u2*u));
        c = dih4 * dih *(12.0/5.0-(3.0/2.0*u));
        d = 3.0/2.0*dih4*dih2*dir;
      }
  
      if ((~compare2) & 0x3) {
	a2 = -1.0/15.0*dir + dih *(8.0/5.0-(4.0/3.0*u2)+u2*u - 3.0/10.0*u4 + 1.0/30.0*u4*u);
	b2 = -1.0/15.0*dir2*dir + dih2*dih*(8.0/3.0 - 3.0 * u + 6.0/5.0*u2 -1.0/6.0*u2*u);
	c2 = -1.0/5.0*dir4*dir + 3.0 * dih4 * dir + dih4*dih*(-12.0/5.0 + 1.0/2.0*u);
	d2 = -1.0*dir4*dir2*dir + 3.0 *dih4 * dir2 * dir - 1.0/2.0*dih4*dih2*dir;
      }

      a1 = ((select2 & a) | andnot(select2, a2));
      b1 = ((select2 & b) | andnot(select2, b2));
      c1 = ((select2 & c) | andnot(select2, c2));
      d1 = ((select2 & d) | andnot(select2, d2));
    }
    if ((~compare) & 0x3) {
      a2 = dir;
      b2 = a2* a2*a2;
      c2 = 3.0*b2*a2*a2;
      d2 = 5.0*c2*a2*a2;
    }

    a = ((select & a1) | andnot(select, a2));
    b = ((select & b1) | andnot(select, b2));
    c = ((select & c1) | andnot(select, c2));
    d = ((select & d1) | andnot(select, d2));
    SSEFloat mxx(m.xx);SSEFloat mxy(m.xy);SSEFloat myy(m.yy); SSEFloat myz(m.yz);SSEFloat mzz(m.zz);SSEFloat mxz(m.xz);
    SSEFloat qirx = mxx*rx + mxy*ry + mxz*rz;
    SSEFloat qiry = mxy*rx + myy*ry + myz*rz;
    SSEFloat qirz = mxz*rx + myz*ry + mzz*rz;
    SSEFloat qir = 0.5 * (qirx*rx+qiry*ry+qirz*rz);
    SSEFloat qir3 = m.totalMass * b  + d*qir - c *tr;

    SSEFloat temppotential(input[i]->potential,input[i+1]->potential,input[i+2]->potential,input[i+3]->potential);

    potential = temppotential - (m.totalMass*a + c*qir - b*tr);
    SSEFloat tempaccX(input[i]->treeAcceleration.x,input[i+1]->treeAcceleration.x,input[i+2]->treeAcceleration.x,input[i+3]->treeAcceleration.x);
    accX = tempaccX -( (qir3 *rx) -(c *qirx));

    SSEFloat tempaccY(input[i]->treeAcceleration.y,input[i+1]->treeAcceleration.y,input[i+2]->treeAcceleration.y,input[i+3]->treeAcceleration.y);
    accY = tempaccY -( (qir3 *ry) -(c *qiry));
    SSEFloat tempaccZ(input[i]->treeAcceleration.z,input[i+1]->treeAcceleration.z,input[i+2]->treeAcceleration.z,input[i+3]->treeAcceleration.z);
    accZ = tempaccZ - (qir3 *rz -c *qirz);
    SSEFloat tempmass(input[i]->mass,input[i+1]->mass,input[i+2]->mass,input[i+3]->mass);
    SSEFloat idt2 = (m.totalMass + tempmass) * b;
    SSEFloat dtGrav(input[i]->dtGrav, input[i+1]->dtGrav,input[i+2]->dtGrav, input[i+3]->dtGrav);
    idt2 = max(idt2,dtGrav);

    float h[4];
    float *p=h;
    storeu(p, potential);
    *(&input[i]->potential) = *p;
   
    p++;
    *(&input[i+1]->potential) = *p;
    p++;
    *(&input[i+2]->potential) = *p;
    p++;
    *(&input[i+3]->potential) = *p;
      
    float d[4];
    p=d;
    storeu(p, accX);
    *(&input[i]->treeAcceleration.x) = *p;
    p++;
    *(&input[i+1]->treeAcceleration.x) = *p;
    p++;
    *(&input[i+2]->treeAcceleration.x) = *p;
    p++;
    *(&input[i+3]->treeAcceleration.x) = *p;
    float l[4];
    p=l;
    storeu(p, accY);
    *(&input[i]->treeAcceleration.y) = *p;
    p++;
    *(&input[i+1]->treeAcceleration.y) = *p;
    p++;
    *(&input[i+2]->treeAcceleration.y) = *p;
    p++;
    *(&input[i+3]->treeAcceleration.y) = *p;
    float x[4];
    p=x;
    storeu(p, accZ);
    *(&input[i]->treeAcceleration.z) = *p;
    p++;
    *(&input[i+1]->treeAcceleration.z) = *p;
    p++;
    *(&input[i+2]->treeAcceleration.z) = *p;
    p++;
    *(&input[i+3]->treeAcceleration.z) = *p;
    float y[4];
    p=y;
    storeu(p, idt2);
    *(&input[i]->dtGrav) = *p;
    p++;
    *(&input[i+1]->dtGrav) = *p;
    p++;
    *(&input[i+2]->dtGrav) = *p;
    p++;
    *(&input[i+3]->dtGrav) = *p;
  }
  return computed;
}
#elif  defined(__SSE2__) && !defined(COSMO_FLOAT) && defined(HEXADECAPOLE)
      
inline
int nodeBucketForce(Tree::GenericTreeNode *node, // source of force
		    Tree::GenericTreeNode *req,  // bucket descriptor
		    GravityParticle *particles, // particles in bucket
		    Vector3D<double> offset,    // offset if a periodic replica
		    int activeRung)         // rung (and above) at which to
     // calculate forces
{
    
  int computed = 0;
  SSEDouble rx,ry,rz,select2;
  SSEDouble rsq,twoh,a,b,a1,a2,b1,b2;
  MultipoleMoments m1 = node->moments;
  SSEDouble mtotalMass(m1.totalMass);
  GravityParticle dummyPart = particles[0];
  Vector3D<double> cm(m1.cm + offset);
  GravityParticle **input = (GravityParticle**)alloca((req->lastParticle - req->firstParticle + 2) * sizeof(GravityParticle*));
  int activeParts = 0;
  for (int j=req->firstParticle; j <= req->lastParticle; ++j) {
    if (particles[j].rung >= activeRung) input[activeParts++] = &particles[j];
  }
  input[activeParts] = &dummyPart;

  /*
  int rem = 3-(activeParts % 4);
  for(int cnt =0;cnt<rem;cnt++)
    input[++activeParts] = &dummyPart;
  */

  if(openSoftening(node, req, offset)) {
    for (int i=0; i<activeParts; i+=2) {
      computed+=2;
      SSEDouble accX(0.0),accY(0.0),accZ(0.0),potential(0.0);
      SSEDouble posX(input[i]->position.x, input[i+1]->position.x);
      SSEDouble posY(input[i]->position.y, input[i+1]->position.y);
      SSEDouble posZ(input[i]->position.z, input[i+1]->position.z);
      SSEDouble mcmx(m1.cm.x); SSEDouble mcmy(m1.cm.y); SSEDouble mcmz(m1.cm.z);
      SSEDouble offsetx(offset.x); SSEDouble offsety(offset.y); SSEDouble offsetz(offset.z);
      rx = offsetx - posX + mcmx;
      ry = offsety - posY + mcmy;
      rz = offsetz - posZ + mcmz;
      rsq= rx*rx+(ry*ry+rz*rz);
      SSEDouble k(m1.soft);SSEDouble t(input[i]->soft,input[i+1]->soft);
      twoh = k + t ;
      SSEDouble select0 = rsq > 0.0;
      int compare0 = movemask(select0);

      if (compare0) {
	SSEDouble r = sqrt(rsq);  
	SSEDouble select = r < twoh;
	int compare = movemask(select);
        if(compare) {
	  SSEDouble dih = 2.0/twoh;
	  SSEDouble u = r * dih;
	  select2 = u < 1.0;
	  int compare2 = movemask(select2);
	  SSEDouble u2 = u *u;
	  SSEDouble u4 = u2 * u2;
	  SSEDouble dih2 =dih *dih;
	  SSEDouble dih4 =dih2 *dih2;
	  if (compare2) {
	    a = dih * ( (7.0/5.0) - ((2.0/3.0)* u2) + ((3.0/10.0)*u4) - ((1.0/10.0)* u4 *u));
	    b = dih2*dih*(4.0/3.0-(6.0/5.0*u2) + (1.0/2.0*u2*u));
	  }

	  if ((~compare2) & 0x3) {
	    SSEDouble dir = 1.0 / r;
	    SSEDouble  dir2 = dir*dir;
	    a2 = -1.0/15.0*dir + dih *(8.0/5.0-(4.0/3.0*u2)+u2*u - 3.0/10.0*u4 + 1.0/30.0*u4*u);
	    b2 = -1.0/15.0*dir2*dir + dih2*dih*(8.0/3.0 - 3.0 * u + 6.0/5.0*u2 -1.0/6.0*u2*u);
	  }
	  a1 = ((select2 & a) | andnot(select2, a2));
	  b1 = ((select2 & b) | andnot(select2, b2));
	}
	if ((~compare) & 0x3) {
	  a2 = 1.0 /r ;
	  b2 = a2* a2*a2;
	}
        a = ((select & a1) | andnot(select, a2));
        b = ((select & b1) | andnot(select, b2));

        SSEDouble temppotential(input[i]->potential,input[i+1]->potential);
        potential = temppotential - (m1.totalMass * a);
        SSEDouble tempaccX(input[i]->treeAcceleration.x,input[i+1]->treeAcceleration.x);
        accX =  tempaccX + rx * (m1.totalMass *b);
        SSEDouble tempaccY(input[i]->treeAcceleration.y,input[i+1]->treeAcceleration.y);
        accY =  tempaccY + ry * (m1.totalMass *b) ;
        SSEDouble tempaccZ(input[i]->treeAcceleration.z,input[i+1]->treeAcceleration.z);
        accZ =   tempaccZ + rz * (m1.totalMass *b) ;
	SSEDouble tempmass(input[i]->mass,input[i+1]->mass);
	SSEDouble idt2 = (m1.totalMass + tempmass) * b;
	SSEDouble dtGrav(input[i]->dtGrav, input[i+1]->dtGrav);
	idt2 = max(idt2,dtGrav);
	storel(&input[i]->potential, potential);
	storeh(&input[i+1]->potential, potential);
	storel(&input[i]->treeAcceleration.x, accX);
	storeh(&input[i+1]->treeAcceleration.x, accX);
	storel(&input[i]->treeAcceleration.y, accY);
	storeh(&input[i+1]->treeAcceleration.y, accY);
	storel(&input[i]->treeAcceleration.z, accZ);
	storeh(&input[i+1]->treeAcceleration.z, accZ);
	storel(&input[i]->dtGrav, idt2);
	storeh(&input[i+1]->dtGrav, idt2);
      }
    }
    return computed;
  }

  MOMR *m = &m1.mom;
  for (int i=0; i<activeParts; i+=2) {
    computed+=2;
    SSEDouble accX(0.0),accY(0.0),accZ(0.0),potential(0.0);
    SSEDouble posX(input[i]->position.x, input[i+1]->position.x);
    SSEDouble posY(input[i]->position.y, input[i+1]->position.y);
    SSEDouble posZ(input[i]->position.z, input[i+1]->position.z);
    SSEDouble cmx(cm.x); SSEDouble cmy(cm.y); SSEDouble cmz(cm.z);
    rx = posX - cmx;
    ry = posY - cmy;
    rz = posZ - cmz;
    rsq= rx*rx+(ry*ry+rz*rz);
    SSEDouble dir = 1.0/sqrt(rsq);
    SSEDouble dir0 = dir;

    rx = -1.0 * rx;
    ry = -1.0 * ry;
    rz = -1.0 * rz; 
    const SSEDouble onethird = 1.0/3.0;
    SSEDouble xx,xy,xz,yy,yz,zz;
    SSEDouble xxx,xxy,xxz,xyy,yyy,yyz,xyz;
    SSEDouble tx,ty,tz,dir2,g2,g3,g4;
    dir = -1.0 * dir;
    dir2 = dir*dir;
    g2 = 3.0*dir*dir2*dir2;
    g3 = -5.0*g2*dir2;
    g4 = -7.0*g3*dir2;
    /**** Calculate the funky distance terms.  */
    xx = 0.5*rx*rx;
    xy = rx*ry;
    xz = rx*rz;
    yy = 0.5*ry*ry;
    yz = ry*rz;
    zz = 0.5*rz*rz;
    xxx = rx*(onethird*xx - zz);
    xxz = rz*(xx - onethird*zz);
    yyy = ry*(onethird*yy - zz);
    yyz = rz*(yy - onethird*zz);
    xx = xx - zz;
    yy = yy - zz;
    xxy = ry*xx;
    xyy = rx*yy;
    xyz = xy*rz;

    /* *  ** Now calculate the interaction up to Hexadecapole order.      */

    tx = g4*(m->xxxx*xxx + m->xyyy*yyy + m->xxxy*xxy + m->xxxz*xxz + m->xxyy*xyy + m->xxyz*xyz + m->xyyz*yyz);
    ty = g4*(m->xyyy*xyy + m->xxxy*xxx + m->yyyy*yyy + m->yyyz*yyz + m->xxyy*xxy + m->xxyz*xxz + m->xyyz*xyz);
    tz = g4*(-m->xxxx*xxz - (m->xyyy + m->xxxy)*xyz - m->yyyy*yyz + m->xxxz*xxx + m->yyyz*yyy - m->xxyy*(xxz + yyz) + m->xxyz*xxy + m->xyyz*xyy);
    g4 = 0.25*(tx*rx + ty*ry + tz*rz);
    xxx = g3*(m->xxx*xx + m->xyy*yy + m->xxy*xy + m->xxz*xz + m->xyz*yz);
    xxy = g3*(m->xyy*xy + m->xxy*xx + m->yyy*yy + m->yyz*yz + m->xyz*xz);
    xxz = g3*(-(m->xxx + m->xyy)*xz - (m->xxy + m->yyy)*yz + m->xxz*xx + m->yyz*yy + m->xyz*xy);
    g3 = onethird*(xxx*rx + xxy*ry + xxz*rz);
    xx = g2*(m->xx*rx + m->xy*ry + m->xz*rz);
    xy = g2*(m->yy*ry + m->xy*rx + m->yz*rz);
    xz = g2*(-(m->xx + m->yy)*rz + m->xz*rx + m->yz*ry);
    g2 = 0.5*(xx*rx + xy*ry + xz*rz);
    dir = m->m * dir;
    dir2 =  - 1.0 * (dir + 5.0*g2 + 7.0*g3 + 9.0*g4) * dir2;

    SSEDouble temppotential(input[i]->potential,input[i+1]->potential);

    potential = temppotential + dir + g2 + g3 + g4;
    SSEDouble tempaccX(input[i]->treeAcceleration.x,input[i+1]->treeAcceleration.x);
    accX =  tempaccX + xx + xxx + tx + rx*dir2;
    SSEDouble tempaccY(input[i]->treeAcceleration.y,input[i+1]->treeAcceleration.y);
    accY =  tempaccY + xy + xxy + ty + ry*dir2;
    SSEDouble tempaccZ(input[i]->treeAcceleration.z,input[i+1]->treeAcceleration.z); 
    accZ =   tempaccZ + xz + xxz + tz + rz*dir2;
   
    SSEDouble tempmass(input[i]->mass,input[i+1]->mass);
    SSEDouble idt2 = (m1.totalMass + tempmass) * dir0 * dir0 * dir0;
    SSEDouble dtGrav(input[i]->dtGrav, input[i+1]->dtGrav);  
    idt2 = max(idt2,dtGrav);
    storel(&input[i]->potential, potential);
    storeh(&input[i+1]->potential, potential);
    storel(&input[i]->treeAcceleration.x, accX);
    storeh(&input[i+1]->treeAcceleration.x, accX);
    storel(&input[i]->treeAcceleration.y, accY);
    storeh(&input[i+1]->treeAcceleration.y, accY);
    storel(&input[i]->treeAcceleration.z, accZ);
    storeh(&input[i+1]->treeAcceleration.z, accZ);
    storel(&input[i]->dtGrav, idt2);
    storeh(&input[i+1]->dtGrav, idt2);
  }
  return computed;
}

#elif  defined(__SSE2__) && defined(COSMO_FLOAT) && defined(HEXADECAPOLE)

inline
int nodeBucketForce(Tree::GenericTreeNode *node, // source of force
		    Tree::GenericTreeNode *req,  // bucket descriptor
		    GravityParticle *particles, // particles in bucket
		    Vector3D<double> offset,    // offset if a periodic replica
		    int activeRung)         // rung (and above) at which to
     // calculate forces
{
  int computed = 0;
  SSEFloat rx,ry,rz,select2;
  SSEFloat rsq,twoh,a,b,a1,a2,b1,b2;
  MultipoleMoments m1 = node->moments;
  SSEFloat mtotalMass(m1.totalMass);
  GravityParticle dummyPart = particles[0];
  Vector3D<double> cm(m1.cm + offset);
                     
  GravityParticle **input = (GravityParticle**)alloca((req->lastParticle - req->firstParticle + 5) * sizeof(GravityParticle*));
                      
  int activeParts = 0;
  for (int j=req->firstParticle; j <= req->lastParticle; ++j) {
                        
    if (particles[j].rung >= activeRung) input[activeParts++] = &particles[j];
  }
  input[activeParts] = &dummyPart;
                         
  int rem = 3-(activeParts % 4);

  for(int cnt =0;cnt<rem;cnt++)
    input[++activeParts] = &dummyPart;

  if(openSoftening(node, req, offset)) {

    for (int i=0; i<activeParts; i+=4) {
      computed+=4;

      SSEFloat accX(0.0),accY(0.0),accZ(0.0),potential(0.0);

      SSEFloat posX(input[i]->position.x, input[i+1]->position.x,input[i+2]->position.x, input[i+3]->position.x);
      SSEFloat posY(input[i]->position.y, input[i+1]->position.y,input[i+2]->position.y, input[i+3]->position.y);
      SSEFloat posZ(input[i]->position.z, input[i+1]->position.z,input[i+2]->position.z, input[i+3]->position.z);
      SSEFloat mcmx(m1.cm.x); SSEFloat mcmy(m1.cm.y); SSEFloat mcmz(m1.cm.z);
      SSEFloat offsetx(offset.x); SSEFloat offsety(offset.y); SSEFloat offsetz(offset.z);
      rx = offsetx - posX + mcmx;
      ry = offsety - posY + mcmy;
      rz = offsetz - posZ + mcmz;
      rsq= rx*rx+(ry*ry+rz*rz);
      SSEFloat k(m1.soft);SSEFloat t(input[i]->soft,input[i+1]->soft,input[i+2]->soft,input[i+3]->soft);
      twoh = k + t ;
      SSEFloat select0 = rsq > 0.0;
      int compare0 = movemask(select0);

      if (compare0) {
	SSEFloat r = sqrt(rsq);
	SSEFloat select = r < twoh;
	int compare = movemask(select);

        if(compare) {
	  SSEFloat dih = 2.0/twoh;
	  SSEFloat u = r * dih;
	  select2 = u < 1.0;
	  int compare2 = movemask(select2);
	  SSEFloat u2 = u *u;
	  SSEFloat u4 = u2 * u2;
	  SSEFloat dih2 =dih *dih;
	  SSEFloat dih4 =dih2 *dih2;

	  if (compare2) {
	    a = dih * ( (7.0/5.0) - ((2.0/3.0)* u2) + ((3.0/10.0)*u4) - ((1.0/10.0)* u4 *u));
	    b = dih2*dih*(4.0/3.0-(6.0/5.0*u2) + (1.0/2.0*u2*u));
	  }

	  if ((~compare2) & 0x3) {
	    SSEFloat dir = 1.0 / r;
	    SSEFloat  dir2 = dir*dir;
	    a2 = -1.0/15.0*dir + dih *(8.0/5.0-(4.0/3.0*u2)+u2*u - 3.0/10.0*u4 + 1.0/30.0*u4*u);
	    b2 = -1.0/15.0*dir2*dir + dih2*dih*(8.0/3.0 - 3.0 * u + 6.0/5.0*u2 -1.0/6.0*u2*u);
	  }
	  a1 = ((select2 & a) | andnot(select2, a2));
	  b1 = ((select2 & b) | andnot(select2, b2));
	}
	if ((~compare) & 0x3) {
	  a2 = 1.0 /r ;
	  b2 = a2* a2*a2;
	}
        a = ((select & a1) | andnot(select, a2));
        b = ((select & b1) | andnot(select, b2));
        SSEFloat temppotential(input[i]->potential,input[i+1]->potential,input[i+2]->potential,input[i+3]->potential);
        potential = temppotential - (m1.totalMass * a);
        SSEFloat tempaccX(input[i]->treeAcceleration.x,input[i+1]->treeAcceleration.x,input[i+2]->treeAcceleration.x,input[i+3]->treeAcceleration.x);
        accX =  tempaccX + rx * (m1.totalMass *b);
        SSEFloat tempaccY(input[i]->treeAcceleration.y,input[i+1]->treeAcceleration.y,input[i+2]->treeAcceleration.y,input[i+3]->treeAcceleration.y);
        accY =  tempaccY + ry * (m1.totalMass *b) ;
        SSEFloat tempaccZ(input[i]->treeAcceleration.z,input[i+1]->treeAcceleration.z,input[i+2]->treeAcceleration.z,input[i+3]->treeAcceleration.z);
        accZ =   tempaccZ + rz * (m1.totalMass *b) ;

	SSEFloat tempmass(input[i]->mass,input[i+1]->mass,input[i+2]->mass,input[i+3]->mass);
	SSEFloat idt2 = (m1.totalMass + tempmass) * b;
	SSEFloat dtGrav(input[i]->dtGrav, input[i+1]->dtGrav,input[i+2]->dtGrav, input[i+3]->dtGrav);
	idt2 = max(idt2,dtGrav);
	float h[4];
	float *p=h;
	storeu(p, potential);
	*(&input[i]->potential) = *p;
	p++;
	*(&input[i+1]->potential) = *p;
        p++;
	*(&input[i+2]->potential) = *p;
	p++;
	*(&input[i+3]->potential) = *p;
	float d[4];
	p=d;
	storeu(p, accX);
	*(&input[i]->treeAcceleration.x) = *p;
	p++;
	*(&input[i+1]->treeAcceleration.x) = *p;
	p++;
	*(&input[i+2]->treeAcceleration.x) = *p;
	p++;
	*(&input[i+3]->treeAcceleration.x) = *p;
	float l[4];
	p=l;
	storeu(p, accY);
	*(&input[i]->treeAcceleration.y) = *p;
	p++;
	*(&input[i+1]->treeAcceleration.y) = *p;
	p++;
	*(&input[i+2]->treeAcceleration.y) = *p;
	p++;
	*(&input[i+3]->treeAcceleration.y) = *p;
	float x[4];
	p=x;
	storeu(p, accZ);
	*(&input[i]->treeAcceleration.z) = *p;
	p++;
	*(&input[i+1]->treeAcceleration.z) = *p;
	p++;
	*(&input[i+2]->treeAcceleration.z) = *p;
	p++;
	*(&input[i+3]->treeAcceleration.z) = *p;
	float y[4];
	p=y;
	storeu(p, idt2);
	*(&input[i]->dtGrav) = *p;
	p++;
        *(&input[i+1]->dtGrav) = *p;
	p++;
        *(&input[i+2]->dtGrav) = *p;
	p++;
        *(&input[i+3]->dtGrav) = *p;
      }
    }

    return computed;
  }

  MOMR *m = &m1.mom;
  for (int i=0; i<activeParts; i+=4) {
    computed+=4;

    SSEFloat accX(0.0),accY(0.0),accZ(0.0),potential(0.0);

    SSEFloat posX(input[i]->position.x, input[i+1]->position.x,input[i+2]->position.x, input[i+3]->position.x);
    SSEFloat posY(input[i]->position.y, input[i+1]->position.y,input[i+2]->position.y, input[i+3]->position.y);
    SSEFloat posZ(input[i]->position.z, input[i+1]->position.z,input[i+2]->position.z, input[i+3]->position.z);
    SSEFloat cmx(cm.x); SSEFloat cmy(cm.y); SSEFloat cmz(cm.z);

    rx = posX - cmx;
    ry = posY - cmy;
    rz = posZ - cmz;
    rsq= rx*rx+(ry*ry+rz*rz);
    SSEFloat dir = 1.0/sqrt(rsq);
    SSEFloat dir0 = dir;
    rx = -1.0 * rx;
    ry = -1.0 * ry;
    rz = -1.0 * rz;

    const SSEFloat onethird = 1.0/3.0;
    SSEFloat xx,xy,xz,yy,yz,zz;
    SSEFloat xxx,xxy,xxz,xyy,yyy,yyz,xyz;
    SSEFloat tx,ty,tz,dir2,g2,g3,g4;

    dir = -1.0 * dir;
    dir2 = dir*dir;
    g2 = 3.0*dir*dir2*dir2;
    g3 = -5.0*g2*dir2;
    g4 = -7.0*g3*dir2;

    /**** Calculate the funky distance terms.  */
    xx = 0.5*rx*rx;
    xy = rx*ry;
    xz = rx*rz;
    yy = 0.5*ry*ry;
    yz = ry*rz;
    zz = 0.5*rz*rz;
    xxx = rx*(onethird*xx - zz);
    xxz = rz*(xx - onethird*zz);
    yyy = ry*(onethird*yy - zz);
    yyz = rz*(yy - onethird*zz);
    xx = xx - zz;
    yy = yy - zz;
    xxy = ry*xx;
    xyy = rx*yy;
    xyz = xy*rz;

    /* *  ** Now calculate the interaction up to Hexadecapole order.      */

    tx = g4*(m->xxxx*xxx + m->xyyy*yyy + m->xxxy*xxy + m->xxxz*xxz + m->xxyy*xyy + m->xxyz*xyz + m->xyyz*yyz);
    ty = g4*(m->xyyy*xyy + m->xxxy*xxx + m->yyyy*yyy + m->yyyz*yyz + m->xxyy*xxy + m->xxyz*xxz + m->xyyz*xyz);
    tz = g4*(-m->xxxx*xxz - (m->xyyy + m->xxxy)*xyz - m->yyyy*yyz + m->xxxz*xxx + m->yyyz*yyy - m->xxyy*(xxz + yyz) + m->xxyz*xxy + m->xyyz*xyy);
    g4 = 0.25*(tx*rx + ty*ry + tz*rz);
    xxx = g3*(m->xxx*xx + m->xyy*yy + m->xxy*xy + m->xxz*xz + m->xyz*yz);
    xxy = g3*(m->xyy*xy + m->xxy*xx + m->yyy*yy + m->yyz*yz + m->xyz*xz);
    xxz = g3*(-(m->xxx + m->xyy)*xz - (m->xxy + m->yyy)*yz + m->xxz*xx + m->yyz*yy + m->xyz*xy);
    g3 = onethird*(xxx*rx + xxy*ry + xxz*rz);
    xx = g2*(m->xx*rx + m->xy*ry + m->xz*rz);
    xy = g2*(m->yy*ry + m->xy*rx + m->yz*rz);
    xz = g2*(-(m->xx + m->yy)*rz + m->xz*rx + m->yz*ry);
    g2 = 0.5*(xx*rx + xy*ry + xz*rz);
    dir = m->m * dir;
    dir2 =  - 1.0 * (dir + 5.0*g2 + 7.0*g3 + 9.0*g4) * dir2;

    SSEFloat temppotential(input[i]->potential,input[i+1]->potential,input[i+2]->potential,input[i+3]->potential);
    potential = temppotential + dir + g2 + g3 + g4;
    SSEFloat tempaccX(input[i]->treeAcceleration.x,input[i+1]->treeAcceleration.x,input[i+2]->treeAcceleration.x,input[i+3]->treeAcceleration.x);
    accX =  tempaccX + xx + xxx + tx + rx*dir2;
    SSEFloat tempaccY(input[i]->treeAcceleration.y,input[i+1]->treeAcceleration.y,input[i+2]->treeAcceleration.y,input[i+3]->treeAcceleration.y);
    accY =  tempaccY + xy + xxy + ty + ry*dir2;
    SSEFloat tempaccZ(input[i]->treeAcceleration.z,input[i+1]->treeAcceleration.z,input[i+2]->treeAcceleration.z,input[i+3]->treeAcceleration.z);
    accZ =   tempaccZ + xz + xxz + tz + rz*dir2;
    SSEFloat tempmass(input[i]->mass,input[i+1]->mass,input[i+2]->mass,input[i+3]->mass);
    SSEFloat idt2 = (m1.totalMass + tempmass) * dir0 * dir0 * dir0;
    SSEFloat dtGrav(input[i]->dtGrav, input[i+1]->dtGrav,input[i+2]->dtGrav, input[i+3]->dtGrav);
    idt2 = max(idt2,dtGrav);
    float h[4];
    float *p=h;
    storeu(p, potential);
    *(&input[i]->potential) = *p;
    p++;
    *(&input[i+1]->potential) = *p;
    p++;
    *(&input[i+2]->potential) = *p;
    p++;
    *(&input[i+3]->potential) = *p;
    float d[4];
    p=d;
    storeu(p, accX);
    *(&input[i]->treeAcceleration.x) = *p;
    p++;
    *(&input[i+1]->treeAcceleration.x) = *p;
    p++;
    *(&input[i+2]->treeAcceleration.x) = *p;
    p++;
    *(&input[i+3]->treeAcceleration.x) = *p;
    float l[4];
    p=l;
    storeu(p, accY);
    *(&input[i]->treeAcceleration.y) = *p;
    p++;
    *(&input[i+1]->treeAcceleration.y) = *p;
    p++;
    *(&input[i+2]->treeAcceleration.y) = *p;
    p++;
    *(&input[i+3]->treeAcceleration.y) = *p;
    float x[4];
    p=x;
    storeu(p, accZ);
    *(&input[i]->treeAcceleration.z) = *p;
    p++;
    *(&input[i+1]->treeAcceleration.z) = *p;
    p++;
    *(&input[i+2]->treeAcceleration.z) = *p;
    p++;
    *(&input[i+3]->treeAcceleration.z) = *p;
    float y[4];
    p=y;
    storeu(p, idt2);
    *(&input[i]->dtGrav) = *p;
    p++;
    *(&input[i+1]->dtGrav) = *p;
    p++;
    *(&input[i+2]->dtGrav) = *p;
    p++;
    *(&input[i+3]->dtGrav) = *p;
  }

  return computed;
}
   
#endif




// Return true if the node's opening radius intersects the
// boundingBox, i.e. the node needs to be opened.

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

// return -1 if there is an intersection
//    0 if no intersection
//    +1 if completely contained

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
        return 1;
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
            return 1;
            }
#else
        return 0;
#endif
    }
}

#endif
