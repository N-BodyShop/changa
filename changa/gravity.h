#ifndef __GRAVITY_H__
#define __GRAVITY_H__

#include "TreeNode.h"
#include "GenericTreeNode.h"
#include "Space.h"

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
#ifdef COSMO_FLOAT
inline
void SPLINEQ(float invr,float r2,float twoh,float& a,float& b,
	     float& c,float& d)
{
	float u,dih,dir=(invr);
	if ((r2) < (twoh)*(twoh)) {
		dih = 2.0f/(twoh);
		u = dih/dir;
		if (u < 1.0f) {
			a = dih*(7.0f/5.0f - 2.0f/3.0f*u*u + 3.0f/10.0f*u*u*u*u
					 - 1.0f/10.0f*u*u*u*u*u);
			b = dih*dih*dih*(4.0f/3.0f - 6.0f/5.0f*u*u + 1.0f/2.0f*u*u*u);
		    c = dih*dih*dih*dih*dih*(12.0f/5.0f - 3.0f/2.0f*u);
			d = 3.0f/2.0f*dih*dih*dih*dih*dih*dih*dir;
			}
		else {
			a = -1.0f/15.0f*dir + dih*(8.0f/5.0f - 4.0f/3.0f*u*u + u*u*u
			              - 3.0f/10.0f*u*u*u*u + 1.0f/30.0f*u*u*u*u*u);
			b = -1.0f/15.0f*dir*dir*dir + dih*dih*dih*(8.0f/3.0f - 3.0f*u + 6.0f/5.0f*u*u - 1.0f/6.0f*u*u*u);
			c = -1.0f/5.0f*dir*dir*dir*dir*dir + 3.0f*dih*dih*dih*dih*dir
				+ dih*dih*dih*dih*dih*(-12.0f/5.0f + 1.0f/2.0f*u);
			d = -dir*dir*dir*dir*dir*dir*dir
				+ 3.0f*dih*dih*dih*dih*dir*dir*dir
					- 1.0f/2.0f*dih*dih*dih*dih*dih*dih*dir;
			}
		}
	else {
		a = dir;
		b = a*a*a;
		c = 3.0f*b*a*a;
		d = 5.0f*c*a*a;
		}
	}
#else
inline
void SPLINEQ(double invr,double r2,double twoh,double& a,double& b,
	     double& c,double& d)
{
	double u,dih,dir=(invr);
	if ((r2) < (twoh)*(twoh)) {
		dih = 2.0/(twoh);
		u = dih/dir;
		if (u < 1.0) {
			a = dih*(7.0/5.0 - 2.0/3.0*u*u + 3.0/10.0*u*u*u*u
					 - 1.0/10.0*u*u*u*u*u);
			b = dih*dih*dih*(4.0/3.0 - 6.0/5.0*u*u + 1.0/2.0*u*u*u);
		    c = dih*dih*dih*dih*dih*(12.0/5.0 - 3.0/2.0*u);
			d = 3.0/2.0*dih*dih*dih*dih*dih*dih*dir;
			}
		else {
			a = -1.0/15.0*dir + dih*(8.0/5.0 - 4.0/3.0*u*u + u*u*u
			              - 3.0/10.0*u*u*u*u + 1.0/30.0*u*u*u*u*u);
			b = -1.0/15.0*dir*dir*dir + dih*dih*dih*(8.0/3.0 - 3.0*u + 6.0/5.0*u*u - 1.0/6.0*u*u*u);
			c = -1.0/5.0*dir*dir*dir*dir*dir + 3.0*dih*dih*dih*dih*dir
				+ dih*dih*dih*dih*dih*(-12.0/5.0 + 1.0/2.0*u);
			d = -dir*dir*dir*dir*dir*dir*dir
				+ 3.0*dih*dih*dih*dih*dir*dir*dir
					- 1.0/2.0*dih*dih*dih*dih*dih*dih*dir;
			}
		}
	else {
		a = dir;
		b = a*a*a;
		c = 3.0*b*a*a;
		d = 5.0*c*a*a;
		}
	}
#endif

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
			b = -1.0/15.0*dir*dir*dir + dih*dih*dih*(8.0/3.0 - 3.0*u + 6.0/5.0*u*u - 1.0/6.0*u*u*u);
			}
		}
	else {
		a = dir;
		b = a*a*a;
		}
	}


inline
void SPLINE(double r2, double twoh, double &a, double &b)
{
	double r, u,dih,dir;
	r = sqrt(r2);
	if (r < (twoh)) {
		dih = 2.0/(twoh);
		u = r*dih;
		if (u < 1.0) {
			a = dih*(7.0/5.0 - 2.0/3.0*u*u + 3.0/10.0*u*u*u*u
					 - 1.0/10.0*u*u*u*u*u);
			b = dih*dih*dih*(4.0/3.0 - 6.0/5.0*u*u + 1.0/2.0*u*u*u);
			}
		else {
			dir = 1.0/r;
			a = -1.0/15.0*dir + dih*(8.0/5.0 - 4.0/3.0*u*u + u*u*u
			              - 3.0/10.0*u*u*u*u + 1.0/30.0*u*u*u*u*u);
			b = -1.0/15.0*dir*dir*dir + dih*dih*dih*(8.0/3.0 - 3.0*u + 6.0/5.0*u*u - 1.0/6.0*u*u*u);
			}
		}
	else {
		a = 1.0/r;
		b = a*a*a;
		}
	}

//
// Return true if the soften nodes overlap, i.e. the forces involve softening
// @param node
// @param myNode
// @param offset Periodic offset to applied to node
//
inline
int openSoftening(Tree::GenericTreeNode *node, Tree::GenericTreeNode *myNode,
          Vector3D<double> offset)
{
    Sphere<double> s(node->moments.cm + offset, 2.0*node->moments.soft);
    Sphere<double> myS(myNode->moments.cm, 2.0*myNode->moments.soft);
    return Space::intersect(myS, s);
    }

#ifdef CMK_VERSION_BLUEGENE
static int forProgress = 0;
#endif

inline int partBucketForce(ExternalGravityParticle *part, Tree::GenericTreeNode *req, GravityParticle *particles, Vector3D<double> offset, int activeRung) {
  int computed = 0;
#ifndef BENCHMARK_NO_WORK
  Vector3D<double> r;
  double rsq;
  double twoh, a, b;

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
    double idt2 = (particles[j].mass + part->mass)*b; // (timescale)^-2
                              // of interaction
        particles[j].treeAcceleration += r * (b * part->mass);
        particles[j].potential -= part->mass * a;
    if(idt2 > particles[j].dtGrav)
        particles[j].dtGrav = idt2;
      }
    }
  }
#endif
  return computed;
}

#if defined(__SSE2__) && !defined(HEXADECAPOLE)
#include <emmintrin.h>

inline
int nodeBucketForce(Tree::GenericTreeNode *node, // source of force
                    Tree::GenericTreeNode *req,  // bucket descriptor
                    GravityParticle *particles, // particles in bucket
                    Vector3D<double> offset,    // offset if a periodic replica
                    int activeRung)         // rung (and above) at which to
// calculate forces
{
  int computed = 0;
#ifndef BENCHMARK_NO_WORK
  __m128d rx, ry, rz;
  __m128d rsq;
  __m128d twoh, a, b, c, d, a1, a2, b1, b2, c1, c2, d1, d2;
  MultipoleMoments m = node->moments;
  __m128d posX, posY, posZ, accX, accY, accZ, potential;

  GravityParticle dummyPart = particles[0];
  Vector3D<double> cm(m.cm + offset);
  __m128d tr = _mm_set1_pd(0.5*(m.xx + m.yy + m.zz));

  GravityParticle **input = (GravityParticle**)alloca((req->lastParticle - req->firstParticle + 2) * sizeof(GravityParticle*));

  //int printbool = 0;
  //for (int j=req->firstParticle; j <= req->lastParticle; ++j) {
  //  if (particles[j].iOrder < 10) printbool = 1;
  //}

  int activeParts = 0;
  for (int j=req->firstParticle; j <= req->lastParticle; ++j) {
    //if (printbool) CkPrintf("parts: %d, iorder: %d\n",req->lastParticle-req->firstParticle+1,particles[j].iOrder);
    if (particles[j].rung >= activeRung) input[activeParts++] = &particles[j];
  }
  input[activeParts] = &dummyPart;

  for (int i=0; i<activeParts; i+=2) {
      computed+=2;
      //if (printbool) CkPrintf("values: %d-%d %d-%d (%d-%d)\n",input[i]-particles,input[i]->iOrder,input[i+1]-particles,input[i+1]->iOrder,input[i+2]-particles,input[i+2]->iOrder);
      accX = accY = accZ = potential = _mm_setzero_pd();
      posX = _mm_setr_pd(input[i]->position.x, input[i+1]->position.x);
      posY = _mm_setr_pd(input[i]->position.y, input[i+1]->position.y);
      posZ = _mm_setr_pd(input[i]->position.z, input[i+1]->position.z);
      rx = _mm_sub_pd(posX, _mm_set1_pd(cm.x));
      ry = _mm_sub_pd(posY, _mm_set1_pd(cm.y));
      rz = _mm_sub_pd(posZ, _mm_set1_pd(cm.z));
      rsq = _mm_add_pd(_mm_mul_pd(rx,rx), _mm_add_pd(_mm_mul_pd(ry,ry), _mm_mul_pd(rz,rz)));
      twoh = _mm_add_pd(_mm_set1_pd(m.soft), _mm_setr_pd(input[i]->soft, input[i+1]->soft));
      __m128d dir = _mm_div_pd(_mm_set1_pd(1.0), _mm_sqrt_pd(rsq));
      //__m128d myvec = _mm_set_pd(rsq, r.x);
      //__m128d myvec2 = _mm_div_pd(_mm_set1_pd(1.0), _mm_sqrt_pd(myvec));
      //double dir;
      //_mm_store_sd(&dir, myvec2);
      //      SPLINEQ(dir, rsq, twoh, a, b, c, d);
      __m128d u,dih;
      __m128d select = _mm_cmplt_pd(rsq, _mm_mul_pd(twoh,twoh));
      int compare = _mm_movemask_pd(select);
      if (compare) {
        dih = _mm_div_pd(_mm_set1_pd(2.0), twoh);
        u = _mm_div_pd(dih, dir);
        __m128d select2 = _mm_cmplt_pd(u, _mm_set1_pd(1.0));
        int compare2 = _mm_movemask_pd(select2);
        __m128d u2 = _mm_mul_pd(u,u);
        __m128d u4 = _mm_mul_pd(u2,u2);
        __m128d dih2 = _mm_mul_pd(dih,dih);
        __m128d dih4 = _mm_mul_pd(dih2,dih2);
        if (compare2) {
          a = _mm_mul_pd(dih, _mm_add_pd(_mm_sub_pd(_mm_set1_pd(7.0/5.0), _mm_mul_pd(_mm_set1_pd(2.0/3.0),u2)), _mm_sub_pd(_mm_mul_pd(_mm_set1_pd(3.0/10.0),u4), _mm_mul_pd(_mm_set1_pd(1.0/10.0),_mm_mul_pd(u4,u)))));
          b = _mm_mul_pd(_mm_mul_pd(dih2,dih), _mm_add_pd(_mm_sub_pd(_mm_set1_pd(4.0/3.0), _mm_mul_pd(_mm_set1_pd(6.0/5.0),u2)), _mm_mul_pd(_mm_set1_pd(1.0/2.0),_mm_mul_pd(u2,u))));
          c = _mm_mul_pd(_mm_mul_pd(dih4,dih),_mm_sub_pd(_mm_set1_pd(12.0/5.0), _mm_mul_pd(_mm_set1_pd(3.0/2.0),u)));
          d = _mm_mul_pd(_mm_mul_pd(_mm_set1_pd(3.0/2.0),dih4), _mm_mul_pd(dih2,dir));
        }
        if ((~compare2) & 0x3) {
          a2 = _mm_add_pd(_mm_mul_pd(_mm_set1_pd(-1.0/15.0),dir), _mm_mul_pd(dih, _mm_add_pd(_mm_sub_pd(_mm_set1_pd(8.0/5.0), _mm_mul_pd(_mm_set1_pd(4.0/3.0),u2)), _mm_add_pd(_mm_sub_pd(_mm_mul_pd(u2,u), _mm_mul_pd(_mm_set1_pd(3.0/10.0),u4)), _mm_mul_pd(_mm_set1_pd(1.0/30.0),_mm_mul_pd(u4,u))))));
          b2 = _mm_add_pd(_mm_mul_pd(_mm_mul_pd(_mm_set1_pd(-1.0/15.0),dir),_mm_mul_pd(dir,dir)), _mm_mul_pd(_mm_mul_pd(dih2,dih),_mm_add_pd(_mm_sub_pd(_mm_set1_pd(8.0/3.0), _mm_mul_pd(_mm_set1_pd(3.0),u)), _mm_sub_pd(_mm_mul_pd(_mm_set1_pd(6.0/5.0),u2), _mm_mul_pd(_mm_set1_pd(1.0/6.0),u)))));
          c2 = _mm_add_pd(_mm_mul_pd(_mm_mul_pd(_mm_set1_pd(-1.0/5.0),dir),_mm_mul_pd(_mm_mul_pd(dir,dir),_mm_mul_pd(dir,dir))), _mm_add_pd(_mm_mul_pd(_mm_set1_pd(3.0),_mm_mul_pd(dih4,dir)), _mm_mul_pd(_mm_mul_pd(dih4,dih),_mm_add_pd(_mm_set1_pd(-12.0/5.0), _mm_mul_pd(_mm_set1_pd(1.0/2.0),u)))));
          d2 = _mm_sub_pd(_mm_sub_pd(_mm_mul_pd(_mm_set1_pd(3.0),_mm_mul_pd(_mm_mul_pd(dih4,dir),_mm_mul_pd(dir,dir))), _mm_mul_pd(_mm_mul_pd(dir,_mm_mul_pd(dir,dir)),_mm_mul_pd(_mm_mul_pd(dir,dir),_mm_mul_pd(dir,dir)))), _mm_mul_pd(_mm_mul_pd(_mm_set1_pd(1.0/2.0),dih4),_mm_mul_pd(dih2,dir)));
        }
        a1 = _mm_or_pd(_mm_and_pd(select2, a), _mm_andnot_pd(select2, a2));
        b1 = _mm_or_pd(_mm_and_pd(select2, b), _mm_andnot_pd(select2, b2));
        c1 = _mm_or_pd(_mm_and_pd(select2, c), _mm_andnot_pd(select2, c2));
        d1 = _mm_or_pd(_mm_and_pd(select2, d), _mm_andnot_pd(select2, d2));
      }
      if ((~compare) & 0x3) {
        a2 = dir;
        b2 = _mm_mul_pd(a2, _mm_mul_pd(a2,a2));
        c2 = _mm_mul_pd(_mm_mul_pd(_mm_set1_pd(3.0),b2), _mm_mul_pd(a2,a2));
        d2 = _mm_mul_pd(_mm_mul_pd(_mm_set1_pd(5.0),c2), _mm_mul_pd(a2,a2));
      }
      a = _mm_or_pd(_mm_and_pd(select, a1), _mm_andnot_pd(select, a2));
      b = _mm_or_pd(_mm_and_pd(select, b1), _mm_andnot_pd(select, b2));
      c = _mm_or_pd(_mm_and_pd(select, c1), _mm_andnot_pd(select, c2));
      d = _mm_or_pd(_mm_and_pd(select, d1), _mm_andnot_pd(select, d2));

      __m128d qirx = _mm_add_pd(_mm_mul_pd(_mm_set1_pd(m.xx),rx),
                     _mm_add_pd(_mm_mul_pd(_mm_set1_pd(m.xy),ry),
                                _mm_mul_pd(_mm_set1_pd(m.xz),rz)));
      __m128d qiry = _mm_add_pd(_mm_mul_pd(_mm_set1_pd(m.xy),rx),
                     _mm_add_pd(_mm_mul_pd(_mm_set1_pd(m.yy),ry),
                                _mm_mul_pd(_mm_set1_pd(m.yz),rz)));
      __m128d qirz = _mm_add_pd(_mm_mul_pd(_mm_set1_pd(m.xz),rx),
                     _mm_add_pd(_mm_mul_pd(_mm_set1_pd(m.yz),ry),
                                _mm_mul_pd(_mm_set1_pd(m.zz),rz)));
      __m128d qir = _mm_mul_pd(_mm_set1_pd(0.5),
                               _mm_add_pd(_mm_mul_pd(qirx,rx),
                                   _mm_add_pd(_mm_mul_pd(qiry,ry), _mm_mul_pd(qirz,rz))));
      __m128d qir3 = _mm_add_pd(_mm_mul_pd(b,_mm_set1_pd(m.totalMass)), _mm_sub_pd(_mm_mul_pd(d,qir), _mm_mul_pd(c,tr)));
      potential = _mm_sub_pd(_mm_setr_pd(input[i]->potential,input[i+1]->potential), _mm_add_pd(_mm_mul_pd(_mm_set1_pd(m.totalMass), a), _mm_sub_pd(_mm_mul_pd(c,qir), _mm_mul_pd(b,tr))));
      accX = _mm_sub_pd(_mm_setr_pd(input[i]->treeAcceleration.x,input[i+1]->treeAcceleration.x), _mm_sub_pd(_mm_mul_pd(qir3,rx), _mm_mul_pd(c,qirx)));
      accY = _mm_sub_pd(_mm_setr_pd(input[i]->treeAcceleration.y,input[i+1]->treeAcceleration.y), _mm_sub_pd(_mm_mul_pd(qir3,ry), _mm_mul_pd(c,qiry)));
      accZ = _mm_sub_pd(_mm_setr_pd(input[i]->treeAcceleration.z,input[i+1]->treeAcceleration.z), _mm_sub_pd(_mm_mul_pd(qir3,rz), _mm_mul_pd(c,qirz)));
      __m128d idt2 = _mm_mul_pd(_mm_add_pd(_mm_setr_pd(input[i]->mass,input[i+1]->mass), _mm_set1_pd(m.totalMass)), b);
      idt2 = _mm_max_pd(idt2, _mm_setr_pd(input[i]->dtGrav, input[i+1]->dtGrav));
      _mm_storel_pd(&input[i]->potential, potential);
      _mm_storeh_pd(&input[i+1]->potential, potential);
      _mm_storel_pd(&input[i]->treeAcceleration.x, accX);
      _mm_storeh_pd(&input[i+1]->treeAcceleration.x, accX);
      _mm_storel_pd(&input[i]->treeAcceleration.y, accY);
      _mm_storeh_pd(&input[i+1]->treeAcceleration.y, accY);
      _mm_storel_pd(&input[i]->treeAcceleration.z, accZ);
      _mm_storeh_pd(&input[i+1]->treeAcceleration.z, accZ);
      _mm_storel_pd(&input[i]->dtGrav, idt2);
      _mm_storeh_pd(&input[i+1]->dtGrav, idt2);
      //if(idt2 > input[i]->dtGrav)
      //  input[i]->dtGrav = idt2;
  }
#endif
  return computed;
}
#else
#ifdef COSMO_FLOAT
inline
int nodeBucketForce(Tree::GenericTreeNode *node, // source of force
            Tree::GenericTreeNode *req,  // bucket descriptor
            GravityParticle *particles, // particles in bucket
            Vector3D<double> offset,    // offset if a periodic replica
            int activeRung)         // rung (and above) at which to
                        // calculate forces
{
  int computed = 0;
  Vector3D<float> r;
  float rsq;
  float twoh, a, b, c, d;
  MultipoleMoments m = node->moments;
    
  Vector3D<float> cm(m.cm + offset);

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
      r = Vector3D<float>(particles[j].position) - cm;
      rsq = r.lengthSquared();
#ifdef HEXADECAPOLE
      // Here we assume that the separation is larger than the softening.
      float dir = 1.0/sqrt(rsq);
      momEvalMomr(&m.mom, dir, -r.x, -r.y, -r.z, &particles[j].potential,
          &particles[j].treeAcceleration.x,
          &particles[j].treeAcceleration.y,
          &particles[j].treeAcceleration.z);
      
      double idt2 = (particles[j].mass + m.totalMass)*dir*dir*dir;
#else
      twoh = (float)(m.soft + particles[j].soft);
        float dir = 1.0f/sqrt(rsq);
	//__m128d myvec = _mm_set_pd(rsq, r.x);
	//__m128d myvec2 = _mm_div_pd(_mm_set1_pd(1.0), _mm_sqrt_pd(myvec));
	//double dir;
	//_mm_store_sd(&dir, myvec2);
        SPLINEQ(dir, rsq, twoh, a, b, c, d);
        float qirx = (float)m.xx*r.x + (float)m.xy*r.y + (float)m.xz*r.z;
        float qiry = (float)m.xy*r.x + (float)m.yy*r.y + (float)m.yz*r.z;
        float qirz = (float)m.xz*r.x + (float)m.yz*r.y + (float)m.zz*r.z;
        float qir = 0.5f*(qirx*r.x + qiry*r.y + qirz*r.z);
        float tr = 0.5f*((float)m.xx + (float)m.yy + (float)m.zz);
        float qir3 = b*(float)m.totalMass + d*qir - c*tr;
        particles[j].potential -= (float)m.totalMass * a + c*qir - b*tr;
        particles[j].treeAcceleration.x -= qir3*r.x - c*qirx;
        particles[j].treeAcceleration.y -= qir3*r.y - c*qiry;
        particles[j].treeAcceleration.z -= qir3*r.z - c*qirz;
    float idt2 = ((float)particles[j].mass + (float)m.totalMass)*b;
#endif
    if(idt2 > particles[j].dtGrav)
        particles[j].dtGrav = idt2;
    }
  }
  return computed;
}
#else
inline
int nodeBucketForce(Tree::GenericTreeNode *node, // source of force
            Tree::GenericTreeNode *req,  // bucket descriptor
            GravityParticle *particles, // particles in bucket
            Vector3D<double> offset,    // offset if a periodic replica
            int activeRung)         // rung (and above) at which to
                        // calculate forces
{
  int computed = 0;
  Vector3D<double> r;
  double rsq;
  double twoh, a, b, c, d;
  MultipoleMoments m = node->moments;
    
  Vector3D<double> cm(m.cm + offset);

  //CkPrintf("nodeBucketForce: ext %f %f %f (%f %f %f) %f %f %f %f %f %f\n",node->moments.radius,node->moments.soft,node->moments.totalMass,cm.x,cm.y,cm.z,node->moments.xx,node->moments.xy,node->moments.xz,node->moments.yy,node->moments.yz,node->moments.zz);
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
      r = particles[j].position - cm;
      rsq = r.lengthSquared();
#ifdef HEXADECAPOLE
      // Here we assume that the separation is larger than the softening.
      double dir = 1.0/sqrt(rsq);
      momEvalMomr(&m.mom, dir, -r.x, -r.y, -r.z, &particles[j].potential,
          &particles[j].treeAcceleration.x,
          &particles[j].treeAcceleration.y,
          &particles[j].treeAcceleration.z);
      
      double idt2 = (particles[j].mass + m.totalMass)*dir*dir*dir;
#else
      twoh = m.soft + particles[j].soft;
        double dir = 1.0/sqrt(rsq);
        SPLINEQ(dir, rsq, twoh, a, b, c, d);
        double qirx = m.xx*r.x + m.xy*r.y + m.xz*r.z;
        double qiry = m.xy*r.x + m.yy*r.y + m.yz*r.z;
        double qirz = m.xz*r.x + m.yz*r.y + m.zz*r.z;
        double qir = 0.5*(qirx*r.x + qiry*r.y + qirz*r.z);
        double tr = 0.5*(m.xx + m.yy + m.zz);
        double qir3 = b*m.totalMass + d*qir - c*tr;
        particles[j].potential -= m.totalMass * a + c*qir - b*tr;
        particles[j].treeAcceleration.x -= qir3*r.x - c*qirx;
        particles[j].treeAcceleration.y -= qir3*r.y - c*qiry;
        particles[j].treeAcceleration.z -= qir3*r.z - c*qirz;
    double idt2 = (particles[j].mass + m.totalMass)*b;
#endif
    if(idt2 > particles[j].dtGrav)
        particles[j].dtGrav = idt2;
    //CkPrintf("nodeBucketForce of %d (%d): %f (%f %f %f) %f / (%f %f %f) %f %f\n",j,particles[j].iOrder,particles[j].mass,particles[j].position.x,particles[j].position.y,particles[j].position.z,particles[j].soft,particles[j].treeAcceleration.x,particles[j].treeAcceleration.y,particles[j].treeAcceleration.z,particles[j].potential,particles[j].dtGrav);
    }
  }
  return computed;
}
#endif
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
  /*
  float rsq = radius * radius;
  __m128 zero = _mm_setzero_ps();
  __m128 lesser = _mm_loadu_ps(&bucketNode->boundingBox.lesser_corner.x);
  __m128 greater = _mm_loadu_ps(&bucketNode->boundingBox.greater_corner.x);
  __m128 cm = _mm_loadu_ps(&node->moments.cm.x);
  __m128 off = _mm_loadu_ps(&offset.x);
  __m128 origin = _mm_add_ps(cm, off);
  __m128 delta1 = _mm_max_ps(_mm_sub_ps(lesser, origin), zero);
  __m128 delta2 = _mm_max_ps(_mm_sub_ps(origin, greater), zero);
  __m128 dsq = _mm_add_ps(_mm_mul_ps(delta1, delta1), _mm_mul_ps(delta2, delta2));
  float dsqf = _mm_cvtss_f32(_mm_add_ps(_mm_add_ps(dsq, _mm_shuffle_ps(dsq,dsq,1)), _mm_unpackhi_ps(dsq,dsq)));
  return dsqf <= rsq;
*/
  return Space::intersect(bucketNode->boundingBox, s);
#endif
}

// return 1 if there is an intersection
//    0 if no intersection
//    -1 if completely contained

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
