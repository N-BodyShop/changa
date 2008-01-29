/*
 ** see (A1) and (A2) of TREESPH: A UNIFICATION OF SPH WITH THE 
 ** HIERARCHICAL TREE METHOD by Lars Hernquist and Neal Katz.
 ** APJ Supplemant Series 70:416-446, 1989
 ** 
 ** Higher derivative terms c and d for use with quadrupole spline
 ** softening (Joachim Stadel, Dec. 94).
 */
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
int openSoftening(GenericTreeNode *node, GenericTreeNode *myNode,
          Vector3D<double> offset)
{
    Sphere<double> s(node->moments.cm + offset, 2.0*node->moments.soft);
    Sphere<double> myS(myNode->moments.cm, 2.0*myNode->moments.soft);
    return Space::intersect(myS, s);
    }

#ifdef CMK_VERSION_BLUEGENE
static int forProgress = 0;
#endif

inline int partBucketForce(ExternalGravityParticle *part, GenericTreeNode *req, GravityParticle *particles, Vector3D<double> offset, int activeRung) {
  int computed = 0;
  Vector3D<double> r;
  double rsq;
  double twoh, a, b;

  //CkPrintf("partBucketForce: ext %f (%f %f %f) %f\n",part->mass,part->position.x,part->position.y,part->position.z,part->soft);
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
      //CkPrintf("partBucketForce of %d (%d): %f (%f %f %f) %f / (%f %f %f) %f %f\n",j,particles[j].iOrder,particles[j].mass,particles[j].position.x,particles[j].position.y,particles[j].position.z,particles[j].soft,particles[j].treeAcceleration.x,particles[j].treeAcceleration.y,particles[j].treeAcceleration.z,particles[j].potential,particles[j].dtGrav);
    }
  }
  return computed;
}

inline
int nodeBucketForce(GenericTreeNode *node, // source of force
            GenericTreeNode *req,  // bucket descriptor
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

// Return true if the node's opening radius intersects the
// boundingBox, i.e. the node needs to be opened.

inline bool
openCriterionBucket(GenericTreeNode *node,
                   GenericTreeNode *bucketNode,
                   Vector3D<double> offset, // Offset of node
                   int localIndex // requesting TreePiece
                   ) {
  // mark the node as used by the requesting TreePiece
  node->markUsedBy(localIndex);

#if COSMO_STATS > 0
  node->used = true;
#endif
  // Note that some of this could be pre-calculated into an "opening radius"
  double radius = opening_geometry_factor * node->moments.radius / theta;
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
      radius = opening_geometry_factor*node->moments.radius/thetaMono;
      Sphere<double> sM(node->moments.cm + offset, radius);
      return Space::intersect(bucketNode->boundingBox, sM);
      }
      }
  return true;
#else
  return Space::intersect(bucketNode->boundingBox, s);
#endif
}

// return 1 if there is an intersection
//    0 if no intersection
//    -1 if completely contained

inline int openCriterionNode(GenericTreeNode *node,
                    GenericTreeNode *myNode,
                    Vector3D<double> offset,
                    int localIndex // requesting TreePiece
                    ) {
  // mark the node as used by this TreePiece
  node->markUsedBy(localIndex);

#if COSMO_STATS > 0
  node->used = true;
#endif
  // Note that some of this could be pre-calculated into an "opening radius"
  double radius = opening_geometry_factor * node->moments.radius / theta;
  if(radius < node->moments.radius)
      radius = node->moments.radius;

  Sphere<double> s(node->moments.cm + offset, radius);

  if(myNode->getType()==Bucket || myNode->getType()==CachedBucket || myNode->getType()==NonLocalBucket){
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
            radius = opening_geometry_factor*node->moments.radius/thetaMono;
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
                radius = opening_geometry_factor*node->moments.radius/thetaMono;
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

