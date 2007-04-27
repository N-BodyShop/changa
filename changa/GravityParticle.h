#ifndef GRAVITYPARTICLE_H
#define GRAVITYPARTICLE_H

#include "SFC.h"
#include <vector>

// Object to bookkeep a Bucket Walk.  It has accumulators for the
// acclerations and potentials of the particles in the bucket.

class BucketGravityRequest {
public:
		
	unsigned int numAdditionalRequests;
	int finished;

	BucketGravityRequest(unsigned int bucketSize = 0) :
	  //identifier(0), numParticlesInBucket(bucketSize),
	  numAdditionalRequests(0), finished(0) {
	}
	
};

// Information needed to calculate gravity

class ExternalGravityParticle {
 public:
  double mass;
  double soft;
  Vector3D<double> position;

  void pup(PUP::er &p) {
    p | position;
    p | mass;
    p | soft;
  }
};

class GravityParticle : public ExternalGravityParticle {
public:

	SFC::Key key;
	Vector3D<double> velocity;
	Vector3D<double> treeAcceleration;
	double potential;
	double dtGrav;
	int iOrder;		/* input order of particles */
        int rung;  ///< the current rung (greater means faster)
	//unsigned int numCellInteractions;
	//unsigned int numParticleInteractions;
	//unsigned int numMACChecks;
	//unsigned int numEntryCalls;

#if COSMO_STATS > 1
	double intcellmass;
	double intpartmass;
	double extcellmass;
	double extpartmass;
#endif
	
	GravityParticle(SFC::Key k = 0) : ExternalGravityParticle() {
          key = k;
          rung = 0;
        }

	inline bool operator<(const GravityParticle& p) const {
		return key < p.key;
	}

	void pup(PUP::er &p) {
          ExternalGravityParticle::pup(p);
          p | key;
          p | velocity;
          p | iOrder;
          p | rung;
        }
};

#endif
