#ifndef GRAVITYPARTICLE_H
#define GRAVITYPARTICLE_H

#include "SFC.h"
#include <vector>

// Object to bookkeep a Bucket Walk.  It has accumulators for the
// acclerations and potentials of the particles in the bucket.

class BucketGravityRequest {
public:
		
	SFC::Key startingNode;
	unsigned int identifier;
	unsigned int requestingPieceIndex;
	OrientedBox<double> boundingBox;
	unsigned int numParticlesInBucket;
	double *softs;
	Vector3D<double>* positions;
	Vector3D<double>* accelerations;
	double *potentials;
	unsigned int numAdditionalRequests;
	int finished;

        /*#if COSMO_DEBUG > 1
          std::vector<u_int64_t> requestedNodes;
          #endif*/
  
	BucketGravityRequest(unsigned int bucketSize = 0) : identifier(0),
	    numParticlesInBucket(bucketSize), numAdditionalRequests(0),
	    finished(0) {
		if(numParticlesInBucket) {
			softs = new double[numParticlesInBucket];
			positions = new Vector3D<double>[numParticlesInBucket];
			accelerations = new Vector3D<double>[numParticlesInBucket];
			potentials = new double[numParticlesInBucket];
		    } else {
			positions = accelerations = 0;
			softs = potentials = 0;
			}
                /*#if COSMO_DEBUG > 1
                  requestedNodes.clear();
                  #endif*/
	}
	
	BucketGravityRequest(const BucketGravityRequest& req) {
		startingNode = req.startingNode;
		identifier = req.identifier;
		requestingPieceIndex = req.requestingPieceIndex;
		boundingBox = req.boundingBox;
		numParticlesInBucket = req.numParticlesInBucket;
		finished = req.finished;
		if(numParticlesInBucket) {
			softs = new double[numParticlesInBucket];
			positions = new Vector3D<double>[numParticlesInBucket];
			accelerations = new Vector3D<double>[numParticlesInBucket];
			potentials = new double[numParticlesInBucket];
			for(unsigned int i = 0; i < numParticlesInBucket; ++i) {
				softs[i] = req.softs[i];
				positions[i] = req.positions[i];
				accelerations[i] = req.accelerations[i];
				potentials[i] = req.potentials[i];
			}
		    } else {
			positions = accelerations = 0;
			softs = potentials = 0;
			}
		numAdditionalRequests = req.numAdditionalRequests;
	}
	
	BucketGravityRequest& operator=(const BucketGravityRequest& req) {
		startingNode = req.startingNode;
		identifier = req.identifier;
		requestingPieceIndex = req.requestingPieceIndex;
		boundingBox = req.boundingBox;
		numParticlesInBucket = req.numParticlesInBucket;
		finished = req.finished;
		delete[] positions;
		delete[] accelerations;
		delete[] softs;
		delete[] potentials;
		if(numParticlesInBucket) {
			softs = new double[numParticlesInBucket];
			positions = new Vector3D<double>[numParticlesInBucket];
			accelerations = new Vector3D<double>[numParticlesInBucket];
			potentials = new double[numParticlesInBucket];
			for(unsigned int i = 0; i < numParticlesInBucket; ++i) {
				softs[i] = req.softs[i];
				positions[i] = req.positions[i];
				accelerations[i] = req.accelerations[i];
				potentials[i] = req.potentials[i];
			}
		    } else {
			positions = accelerations = 0;
			softs = potentials = 0;
			}
		numAdditionalRequests = req.numAdditionalRequests;
                /*#if COSMO_DEBUG > 1
                  for(unsigned int i=0;i<req.requestedNodes.size();i++)
                  requestedNodes.push_back(req.requestedNodes[i]);
                  #endif*/
		return *this;
	}
	
	~BucketGravityRequest() {
		delete[] softs;
		delete[] positions;
		delete[] accelerations;
		delete[] potentials;
	}
	
	void merge(const BucketGravityRequest& req) {
	    for(unsigned int i = 0; i < numParticlesInBucket; ++i) {
		    accelerations[i] += req.accelerations[i];
		    potentials[i] += req.potentials[i];
	    }
	}
	
	void pup(PUP::er &p) {
		p | startingNode;
		p | identifier;
		p | requestingPieceIndex;
		p | boundingBox;
		p | numParticlesInBucket;
		if(p.isUnpacking()) {
			if(numParticlesInBucket) {
				softs = new double[numParticlesInBucket];
				positions = new Vector3D<double>[numParticlesInBucket];
				accelerations = new Vector3D<double>[numParticlesInBucket];
				potentials = new double[numParticlesInBucket];
			}
		}
		for(unsigned int i = 0; i < numParticlesInBucket; ++i) {
			p | softs[i];
			p | positions[i];
			p | accelerations[i];
			p | potentials[i];
		}	
		p | numAdditionalRequests;
	}	
};

/*
class GravityRequest {
public:
	SFC::Key startingNode;
	unsigned int identifier;
	unsigned int requestingPieceIndex;
	double soft;
	Vector3D<double> position;
	Vector3D<double> acceleration;
	double potential;
	unsigned int numAdditionalRequests;
	unsigned int numCellInteractions;
	unsigned int numParticleInteractions;
	unsigned int numMACChecks;
	unsigned int numEntryCalls;
	
	GravityRequest() : identifier(0), requestingPieceIndex(0), numAdditionalRequests(0), numCellInteractions(0), numParticleInteractions(0), numMACChecks(0), numEntryCalls(0) { }
	
	void merge(const GravityRequest& req) {
		acceleration += req.acceleration;
		potential += req.potential;
		numCellInteractions += req.numCellInteractions;
		numParticleInteractions += req.numParticleInteractions;
		numMACChecks += req.numMACChecks;
		numEntryCalls += req.numEntryCalls;
	}
	
	void pup(PUP::er &p) {
		p | startingNode;
		p | identifier;
		p | requestingPieceIndex;
		p | soft;
		p | position;
		p | acceleration;
		p | potential;
		p | numAdditionalRequests;
		p | numCellInteractions;
		p | numParticleInteractions;
		p | numMACChecks;
		p | numEntryCalls;
	}
};
*/

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
	//double mass;
	//double soft;
	//Vector3D<double> position;
	Vector3D<double> velocity;
	//Vector3D<double> acceleration;
	Vector3D<double> treeAcceleration;
	double potential;
	int iOrder;		/* input order of particles */
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
        }

        /*
	GravityParticle(float m, Vector3D<float> p, Vector3D<float> a) : mass(m), position(p), acceleration(a)
          //, numCellInteractions(0), numParticleInteractions(0), numMACChecks(0), numEntryCalls(0)
          { }
        */

	inline bool operator<(const GravityParticle& p) const {
		return key < p.key;
	}

        /*	
	void update(const GravityRequest& req) {
		treeAcceleration += req.acceleration;
		potential += req.potential;
		numCellInteractions += req.numCellInteractions;
		numParticleInteractions += req.numParticleInteractions;
		numMACChecks += req.numMACChecks;
		numEntryCalls += req.numEntryCalls;
	}
        */

	void pup(PUP::er &p) {
          ExternalGravityParticle::pup(p);
          p | key;
          //p | position;
          p | velocity;
          //p | mass;
          //p | soft;
          p | iOrder;
        }
};

#endif
