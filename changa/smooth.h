/*
 * structure for the priority queue.  Also contains information for
 * smooth calculation.
 */
#include <queue>
#include "Compute.h"

class pqSmoothNode
{
 public:
    double fKey;	// distance squared -> place in priority queue
    Vector3D<double> dx; // displacement of this particle
    // XXX to be replaced with "ExternalSmoothParticle"
    ExternalGravityParticle *p; // pointer to rest of particle data
    
    inline bool operator<(const pqSmoothNode& n) const {
	return fKey < n.fKey;
	}
    };

	
// Object to bookkeep a Bucket Smooth Walk.

class NearNeighborState: public State {
public:
	priority_queue <pqSmoothNode> *Qs;
    NearNeighborState(int nParts) 
	{
	    Qs = new (priority_queue<pqSmoothNode>[nParts+2]);
	    }
    

};

class SmoothCompute : public Compute 
{
    int nSmooth;
    void (*fcnSmooth)(GravityParticle *p, int nSmooth, pqSmoothNode *nList);
    void bucketCompare(TreePiece *tp,
		       ExternalGravityParticle *p,  // Particle to test
		       GenericTreeNode *node, // bucket
		       GravityParticle *particles, // local particle data
		       Vector3D<double> offset
		       ) ;
    
public:
 SmoothCompute(void (*fcn)(GravityParticle *p, int nSmooth,
			   pqSmoothNode *nList), int nSm)
     : Compute(Smooth){
	nSmooth = nSm;
	fcnSmooth = fcn;
	}
    bool openCriterion(TreePiece *ownerTP, GenericTreeNode *node, int reqID);
    
    int doWork(GenericTreeNode *node,
	       TreeWalk *tw,
	       State *state,
	       int chunk,
	       int reqID,
	       bool isRoot, 
	       bool &didcomp, int awi);
    int startNodeProcessEvent(TreePiece *owner){}
    int finishNodeProcessEvent(TreePiece *owner){}
    int nodeMissedEvent(TreePiece *owner, int chunk);
    State *getResumeState(int bucketIdx);
    void walkDone() ;
    };

void Density(GravityParticle *p,int nSmooth,pqSmoothNode *nnList);
