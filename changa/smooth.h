/*
 * structure for the priority queue.  Also contains information for
 * smooth calculation.
 */
class pqSmoothNode
{
 public:
    double fKey;	// distance squared -> place in priority queue
    Vector3D<double> dx; // displacement of this particle
    ExternalSmoothParticle *p; // pointer to rest of particle data
    
    inline bool operator<(const pqSmoothNode& n) const {
	return fKey < n.fKey;
	}
    };

	
// Object to bookkeep a Bucket Smooth Walk.

class NearNeighborState: State {
public:
	priority_queue<pqSmoothNode> *Qs;

};

class SmoothCompute : Compute 
{
    int doWork(GenericTreeNode *node,
	       TreeWalk *tw,
	       State *state,
	       int chunk,
	       bool isRoot, 
	       bool &didcomp)
    };

