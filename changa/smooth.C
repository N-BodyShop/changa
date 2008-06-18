/*
 * Implementation of the "smooth" operations needed for SPH.
 * 
 * "smooth" itself is a search for the kth nearest neighbors.  No
 * calculations can be done until all the neighbors are found since
 * the operations will depend on the distance to the kth nearest neighbor.
 *
 * Proposed strategy for bucketwise smooth (this is different than
 * PKDGRAV, which does a particle at a time): each particle in a bucket
 * maintains its own priority queue.  Disadvantage: each "openCriterion"
 * results in nBucket tests, but most of those tests have to be performed
 * anyway.
 *
 * Priority queues are kept in the NearNeighborState class.
 *
 * How do I keep track of which particles are already in the queues?
 * 
 * Discussion of strategies:
 *
 * PKDGRAV: primes the priority queue with particles in tree order
 * starting with the bucket particles.
 * Once a particle is finished, the nearest unfinished particle in the
 * priority queue is done next ("the snake").
 * This requires: 1) a "smooth done" bit for each particle.  2) This
 * also requires a hash table to know if a particle is already in the
 * queue.
 * 
 * We would like to search more than one particle at once.  
 * That last bit means we need a hash table per particle, which seems
 * overkill.  We could avoid this by using a walk that 1) goes through
 * the particles in order, and 2) skips the preloaded particles.
 * Another way to avoid the hash table is to start the particle
 * priority queue with a single particle that is known not to be in
 * the final list, e.g. the 33rd closest particle.
 */

#include "ParallelGravity.h"
#include "TreeWalk.h"
#include "DataManager.h"
#include "Opt.h"
#include "smooth.h"
#include "Space.h"

// called after constructor, so tp should be set
State *SmoothCompute::getNewState(){
  NearNeighborState *state = new NearNeighborState(tp->myNumParticles+2);
  state->counterArrays.reserve(1);		   // array to keep track of
  state->counterArrays[0].reserve(tp->numBuckets); // outstanding requests
  for (int j=0; j<tp->numBuckets; ++j) {
    initSmoothPrioQueue(j, state);
    state->counterArrays[0][j] = 1;	// so we know that the local
					// walk is not finished.
  }
  state->started = true;
  state->nParticlesPending = tp->myNumParticles;
  return state;
}

/*
 * Opening criterion for the smoothBucket walk.
 * Return true if we must open the bucket.
 */
bool SmoothCompute::openCriterion(TreePiece *ownerTP, 
				  GenericTreeNode *node, // Node to test
				  int reqID, State *state) {
    GenericTreeNode *myNode = (GenericTreeNode *) computeEntity;
    GravityParticle *particles = ownerTP->getParticles();
    Vector3D<double> offset = ownerTP->decodeOffset(reqID);
    NearNeighborState *nstate = (NearNeighborState *)state;
    
    for(int j = myNode->firstParticle; j <= myNode->lastParticle; ++j) {
	double r = nstate->Qs[j].top().fKey; // Ball radius
	Sphere<double> s(particles[j].position - offset, r);
	if(Space::intersect(node->boundingBox, s)) {
	   return true;
	   }
	}
    return false;
}

/*
 * Test a given particle against all the priority queues in the
 * bucket.
 */
void SmoothCompute::bucketCompare(TreePiece *ownerTP,
				  ExternalGravityParticle *p,  // Particle to test
				  GenericTreeNode *node, // bucket
				  GravityParticle *particles, // local
							      // particle data
				  Vector3D<double> offset,
                                  State *state
				  ) 
{
    NearNeighborState *nstate = (NearNeighborState *)state;
    for(int j = node->firstParticle; j <= node->lastParticle; ++j) {
	std::priority_queue<pqSmoothNode> *Q = &nstate->Qs[j];
	double rOld = Q->top().fKey; // Ball radius
	Vector3D<double> dr = offset + p->position - particles[j].position;
	
	if(rOld*rOld > dr.lengthSquared()) {  // Perform replacement
	    if(Q->size() == nSmooth)
		Q->pop();	// pop if list is full
	    pqSmoothNode pqNew;
	    pqNew.fKey = dr.length();
	    pqNew.dx = dr;
	    pqNew.p = p;
	    Q->push(pqNew);
	    }
	}
    }

/*
 * There is actually not much "work" for the smooth walk.  This method
 * makes the decision about whether to keep a node or not.  If we keep
 * a node, and its a bucket, then we check each particle with
 * bucketCompare().
 */
int SmoothCompute::doWork(GenericTreeNode *node, // Node to test
			  TreeWalk *tw,
			  State *state,
			  int chunk,
			  int reqID, // bucket we are working with
			  bool isRoot, 
			  bool &didcomp, int awi)
{
    // jetley - no need for this, tp was set in constructor
    //TreePiece *tp = tw->getOwnerTP();
    
    // so that we have a quick return in case of empty nodes
    if(node->getType() == Empty || node->getType() == CachedEmpty){
	return DUMP;
	}
    // check opening criterion
    bool open = openCriterion(tp, node, reqID, state);
    // Turn open into an action
    int action = opt->action(open, node);
    if(action == KEEP){
	// Bounds intersect ball; descend to children
	return KEEP;
	}
    else if(action == COMPUTE) {
	CkAssert(0);
	}
    else if(action == KEEP_LOCAL_BUCKET) {
	// Search bucket for contained particles
	GravityParticle *part = node->particlePointer;
	
	for(int i = node->firstParticle; i <= node->lastParticle; i++) {
	    bucketCompare(tp, &part[i-node->firstParticle],
			  (GenericTreeNode *) computeEntity,
			  tp->getParticles(),
			  tp->decodeOffset(reqID),
                          state);
	    }
	return DUMP;
	}
    else if(action == KEEP_REMOTE_BUCKET) {
	ExternalGravityParticle *part;
	part = tp->requestParticles(node->getKey(), 
				    chunk, 
				    node->remoteIndex, 
				    node->firstParticle, 
				    node->lastParticle, 
				    reqID, awi, false);
	if(part) {
	    // Particles available; Search for contained.
	    for(int i = node->firstParticle; i <= node->lastParticle; i++) {
		bucketCompare(tp, &part[i-node->firstParticle],
			      (GenericTreeNode *) computeEntity,
			      tp->getParticles(),
			      tp->decodeOffset(reqID),
                              state);
		}
	    }
	else {
	    // Missed cache; record this
	    state->counterArrays[0][decodeReqID(reqID)]
		+= node->lastParticle - node->firstParticle+1;
	    }
	
	return DUMP;
	}
    else if(action == DUMP) {
	return DUMP;
	}
    }

/*
 * Process particles received from missed Cache request
 */
void SmoothCompute::recvdParticles(ExternalGravityParticle *part,
				   int num, int chunk,int reqID, State *state,
				   TreePiece *tp){

  Vector3D<double> offset = tp->decodeOffset(reqID);
  int reqIDlist = decodeReqID(reqID);
  CkAssert(num > 0);
  state->counterArrays[0][reqIDlist] -= num;

  GenericTreeNode* reqnode = tp->bucketList[reqIDlist];

  int computed;
  for(int i=0;i<num;i++){
    
      bucketCompare(tp, &part[i], reqnode, tp->myParticles, offset, state);
      }
  ((NearNeighborState *)state)->finishBucketSmooth(reqIDlist, tp);
}

int SmoothCompute::nodeMissedEvent(int reqID, int chunk, State *state)
{
    int reqIDlist = decodeReqID(reqID);
    state->counterArrays[0][reqIDlist]++;
    }

int SmoothCompute::nodeRecvdEvent(TreePiece *owner, int chunk, State *state,
				  int reqIDlist){
  state->counterArrays[0][reqIDlist]--;
  ((NearNeighborState *)state)->finishBucketSmooth(reqIDlist, owner);
}

// Start tree walk and smooth calculation

void TreePiece::startIterationSmooth(int am, // the active rung for
					     // multistepping
				     const CkCallback& cb) {

  callback = cb;
  activeRung = am;

  // XXX I don't believe any of the Chunks are used in the smooth walk.
  int oldNumChunks = numChunks;
  dm->getChunks(numChunks, prefetchRoots);
  CkArrayIndexMax idxMax = CkArrayIndex1D(thisIndex);
  streamingCache[CkMyPe()].cacheSync(numChunks, idxMax, localIndex);
  if (oldNumChunks != numChunks && remainingChunk != NULL) {
    // reallocate remaining chunk to the new size
    delete[] remainingChunk;
    remainingChunk = NULL;
    delete[] nodeInterRemote;
    delete[] particleInterRemote;
  }
  if (remainingChunk == NULL) {
    remainingChunk = new int[numChunks];
    nodeInterRemote = new u_int64_t[numChunks];
    particleInterRemote = new u_int64_t[numChunks];
  }
  nodeInterLocal = 0;
  for (int i=0; i<numChunks; ++i) {
    nodeInterRemote[i] = 0;
    particleInterRemote[i] = 0;
  }
  particleInterLocal = 0;
  nActive = 0;
  iterationNo++;

  //CkAssert(localCache != NULL);
  if(verbosity>1)
    CkPrintf("Node: %d, TreePiece %d: I have %d buckets\n", CkMyNode(),
    	     thisIndex,numBuckets);

  currentBucket = 0;
  currentRemoteBucket = 0;

  initBucketsSmooth();

  for (int i=0; i<numChunks; ++i) remainingChunk[i] = myNumParticles;

  // Create objects that are reused by all buckets
  sTopDown = new TopDownTreeWalk;
  sSmooth = new SmoothCompute(this, Density, nSmooth);
  // creates and initializes nearneighborstate object
  sSmoothState = sSmooth->getNewState();
  optSmooth = new SmoothOpt;

  smoothAwi = addActiveWalk(sTopDown,sSmooth,optSmooth,sSmoothState);

  thisProxy[thisIndex].calculateSmoothLocal();
}

// Initialize particles in bucket and bucket bookkeeping

void TreePiece::initBucketsSmooth() {
  for (unsigned int j=0; j<numBuckets; ++j) {
    GenericTreeNode* node = bucketList[j];
    int numParticlesInBucket = node->particleCount;

    CkAssert(numParticlesInBucket <= TreeStuff::maxBucketSize);
    
    // TODO: active bounds may give a performance boost in the
    // multi-timstep regime.
    // node->boundingBox.reset();  // XXXX dangerous should have separate
				// Active bounds
    for(int i = node->firstParticle; i <= node->lastParticle; ++i) {
      if (myParticles[i].rung >= activeRung) {
	nActive++;
        myParticles[i].fDensity = 0;
	// node->boundingBox.grow(myParticles[i].position);
      }
    }
  }
}

// Start the smoothing

void TreePiece::calculateSmoothLocal() {
    dummyMsg *msg = new (8*sizeof(int)) dummyMsg;
    *((int *)CkPriorityPtr(msg)) = numTreePieces * numChunks + thisIndex + 1;
    CkSetQueueing(msg,CK_QUEUEING_IFIFO);
    msg->val=0;
    thisProxy[thisIndex].nextBucketSmooth(msg);
    }

//
// Do the next set of buckets
//
void TreePiece::nextBucketSmooth(dummyMsg *msg){
  unsigned int i=0;

  while(i<_yieldPeriod && currentBucket<numBuckets){
    smoothNextBucket();
    currentBucket++;
    i++;
  }

  if (currentBucket<numBuckets) {	// Queue up the next set
    thisProxy[thisIndex].nextBucketSmooth(msg);
  } else {
    delete msg;
  }
}

void SmoothCompute::initSmoothPrioQueue(int iBucket, State *state) 
{
  // Prime the queues
  GenericTreeNode *myNode = tp->bucketList[iBucket];
  // state is passed in to function now. 
  NearNeighborState *nstate = (NearNeighborState *)state;
  
  //
  // Get nearest nSmooth particles in tree order
  //
  int firstQueue = myNode->firstParticle;
  int lastQueue = firstQueue + nSmooth + 1; // XXX I'm using one more
					    // particle than needed to
					    // ensure it eventually
					    // gets popped off.
  if(lastQueue > tp->myNumParticles) 
      {
	  lastQueue = tp->myNumParticles;
	  firstQueue = lastQueue - (nSmooth + 1);
	  CkAssert(firstQueue > 0);
	  }
	  
  for(int j = myNode->firstParticle; j <= myNode->lastParticle; ++j) {
      std::priority_queue<pqSmoothNode> *Q = &(nstate->Qs[j]);
      //
      // Find maximum of nearest neighbors
      //
      double drMax2 = 0.0;
      int kMax = 0;
      for(int k = firstQueue; k < lastQueue; k++) 
	  {
	      Vector3D<double> dr = tp->myParticles[k].position
		  - tp->myParticles[j].position;
	      if(dr.lengthSquared() > drMax2) {
		  kMax = k;
		  drMax2 = dr.lengthSquared();
		  }
	      }
      pqSmoothNode pqNew;
      pqNew.fKey = sqrt(drMax2);
      pqNew.p = NULL;  // Make sure this never used
      Q->push(pqNew);
      }
    }

void TreePiece::smoothNextBucket() {
  if(currentBucket >= numBuckets)
    return;

  sTopDown->init(sSmooth, this);
  sSmooth->init(bucketList[currentBucket], activeRung, optSmooth);
  State *smoothState = activeWalks[smoothAwi].s;

  // start the tree walk from the tree built in the cache
  if (bucketList[currentBucket]->rungs >= activeRung) {
    for(int cr = 0; cr < numChunks; cr++){
      GenericTreeNode *chunkRoot = dm->chunkRootToNode(prefetchRoots[cr]);
      for(int x = -nReplicas; x <= nReplicas; x++) {
        for(int y = -nReplicas; y <= nReplicas; y++) {
          for(int z = -nReplicas; z <= nReplicas; z++) {
            sTopDown->walk(chunkRoot, smoothState, cr, encodeOffset(currentBucket, x,y,z), smoothAwi);
          }
        }
      }
    }
  }
  smoothState->counterArrays[0][currentBucket]--;
  ((NearNeighborState *)smoothState)->finishBucketSmooth(currentBucket, this);
}

void NearNeighborState::finishBucketSmooth(int iBucket, TreePiece *tp) {
  GenericTreeNode *node = tp->bucketList[iBucket];

  if(counterArrays[0][iBucket] == 0) {
      tp->sSmooth->walkDone(this);
    nParticlesPending -= node->particleCount;
    if(started && nParticlesPending == 0) {
      started = false;
      cacheManagerProxy[CkMyPe()].finishedChunk(0, 0);
      tp->markWalkDone();
      if(verbosity>1)
	CkPrintf("[%d] TreePiece %d finished smooth with bucket %d\n",CkMyPe(),
		 tp->thisIndex,iBucket);
      if(verbosity > 4)
	ckerr << "TreePiece " << tp->thisIndex << ": My particles are done"
	     << endl;
    }
  }
}

void SmoothCompute::walkDone(State *state) {
  GenericTreeNode *node = (GenericTreeNode *) computeEntity;
  GravityParticle *part = node->particlePointer;

  for(int i = node->firstParticle; i <= node->lastParticle; i++) {
      std::priority_queue<pqSmoothNode> *Q = &((NearNeighborState *)state)->Qs[i];
      double h = Q->top().fKey; // Ball radius
      part[i-node->firstParticle].fBall = h;
      pqSmoothNode NN[Q->size()];
      int nCnt = 0;
      while (!Q->empty()) {
	  NN[nCnt] = Q->top();
	  Q->pop();
	  nCnt++;
	  }
      fcnSmooth(&part[i-node->firstParticle], nCnt, NN);
      }
      // XXX jetley - the nearneighborstate allocated for this compute
      // may be deleted after this function is called. 
}

/* Standard M_4 Kernel */
/* return 1/(h_smooth)^2 for a particle */
double invH2(GravityParticle *p) 
{
    return 4.0/(p->fBall*p->fBall);
    }

inline double KERNEL(double ar2) 
{
    double ak;
    ak = 2.0 - sqrt(ar2);
    if (ar2 < 1.0) ak = (1.0 - 0.75*ak*ar2);
    else ak = 0.25*ak*ak*ak;
    return ak;
    }
inline double DKERNEL(double ar2) 
{
    double adk;
    adk = sqrt(ar2);
    if (ar2 < 1.0) {
	adk = -3 + 2.25*adk;
	}
    else {
	adk = -0.75*(2.0-adk)*(2.0-adk)/adk;
	}
    return adk;
    }
/*
 * XXX Place holder from PKDGRAV
 */
void Density(GravityParticle *p,int nSmooth,pqSmoothNode *nnList)
{
	double ih2,r2,rs,fDensity;
	int i;

	ih2 = invH2(p);
	fDensity = 0.0;
	for (i=0;i<nSmooth;++i) {
		double fDist2 = nnList[i].fKey*nnList[i].fKey;
		r2 = fDist2*ih2;
		rs = KERNEL(r2);
		fDensity += rs*nnList[i].p->mass;
		}
	p->fDensity = M_1_PI*sqrt(ih2)*ih2*fDensity; 
	}
