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
 * the final list, e.g. the 33rd closest particle. This is what is
 * done below.
 */

#include "ParallelGravity.h"
#include "TreeWalk.h"
#include "DataManager.h"
#include "Opt.h"
#include "smooth.h"
#include "Space.h"


SmoothParams *globalSmoothParams;

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
    int open = openCriterion(tp, node, reqID, state);
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
	    if(!TYPETest(&part[i-node->firstParticle], params->iType))
		continue;
	    bucketCompare(tp, &part[i-node->firstParticle],
			  (GenericTreeNode *) computeEntity,
			  tp->getParticles(),
			  tp->decodeOffset(reqID),
                          state);
	    }
	return DUMP;
	}
    else if(action == KEEP_REMOTE_BUCKET) {
	GravityParticle *part;
	part = tp->requestSmoothParticles(node->getKey(), 
				    chunk, 
				    node->remoteIndex, 
				    node->firstParticle, 
				    node->lastParticle, 
				    reqID, awi, computeEntity, false);
	if(part) {
	    // Particles available; Search for contained.
	    for(int i = node->firstParticle; i <= node->lastParticle; i++) {
		if(!TYPETest(&part[i-node->firstParticle], params->iType))
		    continue;
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
 * reassociate bucket
 */
void SmoothCompute::reassoc(void *ce, int ar, Opt *o){
  computeEntity = ce; 
  activeRung = ar;
  opt = o;
}

int SmoothCompute::nodeMissedEvent(int reqID, int chunk, State *state, TreePiece *tp)
{
    int reqIDlist = decodeReqID(reqID);
    state->counterArrays[0][reqIDlist]++;
    }

// called after constructor, so tp should be set
State *KNearestSmoothCompute::getNewState(int d1){
  NearNeighborState *state = new NearNeighborState(tp->myNumParticles+2, nSmooth);
  //state->counterArrays.reserve(1);		   // array to keep track of
  state->counterArrays[0] = new int [d1];
  state->counterArrays[1] = 0;
  //state->counterArrays[0].reserve(tp->numBuckets); // outstanding requests
  //for (int j=0; j<tp->numBuckets; ++j) {
  for(int j = 0; j < d1; ++j){
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
 * Return true if we must open the node.
 */
int KNearestSmoothCompute::openCriterion(TreePiece *ownerTP, 
				  GenericTreeNode *node, // Node to test
				  int reqID, State *state) {
    GenericTreeNode *myNode = (GenericTreeNode *) computeEntity;
    GravityParticle *particles = ownerTP->getParticles();
    Vector3D<double> offset = ownerTP->decodeOffset(reqID);
    NearNeighborState *nstate = (NearNeighborState *)state;
    
    for(int j = myNode->firstParticle; j <= myNode->lastParticle; ++j) {
	if(!TYPETest(&particles[j], params->iType))
	    continue;
	double r = nstate->Qs[j][0].fKey; // Ball radius
	Sphere<double> s(particles[j].position - offset, r);
	if(Space::intersect(node->boundingBox, s)) {
	   return 1;
	   }
	}
    return 0;
}

/*
 * Test a given particle against all the priority queues in the
 * bucket.
 */

void KNearestSmoothCompute::bucketCompare(TreePiece *ownerTP,
				  GravityParticle *p,  // Particle to test
				  GenericTreeNode *node, // bucket
				  GravityParticle *particles, // local
							      // particle data
				  Vector3D<double> offset,
                                  State *state
				  ) 
{
    NearNeighborState *nstate = (NearNeighborState *)state;
    Vector3D<double> rp = offset + p->position;
    for(int j = node->firstParticle; j <= node->lastParticle; ++j) {
	if(!TYPETest(&particles[j], params->iType))
	    continue;
	pqSmoothNode *Q = nstate->Qs[j];
	double rOld = Q[0].fKey; // Ball radius
	Vector3D<double> dr = particles[j].position - rp;
	
	if(rOld*rOld > dr.lengthSquared()) {  // Perform replacement
	    if(nstate->heap_sizes[j] == nSmooth) {
	        std::pop_heap(Q + 0, Q + nSmooth);
	        nstate->heap_sizes[j]--;	// pop if list is full
	        }
	    int end = nstate->heap_sizes[j];
	    Q[end].fKey = dr.length();
	    Q[end].dx = dr;
	    Q[end].p = p;
	    std::push_heap(Q + 0, Q + end + 1); 
	    nstate->heap_sizes[j]++;
	    }
	}
    }


/*
 * Process particles received from missed Cache request
 */
void KNearestSmoothCompute::recvdParticlesFull(GravityParticle *part,
				   int num, int chunk,int reqID, State *state,
				   TreePiece *tp, Tree::NodeKey &remoteBucket){

  Vector3D<double> offset = tp->decodeOffset(reqID);
  int reqIDlist = decodeReqID(reqID);
  CkAssert(num > 0);
  state->counterArrays[0][reqIDlist] -= num;

  GenericTreeNode* reqnode = tp->bucketList[reqIDlist];

  for(int i=0;i<num;i++){
      if(!TYPETest(&part[i], params->iType))
	  continue;
      bucketCompare(tp, &part[i], reqnode, tp->myParticles, offset, state);
      }
  ((NearNeighborState *)state)->finishBucketSmooth(reqIDlist, tp);
}

int KNearestSmoothCompute::nodeRecvdEvent(TreePiece *owner, int chunk, State *state,
				  int reqIDlist){
  state->counterArrays[0][reqIDlist]--;
  ((NearNeighborState *)state)->finishBucketSmooth(reqIDlist, owner);
}

void TreePiece::setupSmooth() {

  // XXX I don't believe any of the Chunks are used in the smooth walk.
  int oldNumChunks = numChunks;
  dm->getChunks(numChunks, prefetchRoots);
  CkArrayIndexMax idxMax = CkArrayIndex1D(thisIndex);
  bWalkDonePending = 0;
  streamingCache[CkMyPe()].cacheSync(numChunks, idxMax, localIndex);
  
  if (myNumParticles == 0) {
    // No particles assigned to this TreePiece
    for (int i=0; i< numChunks; ++i) streamingCache[CkMyPe()].finishedChunk(i, 0);
    contribute(0, 0, CkReduction::concat, callback);
    return;
  }
  
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

}

// Start tree walk and smooth calculation

void TreePiece::startIterationSmooth(// type of smoothing and parameters
				     SmoothParams* params,
				     const CkCallback& cb) {

  callback = cb;
  activeRung = params->activeRung;

  setupSmooth();

  for (int i=0; i<numChunks; ++i) remainingChunk[i] = myNumParticles;

  // Create objects that are reused by all buckets
  twSmooth = new BottomUpTreeWalk;
  sSmooth = new KNearestSmoothCompute(this, params, nSmooth);
      
  initBucketsSmooth((SmoothCompute *) sSmooth);

  // creates and initializes nearneighborstate object
  sSmoothState = sSmooth->getNewState(numBuckets);
  optSmooth = new SmoothOpt;

  completedActiveWalks = 0;	// XXX Potential race with Gravity Walk
  //smoothAwi = addActiveWalk(sTopDown,sSmooth,optSmooth,sSmoothState);
  smoothAwi = addActiveWalk(twSmooth,sSmooth,optSmooth,sSmoothState);
#ifdef CHECK_WALK_COMPLETIONS
  CkPrintf("[%d] addActiveWalk smooth (%d)\n", thisIndex, activeWalks.length());
#endif

  thisProxy[thisIndex].calculateSmoothLocal();
}

// Initialize particles in bucket and bucket bookkeeping
// Template for both SmoothCompute and ReSmoothCompute types.

template <class Tsmooth>
void TreePiece::initBucketsSmooth(Tsmooth tSmooth) {
  for (unsigned int j=0; j<numBuckets; ++j) {
    GenericTreeNode* node = bucketList[j];
    int numParticlesInBucket = node->particleCount;

    CkAssert(numParticlesInBucket <= TreeStuff::maxBucketSize);
    
    // TODO: active bounds may give a performance boost in the
    // multi-timstep regime.
    // node->boundingBox.reset();  // XXXX dangerous should have separate
				// Active bounds
    for(int i = node->firstParticle; i <= node->lastParticle; ++i) {
	if (myParticles[i].rung >= activeRung
	    && TYPETest(&myParticles[i], tSmooth->params->iType)) {
	    nActive++;
	    tSmooth->params->initSmoothParticle(&myParticles[i]);
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

void KNearestSmoothCompute::initSmoothPrioQueue(int iBucket, State *state) 
{
  // Prime the queues
  GenericTreeNode *myNode = tp->bucketList[iBucket];
  // state is passed in to function now. 
  NearNeighborState *nstate = (NearNeighborState *)state;
  
  //
  // Get nearest nSmooth particles in tree order
  //
  int firstQueue = myNode->firstParticle;
  int iCount = 0;
  int lastQueue = firstQueue;
  // Find only particles of interest
  for(lastQueue = firstQueue;
      iCount <= nSmooth && lastQueue <= tp->myNumParticles;
      lastQueue++) {
      if(TYPETest(&tp->myParticles[lastQueue], params->iType))
	  iCount++;
      }
  
  if(lastQueue > tp->myNumParticles) 
      {
	  lastQueue = tp->myNumParticles;
	  firstQueue = myNode->firstParticle - 1;
	  for(; iCount <= nSmooth; firstQueue--) {
	      if(TYPETest(&tp->myParticles[firstQueue], params->iType))
		  iCount++;
	      }
	  CkAssert(firstQueue > 0);
	  }
	  
  for(int j = myNode->firstParticle; j <= myNode->lastParticle; ++j) {
      if(!TYPETest(&tp->myParticles[j], params->iType))
	  continue;
      pqSmoothNode *Q = nstate->Qs[j]; 
      //
      // Find maximum of nearest neighbors
      //
      double drMax2 = 0.0;
      int kMax = 0;
      for(int k = firstQueue; k < lastQueue; k++) 
	  {
	      if(!TYPETest(&tp->myParticles[k], params->iType))
		  continue;
	      Vector3D<double> dr = tp->myParticles[k].position
		  - tp->myParticles[j].position;
	      if(dr.lengthSquared() > drMax2) {
		  kMax = k;
		  drMax2 = dr.lengthSquared();
		  }
	      }
      int end = nstate->heap_sizes[j];
      Q[end].fKey = sqrt(drMax2);
      Q[end].p = NULL; 
      std::push_heap(Q + 0, Q + end + 1); 
      nstate->heap_sizes[j]++; 
      }
    }

void TreePiece::smoothBucketComputation() {
  twSmooth->init(sSmooth, this);
  sSmooth->init(bucketList[currentBucket], activeRung, optSmooth);
  State *smoothState = activeWalks[smoothAwi].s;

  // start the tree walk from the tree built in the cache
  //  if (bucketList[currentBucket]->rungs >= activeRung) {
    for(int cr = 0; cr < numChunks; cr++){
      GenericTreeNode *chunkRoot = dm->chunkRootToNode(prefetchRoots[cr]);
      if(!chunkRoot){
        continue;
      }
      twSmooth->walk(chunkRoot, smoothState, cr,
		     encodeOffset(currentBucket, 0,0,0), smoothAwi);
      for(int x = -nReplicas; x <= nReplicas; x++) {
        for(int y = -nReplicas; y <= nReplicas; y++) {
          for(int z = -nReplicas; z <= nReplicas; z++) {
	      if(x || y || z)
		  twSmooth->walk(chunkRoot, smoothState, cr,
				 encodeOffset(currentBucket, x,y,z), smoothAwi);
          }
        }
      }
    }
  smoothState->counterArrays[0][currentBucket]--;
} 

void TreePiece::smoothNextBucket() {
  if(currentBucket >= numBuckets)
    return;
  State *smoothState = activeWalks[smoothAwi].s;
  smoothBucketComputation();
  ((NearNeighborState *)smoothState)->finishBucketSmooth(currentBucket, this);
}

void NearNeighborState::finishBucketSmooth(int iBucket, TreePiece *tp) {
  GenericTreeNode *node = tp->bucketList[iBucket];

  if(counterArrays[0][iBucket] == 0) {
      tp->sSmooth->walkDone(this);
  if(verbosity>1)
	CkPrintf("[%d] TreePiece %d finished smooth with bucket %d\n",CkMyPe(),
		 tp->thisIndex,iBucket);
    nParticlesPending -= node->particleCount;
    if(started && nParticlesPending == 0) {
      started = false;
      cacheManagerProxy[CkMyPe()].finishedChunk(0, 0);
#ifdef CHECK_WALK_COMPLETIONS
      CkPrintf("[%d] markWalkDone NearNeighborState\n", tp->getIndex());
#endif
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

void KNearestSmoothCompute::walkDone(State *state) {
  GenericTreeNode *node = (GenericTreeNode *) computeEntity;
  GravityParticle *part = node->particlePointer;

  for(int i = node->firstParticle; i <= node->lastParticle; i++) {
      if(!TYPETest(&part[i-node->firstParticle], params->iType))
	  continue;
      NearNeighborState *nstate = (NearNeighborState *)state;
      pqSmoothNode *Q = nstate->Qs[i];
      double h = Q[0].fKey; // Ball radius
      part[i-node->firstParticle].fBall = h;
      int nCnt = nstate->heap_sizes[i];
      params->fcnSmooth(&part[i-node->firstParticle], nCnt, Q);
      }
      // XXX jetley - the nearneighborstate allocated for this compute
      // may be deleted after this function is called. 
}

// From here down are "ReSmooth" methods.

// called after constructor, so tp should be set
State *ReSmoothCompute::getNewState(int nBucket){
  ReNearNeighborState *state = new ReNearNeighborState(tp->myNumParticles+2);
  //state->counterArrays.reserve(1);		   // array to keep track of
  state->counterArrays[0] = new int [nBucket];
  state->counterArrays[1] = 0;
  //state->counterArrays[0].reserve(tp->numBuckets); // outstanding requests
  for (int j = 0; j < nBucket; ++j) {
    state->counterArrays[0][j] = 1;	// so we know that the local
					// walk is not finished.
    }
  state->started = true;
  state->nParticlesPending = tp->myNumParticles;
  return state;
}

/*
 * Opening criterion for the reSmoothBucket walk.
 * Return true if we must open the node.
 */
int ReSmoothCompute::openCriterion(TreePiece *ownerTP, 
				  GenericTreeNode *node, // Node to test
				  int reqID, State *state) {
    GenericTreeNode *myNode = (GenericTreeNode *) computeEntity;
    GravityParticle *particles = ownerTP->getParticles();
    Vector3D<double> offset = ownerTP->decodeOffset(reqID);
    
    for(int j = myNode->firstParticle; j <= myNode->lastParticle; ++j) {
	if(!TYPETest(&particles[j], params->iType))
	    continue;
	double r = particles[j].fBall; // Ball radius
	Sphere<double> s(particles[j].position - offset, r);
	if(Space::intersect(node->boundingBox, s)) {
	   return 1;
	   }
	}
    return 0;
}

/*
 * Test a given particle against all the priority queues in the
 * bucket.
 */

void ReSmoothCompute::bucketCompare(TreePiece *ownerTP,
				  GravityParticle *p,  // Particle to test
				  GenericTreeNode *node, // bucket
				  GravityParticle *particles, // local
							      // particle data
				  Vector3D<double> offset,
                                  State *state
				  ) 
{
    ReNearNeighborState *nstate = (ReNearNeighborState *)state;
    Vector3D<double> rp = offset + p->position;
    for(int j = node->firstParticle; j <= node->lastParticle; ++j) {
	if(!TYPETest(&particles[j], params->iType))
	    continue;
	std::vector<pqSmoothNode> *Q = &nstate->Qs[j];
	double rOld = particles[j].fBall; // Ball radius
	Vector3D<double> dr = particles[j].position - rp;
	
	if(rOld*rOld >= dr.lengthSquared()) {  // Add to list
	    pqSmoothNode pqNew;
	    pqNew.fKey = dr.length();
	    pqNew.dx = dr;
	    pqNew.p = p;
	    Q->push_back(pqNew);
	    }
	}
    }

/*
 * Process particles received from missed Cache request
 */
void ReSmoothCompute::recvdParticlesFull(GravityParticle *part,
				   int num, int chunk,int reqID, State *state,
				   TreePiece *tp, Tree::NodeKey &remoteBucket){

  Vector3D<double> offset = tp->decodeOffset(reqID);
  int reqIDlist = decodeReqID(reqID);
  CkAssert(num > 0);
  state->counterArrays[0][reqIDlist] -= num;

  GenericTreeNode* reqnode = tp->bucketList[reqIDlist];

  for(int i=0;i<num;i++){
      if(!TYPETest(&part[i], params->iType))
	  continue;
      bucketCompare(tp, &part[i], reqnode, tp->myParticles, offset, state);
      }
  ((ReNearNeighborState *)state)->finishBucketSmooth(reqIDlist, tp);
}

int ReSmoothCompute::nodeRecvdEvent(TreePiece *owner, int chunk, State *state,
				  int reqIDlist){
  state->counterArrays[0][reqIDlist]--;
  ((ReNearNeighborState *)state)->finishBucketSmooth(reqIDlist, owner);
}

// Start tree walk and smooth calculation

void TreePiece::startIterationReSmooth(SmoothParams* params,
				       const CkCallback& cb) {

  callback = cb;
  activeRung = params->activeRung;

  setupSmooth();

  for (int i=0; i<numChunks; ++i) remainingChunk[i] = myNumParticles;

  // Create objects that are reused by all buckets
  twSmooth = new TopDownTreeWalk;
  sSmooth = new ReSmoothCompute(this, params);

  initBucketsSmooth((ReSmoothCompute *) sSmooth);

  // creates and initializes nearneighborstate object
  sSmoothState = sSmooth->getNewState(numBuckets);
  optSmooth = new SmoothOpt;

  completedActiveWalks = 0;	// XXX Potential race with Gravity Walk
  smoothAwi = addActiveWalk(twSmooth,sSmooth,optSmooth,sSmoothState);
#ifdef CHECK_WALK_COMPLETIONS
  CkPrintf("[%d] addActiveWalk reSmooth (%d)\n", thisIndex, activeWalks.length());
#endif

  thisProxy[thisIndex].calculateReSmoothLocal();
}

// Start the smoothing

void TreePiece::calculateReSmoothLocal() {
    dummyMsg *msg = new (8*sizeof(int)) dummyMsg;
    *((int *)CkPriorityPtr(msg)) = numTreePieces * numChunks + thisIndex + 1;
    CkSetQueueing(msg,CK_QUEUEING_IFIFO);
    msg->val=0;
    thisProxy[thisIndex].nextBucketReSmooth(msg);
    }

//
// Do the next set of buckets
//
void TreePiece::nextBucketReSmooth(dummyMsg *msg){
  unsigned int i=0;

  while(i<_yieldPeriod && currentBucket<numBuckets){
    reSmoothNextBucket();
    currentBucket++;
    i++;
  }

  if (currentBucket<numBuckets) {	// Queue up the next set
    thisProxy[thisIndex].nextBucketReSmooth(msg);
  } else {
    delete msg;
  }
}

void TreePiece::reSmoothNextBucket() {
  if(currentBucket >= numBuckets)
    return;
  State *smoothState = activeWalks[smoothAwi].s;
  smoothBucketComputation();
  ((ReNearNeighborState *)smoothState)->finishBucketSmooth(currentBucket, this);
}

void ReNearNeighborState::finishBucketSmooth(int iBucket, TreePiece *tp) {
  GenericTreeNode *node = tp->bucketList[iBucket];

  if(counterArrays[0][iBucket] == 0) {
      tp->sSmooth->walkDone(this);
  if(verbosity>1)
	CkPrintf("[%d] TreePiece %d finished smooth with bucket %d\n",CkMyPe(),
		 tp->thisIndex,iBucket);
    nParticlesPending -= node->particleCount;
    if(started && nParticlesPending == 0) {
      started = false;
      cacheManagerProxy[CkMyPe()].finishedChunk(0, 0);
#ifdef CHECK_WALK_COMPLETIONS
      CkPrintf("[%d] markWalkDone ReNearNeighborState\n", tp->getIndex());
#endif
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

void ReSmoothCompute::walkDone(State *state) {
  GenericTreeNode *node = (GenericTreeNode *) computeEntity;
  GravityParticle *part = node->particlePointer;

  for(int i = node->firstParticle; i <= node->lastParticle; i++) {
      if(!TYPETest(&part[i-node->firstParticle], params->iType))
	  continue;
      std::vector<pqSmoothNode> *Q = &((ReNearNeighborState *)state)->Qs[i];
      pqSmoothNode *NN = &((*Q)[0]);
      int nCnt = Q->size();
      params->fcnSmooth(&part[i-node->firstParticle], nCnt, NN);
      Q->clear();
      }
}

/*
 * Functions from PKDGRAV
 */
void Density(GravityParticle *p,int nSmooth, pqSmoothNode *nnList)
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

void DensitySmoothParams::initSmoothParticle(GravityParticle *p) 
{
    p->fDensity = 0.0;
    }

void DensitySmoothParams::initSmoothCache(GravityParticle *p) 
{
    p->fDensity = 0.0;
    }

void DensitySmoothParams::combSmoothCache(GravityParticle *p1,
					  ExternalSmoothParticle *p2) {
    p1->fDensity += p2->fDensity;
}

void DensitySmoothParams::fcnSmooth(GravityParticle *p,int nSmooth,
				    pqSmoothNode *nnList)
{
	GravityParticle *q;
	double fNorm,ih2,r2,rs;
	int i;

	ih2 = invH2(p);
	fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
	for (i=0;i<nSmooth;++i) {
		double fDist2 = nnList[i].fKey*nnList[i].fKey;
		r2 = fDist2*ih2;
		rs = KERNEL(r2);
		rs *= fNorm;
		q = nnList[i].p;
		p->fDensity += rs*q->mass;
		q->fDensity += rs*p->mass;
		}
	}
