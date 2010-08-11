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
    // so that we have a quick return in case of empty nodes
    if(node->getType() == Empty || node->getType() == CachedEmpty){
	return DUMP;
	}
    if(!(node->iParticleTypes & params->iType)) {
	return DUMP; // no particles of appropriate type
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
	    // particle is part of search tree.
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
		// particle is part of search tree.
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
    CkAbort("SmoothCompute: bad walk");
    return -1;
    }

/*
 * reassociate bucket
 */
void SmoothCompute::reassoc(void *ce, int ar, Opt *o){
  computeEntity = ce; 
  activeRung = ar;
  opt = o;
}

void SmoothCompute::nodeMissedEvent(int reqID, int chunk, State *state, TreePiece *tp)
{
    int reqIDlist = decodeReqID(reqID);
    state->counterArrays[0][reqIDlist]++;
    }

// called after constructor, so tp should be set
State *KNearestSmoothCompute::getNewState(int nBuckets){
  NearNeighborState *state = new NearNeighborState(tp->myNumParticles+2, nSmooth);
  // array to keep track of outstanding requests
  state->counterArrays[0] = new int [nBuckets];
  state->counterArrays[1] = 0;

  for(int j = 0; j < nBuckets; ++j){
    initSmoothPrioQueue(j, state);
    state->counterArrays[0][j] = 1;	// so we know that the local
					// walk is not finished.
  }
  state->started = true;
  state->nParticlesPending = tp->myNumParticles;
  state->currentBucket = 0;
  state->bWalkDonePending = 0;
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
	if(!params->isSmoothActive(&particles[j]))
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
	if(!params->isSmoothActive(&particles[j]))
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

void KNearestSmoothCompute::nodeRecvdEvent(TreePiece *owner, int chunk, State *state,
				  int reqIDlist){
  state->counterArrays[0][reqIDlist]--;
  ((NearNeighborState *)state)->finishBucketSmooth(reqIDlist, owner);
}

void TreePiece::setupSmooth() {

  // XXX I don't believe any of the Chunks are used in the smooth walk.
  dm->getChunks(numChunks, prefetchRoots);
  CkArrayIndexMax idxMax = CkArrayIndex1D(thisIndex);
  if (numChunks == 0 && myNumParticles == 0) numChunks = 1;
  cacheSmoothPart[CkMyPe()].cacheSync(numChunks, idxMax, localIndex);
  cacheNode[CkMyPe()].cacheSync(numChunks, idxMax, localIndex);
  
  // The following if is necessary to prevent nodes containing only TreePieces
  // without particles from getting stuck and crashing.
  if (myNumParticles == 0) {
    // No particles assigned to this TreePiece
      for (int i=0; i< numChunks; ++i) {
	  cacheNode[CkMyPe()].finishedChunk(i, 0);
	  cacheSmoothPart[CkMyPe()].finishedChunk(i, 0);
	  }
      return;
  }
  
}

// Start tree walk and smooth calculation

void TreePiece::startIterationSmooth(// type of smoothing and parameters
				     SmoothParams* params,
				     const CkCallback& cb) {

  cbSmooth = cb;
  activeRung = params->activeRung;

  setupSmooth();
  if (myNumParticles == 0) {
      markSmoothWalkDone();
      return;
      }

  // It is assumed that there are at least nSmooth particles on a
  // treepiece in order to prime the search queue.
  CkAssert(myNumParticles > nSmooth);
  
  //for (int i=0; i<numChunks; ++i) remaining Chunk[i] = myNumParticles;

  // Create objects that are reused by all buckets
  twSmooth = new BottomUpTreeWalk;
  sSmooth = new KNearestSmoothCompute(this, params, nSmooth);
      
  initBucketsSmooth((SmoothCompute *) sSmooth);

  // creates and initializes nearneighborstate object
  sSmoothState = sSmooth->getNewState(numBuckets);
  optSmooth = new SmoothOpt;
  processReqSmoothParticles();

  activeWalks.reserve(maxAwi);
  addActiveWalk(smoothAwi, twSmooth,sSmooth,optSmooth,sSmoothState);
#ifdef CHECK_WALK_COMPLETIONS
  CkPrintf("[%d] addActiveWalk smooth (%d)\n", thisIndex, activeWalks.length());
#endif

  thisProxy[thisIndex].calculateSmoothLocal();
}

// Initialize particles in bucket and bucket bookkeeping
// Template for both SmoothCompute and ReSmoothCompute types.

template <class Tsmooth>
void TreePiece::initBucketsSmooth(Tsmooth tSmooth) {
  tSmooth->nActive = 0;
  for (unsigned int j=0; j<numBuckets; ++j) {
    GenericTreeNode* node = bucketList[j];
    int numParticlesInBucket = node->particleCount;

    CkAssert(numParticlesInBucket <= TreeStuff::maxBucketSize);
    
    // TODO: active bounds may give a performance boost in the
    // multi-timstep regime.
    // node->boundingBox.reset();  // XXXX dangerous should have separate
				// Active bounds
    for(int i = node->firstParticle; i <= node->lastParticle; ++i) {
	if(tSmooth->params->isSmoothActive(&myParticles[i])) {
	    tSmooth->nActive++;
	    tSmooth->params->initSmoothParticle(&myParticles[i]);
	    // node->boundingBox.grow(myParticles[i].position);
	    }
	else if (TYPETest(&myParticles[i], tSmooth->params->iType)) {
	    tSmooth->params->initTreeParticle(&myParticles[i]);
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

  int currentBucket = sSmoothState->currentBucket;
  
  while(i<_yieldPeriod && currentBucket<numBuckets){
    smoothNextBucket();
    currentBucket++;
    sSmoothState->currentBucket++;
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
  // First search from the start of the bucket to the end of this treepiece
  for(lastQueue = firstQueue;
      iCount <= nSmooth && lastQueue <= tp->myNumParticles;
      lastQueue++) {
      if(TYPETest(&tp->myParticles[lastQueue], params->iType))
	  iCount++;
      }
  
  // We still don't have enough particles.  Search back to the start of the
  // treepiece
  int bEnough = 1;	// Do we have enough particles on piece to get a limit?
  if(lastQueue > tp->myNumParticles) 
      {
	  lastQueue = tp->myNumParticles;
	  firstQueue = myNode->firstParticle - 1;
	  for(; iCount <= nSmooth; firstQueue--) {
	      if(firstQueue == 0) {
		  bEnough = 0; // Ran out of particles
		  break;
		  }
	      if(TYPETest(&tp->myParticles[firstQueue], params->iType))
		  iCount++;
	      }
	  // CkAssert(firstQueue > 0);
	  }
	  
  for(int j = myNode->firstParticle; j <= myNode->lastParticle; ++j) {
      if(!params->isSmoothActive(&tp->myParticles[j]))
	  continue;
      pqSmoothNode *Q = nstate->Qs[j] = new pqSmoothNode[nSmooth];
      //
      // Find maximum of nearest neighbors
      //
      double drMax2 = 0.0;
      for(int k = firstQueue; k < lastQueue; k++) 
	  {
	      if(!TYPETest(&tp->myParticles[k], params->iType))
		  continue;
	      Vector3D<double> dr = tp->myParticles[k].position
		  - tp->myParticles[j].position;
	      if(dr.lengthSquared() > drMax2) {
		  drMax2 = dr.lengthSquared();
		  }
	      }
      int end = nstate->heap_sizes[j];
      if(bEnough)
	  Q[end].fKey = sqrt(drMax2);
      else
	  Q[end].fKey = HUGE;
      Q[end].p = NULL; 
      std::push_heap(Q + 0, Q + end + 1); 
      nstate->heap_sizes[j]++; 
      }
    }

void TreePiece::smoothBucketComputation() {
  twSmooth->init(sSmooth, this);
  int currentBucket = sSmoothState->currentBucket;
  sSmooth->init(bucketList[currentBucket], activeRung, optSmooth);

  // start the tree walk from the tree built in the cache
  //  if (bucketList[currentBucket]->rungs >= activeRung) {
    for(int cr = 0; cr < numChunks; cr++){
      GenericTreeNode *chunkRoot = dm->chunkRootToNode(prefetchRoots[cr]);
      if(!chunkRoot){
        continue;
      }
      twSmooth->walk(chunkRoot, sSmoothState, cr,
		     encodeOffset(currentBucket, 0,0,0), smoothAwi);
      for(int x = -nReplicas; x <= nReplicas; x++) {
        for(int y = -nReplicas; y <= nReplicas; y++) {
          for(int z = -nReplicas; z <= nReplicas; z++) {
	      if(x || y || z)
		  twSmooth->walk(chunkRoot, sSmoothState, cr,
				 encodeOffset(currentBucket, x,y,z), smoothAwi);
          }
        }
      }
    }
  sSmoothState->counterArrays[0][currentBucket]--;
} 

void TreePiece::smoothNextBucket() {
  int currentBucket = sSmoothState->currentBucket;
  if(currentBucket >= numBuckets)
    return;
  smoothBucketComputation();
  ((NearNeighborState *)sSmoothState)->finishBucketSmooth(currentBucket, this);
}

void NearNeighborState::finishBucketSmooth(int iBucket, TreePiece *tp) {
  GenericTreeNode *node = tp->bucketList[iBucket];

  if(counterArrays[0][iBucket] == 0) {
    tp->sSmooth->walkDone(this);
    if(verbosity>4)
	CkPrintf("[%d] TreePiece %d finished smooth with bucket %d\n",CkMyPe(),
		 tp->thisIndex,iBucket);
    nParticlesPending -= node->particleCount;
    if(started && nParticlesPending == 0) {
      started = false;
      cacheNode[CkMyPe()].finishedChunk(0, 0);
      cacheSmoothPart[CkMyPe()].finishedChunk(0, 0);
#ifdef CHECK_WALK_COMPLETIONS
      CkPrintf("[%d] markWalkDone NearNeighborState\n", tp->getIndex());
#endif
      tp->markSmoothWalkDone();
      if(verbosity > 1)
	ckerr << "TreePiece " << tp->thisIndex << ": My particles are done"
	     << endl;
    }
  }
}

void TreePiece::markSmoothWalkDone() 
{
	CkCallback cb = CkCallback(CkIndex_TreePiece::finishSmoothWalk(),
				   pieces);
	// Use shadow array to avoid reduction conflict
	smoothProxy[thisIndex].ckLocal()->contribute(0, 0,
						     CkReduction::concat, cb);
    }

void TreePiece::finishSmoothWalk()
{
  // Furthermore, we need to wait for outstanding flushes to be
  // processed for the combiner cache.
  if(sSmooth && nCacheAccesses > 0) {
    sSmoothState->bWalkDonePending = 1;
    return;
  }

  nCacheAccesses = 0; // reset for next walk.

  if(myNumParticles != 0) {
      sSmooth->freeState(sSmoothState);
      delete sSmooth;
      delete optSmooth;
      delete twSmooth;
      sSmooth = 0;
      }
#ifdef CHECK_WALK_COMPLETIONS
  CkPrintf("[%d] inside finishSmoothWalk contrib callback\n", thisIndex);
#endif
  smoothProxy[thisIndex].ckLocal()->contribute(0, 0, CkReduction::concat,
					       cbSmooth);
}

void KNearestSmoothCompute::walkDone(State *state) {
  GenericTreeNode *node = (GenericTreeNode *) computeEntity;
  GravityParticle *part = node->particlePointer;

  for(int i = node->firstParticle; i <= node->lastParticle; i++) {
      if(!params->isSmoothActive(&part[i-node->firstParticle]))
	  continue;
      NearNeighborState *nstate = (NearNeighborState *)state;
      pqSmoothNode *Q = nstate->Qs[i];
      double h = Q[0].fKey; // Ball radius
      part[i-node->firstParticle].fBall = h;
      int nCnt = nstate->heap_sizes[i];
      params->fcnSmooth(&part[i-node->firstParticle], nCnt, Q);
      delete [] Q;
      nstate->Qs[i] = NULL;
      }
      // XXX jetley - the nearneighborstate allocated for this compute
      // may be deleted after this function is called. 
}

// From here down are "ReSmooth" methods.

// called after constructor, so tp should be set
State *ReSmoothCompute::getNewState(int nBucket){
  ReNearNeighborState *state = new ReNearNeighborState(tp->myNumParticles+2);
  // array to keep track of outstanding requests
  state->counterArrays[0] = new int [nBucket];
  state->counterArrays[1] = 0;
  for (int j = 0; j < nBucket; ++j) {
    state->counterArrays[0][j] = 1;	// so we know that the local
					// walk is not finished.
    }
  state->started = true;
  state->nParticlesPending = tp->myNumParticles;
  state->currentBucket = 0;
  state->bWalkDonePending = 0;
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
	if(!params->isSmoothActive(&particles[j]))
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
	if(!params->isSmoothActive(&particles[j]))
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

void ReSmoothCompute::nodeRecvdEvent(TreePiece *owner, int chunk, State *state,
				  int reqIDlist){
  state->counterArrays[0][reqIDlist]--;
  ((ReNearNeighborState *)state)->finishBucketSmooth(reqIDlist, owner);
}

// Start tree walk and smooth calculation

void TreePiece::startIterationReSmooth(SmoothParams* params,
				       const CkCallback& cb) {

  cbSmooth = cb;
  activeRung = params->activeRung;

  setupSmooth();
  if (myNumParticles == 0) {
      markSmoothWalkDone();
      return;
      }

  //for (int i=0; i<numChunks; ++i) remaining Chunk[i] = myNumParticles;

  // Create objects that are reused by all buckets
  twSmooth = new TopDownTreeWalk;
  sSmooth = new ReSmoothCompute(this, params);

  initBucketsSmooth((ReSmoothCompute *) sSmooth);

  // creates and initializes nearneighborstate object
  sSmoothState = sSmooth->getNewState(numBuckets);
  optSmooth = new SmoothOpt;
  processReqSmoothParticles();

  activeWalks.reserve(maxAwi);
  addActiveWalk(smoothAwi, twSmooth,sSmooth,optSmooth,sSmoothState);
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
  int currentBucket = sSmoothState->currentBucket;
  
  while(i<_yieldPeriod && currentBucket<numBuckets){
    reSmoothNextBucket();
    currentBucket++;
    sSmoothState->currentBucket++;
    i++;
  }

  if (currentBucket<numBuckets) {	// Queue up the next set
    thisProxy[thisIndex].nextBucketReSmooth(msg);
  } else {
    delete msg;
  }
}

void TreePiece::reSmoothNextBucket() {
  int currentBucket = sSmoothState->currentBucket;
  if(currentBucket >= numBuckets)
    return;
  smoothBucketComputation();
  ((ReNearNeighborState *)sSmoothState)->finishBucketSmooth(currentBucket, this);
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
      cacheNode[CkMyPe()].finishedChunk(0, 0);
      cacheSmoothPart[CkMyPe()].finishedChunk(0, 0);
#ifdef CHECK_WALK_COMPLETIONS
      CkPrintf("[%d] markWalkDone ReNearNeighborState\n", tp->getIndex());
#endif
      tp->markSmoothWalkDone();
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
      if(!params->isSmoothActive(&part[i-node->firstParticle]))
	  continue;
      std::vector<pqSmoothNode> *Q = &((ReNearNeighborState *)state)->Qs[i];
      pqSmoothNode *NN = &((*Q)[0]);
      int nCnt = Q->size();
      params->fcnSmooth(&part[i-node->firstParticle], nCnt, NN);
      Q->clear();
      }
}

// called after constructor, so tp should be set
State *MarkSmoothCompute::getNewState(int nBucket){
  MarkNeighborState *state = new MarkNeighborState(tp->myNumParticles+2);
  state->counterArrays[0] = new int [nBucket];
  state->counterArrays[1] = 0;
  for (int j = 0; j < nBucket; ++j) {
    state->counterArrays[0][j] = 1;	// so we know that the local
					// walk is not finished.
    }
  state->started = true;
  state->nParticlesPending = tp->myNumParticles;
  state->currentBucket = 0;
  state->bWalkDonePending = 0;
  return state;
}

/*
 * Opening criterion for the MarkSmoothBucket walk.
 * Return true if we must open the node.
 */
int MarkSmoothCompute::openCriterion(TreePiece *ownerTP, 
				  GenericTreeNode *node, // Node to test
				  int reqID, State *state) {
    GenericTreeNode *myNode = (GenericTreeNode *) computeEntity; // my bucket
    GravityParticle *particles = ownerTP->getParticles();
    Vector3D<double> offset = ownerTP->decodeOffset(reqID);
    
    for(int j = myNode->firstParticle; j <= myNode->lastParticle; ++j) {
	if(!params->isSmoothActive(&particles[j]))
	    continue;
	Vector3D<double> p(particles[j].position - offset);
	if(Space::contains(node->bndBoxBall, p)) {
	   return 1;
	   }
	}
    return 0;
}

/*
 * Test a given particle against all the priority queues in the
 * bucket.
 */

void MarkSmoothCompute::bucketCompare(TreePiece *ownerTP,
				  GravityParticle *p,  // Particle to test
				  GenericTreeNode *node, // bucket
				  GravityParticle *particles, // local
							      // particle data
				  Vector3D<double> offset,
                                  State *state
				  ) 
{
    Vector3D<double> rp = offset + p->position;
    for(int j = node->firstParticle; j <= node->lastParticle; ++j) {
	if(!params->isSmoothActive(&particles[j]))
	    continue;
	double rOld = p->fBallMax(); // Ball radius
	Vector3D<double> dr = particles[j].position - rp;
	
	if(rOld*rOld >= dr.lengthSquared()) {  // Mark as Neighbor
	    TYPESet(p, TYPE_NbrOfACTIVE);
	    }
	}
    }

/*
 * Process particles received from missed Cache request
 */
void MarkSmoothCompute::recvdParticlesFull(GravityParticle *part,
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
  ((MarkNeighborState *)state)->finishBucketSmooth(reqIDlist, tp);
}

void MarkSmoothCompute::nodeRecvdEvent(TreePiece *owner, int chunk,
				      State *state, int reqIDlist){
  state->counterArrays[0][reqIDlist]--;
  ((MarkNeighborState *)state)->finishBucketSmooth(reqIDlist, owner);
}

// Start marksmooth walk

void TreePiece::startIterationMarkSmooth(SmoothParams* params,
				       const CkCallback& cb) {

  cbSmooth = cb;
  activeRung = params->activeRung;

  setupSmooth();
  if (myNumParticles == 0) {
      markSmoothWalkDone();
      return;
      }

  //for (int i=0; i<numChunks; ++i) remaining Chunk[i] = myNumParticles;

  // Create objects that are reused by all buckets
  twSmooth = new TopDownTreeWalk;
  sSmooth = new MarkSmoothCompute(this, params);

  initBucketsSmooth((MarkSmoothCompute *) sSmooth);

  // creates and initializes nearneighborstate object
  sSmoothState = sSmooth->getNewState(numBuckets);
  optSmooth = new SmoothOpt;
  processReqSmoothParticles();

  addActiveWalk(smoothAwi, twSmooth,sSmooth,optSmooth,sSmoothState);
#ifdef CHECK_WALK_COMPLETIONS
  CkPrintf("[%d] addActiveWalk reSmooth (%d)\n", thisIndex, activeWalks.length());
#endif

  thisProxy[thisIndex].calculateMarkSmoothLocal();
}

// Start the smoothing

void TreePiece::calculateMarkSmoothLocal() {
    dummyMsg *msg = new (8*sizeof(int)) dummyMsg;
    *((int *)CkPriorityPtr(msg)) = numTreePieces * numChunks + thisIndex + 1;
    CkSetQueueing(msg,CK_QUEUEING_IFIFO);
    msg->val=0;
    thisProxy[thisIndex].nextBucketMarkSmooth(msg);
    }

//
// Do the next set of buckets
//
void TreePiece::nextBucketMarkSmooth(dummyMsg *msg){
  unsigned int i=0;
  int currentBucket = sSmoothState->currentBucket;
  
  while(i<_yieldPeriod && currentBucket<numBuckets){
    markSmoothNextBucket();
    currentBucket++;
    sSmoothState->currentBucket++;
    i++;
  }

  if (currentBucket<numBuckets) {	// Queue up the next set
    thisProxy[thisIndex].nextBucketMarkSmooth(msg);
  } else {
    delete msg;
  }
}

void TreePiece::markSmoothNextBucket() {
  int currentBucket = sSmoothState->currentBucket;
  if(currentBucket >= numBuckets)
    return;
  smoothBucketComputation();
  ((MarkNeighborState *)sSmoothState)->finishBucketSmooth(currentBucket, this);
}

void MarkNeighborState::finishBucketSmooth(int iBucket, TreePiece *tp) {
  GenericTreeNode *node = tp->bucketList[iBucket];

  if(counterArrays[0][iBucket] == 0) {
      tp->sSmooth->walkDone(this);
  if(verbosity>1)
	CkPrintf("[%d] TreePiece %d finished smooth with bucket %d\n",CkMyPe(),
		 tp->thisIndex,iBucket);
    nParticlesPending -= node->particleCount;
    if(started && nParticlesPending == 0) {
      started = false;
      cacheNode[CkMyPe()].finishedChunk(0, 0);
      cacheSmoothPart[CkMyPe()].finishedChunk(0, 0);
#ifdef CHECK_WALK_COMPLETIONS
      CkPrintf("[%d] markWalkDone ReNearNeighborState\n", tp->getIndex());
#endif
      tp->markSmoothWalkDone();
      if(verbosity>1)
	CkPrintf("[%d] TreePiece %d finished smooth with bucket %d\n",CkMyPe(),
		 tp->thisIndex,iBucket);
      if(verbosity > 4)
	ckerr << "TreePiece " << tp->thisIndex << ": My particles are done"
	     << endl;
    }
  }
}

void MarkSmoothCompute::walkDone(State *state) {
    // Nothing to do
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

int DensitySmoothParams::isSmoothActive(GravityParticle *p) 
{
    return (p->rung >= activeRung && TYPETest(p, iType));
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
