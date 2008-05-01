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
 * Priority queues are kept in the "BucketSmoothRequest" structure.
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
 * 
 */

/*
 * Opening criterion for the smoothBucket walk.
 */
bool SmoothCompute::openCriterion(TreePiece *ownerTP, 
				  GenericTreeNode *node, // Node to test
				  int reqID) {
    GenericTreeNode *myNode = (GenericTreeNode *) computeEntity;
    GravityParticle *particles = ownerTP->getParticles();
    int listReqID = ownerTP->decodeReqID(reqID);
    Vector3D<double> offset = ownerTP->decodeOffset(reqID);
    BucketSmoothRequest *req = &ownerTP->bucketSmoothReqs[listReqID];
    
    for(int j = myNode->firstParticle; j <= myNode->lastParticle; ++j) {
	// XXX Qs shouldn't be stored in reqs; Probably needs to go in
	// TreePiece somewhere.
	double r = req->Qs[j-myNode->firstParticle].top().fKey; // Ball radius
	Sphere s(particles[j]->pos - offset, r);
	if(Space::intersect(node->boundingBox, s)
	   return KEEP;
	   }
	}
    return DUMP;
}

void SmoothCompute::bucketCompare(ExternalGravityParticle *p,  // Particle to test
				  GenericTreeNode node, // bucket
				  GravityParticle *particles, // local
							      // particle data
				  Vector3D<double> offset,
				  SmoothBucketReq *req	// request data
				  ) 
{
    for(int j = node->firstParticle; j <= node->lastParticle; ++j) {
	priority_queue<pqSmoothNode> *Q = &req->Qs[j-node->firstParticle];
	double rOld = Q->top().fKey; // Ball radius
	Vector3d<double> dr = offset + part->position - particles[j].position;
	
	if(rOld*rOld > dr.lengthSquared()) {  // Perform replacement
	    Q->pop();
	    pqSmoothNode pqNew;
	    pqNew.fKey = dr.length();
	    pqNew.dx = dr;
	    pqNew.p = part;
	    Q->push(pqNew);
	    }
	}
    }

int SmoothCompute::doWork(GenericTreeNode *node,
			  TreeWalk *tw,
			  State *state,
			  int chunk,
			  bool isRoot, 
			  bool &didcomp)
{
    
    // so that we have a quick return in case of empty nodes
    if(node->getType() == Empty || node->getType() == CachedEmpty){
	return DUMP;
	}
    // check opening criterion
    bool open = openCriterion(tp, node, reqID);
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
	    bucketCompare(&part[i-node->firstParticle],
			  (GenericTreeNode *) computeEntity,
			  tp->getParticles(),
			  tp->decodeOffset(reqID),
			  &tp->bucketSmoothReqs[tp->decodeReqID(reqID)]);
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
				    reqID, false);
	if(part) {
	    // Particles available; Search for contained.
	    for(int i = node->firstParticle; i <= node->lastParticle; i++) {
		bucketCompare(&part[i-node->firstParticle],
			      (GenericTreeNode *) computeEntity,
			      tp->getParticles(),
			      tp->decodeOffset(reqID),
			      &tp->bucketSmoothReqs[tp->decodeReqID(reqID)]);
		}
	    }
	else {
	    // Missed cache; record this
	    if(getOptType() == Remote){
		tp->addToRemainingChunk(chunk, node->lastParticle
					- node->firstParticle+1);
		}
	    }
	
	return DUMP;
	}
    else if(action == DUMP) {
	return DUMP;
	}
    // Check for no more requests here and call smooth function.
    }


void TreePiece::smoothNextBucket() {
  if(currentBucket >= numBuckets)
    return;

  TreeWalk *tw = new TopDownTreeWalk();
  Compute *gc = new SmoothCompute();
  State *nullstate = new NullState();
  Opt *opt = new LocalOpt();

  gc->init(bucketList[currentBucket], activeRung, opt);
  tw->init(gc, this);

  // start the tree walk from the tree built in the cache
  if (bucketList[currentBucket]->rungs >= activeRung) {
    for(int cr = 0; cr < numChunks; cr++){
#ifdef CACHE_TREE
            GenericTreeNode *chunkRoot = localCache->chunkRootToNode(prefetchRoots[cr]);
#else
            GenericTreeNode *chunkRoot = keyToNode(prefetchRoots[cr]);
#endif
      for(int x = -nReplicas; x <= nReplicas; x++) {
        for(int y = -nReplicas; y <= nReplicas; y++) {
          for(int z = -nReplicas; z <= nReplicas; z++) {
            tw->walk(chunkRoot, nullstate, -1, encodeOffset(currentBucket, x,y,z));
          }
        }
      }
    }
  }
  bucketReqs[currentBucket].numAdditionalRequests --;
  finishBucket(currentBucket);
  delete gc;
  delete nullstate; 
  delete opt;
  delete tw;
}

// This will go away.  We'll do the function call elsewhere.
void Compute::finishSmooth(int iBucket) {
  BucketSmoothRequest *req = &bucketReqs[iBucket];
  GenericTreeNode *node = bucketList[iBucket];
  GravityParticle *part = node->particlePointer;

  for(int i = node->firstParticle; i <= node->lastParticle; i++) {
      priority_queue<pqSmoothNode> *Q = &req->Qs[i-node->firstParticle];
      double h = Q->top().fKey; // Ball radius
      part[i].fBall = h;
      pqSmoothNode NN[nSmooth];
      int nCnt = 0;
      while (!Q->empty()) {
	  NN[nCnt] = Q->top();
	  Q->pop();
	  nCnt++;
	  }
      fcnSmooth(&part[i], nCnt, NN);
      }
}

/* Standard M_4 Kernel */
/* return 1/(h_smooth)^2 for a particle */
inline double invH2(SmoothParticle &p) 
{
    return 4.0/p.fBall2;
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
    }
/*
 * XXX Place holder from PKDGRAV
 */
void Density(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	FLOAT ih2,r2,rs,fDensity;
	int i;

	ih2 = invH2(p);
	fDensity = 0.0;
	for (i=0;i<nSmooth;++i) {
		r2 = nnList[i].fDist2*ih2;
		rs = KERNEL(r2);
		fDensity += rs*nnList[i].pPart->fMass;
		}
	p->fDensity = M_1_PI*sqrt(ih2)*ih2*fDensity; 
	}
