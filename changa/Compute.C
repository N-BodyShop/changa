#include "ParallelGravity.h"
#include "GenericTreeNode.h"
//#include "codes.h"

#include "Opt.h"
#include "Compute.h"
#include "TreeWalk.h"
//#include "State.h"
#include "Space.h"
#include "gravity.h" 

extern int decodeReqID(int reqID);

void Compute::setOpt(Opt *_opt){
  opt = _opt;
}

OptType Compute::getOptType(){
  return opt->getSelfType();
}

void Compute::init(void *buck, int ar, Opt *o){
  computeEntity = buck;
  activeRung = ar;
  opt = o;
}

State *Compute::getNewState(int dim1, int dim2){
  State *s = new State();
  // 2 arrays of counters 
  // 0. numAdditionalRequests[] - sized numBuckets, init to numChunks
  // 1. remainingChunk[] - sized numChunks
  s->counterArrays[0] = new int [dim1];
  s->counterArrays[1] = new int [dim2];

  return s;
}

State *Compute::getNewState(int dim1){
  // 0. local component of numAdditionalRequests, init to 1
  State *s = new State();
  s->counterArrays[0] = new int [dim1];
  s->counterArrays[1] = 0;
  return s;
}

State *Compute::getNewState(){
  return 0;
}

void Compute::freeState(State *s){
  if(s->counterArrays[0]){
    delete [] s->counterArrays[0];
    s->counterArrays[0] = 0;
  }
  if(s->counterArrays[1]){
    delete [] s->counterArrays[1];
    s->counterArrays[1] = 0;
  }
  delete s;
}

#if INTERLIST_VER > 0
void ListCompute::freeState(State *s){
  freeDoubleWalkState((DoubleWalkState *)s); 
  Compute::freeState(s);
}

void ListCompute::freeDoubleWalkState(DoubleWalkState *state){
  for(int i = 0; i < INTERLIST_LEVELS; i++){
    state->undlists[i].free();
    state->clists[i].free();
  }
  delete [] state->chklists;
  state->undlists.free();
  state->clists.free();

  if(getOptType() == Remote){
    for(int i = 0; i < INTERLIST_LEVELS; i++){
      state->rplists[i].free();
    }
    state->rplists.free();
  }
  else if(getOptType() == Local){
    for(int i = 0; i < INTERLIST_LEVELS; i++){
      state->lplists[i].free();
    }
    state->lplists.free();
  }

}

DoubleWalkState *ListCompute::allocDoubleWalkState(){
  DoubleWalkState *s = new DoubleWalkState;
  s->level = 0;
  s->chklists = new CheckList[INTERLIST_LEVELS]; 
  s->undlists.resize(INTERLIST_LEVELS);
  s->clists.reserve(INTERLIST_LEVELS);

  if(getOptType() == Remote){
    s->rplists.resize(INTERLIST_LEVELS);
  }
  else if(getOptType() == Local){
    s->lplists.resize(INTERLIST_LEVELS);
  }
  
  return s;
}

State *ListCompute::getNewState(int d1, int d2){
  DoubleWalkState *s = allocDoubleWalkState();
  s->counterArrays[0] = new int [d1];
  s->counterArrays[1] = new int [d2];
  return s;
}

State *ListCompute::getNewState(int d1){
  DoubleWalkState *s = allocDoubleWalkState();
  s->counterArrays[0] = new int [d1];
  s->counterArrays[1] = 0;
  return s;
}

State *ListCompute::getNewState(){
  DoubleWalkState *s = allocDoubleWalkState();
  s->counterArrays[0] = 0;
  s->counterArrays[1] = 0;
  return s;
}

void ListCompute::initState(State *state){
  
  DoubleWalkState *s = (DoubleWalkState *)state;
  int level = s->level;
  UndecidedList &myUndlist = s->undlists[level];

  // *my* undecided list:
  myUndlist.length() = 0;
  // interaction lists:
  s->clists[level].length() = 0;
  if(getOptType() == Local){
    s->lplists[level].length() = 0;
  }
  else if(getOptType() == Remote){
    s->rplists[level].length() = 0;
  }
  else{
    CkAbort("Invalid Opt type for ListCompute");
  }
}
#endif

void GravityCompute::reassoc(void *ce, int ar, Opt *o){
  computeEntity = ce; 
  activeRung = ar;
  opt = o;
}

#if INTERLIST_VER > 0
void ListCompute::reassoc(void *ce, int ar, Opt *o){
  computeEntity = ce; 
  activeRung = ar;
  opt = o;
}
#endif
/*
void PrefetchCompute::init(void *buck, int ar, Opt *o){

  // buck is actually the TreePiece this PrefetchCompute
  // will query when it needs a prefetchRequestStruct
  computeEntity = buck;
  activeRung = ar;
  opt = o;
}
*/
int GravityCompute::nodeMissedEvent(int reqID, int chunk, State *state, TreePiece *tp){
  if(getOptType() == Remote){
    state->counterArrays[0][decodeReqID(reqID)]++;
    state->counterArrays[1][chunk]++;
  }
}

int PrefetchCompute::startNodeProcessEvent(State *state){
  //return owner->incPrefetchWaiting();
  state->counterArrays[0][0]++;
  return state->counterArrays[0][0];
}

int PrefetchCompute::finishNodeProcessEvent(TreePiece *owner, State *state){
  //int save = owner->decPrefetchWaiting();
  int save = --state->counterArrays[0][0];
  if(save == 0){
    owner->startRemoteChunk();
  }    
  return save;
}

#if INTERLIST_VER > 0
int ListCompute::nodeMissedEvent(int reqID, int chunk, State *state, TreePiece *tp){
  if(getOptType() == Remote){
#ifdef CHANGA_REFACTOR_MEMCHECK
    CkPrintf("memcheck before nodemissed\n");
    CmiMemoryCheck();
#endif
    int startBucket;
    int end;
    GenericTreeNode *source = (GenericTreeNode *)computeEntity;
    tp->getBucketsBeneathBounds(source, startBucket, end);
    tp->updateUnfinishedBucketState(startBucket, end, 1, chunk, state);
#ifdef CHANGA_REFACTOR_MEMCHECK
    CkPrintf("memcheck after nodemissed\n");
    CmiMemoryCheck();
#endif
  }
}
#endif

int GravityCompute::openCriterion(TreePiece *ownerTP, 
                          GenericTreeNode *node, int reqID, State *state){
  return 
    openCriterionBucket(node,(GenericTreeNode *)computeEntity,ownerTP->decodeOffset(reqID), ownerTP->getIndex());

}

extern int partBucketForce(ExternalGravityParticle *part, GenericTreeNode *req, GravityParticle *particles, Vector3D<double> offset, int activeRung);

void GravityCompute::recvdParticles(ExternalGravityParticle *part,int num,int chunk,int reqID,State *state,TreePiece *tp, Tree::NodeKey &remoteBucket){
  //TreePiece *tp = tw->getOwnerTP();

  Vector3D<double> offset = tp->decodeOffset(reqID);
  int reqIDlist = decodeReqID(reqID);
  CkAssert(num > 0);
  //tp->bucketReqs[reqIDlist].numAdditionalRequests -= num;
  state->counterArrays[0][reqIDlist] -= 1;
  //tp->remainingChunk[chunk] -= num;
  state->counterArrays[1][chunk] -= 1;

  GenericTreeNode* reqnode = tp->bucketList[reqIDlist];

  int computed;
#ifdef BENCHMARK_TIME_COMPUTE
  double startTime = CmiWallTimer();
#endif
  for(int i=0;i<num;i++){
    
#if COSMO_STATS > 1
    for(int j = reqnode->firstParticle; j <= reqnode->lastParticle; ++j) {
      tp->myParticles[j].extpartmass += part[i].mass;
    }
#endif
    
#ifdef COSMO_EVENTS
    double startTimer = CmiWallTimer();
#endif
#ifdef HPM_COUNTER
    hpmStart(2,"particle force");
#endif
    computed = partBucketForce(&part[i], reqnode, tp->myParticles, offset, activeRung);
#ifdef HPM_COUNTER
    hpmStop(2);
#endif
#ifdef COSMO_EVENTS
    traceUserBracketEvent(partForceUE, startTimer, CmiWallTimer());
#endif
  }
#ifdef BENCHMARK_TIME_COMPUTE
  computeTimePart += CmiWallTimer() - startTime;
#endif
  tp->particleInterRemote[chunk] += computed * num;
#if COSMO_DEBUG > 1 || defined CHANGA_REFACTOR_WALKCHECK
  tp->bucketcheckList[reqIDlist].insert(remoteBucketID);
  tp->combineKeys(remoteBucketID,reqIDlist);
#endif
  tp->finishBucket(reqIDlist);
  CkAssert(state->counterArrays[1][chunk] >= 0);
  if (state->counterArrays[1][chunk] == 0) {
    cacheManagerProxy[CkMyPe()].finishedChunk(chunk, tp->nodeInterRemote[chunk]+tp->particleInterRemote[chunk]);
    if (chunk == tp->numChunks-1) tp->markWalkDone();
  }

}

void PrefetchCompute::recvdParticles(ExternalGravityParticle *egp,int num,int chunk,int reqID,State *state, TreePiece *tp, Tree::NodeKey &remoteBucket){
  // when we receive missed particles, we do the same book-keeping
  // as we did when we received a missed node.
  finishNodeProcessEvent(tp, state);
}

int GravityCompute::nodeRecvdEvent(TreePiece *owner, int chunk, State *state, int reqIDlist){	
	
	state->counterArrays[0][reqIDlist]--;
	owner->finishBucket(reqIDlist);

		CkAssert(chunk >= 0); 
		state->counterArrays[1][chunk] --;
		CkAssert(state->counterArrays[1][chunk] >= 0);
		if (state->counterArrays[1][chunk] == 0) {
			streamingCache[CkMyPe()].finishedChunk(chunk, owner->nodeInterRemote[chunk]+owner->particleInterRemote[chunk]);
			if(chunk == owner->numChunks-1) owner->markWalkDone();
		}// end if finished with chunk
}

#if INTERLIST_VER > 0
int ListCompute::nodeRecvdEvent(TreePiece *owner, int chunk, State *state, int reqIDlist){	
  int start, end;
  GenericTreeNode *source = (GenericTreeNode *)computeEntity;
  owner->getBucketsBeneathBounds(source, start, end); 
  owner->updateBucketState(start, end, 1, chunk, state);

  CkAssert(chunk >= 0); 
  int remainingChunk;
  remainingChunk = state->counterArrays[1][chunk];
#ifdef CHANGA_REFACTOR_MEMCHECK
  CkPrintf("memcheck after noderecvd\n");
  CmiMemoryCheck();
#endif
  CkAssert(remainingChunk >= 0);
  if (remainingChunk == 0) {
    streamingCache[CkMyPe()].finishedChunk(chunk, owner->nodeInterRemote[chunk]+owner->particleInterRemote[chunk]);
    if(chunk == owner->numChunks-1) owner->markWalkDone();
  }// end if finished with chunk
}
#endif


/*
void PrefetchCompute::walkDone(){
  deleteComputeEntity();
}
*/

int GravityCompute::doWork(GenericTreeNode *node, TreeWalk *tw, 
                              State *state, int chunk, int reqID, bool isRoot, bool &didcomp, int awi){
  // ignores state
  
  TreePiece *tp = tw->getOwnerTP();
  if(node->getType() == Empty || node->getType() == CachedEmpty){
#ifdef CHANGA_REFACTOR_WALKCHECK
    if(node->parent->getType() != Boundary || getOptType() == Local){
      int bucketIndex = decodeReqID(reqID);
      tp->addToBucketChecklist(bucketIndex, node->getKey());
      tp->combineKeys(node->getKey(), bucketIndex);
    }
#endif
    return DUMP;
  }
  int open;

  open = openCriterion(tp, node, reqID, state);
  if(opt == NULL){       
    ckerr << "GravityCompute reqID("<<reqID<<"), isRoot("<<isRoot<<") has NULL opt" << endl;
    CkAbort("aborting");
  }
  
  int action = opt->action(open, node);
  if(action == KEEP){ // keep node
    return KEEP;
  }
  else if(action == COMPUTE){
    didcomp = true;
#ifdef BENCHMARK_TIME_COMPUTE
    double startTime = CmiWallTimer();
#endif
    int computed = nodeBucketForce(node, 
                    (GenericTreeNode *)computeEntity, 
                    tp->getParticles(), 
                    tp->decodeOffset(reqID), 
                    activeRung);
#ifdef BENCHMARK_TIME_COMPUTE
    computeTimeNode += CmiWallTimer() - startTime;
#endif
    if(getOptType() == Remote){
      tp->addToNodeInterRemote(chunk, computed);
    }
    else if(getOptType() == Local){
      tp->addToNodeInterLocal(computed);
    }
#ifdef CHANGA_REFACTOR_WALKCHECK
    int bucketIndex = decodeReqID(reqID);
    tp->addToBucketChecklist(bucketIndex, node->getKey());
    tp->combineKeys(node->getKey(), bucketIndex);
#endif

    return DUMP;
  }
  else if(action == KEEP_LOCAL_BUCKET){
    didcomp = true;
#if CHANGA_REFACTOR_DEBUG > 2
    CkAssert(node->getType() == Bucket);
    CkPrintf("[%d] GravityCompute told to KEEP_LOCAL_BUCKET, chunk=%d, remoteIndex=%d, first=%d, last=%d, reqID=%d\n", tp->getIndex(),
                         chunk, node->remoteIndex, 
                         node->firstParticle, 
                         node->lastParticle,
                         reqID);
#endif
      // since this is a local bucket, we should have the particles at hand
      GravityParticle *part = node->particlePointer;
      CkAssert(part);
      //ckerr << "(keep_local_bucket) particlePointer[0] - mass: " << part[0].mass << endl;
      int computed = 0;
#ifdef BENCHMARK_TIME_COMPUTE
      double startTime = CmiWallTimer();
#endif
      for(int i = node->firstParticle; i <= node->lastParticle; i++){
        computed += partBucketForce(
                                  &part[i-node->firstParticle],
                                  (GenericTreeNode *)computeEntity, 
                                  tp->getParticles(), 
                                  tp->decodeOffset(reqID), 
                                  activeRung);
      }
#ifdef BENCHMARK_TIME_COMPUTE
      computeTimePart += CmiWallTimer() - startTime;
#endif
      // we could have done the following, because this is a KEEP_LOCAL_BUCKET
      //
      // tp->addToParticleInterLocal(computed);
      //
      // to be sure, though, we do instead:
      if(getOptType() == Remote){
        tp->addToParticleInterRemote(chunk, computed);
      }
      else if(getOptType() == Local){
        tp->addToParticleInterLocal(computed);
      }
#ifdef CHANGA_REFACTOR_WALKCHECK
      int bucketIndex = decodeReqID(reqID);
      tp->addToBucketChecklist(bucketIndex, node->getKey());
      tp->combineKeys(node->getKey(), bucketIndex);
#endif

      return DUMP;
  }
  else if(action == KEEP_REMOTE_BUCKET){
    didcomp = true;
  // fetch particles and compute.
#if CHANGA_REFACTOR_DEBUG > 2
    CkPrintf("[%d] GravityCompute told to KEEP_REMOTE_BUCKET, chunk=%d, remoteIndex=%d, first=%d, last=%d, reqID=%d\n", tp->getIndex(),
                         chunk, node->remoteIndex, 
                         node->firstParticle, 
                         node->lastParticle,
                         reqID);
#endif
    Tree::NodeKey keyref = node->getKey();
    ExternalGravityParticle *part;
    part = 
        tp->particlesMissed(keyref, 
                                       chunk, 
                                       node->remoteIndex, 
                                       node->firstParticle, 
                                       node->lastParticle, 
                                       reqID, false, awi, computeEntity);
    if(part){
#if CHANGA_REFACTOR_DEBUG > 2
      CkPrintf("Particles found in cache\n");
#endif
      int computed = computeParticleForces(tp, node, part, reqID);
      // jetley
      //CkAssert(chunk >= 0);
      //tp->addToParticleInterRemote(chunk, computed);
      if(getOptType() == Remote){
        tp->addToParticleInterRemote(chunk, computed);
      }
      else if(getOptType() == Local){
        tp->addToParticleInterLocal(computed);
      }

    }
    else{
#if CHANGA_REFACTOR_DEBUG > 2
      CkPrintf("Particles not found in cache\n");
#endif
      if(getOptType() == Remote){
        state->counterArrays[0][decodeReqID(reqID)] += 1;
        state->counterArrays[1][chunk] += 1;
      }
    }   
    return DUMP;
  }
  else if(action == DUMP || action == NOP){
    return DUMP;
  }
}

// Source of force is a group of particles
int GravityCompute::computeParticleForces(TreePiece *ownerTP, GenericTreeNode *node, ExternalGravityParticle *part, int reqID){
  int computed = 0;
#ifdef BENCHMARK_TIME_COMPUTE
  double startTime = CmiWallTimer();
#endif
  for(int i = node->firstParticle; i <= node->lastParticle; i++){
    computed += partBucketForce(
                                  &part[i-node->firstParticle],
                                  (GenericTreeNode *)computeEntity, 
                                  ownerTP->getParticles(), 
                                  ownerTP->decodeOffset(reqID), 
                                  activeRung);
  }
#ifdef BENCHMARK_TIME_COMPUTE
  computeTimePart += CmiWallTimer() - startTime;
#endif
#ifdef CHANGA_REFACTOR_WALKCHECK
  int bucketIndex = decodeReqID(reqID);
  ownerTP->addToBucketChecklist(bucketIndex, node->getKey());
  ownerTP->combineKeys(node->getKey(), bucketIndex);
#endif
  return computed;
}

// Source of force is a node 
int GravityCompute::computeNodeForces(TreePiece *ownerTP, GenericTreeNode *node, int reqID){
  return -1;
}

int PrefetchCompute::openCriterion(TreePiece *ownerTP, 
                          GenericTreeNode *node, int reqID, State *state){
  TreePiece *tp = ownerTP;
  PrefetchRequestStruct prs(tp->prefetchReq, tp->numPrefetchReq);
  Vector3D<double> offset = ownerTP->decodeOffset(reqID);
  for(int i = 0; i < prs.numPrefetchReq; i++){
    BinaryTreeNode testNode;
    testNode.boundingBox = prs.prefetchReq[i];
    //testNode.moments.softBound = 0.0;

    if(openCriterionBucket(node, &testNode, offset, ownerTP->getIndex()))
      return 1;
  }
  return 0;  
}

int PrefetchCompute::doWork(GenericTreeNode *node, TreeWalk *tw, State *state, int chunk, int reqID, bool isRoot, bool &didcomp, int awi){
  TreePiece *tp = tw->getOwnerTP();
  // ignores state 
  if(node == NULL){
    CkAbort("PrefetchComputedoWork() given NULL node");
  }
  int open = 0;
  open = openCriterion(tp, node, reqID, state);

  int decision = opt->action(open, node);
  if (decision == DUMP || decision == KEEP){
    return decision;
  }
  else if(decision == KEEP_REMOTE_BUCKET){
#if CHANGA_REFACTOR_DEBUG > 2
    CkPrintf("[%d] PrefetchCompute told to KEEP_REMOTE_BUCKET, chunk=%d, remoteIndex=%d, first=%d, last=%d, reqID=%d\n", tp->getIndex(),
                         chunk, node->remoteIndex, 
                         node->firstParticle, 
                         node->lastParticle,
                         reqID);
#endif
    Tree::NodeKey keyref = node->getKey();
    ExternalGravityParticle *part = 
                      tp->particlesMissed(keyref, 
                                                         chunk, 
                                                         node->remoteIndex, 
                                                         node->firstParticle, 
                                                         node->lastParticle, 
                                                         reqID,
                                                         true, awi, (void *)0);
    
    if(part == NULL){
      //tp->incPrefetchWaiting();
      state->counterArrays[0][0]++;
#if CHANGA_REFACTOR_DEBUG > 2
      CkPrintf("[%d] Particles not found in cache\n", tp->getIndex());
#endif
    }
    else{
#if CHANGA_REFACTOR_DEBUG > 2
      CkPrintf("[%d] Particles found in cache\n", tp->getIndex());
#endif
    }
    return DUMP;
  }
  else if(decision == KEEP_LOCAL_BUCKET){
#if CHANGA_REFACTOR_DEBUG > 2
    CkPrintf("[%d] PrefetchCompute told to KEEP_LOCAL_BUCKET, chunk=%d, remoteIndex=%d, first=%d, last=%d, reqID=%d\n", tp->getIndex(),
                         chunk, node->remoteIndex, 
                         node->firstParticle, 
                         node->lastParticle,
                         reqID);
#endif
    CkAbort("PrefetchOpt told PrefetchCompute to KEEP_LOCAL_BUCKET. This shouldn't happen.\n");
  }
}

#if INTERLIST_VER > 0
/* List compute object
 * requires state object
 * to store lists of particles/nodes in.
 * */

/// source (computeEntity) is the current local node
/// target (encoded in reqID) is the target bucket (only nodeMissed requests during remote walks make sense)
/// node is the global node being processed
int ListCompute::doWork(GenericTreeNode *node, TreeWalk *tw, State *state, int chunk, int reqID, bool isRoot, bool &didcomp, int awi){

  DoubleWalkState *s = (DoubleWalkState *)state;
  int level = s->level;
  CheckList &chklist = s->chklists[level];
  UndecidedList &undlist = s->undlists[level];

  TreePiece *tp = tw->getOwnerTP();
  Vector3D<double> offset = tp->decodeOffset(reqID);

  if(node->getType() == Empty || node->getType() == CachedEmpty){
#ifdef CHANGA_REFACTOR_WALKCHECK_INTERLIST
    if(node->parent->getType() != Boundary || getOptType() == Local){
      ddNodeToInt(node, reqID, s);
    }
#endif
    return DUMP;
  }
  // check opening criterion
  int open;

  open = openCriterion(tp, node, reqID, state);
  //CkPrintf("[%d] open: %d\n", tp->getIndex(), open);
  
  int fakeOpen;
  // in the interlist version, there are three possible return 
  // values for the opencriterin function, whereas the Opt object
  // only knows about two (true, open and false, no need to open).
  // convert -1 (ancestor fully contained) to true (open) 
  if(open == INTERSECT || open == CONTAIN){
    fakeOpen = 1;
  }
  else{
    fakeOpen = 0;
  }

  int action = opt->action(fakeOpen, node);
  if(action == KEEP){
    if(open == CONTAIN)
    {
      // ancestor is a node and is contained 
      // only enqueue nodes if they are available
      addChildrenToCheckList(node, reqID, chunk, awi, s, chklist, tp);
      // DUMP, because 
      // we've just added the children of the glbl node to the chklist,
      // so it won't be empty when dowork returns to dft
      // if the local node needs to be opened by someone in the chklist, 
      // it will be, eventually. We should return KEEP only when we have
      // added things to the undecided list,
      // i.e. when we are emptying the chklist and not adding nodes to it
      // Technically, this should be a NOP
      return DUMP;
    }// end if contain
    else if(open == INTERSECT)
    {
      GenericTreeNode *localNode = (GenericTreeNode *)computeEntity; 
      if(localNode->getType() == Bucket){
        // if the local node is a bucket, we shouldn't descend further locally
        // instead, add the children of the glblnode to its checklist 
        addChildrenToCheckList(node, reqID, chunk, awi, s, chklist, tp);
        //Vector3D<double> vec = tp->decodeOffset(reqID);
        //CkPrintf("level %d: %d (%f, %f, %f)\n", s->level, node->getKey(), vec.x, vec.y, vec.z);
        return DUMP;
      }
      else{
        // ancestor is a node but only intersects, 
        // modify undecided list and send KEEP to tw:
        OffsetNode on;
        on.node = node;
        on.offsetID = reqID;
        undlist.push_back(on);
        return KEEP;
      }
    }
  }
  else if(action == COMPUTE){
    // only nodes can be COMPUTEd
    // particles must be KEEP_*_BUCKETed
    // so we only need to add node to clist here
    didcomp = true;
    int computed;
    
    // add to list
    /*
    Vector3D<double> v = tp->decodeOffset(reqID);
    CkPrintf("[%d] added node %ld (%1.0f,%1.0f,%1.0f) to intlist\n", tp->getIndex(),
                                                                   node->getKey(),
                                                                   v.x, v.y, v.z);
    */
    addNodeToInt(node, reqID, s);
    // all particles beneath this node have been 
    // scheduled for computation
    computed = node->lastParticle-node->firstParticle+1;
    if(getOptType() == Remote){
      tp->addToNodeInterRemote(chunk, computed);
    }
    else if(getOptType() == Local){
      tp->addToNodeInterLocal(computed);
    }
    return DUMP;
  }
  else if(action == KEEP_LOCAL_BUCKET){
    didcomp = true;
#if CHANGA_REFACTOR_DEBUG > 2
    CkAssert(node->getType() == Bucket);
    CkPrintf("[%d] ListCompute told to KEEP_LOCAL_BUCKET, chunk=%d, remoteIndex=%d, first=%d, last=%d, reqID=%d\n", tp->getIndex(),
                         chunk, node->remoteIndex, 
                         node->firstParticle, 
                         node->lastParticle,
                         reqID);
#endif
    // since this is a local bucket, we should have the particles at hand
    GravityParticle *part = node->particlePointer;
    CkAssert(part);
    int computed = node->lastParticle-node->firstParticle+1;
#if CHANGA_REFACTOR_PRINT_INTERACTIONS
    NodeKey key = node->getKey();
    addLocalParticlesToInt(part, computed, offset, s, key);
#else
    addLocalParticlesToInt(part, computed, offset, s);
#endif

    if(getOptType() == Remote){
      tp->addToParticleInterRemote(chunk, computed);
    }
    else if(getOptType() == Local){
      tp->addToParticleInterLocal(computed);
    }
    return DUMP;
  }
  else if(action == KEEP_REMOTE_BUCKET){
    didcomp = true;
  // fetch particles and compute.
#if CHANGA_REFACTOR_DEBUG > 2
    CkPrintf("[%d] ListCompute told to KEEP_REMOTE_BUCKET, chunk=%d, remoteIndex=%d, first=%d, last=%d, reqID=%d\n", tp->getIndex(),
                         chunk, node->remoteIndex, 
                         node->firstParticle, 
                         node->lastParticle,
                         reqID);
#endif
    Tree::NodeKey keyref = node->getKey();
    ExternalGravityParticle *part;
    part = 
        tp->particlesMissed(keyref, 
                                       chunk, 
                                       node->remoteIndex, 
                                       node->firstParticle, 
                                       node->lastParticle, 
                                       reqID, false, awi, computeEntity);
    if(part){
#if CHANGA_REFACTOR_DEBUG > 2
      CkPrintf("Particles found in cache\n");
#endif
      int computed = node->lastParticle-node->firstParticle+1;

#if CHANGA_REFACTOR_PRINT_INTERACTIONS
      NodeKey key = node->getKey();
      addRemoteParticlesToInt(part, computed, offset, s, key);
#else
      addRemoteParticlesToInt(part, computed, offset, s);
#endif

      if(getOptType() == Remote){
        tp->addToParticleInterRemote(chunk, computed);
      }
      else if(getOptType() == Local){
        tp->addToParticleInterLocal(computed);
      }
    }
    else{
#if CHANGA_REFACTOR_DEBUG > 2
    	CkPrintf("Particles not found in cache\n");
#endif
    	if(getOptType() == Remote){
                // particles missed
                int start, end;
    		GenericTreeNode *source = (GenericTreeNode *)computeEntity;
                tp->getBucketsBeneathBounds(source, start, end);
                tp->updateUnfinishedBucketState(start, end, 1, chunk, state);
#ifdef CHANGA_REFACTOR_MEMCHECK
  CkPrintf("memcheck after particlesmissed (%ld)\n", keyref);
  CmiMemoryCheck();
#endif
    	}
    }   
    return DUMP;
  }
  else if(action == DUMP || action == NOP){
    return DUMP;
  }
}

void ListCompute::recvdParticles(ExternalGravityParticle *part,int num,int chunk,int reqID,State *state_,TreePiece *tp, Tree::NodeKey &remoteBucket){

  Vector3D<double> offset = tp->decodeOffset(reqID);
  int targetBucket = decodeReqID(reqID);
  CkAssert(num > 0);
  GenericTreeNode *source = (GenericTreeNode *)computeEntity;
  
  GenericTreeNode* reqnode;
  
  int startBucket;
  int end;

  DoubleWalkState *state  = (DoubleWalkState *)state_;
  tp->getBucketsBeneathBounds(source, startBucket, end);

  // init state
  bool remoteLists = state->rplists.length() > 0;

  int level = source->getLevel(source->getKey());
  for(int i = 0; i <= level; i++){
    state->clists[i].length() = 0;
  }

  if(remoteLists)
    for(int i = 0; i <= level; i++){
      state->rplists[i].length() = 0;
    }

  // put particles in list at correct level
  // (key) here.
  state->level = level;
#if CHANGA_REFACTOR_PRINT_INTERACTIONS
  NodeKey key = node->getKey();
  addRemoteParticlesToInt(part, computed, offset, s, key);
#else
  addRemoteParticlesToInt(part, num, offset, state);
#endif
  
  state->lowestNode = source;

  stateReady(state, tp, chunk, startBucket, end);

  tp->updateBucketState(startBucket, end, 1, chunk, state);

  int remainingChunk;
#ifdef CHANGA_REFACTOR_MEMCHECK
  CkPrintf("memcheck after particlesrecvd\n");
  CmiMemoryCheck();
#endif
  remainingChunk = state->counterArrays[1][chunk];
  CkAssert(remainingChunk >= 0);
  if (remainingChunk == 0) {
    //CkPrintf("Finished chunk %d from particlerecvd\n", chunk);
    cacheManagerProxy[CkMyPe()].finishedChunk(chunk, tp->nodeInterRemote[chunk]+tp->particleInterRemote[chunk]);
    if (chunk == tp->numChunks-1) tp->markWalkDone();
  }
}


// CUDA: using double instead of float for vector
int ListCompute::openCriterion(TreePiece *ownerTP, 
                          GenericTreeNode *node, int reqID, State *state){
  return 
    openCriterionNode(node,(GenericTreeNode *)computeEntity, ownerTP->decodeOffset(reqID), ownerTP->getIndex());
}

#if defined CHANGA_REFACTOR_PRINT_INTERACTIONS
void ListCompute::addRemoteParticlesToInt(ExternalGravityParticle *parts, int n, Vector3D<double> &offset, DoubleWalkState *s, NodeKey key){
#else
void ListCompute::addRemoteParticlesToInt(ExternalGravityParticle *parts, int n, Vector3D<double> &offset, DoubleWalkState *s){
#endif
  RemotePartInfo rpi;
  int level = s->level;

  rpi.particles = parts;
  rpi.numParticles = n;
  rpi.offset = offset;
#ifdef CHANGA_REFACTOR_PRINT_INTERACTIONS
  rpi.key = key;
#endif

  s->rplists[level].push_back(rpi);
}

#if defined CHANGA_REFACTOR_PRINT_INTERACTIONS
void ListCompute::addLocalParticlesToInt(GravityParticle *parts, int n, Vector3D<double> &offset, DoubleWalkState *s, NodeKey key){
#else
void ListCompute::addLocalParticlesToInt(GravityParticle *parts, int n, Vector3D<double> &offset, DoubleWalkState *s){
#endif
  LocalPartInfo lpi;
  int level = s->level;

  lpi.particles = parts;
  lpi.numParticles = n;
  lpi.offset = offset;
#ifdef CHANGA_REFACTOR_PRINT_INTERACTIONS
  lpi.key = key;
#endif

  s->lplists[level].push_back(lpi);

}

void ListCompute::addNodeToInt(GenericTreeNode *node, int offsetID, DoubleWalkState *s){
  OffsetNode nd;
  int level = s->level;
  nd.node = node;
  nd.offsetID = offsetID;

  s->clists[level].push_back(nd);
}

void ListCompute::stateReady(State *state_, TreePiece *tp, int chunk, int start, int end){
  int thisIndex = tp->getIndex();
  GravityParticle *particles = tp->getParticles();

  DoubleWalkState *state = (DoubleWalkState *)state_;
  GenericTreeNode *lowestNode = state->lowestNode;
  int maxlevel = lowestNode->getLevel(lowestNode->getKey());
  
  bool hasRemoteLists = state->rplists.length() > 0 ? true : false;
  bool hasLocalLists = state->lplists.length() > 0 ? true : false;

  for(int b = start; b < end; b++){
    if(tp->bucketList[b]->rungs >= activeRung){

      for(int level = 0; level <= maxlevel; level++){

        CkVec<OffsetNode> &clist = state->clists[level];
        for(int i = 0; i < clist.length(); i++){

#ifdef CHANGA_REFACTOR_WALKCHECK_INTERLIST
          GenericTreeNode *node = clist[i].node;
          tp->addToBucketChecklist(b, node->getKey());
          tp->combineKeys(node->getKey(), b);
          
          // don't try to compute with empty nodes or 
          if(node->getType() == Empty 
            || node->getType() == CachedEmpty){
            continue;
          }
#endif
          int computed =  nodeBucketForce(clist[i].node, 
                          tp->getBucket(b), 
                          particles, 
                          tp->decodeOffset(clist[i].offsetID), activeRung);
          if(getOptType() == Remote){
            tp->addToNodeInterRemote(chunk, computed);
          }
          else if(getOptType() == Local){
            tp->addToNodeInterLocal(computed);
          }
          Vector3D<double> vec = tp->decodeOffset(clist[i].offsetID);
#ifdef CHANGA_REFACTOR_PRINT_INTERACTIONS
          CkPrintf("[%d]: bucket %d with %ld (%1.0f,%1.0f,%1.0f)\n", thisIndex,  b, clist[i].node->getKey(), vec.x, vec.y, vec.z);
#endif

        }

        // remote particles
        if(hasRemoteLists){
          CkVec<RemotePartInfo> &rpilist = state->rplists[level];
          // for each bunch of particles in list
          for(int i = 0; i < rpilist.length(); i++){
            RemotePartInfo &rpi = rpilist[i];
#ifdef CHANGA_REFACTOR_WALKCHECK_INTERLIST
            NodeKey key = rpi.key;
            tp->addToBucketChecklist(b, key);
            tp->combineKeys(key, b);
#endif
            // for each particle in a bunch
            for(int j = 0; j < rpi.numParticles; j++){
	      int computed = partBucketForce(&rpi.particles[j], 
	                      tp->getBucket(b),
	                      particles,
	                      rpi.offset, activeRung);
              if(getOptType() == Remote){// don't really have to perform this check
	        tp->addToParticleInterRemote(chunk, computed);
              }
            }
            Vector3D<double> &vec = rpi.offset;
#ifdef CHANGA_REFACTOR_PRINT_INTERACTIONS
            CkPrintf("[%d]: bucket %d with %ld (%1.0f,%1.0f,%1.0f)\n", thisIndex, b, rpi.key, vec.x, vec.y, vec.z);
#endif
          }
        }

        // local particles
        if(hasLocalLists){
          CkVec<LocalPartInfo> &lpilist = state->lplists[level];
          for(int i = 0; i < lpilist.length(); i++){
            LocalPartInfo &lpi = lpilist[i];
#ifdef CHANGA_REFACTOR_WALKCHECK_INTERLIST
            NodeKey key = lpi.key;
            tp->addToBucketChecklist(b, key);
            tp->combineKeys(key, b);
#endif
            // for each particle in a bunch
            for(int j = 0; j < lpi.numParticles; j++){
	      int computed = partBucketForce(&lpi.particles[j], 
	                      tp->getBucket(b),
	                      particles,
	                      lpi.offset, activeRung);
	      tp->addToParticleInterLocal(computed);
            }
            Vector3D<double> &vec = lpi.offset;
#ifdef CHANGA_REFACTOR_PRINT_INTERACTIONS
            CkPrintf("[%d]: bucket %d with %ld (%1.0f,%1.0f,%1.0f)\n", thisIndex, b, lpi.key, vec.x, vec.y, vec.z);
#endif
          }
        }

      }// level
    }// active
    //if(b == TEST_BUCKET)
    //CkPrintf("-----------------\n");
  }// bucket

/*
  for(int l = 0; l < maxlevel; l++){// for each level
    // nodes first
    CkVec<OffsetNode> &clist = state->clists[l];

    for(int j = 0; j < clist.length(); j++){
      for(int b = start; b < end; b++){ 
        if(tp->bucketList[b]->rungs >= activeRung){
          int computed =  nodeBucketForce(clist[j].node, 
                          tp->getBucket(b), 
                          particles, 
                          tp->decodeOffset(clist[j].offsetID), activeRung);
          if(getOptType() == Remote){
            tp->addToNodeInterRemote(chunk, computed);
          }
          else if(getOptType() == Local){
            tp->addToNodeInterLocal(computed);
          }
          Vector3D<double> vec = tp->decodeOffset(clist[j].offsetID);
          CkPrintf("%d: %d (%f, %f, %f)\n", b, clist[j].node->getKey(), vec.x, vec.y, vec.z);
        }// end if active
      }// end fo each bucket
    }// end for each item in clist
    int computed = 0;
    // remote particles 
    if(hasRemoteLists){
	    CkVec<RemotePartInfo> &rplist = state->rplists[l];
	    for(int j = 0; j < rplist.length(); j++){// for each entry in this list
	      RemotePartInfo &rpi = rplist[j];// entry
	      for(int part = 0; part < rpi.numParticles; part++){ // for each particle in entry
	        for(int b = start; b < end; b++){ // for each bucket
                  
                  if(tp->bucketList[b]->rungs >= activeRung){
	            computed = partBucketForce(&rpi.particles[part], 
	                            tp->getBucket(b),
	                            particles,
	                            rpi.offset, activeRung);
	            tp->addToParticleInterRemote(chunk, computed);
                    Vector3D<double> &vec = rpi.offset;
                    CkPrintf("%d: %d (%f, %f, %f)\n", b, rpi.nd->getKey(), vec.x, vec.y, vec.z);
                  }// end if bucket active
	        }// end for each bucket
	      }// end for each particle in listitem
	    }// end for each listitem in list
    }// end has remote particle lists
    // local particles 
    if(hasLocalLists){
	    CkVec<LocalPartInfo> &lplist = state->lplists[l];
	    for(int j = 0; j < lplist.length(); j++){// for each entry in this list
	      LocalPartInfo &lpi = lplist[j];// entry
	      for(int part = 0; part < lpi.numParticles; part++){ // for each particle in entry
	        for(int b = start; b < end; b++){ // for each bucket
                  
                  if(tp->bucketList[b]->rungs >= activeRung){
	            computed = partBucketForce(&lpi.particles[part], 
	                            tp->getBucket(b),
	                            particles,
	                            lpi.offset, activeRung);
	            tp->addToParticleInterLocal(computed);
                    Vector3D<double> &vec = lpi.offset;
                    CkPrintf("%d: %d (%f, %f, %f)\n", b, lpi.nd->getKey(), vec.x, vec.y, vec.z);
                  }// end if bucket active
	        }// end for each bucket 
	      }// end for each particle in listitem 
	    }// end for each listitem in list
    }// end has local particle lists
  }
  */
}

void ListCompute::addChildrenToCheckList(GenericTreeNode *node, int reqID, int chunk, int awi, State *s, CheckList &chklist, TreePiece *tp){

  Vector3D<double> vec = tp->decodeOffset(reqID);

  for(int i = 0; i < node->numChildren(); i++)
  {
    GenericTreeNode *child = node->getChildren(i);
    if(!child)
    {
      Tree::NodeKey childKey = node->getChildKey(i);
      // child not in my tree, look in cache
      child = tp->nodeMissed(reqID, 
                                node->remoteIndex, 
                                childKey,
                                chunk,
                                false, 
                                awi, computeEntity);
      if(!child)
      {
        nodeMissedEvent(reqID, chunk, s, tp);
#ifdef CHANGA_REFACTOR_INTERLIST_PRINT_LIST_STATE
        CkPrintf("[%d] missed child %ld (%1.0f,%1.0f,%1.0f) of %ld\n", tp->getIndex(), 
                                                                childKey,
                                                                vec.x, vec.y, vec.z,
                                                                node->getKey());
#endif
        continue;
      }
    } 
#ifdef CHANGA_REFACTOR_INTERLIST_PRINT_LIST_STATE
    CkPrintf("[%d] enq avail child %ld (%1.0f,%1.0f,%1.0f) of %ld\n", tp->getIndex(), 
                                                                child->getKey(), 
                                                                vec.x, vec.y, vec.z,
                                                                node->getKey());
#endif
    // child available, enqueue
    OffsetNode on;
    on.node = child;
    on.offsetID = reqID;
    chklist.enq(on);
  }
}

// to be used when debugging
#ifdef CHANGA_REFACTOR_INTERLIST_PRINT_LIST_STATE
void printClist(DoubleWalkState *state, int level, TreePiece *tp){
  
  {
    CkVec<OffsetNode> &list = state->clists[level];
    Vector3D<double> vec;
    for(int i = 0; i < list.length(); i++){
      vec = tp->decodeOffset(list[i].offsetID);
      CkPrintf("%ld(%1.0f,%1.0f,%1.0f)\n", list[i].node->getKey(), vec.x, vec.y, vec.z);
    }
  }

  bool hasRemoteLists = state->rplists.length() > 0 ? true : false;
  bool hasLocalLists = state->lplists.length() > 0 ? true : false;

  if(hasLocalLists){
    CkPrintf("----------------\n");
    CkVec<LocalPartInfo> &list = state->lplists[level];
    Vector3D<double> vec;
    for(int i = 0; i < list.length(); i++){
      vec = list[i].offset;
      CkPrintf("%ld(%1.0f,%1.0f,%1.0f)\n", list[i].key, vec.x, vec.y, vec.z);
    }
  }

  if(hasRemoteLists){
    CkPrintf("----------------\n");
    CkVec<RemotePartInfo> &list = state->rplists[level];
    Vector3D<double> vec;
    for(int i = 0; i < list.length()-1; i++){
      vec = list[i].offset;
      CkPrintf("%ld(%1.0f,%1.0f,%1.0f)\n", list[i].key, vec.x, vec.y, vec.z);
    }
  }
}

void printUndlist(DoubleWalkState *state, int level, TreePiece *tp){
    CkPrintf("----------------\n");
    CkVec<OffsetNode> &list = state->undlists[level];
    Vector3D<double> vec;
    for(int i = 0; i < list.length()-1; i++){
      vec = tp->decodeOffset(list[i].offsetID);
      CkPrintf("%ld(%1.0f,%1.0f,%1.0f)\n", list[i].node->getKey(), vec.x, vec.y, vec.z);
    }
}
#endif

#endif // INTERLIST_VER > 0
