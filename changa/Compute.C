#include "ParallelGravity.h"
#include "GenericTreeNode.h"
//#include "codes.h"

#include "Opt.h"
#include "Compute.h"
#include "TreeWalk.h"
//#include "State.h"
#include "Space.h"
#include "gravity.h"

#ifdef CELL
#include "spert_ppu.h"
#include "cell_typedef.h"

extern void cellSPE_callback(void *data);
extern void cellSPE_single(void *data);
extern int workRequestOut;
#endif

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
#ifdef CUDA
        // wait for prefetched particles as well
        // this way, all nodes/parts not missed will be handled by RNR
        // and all those missed, by RR
        // if we didn't wait for particles
        // it could transpire that we don't miss on these particles during RNR, but they aren't present on the gpu
        // and so we'd have to have a separate array of missed particles for the RNR (in much the same way that the RR has
        // separate arrays for missed nodes and particles) 
  finishNodeProcessEvent(tp, state);
#endif
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
    //CkPrintf("[%d] Finished chunk %d from nodeRecvdEvent\n",owner->getIndex(),chunk);
#ifdef CUDA
    // no more nodes/particles are going to be delivered by the cache
    // flush the interactions remaining in the state
    DoubleWalkState *ds = (DoubleWalkState *)state;
    CkAssert(ds->resume);
    // at this point,
    sendNodeInteractionsToGpu(ds, owner, true);
    resetCudaNodeState(ds);
#endif
    streamingCache[CkMyPe()].finishedChunk(chunk, owner->nodeInterRemote[chunk]+owner->particleInterRemote[chunk]);
    if(chunk == owner->numChunks-1){
      owner->markWalkDone();
    }
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
#ifdef CUDA
      // waiting for particles for reasons discussed 
      // in comments within recvdParticles
      state->counterArrays[0][0]++;
#endif
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
#if defined CHANGA_REFACTOR_PRINT_INTERACTIONS || defined CUDA
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

#if defined CHANGA_REFACTOR_PRINT_INTERACTIONS || defined CUDA
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
#ifdef CUDA
  state->numInteractions = 0;
#endif
#if defined CHANGA_REFACTOR_PRINT_INTERACTIONS || defined CUDA
  NodeKey key = remoteBucket;
  addRemoteParticlesToInt(part, num, offset, state, key);
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
    //CkPrintf("[%d] Finished chunk %d from recvdParticles\n",tp->getIndex(),chunk);
#ifdef CUDA
    // XXX - since no more chunks or pending
    // interactions remain, we must flush the particle lists
    CkAssert(state->resume);
    // at this point,
    sendPartInteractionsToGpu(state, tp, true);
    resetCudaPartState(state);
#endif
    cacheManagerProxy[CkMyPe()].finishedChunk(chunk, tp->nodeInterRemote[chunk]+tp->particleInterRemote[chunk]);
    if (chunk == tp->numChunks-1){
      tp->markWalkDone();
    }
  }
}


// CUDA: using double instead of float for vector
int ListCompute::openCriterion(TreePiece *ownerTP,
                          GenericTreeNode *node, int reqID, State *state){
  return
    openCriterionNode(node,(GenericTreeNode *)computeEntity, ownerTP->decodeOffset(reqID), ownerTP->getIndex());
}

#if defined CHANGA_REFACTOR_PRINT_INTERACTIONS || defined CUDA
void ListCompute::addRemoteParticlesToInt(ExternalGravityParticle *parts, int n, Vector3D<double> &offset, DoubleWalkState *s, NodeKey key){
#else
void ListCompute::addRemoteParticlesToInt(ExternalGravityParticle *parts, int n, Vector3D<double> &offset, DoubleWalkState *s){
#endif
  RemotePartInfo rpi;
  int level = s->level;

  rpi.particles = parts;
  rpi.numParticles = n;
  rpi.offset = offset;
#if defined CHANGA_REFACTOR_PRINT_INTERACTIONS || defined CUDA
  rpi.key = key;
#endif

#ifdef CUDA
  s->numInteractions++;
#endif

  s->rplists[level].push_back(rpi);
}

#if defined CHANGA_REFACTOR_PRINT_INTERACTIONS || defined CUDA
void ListCompute::addLocalParticlesToInt(GravityParticle *parts, int n, Vector3D<double> &offset, DoubleWalkState *s, NodeKey key){
#else
void ListCompute::addLocalParticlesToInt(GravityParticle *parts, int n, Vector3D<double> &offset, DoubleWalkState *s){
#endif
  LocalPartInfo lpi;
  int level = s->level;

  lpi.particles = parts;
  lpi.numParticles = n;
  lpi.offset = offset;
#if defined CHANGA_REFACTOR_PRINT_INTERACTIONS || defined CUDA
  lpi.key = key;
#endif

#ifdef CUDA
  s->numInteractions++;
#endif

  s->lplists[level].push_back(lpi);

}

void ListCompute::addNodeToInt(GenericTreeNode *node, int offsetID, DoubleWalkState *s){
  OffsetNode nd;
  int level = s->level;
  nd.node = node;
  nd.offsetID = offsetID;

  #ifdef CUDA
  s->numInteractions++;
#endif

  s->clists[level].push_back(nd);
}

#ifdef CUDA
extern int encodeOffset(int reqID, int x, int y, int z);
#include "DataManager.h"
#endif

#ifdef CELL
#define CELLBUFFERSIZE 16*1024
#endif

void ListCompute::stateReady(State *state_, TreePiece *tp, int chunk, int start, int end){
  int thisIndex = tp->getIndex();
  GravityParticle *particles = tp->getParticles();

  DoubleWalkState *state = (DoubleWalkState *)state_;

#ifdef CUDA
  int cudaNodeThreshold = state->nodeThreshold;
  int cudaPartThreshold = state->partThreshold;
  bool resume = state->resume;
  int numActiveBuckets = tp->numActiveBuckets;

  // for local particles
  std::map<NodeKey, int> &lpref = tp->dm->getLocalPartsOnGpuTable();
  // for cached particles
  std::map<NodeKey, int> &cpref = tp->dm->getCachedPartsOnGpuTable();
#endif

  GenericTreeNode *lowestNode = state->lowestNode;
  int maxlevel = lowestNode->getLevel(lowestNode->getKey());

  bool hasRemoteLists = state->rplists.length() > 0 ? true : false;
  bool hasLocalLists = state->lplists.length() > 0 ? true : false;

#ifdef CUDA
  CkPrintf("[%d]: stateReady(%d-%d) %s, %d\n",  tp->getIndex(), start, end, getOptType() == Remote ? "Remote" : "Local", state->resume);
  int numInteractions = state->numInteractions;
#endif

#ifndef CUDA
  for(int b = start; b < end; b++){
    if(tp->bucketList[b]->rungs >= activeRung){
#if defined CELL
      //enableTrace();
      // Combine all the computation in a request to be sent to the cell SPE
      int activePart=0;
      for (int k=tp->bucketList[b]->firstParticle; k<=tp->bucketList[b]->lastParticle; ++k) {
        if (tp->myParticles[k].rung >= activeRung) activePart++;
      }

      int indexActivePart=0;
      int activePartDataSize = ROUNDUP_128((activePart+3)*sizeof(CellGravityParticle));
      CellGravityParticle *activePartData = (CellGravityParticle*)malloc_aligned(activePartDataSize,128);
      GravityParticle **partList = new GravityParticle*[activePart];
      CellGroupRequest *cgr = new CellGroupRequest(tp, b, partList, state);
      WRGroupHandle wrh = createWRGroup(cgr, cellSPE_callback);
      for (int k=tp->bucketList[b]->firstParticle; k<=tp->bucketList[b]->lastParticle; ++k) {
        if (tp->myParticles[k].rung >= activeRung) {
          activePartData[indexActivePart] = tp->myParticles[k];
          //CkPrintf("[%d] particle %d: %d on bucket %d, acc: %f %f %f\n",thisIndex,k,tp->myParticles[k].iOrder,b,activePartData[indexActivePart].treeAcceleration.x,activePartData[indexActivePart].treeAcceleration.y,activePartData[indexActivePart].treeAcceleration.z);
          //CkPrintf("[%d] partList %d: particle %d (%p)\n",thisIndex,indexActivePart,k,&tp->myParticles[k]);
          partList[indexActivePart++] = &tp->myParticles[k];
        }
      }
      /*int numPart=0, numNodes=0;
        for(int k=0;k<=curLevelLocal;k++) {
        numNodes += cellListLocal[k].length();
        for(int kk=0; kk<particleListLocal[k].length(); kk++) {
        numPart += particleListLocal[k][kk].numParticles;
        }
        }*/
      int particlesPerRequest = (CELLBUFFERSIZE - 2*sizeof(int)) / sizeof(CellExternalGravityParticle);
      int nodesPerRequest = (CELLBUFFERSIZE - 2*sizeof(int)) / sizeof(CellMultipoleMoments);
      CellContainer *particleContainer = (CellContainer*)malloc_aligned(CELLBUFFERSIZE,128);
      CellExternalGravityParticle *particleData = (CellExternalGravityParticle*)&particleContainer->data;
      CellContainer *nodesContainer = (CellContainer*)malloc_aligned(CELLBUFFERSIZE,128);
      CellMultipoleMoments *nodesData = (CellMultipoleMoments*)&nodesContainer->data;
      int indexNodes=0, indexPart=0;
#endif

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
          int computed = 0;
#ifdef CELL_NODE
          nodesData[indexNodes] = clist[i].node->moments;
          Vector3D<double> tmpoffsetID = tp->decodeOffset(clist[i].offsetID);
          nodesData[indexNodes].cm += tmpoffsetID;
          indexNodes++;
          if (indexNodes == nodesPerRequest) {
            // send off request
            void *activeData = malloc_aligned(activePartDataSize,128);
            memcpy(activeData, activePartData, activePart*sizeof(CellGravityParticle));
            CellRequest *userData = new CellRequest((CellGravityParticle*)activeData, activePart, nodesContainer, partList, tp);
            nodesContainer->numInt = activePart;
            nodesContainer->numExt = indexNodes;
            //CkPrintf("[%d] sending request 1 %p+%d, %p+%d\n",thisIndex,activeData, activePartDataSize, nodesContainer, CELLBUFFERSIZE);
            sendWorkRequest (1, activeData, activePartDataSize, nodesContainer, CELLBUFFERSIZE, NULL, 0,
                (void*)userData, WORK_REQUEST_FLAGS_BOTH_CALLBACKS, cellSPE_single, wrh);
            workRequestOut ++;
            computed += activePart * indexNodes;
            nodesContainer = (CellContainer*)malloc_aligned(CELLBUFFERSIZE,128);
            nodesData = (CellMultipoleMoments*)&nodesContainer->data;
            indexNodes = 0;
          }
#else
          computed =  nodeBucketForce(clist[i].node,
              tp->getBucket(b),
              particles,
              tp->decodeOffset(clist[i].offsetID), activeRung);
#endif
          if(getOptType() == Remote){
            tp->addToNodeInterRemote(chunk, computed);
          }
          else if(getOptType() == Local){
            tp->addToNodeInterLocal(computed);
          }
#ifdef CHANGA_REFACTOR_PRINT_INTERACTIONS
          Vector3D<double> vec = tp->decodeOffset(clist[i].offsetID);
          CkPrintf("[%d]: bucket %d with %ld (%1.0f,%1.0f,%1.0f)\n", thisIndex,  b, clist[i].node->getKey(), vec.x, vec.y, vec.z);
#endif

        }
#ifdef CELL
        if (indexNodes > 0 && level==maxlevel) {
#ifdef CELL_NODE
          void *activeData = malloc_aligned(activePartDataSize,128);
          memcpy(activeData, activePartData, activePart*sizeof(CellGravityParticle));
          CellRequest *userData = new CellRequest((CellGravityParticle*)activeData, activePart, nodesContainer, partList, tp);
          nodesContainer->numInt = activePart;
          nodesContainer->numExt = indexNodes;
          //CkPrintf("[%d] sending request 2 %p+%d, %p+%d\n",thisIndex,activeData, activePartDataSize, nodesContainer, ROUNDUP_128(indexNodes*sizeof(CellMultipoleMoments)));
          sendWorkRequest (1, activeData, activePartDataSize, nodesContainer, ROUNDUP_128(indexNodes*sizeof(CellMultipoleMoments)+2*sizeof(int)),
              NULL, 0, (void*)userData, WORK_REQUEST_FLAGS_BOTH_CALLBACKS, cellSPE_single, wrh);
          workRequestOut ++;
          if(getOptType() == Remote){
            tp->addToNodeInterRemote(chunk, activePart * indexNodes);
          }
          else if(getOptType() == Local){
            tp->addToNodeInterLocal(activePart * indexNodes);
          }
#endif
        } else if (level==maxlevel) {
          free_aligned(nodesContainer);
        }
#endif

        // remote particles
        if(hasRemoteLists){
          CkVec<RemotePartInfo> &rpilist = state->rplists[level];
          // for each bunch of particles in list
          for(int i = 0; i < rpilist.length(); i++){
            RemotePartInfo &rpi = rpilist[i];

#if defined CHANGA_REFACTOR_WALKCHECK_INTERLIST 
            NodeKey key = rpi.key;
            tp->addToBucketChecklist(b, key);
            tp->combineKeys(key, b);
#endif

            // for each particle in a bunch
            for(int j = 0; j < rpi.numParticles; j++){
              int computed = 0;
#ifdef CELL_PART
              particleData[indexPart] = rpi.particles[j];
              particleData[indexPart].position += rpi.offset;
              indexPart++;
              if (indexPart == particlesPerRequest) {
                // send off request
                void *activeData = malloc_aligned(activePartDataSize,128);
                memcpy(activeData, activePartData, activePart*sizeof(CellGravityParticle));
                CellRequest *userData = new CellRequest((CellGravityParticle*)activeData, activePart, particleContainer, partList, tp);
                particleContainer->numInt = activePart;
                particleContainer->numExt = indexPart;
                //CkPrintf("[%d] sending request 3r %p+%d, %p+%d\n",thisIndex,activeData, activePartDataSize, particleContainer, CELLBUFFERSIZE);
                sendWorkRequest (2, activeData, activePartDataSize, particleContainer, CELLBUFFERSIZE, NULL, 0,
                    (void*)userData, WORK_REQUEST_FLAGS_BOTH_CALLBACKS, cellSPE_single, wrh);
                workRequestOut ++;
                tp->addToParticleInterRemote(chunk, activePart * indexPart);
                particleContainer = (CellContainer*)malloc_aligned(CELLBUFFERSIZE,128);
                particleData = (CellExternalGravityParticle*)&particleContainer->data;
                indexPart = 0;
              }
#else
              computed = partBucketForce(&rpi.particles[j],
                  tp->getBucket(b),
                  particles,
                  rpi.offset, activeRung);
#endif // CELL_PART
              if(getOptType() == Remote){// don't really have to perform this check
                tp->addToParticleInterRemote(chunk, computed);
              }
            }

#ifdef CHANGA_REFACTOR_PRINT_INTERACTIONS
            Vector3D<double> &vec = rpi.offset;
            CkPrintf("[%d]: bucket %d with %ld (%1.0f,%1.0f,%1.0f)\n", thisIndex, b, rpi.key, vec.x, vec.y, vec.z);
#endif
          }
        }

        // local particles
        if(hasLocalLists){
          CkVec<LocalPartInfo> &lpilist = state->lplists[level];
          for(int i = 0; i < lpilist.length(); i++){
            LocalPartInfo &lpi = lpilist[i];
#if defined CHANGA_REFACTOR_WALKCHECK_INTERLIST 
            NodeKey key = lpi.key;
            tp->addToBucketChecklist(b, key);
            tp->combineKeys(key, b);
#endif

            // for each particle in a bunch
            for(int j = 0; j < lpi.numParticles; j++){
#ifdef CELL_PART
              particleData[indexPart] = lpi.particles[j];
              particleData[indexPart].position += lpi.offset;
              indexPart++;
              if (indexPart == particlesPerRequest) {
                // send off request
                void *activeData = malloc_aligned(activePartDataSize,128);
                memcpy(activeData, activePartData, activePart*sizeof(CellGravityParticle));
                CellRequest *userData = new CellRequest((CellGravityParticle*)activeData, activePart, particleContainer, partList, tp);
                particleContainer->numInt = activePart;
                particleContainer->numExt = indexPart;
                //CkPrintf("[%d] sending request 3l %p+%d, %p+%d\n",thisIndex,activeData, activePartDataSize, particleContainer, CELLBUFFERSIZE);
                sendWorkRequest (2, activeData, activePartDataSize, particleContainer, CELLBUFFERSIZE, NULL, 0,
                    (void*)userData, WORK_REQUEST_FLAGS_BOTH_CALLBACKS, cellSPE_single, wrh);
                workRequestOut ++;
                tp->addToParticleInterLocal(activePart * indexPart);
                particleContainer = (CellContainer*)malloc_aligned(CELLBUFFERSIZE,128);
                particleData = (CellExternalGravityParticle*)&particleContainer->data;
                indexPart = 0;
              }
#else
              int computed = partBucketForce(&lpi.particles[j],
                  tp->getBucket(b),
                  particles,
                  lpi.offset, activeRung);
              tp->addToParticleInterLocal(computed);
#endif// CELL_PART
            }

#ifdef CHANGA_REFACTOR_PRINT_INTERACTIONS
            Vector3D<double> &vec = lpi.offset;
            CkPrintf("[%d]: bucket %d with %ld (%1.0f,%1.0f,%1.0f)\n", thisIndex, b, lpi.key, vec.x, vec.y, vec.z);
#endif
          }
        }
#ifdef CELL
        if (indexPart > 0 && level==maxlevel) {
#ifdef CELL_PART
          //void *activeData = malloc_aligned(activePartDataSize,128);
          //memcpy(activeData, activePartData, activePart*sizeof(GravityParticle));
          CellRequest *userData = new CellRequest(activePartData, activePart, particleContainer, partList, tp);
          particleContainer->numInt = activePart;
          particleContainer->numExt = indexPart;
          //CkPrintf("[%d] sending request 4 %p+%d, %p+%d (%d int %d ext)\n",thisIndex,activePartData, activePartDataSize, particleContainer, ROUNDUP_128(indexPart*sizeof(CellExternalGravityParticle)),activePart,indexPart);
          sendWorkRequest (2, activePartData, activePartDataSize, particleContainer, ROUNDUP_128(indexPart*sizeof(CellExternalGravityParticle)+2*sizeof(int)),
              NULL, 0, (void*)userData, WORK_REQUEST_FLAGS_BOTH_CALLBACKS, cellSPE_single, wrh);
          workRequestOut ++;
          if(getOptType() == Remote){
            tp->addToNodeInterRemote(chunk, activePart * indexPart);
          }
          else if(getOptType() == Local){
            tp->addToNodeInterLocal(activePart * indexPart);
          }
#endif
        } else if (level==maxlevel) {
          free_aligned(activePartData);
          free_aligned(particleContainer);
        }
#endif

      }// level

#ifdef CELL
      // Now all the requests have been made
      completeWRGroup(wrh);
      OffloadAPIProgress();
#endif
    }// active
  }// bucket

#else // else part of ifndef CUDA
  // *********************
  for(int b = start; b < end; b++){
    if(tp->bucketList[b]->rungs >= activeRung){
      
      int activePart=0;
      for (int k=tp->bucketList[b]->firstParticle; k<=tp->bucketList[b]->lastParticle; ++k) {
        if (tp->myParticles[k].rung >= activeRung) activePart++;
      }

      if(!resume){
        // another bucket has been finished in the local/rnr mode
        numFinLocalRNRBuckets++;
      }

      int bucketSize, bucketStartIndex;
      // get bucket parameters: size, particle start
      getBucketParameters(tp, b, bucketStartIndex, bucketSize, lpref);
      for(int level = 0; level <= maxlevel; level++){

        CkVec<OffsetNode> &clist = state->clists[level];
        for(int i = 0; i < clist.length(); i++){
          numInteractions--;

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
          int computed = 0;
          // check whether it is already on the gpu
          GenericTreeNode *node = clist[i].node;
          int index = node->nodeArrayIndex;
          if(index < 0){
            // index is -1, node was missed after prefetch during some walk,
            // and this must be a resumed walk
            CkAssert(state->nodes);
            // redundant check
            CkAssert(state->resume);
            // put node in nodes
            index = state->nodes->push_back_v(CudaMultipoleMoments(node->moments));
          }
          tp->addToNodeInterRemote(chunk, activePart);

          // put index in interaction list
          {
            Vector3D<double> vec = tp->decodeOffset(clist[i].offsetID);
            CkPrintf("[%d]: bucket %d with node %ld (%d), (%1.0f,%1.0f,%1.0f)\n", thisIndex,  b, clist[i].node->getKey(), clist[i].node->nodeArrayIndex, vec.x, vec.y, vec.z);
          }
          state->cellLists.push_back(ILCell(index, clist[i].offsetID));
          
          // insert a new bucket in the list of targets only when the following condition is met:
          if(state->bucketsNodes.length() == 0 || state->bucketsNodes[state->bucketsNodes.length()-1] != b){
            insertBucketNode(state, b, bucketStartIndex, bucketSize);
          }
          if(state->cellLists.length() >= cudaNodeThreshold){
            // enough node interactions to offload
            sendNodeInteractionsToGpu(state, tp, numInteractions == 0);
            resetCudaNodeState(state);
          }
          if(getOptType() == Remote){
            tp->addToNodeInterRemote(chunk, computed);
          }
          else if(getOptType() == Local){
            tp->addToNodeInterLocal(computed);
          }
#ifdef CHANGA_REFACTOR_PRINT_INTERACTIONS
          Vector3D<double> vec = tp->decodeOffset(clist[i].offsetID);
          CkPrintf("[%d]: bucket %d with %ld (%1.0f,%1.0f,%1.0f)\n", thisIndex,  b, clist[i].node->getKey(), vec.x, vec.y, vec.z);
#endif

        }

        // remote particles
        if(hasRemoteLists){
          CkVec<RemotePartInfo> &rpilist = state->rplists[level];
          // for each bunch of particles in list
          for(int i = 0; i < rpilist.length(); i++){
            numInteractions--;
            RemotePartInfo &rpi = rpilist[i];
            NodeKey key = rpi.key;
#ifdef CHANGA_REFACTOR_WALKCHECK_INTERLIST
            tp->addToBucketChecklist(b, key);
            tp->combineKeys(key, b);
#endif

            int gpuIndex = -1;
            key <<= 1;
            std::map<NodeKey, int>::iterator q = cpref.find(key);
            if(q != cpref.end()){
              //cachedPartsOnGpu = true;
              gpuIndex = q->second;
            }

            if(gpuIndex < 0){
              // couldn't find it on the gpu (must be from a missed bucket of particles)
              CkAssert(state->particles && state->resume);
              gpuIndex = state->particles->length();
              for(int j = 0; j < rpi.numParticles; j++){
                state->particles->push_back(CompactPartData(rpi.particles[j]));
              }
            }

            tp->addToParticleInterRemote(chunk, activePart*rpi.numParticles);

            // put bucket in interaction list
            Vector3D<double> off = rpi.offset;
            {
              Vector3D<double> &vec = rpi.offset;
              CkPrintf("[%d]: bucket %d with part %ld (%d), (%1.0f,%1.0f,%1.0f)\n", thisIndex, b, rpi.key, gpuIndex, vec.x, vec.y, vec.z);
            }
            state->partLists.push_back(ILPart(gpuIndex, encodeOffset(0, off.x, off.y, off.z), rpi.numParticles));
            // insert a new bucket in the list of targets only when the following condition is met:
            if(state->bucketsParts.length() == 0 || state->bucketsParts[state->bucketsParts.length()-1] != b){
              insertBucketPart(state, b, bucketStartIndex, bucketSize);
            }
            if(state->partLists.length() >= cudaPartThreshold){
              // enough nodes to offload
              sendPartInteractionsToGpu(state, tp, numInteractions == 0);
              resetCudaPartState(state);
            }

#ifdef CHANGA_REFACTOR_PRINT_INTERACTIONS
            Vector3D<double> &vec = rpi.offset;
            CkPrintf("[%d]: bucket %d with %ld (%1.0f,%1.0f,%1.0f)\n", thisIndex, b, rpi.key, vec.x, vec.y, vec.z);
#endif
          }
        }

        // local particles
        if(hasLocalLists){
          CkVec<LocalPartInfo> &lpilist = state->lplists[level];
          for(int i = 0; i < lpilist.length(); i++){
            numInteractions--;
            LocalPartInfo &lpi = lpilist[i];
            NodeKey key = lpi.key;
#if defined CHANGA_REFACTOR_WALKCHECK_INTERLIST
            tp->addToBucketChecklist(b, key);
            tp->combineKeys(key, b);
#endif

            int gpuIndex = -1;
            key <<= 1;
            std::map<NodeKey, int>::iterator p = lpref.find(key);
            if(p != lpref.end()){
              gpuIndex = p->second;
            }
            if(gpuIndex < 0){
              // couldn't find it on the gpu (must be from a missed bucket of particles)
              CkAssert(state->particles && state->resume);
              gpuIndex = state->particles->length();
              for(int j = 0; j < lpi.numParticles; j++){
                state->particles->push_back(CompactPartData(lpi.particles[j]));
              }
            }

            tp->addToParticleInterLocal(activePart*lpi.numParticles);
            // put bucket in interaction list
            Vector3D<double> &off = lpi.offset;
            {
              Vector3D<double> &vec = lpi.offset;
              CkPrintf("[%d]: bucket %d with part %ld (%d), (%1.0f,%1.0f,%1.0f)\n", thisIndex, b, lpi.key, gpuIndex, vec.x, vec.y, vec.z);
            }
            state->partLists.push_back(ILPart(gpuIndex, encodeOffset(0, off.x, off.y, off.z), lpi.numParticles));
            // insert a new bucket in the list of targets only when the following condition is met:
            if(state->bucketsParts.length() == 0 || state->bucketsParts[state->bucketsParts.length()-1] != b){
              insertBucketPart(state, b, bucketStartIndex, bucketSize);
            }
            if(state->partLists.length() >= cudaPartThreshold){
              // enough nodes to offload
              sendPartInteractionsToGpu(state, tp, numInteractions == 0);
              resetCudaPartState(state);
            }

#ifdef CHANGA_REFACTOR_PRINT_INTERACTIONS
            Vector3D<double> &vec = lpi.offset;
            CkPrintf("[%d]: bucket %d with %ld (%1.0f,%1.0f,%1.0f)\n", thisIndex, b, lpi.key, vec.x, vec.y, vec.z);
#endif
          }
        }
      }// level

    }// active
  }// bucket
  // here, we have processed all requested buckets 
  // now, check whether all active buckets have been processed,
  // and if so, flush remaining lists to GPU since no more 
  // local/rnr work will be received.
  if(!resume && numFinLocalRNRBuckets == 2*numActiveBuckets){
    // either a local or remote-no-resume walk
    // must send out remaining node/particle interactions here.
    numFinLocalRNRBuckets = 0;
    sendNodeInteractionsToGpu(state, tp, numInteractions == 0);
    resetCudaNodeState(state);
    sendPartInteractionsToGpu(state, tp, numInteractions == 0);
    resetCudaPartState(state);
  }
  // if this is a remote-resume walk, leftovers will be
  // shipped off when there are no more outstanding
  // missed particles/nodes


  // *********************
#endif // ifndef CUDA
}

#ifdef CUDA
#include "HostCUDA.h"


void DeleteHostMoments(CudaMultipoleMoments *array){
  delete [] array;
}

void DeleteHostParticles(CompactPartData *array){
  delete [] array;
}

void ListCompute::initCudaState(DoubleWalkState *state, int numBuckets, int nodeThreshold, int partThreshold, bool resume){
	state->nodeThreshold = nodeThreshold;
	state->partThreshold = partThreshold;
	state->cellLists.reserve(nodeThreshold);
	state->partLists.reserve(partThreshold);

	state->nodeInteractionBucketMarkers.reserve(numBuckets+1);
	state->partInteractionBucketMarkers.reserve(numBuckets+1);

	state->bucketsNodes.reserve(numBuckets+1);
	state->bucketStartMarkersNodes.reserve(numBuckets+1);
	state->bucketSizesNodes.reserve(numBuckets+1);

	state->bucketsParts.reserve(numBuckets+1);
	state->bucketStartMarkersParts.reserve(numBuckets+1);
	state->bucketSizesParts.reserve(numBuckets+1);

	resetCudaNodeState(state);
	resetCudaPartState(state);

	state->resume = resume;
}

void ListCompute::resetCudaNodeState(DoubleWalkState *state){
	state->cellLists.length() = 0;
	state->nodeInteractionBucketMarkers.length() = 0;
	state->bucketsNodes.length() = 0;
	state->bucketStartMarkersNodes.length() = 0;
	state->bucketSizesNodes.length() = 0;
}

void ListCompute::resetCudaPartState(DoubleWalkState *state){
	state->partLists.length() = 0;
	state->partInteractionBucketMarkers.length() = 0;
	state->bucketsParts.length() = 0;
	state->bucketStartMarkersParts.length() = 0;
	state->bucketSizesParts.length() = 0;
}

void ListCompute::getBucketParameters(TreePiece *tp, int bucket, int &bucketStart, int &bucketSize, std::map<NodeKey, int>&lpref){
	// bucket b is listed in this offload
	GenericTreeNode *bucketNode = tp->bucketList[bucket];

	bucketSize = bucketNode->lastParticle - bucketNode->firstParticle + 1;
	NodeKey partKey = bucketNode->getKey();
	partKey <<= 1;

	std::map<NodeKey, int>::iterator iter = lpref.find(partKey);
	CkAssert(iter != lpref.end());
	bucketStart = iter->second;
	//CkPrintf("[%d] bucket %d (key %ld) start: %d\n", tp->getIndex(), bucket, partKey >> 1, bucketStart);
	CkAssert(bucketStart >= 0);
}

void ListCompute::insertBucketNode(DoubleWalkState *state, int b, int start, int size){
	state->bucketsNodes.push_back(b);
	state->bucketSizesNodes.push_back(size);
	state->bucketStartMarkersNodes.push_back(start);
	int currentIntListLength = state->cellLists.length();
	// an element has already been added to the list by the time we get here
	state->nodeInteractionBucketMarkers.push_back(currentIntListLength-1);
}

void ListCompute::insertBucketPart(DoubleWalkState *state, int b, int start, int size){
	state->bucketsParts.push_back(b);
	state->bucketSizesParts.push_back(size);
	state->bucketStartMarkersParts.push_back(start);
	int currentIntListLength = state->partLists.length();
	// an element has already been added to the list by the time we get here
	state->partInteractionBucketMarkers.push_back(currentIntListLength-1);
}

void cudaCallback(void *param, void *msg){
  CudaRequest *data = (CudaRequest *)param;
  int *affectedBuckets = data->affectedBuckets;
  TreePiece *tp = (TreePiece*)data->tp;
  DoubleWalkState *state = (DoubleWalkState *)data->state;
  int bucket;

  int numBucketsDone = data->numBucketsPlusOne-1;
  if(!data->lastBucketComplete){
    numBucketsDone--;
  }
  
  // bucket computations finished
  for(int i = 0; i < numBucketsDone; i++){
    bucket = affectedBuckets[i];
    state->counterArrays[0][bucket]--;
    tp->finishBucket(bucket);
  }

  // free data structures 
  free(data->affectedBuckets);
  delete ((CkCallback *)data->cb);
  delete data; 
}

void ListCompute::sendNodeInteractionsToGpu(DoubleWalkState *state, TreePiece *tp, bool lastBucketComplete){
  // all interactions to be shipped have been recorded at this point
  // need to cap the end of state->nodeInteractionBucketMarkers
  int length = state->cellLists.length();
  state->nodeInteractionBucketMarkers.push_back(length);
  if(length > 0){
    /*
       int index = tp->getIndex();
       CkPrintf("[%d] %d cellLists: ", index, getOptType());
       for(int i = 0; i < state->cellLists.length(); i++){
       Vector3D<double> v = tp->decodeOffset(state->cellLists[i].offsetID);
       CkPrintf("%d (%1.0f,%1.0f,%1.0f), ", state->cellLists[i].index, v.x, v.y, v.z);
       }
       CkPrintf("\nnodeInteractionBucketMarkers: ");
    // state->nodeInteractionBucketMarkers
    for(int i = 0; i < state->nodeInteractionBucketMarkers.length(); i++){
    CkPrintf("%d, ", state->nodeInteractionBucketMarkers[i]);
    }

    CkPrintf("\nbucketSizes: ");
    for(int i = 0; i < state->bucketSizesNodes.length(); i++){
    CkPrintf("%d, ", state->bucketSizesNodes[i]);
    }
    CkPrintf("\nbucketStartMarkersNodes: ");
    for(int i = 0; i < state->bucketStartMarkersNodes.length(); i++){
    CkPrintf("%d, ", state->bucketStartMarkersNodes[i]);
    }
    CkPrintf("\n");
    */


    ILCell *cellList = state->cellLists.getVec();
    int *cellListBucketMarkers = state->nodeInteractionBucketMarkers.getVec();
    int *bucketStarts = state->bucketStartMarkersNodes.getVec();
    int *bucketSizes = state->bucketSizesNodes.getVec();
    int numBuckets = state->bucketSizesNodes.length();
    int numInteractions = state->cellLists.length();

    int total = lastBucketComplete ? numBuckets : numBuckets-1;  
    int bucket;
    CudaRequest *data = new CudaRequest;
    // offloading computation 
    for(int i = 0; i < total; i++){
      bucket = data->affectedBuckets[i];
      state->counterArrays[0][bucket]++;
    }

    data->list = cellList;
    data->bucketMarkers = cellListBucketMarkers;
    data->bucketStarts = bucketStarts;
    data->bucketSizes = bucketSizes;
    data->numInteractions = numInteractions;
    data->numBucketsPlusOne = numBuckets+1;
    data->tp = (void *)tp;
    data->affectedBuckets = state->bucketsNodes.getVec();
    data->lastBucketComplete = lastBucketComplete;
    data->state = state;

    OptType type = getOptType();
    data->cb = new CkCallback(cudaCallback, data);
    if(type == Local){
      TreePieceCellListDataTransferLocal(data);
    }
    else if(!state->resume){
      TreePieceCellListDataTransferRemote(data);
    }
    else{
      CudaMultipoleMoments *missedNodes = state->nodes->getVec();
      int len = state->nodes->length();
      CkAssert(missedNodes);
      TreePieceCellListDataTransferRemoteResume(data, missedNodes, len);
    }

  }
}

void ListCompute::sendPartInteractionsToGpu(DoubleWalkState *state, TreePiece *tp, bool lastBucketComplete){

  int index = tp->getIndex();
  // all interactions to be shipped have been recorded at this point
  // need to cap the end of state->partInteractionBucketMarkers
  int length = state->partLists.length();
  state->partInteractionBucketMarkers.push_back(length);

  if(length > 0){
    /*
       CkPrintf("[%d] %d partLists: ", index, getOptType());
       for(int i = 0; i < state->partLists.length(); i++){
       Vector3D<double> v = tp->decodeOffset(state->partLists[i].off);
       CkPrintf("%d (%1.0f,%1.0f,%1.0f), ", state->partLists[i].index, v.x, v.y, v.z);
       }
       CkPrintf("\npartInteractionBucketMarkers: ");
    // state->partInteractionBucketMarkers
    for(int i = 0; i < state->partInteractionBucketMarkers.length(); i++){
    CkPrintf("%d, ", state->partInteractionBucketMarkers[i]);
    }

    CkPrintf("\nbucketSizes: ");
    for(int i = 0; i < state->bucketSizesParts.length(); i++){
    CkPrintf("%d, ", state->bucketSizesParts[i]);
    }
    CkPrintf("\nbucketStartMarkersParts: ");
    for(int i = 0; i < state->bucketStartMarkersParts.length(); i++){
    CkPrintf("%d, ", state->bucketStartMarkersParts[i]);
    }
    CkPrintf("\n");
    */
    ILPart *partList = state->partLists.getVec();
    int numInteractions = state->partLists.length();

    int *partListBucketMarkers = state->partInteractionBucketMarkers.getVec();
    int *bucketStarts = state->bucketStartMarkersParts.getVec();
    int *bucketSizes = state->bucketSizesParts.getVec();
    int numBuckets= state->bucketSizesParts.length();

    int total = lastBucketComplete ? numBuckets : numBuckets-1;  
    int bucket;
    CudaRequest *data = new CudaRequest;
    // offloading computation 
    for(int i = 0; i < total; i++){
      bucket = data->affectedBuckets[i];
      state->counterArrays[0][bucket]++;
    }

    data->list = partList;
    data->bucketMarkers = partListBucketMarkers;
    data->bucketStarts = bucketStarts;
    data->bucketSizes = bucketSizes;
    data->numInteractions = numInteractions;
    data->numBucketsPlusOne = numBuckets+1;
    data->tp = tp;
    data->affectedBuckets = state->bucketsParts.getVec();
    data->lastBucketComplete = lastBucketComplete;
    data->state = state;

    OptType type = getOptType();
    data->cb = new CkCallback(cudaCallback, data);
    if(type == Local){
      TreePiecePartListDataTransferLocal(data);
    }
    else if(!state->resume){
      TreePiecePartListDataTransferRemote(data);
    }
    else{
      CompactPartData *missedParts = state->particles->getVec();
      int len = state->particles->length();
      CkAssert(missedParts);
      TreePiecePartListDataTransferRemoteResume(data, missedParts, len);
    }
    /*
       OptType type = getOptType();
       if(type == Local){
       TreePiecePartListDataTransferLocal(data);
       }
       else if(type == Remote && !state->resume){
       CompactPartData *missedParts = state->particles->getVec();
       CkAssert(missedParts);
       int len = state->particles->length();
       TreePiecePartListDataTransferRemote(data, missedParts, len);
       }
       else if(type == Remote && state->resume){
       CompactPartData *missedParts = state->particles->getVec();
       CkAssert(missedParts);
       int len = state->particles->length();
       TreePiecePartListDataTransferRemoteResume(data, missedParts, len);
       }
       */
  }
}

#endif

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
