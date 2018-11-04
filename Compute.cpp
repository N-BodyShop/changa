#include "ParallelGravity.h"
#include "GenericTreeNode.h"
//#include "codes.h"

#include "Opt.h"
#include "Compute.h"
#include "TreeWalk.h"
#include "State.h"
#include "Space.h"
#include "gravity.h"

int decodeReqID(int reqID);

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
  // 1. remaining Chunk[] - sized numChunks
  s->counterArrays[0] = new int [dim1];
  s->counterArrays[1] = new int [dim2];
  s->currentBucket = 0;
  s->bWalkDonePending = 0;
  
  // this variable shouldn't be used at all in the remote walk
  s->myNumParticlesPending = -1;
  return s;
}

State *Compute::getNewState(int dim1){
  // 0. local component of numAdditionalRequests, init to 1
  State *s = new State();
  s->counterArrays[0] = new int [dim1];
  s->counterArrays[1] = 0;
  s->currentBucket = 0;
  s->bWalkDonePending = 0;

  // this is used by local walks 
  // not prefetch ones, even though
  // prefetch computes use this version
  // of gtNewState
  s->myNumParticlesPending = dim1;
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
/// @brief Version that frees a DoubleWalkState
void ListCompute::freeState(State *s){
  freeDoubleWalkState((DoubleWalkState *)s);
  Compute::freeState(s);
}

/// @brief Free up all the lists.
void ListCompute::freeDoubleWalkState(DoubleWalkState *state){
  for(int i = 0; i < INTERLIST_LEVELS; i++){
    state->undlists[i].free();
    state->clists[i].free();
  }
  delete [] state->chklists;
  state->undlists.free();
  state->clists.free();

  if(state->rplists.length() > 0){
    for(int i = 0; i < INTERLIST_LEVELS; i++){
      state->rplists[i].free();
    }
    state->rplists.free();
  }
  else if(state->lplists.length() > 0){
    for(int i = 0; i < INTERLIST_LEVELS; i++){
      state->lplists[i].free();
    }
    state->lplists.free();
  }

#ifdef CUDA
  state->nodeLists.free();
  state->particleLists.free();
#endif

  if(state->placedRoots){
    delete [] state->placedRoots;
    state->placedRoots = 0;
  }

}

DoubleWalkState *ListCompute::allocDoubleWalkState(){
  DoubleWalkState *s = new DoubleWalkState;
  s->level = 0;
  s->chklists = new CheckList[INTERLIST_LEVELS];
  s->undlists.resize(INTERLIST_LEVELS);
  s->clists.resize(INTERLIST_LEVELS);

  if(getOptType() == Remote){
    s->rplists.resize(INTERLIST_LEVELS);
  }
  else if(getOptType() == Local){
    s->lplists.resize(INTERLIST_LEVELS);
  }

  return s;
}

/// @brief Version that allocates a DoubleWalkState
State *ListCompute::getNewState(int d1, int d2){
  DoubleWalkState *s = allocDoubleWalkState();
  s->counterArrays[0] = new int [d1];
  s->counterArrays[1] = new int [d2];
  // one boolean for each chunk
  s->placedRoots = new bool [d2];
  s->currentBucket = 0;
  s->bWalkDonePending = 0;

  s->myNumParticlesPending = -1;
  return s;
}

/// @brief Version that allocates a DoubleWalkState
State *ListCompute::getNewState(int d1){
  DoubleWalkState *s = allocDoubleWalkState();
  s->counterArrays[0] = new int [d1];
  s->counterArrays[1] = 0;
  // no concept of chunks in local computation
  s->placedRoots = new bool [1];
  s->currentBucket = 0;
  s->bWalkDonePending = 0;

  s->myNumParticlesPending = d1;
  return s;
}

/// @brief Version that allocates a DoubleWalkState
State *ListCompute::getNewState(){
  DoubleWalkState *s = allocDoubleWalkState();
  s->counterArrays[0] = 0;
  s->counterArrays[1] = 0;
  // this function used for remote-resume states, 
  // no placedRoots vars required
  s->placedRoots = 0;
  s->currentBucket = 0;
  s->bWalkDonePending = 0;

  s->myNumParticlesPending = -1;
  return s;
}

/// @brief Clear lists.
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
/// @brief Reassociate the target node.
void ListCompute::reassoc(void *ce, int ar, Opt *o){
  computeEntity = ce;
  activeRung = ar;
  opt = o;
}
#endif

void GravityCompute::nodeMissedEvent(int reqID, int chunk, State *state, TreePiece *tp){
  if(getOptType() == Remote){
    state->counterArrays[0][decodeReqID(reqID)]++;
    state->counterArrays[1][chunk]++;
  }
}

void PrefetchCompute::startNodeProcessEvent(State *state){
  state->counterArrays[0][0]++;
}

void PrefetchCompute::finishNodeProcessEvent(TreePiece *owner, State *state){
  int save = --state->counterArrays[0][0];
  if(save == 0){
    owner->startRemoteChunk();
  }
}

#if INTERLIST_VER > 0
/// @brief Update state on a node miss.
/// @param reqID unused.
void ListCompute::nodeMissedEvent(int reqID, int chunk, State *state, TreePiece *tp){
  CkAssert(getOptType() == Remote);
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
#endif

int GravityCompute::openCriterion(TreePiece *ownerTP,
                          GenericTreeNode *node, int reqID, State *state){
  return
    openCriterionBucket(node,(GenericTreeNode *)computeEntity,ownerTP->decodeOffset(reqID));

}

void GravityCompute::recvdParticles(ExternalGravityParticle *part,int num,int chunk,int reqID,State *state,TreePiece *tp, Tree::NodeKey &remoteBucket){
  //TreePiece *tp = tw->getOwnerTP();

  Vector3D<cosmoType> offset = tp->decodeOffset(reqID);
  int reqIDlist = decodeReqID(reqID);
  CkAssert(num > 0);
  state->counterArrays[0][reqIDlist] -= 1;
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
  tp->finishBucket(reqIDlist);
  CkAssert(state->counterArrays[1][chunk] >= 0);
  if (state->counterArrays[1][chunk] == 0) {
    cacheGravPart[CkMyPe()].finishedChunk(chunk, tp->particleInterRemote[chunk]);
#ifdef CHECK_WALK_COMPLETIONS
    CkPrintf("[%d] finishedChunk %d GravityCompute::recvdParticles\n", tp->getIndex(), chunk);
#endif
    tp->finishedChunk(chunk);
  }

}

void PrefetchCompute::recvdParticles(ExternalGravityParticle *egp,int num,int chunk,int reqID,State *state, TreePiece *tp, Tree::NodeKey &remoteBucket){
        // wait for prefetched particles as well
        // this way, all nodes/parts not missed will be handled by RNR
        // and all those missed, by RR
        // if we didn't wait for particles
        // it could transpire that we don't miss on these particles during RNR, but they aren't present on the gpu
        // and so we'd have to have a separate array of missed particles for the RNR (in much the same way that the RR has
        // separate arrays for missed nodes and particles) 
  finishNodeProcessEvent(tp, state);
}

void GravityCompute::nodeRecvdEvent(TreePiece *owner, int chunk, State *state, int reqIDlist){

  state->counterArrays[0][reqIDlist]--;
  owner->finishBucket(reqIDlist);

  CkAssert(chunk >= 0);
  state->counterArrays[1][chunk] --;
  CkAssert(state->counterArrays[1][chunk] >= 0);
  if (state->counterArrays[1][chunk] == 0) {
    cacheGravPart[CkMyPe()].finishedChunk(chunk, owner->particleInterRemote[chunk]);
#ifdef CHECK_WALK_COMPLETIONS
    CkPrintf("[%d] finishedChunk %d GravityCompute::nodeRecvdEvent\n", owner->getIndex(), chunk);
#endif
    owner->finishedChunk(chunk);
  }// end if finished with chunk
}

#if INTERLIST_VER > 0
/// @brief Update state upon receiving a remote node.
///
/// Calls TreePiece::finishedChunk() if all outstanding requests are satisfied.
void ListCompute::nodeRecvdEvent(TreePiece *owner, int chunk, State *state, int reqIDlist){
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
#if COSMO_PRINT_BK > 1
  CkPrintf("[%d] nodeRecvdEvent chunk: %d remainingChunk: %d\n", owner->getIndex(), chunk, remainingChunk);
#endif
  if (remainingChunk == 0) {
#ifdef CUDA
    // no more nodes/particles are going to be delivered by the cache
    // flush the interactions remaining in the state
    DoubleWalkState *ds = (DoubleWalkState *)state;

    if(ds->nodeLists.totalNumInteractions > 0){
      sendNodeInteractionsToGpu(ds, owner);
      resetCudaNodeState(ds);
    }
    if(ds->particleLists.totalNumInteractions > 0){
      sendPartInteractionsToGpu(ds, owner);
      resetCudaPartState(ds);
    }

#endif
#if COSMO_PRINT_BK > 1
    CkPrintf("[%d] FINISHED CHUNK %d from nodeRecvdEvent\n", owner->getIndex(), chunk);
#endif
    cacheGravPart[CkMyPe()].finishedChunk(chunk, owner->particleInterRemote[chunk]);
#ifdef CHECK_WALK_COMPLETIONS
    CkPrintf("[%d] finishedChunk %d ListCompute::nodeRecvdEvent\n", owner->getIndex(), chunk);
#endif
    owner->finishedChunk(chunk);
  }// end if finished with chunk
}
#endif

#include "TreeNode.h"
using namespace TreeStuff;

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
  CkAssert(opt != NULL);

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
    
    GenericTreeNode *b = (GenericTreeNode *)computeEntity;
    updateInterMass(b->particlePointer,b->firstParticle,b->lastParticle,node->moments.totalMass);

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
      int computed = 0;
#ifdef BENCHMARK_TIME_COMPUTE
      double startTime = CmiWallTimer();
#endif
      Vector3D<cosmoType> offset = tp->decodeOffset(reqID);
      for(int i = node->firstParticle; i <= node->lastParticle; i++){
        computed += partBucketForce(
                                  &part[i-node->firstParticle],
                                  (GenericTreeNode *)computeEntity,
                                  tp->getParticles(),
                                  offset,
                                  activeRung);

        GenericTreeNode *b = (GenericTreeNode *)computeEntity;
        updateInterMass(b->particlePointer,b->firstParticle,b->lastParticle,&part[i-node->firstParticle],offset);
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
    part = tp->particlesMissed(keyref,
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
  CkAbort("bad walk state");
  return -1;
}

void GravityCompute::updateInterMass(GravityParticle *p, int start, int end, double totalMass){
  for(int j = start; j <= end; j++){
    p[j-start].interMass += totalMass;
  }
}

void GravityCompute::updateInterMass(GravityParticle *p, int start, int end, GravityParticle *s, Vector3D<cosmoType> &offset){
  Vector3D<cosmoType> r; 
  for(int j = start; j <= end; j++){
    r = offset + s->position - p[j-start].position;
    if(r.lengthSquared() == 0) continue;
    p[j-start].interMass += s->mass;
  }
}

///
/// @brief Calculate forces due to particles in a node
/// @param ownerTP TreePiece of target particles
/// @param node GenericTreeNode containing source particles
/// @param part array of source particles
/// @param reqID bucket number of target particles and offset
///
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

int PrefetchCompute::openCriterion(TreePiece *ownerTP,
                          GenericTreeNode *node, int reqID, State *state){
  TreePiece *tp = ownerTP;
  PrefetchRequestStruct prs(tp->prefetchReq, tp->numPrefetchReq);
  Vector3D<cosmoType> offset = ownerTP->decodeOffset(reqID);
  for(int i = 0; i < prs.numPrefetchReq; i++){
    BinaryTreeNode testNode;
    testNode.boundingBox = prs.prefetchReq[i];

    if(openCriterionBucket(node, &testNode, offset))
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
    ExternalGravityParticle *part;
    part = tp->particlesMissed(keyref,
                               chunk,
                               node->remoteIndex,
                               node->firstParticle,
                               node->lastParticle,
                               reqID,
                               true, awi, (void *)0);

    if(part == NULL){
      // waiting for particles for reasons discussed 
      // in comments within recvdParticles
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
  CkAbort("PrefetchCompute: bad option");
  return -1;
}

#if INTERLIST_VER > 0
/// @brief Process a node.
/// @param node is the global node being processed.
/// @param state contains the lists to be checked.
/// @param chunk chunk we are walking; used in case of a miss
/// @param reqID Encodes offset and bucket
/// @param awi Active walk index; used in case of a miss
/// @return KEEP if we descend further down the tree
int ListCompute::doWork(GenericTreeNode *node, TreeWalk *tw, State *state, int chunk, int reqID, bool isRoot, bool &didcomp, int awi){

  DoubleWalkState *s = (DoubleWalkState *)state;
  int level = s->level;
  CheckList &chklist = s->chklists[level];
  UndecidedList &undlist = s->undlists[level];

  TreePiece *tp = tw->getOwnerTP();
  Vector3D<cosmoType> offset = tp->decodeOffset(reqID);

  if(node->getType() == Empty || node->getType() == CachedEmpty){
#ifdef CHANGA_REFACTOR_WALKCHECK_INTERLIST
    if(node->parent->getType() != Boundary || getOptType() == Local){
      addNodeToInt(node, reqID, s);
    }
#endif
    return DUMP;
  }
  // check opening criterion
  int open;

  open = openCriterion(tp, node, reqID, state);

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

    addNodeToInt(node, reqID, s);
    // all particles beneath this node have been
    // scheduled for computation
    computed = node->lastParticle-node->firstParticle+1;
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
#if defined CHANGA_REFACTOR_PRINT_INTERACTIONS || defined CHANGA_REFACTOR_WALKCHECK_INTERLIST || defined CUDA
    NodeKey key = node->getKey();
    addLocalParticlesToInt(part, computed, offset, s, key, node);
#else
    addLocalParticlesToInt(part, computed, offset, s);
#endif

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
    part = tp->particlesMissed(keyref,
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

#if defined CHANGA_REFACTOR_PRINT_INTERACTIONS || defined CHANGA_REFACTOR_WALKCHECK_INTERLIST || defined CUDA
      NodeKey key = node->getKey();
      addRemoteParticlesToInt(part, computed, offset, s, key);
#else
      addRemoteParticlesToInt(part, computed, offset, s);
#endif
    }
    else{
#if CHANGA_REFACTOR_DEBUG > 2
      CkPrintf("Particles not found in cache\n");
#endif
      CkAssert(getOptType() == Remote);
      // particles missed
      int start, end;
      GenericTreeNode *source = (GenericTreeNode *)computeEntity;
      tp->getBucketsBeneathBounds(source, start, end);
#if COSMO_PRINT_BK > 1
      CkPrintf("[%d] missed parts %ld (chunk %d)\n", tp->getIndex(), keyref << 1, chunk);
#endif
      tp->updateUnfinishedBucketState(start, end, 1, chunk, state);
#ifdef CHANGA_REFACTOR_MEMCHECK
      CkPrintf("memcheck after particlesmissed (%ld)\n", keyref);
      CmiMemoryCheck();
#endif
    }
    return DUMP;
  }
  else if(action == DUMP || action == NOP){
    return DUMP;
  }
  CkAbort("ListCompute: bad walk state");
  return -1;
}

/// @brief Process received remote particles
/// @param part Array of particles received.
/// @param num Number of particles
///
/// Update the state bookkeeping, add the particles to the interaction
/// list and call stateReady() to compute their interactions.  Call
/// TreePiece::finishedChunk() if all outstanding requests are satisfied.
void ListCompute::recvdParticles(ExternalGravityParticle *part,int num,int chunk,int reqID,State *state_,TreePiece *tp, Tree::NodeKey &remoteBucket){

  Vector3D<cosmoType> offset = tp->decodeOffset(reqID);
  CkAssert(num > 0);
  GenericTreeNode *source = (GenericTreeNode *)computeEntity;

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
#if defined CHANGA_REFACTOR_PRINT_INTERACTIONS || defined CHANGA_REFACTOR_WALKCHECK_INTERLIST || defined CUDA
  NodeKey key = remoteBucket >> 1;
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
#if COSMO_PRINT_BK > 1
  CkPrintf("[%d] recvdParticles chunk: %d remainingChunk: %d\n", tp->getIndex(), chunk, remainingChunk);
#endif
  CkAssert(remainingChunk >= 0);
  if (remainingChunk == 0) {
#ifdef CUDA
    if(state->nodeLists.totalNumInteractions > 0){
      sendNodeInteractionsToGpu(state, tp);
      resetCudaNodeState(state);
    }
    if(state->particleLists.totalNumInteractions > 0){
      sendPartInteractionsToGpu(state, tp);
      resetCudaPartState(state);
    }
#endif
#if COSMO_PRINT_BK > 1
    CkPrintf("[%d] FINISHED CHUNK %d from recvdParticles\n", tp->getIndex(), chunk);
#endif
    cacheGravPart[CkMyPe()].finishedChunk(chunk, tp->particleInterRemote[chunk]);
#ifdef CHECK_WALK_COMPLETIONS
    CkPrintf("[%d] finishedChunk %d ListCompute::recvdParticles\n", tp->getIndex(), chunk);
#endif
    tp->finishedChunk(chunk);
  }
}

/// @brief apply node opening criterion to this node.
int ListCompute::openCriterion(TreePiece *ownerTP,
                          GenericTreeNode *node, int reqID, State *state){
  return
    openCriterionNode(node,(GenericTreeNode *)computeEntity, ownerTP->decodeOffset(reqID));
}

#if defined CHANGA_REFACTOR_PRINT_INTERACTIONS || defined CHANGA_REFACTOR_WALKCHECK_INTERLIST || defined CUDA
void ListCompute::addRemoteParticlesToInt(ExternalGravityParticle *parts, int n, Vector3D<cosmoType> &offset, DoubleWalkState *s, NodeKey key){
#else
void ListCompute::addRemoteParticlesToInt(ExternalGravityParticle *parts, int n, Vector3D<cosmoType> &offset, DoubleWalkState *s){
#endif
  RemotePartInfo rpi;
  int level = s->level;

  rpi.particles = parts;
  rpi.numParticles = n;
  rpi.offset = offset;
#if defined CHANGA_REFACTOR_PRINT_INTERACTIONS || defined CHANGA_REFACTOR_WALKCHECK_INTERLIST || defined CUDA
  rpi.key = key;
#endif

  s->rplists[level].push_back(rpi);
}

#if defined CHANGA_REFACTOR_PRINT_INTERACTIONS || defined CHANGA_REFACTOR_WALKCHECK_INTERLIST  || defined CUDA
void ListCompute::addLocalParticlesToInt(GravityParticle *parts, int n, Vector3D<cosmoType> &offset, DoubleWalkState *s, NodeKey key, GenericTreeNode *gtn){
#else
void ListCompute::addLocalParticlesToInt(GravityParticle *parts, int n, Vector3D<cosmoType> &offset, DoubleWalkState *s){
#endif
  LocalPartInfo lpi;
  int level = s->level;

  lpi.particles = parts;
  lpi.numParticles = n;
  lpi.offset = offset;
#if defined CHANGA_REFACTOR_PRINT_INTERACTIONS || defined CHANGA_REFACTOR_WALKCHECK_INTERLIST || defined CUDA
  lpi.key = key;
  lpi.nd = gtn;
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

#ifdef CUDA
#include "DataManager.h"
#endif

#ifdef CUDA
void allocatePinnedHostMemory(void **, size_t);

template<typename T>
CudaRequest *GenericList<T>::serialize(TreePiece *tp){
    // get count of buckets with interactions first
    int numFilledBuckets = 0;
    int listpos = 0;
    int curbucket = 0;

#ifdef HAPI_TRACE
    double starttime = CmiWallTimer();
#endif
    for(int i = 0; i < lists.length(); i++){
      if(lists[i].length() > 0){
        numFilledBuckets++;
      }
    }

    // create flat lists and associated data structures
    // allocate memory and flatten lists only if there
    // are interactions to transfer. 
    T *flatlists = NULL;
    int *markers = NULL;
    int *starts = NULL;
    int *sizes = NULL;
    int *affectedBuckets = NULL;

    if(totalNumInteractions > 0){
#ifdef HAPI_USE_CUDAMALLOCHOST
      allocatePinnedHostMemory((void **)&flatlists, totalNumInteractions*sizeof(T));
      allocatePinnedHostMemory((void **)&markers, (numFilledBuckets+1)*sizeof(int));
      allocatePinnedHostMemory((void **)&starts, (numFilledBuckets)*sizeof(int));
      allocatePinnedHostMemory((void **)&sizes, (numFilledBuckets)*sizeof(int));
#else
      flatlists = (T *) malloc(totalNumInteractions*sizeof(T));
      markers = (int *) malloc((numFilledBuckets+1)*sizeof(int));
      starts = (int *) malloc(numFilledBuckets*sizeof(int));
      sizes = (int *) malloc(numFilledBuckets*sizeof(int));
#endif
      affectedBuckets = new int[numFilledBuckets];

      // populate flat lists
      int listslen = lists.length();
      for(int i = 0; i < listslen; i++){
        int listilen = lists[i].length();
        if(listilen > 0){
          memcpy(&flatlists[listpos], lists[i].getVec(), listilen*sizeof(T));
          markers[curbucket] = listpos;
          if(tp->largePhase()){
            getBucketParameters(tp, i, starts[curbucket], sizes[curbucket]);
          }
          else{
            getActiveBucketParameters(tp, i, starts[curbucket], sizes[curbucket]);
          }
          affectedBuckets[curbucket] = i;
          listpos += listilen;
          curbucket++;
        }
      }
      markers[numFilledBuckets] = listpos;
      CkAssert(listpos == totalNumInteractions);
    }

    CudaRequest *request = new CudaRequest;
    request->list = (void *)flatlists;
    request->bucketMarkers = markers;
    request->bucketStarts = starts;
    request->bucketSizes = sizes;
    request->numInteractions = totalNumInteractions;
    request->numBucketsPlusOne = numFilledBuckets+1;
    request->affectedBuckets = affectedBuckets;
    request->tp = (void *)tp;
    request->fperiod = tp->fPeriod.x;

#ifdef HAPI_TRACE
    traceUserBracketEvent(CUDA_SER_LIST, starttime, CmiWallTimer());
#endif

    return request;
  }

// specialization for particle interlists
// need to send interactions with particles instead of those with buckets
// so that we can allot interactions within a bucket-bucket interaction
// to different threads 

template<>
CudaRequest *GenericList<ILPart>::serialize(TreePiece *tp){
    // get count of buckets with interactions first
    int numFilledBuckets = 0;
    int listpos = 0;
    int curbucket = 0;

#ifdef HAPI_TRACE
    double starttime = CmiWallTimer();
#endif
    int numParticleInteractions = 0;
    for(int i = 0; i < lists.length(); i++){
      if(lists[i].length() > 0){
        numFilledBuckets++;
        for(int j = 0; j < lists[i].length(); j++){
          numParticleInteractions += lists[i][j].num;
        }
      }
    }

    // create flat lists and associated data structures
    // allocate memory and flatten lists only if there
    // are interactions to transfer. 
    // use ILCell instead of ILPart because we don't need num field anymore
    // however, we pay the price of replicating the offset for each particle in the bucket 
    ILCell *flatlists = NULL;
    int *markers = NULL;
    int *starts = NULL;
    int *sizes = NULL;
    int *affectedBuckets = NULL;

    if(totalNumInteractions > 0){
#ifdef HAPI_USE_CUDAMALLOCHOST
      allocatePinnedHostMemory((void **)&flatlists, numParticleInteractions*sizeof(ILCell));
      allocatePinnedHostMemory((void **)&markers, (numFilledBuckets+1)*sizeof(int));
      allocatePinnedHostMemory((void **)&starts, (numFilledBuckets)*sizeof(int));
      allocatePinnedHostMemory((void **)&sizes, (numFilledBuckets)*sizeof(int));
#else
      flatlists = (T *) malloc(numParticleInteractions*sizeof(ILCell));
      markers = (int *) malloc((numFilledBuckets+1)*sizeof(int));
      starts = (int *) malloc(numFilledBuckets*sizeof(int));
      sizes = (int *) malloc(numFilledBuckets*sizeof(int));
#endif
      affectedBuckets = new int[numFilledBuckets];

      // populate flat lists
      int listslen = lists.length();
      // for each level i
      for(int i = 0; i < listslen; i++){
        int listilen = lists[i].length();
        if(listilen > 0){
          
          markers[curbucket] = listpos;
          // for each bucket j at level i
          for(int j = 0; j < listilen; j++){

            ILPart &ilp = lists[i][j];
            int bucketLength = ilp.num;
            int bucketStart = ilp.index;
            int offsetID = ilp.off;
            
            // for each particle k in bucket j
            for(int k = 0; k < bucketLength; k++){
              flatlists[listpos].index = bucketStart+k;
              flatlists[listpos].offsetID = offsetID;
              listpos++;
            }
          }
          if(tp->largePhase()){
            getBucketParameters(tp, i, starts[curbucket], sizes[curbucket]);
          }
          else{
            getActiveBucketParameters(tp, i, starts[curbucket], sizes[curbucket]);
          }
          affectedBuckets[curbucket] = i;
          curbucket++;
        }
      }
      markers[numFilledBuckets] = listpos;
      CkAssert(listpos == numParticleInteractions);
    }

    CudaRequest *request = new CudaRequest;
    request->list = (void *)flatlists;
    request->bucketMarkers = markers;
    request->bucketStarts = starts;
    request->bucketSizes = sizes;
    request->numInteractions = numParticleInteractions;
    request->numBucketsPlusOne = numFilledBuckets+1;
    request->affectedBuckets = affectedBuckets;
    request->tp = (void *)tp;
    request->fperiod = tp->fPeriod.x;

#ifdef HAPI_TRACE
    traceUserBracketEvent(CUDA_SER_LIST, starttime, CmiWallTimer());
#endif

    return request;
  }



template<typename T>
void GenericList<T>::push_back(int b, T &ilc, DoubleWalkState *state, TreePiece *tp){
    if(lists[b].length() == 0){
        state->counterArrays[0][b]++;
#if COSMO_PRINT_BK > 1
        CkPrintf("[%d] request out bucket %d numAddReq: %d,%d\n", tp->getIndex(), b, tp->sRemoteGravityState->counterArrays[0][b], tp->sLocalGravityState->counterArrays[0][b]);
#endif
    }
    lists[b].push_back(ilc);
    totalNumInteractions++;
  }
#endif

/// @brief Populate CkLoop Data
///
/// This populates the local, remote lists. Once all the nodes and particles are
/// identified with which force calculations need to be performed, this list is
/// passed on to CkLoop function where it gets parallelized.
void ListCompute::fillLists(State *state_, TreePiece *tp, int chunk, int start,
    int end, CkVec<OffsetNode>& clistforb, CkVec<RemotePartInfo>& rplistforb,
    CkVec<LocalPartInfo>& lplistforb) {
  GravityParticle *particles = tp->getParticles();

  DoubleWalkState *state = (DoubleWalkState *)state_;

  GenericTreeNode *lowestNode = state->lowestNode;
  int maxlevel = lowestNode->getLevel(lowestNode->getKey());

  bool hasRemoteLists = state->rplists.length() > 0 ? true : false;
  bool hasLocalLists = state->lplists.length() > 0 ? true : false;

  int numNodes = 0;
  int numLParticles = 0;
  int numRParticles = 0;
  bool filled = false;

  for(int level = 0; level <= maxlevel; level++){
    // node interactions
    numNodes += state->clists[level].length();

    // remote particle interactions
    if(hasRemoteLists){
      int listlen = state->rplists[level].length();
      for(int part = 0; part < listlen; part++){
        numRParticles += (state->rplists[level])[part].numParticles;
      }
    }
    // local particle interactions
    if(hasLocalLists){
      int listlen = state->lplists[level].length();
      for(int part = 0; part < listlen; part++){
        numLParticles += (state->lplists[level])[part].numParticles;
      }
    }
  }

  for(int b = start; b < end; b++) {
    if(tp->bucketList[b]->rungs >= activeRung) {

      int activePart=0;
      for (int k=tp->bucketList[b]->firstParticle; k<=tp->bucketList[b]->lastParticle; ++k) {
        if (tp->myParticles[k].rung >= activeRung) activePart++;
      }

      OptType type = getOptType();
      if(type == Local){
        tp->addToNodeInterLocal(numNodes*activePart);
        tp->addToParticleInterLocal(numLParticles*activePart);
      }
      else if(type == Remote){
        tp->addToNodeInterRemote(chunk,numNodes*activePart);
        tp->addToParticleInterRemote(chunk,numRParticles*activePart);
      }

      if (!filled) {
        filled = true;
        for(int level = 0; level <= maxlevel; level++){

          CkVec<OffsetNode> &clist = state->clists[level];
          for(unsigned int i = 0; i < clist.length(); i++){
            clistforb.insertAtEnd(clist[i]);
          }

          // remote particles
          if(hasRemoteLists){
            CkVec<RemotePartInfo> &rpilist = state->rplists[level];
            // for each bunch of particles in list
            for(unsigned int i = 0; i < rpilist.length(); i++){
              RemotePartInfo &rpi = rpilist[i];
              rplistforb.insertAtEnd(rpi);
            }
          }

          // local particles
          if(hasLocalLists) {
            CkVec<LocalPartInfo> &lpilist = state->lplists[level];
            for(unsigned int i = 0; i < lpilist.length(); i++) {
              LocalPartInfo &lpi = lpilist[i];
              lplistforb.insertAtEnd(lpi);
            }
          }
        }// level
      }
    }// active
  }// bucket
}

template<class type> int calcNodeForces(TreePiece *tp, int b, int activeRung,
    CkVec<type>& clist) {

  GravityParticle *particles = tp->getParticles();
  int computed = 0;
  for(unsigned int i = 0; i < clist.length(); i++){

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
    computed +=  nodeBucketForce(clist[i].node,
        tp->getBucket(b),
        particles,
        tp->decodeOffset(clist[i].offsetID), activeRung);
#ifdef CHANGA_REFACTOR_PRINT_INTERACTIONS
    if(b == TEST_BUCKET && tp->getIndex() == TEST_TP){
      Vector3D<double> vec = tp->decodeOffset(clist[i].offsetID);
      CkPrintf("[%d]: bucket %d with node %ld (%1.0f,%1.0f,%1.0f)\n", tp->thisIndex,  b, clist[i].node->getKey(), vec.x, vec.y, vec.z);
    }
#endif
  }
  return computed;
}

template<class type> int calcParticleForces(TreePiece *tp, int b, int activeRung,
    CkVec<type>& clist) {

  GravityParticle *particles = tp->getParticles();
  int computed = 0;
  // for each bunch of particles in list
  for(unsigned int i = 0; i < clist.length(); i++){
    type &cli = clist[i];

#ifdef CHANGA_REFACTOR_WALKCHECK_INTERLIST
      NodeKey key = cli.key;
      tp->addToBucketChecklist(b, key);
      tp->combineKeys(key, b);
#endif

    // for each particle in a bunch
    for(int j = 0; j < cli.numParticles; j++){
      computed +=  partBucketForce(&cli.particles[j],
          tp->getBucket(b),
          particles,
          cli.offset, activeRung);

#ifdef CHANGA_REFACTOR_PRINT_INTERACTIONS
      if(b == TEST_BUCKET && tp->getIndex() == TEST_TP){
        Vector3D<double> &vec = cli.offset;
        CkPrintf("[%d]: bucket %d with remote part %ld (%1.0f,%1.0f,%1.0f)\n",
          tp->thisIndex, b, cli.key, vec.x, vec.y, vec.z);
      }
#endif
    }
  }
  return computed;
}

#ifdef GPU_LOCAL_TREE_WALK

void cudaCallbackForAllBuckets(void *param, void *msg) {
  CudaRequest *data = (CudaRequest *)param;
  int *affectedBuckets = data->affectedBuckets;
  TreePiece *tp = (TreePiece*)data->tp;
  DoubleWalkState *state = (DoubleWalkState *)data->state;
  int bucket;

  int numBucketsDone = data->numBucketsPlusOne-1;

  // bucket computations finished
  //
  for(int i = 0; i < numBucketsDone; i++){
    bucket = affectedBuckets[i];
    state->counterArrays[0][bucket]--;
    tp->finishBucket(bucket);
  }

  // free data structures
  if(numBucketsDone > 0){
    delete [] data->affectedBuckets;
  }
  delete ((CkCallback *)data->cb);
  delete data;
}

/**
 * This function is designed to send an ignition signal to the GPU manager.
 * To make minor change to existing ChaNGa code, we mimic a nodeGravityCompute
 * request. This request will eventually call our GPU local tree walk kernel.
 */
void ListCompute::sendLocalTreeWalkTriggerToGpu(State *state, TreePiece *tp,
  int activeRung, int startBucket, int endBucket) {
  int numFilledBuckets = 0;
  for (int i = startBucket; i < endBucket; ++i) {
    if (tp->bucketList[i]->rungs >= activeRung) {
      ++numFilledBuckets;
    }
  }

  // No necessary to call GPU kernel if there is no active bucket in
  // current tree piece
  if (numFilledBuckets == 0) {
    return;
  }

  int *affectedBuckets = new int[numFilledBuckets];

  // Set up a series of dummy parameters to match existing function interfaces
  int dummyTotalNumInteractions = 1;
  int dummyCurBucket = 0;
  ILCell *dummyFlatlists = NULL;
  int *dummyNodeMarkers = NULL;
  int *dummyStarts = NULL;
  int *dummySizes = NULL;

#ifdef HAPI_USE_CUDAMALLOCHOST
  allocatePinnedHostMemory((void **)&dummyFlatlists, dummyTotalNumInteractions *
                                                    sizeof(ILCell));
  allocatePinnedHostMemory((void **)&dummyNodeMarkers, (numFilledBuckets+1) *
                                                      sizeof(int));
  allocatePinnedHostMemory((void **)&dummyStarts, numFilledBuckets *
                                                  sizeof(int));
  allocatePinnedHostMemory((void **)&dummySizes, numFilledBuckets * sizeof(int));
#else
  dummyFlatlists = (ILCell *) malloc(dummy_totalNumInteractions*sizeof(ILCell));
  dummyNodeMarkers = (int *) malloc((numFilledBuckets+1)*sizeof(int));
  dummyStarts = (int *) malloc(numFilledBuckets*sizeof(int));
  dummySizes = (int *) malloc(numFilledBuckets*sizeof(int));
#endif

//  No need to memset the interaction list array since we're not using it at all
//  And, we can't directly memset it.
  ILCell temp_ilc;
  memcpy(&dummyFlatlists[0], &temp_ilc, dummyTotalNumInteractions *
                                        sizeof(ILCell));
  for (int i = startBucket; i < endBucket; ++i) {
    if (tp->bucketList[i]->rungs >= activeRung) {
      ((DoubleWalkState *)state)->counterArrays[0][i] ++;
      dummyNodeMarkers[dummyCurBucket] = tp->bucketList[i]->nodeArrayIndex;
      int tempNum = 0;
      memcpy(&dummyStarts[dummyCurBucket], &tempNum, sizeof(int));
      memcpy(&dummySizes[dummyCurBucket], &tempNum, sizeof(int));
      affectedBuckets[dummyCurBucket] = i;
      dummyCurBucket++;
    }
  }

  CudaRequest *request = new CudaRequest;

  request->numBucketsPlusOne = numFilledBuckets+1;
  request->affectedBuckets = affectedBuckets;
  request->tp = (void *)tp;
  request->state = (void *)state;
  request->node = true;
  request->remote = false;
  request->firstParticle = tp->FirstGPUParticleIndex;
  request->lastParticle = tp->LastGPUParticleIndex;
  // In DataManager serializes the local tree so that the root of the local tree
  // will always be the 0th element in the moments array.
  request->rootIdx = 0;
  request->theta = theta;
  request->thetaMono = thetaMono;
  request->nReplicas = tp->nReplicas;
  request->fperiod = tp->fPeriod.x;
  request->fperiodY = tp->fPeriod.y;
  request->fperiodZ = tp->fPeriod.z;
  request->cb = new CkCallback(cudaCallbackForAllBuckets, request);

  request->list = (void *)dummyFlatlists;
  request->bucketMarkers = dummyNodeMarkers;
  request->bucketStarts = dummyStarts;
  request->bucketSizes = dummySizes;
  request->numInteractions = dummyTotalNumInteractions;

  TreePieceCellListDataTransferLocal(request);
}
#endif //GPU_LOCAL_TREE_WALK

/// @brief Check for computation
/// Computation can be done on buckets indexed from start to end
/// @param state_ State to be checked
/// @param start first bucket to apply computations
/// @param end last + 1 bucket to apply computations
///
/// state->lowestNode is used to determine how deep the walk got.
/// Cells and particles from that level and above are processed
/// on these buckets.
void ListCompute::stateReady(State *state_, TreePiece *tp, int chunk, int start, int end){
  int thisIndex = tp->getIndex();
  GravityParticle *particles = tp->getParticles();

  DoubleWalkState *state = (DoubleWalkState *)state_;

#ifdef CUDA
  bool resume = state->resume;
  int numActiveBuckets = tp->numActiveBuckets;

  // for cached particles
  std::map<NodeKey, int> &cpref = tp->dm->getCachedPartsOnGpuTable();
#endif

  GenericTreeNode *lowestNode = state->lowestNode;
  int maxlevel = lowestNode->getLevel(lowestNode->getKey());

  bool hasRemoteLists = state->rplists.length() > 0 ? true : false;
  bool hasLocalLists = state->lplists.length() > 0 ? true : false;

#if defined CUDA && COSMO_PRINT_BK > 1
  CkPrintf("[%d]: stateReady(%d-%d) %s, resume: %d\n",  tp->getIndex(), start, end, getOptType() == Remote ? "Remote" : "Local", state->resume);
#endif

#ifdef CHANGA_REFACTOR_MEMCHECK
    CkPrintf("memcheck before list iteration\n");
    CmiMemoryCheck();
#endif
#ifndef CUDA
  for(int b = start; b < end; b++){
    if(tp->bucketList[b]->rungs >= activeRung){

      for(int level = 0; level <= maxlevel; level++){

        CkVec<OffsetNode> &clist = state->clists[level];
        int computed;
        computed = calcNodeForces(tp, b, activeRung, clist);
        if(getOptType() == Remote){
          tp->addToNodeInterRemote(chunk, computed);
        } else if(getOptType() == Local){
          tp->addToNodeInterLocal(computed);
        }

        // remote particles
        if(hasRemoteLists){
          CkVec<RemotePartInfo> &rpilist = state->rplists[level];
          computed = calcParticleForces(tp, b, activeRung, rpilist);
          if(getOptType() == Remote){// don't really have to perform this check
            tp->addToParticleInterRemote(chunk, computed);
          }
        }

        // local particles
        if(hasLocalLists){
          CkVec<LocalPartInfo> &lpilist = state->lplists[level];
          computed = calcParticleForces(tp, b, activeRung, lpilist);
          tp->addToParticleInterLocal(computed);
        }
      }// level

    }// active
  }// bucket

#else // else part of ifndef CUDA
  // *************************************************************
  // calculate number of interactions first
  int numNodes = 0;
  int numLParticles = 0;
  int numRParticles = 0;

  for(int level = 0; level <= maxlevel; level++){
    // node interactions
    numNodes += state->clists[level].length();

    // remote particle interactions
    if(hasRemoteLists){
      int listlen = state->rplists[level].length();
      for(int part = 0; part < listlen; part++){
        numRParticles += (state->rplists[level])[part].numParticles;
      }
    }
    // local particle interactions
    if(hasLocalLists){
      int listlen = state->lplists[level].length();
      for(int part = 0; part < listlen; part++){
        numLParticles += (state->lplists[level])[part].numParticles;
      }
    }
  }

  for(int b = start; b < end; b++){
    if(tp->bucketList[b]->rungs >= activeRung){
      
      int activePart=0;
      for (int k=tp->bucketList[b]->firstParticle; k<=tp->bucketList[b]->lastParticle; ++k) {
        if (tp->myParticles[k].rung >= activeRung) activePart++;
      }

      OptType type = getOptType();
      if(type == Local){
        tp->addToNodeInterLocal(numNodes*activePart);
        tp->addToParticleInterLocal(numLParticles*activePart);
      }
      else if(type == Remote){
        tp->addToNodeInterRemote(chunk,numNodes*activePart);
        tp->addToParticleInterRemote(chunk,numRParticles*activePart);
      }

      for(int level = 0; level <= maxlevel; level++){

        CkVec<OffsetNode> &clist = state->clists[level];
        for(int i = 0; i < clist.length(); i++){
          GenericTreeNode *node = clist[i].node;
#ifdef CHANGA_REFACTOR_WALKCHECK_INTERLIST
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
          int index = node->nodeArrayIndex;

#ifdef CHANGA_REFACTOR_PRINT_INTERACTIONS
          {

          if(b == TEST_BUCKET && tp->getIndex() == TEST_TP){
            Vector3D<cosmoType> vec = tp->decodeOffset(clist[i].offsetID);
            CkPrintf("[%d]: remote: %d resume: %d bucket %d with node %ld (%d), (%1.0f,%1.0f,%1.0f)\n", thisIndex, getOptType() == Remote, state->resume, b, clist[i].node->getKey(), clist[i].node->nodeArrayIndex, vec.x, vec.y, vec.z);
          }
          
#endif

#ifdef HEXADECAPOLE
          GenericTreeNode *bucketNode = tp->bucketList[b];
          if (openSoftening(node, bucketNode,
                            tp->decodeOffset(clist[i].offsetID))) {

            ExternalGravityParticle tmpPart;
            MultipoleMoments &m = node->moments;
            tmpPart.mass = m.totalMass;
            tmpPart.soft = m.soft;
            tmpPart.position = m.cm;

            partBucketForce(&tmpPart, bucketNode, particles,
                            tp->decodeOffset(clist[i].offsetID), activeRung);
            continue;
          }
#endif
          DoubleWalkState *rrState;
          if(state->resume || (!state->resume && index < 0)){
            CkAssert(getOptType() == Remote);
            rrState = (DoubleWalkState*) tp->sInterListStateRemoteResume;
            std::unordered_map<NodeKey,int>::iterator it = rrState->nodeMap.find(node->getKey());
            if(it == rrState->nodeMap.end()){
                index = rrState->nodes->push_back_v(CudaMultipoleMoments(node->moments));
                rrState->nodeMap[node->getKey()] = index;
                }
            else{
                index = it->second;
                }
          }
          else{
            rrState = state;
          }
#ifdef CHANGA_REFACTOR_PRINT_INTERACTIONS
          if(b == TEST_BUCKET && tp->getIndex() == TEST_TP){
            CkPrintf("[%d]: pushing node 0x%x (%f,%f,%f,%f,%f,%f)\n", thisIndex, node, node->moments.soft, node->moments.totalMass, node->moments.radius, node->moments.cm.x, node->moments.cm.y, node->moments.cm.z);
          }
#endif
              
          // now to add the node index to the list of interactions
          CkAssert(index >= 0);
#ifdef CHANGA_REFACTOR_PRINT_INTERACTIONS
          {

          if(b == TEST_BUCKET && tp->getIndex() == TEST_TP){
            CkPrintf("[%d]: (-1) -> (%d)\n", thisIndex, index);
          }
          }
#endif

          ILCell tilc(index, clist[i].offsetID);
          rrState->nodeLists.push_back(b, tilc, rrState, tp);
          if(rrState->nodeOffloadReady()){
            // enough node interactions to offload
            sendNodeInteractionsToGpu(rrState, tp);
            resetCudaNodeState(rrState);
          }
        }// length
      }// level
#ifdef CHANGA_REFACTOR_MEMCHECK
    CkPrintf("memcheck after node calc, bucket %d\n", b);
    CmiMemoryCheck();
#endif
      
    if(hasRemoteLists){
      for(int level = 0; level <= maxlevel; level++){
        // remote particles
        CkVec<RemotePartInfo> &rpilist = state->rplists[level];
        // for each bunch of particles in list
        for(int i = 0; i < rpilist.length(); i++){

          RemotePartInfo &rpi = rpilist[i];
          NodeKey key = rpi.key;
#ifdef CHANGA_REFACTOR_WALKCHECK_INTERLIST
          tp->addToBucketChecklist(b, key);
          tp->combineKeys(key, b);
#endif

          int gpuIndex = -1;
          // N.B. Keys for the particles are shifted to the left to
          // distinguish from the node.
          key <<= 1;
          std::map<NodeKey, int>::iterator q = cpref.find(key);
          if(q != cpref.end()){
            //cachedPartsOnGpu = true;
            gpuIndex = q->second;
          }
#ifdef CHANGA_REFACTOR_PRINT_INTERACTIONS
          {

          if(b == TEST_BUCKET && tp->getIndex() == TEST_TP){
            Vector3D<cosmoType> &vec = rpi.offset;
            CkPrintf("[%d]: remote: %d resume: %d bucket %d with remote part %ld (%d), (%1.0f,%1.0f,%1.0f)\n", thisIndex, getOptType() == Remote, state->resume, b, key, gpuIndex, vec.x, vec.y, vec.z);
          }
          }
#endif
          DoubleWalkState *rrState;
          if(state->resume || (!state->resume && gpuIndex < 0)){
            rrState = (DoubleWalkState*) tp->sInterListStateRemoteResume;
            CkAssert(rrState->particles);
            std::unordered_map<NodeKey,int>::iterator it = rrState->partMap.find(key);
            if(it == rrState->partMap.end()){
              gpuIndex = rrState->particles->length();
              rrState->partMap[key] = gpuIndex;
              for(int j = 0; j < rpi.numParticles; j++){
                rrState->particles->push_back(CompactPartData(rpi.particles[j]));
              }
            }
            else{
              gpuIndex = it->second;
            }
          }
          else{ // index >= 0
            rrState = state;
          }

          // now to add the node index to the list of interactions
          CkAssert(gpuIndex >= 0);

          // put bucket in interaction list
          Vector3D<cosmoType> off = rpi.offset;
          ILPart tilp (gpuIndex, encodeOffset(0, off.x, off.y, off.z), rpi.numParticles);
          rrState->particleLists.push_back(b, tilp, rrState, tp);
          if(rrState->partOffloadReady()){
            // enough nodes to offload
            sendPartInteractionsToGpu(rrState, tp);
            resetCudaPartState(rrState);
          }

        }// length
      }// level
    }// hasRemoteLists
#ifdef CHANGA_REFACTOR_MEMCHECK
    CkPrintf("memcheck after rp calc, bucket %d\n", b);
    CmiMemoryCheck();
#endif

      // local particles
      if(hasLocalLists){
        if(tp->largePhase()){
          for(int level = 0; level <= maxlevel; level++){
            CkVec<LocalPartInfo> &lpilist = state->lplists[level];
            for(int i = 0; i < lpilist.length(); i++){
              LocalPartInfo &lpi = lpilist[i];
#if defined CHANGA_REFACTOR_WALKCHECK_INTERLIST || defined CHANGA_REFACTOR_PRINT_INTERACTIONS
              NodeKey key = lpi.key;
#endif
#if defined CHANGA_REFACTOR_WALKCHECK_INTERLIST
              tp->addToBucketChecklist(b, key);
              tp->combineKeys(key, b);
#endif

              int gpuIndex = lpi.nd->bucketArrayIndex;

#ifdef CHANGA_REFACTOR_PRINT_INTERACTIONS
              {
                if(b == TEST_BUCKET && tp->getIndex() == TEST_TP){
                  Vector3D<cosmoType> &vec = lpi.offset;
                  CkPrintf("[%d]: remote: %d resume: %d bucket %d with local part %ld (%d), (%1.0f,%1.0f,%1.0f)\n", thisIndex, getOptType() == Remote, state->resume, b, key, gpuIndex, vec.x, vec.y, vec.z);
                }
              }
#endif
          
              if(gpuIndex >= 0){

              // put bucket in interaction list
                Vector3D<cosmoType> &off = lpi.offset;
                ILPart tilp(gpuIndex, encodeOffset(0, off.x, off.y, off.z), lpi.numParticles);
                state->particleLists.push_back(b, tilp, state, tp);
                if(state->partOffloadReady()){
                  // enough nodes to offload
                  sendPartInteractionsToGpu(state, tp);
                  resetCudaPartState(state);
                }
              }
            }
          }
        }
        else{ // small phase, need to attach particle data to state->particles
          for(int level = 0; level <= maxlevel; level++){
            CkVec<LocalPartInfo> &lpilist = state->lplists[level];
            for(int i = 0; i < lpilist.length(); i++){
              LocalPartInfo &lpi = lpilist[i];
#if defined CHANGA_REFACTOR_WALKCHECK_INTERLIST || defined CHANGA_REFACTOR_PRINT_INTERACTIONS
              NodeKey key = lpi.key;
#endif
#if defined CHANGA_REFACTOR_WALKCHECK_INTERLIST
              tp->addToBucketChecklist(b, key);
              tp->combineKeys(key, b);
#endif

              int gpuIndex = lpi.nd->bucketArrayIndex;
              if(gpuIndex < 0){
                CkAssert(state->particles != NULL);
                gpuIndex = state->particles->length();
                lpi.nd->bucketArrayIndex = gpuIndex;
                state->markedBuckets.push_back(lpi.nd);
                for(int j = 0; j < lpi.numParticles; j++){
                  state->particles->push_back(CompactPartData(lpi.particles[j]));
                }
              }

#ifdef CHANGA_REFACTOR_PRINT_INTERACTIONS
              {
                if(b == TEST_BUCKET && tp->getIndex() == TEST_TP){
                  Vector3D<cosmoType> &vec = lpi.offset;
                  CkPrintf("[%d]: remote: %d resume: %d bucket %d with local part %ld (%d), (%1.0f,%1.0f,%1.0f)\n", thisIndex, getOptType() == Remote, state->resume, b, key, gpuIndex, vec.x, vec.y, vec.z);
                }
              }
#endif
              CkAssert(gpuIndex >= 0);

              // put bucket in interaction list
              Vector3D<cosmoType> &off = lpi.offset;
              ILPart tilp(gpuIndex, encodeOffset(0, off.x, off.y, off.z), lpi.numParticles);
              state->particleLists.push_back(b, tilp, state, tp);
              if(state->partOffloadReady()){
                // enough nodes to offload
                sendPartInteractionsToGpu(state, tp);
                resetCudaPartState(state);
              }
            }
          }
        }
      }// level
#ifdef CHANGA_REFACTOR_MEMCHECK
    CkPrintf("memcheck after lp calc, bucket %d\n", b);
    CmiMemoryCheck();
#endif

    }// active
  }// bucket
  // *********************
#endif // ifndef CUDA
}

/// @brief CkLoop version of stateReady()
///
/// This function is called in the ckloop so multiple threads are concurrently
/// updating forces on its set of particles. Do not modify state variables of
/// TreePiece without locking.
void ListCompute::stateReadyPar(TreePiece *tp, int start, int end,
    CkVec<OffsetNode>& clist, CkVec<RemotePartInfo>& rpilist,
    CkVec<LocalPartInfo>& lpilist) {

#ifndef CUDA
  GravityParticle *particles = tp->getParticles();

  bool hasRemoteLists = rpilist.length() > 0 ? true : false;
  bool hasLocalLists = lpilist.length() > 0 ? true : false;

  for(int b = start; b < end; b++){
    if(tp->bucketList[b]->rungs >= activeRung){
      calcNodeForces(tp, b, activeRung, clist);

      // remote particles
      if(hasRemoteLists){
        calcParticleForces(tp, b, activeRung, rpilist);
      }

      // local particles
      if(hasLocalLists){
        calcParticleForces(tp, b, activeRung, lpilist);
      }
    }// active
  }// bucket
#else
  CkAbort("bUseCkLoopPar option does not work with CUDA. Turn bUseCkLoopPar off!\n");
#endif
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
	state->resume = resume;

        state->nodeLists.init(numBuckets, NUM_INIT_MOMENT_INTERACTIONS_PER_BUCKET);
        state->particleLists.init(numBuckets, NUM_INIT_PARTICLE_INTERACTIONS_PER_BUCKET);

}

/// @brief Reset node array after interactions have been sent to the GPU.
void ListCompute::resetCudaNodeState(DoubleWalkState *state){
  GenericTreeNode *tmp;
  state->nodeLists.reset();
  if(state->nodes){
    state->nodes->length() = 0;
    state->nodeMap.clear();
  }
  
}

void ListCompute::resetCudaPartState(DoubleWalkState *state){
  state->particleLists.reset();
  if(state->particles){
    state->particles->length() = 0;
    state->partMap.clear();
  }
}

void cudaCallback(void *param, void *msg){
  CudaRequest *data = (CudaRequest *)param;
  // Paranoid about data corruption here.
  CkAssert(((CkCallback *)data->cb)->type == CkCallback::callbackType::callCFn);
  CkAssert(((CkCallback *)data->cb)->d.cfn.fn == cudaCallback);
  CkAssert(data->numBucketsPlusOne > 0);
  CkAssert(data->numInteractions >= 0);

  int *affectedBuckets = data->affectedBuckets;
  TreePiece *tp = (TreePiece*)data->tp;
  DoubleWalkState *state = (DoubleWalkState *)data->state;
  int bucket;

  int numBucketsDone = data->numBucketsPlusOne-1;
  
#if COSMO_PRINT_BK > 1
  CkPrintf("[%d] CUDACALLBACK(node/part: %d, remote: %d, resume: %d), numBucketsDone: %d\n", 
                                  tp->getIndex(),
                                  data->node, data->remote, state->resume,
                                  numBucketsDone);
#endif
  // bucket computations finished
  //
#ifdef CHANGA_REFACTOR_MEMCHECK
  CkPrintf("memcheck before cudaCallback\n");
  CmiMemoryCheck();
#endif
  for(int i = 0; i < numBucketsDone; i++){
    bucket = affectedBuckets[i];
    state->counterArrays[0][bucket]--;
#if COSMO_PRINT_BK > 1
    CkPrintf("[%d] bucket %d numAddReq: %d,%d\n", tp->getIndex(), bucket, tp->getSRemoteGravityState()->counterArrays[0][bucket], tp->getSLocalGravityState()->counterArrays[0][bucket]);
#endif
    tp->finishBucket(bucket);
  }

#ifdef CHANGA_REFACTOR_MEMCHECK
  CkPrintf("memcheck before cudaCallback freeing\n");
  CmiMemoryCheck();
#endif

  // free data structures 
  if(numBucketsDone > 0){
    delete [] data->affectedBuckets;
  }
  delete ((CkCallback *)data->cb);
  delete data; 
#ifdef CHANGA_REFACTOR_MEMCHECK
  CkPrintf("memcheck after cudaCallback freeing\n");
  CmiMemoryCheck();
#endif
}

void ListCompute::sendNodeInteractionsToGpu(DoubleWalkState *state,
                                            TreePiece *tp){

  int thisIndex = tp->getIndex();

#ifdef CHANGA_REFACTOR_MEMCHECK
  CkPrintf("memcheck before sendNodeInteractionsToGpu\n");
  CmiMemoryCheck();
#endif

  CudaRequest *data = state->nodeLists.serialize(tp);

  data->state = (void *)state;
  data->node = true;
  data->remote = (getOptType() == Remote);

#ifdef HAPI_INSTRUMENT_WRS
  data->tpIndex = tp->getInstrumentId();
  data->phase = tp->getActiveRung();
#endif

#ifdef CUDA_PRINT_TRANSFERRED_INTERACTIONS
  CkPrintf("*************\n");
  CkPrintf("[%d] %d cellLists:\n", thisIndex, getOptType());
  CkPrintf("*************\n");

  ILCell *cellLists = (ILCell *)data->list;

  for(int i = 0; i < data->numInteractions; i++){
    Vector3D<cosmoType> v = tp->decodeOffset(cellLists[i].offsetID);
    CkPrintf("[%d] %d (%1.0f,%1.0f,%1.0f)\n", thisIndex, cellLists[i].index, v.x, v.y, v.z);
  }
  CkPrintf("\nnodeInteractionBucketMarkers:\n");
  for(int i = 0; i < data->numBucketsPlusOne; i++){
    CkPrintf("[%d] %d\n", thisIndex, data->bucketMarkers[i]);
  }

  CkPrintf("\nbucketSizes:\n");
  for(int i = 0; i < data->numBucketsPlusOne-1; i++){
    CkPrintf("[%d] %d\n", thisIndex, data->bucketSizes[i]);
  }
  CkPrintf("\nbucketStartMarkersNodes:\n");
  for(int i = 0; i < data->numBucketsPlusOne-1; i++){
    CkPrintf("[%d] %d\n", thisIndex, data->bucketStarts[i]);
  }
  CkPrintf("*************\n");
#endif



#if COSMO_PRINT_BK > 1
  CkPrintf("[%d] OFFLOAD (remote: %d, node: %d, resume: %d) nodeInteractions: %d, numBuckets: %d\n", thisIndex, data->remote, true, state->resume, data->numInteractions, data->numBucketsPlusOne-1);
#endif

  OptType type = getOptType();
  data->cb = new CkCallback(cudaCallback, data);

  if(data->numBucketsPlusOne == 1) { // Nothing for the GPU to do
      cudaCallback(data, NULL);     // Clean up
      return;
  }

#ifdef HAPI_INSTRUMENT_WRS
  double time = state->nodeListConstructionTimeStop();
#endif
  if(type == Local){
#ifdef CUDA_STATS
    tp->localNodeInteractions += state->nodeLists.totalNumInteractions;
#endif
    TreePieceCellListDataTransferLocal(data);
#ifdef HAPI_INSTRUMENT_WRS
    tp->localNodeListConstructionTime += time;
    tp->nLocalNodeReqs++;
#endif
  }
  else if(type == Remote && !state->resume){
#ifdef CUDA_STATS
    tp->remoteNodeInteractions += state->nodeLists.totalNumInteractions;
#endif
    TreePieceCellListDataTransferRemote(data);
#ifdef HAPI_INSTRUMENT_WRS
    tp->remoteNodeListConstructionTime += time;
    tp->nRemoteNodeReqs++;
#endif
  }
  else if(type == Remote && state->resume){
    CudaMultipoleMoments *missedNodes = state->nodes->getVec();
    int len = state->nodes->length();
    CkAssert(missedNodes);
#ifdef CUDA_STATS
    tp->remoteResumeNodeInteractions += state->nodeLists.totalNumInteractions;
#endif
    TreePieceCellListDataTransferRemoteResume(data, missedNodes, len);
#ifdef HAPI_INSTRUMENT_WRS
    tp->remoteResumeNodeListConstructionTime += time;
    tp->nRemoteResumeNodeReqs++;
#endif
  }
  else {CkAssert(0);}
#ifdef CHANGA_REFACTOR_MEMCHECK
  CkPrintf("memcheck after sendNodeInteractionsToGpu\n");
  CmiMemoryCheck();
#endif
#ifdef HAPI_INSTRUMENT_WRS
  state->nodeListConstructionTimeStart();
#endif
}

void ListCompute::sendPartInteractionsToGpu(DoubleWalkState *state,
                                            TreePiece *tp){

  int thisIndex = tp->getIndex();

#ifdef CHANGA_REFACTOR_MEMCHECK
  CkPrintf("memcheck before sendPartInteractionsToGpu\n");
  CmiMemoryCheck();
#endif

  CudaRequest *data = state->particleLists.serialize(tp);

  data->state = (void *)state;
  data->node = false;
  data->remote = (getOptType() == Remote);

#ifdef HAPI_INSTRUMENT_WRS
  data->tpIndex = tp->getInstrumentId();
  data->phase = tp->getActiveRung();
#endif

#ifdef CUDA_PRINT_TRANSFERRED_INTERACTIONS
  CkPrintf("*************\n");
  CkPrintf("[%d] partLists: ", thisIndex);
  CkPrintf("*************\n");

  ILPart *partLists = (ILPart *)data->list;

  for(int i = 0; i < data->numInteractions; i++){
    Vector3D<cosmoType> v = tp->decodeOffset(partLists[i].off);
    CkPrintf("[%d] %d (%1.0f,%1.0f,%1.0f)\n", thisIndex, partLists[i].index, v.x, v.y, v.z);
  }
  CkPrintf("\npartInteractionBucketMarkers:\n");
  for(int i = 0; i < data->numBucketsPlusOne; i++){
    CkPrintf("[%d] %d\n", thisIndex, data->bucketMarkers[i]);
  }

  CkPrintf("\nbucketSizes:\n");
  for(int i = 0; i < data->numBucketsPlusOne-1; i++){
    CkPrintf("[%d] %d\n", thisIndex, data->bucketSizes[i]);
  }
  CkPrintf("\nbucketStartMarkersParts:\n");
  for(int i = 0; i < data->numBucketsPlusOne-1; i++){
    CkPrintf("[%d] %d\n", thisIndex, data->bucketStarts[i]);
  }
  CkPrintf("*************\n");
#endif


#if COSMO_PRINT_BK > 1
  CkPrintf("[%d] OFFLOAD (remote: %d, node: %d, resume: %d) partInteractions: %d, numBuckets: %d\n", thisIndex, data->remote, false, state->resume, data->numInteractions, data->numBucketsPlusOne-1);
#endif

  OptType type = getOptType();
  data->cb = new CkCallback(cudaCallback, data);

  if(data->numBucketsPlusOne == 1) { // Nothing for the GPU to do
      cudaCallback(data, NULL);     // Clean up
      return;
  }

#ifdef HAPI_INSTRUMENT_WRS
  double time = state->partListConstructionTimeStop();
#endif

  if(type == Local){
#ifdef CUDA_STATS
    tp->localPartInteractions += state->particleLists.totalNumInteractions;
#endif
    if(tp->largePhase()){
      TreePiecePartListDataTransferLocal(data);
    }
    else{
      CompactPartData *parts = state->particles->getVec();
      int leng = state->particles->length();
      TreePiecePartListDataTransferLocalSmallPhase(data, parts, leng);
      tp->clearMarkedBuckets(state->markedBuckets);
    }
#ifdef HAPI_INSTRUMENT_WRS
    tp->localPartListConstructionTime += time;
    tp->nLocalPartReqs++;
#endif
  }
  else if(type == Remote && !state->resume){
#ifdef CUDA_STATS
    tp->remotePartInteractions += state->particleLists.totalNumInteractions;
#endif
    TreePiecePartListDataTransferRemote(data);
#ifdef HAPI_INSTRUMENT_WRS
    tp->remotePartListConstructionTime += time;
    tp->nRemotePartReqs++;
#endif
  }
  else if(type == Remote && state->resume){
    CompactPartData *missedParts = state->particles->getVec();
    int len = state->particles->length();
    CkAssert(missedParts);
#ifdef CUDA_STATS
    tp->remoteResumePartInteractions += state->particleLists.totalNumInteractions;
#endif
    TreePiecePartListDataTransferRemoteResume(data, missedParts, len);
#ifdef HAPI_INSTRUMENT_WRS
    tp->remoteResumePartListConstructionTime += time;
    tp->nRemoteResumePartReqs++;
#endif
  }
  else {CkAssert(0);}
#ifdef CHANGA_REFACTOR_MEMCHECK
  CkPrintf("memcheck after sendPartInteractionsToGpu\n");
  CmiMemoryCheck();
#endif
#ifdef HAPI_INSTRUMENT_WRS
  state->partListConstructionTimeStart();
#endif
}

#endif

void ListCompute::addChildrenToCheckList(GenericTreeNode *node, int reqID, int chunk, int awi, State *s, CheckList &chklist, TreePiece *tp){

  Vector3D<cosmoType> vec = tp->decodeOffset(reqID);
#ifdef CHANGA_REFACTOR_INTERLIST_PRINT_LIST_STATE
  int bucket = decodeReqID(reqID);
#endif

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
#if COSMO_PRINT_BK > 1
        CkPrintf("[%d] missed node %ld (chunk %d)\n", tp->getIndex(), childKey, chunk);
#endif

        nodeMissedEvent(reqID, chunk, s, tp);
#ifdef CHANGA_REFACTOR_INTERLIST_PRINT_LIST_STATE
        if(bucket == TEST_BUCKET && tp->getIndex() == TEST_TP){
          CkPrintf("[%d] missed child %ld (%1.0f,%1.0f,%1.0f) of %ld\n", tp->getIndex(),
              childKey,
              vec.x, vec.y, vec.z,
              node->getKey());
        }
#endif
        continue;
      }
    }
#ifdef CHANGA_REFACTOR_INTERLIST_PRINT_LIST_STATE
    if(bucket == TEST_BUCKET && tp->getIndex() == TEST_TP){
      CkPrintf("[%d] enq avail child %ld (%1.0f,%1.0f,%1.0f) of %ld\n", tp->getIndex(),
          child->getKey(),
          vec.x, vec.y, vec.z,
          node->getKey());
    }
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
    Vector3D<cosmoType> vec;
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
    Vector3D<cosmoType> vec;
    for(int i = 0; i < list.length(); i++){
      vec = list[i].offset;
      CkPrintf("%ld(%1.0f,%1.0f,%1.0f)\n", list[i].key, vec.x, vec.y, vec.z);
    }
  }

  if(hasRemoteLists){
    CkPrintf("----------------\n");
    CkVec<RemotePartInfo> &list = state->rplists[level];
    Vector3D<cosmoType> vec;
    for(int i = 0; i < list.length()-1; i++){
      vec = list[i].offset;
      CkPrintf("%ld(%1.0f,%1.0f,%1.0f)\n", list[i].key, vec.x, vec.y, vec.z);
    }
  }
}

void printUndlist(DoubleWalkState *state, int level, TreePiece *tp){
    CkPrintf("----------------\n");
    CkVec<OffsetNode> &list = state->undlists[level];
    Vector3D<cosmoType> vec;
    for(int i = 0; i < list.length()-1; i++){
      vec = tp->decodeOffset(list[i].offsetID);
      CkPrintf("%ld(%1.0f,%1.0f,%1.0f)\n", list[i].node->getKey(), vec.x, vec.y, vec.z);
    }
}
#endif

#endif // INTERLIST_VER > 0

void RemoteTreeBuilder::registerNode(GenericTreeNode *node){
  tp->nodeLookupTable[node->getKey()] = node;
}

void LocalTreeBuilder::registerNode(GenericTreeNode *node){
  tp->nodeLookupTable[node->getKey()] = node;
}

bool RemoteTreeBuilder::work(GenericTreeNode *node, int level){
  CkAssert(node != NULL);
  CkAssert(node->isValid() && !node->isCached());

  Tree::NodeType type = node->getType();

  switch(type){
    case NonLocalBucket:
    case NonLocal:
    {
      // find a remote index for the node
      int first, last;
      bool isShared = tp->nodeOwnership(node->getKey(), first, last);
      if (last < first) {
        // the node is really empty because falling between two TreePieces
        node->makeEmpty();
        node->remoteIndex = tp->thisIndex;
      } 
      else{
        // Choose a piece from among the owners from which to
        // request moments in such a way that if I am a piece with a
        // higher index, I request from a higher indexed treepiece.
        node->remoteIndex = tp->getResponsibleIndex(first,last);
        // request the remote chare to fill this node with the Moments
        if(requestNonLocalMoments){
          CkEntryOptions opts;
          opts.setPriority((unsigned int) -110000000);
          streamingProxy[node->remoteIndex].requestRemoteMoments(node->getKey(), tp->thisIndex, &opts);
        }
      }
      registerNode(node);
      return false;
    }

    case Empty:
    case Internal:
      return false;

    case Boundary:
      node->makeOctChildren(tp->myParticles,tp->myNumParticles,level,
			    tp->pTreeNodes);
      // The boundingBox was used above to determine the spatially equal
      // split between the children.  Now reset it so it can be calculated
      // from the particle positions.
      node->boundingBox.reset();

      registerNode(node);
      return true;

    default:
      CkAbort("Bad node type in RemoteTreeBuilder\n");
      return false;
  }
}

void RemoteTreeBuilder::doneChildren(GenericTreeNode *node, int level){
  CkAssert(node->getType() == Boundary ||
           node->getType() == Internal);

  if(node->getType() == Boundary){
    for(int i = 0; i < node->numChildren(); i++){
      GenericTreeNode *child = node->getChildren(i);
      if(child->getType() == NonLocal ||
         child->getType() == Boundary) node->remoteIndex--;
    }
  }
}

bool LocalTreeBuilder::work(GenericTreeNode *node, int level){
  if (level == NodeKeyBits-2) {
    ckerr << tp->thisIndex << ": TreePiece: This piece of tree has exhausted all the bits in the keys.  Super double-plus ungood!" << endl;
    ckerr << "Left particle: " << (node->firstParticle) << " Right particle: " << (node->lastParticle) << endl;
    ckerr << "Left key : " << keyBits((tp->myParticles[node->firstParticle]).key, KeyBits).c_str() << endl;
    ckerr << "Right key: " << keyBits((tp->myParticles[node->lastParticle]).key, KeyBits).c_str() << endl;
    ckerr << "Node type: " << node->getType() << endl;
    ckerr << "myNumParticles: " << tp->myNumParticles << endl;
    CkAbort("Tree is too deep!");
    return false;
  }

#if INTERLIST_VER > 0
  node->startBucket = tp->numBuckets; 
#endif

  if(node->getType() == NonLocal || node->getType() == NonLocalBucket){
    // don't deliver moments of this node to clients
    // because:
    // 1. either this is a NonLocal node which my TP 
    // requested from a remote source, in which
    // case we couldn't have received the response
    // yet, since the local tree build is still on,
    // and we haven't yielded the processor since
    // the remote requests were sent out. That is,
    // the moments of this node aren't ready yet.

    // 2. or this is a NonLocal node whose moments
    // can be supplied by a TreePiece on this PE.
    // In this case, my TP cannot be the supplier
    // of this node to any other TP. That is, my
    // TP doesn't have to deliver this node to anyone
    return false;
  }
  else if(node->getType() == Internal &&
          (node->lastParticle - node->firstParticle < maxBucketSize || level >= NodeKeyBits-3)){

       if(level >= NodeKeyBits-3
	  && node->lastParticle - node->firstParticle >= maxBucketSize)
           ckerr << "Truncated tree with "
                 << node->lastParticle - node->firstParticle
                 << " particle bucket" << endl;
       
      CkAssert(node->firstParticle != 0 && node->lastParticle != tp->myNumParticles+1);
      node->remoteIndex = tp->thisIndex;
      node->makeBucket(tp->myParticles);
      tp->bucketList.push_back(node);
      tp->numBuckets++;

      registerNode(node);
      // deliver moments, since doneChildren() will
      // never be called with a bucket
      tp->deliverMomentsToClients(node);
      return false;
  }
  else if(node->getType() == Empty){
    node->remoteIndex = tp->thisIndex;

    registerNode(node);
    // deliver MomentsToClients, since remote pieces don't know its
    // an empty node.
    tp->deliverMomentsToClients(node);
    return false;
  }
  else if(node->getType() == Internal){
    node->makeOctChildren(tp->myParticles,tp->myNumParticles,level,
	                  tp->pTreeNodes);
    // The boundingBox was used above to determine the spacially equal
    // split between the children.  Now reset it so it can be calculated
    // from the particle positions.
    node->boundingBox.reset();
    node->remoteIndex = tp->thisIndex;

    registerNode(node);
    // don't deliver MomentsToClients, since we
    // haven't yet computed the moments for this 
    // node: this will happen only after the moments
    // of both children of this node have been
    // completed (i.e. childrenDone() is called with
    // this node.
    return true;
  }
  else if(node->getType() == Boundary){
    CkAssert(node->getChildren(0) != NULL);
    // don't deliver MomentsToClients yet: same 
    // explanation as above
    return true;
  }
  CkAbort("Unhandled node type");
  return false;
}

void LocalTreeBuilder::doneChildren(GenericTreeNode *node, int level){
  CkAssert(node->getType() == Boundary ||
           node->getType() == Internal);
  node->rungs = 0;
#if INTERLIST_VER > 0
  node->numBucketsBeneath = 0;
#endif

  if(node->getType() == Boundary){
    for(int i = 0; i < node->numChildren(); i++){
      GenericTreeNode *child = node->getChildren(i);
      if((child->getType() != NonLocal && child->getType() != NonLocalBucket
	  && child->getType() != Empty)
	 && child->rungs > node->rungs)
	  node->rungs = child->rungs;
#if INTERLIST_VER > 0
      node->numBucketsBeneath += child->numBucketsBeneath;
#endif
    }

    if(node->remoteIndex == 0){
      tp->boundaryParentReady(node);
      // call deliver MomentsToClients, since the 
      // moments for this node have been computed
      // in boundaryParent Ready
      tp->deliverMomentsToClients(node);
    }
  }
  else{
    for(int i = 0; i < node->numChildren(); i++){
      GenericTreeNode *child = node->getChildren(i);
      if(child->getType() != Empty) {
	  tp->accumulateMomentsFromChild(node,child); 

	  if(child->rungs > node->rungs) node->rungs = child->rungs;
#if INTERLIST_VER > 0
	  node->numBucketsBeneath += child->numBucketsBeneath;
#endif
	  }
    }

    calculateRadiusFarthestCorner(node->moments, node->boundingBox);

    tp->deliverMomentsToClients(node);
  }
}

const char *typeString(NodeType type);
bool LocalTreePrinter::work(GenericTreeNode *node, int level){
  CkAssert(node != NULL);

#ifdef BIGKEYS
  uint64_t upper = node->getKey() >> 64;
  uint64_t lower = node->getKey();
  file << upper << lower
       << "[label=\"" << upper << lower
#else  
  file << node->getKey() 
       << "[label=\"" << node->getKey() 
#endif
       << " (" << typeString(node->getType()) << ")\\n"
       << node->remoteIndex
       << "\"]" << std::endl;

  if(node->getType() == Internal ||
     node->getType() == Boundary)
    return true;
  else 
    return false;
}

void LocalTreePrinter::doneChildren(GenericTreeNode *node, int level){
#ifdef BIGKEYS
  uint64_t upper = node->getKey() >> 64;
  uint64_t lower = node->getKey();
#endif
  for(int i = 0; i < node->numChildren(); i++){
    CkAssert(node->getChildren(i) != NULL);
#ifdef BIGKEYS
    uint64_t ch_upper = node->getChildren(i)->getKey() >> 64;
    uint64_t ch_lower = node->getChildren(i)->getKey();
    file << upper << lower << " -> " << ch_upper << ch_lower << std::endl;
#else
    file << node->getKey() << " -> " << node->getChildren(i)->getKey() << std::endl;
#endif
  }
}

void LocalTreePrinter::openFile(){
  std::ostringstream oss;
  oss << description << ".tree." << index << ".dot";
  file.open(oss.str().c_str());
  CkAssert(file.is_open());
  file << "digraph " << description << index << "{" << std::endl;
}
