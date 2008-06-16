#include "ParallelGravity.h"
#include "GenericTreeNode.h"
//#include "codes.h"

#include "Opt.h"
#include "Compute.h"
#include "TreeWalk.h"
//#include "State.h"
#include "Space.h"
#include "gravity.h"  // for openCriterion(..)

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

State *Compute::getResumeState(int bucketIdx){
  return 0;
}

State *Compute::getNewState(){
  return 0;
}

State *GravityCompute::getNewState(){
  State *s = new State();
  if(getOptType() == Remote){
    // 2 arrays of counters 
    // 0. numAdditionalRequests[] - sized numBuckets, init to numChunks
    // 1. remainingChunk[] - sized numChunks
    s->counterArrays.reserve(2);
  }
  else if(getOptType() == Local){
    // 0. local component of numAdditionalBuckets, init to 1
    s->counterArrays.reserve(1);
  }
  return s;
}

State *PrefetchCompute::getNewState(){
  State *s = new State();
  // prefetchWaiting
  s->counterArrays.reserve(1);
  s->counterArrays[0].reserve(1);
  // no per-bucket counters
  return s;
}

void GravityCompute::reassoc(void *ce, int ar, Opt *o){
  computeEntity = ce; 
  activeRung = ar;
  opt = o;
}

/*
void PrefetchCompute::init(void *buck, int ar, Opt *o){

  // buck is actually the TreePiece this PrefetchCompute
  // will query when it needs a prefetchRequestStruct
  computeEntity = buck;
  activeRung = ar;
  opt = o;
}
*/
int GravityCompute::nodeMissedEvent(int chunk, State *state){
  if(getOptType() == Remote){
    //owner->addToRemainingChunk(chunk, +1);
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

int ListCompute::nodeMissedEvent(int chunk, State *state){
  if(getOptType() == Remote){
    //owner->addToRemainingChunk(chunk, +1);
    state->counterArrays[1][chunk]++;
  }
}

bool GravityCompute::openCriterion(TreePiece *ownerTP, 
                          GenericTreeNode *node, int reqID, State *state){
  return 
    openCriterionBucket(node,(GenericTreeNode *)computeEntity,ownerTP->decodeOffset(reqID), ownerTP->getIndex());
}

extern int partBucketForce(ExternalGravityParticle *part, GenericTreeNode *req, GravityParticle *particles, Vector3D<double> offset, int activeRung);

void GravityCompute::recvdParticles(ExternalGravityParticle *part,int num,int chunk,int reqID,State *state,TreePiece *tp){
  //TreePiece *tp = tw->getOwnerTP();

  Vector3D<double> offset = tp->decodeOffset(reqID);
  int reqIDlist = decodeReqID(reqID);
  CkAssert(num > 0);
  //tp->bucketReqs[reqIDlist].numAdditionalRequests -= num;
  state->counterArrays[0][reqIDlist] -= num;
  //tp->remainingChunk[chunk] -= num;
  state->counterArrays[1][chunk] -= num;

  GenericTreeNode* reqnode = tp->bucketList[reqIDlist];

  int computed;
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

void PrefetchCompute::recvdParticles(ExternalGravityParticle *egp,int num,int chunk,int reqID,State *state, TreePiece *tp){
  // when we receive missed particles, we do the same book-keeping
  // as we did when we received a missed node.
  finishNodeProcessEvent(tp, state);
}

int GravityCompute::nodeRecvdEvent(TreePiece *owner, int chunk, State *state, int reqIDlist){
  //owner->bucketReqs[reqIDlist].numAdditionalRequests--;
  state->counterArrays[0][reqIDlist]--;
  owner->finishBucket(reqIDlist);

  // needn't perform check, because local gravity walks never enter this function at all
  if(getOptType() == Remote){ 
    CkAssert(chunk >= 0); 
    state->counterArrays[1][chunk]--;
    CkAssert(state->counterArrays[1][chunk] >= 0);
    if (state->counterArrays[1][chunk] == 0) {
      streamingCache[CkMyPe()].finishedChunk(chunk, owner->nodeInterRemote[chunk]+owner->particleInterRemote[chunk]);
      if(chunk == owner->numChunks-1) owner->markWalkDone();
    }// end if finished with chunk
  }// end if remote
}


/*
void PrefetchCompute::walkDone(){
  deleteComputeEntity();
}
*/

int GravityCompute::doWork(GenericTreeNode *node, TreeWalk *tw, 
                              State *state, int chunk, int reqID, bool isRoot, bool &didcomp, int awi){
  // ignores state
  
  TreePiece *tp = tw->getOwnerTP();
  // FIXME - better way to do this? 
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
  // check opening criterion
  bool open = openCriterion(tp, node, reqID, state);
  if(opt == NULL){       
    ckerr << "GravityCompute reqID("<<reqID<<"), isRoot("<<isRoot<<") has NULL opt" << endl;
    CkAbort("aborting");
  }
  int action = opt->action(open, node);
  if(action == KEEP){
    return KEEP;
  }
  else if(action == COMPUTE){
    didcomp = true;
    int computed = nodeBucketForce(node, 
                    (GenericTreeNode *)computeEntity, 
                    tp->getParticles(), 
                    tp->decodeOffset(reqID), 
                    activeRung);
    
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
      for(int i = node->firstParticle; i <= node->lastParticle; i++){
        computed += partBucketForce(
                                  &part[i-node->firstParticle],
                                  (GenericTreeNode *)computeEntity, 
                                  tp->getParticles(), 
                                  tp->decodeOffset(reqID), 
                                  activeRung);
      }
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
                                       reqID, false, awi);
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
        //tp->addToRemainingChunk(chunk, node->lastParticle-node->firstParticle+1);
        state->counterArrays[1][chunk] += node->lastParticle-node->firstParticle+1;
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
  for(int i = node->firstParticle; i <= node->lastParticle; i++){
    computed += partBucketForce(
                                  &part[i-node->firstParticle],
                                  (GenericTreeNode *)computeEntity, 
                                  ownerTP->getParticles(), 
                                  ownerTP->decodeOffset(reqID), 
                                  activeRung);
  }
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

bool PrefetchCompute::openCriterion(TreePiece *ownerTP, 
                          GenericTreeNode *node, int reqID, State *state){
  // redundant in this case, because we already have a pointer 
  // to the ownerTP.
  TreePiece *tp = (TreePiece *)computeEntity;
  PrefetchRequestStruct prs(tp->prefetchReq, tp->numPrefetchReq);
  Vector3D<double> offset = ownerTP->decodeOffset(reqID);
  for(int i = 0; i < prs.numPrefetchReq; i++){
    BinaryTreeNode testNode;
    testNode.boundingBox = prs.prefetchReq[i];
    //testNode.moments.softBound = 0.0;

    if(openCriterionBucket(node, &testNode, offset, ownerTP->getIndex()))
      return true;
  }
  return false;  
}

int PrefetchCompute::doWork(GenericTreeNode *node, TreeWalk *tw, State *state, int chunk, int reqID, bool isRoot, bool &didcomp, int awi){
  TreePiece *tp = tw->getOwnerTP();
  // ignores state 
  if(node == NULL){
    CkAbort("PrefetchCompute::doWork() given NULL node");
  }
  bool open = false;
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
                                                         true, awi);
    
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

/* List compute object
 * requires liststate object
 * on KEEP, descends further down tree
 * on DUMP, dumps subtree rooted at node in question
 * on COMPUTE, appends object to list
 * on KEEP_LOCAL_BUCKET, appends local bucket to list
 * on KEEP_REMOTE_BUCKET, appends remote bucket to list
 * *
 * */

int ListCompute::doWork(GenericTreeNode *node, TreeWalk *tw, State *state, int chunk, int reqID, bool isRoot, bool &didcomp, int awi){

  TreePiece *tp = tw->getOwnerTP();
  // FIXME - better way to do this? 
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
  // check opening criterion
  bool open = openCriterion(tp, node, reqID, state);
  if(opt == NULL){       
    ckerr << "ListCompute reqID("<<reqID<<"), isRoot("<<isRoot<<") has NULL opt" << endl;
    CkAbort("aborting");
  }
  int action = opt->action(open, node);
  if(action == KEEP){
    return KEEP;
  }
  else if(action == COMPUTE){
    didcomp = true;
    int computed;
    // FIXME - add node to node compute list, update computed
    
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
    CkPrintf("[%d] ListCompute told to KEEP_LOCAL_BUCKET, chunk=%d, remoteIndex=%d, first=%d, last=%d, reqID=%d\n", tp->getIndex(),
                         chunk, node->remoteIndex, 
                         node->firstParticle, 
                         node->lastParticle,
                         reqID);
#endif
      // since this is a local bucket, we should have the particles at hand
      GravityParticle *part = node->particlePointer;
      CkAssert(part);
      int computed = 0;
      for(int i = node->firstParticle; i <= node->lastParticle; i++){
        // FIXME - add particles to paricle list, update computed
      }
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
                                       reqID, false, awi);
    if(part){
#if CHANGA_REFACTOR_DEBUG > 2
      CkPrintf("Particles found in cache\n");
#endif
      // FIXME - add particles to particle list, update computed
      int computed;
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
        //tp->addToRemainingChunk(chunk, node->lastParticle-node->firstParticle+1);
        state->counterArrays[1][chunk] += node->lastParticle-node->firstParticle+1;
      }
    }   
    return DUMP;
  }
  else if(action == DUMP || action == NOP){
    return DUMP;
  }
}

bool ListCompute::openCriterion(TreePiece *ownerTP, 
                          GenericTreeNode *node, int reqID, State *state){
  return 
    openCriterionBucket(node,(GenericTreeNode *)computeEntity,ownerTP->decodeOffset(reqID), ownerTP->getIndex());
}

void ListCompute::addParticlesToList(ExternalGravityParticle *parts, int n, ListState *state){
  /*
  RemotePartInfo pinfo;
  pinfo.particles = parts;
  pinfo.numParticles = n;
  pinfo.offset = 
  state->remoteParticleList.push_back();
  */
}

void ListCompute::addParticlesToList(GravityParticle *, int n, ListState *state){
}
