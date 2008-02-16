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

bool GravityCompute::openCriterion(TreePiece *ownerTP, 
                          GenericTreeNode *node, int reqID){
  return 
    openCriterionBucket(node,(GenericTreeNode *)computeEntity,ownerTP->decodeOffset(reqID), ownerTP->getIndex());
}

int GravityCompute::doWork(GenericTreeNode *node, TreeWalk *tw, 
                              State *state, int chunk, int reqID, bool isRoot, bool &didcomp){
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
  bool open = openCriterion(tp, node, reqID);
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
    
    //tp->addToNodeInterRemote(chunk, computed);
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
    ExternalGravityParticle *part;
    part = 
        tp->requestParticles(node->getKey(), 
                                       chunk, 
                                       node->remoteIndex, 
                                       node->firstParticle, 
                                       node->lastParticle, 
                                       reqID, false);
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
        tp->addToRemainingChunk(chunk, node->lastParticle-node->firstParticle+1);
      }
    }   
    return DUMP;
  }
  else if(action == NOP || action == DUMP){
    return DUMP;
  }
  // can do without DEFER if local walk begins from prefetchroots
  /*
  else if(action == DEFER){	//ask tw to check node ancestors 
    bool unopenedAncestors = false;
    if(isRoot){       // this is the first node in this walk: must check
      unopenedAncestors = tw->ancestorCheck(node, reqID);
    } 
    // unopenedAncestors will be set to true here if some ancestor of node
    // remained unopened and so, was or will be used in another computation
    if(unopenedAncestors){
      return DUMP;
    }
    else{       // no unopenedAncestors
      if(!open){ // no need to open node
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
      }// end if(!open)
      else{     // had to open node, check whether it is a bucket
        if(node->getType()==NonLocalBucket || node->getType()==CachedBucket){
#if CHANGA_REFACTOR_DEBUG > 2
          CkPrintf("[%d] GravityCompute told to DEFER, chunk=%d, remoteIndex=%d, first=%d, last=%d, reqID=%d\n", tp->getIndex(),
                         chunk, node->remoteIndex, 
                         node->firstParticle, 
                         node->lastParticle,
                         reqID);
#endif
          ExternalGravityParticle *part = 
                      tp->requestParticles(node->getKey(), 
                                                         chunk, 
                                                         node->remoteIndex, 
                                                         node->firstParticle, 
                                                         node->lastParticle, 
                                                         reqID, false);
                      didcomp = true;
          if(part){     // particles immediately available
            int computed = computeParticleForces(tp, node, part, reqID);
            //tp->addToParticleInterRemote(chunk, computed);
            if(getOptType() == Remote){
              tp->addToParticleInterRemote(chunk, computed);
            }
            else if(getOptType() == Local){
              tp->addToParticleInterLocal(computed);
            }

          }
          else if(getOptType() == Remote){
            tp->addToRemainingChunk(chunk, node->lastParticle-node->firstParticle+1);
          }
          return DUMP;
        }// end check whether node is a bucket
        else if(node->getType()==NonLocal || node->getType()==Cached){
          return KEEP;
        }// end check whether node is nonlocal or cached
        else{   // this shouldn't happen
          ckerr << "GravityCompute::doWork() : DEFERed node was of type "<< node->getType() << endl;
          CkAbort("");
        }// end error condition - wrong type of node DEFERed
      }// end node needs to be opened
    }// end there were no unopenedAncestors
  }// end action DEFERed
  */
}

// Source of force is a group of particles
int GravityCompute::computeParticleForces(TreePiece *ownerTP, GenericTreeNode *node, ExternalGravityParticle *part, int reqID){
  // jetley 
  /*
    for(int i = node->firstParticle; i <= node->lastParticle; i++){
      if(isnan(part[i-node->firstParticle].mass)){
        ostringstream oss;
        oss << "part[" << i-node->firstParticle << "].mass = nan" << endl;
        CkAbort(oss.str().data());
      }
    }
    */
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

/*
void PrefetchCompute::init(void *prs, int ar, Opt *o){
  computeEntity = prs;
  activeRung = ar;
  opt = o;
}
*/

// Called by TreeWalk when a node is requested on the Compute's behalf
// Also called by PrefetchCompute itself when it requests the prefetch of a
// bucket of particles

/*
void PrefetchCompute::nodeRequestSent(TreePiece *owner){
  owner->incPrefetchWaiting();
}
*/

bool PrefetchCompute::openCriterion(TreePiece *ownerTP, 
                          GenericTreeNode *node, int reqID){
  PrefetchRequestStruct * prs = (PrefetchRequestStruct *)computeEntity;
  Vector3D<double> offset = ownerTP->decodeOffset(reqID);
  for(int i = 0; i < prs->numPrefetchReq; i++){
    BinaryTreeNode testNode;
    testNode.boundingBox = prs->prefetchReq[i];
    //testNode.moments.softBound = 0.0;

    if(openCriterionBucket(node, &testNode, offset, ownerTP->getIndex()))
      return true;
  }
  return false;  
}

int PrefetchCompute::doWork(GenericTreeNode *node, TreeWalk *tw, State *state, int chunk, int reqID, bool isRoot, bool &didcomp){
  TreePiece *tp = tw->getOwnerTP();
  // ignores state 
  if(node == NULL){
    CkAbort("PrefetchCompute::doWork() given NULL node");
  }
  bool open = false;
  open = openCriterion(tp, node, reqID);

  int decision = opt->action(open, node);
  // the local node will become a CachedNonLocal[...] from a NonLocal[...]
  // tp->requestNode(node->remoteIndex, node->getKey(), chunk, reqID, true);
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
    ExternalGravityParticle *part = 
                      tp->requestParticles(node->getKey(), 
                                                         chunk, 
                                                         node->remoteIndex, 
                                                         node->firstParticle, 
                                                         node->lastParticle, 
                                                         reqID, true);
    
    if(part == NULL){
      tp->incPrefetchWaiting();
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


