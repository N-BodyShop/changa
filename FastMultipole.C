/*
 * FMM implementation based on Joachim's implementation in PKDGRAV2M
 *
 * Proposed strategy:
 * 
 * Start walk on the "owner" of the root cell.  For the first "shifts"
 * call a remote entry to start checking on another TreePiece.
 * 
 */

#include "ParallelGravity.h"
#include "GenericTreeNode.h"
#include "State.h"
#include "TreeWalk.h"
#include "Opt.h"
#include "moments.h"
#include "FastMultipole.h"

enum FMMOpenType {UNDECIDED, ACCEPT_PARTICLES, OPEN_BUCKET, OPEN_CHECK,
                  ACCEPT_MULTIPOLE, ACCEPT_MONO, ACCEPT_LOCAL_PARTICLES,
                  ACCEPT_SOFTENED_LOCAL_PARTICLES, ACCEPT_LOCAL,
                  ACCEPT_LOCAL_MONO, IGNORE};
        
/// @brief Calculate forces from local expansion on particles in a bucket
inline
int localBucketForce(LocalMoments &L,
		    Tree::GenericTreeNode *req,  
		    GravityParticle *particles, 
		    int activeRung)

{
    int computed = 0;
    Vector3D<cosmoType> r;

    for(int j = req->firstParticle; j <= req->lastParticle; ++j) {
        if (particles[j].rung >= activeRung) {
            particles[j].interMass += L.totalMass;
            computed++;
            r.x = particles[j].position.x - req->moments.cm.x;
            r.y = particles[j].position.y - req->moments.cm.y;
            r.z = particles[j].position.z - req->moments.cm.z;
            momEvalLocr(&L.mom, r.x, r.y, r.z, &particles[j].potential,
                        &particles[j].treeAcceleration.x,
                        &particles[j].treeAcceleration.y,
                        &particles[j].treeAcceleration.z);
            }
        }
    return computed;
}
/*
 * Summary of PKDGRAV2M walk:
 *
 * Updated to PKDGRAV2 version of Sept. 2013
 *
 * Given Check list (includes periodic replicas) and current cell:
 *
 * Process the check list:
 * 	For each cell on the list:
 *		if local expansion is valid and softening is valid
 *				Add cell moments to current local expansion
 *              else if check cell is larger
 *			if not a bucket
 *				add children to checklist
 *			else
 *                              add particles to checklist
 *              if cell fails far field opening or too few particles
 *                      if bucket add particles to p-p list
 *                      if cell add children to checklist
 *              else if not softened
 *                      Add to cell interaction list
 *              else if satisfies monopole opening
 *			Add to particle list
 *              else
 *                      if bucket add particles to p-p list
 *                      if cell add children to checklist
 *
 *       if current cell is not a bucket
 *		current cell is first child
 *		add sibling to checklist
 *		push particle list, cell list, check list, local
 * 							expansion on to stack
 *		shift local expansion
 *	 else
 *		calculate forces on all particles
 *		current cell is "next cell"
 *		pop stack for particle, cell, check, and local expansion
 */

/*
 * Below is an implementation of openCriterion() that checks all the
 * above conditions.
 * There are now several possible "returns" from an opencriterion
 * test.   Here is what doWork() should do based on these returns:
 *
 * Updated 9/12/13 to co-incide with 2013 version of pkdgrav2.  A big
 * difference is that there is no longer "open sink cell" conditions:
 * all cells on the check list are processed.
 *
 * UNDECIDED (0): source cell ("check") stays on the checklist
 * ACCEPT_PARTICLES (1): the source's particles are placed on the particle
 * interaction list.
 * OPEN_BUCKET (2): replace check cell with its particles;
 * 		continue processing check list.
 * OPEN_CHECK (3): replace check cell with children on check list;
 * 		continue processing check list.
 * ACCEPT_MULTIPOLE (4): accept as remote expansion for destination cell;
 * 	place on node interaction list.
 * ACCEPT_MONO (5): accept source cell as softened monopole for
 * 	destination cell, place on particle interaction list.
 * ACCEPT_LOCAL_PARTICLES (6): accept source, which are particles, as
 *      local expansion.  Can't happen.
 * ACCEPT_SOFTENED_LOCAL_PARTICLES (7): accept source, which are
 *      softened particles, as local expansion.  Can't happen.
 * ACCEPT_LOCAL (8): accept source (i.e. "check") cell as local expansion
 * 	for destination cell: use momLocrAddMomr5() to contribute to
 * 	local expansion.
 * ACCEPT_LOCAL_MONO (9): accept source (i.e. "check") particle as
 *      local expansion; Can't happen.
 * IGNORE (10): cell has no mass.
 */

// Below "myNode" is the node for which we are examining the check
// list, accumulating the local expansion and interaction terms.
// "node" is the node under consideration for acceptance.

FMMOpenType openCriterionFMM(Tree::GenericTreeNode *node, // "source" node
			 Tree::GenericTreeNode *myNode,	// "destination" node
                    Vector3D<double> offset,
                    int localIndex // requesting TreePiece
                    ) {
  // Note that some of this could be pre-calculated into an "opening radius"
  double radius = TreeStuff::opening_geometry_factor * node->moments.radius / theta;
  if(radius < node->moments.radius)
      radius = node->moments.radius;

  // Enhance my radius by a factor of 1.5 because of a far field -
  // local expansion assymetry.
  double myRadius = 1.5*TreeStuff::opening_geometry_factor * myNode->moments.radius / theta;
  if(myRadius < myNode->moments.radius)
      myRadius = myNode->moments.radius;

  Sphere<double> s(node->moments.cm + offset, radius);
  Sphere<double> myS(myNode->moments.cm, myRadius);

  if(!Space::intersect(s, myS)) {
      /// XXX check this
      if(!openSoftening(node, myNode, offset))
	  return ACCEPT_LOCAL;
      }
  else if(radius > myRadius) {	// Check Cell is larger
      if(node->isBucket()) {
          return OPEN_BUCKET;
          }
      else { // not a bucket
	  return OPEN_CHECK;
	  }
      }
  
  FMMOpenType iOpenA, iOpenB;
  if(node->isBucket())
      iOpenA = ACCEPT_PARTICLES;
  else
      iOpenA = OPEN_CHECK;
  
  // Always open node if this many particles or fewer.
  const int nMinParticleNode = 6;
  if(node->particleCount <= nMinParticleNode
     || Space::intersect(myNode->boundingBox, s)) {
      iOpenB = iOpenA;
      }
  else if(!openSoftening(node, myNode, offset)) {
      iOpenB = ACCEPT_MULTIPOLE;
      }
  else {
      double monoRadius = TreeStuff::opening_geometry_factor*node->moments.radius/thetaMono;
      Sphere<double> sM(node->moments.cm + offset, monoRadius);
      if(!Space::intersect(myNode->boundingBox, sM))
          iOpenB = ACCEPT_MONO;
      else {
          iOpenB = iOpenA;
          }
      }
      if(!myNode->isBucket())
          return UNDECIDED;
      else
          return iOpenB;
    }

int FMMCompute::openCriterion(TreePiece *ownerTP,
                          GenericTreeNode *node, int reqID, State *state){
  return openCriterionFMM(node,(GenericTreeNode *)computeEntity,
                          ownerTP->decodeOffset(reqID), ownerTP->getIndex());
}

/// @brief FMM version of ListCompute
int FMMCompute::doWork(GenericTreeNode *node,
			  TreeWalk *tw,
			  State *state,
			  int chunk,
			  int reqID,
			  bool isRoot, 
			  bool &didcomp, int awi)
{
    DoubleWalkState *s = (DoubleWalkState *)state;
    int level = s->level;
    CheckList &chklist = s->chklists[level];
    UndecidedList &undlist = s->undlists[level];
    LocalMoments &L = s->momLocal[level];

    TreePiece *tp = tw->getOwnerTP();
    Vector3D<double> offset = tp->decodeOffset(reqID);

    // so that we have a quick return in case of empty nodes
    if(node->getType() == Empty || node->getType() == CachedEmpty){
	return DUMP;
	}
    int open = openCriterion(tp, node, reqID, state);
    int fakeOpen;
    // in the FMM version, there are many possible return
    // values for the opencriterin function, whereas the Opt object
    // only knows about two (true, open and false, no need to open).
    // Catch the ones that involve opening.
    if(open == ACCEPT_PARTICLES || open == OPEN_BUCKET
       || open == OPEN_CHECK || open == UNDECIDED) {
        fakeOpen = 1;
        }
    else{
        fakeOpen = 0;
        }

    int action = opt->action(fakeOpen, node);
    switch(open) {
    case UNDECIDED:
        if(action == KEEP || action == KEEP_LOCAL_BUCKET
           || action == KEEP_REMOTE_BUCKET) {
            OffsetNode on;
            on.node = node;
            on.offsetID = reqID;
            undlist.push_back(on);                /// Is this the right list?
            return DUMP;
            }
        break;
    case ACCEPT_PARTICLES:
        if(action == KEEP_LOCAL_BUCKET){
            didcomp = true;
            // since this is a local bucket, we should have the
            // particles at hand
            GravityParticle *part = node->particlePointer;
            CkAssert(part);
            int computed = node->lastParticle-node->firstParticle+1;
            addLocalParticlesToInt(part, computed, offset, s);
            return DUMP;
            }
        else if(action == KEEP_REMOTE_BUCKET){
            didcomp = true;
            // fetch particles and compute.
            Tree::NodeKey keyref = node->getKey();
            ExternalGravityParticle *part;
            part = tp->particlesMissed(keyref,
                                       chunk,
                                       node->remoteIndex,
                                       node->firstParticle,
                                       node->lastParticle, reqID, false, awi,
                                       computeEntity);
            if(part){
                int computed = node->lastParticle-node->firstParticle+1;

                addRemoteParticlesToInt(part, computed, offset, s);
                }
            else{
                CkAssert(getOptType() == Remote);
                // particles missed
                int start, end;
                GenericTreeNode *source = (GenericTreeNode *)computeEntity;
                tp->getBucketsBeneathBounds(source, start, end);
                tp->updateUnfinishedBucketState(start, end, 1, chunk, state);
                }
            return DUMP;
            }
        break;
    case OPEN_BUCKET:
        if(action == KEEP_LOCAL_BUCKET){
            didcomp = true;
            // since this is a local bucket, we should have the
            // particles at hand
            GravityParticle *part = node->particlePointer;
            CkAssert(part);
            int computed = node->lastParticle-node->firstParticle+1;
            // XXX the following is dangerous
            addParticlesToChkList(part, computed,
                                  reqID, chklist, tp);
            return DUMP;
            }
        else if(action == KEEP_REMOTE_BUCKET){
            didcomp = true;
            // fetch particles and compute.
            Tree::NodeKey keyref = node->getKey();
            ExternalGravityParticle *part;
            part = tp->particlesMissed(keyref,
                                       chunk,
                                       node->remoteIndex,
                                       node->firstParticle,
                                       node->lastParticle, reqID, false, awi,
                                       computeEntity);
            if(part) {
                int computed = node->lastParticle-node->firstParticle+1;
                addParticlesToChkList<ExternalGravityParticle>(part, computed,
                                                               reqID, chklist,
                                                               tp);
                }
            else {
                CkAssert(getOptType() == Remote);
                // particles missed
                int start, end;
                GenericTreeNode *source = (GenericTreeNode *)computeEntity;
                tp->getBucketsBeneathBounds(source, start, end);
                tp->updateUnfinishedBucketState(start, end, 1, chunk, state);
                }
            return DUMP;
            }
    case OPEN_CHECK:
        if(action == KEEP) {
            addChildrenToCheckList(node, reqID, chunk, awi, s, chklist, tp);
            return DUMP;
            }
        break;
    case ACCEPT_MULTIPOLE:
    case ACCEPT_MONO:  /// The node bucket force detects whether this
                       /// is softend or not.
        if(action == COMPUTE) {
            didcomp = true;
            addNodeToInt(node, reqID, s);
            // all particles beneath this node have been
            // scheduled for computation
            return DUMP;
            }
        break;
    case ACCEPT_LOCAL:
        if(action == COMPUTE) {
            didcomp = true;
            MultipoleMoments &m = node->moments;
            Vector3D<cosmoType> cm(m.cm + offset);
            GenericTreeNode *myNode = (GenericTreeNode *)computeEntity;
            Vector3D<cosmoType> r = myNode->moments.cm - cm;
            cosmoType rsq = r.lengthSquared();
            cosmoType dir = COSMO_CONST(1.0)/sqrt(rsq);
            double tax, tay, taz;   // Used for gravstep.
            momLocrAddMomr5(&L.mom,&node->moments.mom,dir,r.x,r.y,r.z,
                             &tax,&tay,&taz);
            L.totalMass += node->moments.totalMass;
            return DUMP;
            }
        break;
    default:
        CkAbort("Bad openCriterion return");
        }
    if(action == DUMP || action == NOP){
        return DUMP;
        }
    CkAbort("ListCompute: bad walk state");
    return -1;
    }

/// @brief Turn an array of particles into cells and add them to the
/// check list.
/// @param part pointer to array of particles
/// @param nPart number of particles
/// @param reqID identifier for destination node
/// @param chklist Check list
/// @param tp TreePiece
template <typename ParticleT>
void FMMCompute::addParticlesToChkList(ParticleT *part,
                                         int nPart, int reqID,
                                         CheckList &chklist, TreePiece *tp) 
{
    for(int i = 0; i < nPart; i++) {
        // Create a bucket with one particle.
        GenericTreeNode *node = tp->pTreeNodes->alloc_one(0, Bucket, 0, 0, NULL);
        node->makeEmpty();
        // XXX This should be moved into a special creation interface
        // The following cast is dangerous, but we are only using the
        // External parts of Gravity particle.
        node->particlePointer = (GravityParticle *) part;
        node->moments += *part;
        node->boundingBox.grow(part->position);
        node->moments.radius = 0.0;
        
        OffsetNode on;
        on.node = node;
        on.offsetID = reqID;
        chklist.enq(on);
        }
    }

void FMMCompute::addRemoteParticlesToInt(ExternalGravityParticle *parts, int n, Vector3D<double> &offset, DoubleWalkState *s){
  RemotePartInfo rpi;
  int level = s->level;

  rpi.particles = parts;
  rpi.numParticles = n;
  rpi.offset = offset;

  s->rplists[level].push_back(rpi);
}

void FMMCompute::addLocalParticlesToInt(GravityParticle *parts, int n, Vector3D<double> &offset, DoubleWalkState *s){
  LocalPartInfo lpi;
  int level = s->level;

  lpi.particles = parts;
  lpi.numParticles = n;
  lpi.offset = offset;

  s->lplists[level].push_back(lpi);

}

void FMMCompute::addNodeToInt(GenericTreeNode *node, int offsetID, DoubleWalkState *s){
  OffsetNode nd;
  int level = s->level;
  nd.node = node;
  nd.offsetID = offsetID;

  s->clists[level].push_back(nd);
}

void FMMCompute::addChildrenToCheckList(GenericTreeNode *node, int reqID, int chunk, int awi, State *s, CheckList &chklist, TreePiece *tp){

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
        continue;
      }
    }
    // child available, enqueue
    OffsetNode on;
    on.node = child;
    on.offsetID = reqID;
    chklist.enq(on);
  }
}

void FMMCompute::stateReady(State *state_, TreePiece *tp, int chunk, int start, int end){
  int thisIndex = tp->getIndex();
  GravityParticle *particles = tp->getParticles();

  DoubleWalkState *state = (DoubleWalkState *)state_;

  GenericTreeNode *lowestNode = state->lowestNode;
  int maxlevel = lowestNode->getLevel(lowestNode->getKey());

  bool hasRemoteLists = state->rplists.length() > 0 ? true : false;
  bool hasLocalLists = state->lplists.length() > 0 ? true : false;

  for(int b = start; b < end; b++){
    if(tp->bucketList[b]->rungs >= activeRung){
      if(start + 1 == end)
          localBucketForce(state->momLocal[maxlevel], tp->getBucket(b),
                         particles, activeRung);
        
      for(int level = 0; level <= maxlevel; level++){
        CkVec<OffsetNode> &clist = state->clists[level];
        for(unsigned int i = 0; i < clist.length(); i++){

          int computed = 0;
          computed =  nodeBucketForce(clist[i].node,
              tp->getBucket(b),
              particles,
              tp->decodeOffset(clist[i].offsetID), activeRung);
          if(getOptType() == Remote){
            tp->addToNodeInterRemote(chunk, computed);
          }
          else if(getOptType() == Local){
            tp->addToNodeInterLocal(computed);
          }

        }

        // remote particles
        if(hasRemoteLists){
          CkVec<RemotePartInfo> &rpilist = state->rplists[level];
          // for each bunch of particles in list
          for(unsigned int i = 0; i < rpilist.length(); i++){
            RemotePartInfo &rpi = rpilist[i];

            // for each particle in a bunch
            for(int j = 0; j < rpi.numParticles; j++){
              int computed = 0;
              computed = partBucketForce(&rpi.particles[j],
                  tp->getBucket(b),
                  particles,
                  rpi.offset, activeRung);
              if(getOptType() == Remote){// don't really have to perform this check
                tp->addToParticleInterRemote(chunk, computed);
              }
            }

          }
        }

        // local particles
        if(hasLocalLists){
          CkVec<LocalPartInfo> &lpilist = state->lplists[level];
          for(unsigned int i = 0; i < lpilist.length(); i++){
            LocalPartInfo &lpi = lpilist[i];

            // for each particle in a bunch
            for(int j = 0; j < lpi.numParticles; j++){
              int computed = partBucketForce(&lpi.particles[j],
                  tp->getBucket(b),
                  particles,
                  lpi.offset, activeRung);
              tp->addToParticleInterLocal(computed);
            }

          }
        }
      }// level

    }// active
  }// bucket
}

void FMMCompute::reassoc(void *ce, int ar, Opt *o){
  computeEntity = ce;
  activeRung = ar;
  opt = o;
}

void FMMCompute::initState(State *state){

  DoubleWalkState *s = (DoubleWalkState *)state;
  int level = s->level;
  UndecidedList &myUndlist = s->undlists[level];

  // *my* undecided list:
  myUndlist.length() = 0;
  // interaction lists:
  s->clists[level].length() = 0;
  s->clists[level].reserve(1000);
  s->momLocal[level].totalMass = 0.0;
  momClearLocr(&s->momLocal[level].mom);
  if(getOptType() == Local){
    s->lplists[level].length() = 0;
    s->lplists[level].reserve(100);
  }
  else if(getOptType() == Remote){
    s->rplists[level].length() = 0;
    s->rplists[level].reserve(100);
  }
  else{
    CkAbort("Invalid Opt type for ListCompute");
  }
}

void FMMCompute::nodeMissedEvent(int reqID, int chunk, State *state, TreePiece *tp){
  CkAssert(getOptType() == Remote);
  int startBucket;
  int end;
  GenericTreeNode *source = (GenericTreeNode *)computeEntity;
  tp->getBucketsBeneathBounds(source, startBucket, end);
  tp->updateUnfinishedBucketState(startBucket, end, 1, chunk, state);
}

void FMMCompute::nodeRecvdEvent(TreePiece *owner, int chunk, State *state, int reqIDlist){
  int start, end;
  GenericTreeNode *source = (GenericTreeNode *)computeEntity;
  owner->getBucketsBeneathBounds(source, start, end);
  owner->updateBucketState(start, end, 1, chunk, state);

  CkAssert(chunk >= 0);
  int remainingChunk;
  remainingChunk = state->counterArrays[1][chunk];
  CkAssert(remainingChunk >= 0);
  if (remainingChunk == 0) {
    cacheGravPart[CkMyPe()].finishedChunk(chunk, owner->particleInterRemote[chunk]);
    owner->finishedChunk(chunk);
  }// end if finished with chunk
}

/// XXX here need to distinguish between interaction particles and
/// "multipole particles".
void FMMCompute::recvdParticles(ExternalGravityParticle *part,int num,int chunk,int reqID,State *state_,TreePiece *tp, Tree::NodeKey &remoteBucket){

  Vector3D<double> offset = tp->decodeOffset(reqID);
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
      state->momLocal[i].totalMass = 0.0;
      momClearLocr(&state->momLocal[i].mom);
      }

  if(remoteLists)
    for(int i = 0; i <= level; i++){
        state->rplists[i].length() = 0;
        }

  // put particles in list at correct level
  // (key) here.
  state->level = level;
  addRemoteParticlesToInt(part, num, offset, state);

  state->lowestNode = source;

  stateReady(state, tp, chunk, startBucket, end);

  tp->updateBucketState(startBucket, end, 1, chunk, state);

  int remainingChunk;
  remainingChunk = state->counterArrays[1][chunk];
  CkAssert(remainingChunk >= 0);
  if (remainingChunk == 0) {
    cacheGravPart[CkMyPe()].finishedChunk(chunk, tp->particleInterRemote[chunk]);
    tp->finishedChunk(chunk);
  }
}

State *FMMCompute::getNewState(){
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

void FMMCompute::freeState(State *s){
  freeDoubleWalkState((DoubleWalkState *)s);
  Compute::freeState(s);
}

void FMMCompute::freeDoubleWalkState(DoubleWalkState *state){
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
  if(state->momLocal.length() > 0) {
      state->momLocal.free();
      }
  
  if(state->placedRoots){
    delete [] state->placedRoots;
    state->placedRoots = 0;
  }
}


DoubleWalkState *FMMCompute::allocDoubleWalkState(){
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

  s->momLocal.resize(INTERLIST_LEVELS);
  return s;
}

State *FMMCompute::getNewState(int d1, int d2){
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

State *FMMCompute::getNewState(int d1){
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

