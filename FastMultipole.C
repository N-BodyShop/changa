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

enum FMMOpenType {UNDECIDED, ACCEPT_PARTICLES, OPEN_BUCKET, OPEN_CHECK,
                  ACCEPT_MULTIPOLE, ACCEPT_MONO, ACCEPT_LOCAL_PARTICLES,
                  ACCEPT_SOFTENED_LOCAL_PARTICLES, ACCEPT_LOCAL,
                  ACCEPT_LOCAL_MONO, IGNORE};
        
class FMMCompute : public Compute {
  public:
  FMMCompute() : Compute(FMM) {}

  int doWork(GenericTreeNode *, TreeWalk *tw, State *state, int chunk, int reqID, bool isRoot, bool &didcomp, int awi);

  int openCriterion(TreePiece *ownerTP, GenericTreeNode *node, int reqID, State *state);

  // book keeping on notifications
  void nodeMissedEvent(int reqID, int chunk, State *state, TreePiece *tp);
  void nodeRecvdEvent(TreePiece *owner, int chunk, State *state, int bucket);
  void recvdParticles(ExternalGravityParticle *egp,int num,int chunk,int reqID,State *state, TreePiece *tp, Tree::NodeKey &remoteBucket);

  void initState(State *state);
  void stateReady(State *, TreePiece *, int chunk, int start, int end);
//  void printUndlist(DoubleWalkState *state, int level, TreePiece *tp);
//  void printClist(DoubleWalkState *state, int level, TreePiece *tp);
  void reassoc(void *cE, int activeRung, Opt *o);
  State *getNewState(int d1, int d2);
  State *getNewState(int d1);
  State *getNewState();
  void freeState(State *state);
  void freeDoubleWalkState(DoubleWalkState *state);

  private:

  void addChildrenToCheckList(GenericTreeNode *node, int reqID, int chunk, int awi, State *s, CheckList &chklist, TreePiece *tp);
  void addNodeToInt(GenericTreeNode *node, int offsetID, DoubleWalkState *s);
  void addLocalParticlesToInt(GravityParticle *parts, int n, Vector3D<double> &offset, DoubleWalkState *s);
  void addRemoteParticlesToInt(ExternalGravityParticle *parts, int n, Vector3D<double> &offset, DoubleWalkState *s);
  template <typename ParticleT>
  void addParticlesToChkList(ParticleT *part,
                               int nPart, int reqID,
                               CheckList &chklist, TreePiece *tp);
  DoubleWalkState *allocDoubleWalkState();
};

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
      if(myNode->isBucket())
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
    LOCR &L = s->momLocal;

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
       || open == OPEN_CHECK) {
        fakeOpen = 1;
        }
    else{
        fakeOpen = 0;
        }

    int action = opt->action(fakeOpen, node);
    switch(open) {
    case UNDECIDED:
        OffsetNode on;
        on.node = node;
        on.offsetID = reqID;
        undlist.push_back(on);                /// Is this the right list?
        return KEEP;
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
            addParticlesToChkList<GravityParticle>(part, computed, reqID,
                                                   chklist, tp);
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
            momLocrAddMomr5(&L,node->moments.mom,dir,-r.x,-r.y,-r.z,
                             &tax,&tay,&taz);
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
        node->makeBucket(&part[i]);
        OffsetNode on;
        on.node = node;
        on.offsetID = reqID;
        chklist.enq(on);
        }
    }

