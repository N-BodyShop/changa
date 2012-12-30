#include "config.h"
#include "GenericTreeNode.h"
#include "ParallelGravity.h"
#include "TreeWalk.h"
#include "State.h"
#include "Compute.h"
#include <stack>

#ifdef CHANGA_REFACTOR_WALKCHECK
#include <string>
#endif

const char *typeString(NodeType type);

void TreeWalk::init(Compute *c, TreePiece *owner){
  comp = c;
  ownerTP = owner;
}

void TreeWalk::reassoc(Compute *c){
  comp = c;
}

void TopDownTreeWalk::walk(GenericTreeNode *startNode, State *state, int chunk, int reqID, int awi){
#ifdef BENCHMARK_TIME_WALK
  double startTime = CmiWallTimer();
#endif
#ifndef CHANGA_REFACTOR_WALKCHECK
#ifdef TREE_BREADTH_FIRST
  bft(startNode, state, chunk, reqID, true, awi);       // isRoot
#else
  dft(startNode, state, chunk, reqID, true, awi);       // isRoot
#endif
#else
  bool doprint = false;
  if(comp->getSelfType() == Gravity && comp->getOptType() == Remote){
    // FIXME - TP::getCurrentRemote Bucket method is no more
    doprint = (ownerTP->/*getCurrentRemote Bucket()*/ == CHECK_BUCKET && ownerTP->getIndex() == CHECK_INDEX);
  }
  else if(comp->getSelfType() == Gravity && comp->getOptType() == Local){
    // FIXME - TP::getCurrent Bucket method is no more
    doprint = (ownerTP->/*getCurrent Bucket()*/ == CHECK_BUCKET && ownerTP->getIndex() == CHECK_INDEX);
  }
  if(doprint){
    CkPrintf("Printing walk (%d: type %c)\n", ownerTP->getIndex(), comp->getOptType() == Local ? 'L' : 'R');
  }
  dft(startNode, state, chunk, reqID, true, 0, doprint);       // isRoot
#endif
#ifdef BENCHMARK_TIME_WALK
  walkTime += CmiWallTimer() - startTime;
#endif
}

#ifndef CHANGA_REFACTOR_WALKCHECK
void TopDownTreeWalk::dft(GenericTreeNode *node, State *state, int chunk, int reqID, bool isRoot, int awi){
#else
void TopDownTreeWalk::dft(GenericTreeNode *node, State *state, int chunk, int reqID, bool isRoot, int shift, bool doprint){
#endif
  int ret;
  if(node == NULL){   // something went wrong here
    ckerr << "TopDownTreeWalk recvd. null node - chunk("<<chunk<<"), reqID("<<reqID<<"), isRoot("<<isRoot<<")" << endl;
    CkAbort("Abort");
  }
#ifdef BENCHMARK_TIME_WALK
  double start1 = CmiWallTimer();
#endif
  NodeKey globalKey = node->getKey();
  // process this node
  bool didcomp = false;
  ret = comp->doWork(node, this, state, chunk, reqID, isRoot, didcomp, awi);
#ifdef BENCHMARK_TIME_WALK
  doWorkTime += CmiWallTimer() - start1;
#endif

#ifdef CHANGA_REFACTOR_WALKCHECK
  if(doprint){
    string s;
    for(int i = 0; i < shift; i++) s = s + " ";
    char *arr = "KD";
    char *c = "NY";
    if(ret == KEEP)
      CkPrintf("%s%ld (%c)\n", s.c_str(),node->getKey(), arr[ret]);
    else if(ret == DUMP)
      CkPrintf("%s%ld (%c,%c)\n", s.c_str(),node->getKey(), arr[ret], c[didcomp]);

  }
#endif

#if CHANGA_REFACTOR_DEBUG > 2
  CkPrintf("[%d]: TopDownTreeWalk - chunk: %d, ComputeType: %d, isPrefetch: %d, reqID: %d\n", ownerTP->getIndex(), chunk, comp->getSelfType(), comp->getOptType() == Pref, reqID);
#endif
  if(ret == KEEP){      // descend further down tree
    for(int i = 0; i < node->numChildren(); i++){
      GenericTreeNode *child = node->getChildren(i);
      globalKey = node->getChildKey(i);

      comp->startNodeProcessEvent(state);

      // check whether child is NULL and get from cache if necessary/possible
      if(child == NULL){
        // needed to descend, but couldn't because node wasn't available
#if CHANGA_REFACTOR_DEBUG > 2
        CkPrintf("[%d]: TopDownTreeWalk missed a child: chunk=%d, TreeWalkType=%d, ComputeType=%d, OptType=%d, isPrefetch=%d\n",
            ownerTP->getIndex(),
            chunk,
            getSelfType(),
            comp->getSelfType(),
            comp->getOptType(),
            comp->getSelfType()==Prefetch);
#endif
#ifdef CHANGA_REFACTOR_WALKCHECK
        if(doprint){
          string s;
          for(int j = 0; j < shift+2; j++) s = s + " ";
          CkPrintf("%s%ld SM\n", s.c_str(), node->getChildKey(i));
        }
#endif
        child = ownerTP->nodeMissed(reqID,
                                    node->remoteIndex,
                                    globalKey,
                                    chunk, comp->getSelfType() == Prefetch,
                                    awi, comp->getComputeEntity());
        if(child == NULL){     // missed in cache, skip node for now
#if CHANGA_REFACTOR_DEBUG > 2
          CkPrintf("[%d]: child not found in cache\n", ownerTP->getIndex());
#endif
          comp->nodeMissedEvent(reqID, chunk, state, ownerTP);
#ifdef CHANGA_REFACTOR_WALKCHECK
          if(doprint){
            string s;
            for(int j = 0; j < shift+2; j++) s = s + " ";
            CkPrintf("%s%ld HM\n", s.c_str(), node->getChildKey(i));
          }
#endif
          continue;
        }
        else{
#if CHANGA_REFACTOR_DEBUG > 2
          CkPrintf("[%d]: child found in cache\n", ownerTP->getIndex());
#endif
        }
      }// end check NULL node

      // process children recursively
      // the next can't be the first node we are processing, so isRoot = false
#ifndef CHANGA_REFACTOR_WALKCHECK
      dft(child, state, chunk, reqID, false, awi);
#else
      dft(child, state, chunk, reqID, false, shift+2, doprint);
#endif

    }// for each child
  }// if KEEP
  //else // don't need the node anymore, return up the tree
  comp->finishNodeProcessEvent(ownerTP, state);
  return;
}

void TopDownTreeWalk::bft(GenericTreeNode *node, State *state, int chunk, int reqID, bool isRoot, int awi){
  int ret;
  if(node == NULL){   // something went wrong here
    ckerr << "TopDownTreeWalk recvd. null node - chunk("<<chunk<<"), reqID("<<reqID<<"), isRoot("<<isRoot<<")" << endl;
    CkAbort("Abort");
  }

  CkQ<GenericTreeNode*> queue(1024);
  queue.enq(node);

  while (node = queue.deq()) {

    NodeKey globalKey = node->getKey();
    // process this node
    bool didcomp = false;
    ret = comp->doWork(node, this, state, chunk, reqID, isRoot, didcomp, awi);

    if(ret == KEEP){      // descend further down tree
      for(unsigned int i = 0; i < node->numChildren(); i++){
        GenericTreeNode *child = node->getChildren(i);
        globalKey = node->getChildKey(i);

        comp->startNodeProcessEvent(state);

        // check whether child is NULL and get from cache if necessary/possible
        if(child == NULL){
          // needed to descend, but couldn't because node wasn't available
          child = ownerTP->nodeMissed(reqID, node->remoteIndex, globalKey, chunk, comp->getSelfType() == Prefetch, awi, comp->getComputeEntity());
          if(child == NULL){     // missed in cache, skip node for now
            comp->nodeMissedEvent(reqID, chunk, state, ownerTP);
            continue;
          }
          else{
          }
        }// end check NULL node

        // process children recursively
        // the next can't be the first node we are processing, so isRoot = false
        queue.enq(child);

      }// for each child
    }// if KEEP
    //else // don't need the node anymore, return up the tree
    comp->finishNodeProcessEvent(ownerTP, state);

  }
}


bool bIsReplica(int reqID);

///
/// Bottom up treewalk for efficient smooth:
/// If startNode is the root, and not a periodic,
/// then go down to the bucket being walked, pushing siblings onto a
/// stack in the local frame.
/// Once the stack is processed, then all walks are done in the
/// standard "top down" way.
///
void BottomUpTreeWalk::walk(GenericTreeNode *startNode, State *state,
			    int chunk, int reqID, int awi){
    int reqIDlist = decodeReqID(reqID);
    GenericTreeNode *reqnode = ownerTP->bucketList[reqIDlist];
    int ret;
    bool didcomp = false;
    bool isRoot = false;
    GenericTreeNode *node;

    if(reqnode->cpStart != NULL) return;
    
    node = startNode;
    if(bIsReplica(reqID) || node->getKey() != ownerTP->root->getKey()) {
	// Do standard top-down walk.
	if(node == NULL){   // something went wrong here
	    ckerr << "BottomUpTreeWalk recvd. null node - chunk("
		  <<chunk<<"), reqID("<<reqID<<"), isRoot("<<isRoot<<")"
		  << endl;
	    CkAbort("TreeWalk");
	    }
	NodeKey currentGlobalKey = node->getKey();
	// process this node
	ret = comp->doWork(node, this, state, chunk, reqID, isRoot, didcomp,
			   awi);
	if(ret == KEEP){      // descend further down tree
	    for(unsigned int i = 0; i < node->numChildren(); i++){
		GenericTreeNode *child = node->getChildren(i);
		currentGlobalKey = node->getChildKey(i);

		comp->startNodeProcessEvent(state);

		// check whether child is NULL and get from cache if
		// necessary/possible
		if(child == NULL){
		    // needed to descend, but couldn't because node
		    // wasn't available
		    child = ownerTP->nodeMissed(reqID, node->remoteIndex,
						currentGlobalKey, chunk,
						comp->getSelfType() == Prefetch,
						awi, comp->getComputeEntity());
		    if(child == NULL){   // missed in cache, skip node for now
			comp->nodeMissedEvent(reqID, chunk, state, ownerTP);
			continue;
			}
		    } // end check NULL node

		// process children recursively
		// the next can't be the first node we are processing,
		// so isRoot = false
		walk(child, state, chunk, reqID, awi);
		}// for each child
	    }// if KEEP
	// don't need the node anymore, return up the tree
	comp->finishNodeProcessEvent(ownerTP, state);
	return;
	}
    std::stack<GenericTreeNode *> nodeStack;
    std::stack<GenericTreeNode *> parentStack;
    // go down the tree toward the bucket, push all siblings of the
    // ancestor onto the stack.
    while(node != reqnode) {
	int which = node->whichChild(reqnode->getKey());
	for(int iChild = 0; iChild < node->numChildren(); iChild++) {
	    if(iChild != which)
		nodeStack.push(node->getChildren(iChild));
		parentStack.push(node);
	    }
	node = node->getChildren(which);
	}
    // work on the bucket
    ret = comp->doWork(node, this, state, chunk, reqID, isRoot, didcomp, awi);
    CkAssert(ret == DUMP);

    while(!nodeStack.empty()) {
	node = nodeStack.top();
	walk(node, state, chunk, reqID, awi);
	// Sphere containing search radii of all particles in bucket
	Sphere<cosmoType> s(reqnode->centerSm, reqnode->sizeSm + reqnode->fKeyMax);
	// XXX This test should be done inside the KNearestNeighbor
	// class for a cleaner interface.
	if(parentStack.top()->getType() == Internal
	   && Space::contains(parentStack.top()->boundingBox, s)) {
	    reqnode->cpStart = node->parent; // found containing node
	    reqnode->cpStart = parentStack.top();
	    CkAssert(node->parent->numChildren() <= 2);  // Assumes binary tree
	    if(verbosity > 3)
		CkPrintf("[%d] Bucket %d is contained in node %lx\n", ownerTP->thisIndex,
			 reqIDlist, reqnode->cpStart->getKey());
	    // clear stack
	    return;
	    }
	nodeStack.pop();
	parentStack.pop();
	}
}
/**
 * LocalTargetWalk functions.
 */

// This walk interprets what is otherwise the 'reqID' argument as the targetBucketIndex
// It can only be called from the calculateGravityRemote/startNextBucket
#if INTERLIST_VER > 0
void LocalTargetWalk::walk(GenericTreeNode *ancestor, State *state, int chunk, int targetBucketIndex, int awi){

    targetKey = ownerTP->getBucket(targetBucketIndex)->getKey();
    // construct lists
    int ancestorLevel = ancestor->getLevel(ancestor->getKey());
    dft(ancestor, state, chunk, targetBucketIndex, (ancestorLevel == 0), awi, ancestorLevel);
}

int reEncodeOffset(int reqID, int offsetID);

// This walk interprets what is otherwise the 'reqID' argument as the targetBucketIndex
void LocalTargetWalk::dft(GenericTreeNode *localNode, State *state, int chunk, int targetBucketIndex, bool isRoot, int awi, int level){

  bool descend = false;
  // localNode has changed, need to update computeEntity
  comp->setComputeEntity(localNode);

  // get current level's checklist
  DoubleWalkState *s = (DoubleWalkState *)state;
  s->level = level;
  CheckList &chklist = s->chklists[level];
  UndecidedList &myUndlist = s->undlists[level];


  if(!isRoot)
  {
    // we don't clear state for this level in case
    // this is the local root, because its chklist
    // will have been set just before walk was called.

    // this clears the undlist and interaction lists
    // for the current level
    comp->initState(s);

    UndecidedList &parentUndlist = s->undlists[level-1];
    // process parent's undecided nodes first
    // don't modify parentundlist - sibling will need the same list
    CkAssert(chklist.isEmpty());
    for(unsigned int i = 0; i < parentUndlist.length(); i++)
    {
      // get remote node from state
      OffsetNode &glblNode = parentUndlist[i];

      // do work
      bool didcomp = false;
      int reqID = reEncodeOffset(targetBucketIndex, glblNode.offsetID);

      descend = processNode(glblNode.node, s, chunk, reqID, isRoot, didcomp, awi);

#ifdef CHANGA_REFACTOR_INTERLIST_PRINT_LIST_STATE
      int tpindex = ownerTP->getIndex();
      if(targetBucketIndex == TEST_BUCKET && tpindex == TEST_TP){ 
        char arr[2] = {'K', 'D'};
        char arrr[2] = {'N', 'Y'};
        Vector3D<double> vec = ownerTP->decodeOffset(glblNode.offsetID);
        CkPrintf("[%d] undecided level %d: key %ld (%1.0f,%1.0f,%1.0f) - %s, target %d, ret: %c, comp: %c\n", tpindex, level, glblNode.node->getKey(), vec.x, vec.y, vec.z, typeString(glblNode.node->getType()), targetBucketIndex, arr[!descend], arrr[didcomp]);
      }
#endif
    }
  }

  // while there are nodes to process on this level
  while(!chklist.isEmpty())
  {
    // get remote node from state
    const OffsetNode &glblNode = chklist.deq();

    // do work
    bool didcomp = false;
    int reqID = reEncodeOffset(targetBucketIndex, glblNode.offsetID);
    descend = processNode(glblNode.node, s, chunk, reqID, isRoot, didcomp, awi);
#ifdef CHANGA_REFACTOR_INTERLIST_PRINT_LIST_STATE
    char arr[2] = {'K', 'D'};
    char arrr[2] = {'N', 'Y'};
    int tpindex = ownerTP->getIndex();
    if(targetBucketIndex == TEST_BUCKET && tpindex == TEST_TP){ 
      Vector3D<double> vec = ownerTP->decodeOffset(glblNode.offsetID);
      CkPrintf("[%d] chklist level %d: key %ld (%1.0f,%1.0f,%1.0f) - %s, target %d, ret: %c, comp: %c\n", tpindex, level, glblNode.node->getKey(), vec.x, vec.y, vec.z, typeString(glblNode.node->getType()), targetBucketIndex, arr[!descend], arrr[didcomp]);
    }
#endif
  }
  // if the undecided list is non-empty, have to let
  // children decide what to do with those nodes
  // descend into localNode
  bool myUndlistEmpty = myUndlist.length() == 0;
  //CkAssert((myUndlistEmpty && !descend) || (!myUndlistEmpty && descend));

  if(!myUndlistEmpty)
  {
    int which = localNode->whichChild(targetKey);
    GenericTreeNode *child = localNode->getChildren(which);
    CkAssert(child);
    dft(child, s, chunk, targetBucketIndex, false, awi, level+1);
  }
  else{
    // will not descend, must be lowest point on path
    // set lowest node in state
    s->lowestNode = localNode;
  }
  CkAssert(s->lowestNode != 0);
}

bool LocalTargetWalk::processNode(
                   GenericTreeNode *glblNode,
                   State *state,
                   int chunk, int reqID,
                   bool isRoot, bool &didcomp,
                   int awi)
{

  // assume that we won't be descending into the outer tree
  // if we are advised otherwise by the compute, descend
  // is set to true and the walk continued
  bool descend= false;
  // bucket number in reqID below does matter, because
  // we need to keep track of the target bucket.
  // only the remoteWalk is supposed to request
  // nodes; if a local walk started requesting remote nodes,
  // we'd be in trouble - awi is -1 for local walks because
  // they are not supposed to request missing nodes

  int ret = comp->doWork(glblNode, this, state, chunk, reqID, isRoot, didcomp, awi);
  if(ret == KEEP){
    descend = true;
  }
  else if(ret == DUMP){
    // a computation was carried out, need to
    // 'return'
    descend = false;
  }
  return descend;
}

// Here, node is the global node that we missed on earlier
// At this point, the 'source' node, i.e. the local node at which the localTargetWalk is to begin,has been
// set as the computeEntity for the compute. Also, the TW and compute  have been reassociated.
// The target bucket has been obtained from reqID and set as the target of the walk.

// This walk interprets what is otherwise the 'reqID' argument as the targetBucketIndex
void LocalTargetWalk::resumeWalk(GenericTreeNode *node, State *state_, int chunk, int reqID, int activeWalkIndex){
	DoubleWalkState *state = (DoubleWalkState *)state_;

    // initial target
	int targetBucket = decodeReqID(reqID);
	GenericTreeNode *source = (GenericTreeNode *)comp->getComputeEntity();
	// first and last buckets beneath source
	int startBucket, endBucket;
	int prevBucket = -1;

	// so that dummySource is changed if necessary, and source is left
	// untouched.
	GenericTreeNode *dummySource = source;
	ownerTP->getBucketsBeneathBounds(dummySource, startBucket, endBucket);

	// clear all levels up to and including source
	// this is so that we don't include the lists from
	// other resumed walks in the computation for this
	// instance
	// init state
	bool localLists = state->lplists.length() > 0;
	bool remoteLists = state->rplists.length() > 0;

	int level = source->getLevel(source->getKey());
	for(int i = 0; i <= level; i++){
		CheckList &chklist = state->chklists[i];
		while(!chklist.isEmpty()){
			chklist.deq();
		}
		state->clists[i].length() = 0;
		state->undlists[i].length() = 0;
	}

	if(localLists)
		for(int i = 0; i <= level; i++){
			state->lplists[i].length() = 0;
		}

	if(remoteLists)
		for(int i = 0; i <= level; i++){
			state->rplists[i].length() = 0;
		}

	// enqueue
	OffsetNode on;
	on.node = node;
	on.offsetID = reqID;
	state->chklists[level].enq(on);

	int activeRung = comp->getActiveRung();
	while(targetBucket < endBucket){
		GenericTreeNode *targetNode = ownerTP->getBucket(targetBucket);
		if(targetNode->rungs >= activeRung){
			targetKey = targetNode->getKey();
			GenericTreeNode *lca = ownerTP->getStartAncestor(targetBucket, prevBucket, source);

			int lcaLevel = lca->getLevel(lca->getKey());
			dft(lca, state, chunk, targetBucket, (lcaLevel == level), activeWalkIndex, lcaLevel);

			GenericTreeNode *lowestNode = state->lowestNode;
			// start and end buckets beneath lowest node
			int sL, eL;
			ownerTP->getBucketsBeneathBounds(lowestNode, sL, eL);

			CkAssert(targetBucket >= sL);
			comp->stateReady(state, ownerTP, chunk, targetBucket, eL);

			prevBucket = targetBucket;
			targetBucket = eL;
		}
		else{
			if(targetBucket != 0){
				prevBucket = targetBucket;
			}
			while(targetBucket < endBucket && ownerTP->getBucket(targetBucket)->rungs < activeRung){
				targetBucket++;
			}
		}
	}
	comp->setComputeEntity(source);
}
#endif

const char *translations[] = {"",
                                 "Invalid",
                                 "Bucket",
                                 "Internal",
                                 "Boundary",
                                 "NonLocal",
                                 "Empty",
                                 "Top",
                                 "NonLocalBucket",
                                 "Cached",
                                 "CachedBucket",
                                 "CachedEmpty"
                                 };

const char *typeString(NodeType type){
  CkAssert(type > 0 && type <= NUM_NODE_TYPES);
  return translations[type];
}
