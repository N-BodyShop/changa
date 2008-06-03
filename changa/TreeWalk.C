//#include "codes.h"

#include "config.h"
#include "GenericTreeNode.h"
#include "ParallelGravity.h"
//include "TreeWalk.h"
//#include "State.h"
#include "Compute.h"

void TreeWalk::init(Compute *c, TreePiece *owner){
  comp = c;
  ownerTP = owner;
}

void TreeWalk::reassoc(Compute *c){
  comp = c;
}

void TopDownTreeWalk::walk(GenericTreeNode *startNode, State *state, int chunk, int reqID, int awi){
#ifndef CHANGA_REFACTOR_WALKCHECK
  dft(startNode, state, chunk, reqID, true, awi);       // isRoot
#else
  bool doprint = false;
  if(comp->getSelfType() == Gravity && comp->getOptType() == Remote){
    doprint = (ownerTP->getCurrentRemoteBucket() == CHECK_BUCKET && ownerTP->getIndex() == CHECK_INDEX);
  }
  else if(comp->getSelfType() == Gravity && comp->getOptType() == Local){
    doprint = (ownerTP->getCurrentBucket() == CHECK_BUCKET && ownerTP->getIndex() == CHECK_INDEX);
  }
  if(doprint){
    CkPrintf("Printing walk (%d: type %c)\n", ownerTP->getIndex(), comp->getOptType() == Local ? 'L' : 'R');
  }
  dft(startNode, state, chunk, reqID, true, 0, doprint);       // isRoot
#endif
}

#ifndef CHANGA_REFACTOR_WALKCHECK
void TopDownTreeWalk::dft(GenericTreeNode *node, State *state, int chunk, int reqID, bool isRoot, int awi){
#else
void TopDownTreeWalk::dft(GenericTreeNode *node, State *state, int chunk, int reqID, bool isRoot, int shift, bool doprint){
#endif
  int ret;
  GenericTreeNode *tmp;
  if(node == NULL){   // something went wrong here  
    ckerr << "TopDownTreeWalk recvd. null node - chunk("<<chunk<<"), reqID("<<reqID<<"), isRoot("<<isRoot<<")" << endl;
    CkAbort("Abort");
  }
  currentGlobalKey = node->getKey();
  // process this node
  bool didcomp = false;
  ret = comp->doWork(node, this, state, chunk, reqID, isRoot, didcomp, awi);

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
      currentGlobalKey = node->getChildKey(i);

      comp->startNodeProcessEvent(ownerTP);
      /*
      if(comp->getSelfType() == Prefetch){
        ownerTP->incPrefetchWaiting();
      }
      */
      
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
        child = ownerTP->nodeMissed(reqID, node->remoteIndex, currentGlobalKey, chunk, comp->getSelfType() == Prefetch, awi);
        if(child == NULL){     // missed in cache, skip node for now
#if CHANGA_REFACTOR_DEBUG > 2
          CkPrintf("[%d]: child not found in cache\n", ownerTP->getIndex());
#endif
          comp->nodeMissedEvent(ownerTP, chunk);
          /*
          if(comp->getSelfType() == Gravity && comp->getOptType() == Remote){
            ownerTP->addToRemainingChunk(chunk, +1);
          }
          */
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
   // return;
  comp->finishNodeProcessEvent(ownerTP);
  /*
  if(comp->getSelfType() == Prefetch){
    if(ownerTP->decPrefetchWaiting() == 0){
      ownerTP->startRemoteChunk();
    }    
  }
  */
  return;
}

void TopDownTreeWalk::reset(){
  return;
}


// returns true if an ancestor of 'node' remained unopened in a prior tree walk and was
// therefore used in computation - helps avoid the duplication of computations 
// the bucket that is used in this decision belongs to 'comp' and therefore
// needn't be passed to the function
// another reason the bucket isn't passed as an argument to this method
// is that TopDownTreeWalk should constrain itself to work only for compute
// objects representing buckets of particles.
bool TopDownTreeWalk::ancestorCheck(GenericTreeNode *node, int reqID){
  GenericTreeNode *ancestor;
  CkAssert(comp != NULL);

  ancestor = ownerTP->getRoot();
  while(ancestor != NULL){
    // comp will know what to do with this TP and ancestor
    if(!comp->openCriterion(this->getOwnerTP(), ancestor, reqID)){
      return true;
    }
    int which = ancestor->whichChild(node->getKey());
    //CkPrintf("%d: ancestorCheck - getChildren(%d)\n", ownerTP->getIndex(), which);
    ancestor = ancestor->getChildren(which);
  }
  return false;
}

/*
void DoubleWalk::walk(GenericTreeNode *startNode, State *_state, int reqID){
  // no need for state ?
  doubleDft(reqID);
}

void DoubleWalk::doubleDft(GenericTreeNode *ancestor, State *state, int reqID...){
}
*/

/*
void BucketIterator::init(Compute *c, TreePiece *owner){
  comp = c;
  ownerTP = owner;
}

void BucketIterator::reset(){
  return;
}

void BucketIterator::walk(GenericTreeNode *node, State *state, int chunk, int reqID){
  iterate(node, state, chunk, reqID);
}

}
 
void BucketIterator::iterate(GenericTreeNode *node, 
*/
