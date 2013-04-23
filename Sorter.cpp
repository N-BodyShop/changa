/** \file Sorter.cpp
 Implementation of the parallel sort.
 \author originally Graeme Lufkin (gwl@u.washington.edu)
 */

#include <algorithm>

#include "Sorter.h"

using std::vector;
using std::list;
using std::set;

using namespace std;
using namespace SFC;

extern CProxy_TreePiece treeProxy;

/***************ORB Decomposition*****************/
/// @brief Start ORB decomposition
/// @param m message with bounding box of the simulation
///
/// Begins the ORB decomposition by calculating the initial split and
/// broadcasting to TreePiece::evaluateParticleCounts().  Particle
/// counts will be contributed to Sorter::collectORBCounts().
///
void Sorter::doORBDecomposition(CkReductionMsg* m){

  float len=0.0,len2=0.0;
  char dim;
  double pos;
 	
  OrientedBox<float> box = *static_cast<OrientedBox<float> *>(m->getData());
  delete m;

  if(numChares == 1) { // No decomposition to do
      treeProxy[0].initBeforeORBSend(0,0,0,sortingCallback,
				     CkCallback(CkIndex_Sorter::readytoSendORB(0), thishandle));
      return;
      }

  ORBData single;
  single.boundingBox = box;
  
  //Find which is the longest dimension
  len=box.greater_corner.x-box.lesser_corner.x;
  dim=0;
  if(len<0.0) { len = -len; }
  len2=box.greater_corner.y-box.lesser_corner.y;
  if(len2<0.0) { len2 = -len2; }
  if(len2>len) { len = len2; dim=1; }
  len2=box.greater_corner.z-box.lesser_corner.z;
  if(len2<0.0) { len2 = -len2; }
  if(len2>len) { len = len2; dim=2; }
  
  if(box.lesser_corner[dim] <= box.greater_corner[dim]){
    pos = box.lesser_corner[dim] + len/2;
  }
  else{
    pos = box.greater_corner[dim] + len/2;
  }
 
  single.curLow = box.lesser_corner[dim];
  single.curHigh = box.greater_corner[dim];
  single.curDivision = pos;
  single.curDim = dim;
  
  orbData.clear();
  orbData.push_back(single);

	ORBSplittersMsg * splittersMsg = new (1,1) ORBSplittersMsg(1,CkCallback(CkIndex_Sorter::collectORBCounts(0), thishandle));
	splittersMsg->pos[0] = pos;
	splittersMsg->dim[0] = dim;
  treeProxy.evaluateParticleCounts(splittersMsg);
}

/// Calculate candidate divisions for next level in the ORB
/// decomposition tree. If we have enough pieces, proceed to
/// sending particles. 
void Sorter::finishPhase(CkReductionMsg *m){

  float len=0.0,len2=0.0;
  char dim;
  double pos;
  int i,listSize;
  
  delete m;
  std::list<ORBData>::iterator iter;
  std::list<ORBData>::iterator iter2;
  iter = orbData.begin();
	listSize = orbData.size();

  for(int i=0;i<listSize;i++){
    ORBData first,second;
    first.boundingBox.lesser_corner = (*iter).boundingBox.lesser_corner;
    first.boundingBox.greater_corner = (*iter).boundingBox.greater_corner;
    first.boundingBox.greater_corner[(*iter).curDim] = (*iter).curDivision;

    second.boundingBox.lesser_corner = (*iter).boundingBox.lesser_corner;
    second.boundingBox.lesser_corner[(*iter).curDim] = (*iter).curDivision;
    second.boundingBox.greater_corner = (*iter).boundingBox.greater_corner;

    //Find which is the longest dimension
    OrientedBox<float> box1 = first.boundingBox;
    OrientedBox<float> box2 = second.boundingBox;

    len=box1.greater_corner.x-box1.lesser_corner.x;
    dim=0;
    CkAssert(len>=0.0);
    //if(len<0.0) { len = -len; }
    len2=box1.greater_corner.y-box1.lesser_corner.y;
    CkAssert(len2>=0.0);
    //if(len2<0.0) { len2 = -len2; }
    if(len2>len) { len = len2; dim=1; }
    len2=box1.greater_corner.z-box1.lesser_corner.z;
    CkAssert(len2>=0.0);
    //if(len2<0.0) { len2 = -len2; }
    if(len2>len) { len = len2; dim=2; }
    pos = box1.lesser_corner[dim] + len/2;

    first.curLow = box1.lesser_corner[dim];
    first.curHigh = box1.greater_corner[dim];
    first.curDivision = pos;
    first.curDim = dim;

    len=box2.greater_corner.x-box2.lesser_corner.x;
    dim=0;
    CkAssert(len>=0.0);
    //if(len<0.0) { len = -len; }
    len2=box2.greater_corner.y-box2.lesser_corner.y;
    CkAssert(len2>=0.0);
    //if(len2<0.0) { len2 = -len2; }
    if(len2>len) { len = len2; dim=1; }
    len2=box2.greater_corner.z-box2.lesser_corner.z;
    CkAssert(len2>=0.0);
    //if(len2<0.0) { len2 = -len2; }
    if(len2>len) { len = len2; dim=2; }
    pos = box2.lesser_corner[dim] + len/2;

    second.curLow = box2.lesser_corner[dim];
    second.curHigh = box2.greater_corner[dim];
    second.curDivision = pos;
    second.curDim = dim;

    //Insert both the sub-divisions into the list and remove the bigger division
    iter2 = orbData.insert(iter,first);
    iter2 = orbData.insert(iter,second);
    iter = orbData.erase(iter);
  }
 
  if(numChares == orbData.size()){ //Move data around
    for(i=0;i<numChares;i++){
	if(verbosity > 2) {
	    CkPrintf("%d has %d particles\n",i,binCounts[i]);
	    CkPrintf("%d has %d gas particles\n",i,binCountsGas[i]);
	    CkPrintf("%d has %d star particles\n",i,binCountsStar[i]);
	    }
	treeProxy[i].initBeforeORBSend(binCounts[i], binCountsGas[i],
				       binCountsStar[i], sortingCallback,
				       CkCallback(CkIndex_Sorter::readytoSendORB(0), thishandle));
	}
  }
  else{ //Send the next phase of splitters
    ORBSplittersMsg *splittersMsg = new (orbData.size(),orbData.size()) ORBSplittersMsg(orbData.size(),CkCallback(CkIndex_Sorter::collectORBCounts(0), thishandle));
    for(i=0,iter=orbData.begin();iter!=orbData.end();i++,iter++){
      splittersMsg->pos[i] = (*iter).curDivision;
      splittersMsg->dim[i] = (*iter).curDim;
    }
    treeProxy.evaluateParticleCounts(splittersMsg);
  }

}

void Sorter::readytoSendORB(CkReductionMsg* m){
  delete m;

/*
 * Send information to DataManager, then broadcast to send particles.
 */
  dm.acceptResponsibleIndex(&(*chareIDs.begin()), chareIDs.size(),
			    CkCallback(CkIndex_TreePiece::sendORBParticles(),
				       treeProxy));
}

/// @brief Collect particle counts from treepieces and send out new
/// splits.
/// @param m A message with the summed counts for the current ORB
/// splits.
///
/// If the counts are within the tolerances, call
/// TreePiece::finalizeBoundaries().
///
void Sorter::collectORBCounts(CkReductionMsg* m){

  std::list<ORBData>::iterator iter;
  int i;
  
  numCounts = m->getSize() / (3*sizeof(int)); // three separate arrays for
					    // total, gas and stars
  binCounts.resize(numCounts);
  binCountsGas.resize(numCounts);
  binCountsStar.resize(numCounts);
  int* startCounts = static_cast<int *>(m->getData());
  copy(startCounts, startCounts + numCounts, binCounts.begin());
  copy(startCounts + numCounts, startCounts + 2*numCounts,
       binCountsGas.begin());
  copy(startCounts + 2*numCounts, startCounts + 3*numCounts,
       binCountsStar.begin());
  delete m;

  CkAssert(numCounts == 2*orbData.size());
  
  ORBSplittersMsg *splittersMsg;
  float TOLER=0.05;
  int doneCount=0;
  
  for(i=0,iter=orbData.begin(); iter!=orbData.end(); i++,iter++){
      if((binCounts[2*i+1]*(1-TOLER)<=binCounts[2*i]
	  && binCounts[2*i]<=(1+TOLER)*binCounts[2*i+1])
	 || (domainDecomposition == ORB_space_dec)){
      doneCount++;
    }
    else{
      if(binCounts[2*i] > binCounts[2*i+1]){
        (*iter).curHigh = (*iter).curDivision;
        (*iter).curDivision = ((*iter).curLow + (*iter).curHigh)/2;
      }
      else{
        (*iter).curLow = (*iter).curDivision;
        (*iter).curDivision = ((*iter).curLow + (*iter).curHigh)/2;
      }
    }
  }

  
  //Assuming that lesser corner is always smaller than greater corner
  if(doneCount==orbData.size()){
    splittersMsg = new (orbData.size(),orbData.size()) ORBSplittersMsg(orbData.size(),CkCallback(CkIndex_Sorter::finishPhase(0), thishandle));
    for(i=0,iter=orbData.begin();iter!=orbData.end();i++,iter++){
      splittersMsg->pos[i] = (*iter).curDivision;
      splittersMsg->dim[i] = (*iter).curDim;
    }
    //finalize the boundaries in all the Treepieces
    treeProxy.finalizeBoundaries(splittersMsg);
  }
  else{
    splittersMsg = new (orbData.size(),orbData.size()) ORBSplittersMsg(orbData.size(),CkCallback(CkIndex_Sorter::collectORBCounts(0), thishandle));
    for(i=0,iter=orbData.begin();iter!=orbData.end();i++,iter++){
      splittersMsg->pos[i] = (*iter).curDivision;
      splittersMsg->dim[i] = (*iter).curDim;
    }
    treeProxy.evaluateParticleCounts(splittersMsg);
  }
  
}

/**
 * \brief Overall start of domain decomposition
 * @param dataManagerID ID of data manager group
 * @param toler Tolerance within which to have the same number of particles
 * @param decompose Are we still deciding on decomposition?
 */
void Sorter::startSorting(const CkGroupID& dataManagerID,
			  const double toler, bool decompose) {
	numChares = numTreePieces;
	dm = CProxy_DataManager(dataManagerID);
	tolerance = toler;
	sorted = false;
	numIterations = 0;
	
  //Changed for implementing OCT decomposition
  Key delta;
  Key k;
  BinaryTreeNode *rt;

#ifdef REDUCTION_HELPER
  // The reduction helper needs to know the number of pieces on each processor.
  reductionHelperProxy.countTreePieces(CkCallbackResumeThread());
  CProxy_ReductionHelper boundariesTargetProxy = reductionHelperProxy; 
#else
  CProxy_TreePiece boundariesTargetProxy = treeProxy; 
#endif

  decompTime = CmiWallTimer();
    
  switch (domainDecomposition){
    case SFC_dec:
    case SFC_peano_dec:
    case SFC_peano_dec_3D:
    case SFC_peano_dec_2D:
        numKeys = 0;
        if (keyBoundaries.size() == 0) {
	    splitters.clear();
	    int nSplitters = 4*numChares + 1;
	    splitters.reserve(nSplitters);
	    delta = (lastPossibleKey - SFC::firstPossibleKey) / (nSplitters-1);
	    k = firstPossibleKey;
	    for(int i = 0; i < (nSplitters-1); i++, k += delta) {
		if(k != firstPossibleKey)
		    k |= 7L;  // Set bottom bits to avoid trees too deep
		splitters.push_back(k);
		}
	    splitters.push_back(lastPossibleKey);
            }
	else { // reuse the existing splitters from the previous decomposition
	    splitters.assign(keyBoundaries.begin(), keyBoundaries.end());
	    }
        break;

    case Oct_dec:
      {
        refineLevel = octRefineLevel;

        rt = new BinaryTreeNode();
        int numInitialBins = numInitDecompBins;

        if (nodeKeys.size() == 0) {
          nodeKeys.reserve(numInitialBins);
          nodeKeys.resize(numInitialBins>>1, 0);
          NodeKey *tmp = &(*nodeKeys.begin());
          rt->getChunks(nodeKeys.size(),tmp);
        }
        delete rt;

        activeNodes->clear();
        tmpActiveNodes->clear();

        numDecompRoots = nodeKeys.size();
        if (decompRoots == NULL) {
          decompRoots = new OctDecompNode[numDecompRoots];
          for(int i = 0; i < nodeKeys.size(); i++){
            decompRoots[i].key = nodeKeys[i];
            //CkPrintf("init add %llu to activeNodes\n", decompRoots[i].key);
            activeNodes->push_back(&decompRoots[i]);
          }
        }
        else {
          root->getLeafNodes(activeNodes);
        }

        //CkPrintf("Sorter: initially %d keys\n", activeNodes->length());

	if(!joinThreshold) {
	    joinThreshold = particlesPerChare;
	    splitThreshold = (int)(joinThreshold * 1.5);
	    }
        
        //Convert the Node Keys to the splitter keys which will be sent to histogram
        convertNodesToSplitters();
        break;
      }

    case ORB_dec:
    case ORB_space_dec:
	numKeys = 0;
      treeProxy.initORBPieces(CkCallback(CkIndex_Sorter::doORBDecomposition(0), thishandle));
    break;
      default:
    CkAbort("Invalid domain decomposition requested");
  }
  
  if(verbosity >= 3)
    ckout << "Sorter: Initially have " << splitters.size() << " splitters" << endl;
  

  //send out the first guesses to be evaluated
  if((domainDecomposition!=ORB_dec) && (domainDecomposition!=ORB_space_dec)) {

    if (decompose) {
      // XXX: Optimizations available if sort has been done before!
      //create initial evenly distributed guesses for splitter keys
      keyBoundaries.clear();
      keyBoundaries.reserve(numChares + 1);
      keyBoundaries.push_back(firstPossibleKey);      
    } else {
      //send out all the decided keys to get final bin counts
      sorted = true;     
    }

    std::vector<SFC::Key>* keys; 
    if(decompose || domainDecomposition == Oct_dec){
      keys = &splitters; 
    }
    else {
      keys = &keyBoundaries; 
    }

    boundariesTargetProxy.evaluateBoundaries(&(*keys->begin()), keys->size(), 0, CkCallback(CkIndex_Sorter::collectEvaluations(0), thishandle));
  }
}

/**
 * Given "numKeys" node keys ("nodeKeys"), convert these keys into splitters.
 */
void Sorter::convertNodesToSplitters(){
  Key partKey;

  splitters.clear();
  splitters.reserve(nodeKeys.size() + 1);
  //binCounts.reserve(nodeKeys.size());
  const Key mask = Key(1) << KeyBits;
  for(unsigned int i=0;i<nodeKeys.size();i++){
    partKey=Key(nodeKeys[i]);
    while(!(partKey & mask)){
      partKey <<= 1;
    }
    partKey &= ~mask;
    if(partKey != 0) partKey--;
    splitters.push_back(partKey);
  }

  // Note: Splitters are guaranteed to be sorted.
  // Even if nodeKeys are not completely sorted (see getChunks),
  // by construction the left shifting will produce a sorted splitters vector
  splitters.push_back(lastPossibleKey);
}

/**
 * Given "num" node keys, create splitters for these selected nodekeys
 * to refine the information within those nodes: each node in keys is
 * divided a number of times that depends on "refineLevel". Returns a
 * newly allocated Key array containing a concatenation of the
 * splitter keys for each node to be refined for the histogramming phase.
 */
Key * Sorter::convertNodesToSplittersRefine(int num, NodeKey* keys){
  Key partKey = Key(0);

  Key *result = new Key[num * ((1<<refineLevel)+1)];
  Key levelMask = Key(1) << KeyBits;
  levelMask >>= refineLevel;
  int idx = 0;
  const Key mask = Key(1) << KeyBits;
  for(unsigned int i=0;i<num;i++){
    CkAssert(! (partKey & levelMask));
    partKey=Key(keys[i]<<refineLevel);
    //CkPrintf("convertRefine %llx -> ", keys[i]);
    int shift = 0;
    // find how much we need to shift each key (depend on the tree level of the key)
    while (!(partKey<<shift & mask)) {
      ++shift;
    }
    partKey &= ~mask >> shift;
    for (int j=0; j<=(1<<refineLevel); ++j) {
      Key kResult =  ((partKey+j) << shift);
      if(kResult != 0) kResult--;
      result[idx++] = kResult;
      //CkPrintf("%llx,", kResult);
    }
    //CkPrintf("\n");
  }

  // Note: Splitters are guaranteed to be sorted.
  // By construction the left shifting will produce a sorted splitters vector

  return result;
}

void Sorter::collectEvaluations(CkReductionMsg* m) {

  switch (domainDecomposition){
    case SFC_dec:
    case SFC_peano_dec:
    case SFC_peano_dec_3D:
    case SFC_peano_dec_2D:
      collectEvaluationsSFC(m);
      break;
    case Oct_dec:
      collectEvaluationsOct(m);
      break;
    case ORB_dec:
    case ORB_space_dec:
      CkAbort("ORB: We shouldn't have reached here");
    break;
      default:
    CkAbort("Invalid domain decomposition requested");
  }
}

/**
 * Examine the counts for Oct decomposition and determine if further
 * refining is needed.  Call TreePiece::evaluateBoundaries if needed,
 * else send the final keys to DataManager::acceptFinalKeys.
 * Sorter::refineOctSplitting perfoms splits where needed.
 */
void Sorter::collectEvaluationsOct(CkReductionMsg* m) {

  numIterations++;
  numCounts = m->getSize() / sizeof(int);
  int* startCounts = static_cast<int *>(m->getData());

  //call function which will balance the bin counts: define it in GenericTreeNode
  //make it a templated function
  //Pass the bincounts as well as the nodekeys

  if (joinThreshold == 0) {
    int total_particles = std::accumulate(startCounts, startCounts+numCounts, 0);
    joinThreshold = total_particles / (numTreePieces>>1);
    splitThreshold = (int) (joinThreshold * 1.5);
  }
  
  if(verbosity>=3){
    int i=0;
    CkPrintf("Bin Counts in collect eval (%d):",numCounts);
    for ( ; i<numCounts; i++) {
      CkPrintf("%d,",startCounts[i]);
    }
    CkPrintf("\n");
    CkPrintf("Nodekeys:");
    for(int j=0;j<nodeKeys.size();j++)
      CkPrintf("%llx,",nodeKeys[j]);
    CkPrintf("\n");
    if (nodesOpened.size() > 0) {
      CkPrintf("Nodes opened (%d):",nodesOpened.size());
      for (int j=0;j<nodesOpened.size();j++)
        CkPrintf("%llx,",nodesOpened[j]);
      CkPrintf("\n");
    }
  }

  double startTimer = CmiWallTimer();
  bool histogram = refineOctSplitting(numCounts, startCounts);
  traceUserBracketEvent(weightBalanceUE, startTimer, CmiWallTimer());

  //CkPrintf("refineOctSplitting nodesOpened %d took %g s\n", nodesOpened.size(), CmiWallTimer()-startTimer);
  delete m;

  if(verbosity>=3){
    CkPrintf("Nodekeys after (%d):",nodeKeys.size());
    for(int i=0;i<nodeKeys.size();i++)
      CkPrintf("%llx,",nodeKeys[i]);
    CkPrintf("\n");
  }
  if(histogram){
    refineLevel = octRefineLevel;
    int arraySize = (1<<refineLevel)+1;
    startTimer = CmiWallTimer();
    Key *array = convertNodesToSplittersRefine(nodesOpened.size(),nodesOpened.getVec());
    //CkPrintf("convertNodesToSplittersRefine elts %d took %g s\n", nodesOpened.size()*arraySize, CmiWallTimer()-startTimer);
#ifdef REDUCTION_HELPER
    CProxy_ReductionHelper boundariesTargetProxy = reductionHelperProxy; 
#else
    CProxy_TreePiece boundariesTargetProxy = treeProxy; 
#endif
    boundariesTargetProxy.evaluateBoundaries(array, nodesOpened.size()*arraySize, 1<<refineLevel, CkCallback(CkIndex_Sorter::collectEvaluations(0), thishandle));
    delete[] array;
  }
  else{
    sorted=true;
    splitters.clear();

    if (root == NULL) {
      // build a tree out of decompRoots
      // requires numDecompRoots to be a power of OctDecompNode::maxNumChildren
      CkVec<OctDecompNode*> *leaves = new CkVec<OctDecompNode*>;
      root = new OctDecompNode;
      root->key = 1;
      int depth = 0;
      for (int numLeaves = numDecompRoots; numLeaves > OctDecompNode::maxNumChildren; numLeaves /= OctDecompNode::maxNumChildren) {
        depth++;
      }

      root->makeSubTree(depth, leaves);
      // CkPrintf("num leaves: %d numDecompRoots: %d\n", leaves->length(), numDecompRoots);

      for (int i = 0; i < leaves->length(); i++) {
        (*leaves)[i]->nchildren = OctDecompNode::maxNumChildren;
        (*leaves)[i]->children = new OctDecompNode[OctDecompNode::maxNumChildren];
        for (int j = 0; j < OctDecompNode::maxNumChildren; j++) {
          (*leaves)[i]->children[j] = decompRoots[i * OctDecompNode::maxNumChildren + j];
        }
      }

      delete leaves;
    }

    int total = root->buildCounts();
    // CkPrintf("total number of particles: %d\n", total);
    do {
	// Convert Oct domains to splitters, ensuring that we do not exceed
	// the number of available TreePieces.
	nodeKeys.clear();
	binCounts.clear();

        root->combine(joinThreshold, nodeKeys, binCounts);

	if(binCounts.size() > numTreePieces) {
	    CkPrintf("bumping joinThreshold: %d, size: %d\n", joinThreshold,
		     binCounts.size());
	    joinThreshold = (int) (1.1*joinThreshold);
	    }
	}
    while(binCounts.size() > numTreePieces);
    convertNodesToSplitters();

#if 0
    CkPrintf("final splitters: ");
    for(int i = 0; i < splitters.size(); i++){
      CkPrintf("%llx,", splitters[i]);
    }
    CkPrintf("\n");

    CkPrintf("final chares: ");
    for(int i = 0; i < chareIDs.size(); i++){
      CkPrintf("%d,", chareIDs[i]);
    }
    CkPrintf("\n");

    CkPrintf("final counts: ");
    int totalAccum = 0;
    int maxPartCount = 0;
    for(int i = 0; i < binCounts.size(); i++){
      CkPrintf("%d,", binCounts[i]);
      totalAccum += binCounts[i];

      if (binCounts[i] > maxPartCount) {
        maxPartCount = binCounts[i];
      }

    }

    int avePartCount = totalAccum / binCounts.size();
    int numSmallCount = 0;
    for(int i = 0; i < binCounts.size(); i++){
      if (binCounts[i] < 0.1 * avePartCount) {
        numSmallCount++;
      }
    }
    CkPrintf("\ntotal count: %d\n", totalAccum);
    CkPrintf("average number of particles per tree piece: %d\n", avePartCount);
    CkPrintf("number of tree pieces of size less than 10%% of average %d\n", numSmallCount);
    CkPrintf("maximum tree piece size: %d\n", maxPartCount);

#endif

    ckout << "Sorter: Histograms balanced after " << numIterations
          << " iterations. Using " << nodeKeys.size() << " chares." << endl;
    //We have the splitters here because weight balancer didn't change any node keys
    //We also have the final bin counts
		
    if(binCounts.size() > numTreePieces){
      CkPrintf("Need %d tree pieces, available %d\n", binCounts.size(), numTreePieces);
      CkAbort("too few tree pieces\n");
    }

    /*
    CkPrintf("Done decomp\n");
    CkExit();
    return;
    */

    CkPrintf(" histogramming %g sec ... \n", CmiWallTimer()-decompTime);
    
    dm.acceptFinalKeys(&(*splitters.begin()), &(*chareIDs.begin()), &(*binCounts.begin()), splitters.size());
    numIterations = 0;
    sorted = false;
    return;
  }

}

int OctDecompNode::maxNumChildren = 2;
int OctDecompNode::lgMaxNumChildren = 1;

void OctDecompNode::makeSubTree(int refineLevel, CkVec<OctDecompNode*> *active){
  if(refineLevel == 0){
    active->push_back(this);
    //CkPrintf("push %llx into tmpActive\n", key);
    return;
  }

  children = new OctDecompNode[maxNumChildren];
  nchildren = maxNumChildren;

  NodeKey childKey = (key << lgMaxNumChildren);
  for(int i = 0; i < nchildren; i++){
    children[i].key = childKey;
    children[i].makeSubTree(refineLevel-1,active);
    childKey++;
  }
}

int OctDecompNode::buildCounts() {
  if (children == NULL) {
    return nparticles;
  }
  else {
    nparticles = 0;
    for (int i = 0; i < nchildren; i++) {
      nparticles += children[i].buildCounts();
    }
    return nparticles;
  }
}

void OctDecompNode::combine(int joinThreshold, vector<NodeKey> &finalKeys, vector<unsigned int> &counts){
  if(nparticles < joinThreshold || nchildren == 0){
    finalKeys.push_back(key);
    counts.push_back(nparticles);
    deleteBeneath();
    return;
  }

  for(int i = 0; i < nchildren; i++){
    children[i].combine(joinThreshold, finalKeys, counts);
  }

}

void OctDecompNode::deleteBeneath(){
  if(nchildren == 0) return;

  for(int i = 0; i < nchildren; i++){
    children[i].deleteBeneath();
  }
  delete[] children;

  children = NULL;
  nchildren = 0;
}

 void OctDecompNode::getLeafNodes(CkVec<OctDecompNode*> *activeNodes) {
   if (nchildren == 0) {
     activeNodes->push_back(this); 
   }
   else {
     for(int i = 0; i < nchildren; i++){
       children[i].getLeafNodes(activeNodes);
     }
   }
 }

/**
 * This function uses nodeKeys as tree holding the current status of
 * the decomposition. It modifies it to decide when it is time to join two nodes,
 * or to split a node (meaning a TreePiece holding that part of the tree rooted
 * at "node").
 * As a consequence of this modification, "splitters" "chareIDs" and "binCounts"
 * are updated accordingly to reflect the changes. If some node requires extra
 * refinement, "nodesOpened" will be changed to reflect this request for more data.
 * Returns true if more refinement is requested.
 */
bool Sorter::refineOctSplitting(int n, int *count) {
  unsigned int nprocess = 0;
  unsigned int nopen = 0;
  unsigned int njoin = 0;

  CkAssert(activeNodes->length() == n);

  nodesOpened.clear();

  for(int i = 0; i < n; i++){
    OctDecompNode *parent = (*activeNodes)[i];
    parent->nparticles = count[i];
    if(parent->nparticles > splitThreshold){
      // create a subtree of depth 'refineLevel' underneath 'parent'
      // newly created children are pushed into 'tmpActiveNodes'
      // the key of the parent is placed in 'nodesOpened' so that we 
      // can make splitters out of the childrens' keys
      parent->makeSubTree(refineLevel, tmpActiveNodes);
      nodesOpened.push_back(parent->key);
    }
  }

  CkVec<OctDecompNode*> *save = tmpActiveNodes;
  tmpActiveNodes = activeNodes;
  activeNodes = save;

  tmpActiveNodes->length() = 0;

  //CkPrintf("Sorter: refined to get %d keys\n", activeNodes->length());

  return (nodesOpened.size() > 0);

}

/**
 * \brief Collect evaluations for the SFC domain decomposion.
 *
 * Examines the bin counts to see if the iteration has converged.  If
 * yes, a final count is done and then sent to
 * DataManager::acceptFinalKeys.  If not, then Sorter::adjustSplitters
 * is used to refine the search intervals.
 */
void Sorter::collectEvaluationsSFC(CkReductionMsg* m) {
	numIterations++;
	numCounts = m->getSize() / sizeof(int);
	binCounts.resize(numCounts + 1);
	binCounts[0] = 0;
	int* startCounts = static_cast<int *>(m->getData());
	copy(startCounts, startCounts + numCounts, binCounts.begin() + 1);
	delete m;
	
	if(sorted) { //true after the final keys have been binned
		//send out the final splitters and responsibility table
		dm.acceptFinalKeys(&(*keyBoundaries.begin()), &(*chareIDs.begin()), &(*binCounts.begin()) + 1, keyBoundaries.size());
		numIterations = 0;
		sorted = false;
		return;
	}
	
	if(verbosity >= 4)
		ckout << "Sorter: On iteration " << numIterations << endl;
	CkAssert(numIterations < 1000);  // Sorter has not converged.
	
	//sum up the individual bin counts, so each bin has the count of it and all preceding
	partial_sum(binCounts.begin(), binCounts.end(), binCounts.begin());
	
	if(!numKeys) {
		numKeys = binCounts.back();
		int avgValue = numKeys / numChares;
		closeEnough = static_cast<int>(avgValue * tolerance);
		if(closeEnough < 0 || closeEnough >= avgValue) {
			ckerr << "Sorter: Unacceptable tolerance, requiring exact fit." << endl;
			closeEnough = 0;
		}
		
		//each splitter key will split the keys near a goal number of keys
		goals.assign(numChares - 1, avgValue);
		// evenly distribute extra particles.
		int rem = numKeys % numChares;
		int j = 0;
		for(list<int>::iterator i = goals.begin(); j < rem; j++, ++i) {
		    *i = avgValue + 1;
		    }
		partial_sum(goals.begin(), goals.end(), goals.begin());
		
		if(verbosity >= 3)
			ckout << "Sorter: Target keys per chare: " << avgValue << " plus/minus " << (2 * closeEnough) << endl;
	}

	//make adjustments to the splitter keys based on the results of the previous iteration
	adjustSplitters();

	if(verbosity >= 4) {
		ckout << "Sorter: Probing " << splitters.size() << " splitter keys" << endl;
		ckout << "Sorter: Decided on " << (keyBoundaries.size() - 1) << " splitting keys" << endl;
	}

        std::vector<SFC::Key>* keys; 
	//check if we have found all the splitters
	if(sorted) {
          
		if(verbosity)
			ckout << "Sorter: Histograms balanced after " << numIterations << " iterations." << endl;
		
		sort(keyBoundaries.begin() + 1, keyBoundaries.end());
		keyBoundaries.push_back(lastPossibleKey);
		
		//send out all the decided keys to get final bin counts
                keys = &keyBoundaries; 
	} else {
          //send out the new guesses to be evaluated
          keys = &splitters; 
        }

#ifdef REDUCTION_HELPER
        CProxy_ReductionHelper boundariesTargetProxy = reductionHelperProxy; 
#else
        CProxy_TreePiece boundariesTargetProxy = treeProxy; 
#endif
        boundariesTargetProxy.evaluateBoundaries(&(*keys->begin()), keys->size(), 0, CkCallback(CkIndex_Sorter::collectEvaluations(0), thishandle));
}

/** Generate new guesses for splitter keys based on the histograms that came
 back from the last batch.
 We need to find the keys that split a distribution into even piles.
 We find the bracketing splits, and then insert a new guess in the
 middle for a new evaluation. This is a simultaneous bisection search
 for each splitter key not yet found.
 */
void Sorter::adjustSplitters() {
	
	set<Key> newSplitters;
	
	Key leftBound, rightBound;
	vector<unsigned int>::iterator numLeftKey, numRightKey = binCounts.begin();
	
	//for each goal not yet met (each splitter key not yet found)
	for(list<int>::iterator Ngoal = goals.begin(); Ngoal != goals.end(); ) {

		//find the positions that bracket the goal
		numRightKey = lower_bound(numRightKey, binCounts.end(), *Ngoal);
		numLeftKey = numRightKey - 1;
		
		if(numRightKey == binCounts.begin())
			ckerr << "Sorter: Looking for " << *Ngoal << " How could this happen at the beginning?" << endl;
		if(numRightKey == binCounts.end())
			ckerr << "Sorter: Looking for " << *Ngoal << " How could this happen at the end?" << endl;
		
		//translate the positions into the bracketing keys
		leftBound = splitters[numLeftKey - binCounts.begin()];
		rightBound = splitters[numRightKey - binCounts.begin()];
		
		//check if one of the bracketing keys is close enough to the goal
		if(abs((int)*numLeftKey - *Ngoal) <= closeEnough) {
			//add this key to the list of decided splitter keys
			keyBoundaries.push_back(leftBound);
			//the goal has been met, delete it
			list<int>::iterator temp = Ngoal;
			++Ngoal;
			goals.erase(temp);
		} else if(abs((int)*numRightKey - *Ngoal) <= closeEnough) {
			keyBoundaries.push_back(rightBound);
			list<int>::iterator temp = Ngoal;
			++Ngoal;
			goals.erase(temp);
		} else {
			// not close enough yet, add the bracketing keys and
			// the middle to the guesses
			// Set bottom bits to avoid trees to deep.
			newSplitters.insert(leftBound | 7L);
			newSplitters.insert((leftBound / 4 * 3 + rightBound / 4)
					    | 7L);
			newSplitters.insert((leftBound / 2 + rightBound / 2)
					    | 7L);
			newSplitters.insert((leftBound / 4 + rightBound / 4 * 3)
					    | 7L);
			newSplitters.insert(rightBound | 7L);
			++Ngoal;
		}
	}
	
	//if we don't have any new keys to probe, then we're done
	if(newSplitters.empty())
		sorted = true;
	else {
		//evaluate the new set of splitters
		newSplitters.insert(firstPossibleKey);
		newSplitters.insert(lastPossibleKey);
		// The following statement breaks with the PGI
		// compiler.  (version PGCC/x86 Linux/x86-64 6.1-2)
		// The following loop is a work around.
		// splitters.assign(newSplitters.begin(), newSplitters.end());
		splitters.clear();
		if(verbosity >=4 ) CkPrintf("Keys:");
		for(set<Key>::iterator iterNew = newSplitters.begin();
		    iterNew != newSplitters.end(); iterNew++) {
		    if(verbosity >= 4) CkPrintf("%lx,", *iterNew);
		    splitters.push_back(*iterNew);
		    }
		if(verbosity >=4 ) CkPrintf("\n");
	}
}
