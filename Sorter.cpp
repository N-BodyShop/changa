/** \file Sorter.cpp
 Implementation of the parallel sort.
 \author Graeme Lufkin (gwl@u.washington.edu)
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
void Sorter::doORBDecomposition(CkReductionMsg* m){

  float len=0.0,len2=0.0;
  char dim;
  double pos;
 	
  OrientedBox<float> box = *static_cast<OrientedBox<float> *>(m->getData());
  delete m;

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
  //phaseLeader=0;
  //lastPiece = numTreePieces-1;
  //for(int i=0;i<numTreePieces;i++)
    //thisProxy[i].evaluateFirstTime(pos,dim,phaseLeader);
}

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
 
	//CkPrintf("num chares:%d, partitions got:%d\n",numChares,orbData.size());

  if(numChares == orbData.size()){ //Move data around
    for(i=0;i<numChares;i++){
      if(verbosity > 1)
      CkPrintf("%d has %d particles\n",i,binCounts[i]);
			//treeProxy[i].sendORBParticles(binCounts[i],sortingCallback,CkCallback(CkIndex_Sorter::sendBoundingBoxes(0), thishandle));
			treeProxy[i].initBeforeORBSend(binCounts[i],sortingCallback,CkCallback(CkIndex_Sorter::readytoSendORB(0), thishandle));
		}
  }
  else{ //Send the next phase of splitters
    ORBSplittersMsg *splittersMsg = new (orbData.size(),orbData.size()) ORBSplittersMsg(orbData.size(),CkCallback(CkIndex_Sorter::collectORBCounts(0), thishandle));
    for(i=0,iter=orbData.begin();iter!=orbData.end();i++,iter++){
      splittersMsg->pos[i] = (*iter).curDivision;
      splittersMsg->dim[i] = (*iter).curDim;
    }
    /*for(int i=0;i<orbData.size();i++){
      splittersMsg->pos[i] = orbData[i].curDivision;
      splittersMsg->dim[i] = (char) orbData[i].curDim;
    }*/
    treeProxy.evaluateParticleCounts(splittersMsg);
  }

}

void Sorter::readytoSendORB(CkReductionMsg* m){
  delete m;

  treeProxy.sendORBParticles();
}

void Sorter::collectORBCounts(CkReductionMsg* m){

  std::list<ORBData>::iterator iter;
  int i;
  
  numCounts = m->getSize() / sizeof(int);
	binCounts.resize(numCounts);
	//binCounts[0] = 0;
	int* startCounts = static_cast<int *>(m->getData());
	//copy(startCounts, startCounts + numCounts, binCounts.begin() + 1);
	copy(startCounts, startCounts + numCounts, binCounts.begin());
	delete m;

  CkAssert(numCounts == 2*orbData.size());
  
  ORBSplittersMsg *splittersMsg;
  float TOLER=0.05;
  int doneCount=0;
  
  for(i=0,iter=orbData.begin(); iter!=orbData.end(); i++,iter++){
    if(binCounts[2*i+1]*(1-TOLER)<=binCounts[2*i] && binCounts[2*i]<=(1+TOLER)*binCounts[2*i+1]){
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
    /*for(int i=0;i<orbData.size();i++){
      splittersMsg->pos[i] = orbData[i].curDivision;
      splittersMsg->dim[i] = (char) orbData[i].curDim;
    }*/
    //finalize the boundaries in all the Treepieces
    treeProxy.finalizeBoundaries(splittersMsg);
  }
  else{
    splittersMsg = new (orbData.size(),orbData.size()) ORBSplittersMsg(orbData.size(),CkCallback(CkIndex_Sorter::collectORBCounts(0), thishandle));
    for(i=0,iter=orbData.begin();iter!=orbData.end();i++,iter++){
      splittersMsg->pos[i] = (*iter).curDivision;
      splittersMsg->dim[i] = (*iter).curDim;
    }
    /*for(int i=0;i<orbData.size();i++){
      splittersMsg->pos[i] = orbData[i].curDivision;
      splittersMsg->dim[i] = (char) orbData[i].curDim;
    }*/
    treeProxy.evaluateParticleCounts(splittersMsg);
  }
  
}

/*void Sorter::sendBoundingBoxes(CkReductionMsg* m){
  delete m;

  std::list<ORBData>::iterator iter;
  int i;
  
  BoundingBoxes *bounding = new (numChares) BoundingBoxes();

  for(i=0,iter=orbData.begin();iter!=orbData.end();iter++,i++){
    bounding->boxes[i] = (*iter).boundingBox;
  }

  treeProxy.receiveBoundingBoxes(bounding);

}*/

/************************************************/

/*
class WeighedNode {
public:
  Key key;
  u_int64_t count;
  WeighedNode children[2];
  WeighedNode(Key k, u_int64_t c) : key(k), count(c) { }
};
*/

void Sorter::startSorting(const CkGroupID& dataManagerID,
			  const double toler, const CkCallback& cb, bool decompose) {
	numChares = numTreePieces;
	dm = CProxy_DataManager(dataManagerID);
	tolerance = toler;
	sorted = false;
	sortingCallback = cb;
	numIterations = 0;
	
  //Changed for implementing OCT decomposition
  Key delta;
  Key k;
  BinaryTreeNode *rt;
  switch (domainDecomposition){
    case SFC_dec:
    case SFC_peano_dec:
    case SFC_peano_dec_3D:
    case SFC_peano_dec_2D:
        numKeys = 0;
        if (splitters.size() == 0) {
          // reuse the existing splitters from the previous decomposition
	    splitters.clear();
	    splitters.reserve(3 * numChares - 1);
	    delta = (lastPossibleKey - SFC::firstPossibleKey) / (3 * numChares - 2);
	    k = firstPossibleKey;
	    for(int i = 0; i < (3 * numChares - 2); i++, k += delta) {
		if(k != firstPossibleKey)
		    k |= 7L;  // Set bottom bits to avoid trees too deep
		splitters.push_back(k);
		}
	    splitters.push_back(lastPossibleKey);
        }
        break;
    case Oct_dec:
      if (nodeKeys.size() == 0) {
        rt = new BinaryTreeNode();
        //numKeys = numChares;
        //keysSize = (int) (numChares * 1.1);
        nodeKeys.reserve(numChares);
        nodeKeys.resize(numChares>>1, 0);
        chareIDs.resize(nodeKeys.size());
        NodeKey *tmp = &(*nodeKeys.begin());
        rt->getChunks(nodeKeys.size(),tmp);
        delete rt;
        // place the unused chares in the available list
        availableChares.reserve(numChares>>1 + 10);
        for (int i = numChares-1; i>=nodeKeys.size(); --i) {
          availableChares.push_back(i);
        }
        joinThreshold = particlesPerChare;
       splitThreshold = joinThreshold * 1.5;
        /*
        numKeys = 1;
        treeRoot = new WeighedNode(1, 1<<30);
        */
      }
      //wbState = new WeightBalanceState<int>();
      //Convert the Node Keys to the splitter keys which will be sent to histogram
      convertNodesToSplitters();
      break;
    case ORB_dec:
	numKeys = 0;
      treeProxy.initORBPieces(CkCallback(CkIndex_Sorter::doORBDecomposition(0), thishandle));
    break;
      default:
    CkAbort("Invalid domain decomposition requested");
  }

	if(verbosity >= 3)
		cout << "Sorter: Initially have " << splitters.size() << " splitters" << endl;

	//send out the first guesses to be evaluated
  if(domainDecomposition!=ORB_dec){
    if (decompose) {

	// XXX: Optimizations available if sort has been done before!
	//create initial evenly distributed guesses for splitter keys
	keyBoundaries.clear();
	keyBoundaries.reserve(numChares + 1);
	keyBoundaries.push_back(firstPossibleKey);

//      dm.acceptCandidateKeys(&(*splitters.begin()), splitters.size(), 0, CkCallback(CkIndex_Sorter::collectEvaluations(0), thishandle));
	  treeProxy.evaluateBoundaries(&(*splitters.begin()), splitters.size(), 0, CkCallback(CkIndex_Sorter::collectEvaluations(0), thishandle));
    } else {
      //send out all the decided keys to get final bin counts
      sorted = true;
//      dm.acceptCandidateKeys(&(*keyBoundaries.begin()), keyBoundaries.size(), 0, CkCallback(CkIndex_Sorter::collectEvaluations(0), thishandle));
      treeProxy.evaluateBoundaries(&(*keyBoundaries.begin()), keyBoundaries.size(), 0, CkCallback(CkIndex_Sorter::collectEvaluations(0), thishandle));
    }
  }
}

/**
 * Given "numKeys" node keys ("nodeKeys"), convert these keys into splitters.
 */
void Sorter::convertNodesToSplitters(){
  Key partKey;

  splitters.clear();
  splitters.reserve(nodeKeys.size() + 1);
  binCounts.reserve(nodeKeys.size());
  const Key mask = Key(1) << 63;
  for(unsigned int i=0;i<nodeKeys.size();i++){
    partKey=Key(nodeKeys[i]);
    while(!(partKey & mask)){
      partKey <<= 1;
    }
    partKey &= ~mask;
    splitters.push_back(partKey);
  }
  //Sort here to make sure that splitters go sorted to histogramming
  //They might be unsorted here due to sorted or unsorted node keys
  // FILIPPO: no, by construction they must be ordered already!
  //sort(splitters.begin(),splitters.end());
  splitters.push_back(lastPossibleKey);
}

/**
 * Given "num" node keys, create splitters for these selected nodekeys to 
 * refine the information within those nodes: each node is divided
 * "level" times. Returns a newly allocated Key array containing the
 * splitters key for the histogramming phase.
 */
Key * Sorter::convertNodesToSplittersRefine(int num, NodeKey* keys){
  Key partKey = Key(0);

  Key *result = new Key[num * ((2<<refineLevel)+1)];
  int64_t levelMask = int64_t(1) << 63;
  levelMask >>= refineLevel;
  int idx = 0;
  const Key mask = Key(1) << 63;
  for(unsigned int i=0;i<num;i++){
    CkAssert(! (partKey & levelMask));
    partKey=Key(keys[i]<<refineLevel);
    int shift = 0;
    // find how much we need to shift each key (depend on the tree level of the key)
    while (!(partKey<<shift & mask)) {
      ++shift;
    }
    partKey &= ~mask >> shift;
    for (int j=0; j<=(1<<refineLevel); ++j) {
      result[idx++] = (partKey+j) << shift;
    }
  }
  //Sort here to make sure that splitters go sorted to histogramming
  //They might be unsorted here due to sorted or unsorted node keys
  // FILIPPO: no, by construction they must be ordered already!
  //sort(splitters.begin(),splitters.end());
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
      CkAbort("ORB: We shouldn't have reached here");
    break;
      default:
    CkAbort("Invalid domain decomposition requested");
  }
}

void Sorter::collectEvaluationsOct(CkReductionMsg* m) {

  numIterations++;
  numCounts = m->getSize() / sizeof(int);
  //binCounts.resize(numCounts);
  //binCounts[0] = 0;
  int* startCounts = static_cast<int *>(m->getData());
  //copy(startCounts, startCounts + numCounts, binCounts.begin() + 1);
  //copy(startCounts, startCounts + numCounts, binCounts.begin());

  //CkAssert(numCounts==numChares);
  //call function which will balance the bin counts: define it in GenericTreeNode
  //make it a templated function
  //Pass the bincounts as well as the nodekeys
  
  if (joinThreshold == 0) {
    int total_particles = std::accumulate(startCounts, startCounts+numCounts, 0);
    joinThreshold = total_particles / (numTreePieces>>1);
    splitThreshold = joinThreshold * 1.5;
  }
  
  if(verbosity>=3){
    int i=0;
    CkPrintf("Bin Counts in collect eval (%d):",numCounts);
    for ( ; i<numCounts; i++) {
      CkPrintf("%d,",startCounts[i]);
    }
    CkPrintf("\n");
    CkPrintf("Nodekeys:");
    for(int j=0;j<numChares;j++)
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
  //bool histogram=weightBalance<int>(nodeKeys,startCounts,numKeys,keysSize,numChares,&zeros,&nodesOpened,wbState);
  bool histogram = refineOctSplitting(numCounts, startCounts);
  traceUserBracketEvent(weightBalanceUE, startTimer, CmiWallTimer());
  delete m;

  if(verbosity>=3){
    CkPrintf("Nodekeys after (%d):",nodeKeys.size());
    for(int i=0;i<nodeKeys.size();i++)
      CkPrintf("%llx,",nodeKeys[i]);
    CkPrintf("\n");
  }
  if(histogram){
    //convertNodesToSplitters(numKeys,nodeKeys);
    refineLevel = 2;
    int arraySize = (1<<refineLevel)+1;
    Key *array = convertNodesToSplittersRefine(nodesOpened.size(),nodesOpened.getVec());
//    dm.acceptCandidateKeys(&(*splitters.begin()), splitters.size(), 1, CkCallback(CkIndex_Sorter::collectEvaluations(0), thishandle));
    treeProxy.evaluateBoundaries(array, nodesOpened.size()*arraySize, 1<<refineLevel, CkCallback(CkIndex_Sorter::collectEvaluations(0), thishandle));
    delete[] array;
  }
  else{
    sorted=true;
    if(verbosity)
      ckout << "Sorter: Histograms balanced after " << numIterations
	    << " iterations." << endl;
    //We have the splitters here because weight balancer didn't change any node keys
    //We also have the final bin counts
		
    //determine which TreePiece is responsible for each interval
    //vector<int> chareIDs(numChares, 1);
    //chareIDs[0] = 0;
    //partial_sum(chareIDs.begin(), chareIDs.end(), chareIDs.begin());
	
    //send out the final splitters and responsibility table
    //convertNodesToSplittersNoZeros(numKeys,nodeKeys,zeros);
    //convertNodesToSplitters(); // Filippo: not needed anymore since splitters is kept in sync with
                                 // the other arrays by the function refileOctSplitting.
    dm.acceptFinalKeys(&(*splitters.begin()), &(*chareIDs.begin()), &(*binCounts.begin()), splitters.size(), sortingCallback);
    numIterations = 0;
    sorted = false;
    //delete [] nodeKeys;
    //delete wbState;
    return;
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
  int i, idx;
  if (nodesOpened.size() == 0) {
    // This means that we are not refining the counts, but we got a brand new histogramming     
    CkAssert(n == nodeKeys.size());
    binCounts.resize(n);
    copy(count, count+n, binCounts.begin());
    //CkPrintf("Sorter: threshold %d %d\n",splitThreshold,joinThreshold);
    //CkPrintf("Sorter: incoming> {");
    //for (i=0; i<n; ++i) CkPrintf(" (%llx)%d",nodeKeys[i],binCounts[i]);
    //CkPrintf(" }\n");
    // Walk the keys and assign the new counts, while walking also decide if some nodes
    // are to be opened further or joined together
    for (i=0, idx=0; i<n; ++i, ++idx) {
      // Check if the node has too many particles and needs to be split
      if (binCounts[idx] > splitThreshold) {
        nodesOpened.push_back(nodeKeys[idx]);
        //CkPrintf("Sorter: opening %llx (%d)\n",nodeKeys[idx],binCounts[idx]);
      }
      // Check if two nodes can be joined together
      while (idx>0 && (binCounts[idx-1]+binCounts[idx] < joinThreshold) && (nodeKeys[idx-1]>>1 == nodeKeys[idx]>>1)) {
        // Join and repeat the check recursively
        //CkPrintf("Sorter: joining %llx and %llx (%d + %d)\n",nodeKeys[idx-1],nodeKeys[idx],binCounts[idx-1],binCounts[idx]);
        nodeKeys[idx-1] >>= 1;
        nodeKeys.erase(nodeKeys.begin()+idx);
        splitters.erase(splitters.begin()+idx);
        binCounts[idx-1] += binCounts[idx];
        binCounts.erase(binCounts.begin()+idx);
        availableChares.push_back(chareIDs[idx]);
        chareIDs.erase(chareIDs.begin()+idx);
        --idx;
      }
    }
  } else {
    CkVec<NodeKey> tmpOpened;
    int64_t levelMask = int64_t(1) << 63;
    levelMask >>= refineLevel;
    const Key mask = Key(1) << 63;
    for (i=0; i<nodesOpened.size(); ++i) {
      //CkPrintf("Sorter: considering %llx\n",nodesOpened[i]);
      if (availableChares.size() < 1<<refineLevel) {
	CkPrintf("availableChares size is %d, cannot refine further\n", availableChares.size());
        break;
      }
      NodeKey key = nodesOpened[i];
      int index = std::find(nodeKeys.begin(), nodeKeys.end(), key) - nodeKeys.begin();

      CkAssert(! (key & levelMask));
      key <<= refineLevel;
      int shift = 0;
      // find how much we need to shift each key (depend on the tree level of the key)
      while (!(key<<shift & mask)) {
        ++shift;
      }

      // Add all the nodes we just requested for histogram
      nodeKeys[index] <<= refineLevel;
      nodeKeys.insert(nodeKeys.begin()+index+1, (1<<refineLevel)-1, key);
      splitters.insert(splitters.begin()+index+1, (1<<refineLevel)-1, key&(~mask >> shift));
      for (int j=1; j<(1<<refineLevel); ++j) {
        nodeKeys[index+j] += j;
        splitters[index+j] = (splitters[index+j] + j) << shift;
      }
      binCounts[index] = count[i*(1<<refineLevel)];
      binCounts.insert(binCounts.begin()+index+1, &count[i*(1<<refineLevel)+1], &count[(i+1)*(1<<refineLevel)]);
      chareIDs.insert(chareIDs.begin()+index+1, availableChares.end()-(1<<refineLevel)+1, availableChares.end());

      if (verbosity >= 4 ) 
	CkPrintf("Split node index %d, last added chare is %d (refine level = %d), %d available chares left\n", index, availableChares.back(), (1<<refineLevel), availableChares.size()-1);

      
      availableChares.erase(availableChares.end()-(1<<refineLevel)+1, availableChares.end());
     
 
      // Trim down what we over-opened just above
      for (int j=1, idx=0; j<=(1<<refineLevel); ++j, ++idx) {
        if (binCounts[index+idx] > splitThreshold) {
          //CkPrintf("Sorter: further opening %llx (%d)\n",nodeKeys[index+idx],binCounts[index+idx]);
          tmpOpened.push_back(nodeKeys[index+idx]);
        }
        while (idx>0 && (binCounts[index+idx-1]+binCounts[index+idx] <= splitThreshold) && (nodeKeys[index+idx-1]>>1 == nodeKeys[index+idx]>>1)) {
          // Join and repeat the check recursively
          //CkPrintf("Sorter: re-joining %llx and %llx (%d + %d)\n",nodeKeys[index+idx-1],nodeKeys[index+idx],binCounts[index+idx-1],binCounts[index+idx]);
          nodeKeys[index+idx-1] >>= 1;
          nodeKeys.erase(nodeKeys.begin()+index+idx);
          splitters.erase(splitters.begin()+index+idx);
          binCounts[index+idx-1] += binCounts[index+idx];
          binCounts.erase(binCounts.begin()+index+idx);
          availableChares.push_back(chareIDs[index+idx]);
          chareIDs.erase(chareIDs.begin()+index+idx);
          --idx;
        }
      }
    }
    // Copy the new list of nodes to be opened to nodesOpened
    nodesOpened.removeAll();
    for (i=0; i<tmpOpened.size(); ++i) {
      nodesOpened.push_back(tmpOpened[i]);
    }
  }
  //CkPrintf("Sorter: final situation {");
  //for (i=0; i<nodeKeys.size(); ++i) CkPrintf(" %llx(%d)",nodeKeys[i],binCounts[i]);
  //CkPrintf(" }\n");
  return nodesOpened.size() > 0;
}

void Sorter::collectEvaluationsSFC(CkReductionMsg* m) {
	numIterations++;
	numCounts = m->getSize() / sizeof(int);
	binCounts.resize(numCounts + 1);
	binCounts[0] = 0;
	int* startCounts = static_cast<int *>(m->getData());
	copy(startCounts, startCounts + numCounts, binCounts.begin() + 1);
	delete m;
	
	if(sorted) { //true after the final keys have been binned
		//determine which TreePiece is responsible for each interval
		//vector<int> chareIDs(numChares, 1);
		//chareIDs[0] = 0;
		//partial_sum(chareIDs.begin(), chareIDs.end(), chareIDs.begin());
		
		//send out the final splitters and responsibility table
		dm.acceptFinalKeys(&(*keyBoundaries.begin()), &(*chareIDs.begin()), &(*binCounts.begin()) + 1, keyBoundaries.size(), sortingCallback);
		numIterations = 0;
		sorted = false;
		return;
	}
	
	if(verbosity >= 4)
		cout << "Sorter: On iteration " << numIterations << endl;
	CkAssert(numIterations < 1000);  // Sorter has not converged.
	
	//sum up the individual bin counts, so each bin has the count of it and all preceding
	partial_sum(binCounts.begin(), binCounts.end(), binCounts.begin());
	
	if(!numKeys) {
		numKeys = binCounts.back();
		int avgValue = numKeys / numChares;
		closeEnough = static_cast<int>(avgValue * tolerance);
		if(closeEnough < 0 || closeEnough >= avgValue) {
			cerr << "Sorter: Unacceptable tolerance, requiring exact fit." << endl;
			closeEnough = 0;
		}
		
		//each splitter key will split the keys near a goal number of keys
		goals.assign(numChares - 1, avgValue);
		partial_sum(goals.begin(), goals.end(), goals.begin());
		
		if(verbosity >= 3)
			cout << "Sorter: Target keys per chare: " << avgValue << " plus/minus " << (2 * closeEnough) << endl;
	}

	//make adjustments to the splitter keys based on the results of the previous iteration
	adjustSplitters();

	if(verbosity >= 4) {
		cout << "Sorter: Probing " << splitters.size() << " splitter keys" << endl;
		cout << "Sorter: Decided on " << (keyBoundaries.size() - 1) << " splitting keys" << endl;
	}
	
	//check if we have found all the splitters
	if(sorted) {
		if(verbosity)
			cerr << "Sorter: Histograms balanced after " << numIterations << " iterations." << endl;
		
		sort(keyBoundaries.begin() + 1, keyBoundaries.end());
		keyBoundaries.push_back(lastPossibleKey);
		
		//send out all the decided keys to get final bin counts
//		dm.acceptCandidateKeys(&(*keyBoundaries.begin()), keyBoundaries.size(), 0, CkCallback(CkIndex_Sorter::collectEvaluations(0), thishandle));
        treeProxy.evaluateBoundaries(&(*keyBoundaries.begin()), keyBoundaries.size(), 0, CkCallback(CkIndex_Sorter::collectEvaluations(0), thishandle));
	} else //send out the new guesses to be evaluated
//		dm.acceptCandidateKeys(&(*splitters.begin()), splitters.size(), 0, CkCallback(CkIndex_Sorter::collectEvaluations(0), thishandle));
	    treeProxy.evaluateBoundaries(&(*splitters.begin()), splitters.size(), 0, CkCallback(CkIndex_Sorter::collectEvaluations(0), thishandle));
}

/** Generate new guesses for splitter keys based on the histograms that came
 back from the last batch.
 We need to find the keys that split a distribution into well defined piles.
 We send out low, high, and middle guesses for each split.  We then pick the
 left or right side and move into there, sending out for evaluation.  This
 is a simultaneous binary search for each splitter key not yet found.
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
			cerr << "Sorter: Looking for " << *Ngoal << " How could this happen at the beginning?" << endl;
		if(numRightKey == binCounts.end())
			cerr << "Sorter: Looking for " << *Ngoal << " How could this happen at the end?" << endl;
		
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
			newSplitters.insert((leftBound / 2 + rightBound / 2)
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
		for(set<Key>::iterator iterNew = newSplitters.begin();
		    iterNew != newSplitters.end(); iterNew++) {
		    splitters.push_back(*iterNew);
		    }
	}
}
