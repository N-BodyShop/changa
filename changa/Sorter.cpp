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

void Sorter::startSorting(const CkGroupID& dataManagerID, const int nChares,
			  const double toler, const CkCallback& cb, bool decompose) {
	numChares = nChares;
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
	numKeys = 0;
        //if (splitters.size() == 0) {
          // reuse the existing splitters from the previous decomposition
            splitters.clear();
	    splitters.reserve(3 * numChares - 1);
	    delta = (lastPossibleKey - SFC::firstPossibleKey) / (3 * numChares - 2);
	    k = firstPossibleKey;
	    for(int i = 0; i < (3 * numChares - 2); i++, k += delta)
		    splitters.push_back(k);
	    splitters.push_back(lastPossibleKey);
            //}
        break;
    case Oct_dec:
      if (nodeKeys == NULL) {
        rt = new BinaryTreeNode();
        numKeys = numChares;
        keysSize = (int) (numChares * 1.1);
        nodeKeys = new NodeKey[keysSize];
        rt->getChunks(numKeys,nodeKeys);
        delete rt;
      }
      wbState = new WeightBalanceState<int>();
      //Convert the Node Keys to the splitter keys which will be sent to histogram
      convertNodesToSplitters(numKeys,nodeKeys);
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

      dm.acceptCandidateKeys(&(*splitters.begin()), splitters.size(), 0, CkCallback(CkIndex_Sorter::collectEvaluations(0), thishandle));
    } else {
      //send out all the decided keys to get final bin counts
      sorted = true;
      dm.acceptCandidateKeys(&(*keyBoundaries.begin()), keyBoundaries.size(), 0, CkCallback(CkIndex_Sorter::collectEvaluations(0), thishandle));
    }
  }
}

void Sorter::convertNodesToSplitters(int num, NodeKey* nodeKeys){
  Key partKey;

  splitters.clear();
  splitters.reserve(num + 1);
  const Key mask = Key(1) << 63;
  for(unsigned int i=0;i<num;i++){
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

void Sorter::convertNodesToSplittersRefine(int num, NodeKey* nodeKeys){
  Key partKey;
  Key nextKey;
  Key lastKey;

  splitters.clear();
  splitters.reserve(num * 4);
  const Key mask = Key(1) << 63;
  for(unsigned int i=0;i<num;i++){
    partKey=Key(nodeKeys[i]<<1);
    nextKey=partKey+1;
    lastKey=nextKey+1;
    while(!(partKey & mask)){
      partKey <<= 1;
      nextKey <<= 1;
      lastKey <<= 1;
    }
    partKey &= ~mask;
    nextKey &= ~mask;
    lastKey &= ~mask;
    splitters.push_back(partKey);
    splitters.push_back(nextKey);
    splitters.push_back(nextKey);
    splitters.push_back(lastKey);
  }
  //Sort here to make sure that splitters go sorted to histogramming
  //They might be unsorted here due to sorted or unsorted node keys
  // FILIPPO: no, by construction they must be ordered already!
  //sort(splitters.begin(),splitters.end());
}

void Sorter::convertNodesToSplittersNoZeros(int num, NodeKey* nodeKeys, CkVec<int> &zero){
  Key partKey;

  splitters.clear();
  splitters.reserve(num + 1);
  binCounts.clear();
  binCounts.reserve(num + 1);
  const Key mask = Key(1) << 63;
  for(unsigned int i=0;i<num;i++){
    if (zero[i] == 0) continue;
    partKey=Key(nodeKeys[i]);
    while(!(partKey & mask)){
      partKey <<= 1;
    }
    partKey &= ~mask;
    splitters.push_back(partKey);
    binCounts.push_back(zero[i]);
  }
  //Sort here to make sure that splitters go sorted to histogramming
  //They might be unsorted here due to sorted or unsorted node keys
  // FILIPPO: no, by construction they must be ordered already!
  //sort(splitters.begin(),splitters.end());
  splitters.push_back(lastPossibleKey);
}

void Sorter::collectEvaluations(CkReductionMsg* m) {

  switch (domainDecomposition){
    case SFC_dec:
    case SFC_peano_dec:
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
  binCounts.resize(numCounts);
  //binCounts[0] = 0;
  int* startCounts = static_cast<int *>(m->getData());
  //copy(startCounts, startCounts + numCounts, binCounts.begin() + 1);
  //copy(startCounts, startCounts + numCounts, binCounts.begin());

  //CkAssert(numCounts==numChares);
  //call function which will balance the bin counts: define it in GenericTreeNode
  //make it a templated function
  //Pass the bincounts as well as the nodekeys
  
  if(verbosity>=3){
    std::vector<unsigned int>::iterator iter;
    int i=0;
    CkPrintf("Bin Counts in collect eval (%d):",numCounts);
    //for(iter=binCounts.begin();iter!=binCounts.end();iter++,i++){
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
  bool histogram=weightBalance<int>(nodeKeys,startCounts,numKeys,keysSize,numChares,&zeros,&nodesOpened,wbState);
  //bool histogram=weightBalance<int>(nodeKeys,startCounts,numKeys,1);
  traceUserBracketEvent(weightBalanceUE, startTimer, CmiWallTimer());
  delete m;

  if(verbosity>=3){
    CkPrintf("Nodekeys after (%d):",numKeys);
    for(int i=0;i<numKeys;i++)
      CkPrintf("%llx,",nodeKeys[i]);
    CkPrintf("\n");
  }
  if(histogram){
    //convertNodesToSplitters(numKeys,nodeKeys);
    convertNodesToSplittersRefine(nodesOpened.size(),nodesOpened.getVec());
    dm.acceptCandidateKeys(&(*splitters.begin()), splitters.size(), 1, CkCallback(CkIndex_Sorter::collectEvaluations(0), thishandle));
  }
  else{
    sorted=true;
    if(verbosity)
      cerr << "Sorter: Histograms balanced after " << numIterations << " iterations." << endl;
    //We have the splitters here because weight balancer didn't change any node keys
    //We also have the final bin counts
		
    //determine which TreePiece is responsible for each interval
    vector<int> chareIDs(numChares, 1);
    chareIDs[0] = 0;
    partial_sum(chareIDs.begin(), chareIDs.end(), chareIDs.begin());
	
    //send out the final splitters and responsibility table
    convertNodesToSplittersNoZeros(numKeys,nodeKeys,zeros);
    dm.acceptFinalKeys(&(*splitters.begin()), &(*chareIDs.begin()), &(*binCounts.begin()), splitters.size(), sortingCallback);
    numIterations = 0;
    sorted = false;
    //delete [] nodeKeys;
    delete wbState;
    return;
  }
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
		vector<int> chareIDs(numChares, 1);
		chareIDs[0] = 0;
		partial_sum(chareIDs.begin(), chareIDs.end(), chareIDs.begin());
		
		//send out the final splitters and responsibility table
		dm.acceptFinalKeys(&(*keyBoundaries.begin()), &(*chareIDs.begin()), &(*binCounts.begin()) + 1, keyBoundaries.size(), sortingCallback);
		numIterations = 0;
		sorted = false;
		return;
	}
	
	if(verbosity >= 4)
		cout << "Sorter: On iteration " << numIterations << endl;
	
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
		dm.acceptCandidateKeys(&(*keyBoundaries.begin()), keyBoundaries.size(), 0, CkCallback(CkIndex_Sorter::collectEvaluations(0), thishandle));
	} else //send out the new guesses to be evaluated
		dm.acceptCandidateKeys(&(*splitters.begin()), splitters.size(), 0, CkCallback(CkIndex_Sorter::collectEvaluations(0), thishandle));
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
			//not close enough yet, add the bracketing keys and the middle to the guesses
			newSplitters.insert(leftBound);
			newSplitters.insert(leftBound / 2 + rightBound / 2);
			newSplitters.insert(rightBound);
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
