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

void Sorter::startSorting(const CkGroupID& dataManagerID, const int nChares,
			  const double toler, const CkCallback& cb) {	
	numChares = nChares;
	dm = CProxy_DataManager(dataManagerID);
	tolerance = toler;
	sorted = false;
	sortingCallback = cb;
	numIterations = 0;
	numKeys = 0;
	
	// XXX: Optimizations available if sort has been done before!
	//create initial evenly distributed guesses for splitter keys
	keyBoundaries.clear();
	keyBoundaries.reserve(numChares + 1);
	keyBoundaries.push_back(firstPossibleKey);

  //Changed for implementing OCT decomposition
  Key delta;
  Key k;
  BinaryTreeNode *rt;
  switch (domainDecomposition){
    case SFC_dec:
	    splitters.clear();
	    splitters.reserve(3 * numChares - 1);
	    delta = (lastPossibleKey - SFC::firstPossibleKey) / (3 * numChares - 2);
	    k = firstPossibleKey;
	    for(int i = 0; i < (3 * numChares - 2); i++, k += delta)
		    splitters.push_back(k);
	    splitters.push_back(lastPossibleKey);
	    break;
    case Oct_dec:
	    rt = new BinaryTreeNode();
      nodeKeys = new NodeKey[numChares];
      rt->getChunks(numChares,nodeKeys);
      delete rt;
      //Convert the Node Keys to the splitter keys which will be sent to histogram
      convertNodesToSplitters(numChares,nodeKeys);
	    break;
    case ORB_dec:
      CkAbort("ORB domain decomposition not yet implemented");
    break;
      default:
    CkAbort("Invalid domain decomposition requested");
  }

	if(verbosity >= 3)
		cout << "Sorter: Initially have " << splitters.size() << " splitters" << endl;

	//send out the first guesses to be evaluated
	dm.acceptCandidateKeys(&(*splitters.begin()), splitters.size(), CkCallback(CkIndex_Sorter::collectEvaluations(0), thishandle));
}

void Sorter::convertNodesToSplitters(int numChares, NodeKey* nodeKeys){
  Key partKey;

  splitters.clear();
  splitters.reserve(numChares + 1);
  const Key mask = Key(1) << 63;
  for(unsigned int i=0;i<numChares;i++){
    partKey=Key(nodeKeys[i]);
    while(!(partKey & mask)){
      partKey <<= 1;
    }
    partKey &= ~mask;
    splitters.push_back(partKey);
  }
  //Sort here to make sure that splitters go sorted to histogramming
  //They might be unsorted here due to sorted or unsorted node keys
  sort(splitters.begin(),splitters.end());
  splitters.push_back(lastPossibleKey);
}

void Sorter::collectEvaluations(CkReductionMsg* m) {

  switch (domainDecomposition){
    case SFC_dec:
      collectEvaluationsSFC(m);
      break;
    case Oct_dec:
      collectEvaluationsOct(m);
      break;
    case ORB_dec:
      CkAbort("ORB domain decomposition not yet implemented");
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
	copy(startCounts, startCounts + numCounts, binCounts.begin());
	delete m;

	CkAssert(numCounts==numChares);
  //call function which will balance the bin counts: define it in GenericTreeNode
  //make it a templated function
  //Pass the bincounts as well as the nodekeys
  
  if(verbosity>=3){
    std::vector<unsigned int>::iterator iter;
    int i=0;
    CkPrintf("Bin Counts in collect eval:");
    for(iter=binCounts.begin();iter!=binCounts.end();iter++){
      CkPrintf("%u,",*iter);
      i++;
    }
    CkPrintf("\n");
    CkPrintf("Nodekeys:");
    for(int j=0;j<numChares;j++)
      CkPrintf("%llx,",nodeKeys[j]);
    CkPrintf("\n");
  }
  
	bool histogram=weightBalance<unsigned int>(nodeKeys,&(*binCounts.begin()),numChares,1);
  if(verbosity>=3){
    CkPrintf("Nodekeys after:");
    for(int i=0;i<numChares;i++)
      CkPrintf("%llx,",nodeKeys[i]);
    CkPrintf("\n");
  }
	if(histogram){
    convertNodesToSplitters(numChares,nodeKeys);
    dm.acceptCandidateKeys(&(*splitters.begin()), splitters.size(), CkCallback(CkIndex_Sorter::collectEvaluations(0), thishandle));
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
		dm.acceptFinalKeys(&(*splitters.begin()), &(*chareIDs.begin()), &(*binCounts.begin()), splitters.size(), sortingCallback);
		numIterations = 0;
		sorted = false;
    delete [] nodeKeys;
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
		dm.acceptCandidateKeys(&(*keyBoundaries.begin()), keyBoundaries.size(), CkCallback(CkIndex_Sorter::collectEvaluations(0), thishandle));
	} else //send out the new guesses to be evaluated
		dm.acceptCandidateKeys(&(*splitters.begin()), splitters.size(), CkCallback(CkIndex_Sorter::collectEvaluations(0), thishandle));
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
		splitters.assign(newSplitters.begin(), newSplitters.end());
	}
}
