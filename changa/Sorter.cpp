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
      //NodeKey *nodeKeys;
      nodeKeys = new NodeKey[numChares];
      rt->getChunks(numChares,nodeKeys);
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

/*
class compare{
  public:
    compare(){}
    
    bool operator()(NodeKey key1, NodeKey key2){

      NodeKey tmp = NodeKey(1);
      int len1=0, len2=0;
      int cnt=1;
  
      while(tmp<=key1 || tmp<=key2){
        tmp<<=1;
        cnt++;
        if(len1==0 && tmp>key1){
          len1=cnt-1;
        }
        if(len2==0 && tmp>key2){
          len2=cnt-1;
        }
      }
  
      if(len1==len2){
        return key1<key2;
      }
      else if(len1>len2){
        key1>>=(len1-len2);
        return key1<key2;
      }
      else{
        key2>>=(len2-len1);
        return key1<key2;
      }
    }
};

template <class T>
inline bool weightBalance(NodeKey *nodeKeys, T* weights, int num){

  //T can be signed or unsigned
  
	NodeKey curHeaviest;
	map<NodeKey,T,compare>::iterator curLightest;
	T lightestWt= ~T(0);
	T tmpWt=0;

	NodeKey parent,child1,child2;
  int numBalances=0;
  
	//Need to construct a temporary copy of the input data to operate
	//construct a map indexed by the nodekey
	map<NodeKey,T,compare> curNodeWts;
	map<NodeKey,T,compare>::iterator iter;
	map<NodeKey,T,compare>::iterator iter2;
	curNodeWts.clear();
	for(int i=0;i<num;i++){
		curNodeWts[nodeKeys[i]]=weights[i];
	}

	//loop here
  while(1){
    tmpWt=0;
    lightestWt=~T(0);
	  //find the heaviest Node
	  for(iter=curNodeWts.begin();iter!=curNodeWts.end();iter++){
		  if((*iter).second>tmpWt){
			  tmpWt=(*iter).second;
			  curHeaviest=(*iter).first;
		  }
	  }
    if(tmpWt==0) //In case, no-one had weight > 0
      break;
    
	  //find the lightest parent-- implemented only for a binary tree
    iter=curNodeWts.begin();
    iter2=curNodeWts.begin();
    iter2++;
	  for( ;iter2!=curNodeWts.end();iter++,iter2++){
		  if((*iter).second==~T(0) || (*iter2).second==~T(0))//Ignore those which have been opened
        continue;
      if((*iter).first==curHeaviest || (*iter2).first==curHeaviest)
        continue;
      child1=(*iter).first;
		  child2=(*(iter2)).first;
		  child1 >>= 1;
		  child2 >>= 1;
		  if(child1==child2){
			  tmpWt=(*iter).second+(*iter2).second;
			  if(tmpWt<lightestWt || lightestWt==~T(0)){
				  lightestWt=tmpWt;
				  curLightest=iter;
			  }
      }
	  }

	  if((curNodeWts[curHeaviest] > lightestWt) && lightestWt!=~T(0)){
		  numBalances++;
      parent = (*curLightest).first >> 1;
      iter2=curLightest; iter2++; iter2++;
		  //Erase the children and add the lightest parent
      curNodeWts.erase(curLightest,iter2);
		  //curNodeWts[parent]=lightestWt;
      curNodeWts.insert(pair<NodeKey,T>(parent,lightestWt));
      child1=curHeaviest << 1;
      child2=child1 | NodeKey(1);
      //Erase the heaviest and add it's two children
      curNodeWts.erase(curHeaviest);
      curNodeWts[child1]=~T(0);
      curNodeWts[child2]=~T(0);
	  }
	  else //We are done here
		  break;
  }
  //end loop here

  int i=0;
  if(numBalances>0){
    //construct new node key array before returning
	  for(iter=curNodeWts.begin(),i=0;iter!=curNodeWts.end();i++,iter++){
      nodeKeys[i]=(*iter).first;
    }
    CkAssert(i==num);
    return true;
  }
  else { return false; }

}*/

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
    std::vector<int>::iterator iter;
    CkPrintf("Bin Counts in collect eval:");
    for(iter=binCounts.begin();iter!=binCounts.end();iter++){
      CkPrintf("%d,",*iter);
    }
    CkPrintf("\n");
    CkPrintf("nodekeys:");
    for(int i=0;i<numChares;i++)
      CkPrintf("%llx,",nodeKeys[i]);
    CkPrintf("\n");
  }
  
	bool histogram=weightBalance<int>(nodeKeys,&(*binCounts.begin()),numChares);
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
	vector<int>::iterator numLeftKey, numRightKey = binCounts.begin();
	
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
		if(abs(*numLeftKey - *Ngoal) <= closeEnough) {
			//add this key to the list of decided splitter keys
			keyBoundaries.push_back(leftBound);
			//the goal has been met, delete it
			list<int>::iterator temp = Ngoal;
			++Ngoal;
			goals.erase(temp);
		} else if(abs(*numRightKey - *Ngoal) <= closeEnough) {
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
