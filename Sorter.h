/** \file Sorter.h
 This class manages a parallel sort of particles using a DataManager to hold
 suggested splitter keys.
 \author Graeme Lufkin (gwl@u.washington.edu)
*/

#ifndef SORTER_H
#define SORTER_H

#include <list>
#include <vector>
#include <set>
#include <numeric>

#include "ParallelGravity.h"

/**
   \brief Domain decomposition of particles via a parallel sort
*/
/** Domain decomposition is handled differently depending on the type
 of decomposition requested.  SFC (either Morton ordering or
 Peano-Hibert) or Oct (where domains are complete nodes within the
 tree) is described first.
 
 In the SFC and Oct case, the Sorter controls a parallel sorting of particles
 based on their keys.  The goal is to split the particles up into n
 pieces of roughly equal sizes.
 
 The Sorter generates a list of possible splitter keys, and broadcasts this
 list to the TreePieces. The TreePieces evaluate the list, determining
 how many particles they have between each splitter key.  These
 bin counts are then reduced back to the Sorter::collectEvaluations,
 which uses the information to create a new list.  When the
 last iteration of splitter keys is deemed satisfactory, counts for
 the final keys are evaluated, and the keys and counts are broadcast
 to DataManager::acceptFinalKeys(), which coordinates the
 shuffling of the particles to the correct owners.
 */
class Sorter : public Chare {
	/// The total number of keys we're sorting.
	int numKeys;
	/// The number of chares to sort into.
	int numChares;
	// The number of chares currently with assigned data.
	//int numUsedChares;
	/// The indexes of the chares that are responsible for each segment of data.
	std::vector<int> chareIDs;
	/// A list of chare elements to which nothing is assigned
	std::vector<int> availableChares;
	// Total size of the keys allocated (allows a margin to increase)
	//int keysSize;

	/// The percent tolerance to sort keys within.
	double tolerance;
	/// The number of particles on either side of a splitter that corresponds to the requested tolerance.
	int closeEnough;
	/// The number of iterations completed.
	int numIterations;
	/// A flag telling if we're done yet.
	bool sorted;

	std::vector<NodeKey> nodeKeys;
	/// The histogram of counts for the last round of splitter keys.
	std::vector<unsigned int> binCounts;
	std::vector<unsigned int> binCountsGas;
	/// The number of bins in the histogram.
	int numCounts;
	/// The keys I've decided on that divide the objects evenly (within the tolerance).
	std::vector<SFC::Key> keyBoundaries;
	/// The keys I'm sending out to be evaluated.
	std::vector<SFC::Key> splitters;
	/// The list of object number splits not yet met.
	std::list<int> goals;
	
	/// The DataManager I broadcast candidate keys to.
	CProxy_DataManager dm;
	/// The callback for when the sort is complete.
	CkCallback sortingCallback;

	/// Variables to decide when to split or join a TreePiece in the Oct decomposition
	int joinThreshold, splitThreshold;
	/// Specify what is the level of refinement of nodes sent out
	/// for histogramming in Oct decomposition.
	int refineLevel;
	
  ///Variables added for ORB decomposition
  typedef struct DivData{
    OrientedBox<float> boundingBox;
    double curLow;
    double curHigh;
    double curDivision;
    char curDim;
  } ORBData;

  /// Contains the split information for the domain decomposition.
  /// Splits are inserted into the list as the decomposition tree is
  /// built down from the top.
  std::list<ORBData> orbData;

	
        /// The weights for all the keys returned by the weightBalance routine
        CkVec<int> zeros;
        /// The list of nodes opened by the last invocation of weightBalance
        CkVec<NodeKey> nodesOpened;
        /// The transient state used by the weightBalance routine
        WeightBalanceState<int>* wbState;

  Compare comp;

	void adjustSplitters();
	bool refineOctSplitting(int n, int *count);
	
public:
	
	Sorter() {
          chareIDs.resize(numTreePieces, 1);
          chareIDs[0] = 0;
          partial_sum(chareIDs.begin(), chareIDs.end(), chareIDs.begin());
	};
	Sorter(CkMigrateMessage* m) {
          chareIDs.resize(numTreePieces, 1);
          chareIDs[0] = 0;
          partial_sum(chareIDs.begin(), chareIDs.end(), chareIDs.begin());
	};

	/** Sort the particles in an array of TreePieces using a histogram-probing method.
	 The the Sorter sends splitters to the TreePieces for them to evaluate.
	 The evaluation results in a histogram back to the Sorter,
	 which uses the information to generate new guesses for the splitter keys.
	 When all the splitter keys have been found (within the percent tolerance per chare)
	 the TreePieces move the data around to the proper sorted location.  The boundary
	 keys are shared between adjacent TreePieces, as these will be needed in the
	 tree building phase.
	 The callback will receive a CkReductionMsg containing no data.
	 */
	void startSorting(const CkGroupID& dataManagerID, const double toler, const CkCallback& cb, bool decompose);
	void convertNodesToSplitters();
	SFC::Key * convertNodesToSplittersRefine(int num, NodeKey* keys);
	//void convertNodesToSplittersNoZeros(int num, NodeKey* nodeKeys, CkVec<int> &zero);
	void collectEvaluations(CkReductionMsg* m);
	void collectEvaluationsSFC(CkReductionMsg* m);
	void collectEvaluationsOct(CkReductionMsg* m);

  //ORB Decomposition
  void doORBDecomposition(CkReductionMsg* m);
  void finishPhase(CkReductionMsg *m);
  void collectORBCounts(CkReductionMsg* m);
  void readytoSendORB(CkReductionMsg* m);
  //void sendBoundingBoxes(CkReductionMsg* m);
};

#endif //SORTER_H
