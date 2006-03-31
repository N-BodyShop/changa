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

/** The Sorter will control sorting of particles based on their keys.
 The Sorter generates a list of possible splitter keys, and broadcasts this
 list to the DataManager.  When the DataManager has accepted the list, it
 tells its TreePieces to evaluate the list, determining how many particles
 it has between each splitter key.  These bin counts are then reduced back
 to the Sorter, which uses the information to create a new list.  When the
 last iteration of splitter keys is deemed satisfactory, the final keys
 are broadcast, and the TreePieces can shuffle the particles to the correct
 owners.
 */
class Sorter : public Chare {
	/// The total number of keys we're sorting.
	int numKeys;
	/// The number of chares to sort into.
	int numChares;
	
	/// The percent tolerance to sort keys within.
	double tolerance;
	/// The number of particles on either side of a splitter that corresponds to the requested tolerance.
	int closeEnough;
	/// The number of iterations completed.
	int numIterations;
	/// A flag telling if we're done yet.
	bool sorted;

  NodeKey *nodeKeys;
	/// The histogram of counts for the last round of splitter keys.
	std::vector<int> binCounts;
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
	
	void adjustSplitters();
	
public:
	
	Sorter() { };

	/** Sort the particles in an array of TreePieces using a histogram-probing method.
	 The DataManager receives splitter keys from the Sorter, and instructs the TreePieces
	 to evaluate them.  The evaluation results in a histogram back to the Sorter,
	 which uses the information to generate new guesses for the splitter keys.
	 When all the splitter keys have been found (within the percent tolerance per chare)
	 the TreePieces move the data around to the proper sorted location.  The boundary
	 keys are shared between adjacent TreePieces, as these will be needed in the
	 tree building phase.
	 The callback will receive a CkReductionMsg containing no data.
	 */
	void startSorting(const CkGroupID& dataManagerID, const int nChares, const double toler, const CkCallback& cb);
	void convertNodesToSplitters(int numChares, NodeKey* nodeKeys);
	void collectEvaluations(CkReductionMsg* m);
	void collectEvaluationsSFC(CkReductionMsg* m);
	void collectEvaluationsOct(CkReductionMsg* m);
};

#endif //SORTER_H
