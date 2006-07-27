/** @file DataManager.h
 The DataManager holds data that TreePieces need access to, but not
 their own copy.
 @author Graeme Lufkin (gwl@u.washington.edu)
*/

#ifndef DATAMANAGER_H
#define DATAMANAGER_H

#include <vector>
#include <map>
#include <string>
#include "CacheManager.h"

/** The DataManager is used to store information that all TreePieces will need,
 but will not modify.  The first example is the list of splitter keys and the
 responsible chare for each interval.  This data is given to the DataManager by
 the Sorter.  The DataManager then instructs all the TreePieces on its node
 to evaluate the boundary keys.
 */
class DataManager : public CBase_DataManager {
	friend class TreePiece;
	/// The array of TreePieces I hold data for.
	CProxy_TreePiece treePieces;
	/// A pointer to the CacheManager in the local processor
	CacheManager *localCache;

protected:
	
	// An array telling me which TreePieces are actually on my node.
	// USE THE INFO IN THE CACHEMANAGER WHICH IS KEPT UPDATED
	//std::vector<int> myTreePieces;
	
	/// The array of splitter keys for the sort.
	std::vector<SFC::Key> boundaryKeys;
	
	/// An array identifying which chare is responsible for each interval of keys.
	std::vector<int> responsibleIndex;
	
	/// An array with how many particles are held by each TreePiece when sorted.
	std::vector<int> particleCounts;
	
public:
	
	DataManager(const CkArrayID& treePieceID);
	
	~DataManager() { }
	
	/// Called by the Sorter, I ask my TreePieces to evaluate these splitter keys
	void acceptCandidateKeys(const SFC::Key* keys, const int n, int isRefine, const CkCallback& cb);
	/// Called by the Sorter, I save these final keys and the list of which TreePiece is responsible for which interval
	void acceptFinalKeys(const SFC::Key* keys, const int* responsible, unsigned int* bins, const int n, const CkCallback& cb);
	
};

#endif //DATAMANAGER_H
