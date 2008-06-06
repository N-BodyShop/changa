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
#include "GenericTreeNode.h"

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
	// A pointer to the CacheManager in the local processor
	//CacheManager *localCache;

protected:
	
	/// The array of splitter keys for the sort.
	std::vector<SFC::Key> boundaryKeys;
	
	/// An array identifying which chare is responsible for each interval of keys.
	std::vector<int> responsibleIndex;
	
	/// An array with how many particles are held by each TreePiece when sorted.
	std::vector<int> particleCounts;
	
	/// An array with a list of the first and last particle of each TreePiece
	SFC::Key * splitters;
	/// The size of the array splitters
	int numSplitters;
	
	/// A list of roots of the TreePieces in this node
	CkVec<Tree::GenericTreeNode*> registeredChares;
	/// The root of the combined trees
	Tree::GenericTreeNode * root;
	/// The lookup table for nodes created by the DataManager to combine trees
	Tree::NodeLookupType nodeLookupTable;
	
	/// Number of chunks in which the tree was splitted during last combine operation
	int oldNumChunks;
	/// Nodes currently used as roots for remote computation
	Tree::NodeKey *chunkRoots;
    /// Lookup table for the chunkRoots
    Tree::NodeLookupType chunkRootTable;

public:
	
	DataManager(const CkArrayID& treePieceID);
	DataManager(CkMigrateMessage *);
private:
    void init();
public:
	
	~DataManager() { }
	
	/// Collect the boundaries of all TreePieces, and trigger the real treebuild
	void collectSplitters(CkReductionMsg* m);
	/// Called by the Sorter, I ask my TreePieces to evaluate these splitter keys
	//void acceptCandidateKeys(const SFC::Key* keys, const int n, int isRefine, const CkCallback& cb);
	/// Called by the Sorter, I save these final keys and the list of which TreePiece is responsible for which interval
	void acceptFinalKeys(const SFC::Key* keys, const int* responsible, unsigned int* bins, const int n, const CkCallback& cb);
	void pup(PUP::er& p);
	
	// Functions used to create a tree inside the DataManager comprising
	// all the trees in the TreePieces in the local node
private:
	Tree::GenericTreeNode *buildProcessorTree(int n, Tree::GenericTreeNode **gtn);
	int createLookupRoots(Tree::GenericTreeNode *node, Tree::NodeKey *keys);
public:
	void notifyPresence(Tree::GenericTreeNode *root);
	void combineLocalTrees(CkReductionMsg *msg);
        void getChunks(int &num, Tree::NodeKey *&roots);
};

#endif //DATAMANAGER_H
