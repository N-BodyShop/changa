/** @file TreePiece.h
 @author Graeme Lufkin (gwl@u.washington.edu)
 */

#ifndef TREEPIECE_H
#define TREEPIECE_H

#include <vector>
#include <map>

#include "NChilada.h"

class TreeNode;
class DataManager;

/** A piece of tree represents the basic computation unit in our framework.
 We have an array of tree pieces that communicate when performing calculations.
 */
class TreePiece : public ArrayElement1D {
	/// The number of boundaries I've received (0, 1, or 2).
	int numBoundaries;

protected:
	
	/// A proxy to the array I'm a member of.
	CProxy_TreePiece treePieces;
	/// A proxy to the DataManager.
	CProxy_DataManager dataManager;
	/// A local pointer to my DataManager.
	DataManager* dm;
	
	/// The bounding box of all the particles.
	OrientedBox<double> boundingBox;
	/// My index in the responsibility array.
	int myPlace;
	/// A callback to keep around for a while.
	CkCallback callback;
	/// The number of TreePieces in the array I'm in.
	int numTreePieces;

	/** This contains the particles I am given initially, and are un-ordered. */
	std::vector<FullParticle> myParticles;
	/** This contains my set of particles after they have been sorted. */
	std::vector<FullParticle> mySortedParticles;
	
	/// The number of particles I am responsible for.
	int numParticles;
	/// The root of my part of the tree.
	TreeNode* root;
	/// The keys determined by the Sorter that separate me from my neighbors.
	Key leftSplitter, rightSplitter;
	
	FullParticle* leftBoundary;
	FullParticle* rightBoundary;
	
	/** A table of the nodes in my tree, indexed by their keys.
	 @todo XXX: Make this lookup a hash table, so we get O(1) behavior instead of O(log N).
	 */
	std::map<Key, TreeNode *> nodeLookup;
	
	TreeNode* lookupLeftChild(TreeNode* node);
	TreeNode* lookupRightChild(TreeNode* node);
	void buildTree(TreeNode* parent, FullParticle* leftParticle, FullParticle* rightParticle);

public:

	TreePiece();

	/** The migration constructor is empty. */
	TreePiece(CkMigrateMessage* m) { }
	
	~TreePiece();
	
	/** Inform the DataManager of my node that I'm here.
	 The callback will receive a CkReductionMsg containing no data.
	 */
	void registerWithDataManager(const CkGroupID& dataManagerID, const CkCallback& cb);
	
	/** Build my part of the tree.
	 The callback will receive a CkReductionMsg containing no data.
	 */
	void startTreeBuild(const CkCallback& cb);
	
	/** Write a graphviz graph of my piece of the tree and a tipsy box macro to a file. 
	 */
	void report(const CkCallback& cb);
	
	// These functions are conceptually private, but are public because charm entry methods are always public
	
	void receiveParticles(const FullParticle* particles, const int n, const CkCallback& cb);
	void assignKeys(CkReductionMsg* m);
	
	void evaluateBoundaries(const CkCallback& cb);
	void unshuffleParticles(const CkCallback& cb);
	void acceptSortedParticles(const FullParticle* particles, const int n);
	void shareBoundaries(CkReductionMsg* m);
	void acceptBoundaryKey(const Key k);
	
	//void generateImage(liveVizRequestMsg* m);
	//void applyPerParticleFunction(const CkCallback& cb);
	
	//void pup(PUP::er& p);
};


#endif //TREEPIECE_H
