/** @file TreeNode.h
 This header defines the structures used to make up trees.
 @author Graeme Lufkin (gwl@u.washington.edu)
 */

#ifndef TREENODE_H
#define TREENODE_H

#include <sstream>

#include "pup.h"

#include "OrientedBox.h"
#include "SFC.h"
#include "GravityTreeNode.h"

namespace Tree {

using namespace SFC;

class SFCTreeNode : public GravityTreeNode {
protected:
	/** The highest bit set in a 64-bit integer, used with shifts 
	 to identify the level of a key. */
	const static Key lastBitFlag = (static_cast<Key>(1) << 63);

	SFCTreeNode* populateLeftChild(SFCTreeNode* child) {
		child->parent = this;
		child->key = key;
		child->level = level + 1;
		child->boundingBox = cutBoxLeft(boundingBox, level % 3);
		leftChild = child;
		return child;
	}
	
	SFCTreeNode* populateRightChild(SFCTreeNode* child) {
		child->parent = this;
		child->key = key | (static_cast<Key>(1) << (62 - level));
		child->level = level + 1;
		child->boundingBox = cutBoxRight(boundingBox, level % 3);
		rightChild = child;
		return child;
	}
	
public:
		
	Key key;
	unsigned char level;
	unsigned int numOwners;
	
	SFCTreeNode() : key(0), level(0), numOwners(0) { }
	
	SFCTreeNode* createLeftChild() {
		return populateLeftChild(new SFCTreeNode);
	}
	
	SFCTreeNode* createRightChild() {
		return populateRightChild(new SFCTreeNode);
	}
	
	/** The lookup key of this node.
	 A lookup key for a node shifts the key right, leaving only a leading 1 
	 and the \c level bits that uniquely identify this node.
	 */
	inline Key lookupKey() const {
		return (key | lastBitFlag) >> (63 - level);
	}
	
	/// The lookup key of this node's left child.
	inline Key leftChildLookupKey() const {
		return (key | lastBitFlag) >> (63 - level - 1);
	}
	
	/// The lookup key of this node's right child.
	inline Key rightChildLookupKey() const {
		return ((key | lastBitFlag) >> (63 - level - 1)) | static_cast<Key>(1);
	}
	
	/// The key of this node's right child (the left child's key is the same as the parent's)
	inline Key rightChildKey() const {
		return key | (static_cast<Key>(1) << (62 - level));
	}
	
	/// The first key in this node (on the left)
	inline Key leftBoundary() const {
		return key;
	}
	
	/// The last key in this node (on the right)
	inline Key rightBoundary() const {
		return key | ((static_cast<Key>(1) << (63 - level)) - 1);
	}
};

} //close namespace Tree

namespace TreeStuff {

using namespace SFC;

inline std::string keyBits(const Key k, const int numBits) {
	std::ostringstream oss;
	oss << "N";
	for(int i = 0; i < numBits; i++)
		oss << (k & (static_cast<Key>(1) << (62 - i)) ? 1 : 0);
	return oss.str();
}

const double opening_geometry_factor = 2 / sqrt(3.0);

/// The maximum number of particles in a bucket.
extern int maxBucketSize;
	

} //close namespace TreeStuff

#endif //TREENODE_H
