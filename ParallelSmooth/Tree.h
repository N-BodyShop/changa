/** @file Tree.h
 This header defines the structures used to make up trees.
 @author Graeme Lufkin (gwl@u.washington.edu)
 */

#ifndef TREE_H
#define TREE_H

#include <sstream>

#include "pup.h"

#include "OrientedBox.h"
#include "Particle.h"

enum NodeType {
	Invalid,
	Bucket,
	Internal,
	NonLocal,
	Boundary,
	Empty,
	Top
};

template <typename T>
OrientedBox<T> cutBoxDown(const OrientedBox<T>& box, const int axis) {
	OrientedBox<T> newBox = box;
	switch(axis) {
		case 0: //cut x
			newBox.greater_corner.x = (box.greater_corner.x + box.lesser_corner.x) / 2.0;
			break;
		case 1: //cut y
			newBox.greater_corner.y = (box.greater_corner.y + box.lesser_corner.y) / 2.0;				
			break;
		case 2: //cut z
			newBox.greater_corner.z = (box.greater_corner.z + box.lesser_corner.z) / 2.0;
			break;
	}
	return newBox;
}

template <typename T>
OrientedBox<T> cutBoxUp(const OrientedBox<T>& box, const int axis) {
	OrientedBox<T> newBox = box;
	switch(axis) {
		case 0: //cut x
			newBox.lesser_corner.x = (box.greater_corner.x + box.lesser_corner.x) / 2.0;
			break;
		case 1: //cut y
			newBox.lesser_corner.y = (box.greater_corner.y + box.lesser_corner.y) / 2.0;				
			break;
		case 2: //cut z
			newBox.lesser_corner.z = (box.greater_corner.z + box.lesser_corner.z) / 2.0;
			break;
	}
	return newBox;
}

/**
 This class represents a node in the tree of particles.
 The responsibility for building a valid tree using these nodes falls
 to another part of the program.
 If the node is actually a bucket containing particles, then its
 key will have the bucketFlag set.  The beginBucket and endBucket pointers
 will then be valid, pointing to the first particle in the bucket and
 one past the last particle in the bucket, respectively.
 If it is not a bucket, then the leftChild and rightChild pointers are valid.
 You need to check chareID to see if one of the
 children is non-local.  If chareID is not zero, then check its sign.
 If positive, then the rightChild is non-local.  If negative, then the 
 leftChild is non-local.  The value of chareID (making positive, and subtracting
 one) is the index of the chare (not neccessarily the only!) that actually 
 owns (or knows more about) the child node.
 */
class TreeNode : public PUP::able {
private:
	NodeType myType;

public:
	
	/** A special key that is set when a node is really a  bucket containing particles. */
	const static Key bucketFlag = (static_cast<Key>(1) << 63);
	
	/// The maximum number of particles in a bucket.
	const static int maxBucketSize = 12;
	
	/// A pointer to the parent of this node.
	TreeNode* parent;
	
	/// This node's key determines what region of space it contains.
	Key key;
	
	/// This tells you how many bits of the key to look at.
	unsigned char level;
	
	/// This points to a left child node or the first particle in a bucket, depending on \c bucketFlag being set
	union {
		TreeNode* leftChild;
		FullParticle* beginBucket; //only access if getType() == Bucket
		unsigned int chareID; //only access if getType() == NonLocal
	};
	
	/// This points to a right child node or one past the last particle in a bucket, depending on \c bucketFlag being set
	union {
		TreeNode* rightChild;
		FullParticle* endBucket;
	};
	
	/// An axis-oriented box giving the region this node represents.
	OrientedBox<double> box;
	
	/** Given a key, create a vector of floats representing a position.
	 This is almost the inverse of the makeKey() function.  Since the key
	 does not use all the bits, each float generated here will have its last
	 two mantissa bits set zero, regardless of the values when the key was
	 originally generated.
	 */
	static inline Vector3D<float> makeVector(Key k) {
		Vector3D<float> v(1.0, 1.0, 1.0);
		int* ix = reinterpret_cast<int *>(&v.x);
		int* iy = reinterpret_cast<int *>(&v.y);
		int* iz = reinterpret_cast<int *>(&v.z);
		for(int mask = (1 << 2); mask < (1 << 23); mask <<= 1) {
			if(k & 4)
				*ix |= mask;
			if(k & 2)
				*iy |= mask;
			if(k & 1)
				*iz |= mask;
			k >>= 3;	
		}
		return v;
	}
		
	TreeNode() : myType(Invalid), parent(0), key(0), level(0), leftChild(0), rightChild(0) { }
	
	TreeNode* createLeftChild() {
		TreeNode* child = new TreeNode;
		child->parent = this;
		child->key = key;
		child->level = level + 1;
		child->box = cutBoxDown(box, level % 3);
		leftChild = child;
		return child;
	}
	
	TreeNode* createRightChild() {
		TreeNode* child = new TreeNode;
		child->parent = this;
		child->key = key | (static_cast<Key>(1) << (62 - level));
		child->level = level + 1;
		child->box = cutBoxUp(box, level % 3);
		rightChild = child;
		return child;
	}
	
	~TreeNode() {
		if(myType == NonLocal)
			delete rightChild;
		else if(myType != Bucket) {
			delete leftChild;
			delete rightChild;
		}
	}
	
	NodeType getType() const { return myType; }
	void setType(NodeType t) { myType = t; }
	
	/// Is this node really a bucket of particles?
	inline bool isBucket() const {
		return myType == Bucket;
	}
	
	/// Is this node non-local?
	inline bool isNonLocal() const {
		return myType == NonLocal;
	}
	
	/// Is this node completely owned by this piece
	inline bool isInternal() const {
		return myType == Internal || myType == Bucket;
	}
	
	/// Is this node someone's left child?
	inline bool isLeftChild() const {
		return (parent != 0 && parent->leftChild == this);
		//return (parent != 0 && parent->chareID >= 0 && parent->leftChild == this);
	}

	/// Is this node someone's right child?
	inline bool isRightChild() const {
		return (parent != 0 && parent->rightChild == this);
		//return (parent != 0 && parent->chareID <= 0 && parent->rightChild == this);
	}
	
	/** The lookup key of this node.
	 A lookup key for a node shifts the key right, leaving only a leading 1 
	 and the \c level bits that uniquely identify this node.
	 */
	inline Key lookupKey() const {
		return (key | bucketFlag) >> (63 - level);
	}
	
	/// The lookup key of this node's left child.
	inline Key leftChildLookupKey() const {
		return (key | bucketFlag) >> (63 - level - 1);
	}
	
	/// The lookup key of this node's right child.
	inline Key rightChildLookupKey() const {
		return ((key | bucketFlag) >> (63 - level - 1)) | static_cast<Key>(1);
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
	
	PUPable_decl(TreeNode);
	
	TreeNode(CkMigrateMessage* m) : PUP::able(m), myType(Invalid), parent(0), key(0), level(0), leftChild(0), rightChild(0) { }
	/*
	virtual void pup(PUP::er& p) {
		PUP::able::pup(p);
		p(myType);
		p(key);
		p(level);
		
		if(myType == NonLocal)
			p(chareID);
		else if(myType == Bucket) {
			if(p.isUnpacking()) {
				p(beginBucket);
				p(endBucket);
				beginBucket += ???;
				endBucket += ???;
			} else {
				p(beginBucket - ???);
				p(endBucket - ???);
			}
		} else {
			p | leftChild;				
			p | rightChild;
			if(p.isUnpacking()) {
				if(leftChild)
					leftChild->parent = this;
				if(rightChild)
					rightChild->parent = this;
			}
		}
		
		p(box);
	}
	*/
};

#endif //TREE_H
