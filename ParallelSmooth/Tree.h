/** \file Tree.h
 This header defines the structures used to make up trees.
 \author Graeme Lufkin (gwl@u.washington.edu)
 */

#ifndef TREE_H
#define TREE_H

#include <iostream>
#include <bitset>

#include "OrientedBox.h"
#include "Particle.h"

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
class TreeNode {
private:
		

public:
	
	/** A special key that is set when a node is really a  bucket containing particles. */
	const static Key bucketFlag = ((Key) 1 << 63);
	
	/// The maximum number of particles in a bucket.
	const static int maxBucketSize = 12;
	
	/// A pointer to the parent of this node.
	TreeNode* parent;
	
	/** This contains the chare index of a non-local child.
	 If this is negative, change the sign and subtract one to get the index
	 of the non-local left child.  If it is positive, subtract one to get
	 the index of the non-local right child.  If it is zero, then both children
	 are local. */
	int chareID;
	
	/// This node's key determines what region of space it contains.
	Key key;
	
	/// This tells you how many bits of the key to look at.
	int level;
	
	/// This points to a left child node or the first particle in a bucket, depending on \c bucketFlag being set
	union {
		TreeNode* leftChild;
		FullParticle* beginBucket;
	};
	
	/// This points to a right child node or one past the last particle in a bucket, depending on \c bucketFlag being set
	union {
		TreeNode* rightChild;
		FullParticle* endBucket;
	};
	
	/// An axis-oriented box giving the region this node represents.
	OrientedBox<float> box;
	/// A flag that is set if this node and all its children are owned by the same chare.
	bool whollyOwned;
	
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
		
	TreeNode() : parent(0), chareID(0), key(0), leftChild(0), rightChild(0) { }
	
	TreeNode(TreeNode* p, const Key k) : parent(p), chareID(0), key(k), leftChild(0), rightChild(0) {
		// XXX: Optimization available to make bounding box from parent's box
		if(parent == 0)
			level = 0;
		else
			level = parent->level + 1;
		
		Vector3D<float> oneAndAHalf(1.5, 1.5, 1.5);
		box = OrientedBox<float>(makeVector(key) - oneAndAHalf, makeVector(key | ((static_cast<Key>(1) << (63 - level)) - 1)) - oneAndAHalf);
		
		whollyOwned = false;
	}
	
	~TreeNode() { }
	
	/// Is this node really a bucket of particles?
	inline bool isBucket() const {
		return key & bucketFlag;
	}
	
	/// Is this node someone's left child?
	inline bool isLeftChild() const {
		return (parent != 0 && parent->chareID >= 0 && parent->leftChild == this);
	}

	/// Is this node someone's right child?
	inline bool isRightChild() const {
		return (parent != 0 && parent->chareID <= 0 && parent->rightChild == this);
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
	
};

#endif //TREE_H
