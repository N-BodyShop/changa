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
/*
class Moments {
public:
	
	double totalMass;
	/// The center of mass (not divided by the total mass)
	Vector3D<double> cm;
	
	/// Provides a physical size for this multipole expansion
	double radius;
	
	Moments(double m = 0, const Vector3D<double>& com = Vector3D<double>()) : totalMass(m), cm(com) { }
	
	Moments& operator+=(const Moments& m) {
		totalMass += m.totalMass;
		cm += m.cm;
		return *this;
	}
	
	template <typename ParticleType>
	Moments& operator+=(const ParticleType& p) {
		totalMass += p.mass;
		cm += p.mass * p.position;
		return *this;
	}
	
	Moments operator-(const Moments& m) {
		return Moments(totalMass - m.totalMass, cm - m.cm);
	}
	
	void clear() {
		totalMass = 0;
		cm.x = cm.y = cm.z = 0;
		radius = 0;
	}
	
	/// Determine the distance from the cm to the farthest corner of a box
	void calculateRadius(const OrientedBox<double>& box) {
		Vector3D<double> delta1 = cm / totalMass - box.lesser_corner;	
		Vector3D<double> delta2 = box.greater_corner - cm / totalMass;
		delta1.x = max(delta1.x, delta2.x);	
		delta1.y = max(delta1.y, delta2.y);	
		delta1.z = max(delta1.z, delta2.z);
		radius = delta1.length();
	}
};

inline void operator|(PUP::er& p, Moments& m) {
	p | m.totalMass;
	p | m.cm;
	p | m.radius;
}

enum NodeType {
	Invalid,
	Bucket,
	Internal,
	NonLocal,
	Empty,
	Boundary,
	Top
};
*/
/// The maximum number of particles in a bucket.
extern int maxBucketSize;
	
/**
 This class represents a node in the tree of particles.
 The responsibility for building a valid tree using these nodes falls
 to another part of the program.
 If the node is actually a bucket containing particles, then its
 NodeType will be Bucket.  The beginBucket and endBucket indices
 will then be valid, indexing to the first particle in the bucket and
 one past the last particle in the bucket, respectively.
 If it is not a bucket, then the leftChild and rightChild pointers are valid.
 If its NodeType is Boundary, then one of its children is NonLocal.
 If its NodeType is NonLocal, then chareID contains the index in the chare
 array to the piece which owns (or know more about) the node.
 */
 /*
class TreeNode : public PUP::able {
protected:
	NodeType myType;

	/// The highest bit set in a 64-bit integer, used with shifts to identify the level of a key. 
	const static Key lastBitFlag = (static_cast<Key>(1) << 63);

	inline TreeNode* populateLeftChild(TreeNode* child) {
		child->parent = this;
		child->key = key;
		child->level = level + 1;
		child->boundingBox = cutBoxLeft(boundingBox, level % 3);
		leftChild = child;
		return child;
	}
	
	inline TreeNode* populateRightChild(TreeNode* child) {
		child->parent = this;
		child->key = key | (static_cast<Key>(1) << (62 - level));
		child->level = level + 1;
		child->boundingBox = cutBoxRight(boundingBox, level % 3);
		rightChild = child;
		return child;
	}
	
public:
	
	/// A pointer to the parent of this node.
	TreeNode* parent;
	
	/// This node's key determines what region of space it contains.
	Key key;
	
	/// This tells you how many bits of the key to look at.
	unsigned char level;
	
	/// This points to a left child node or the first particle in a bucket, depending on \c bucketFlag being set
	union {
		TreeNode* leftChild;
		unsigned int beginBucket; //only access if getType() == Bucket
		unsigned int chareID; //only access if getType() == NonLocal
	};
	
	/// This points to a right child node or one past the last particle in a bucket, depending on \c bucketFlag being set
	union {
		TreeNode* rightChild;
		unsigned int endBucket;
	};
	
	/// An axis-oriented box giving the region this node represents.
	OrientedBox<double> boundingBox;

	/// The gravitational moments of this node
	Moments moments;
	
	/// The number of particles contained by this node
	u_int64_t numParticles;
	
	/// The number of TreePieces that hold part of this node
	unsigned int numOwners;
		
	TreeNode() : myType(Invalid), parent(0), key(0), level(0), leftChild(0), rightChild(0), numParticles(0), numOwners(0) { }
	
	inline TreeNode* createLeftChild() {
		return populateLeftChild(new TreeNode);
	}
	
	inline TreeNode* createRightChild() {
		return populateRightChild(new TreeNode);
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
	}

	/// Is this node someone's right child?
	inline bool isRightChild() const {
		return (parent != 0 && parent->rightChild == this);
	}
	
	TreeNode* getSibling() const {
		if(parent)
			return (parent->leftChild == this ? parent->rightChild : parent->leftChild);
		else
			return 0;
	}
	
	/** The lookup key of this node.
	 A lookup key for a node shifts the key right, leaving only a leading 1 
	 and the \c level bits that uniquely identify this node.
	 *//*
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
	
	PUPable_decl(TreeNode);	
	TreeNode(CkMigrateMessage* m) : PUP::able(m), myType(Invalid), parent(0), key(0), level(0), leftChild(0), rightChild(0) { }
	
	virtual void pup(PUP::er& p) {
		PUP::able::pup(p);
		
		if(p.isUnpacking()) {
			int myType_int;
			p | myType_int;
			myType = NodeType(myType_int);
		} else
			p | static_cast<int>(myType);
		
		p | key;
		p | level;
		
		if(myType == NonLocal)
			p | chareID;
		else if(myType == Bucket) {
			p | beginBucket;
			p | endBucket;
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
		
		p | boundingBox;
		p | moments;		
		p | numParticles;
		p | numOwners;
	}
};
*/
} //close namespace TreeStuff

#endif //TREENODE_H
