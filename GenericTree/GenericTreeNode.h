/** \file GenericTreeNode.h
 This file defines the a generic tree node class.
 @author Graeme Lufkin (gwl@u.washington.edu), originally Charm folks
 @version 1.0
 */

#ifndef GENERICTREENODE_H
#define GENERICTREENODE_H

#include "OrientedBox.h"

namespace Tree {

/// This enumeration determines the different types of node a GenericTreeNode can be
enum NodeType {
	Invalid = 1,
	Bucket,
	Internal,
	Boundary,
	NonLocal,
	Empty,
	Top
};

class GenericTreeNode {
protected:
	NodeType myType;
	
public:
	
	/// The parent of this node, or null if none
	GenericTreeNode* parent;
	/// The axis-aligned bounding box of this node
	OrientedBox<double> boundingBox;
	/// An index for the first particle contained by this node
	u_int64_t beginParticle;
	/// An index to one past the last particle contained by this node
	u_int64_t endParticle;
	/// An index to the real location of this node
	unsigned int remoteIndex;
	
	GenericTreeNode() : myType(Invalid), parent(0), beginParticle(0), endParticle(0), remoteIndex(0) { }

	virtual ~GenericTreeNode() { }
	
	NodeType getType() const {
		return myType;
	}
	
	void setType(NodeType t) {
		myType = t;
	}
	
	virtual unsigned int numChildren() const = 0;
	
	virtual GenericTreeNode** getChildren() = 0;
	
};


class BinaryTreeNode : public GenericTreeNode {
public:
	BinaryTreeNode* leftChild;
	BinaryTreeNode* rightChild;
		
	BinaryTreeNode() : leftChild(0), rightChild(0) { }

	virtual ~BinaryTreeNode() {
		delete leftChild;
		delete rightChild;
	}
	
	virtual unsigned int numChildren() const {
		return 2;
	}
	
	virtual GenericTreeNode** getChildren() {
		return reinterpret_cast<GenericTreeNode **>(&leftChild);
	}
	
	bool isLeftChild() const {
		return (dynamic_cast<BinaryTreeNode *>(parent) && dynamic_cast<BinaryTreeNode *>(parent)->leftChild == this);
	}

	bool isRightChild() const {
		return (dynamic_cast<BinaryTreeNode *>(parent) && dynamic_cast<BinaryTreeNode *>(parent)->rightChild == this);
	}
	
	BinaryTreeNode* getSibling() const {
		BinaryTreeNode* p = dynamic_cast<BinaryTreeNode *>(parent);
		if(p)
			return (p->leftChild == this ? p->rightChild : p->leftChild);
		else
			return 0;
	}
};

class OctTreeNode : public virtual GenericTreeNode {
protected:
	GenericTreeNode* children[8];
public:
		
	OctTreeNode() {
		children[0] = 0;
		children[1] = 0;
		children[2] = 0;
		children[3] = 0;
		children[4] = 0;
		children[5] = 0;
		children[6] = 0;
		children[7] = 0;
	}

	virtual unsigned int numChildren() const {
		return 8;
	}
	
	virtual GenericTreeNode** getChildren() {
		return children;
	}
};

} //close namespace Tree

#endif //GENERICTREENODE_H
