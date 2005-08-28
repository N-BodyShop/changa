/** \file GenericTreeNode.h
    This file defines the generic tree structures.
    @author Filippo Gioachin (previously Graeme Lufkin and Chao Huang)
    @version 2.0
*/

#ifndef GENERICTREENODE_H
#define GENERICTREENODE_H

#include "pup.h"

#include "OrientedBox.h"
#include "MultipoleMoments.h"

#include "GravityParticle.h"

namespace Tree {

  /** @brief This key is the identification of a node inside the global tree,
      and it is unique for the node. This is used to lookup nodes in any hash or
      cache table.

      The position of the leftmost bit which is 1 represents the depth of the
      node into the tree, all the bits at its right describe the path of this
      node into the tree, and the bits at its left are clearly 0 and unused.
   */
  typedef u_int64_t NodeKey;

  /// This enumeration determines the different types of node a GenericTreeNode can be
  enum NodeType {
    Invalid = 1,
    Bucket,
    Internal,
    Boundary,
    NonLocal,
    Empty,
    Top,
    NonLocalBucket,
    Cached,
    CachedBucket,
    CachedEmpty
  };

  /// A simple counter for how many nodes are empty - statistics
  extern int numEmptyNodes;

  class GenericTreeNode {
  protected:
    NodeType myType;
    NodeKey key;
	
    GenericTreeNode() : myType(Invalid), key(0), parent(0), firstParticle(0), lastParticle(0), remoteIndex(0) { }

  public:

    /// The moments for the gravity computation
    MultipoleMoments moments;
    /// The parent of this node, or null if none
    GenericTreeNode* parent;
    /// The axis-aligned bounding box of this node
    OrientedBox<double> boundingBox;
    /// An index for the first particle contained by this node, 0 means outside the node
    int firstParticle;
    /// An index to the last particle contained by this node, myNumParticles+1 means outside the node
    int lastParticle;
    /// An index to the real location of this node if this node is NonLocal, if
    /// it is Boundary or Internal or Bucket it is equal to thisIndex
    unsigned int remoteIndex;
    /// If this node is partially local, total number of particles contained (across all chares)
    unsigned int particleCount;
	
    //GenericTreeNode() : myType(Invalid), key(0), parent(0), beginParticle(0), endParticle(0), remoteIndex(0) { }

    GenericTreeNode(NodeKey k, NodeType type, int first, int last, GenericTreeNode *p) : myType(type), key(k), parent(p), firstParticle(first), lastParticle(last), remoteIndex(0) { }

    virtual ~GenericTreeNode() { }
    
    inline const NodeType getType() const { return myType; }
    inline void setType(NodeType t) { myType = t; }

    inline const NodeKey getKey() const { return key; }

    /// return the number of children this node has
    virtual unsigned int numChildren() const = 0;
    /// return the pointers to the specified child of this node
    virtual GenericTreeNode* getChildren(int) = 0;
    /// set the specified child of this node to the passed pointer
    virtual void setChildren(int, GenericTreeNode*) = 0;
    /// return the keys for the specified child
    virtual NodeKey getChildKey(int) = 0;
    /// return the key for the parent
    virtual NodeKey getParentKey() = 0;
    /// return an integer with the number of the child reflecting the key
    virtual int whichChild(NodeKey childkey) = 0;

    /// construct the children of the "this" node following the given logical
    /// criteria (Oct/Orb)
    virtual void makeOctChildren(GravityParticle *part, int totalPart, int level) = 0;
    virtual void makeOrbChildren(GravityParticle *part, int totalPart, int level) = 0;

    /// transform an internal node into a bucket
    inline void makeBucket(GravityParticle *part) {
      myType = Bucket;
      for (int i = firstParticle; i <= lastParticle; ++i)
	moments += part[i];
      calculateRadiusFarthestParticle(moments, &part[firstParticle], &part[lastParticle+1]);
    }

    inline void makeEmpty() { myType = Empty; }

    virtual GenericTreeNode *createNew() const = 0;
    virtual GenericTreeNode *clone() const = 0;

    virtual void pup(PUP::er &p, int depth) = 0;
    virtual void pup(PUP::er &p) {
      int iType;
      if(p.isUnpacking()) {
	p | iType;
	myType = (NodeType) iType;
      } else {
	iType = (int) myType;
	p | iType;
      }
      p | key;
      p | moments;
      p | boundingBox;
      p | firstParticle;
      p | lastParticle;
      p | remoteIndex;
      p | particleCount;
    }
  };
  
  /** A table of the nodes in my tree, indexed by their keys.
      @todo XXX: Make this lookup a hash table, so we get O(1) behavior instead of O(log N).
  */
  typedef std::map<NodeKey, GenericTreeNode *> NodeLookupType;

  class BinaryTreeNode : public GenericTreeNode {
  protected:
  public:
    BinaryTreeNode* children[2];

    BinaryTreeNode() : GenericTreeNode() {
      children[0] = 0;
      children[1] = 0;
    }

  public:
    BinaryTreeNode(NodeKey k, NodeType type, int first, int nextlast, BinaryTreeNode *p) : GenericTreeNode(k, type, first, nextlast, p) {
      children[0] = 0;
      children[1] = 0;
    }

    ~BinaryTreeNode() {
      //delete children[0];
      //delete children[1];
    }
    
    unsigned int numChildren() const {
      return 2;
    }
    
    GenericTreeNode* getChildren(int i) {
      CkAssert(i>=0 && i<2);
      return children[i];
    }
    
    void setChildren(int i, GenericTreeNode* node) {
      CkAssert(i>=0 && i<2);
      children[i] = (BinaryTreeNode*)node;
    }
    
    NodeKey getChildKey(int i) {
      CkAssert(i>=0 && i<2);
      return (key<<1) + i;
    }
    
    NodeKey getParentKey() {
      return (key>>1);
    }

    int whichChild(NodeKey child) {
      return (child ^ (key<<1));
    }
    
    bool isLeftChild() const {
      return (dynamic_cast<BinaryTreeNode *>(parent) && dynamic_cast<BinaryTreeNode *>(parent)->children[0] == this);
    }
    
    bool isRightChild() const {
      return (dynamic_cast<BinaryTreeNode *>(parent) && dynamic_cast<BinaryTreeNode *>(parent)->children[1] == this);
    }
    
    BinaryTreeNode* getSibling() const {
      BinaryTreeNode* p = dynamic_cast<BinaryTreeNode *>(parent);
      if(p)
	return (p->children[0] == this ? p->children[1] : p->children[0]);
      else
	return 0;
    }

    void makeOctChildren(GravityParticle *part, int totalPart, int level) {
      children[0] = new BinaryTreeNode();
      children[1] = new BinaryTreeNode();
      children[0]->parent = this;
      children[1]->parent = this;
      children[0]->boundingBox = boundingBox;
      children[1]->boundingBox = boundingBox;
      double split;
      switch (level%3) {
      case 0:
	split = 0.5 * (boundingBox.greater_corner.x + boundingBox.lesser_corner.x);
	children[0]->boundingBox.greater_corner.x = split;
	children[1]->boundingBox.lesser_corner.x = split;
	break;
      case 1:
	split = 0.5 * (boundingBox.greater_corner.y + boundingBox.lesser_corner.y);
	children[0]->boundingBox.greater_corner.y = split;
	children[1]->boundingBox.lesser_corner.y = split;
	break;
      case 2:
	split = 0.5 * (boundingBox.greater_corner.z + boundingBox.lesser_corner.z);
	children[0]->boundingBox.greater_corner.z = split;
	children[1]->boundingBox.lesser_corner.z = split;
	break;
      }
      children[0]->key = key << 1;
      children[1]->key = (key << 1) + 1;
      children[0]->firstParticle = firstParticle;
      children[1]->lastParticle = lastParticle;

      SFC::Key mask = SFC::Key(1) << (62 - level);
      SFC::Key leftBit = part[firstParticle].key & mask;
      SFC::Key rightBit = part[lastParticle].key & mask;

      if (leftBit == rightBit) {
	// we have only one child, the other is either NonLocal or Empty
	if (leftBit) {
	  // left child missing
	  if (firstParticle == 0) {
	    children[0]->myType = NonLocal;
	    children[1]->myType = Boundary;
	  } else {
	    children[0]->myType = Empty;
	    children[1]->myType = lastParticle==totalPart+1 ? Boundary : Internal;
	  }
	  children[0]->lastParticle = firstParticle-1;
	  children[1]->firstParticle = firstParticle;
	  children[1]->particleCount = particleCount;
	} else {
	  // right child missing
	  if (lastParticle == totalPart+1) {
	    children[1]->myType = NonLocal;
	    children[0]->myType = Boundary;
	  } else {
	    children[1]->myType = Empty;
	    children[0]->myType = firstParticle==0 ? Boundary : Internal;
	  }
	  children[1]->firstParticle = lastParticle+1;
	  children[0]->lastParticle = lastParticle;
	  children[0]->particleCount = particleCount;
	}
      } else if (leftBit < rightBit) {
	// both children are present
	if (firstParticle == 0 && leftBit != (part[1].key & mask)) {
	  // the left child is NonLocal
	  children[0]->myType = NonLocal;
	  children[0]->lastParticle = firstParticle;
	  children[1]->myType = lastParticle==totalPart+1 ? Boundary : Internal;
	  children[1]->firstParticle = firstParticle+1;
	  children[1]->particleCount = particleCount;
	} else if (lastParticle == totalPart+1 && rightBit != (part[totalPart].key & mask)) {
	  // the right child is NonLocal
	  children[1]->myType = NonLocal;
	  children[1]->firstParticle = lastParticle;
	  children[0]->myType = firstParticle==0 ? Boundary : Internal;
	  children[0]->lastParticle = lastParticle-1;
	  children[0]->particleCount = particleCount;
	} else {
	  // the splitting point is somewhere in the middle
	  GravityParticle *splitParticle = std::lower_bound(&part[firstParticle], &part[lastParticle+1], part[lastParticle].key & ((~ (SFC::Key)0) << (62-level)));
	  children[0]->lastParticle = splitParticle - part - 1;
	  children[1]->firstParticle = splitParticle - part;
	  children[0]->particleCount = children[0]->lastParticle - firstParticle + 1;
	  children[1]->particleCount = lastParticle - children[1]->firstParticle + 1;
	  children[0]->myType = children[0]->particleCount==0 ? Empty : (firstParticle==0 ? Boundary : Internal);
	  children[1]->myType = children[1]->particleCount==0 ? Empty : (lastParticle==totalPart+1 ? Boundary : Internal);
	  if (firstParticle == 0) children[0]->particleCount --;
	  if (lastParticle == totalPart+1) children[1]->particleCount --;
	}
      } else {
	CkAbort("Ok, the particles should be ordered so, how the hell did we get here?!?");
      }
    }

    void makeOrbChildren(GravityParticle *part, int totalPart, int level) {}

    GenericTreeNode *createNew() const {
      return new BinaryTreeNode();
    }

    GenericTreeNode *clone() const {
      BinaryTreeNode *tmp = new BinaryTreeNode();
      *tmp = *this;
      return tmp;
    }

    void pup(PUP::er &p) { pup(p, -1); }
    void pup(PUP::er &p, int depth) {
      //CkPrintf("Pupper of BinaryTreeNode(%d) called for %s (%d)\n",depth,p.isPacking()?"Packing":p.isUnpacking()?"Unpacking":"Sizing",p.isSizing()?((PUP::sizer*)&p)->size():((PUP::mem*)&p)->size());
      GenericTreeNode::pup(p);
      int isNull;
      for (int i=0; i<2; ++i) {
	isNull = (children[i]==NULL || depth==0) ? 0 : 1;
	p | isNull;
	CkAssert(isNull==0 || isNull==1);
	if (isNull != 0 && depth != 0) {
	  if (p.isUnpacking()) children[i] = new BinaryTreeNode();
	  children[i]->pup(p, depth-1);
	  if (p.isUnpacking()) children[i]->parent = this;
	}
      }
    }
  };
  
  class OctTreeNode : public GenericTreeNode {
  protected:
  public:
    /// the 8 different children
    OctTreeNode* children[8];

    OctTreeNode() : GenericTreeNode() {
      children[0] = 0;
      children[1] = 0;
      children[2] = 0;
      children[3] = 0;
      children[4] = 0;
      children[5] = 0;
      children[6] = 0;
      children[7] = 0;
    }

  public:
    OctTreeNode(NodeKey k, NodeType type, int first, int last, BinaryTreeNode *p) : GenericTreeNode(k, type, first, last, p) {
      children[0] = 0;
      children[1] = 0;
      children[2] = 0;
      children[3] = 0;
      children[4] = 0;
      children[5] = 0;
      children[6] = 0;
      children[7] = 0;
    }
    
    virtual ~OctTreeNode() {
      //delete children[0];
      //delete children[1];
      //delete children[2];
      //delete children[3];
      //delete children[4];
      //delete children[5];
      //delete children[6];
      //delete children[7];
    }

    virtual unsigned int numChildren() const {
      return 8;
    }
    
    virtual GenericTreeNode* getChildren(int i) {
      CkAssert(i>=0 && i<8);
      return children[i];
    }

    virtual void setChildren(int i, GenericTreeNode* node) {
      CkAssert(i>=0 && i<8);
      children[i] = (OctTreeNode*)node;
    }

    NodeKey getChildKey(int i) {
      CkAssert(i>=0 && i<8);
      return (key<<3) + i;
    }

    NodeKey getParentKey() {
      return (key>>3);
    }

    int whichChild(NodeKey child) {
      return (child ^ (key<<3));
    }
    
    void makeOctChildren(GravityParticle *part, int totalPart, int level) {}

    void makeOrbChildren(GravityParticle *part, int totalPart, int level) {}

    GenericTreeNode *createNew() const {
      return new OctTreeNode();
    }

    GenericTreeNode *clone() const {
      OctTreeNode *tmp = new OctTreeNode();
      *tmp = *this;
      return tmp;
    }

    void pup(PUP::er &p);
    void pup(PUP::er &p, int depth);
  };

  /** Possible types of trees: the first word means the physical structure used
      to implement the tree (Binary or Oct), the second work means the logic
      used to decompose the particles of the space amoung the nodes (Oct or
      ORB). Notice that the SFC domain decomposition builds an oct-tree!
   */
  enum GenericTrees {
    Binary_Oct,
    Oct_Oct,
    Binary_ORB
  };

} //close namespace Tree

inline void operator|(PUP::er &p,Tree::NodeType &nt) {
  int nti;
  if (p.isUnpacking()) {
    p | nti;
    nt = (Tree::NodeType)nti;
  } else {
    nti = (int)nt;
    p | nti;
  }
}

inline void operator|(PUP::er &p,Tree::GenericTrees &gt) {
  int gti;
  if (p.isUnpacking()) {
    p | gti;
    gt = (Tree::GenericTrees)gti;
  } else {
    gti = (int)gt;
    p | gti;
  }
}

#endif //GENERICTREENODE_H
