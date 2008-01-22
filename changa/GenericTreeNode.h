/** \file GenericTreeNode.h
    This file defines the generic tree structures.
    @author Filippo Gioachin (previously Graeme Lufkin and Chao Huang)
    @version 2.0
*/

#ifndef GENERICTREENODE_H
#define GENERICTREENODE_H

#include "pup.h"
#include "ckpool.h"

#include <map>
#include <vector>
#include <algorithm>
#include <set>

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
    u_int64_t usedBy;
	
    GenericTreeNode() : myType(Invalid), key(0), parent(0), firstParticle(0), lastParticle(0), remoteIndex(0), usedBy(0) {
#if COSMO_STATS > 0
      used = false;
#endif
#if INTERLIST_VER > 0
      visitedR=false;
      visitedL=false;
      bucketListIndex=-1;
      startBucket=-1;
#endif
    }

  public:
#if COSMO_STATS > 0
    bool used; // FIXME: this field can now be replaced by "usedBy"
#endif

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
    /// Pointer to the first particle in this node
    GravityParticle *particlePointer;

    /// The greatest rung amoung the particles contained in this node (greater
    /// means faster). This information is limited to the nodes in the current
    /// TreePiece, and do not consider non-local data.
    int rungs;

#if INTERLIST_VER > 0
    bool visitedR;
    bool visitedL;
    int bucketListIndex;
    int startBucket;
#endif
    //GenericTreeNode() : myType(Invalid), key(0), parent(0), beginParticle(0), endParticle(0), remoteIndex(0) { }

    GenericTreeNode(NodeKey k, NodeType type, int first, int last, GenericTreeNode *p) : myType(type), key(k), parent(p), firstParticle(first), lastParticle(last), remoteIndex(0), usedBy(0) { 
    #if INTERLIST_VER > 0
      visitedR=false;
      visitedL=false;
			bucketListIndex=-1;
      startBucket=-1;
    #endif
    }

    virtual ~GenericTreeNode() { }
    virtual void fullyDelete() = 0;
    
    inline NodeType getType() const { return myType; }
    inline void setType(NodeType t) { myType = t; }

    inline NodeKey getKey() const { return key; }

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
#if INTERLIST_VER > 0
    virtual bool contains(NodeKey nodekey) = 0;
#endif
    
    // these two functions are used to track the communication between objects:
    // a nodes is marked usedBy when a local TreePiece has touched it
    void markUsedBy(int index) { usedBy |= (((u_int64_t)1) << index); }
    bool isUsedBy(int index) { return (usedBy & (((u_int64_t)1) << index)); }

    /// construct the children of the "this" node following the given logical
    /// criteria (Oct/Orb)
    virtual void makeOctChildren(GravityParticle *part, int totalPart, int level) = 0;
    virtual void makeOrbChildren(GravityParticle *part, int totalPart, int level, int rootsLevel, bool (*compFnPtr[])(GravityParticle, GravityParticle)) = 0;

    // get the number of chunks possible for the given request
    // @return a number greater or iqual to th the request
    //virtual int getNumChunks(int num) = 0;
    /// get the nodes corresponding to a particular number of chunks requested
    virtual void getChunks(int num, NodeKey *&ret) = 0;

    /// transform an internal node into a bucket
    inline void makeBucket(GravityParticle *part) {
      myType = Bucket;
      boundingBox.reset();
      rungs = 0;
      for (int i = firstParticle; i <= lastParticle; ++i) {
        moments += part[i];
        boundingBox.grow(part[i].position);
        if (part[i].rung > rungs) rungs = part[i].rung;
      }
      calculateRadiusFarthestParticle(moments, &part[firstParticle], &part[lastParticle+1]);
    }

    inline void makeEmpty() {
      myType = Empty;
      particleCount = 0;
      moments.clear();
      boundingBox.reset();
    }

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
#if INTERLIST_VER > 0
      p | bucketListIndex;
      p | startBucket;
      p | visitedR;
      p | visitedL;
#endif
    }
  };
  
  /** A table of the nodes in my tree, indexed by their keys.
      @todo XXX: Make this lookup a hash table, so we get O(1) behavior instead of O(log N).
  */
  typedef std::map<NodeKey, GenericTreeNode *> NodeLookupType;

  class BinaryTreeNode : public GenericTreeNode {//, public CkPool<BinaryTreeNode, 32> {
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

    void fullyDelete() {
      if (children[0] != NULL) children[0]->fullyDelete();
      delete children[0];
      if (children[1] != NULL) children[1]->fullyDelete();
      delete children[1];
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

    int getLevel(NodeKey k) {
      int i = 0;
      k >>= 1;
      while (k!=0) {
	k >>= 1;
	++i;
      }
      return i;
    }

    int whichChild(NodeKey child) {
      int thisLevel = getLevel(key);
      int childLevel = getLevel(child);
      return ((child >> (childLevel-thisLevel-1)) ^ (key << 1));
    }
 
#if INTERLIST_VER > 0
    bool contains(NodeKey node) {
      int thisLevel = getLevel(key);
      int nodeLevel = getLevel(node);
      
      if(nodeLevel<thisLevel) return false;
      
      return ((node>>(nodeLevel-thisLevel))==key);
    }
#endif

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
	    children[0]->makeEmpty();
	    children[1]->myType = lastParticle==totalPart+1 ? Boundary : Internal;
	  }
	  children[0]->lastParticle = firstParticle-1;
          children[0]->particleCount = 0;
	  children[1]->firstParticle = firstParticle;
	  children[1]->particleCount = particleCount;
	} else {
	  // right child missing
	  if (lastParticle == totalPart+1) {
	    children[1]->myType = NonLocal;
	    children[0]->myType = Boundary;
	  } else {
	    children[1]->makeEmpty();
	    children[0]->myType = firstParticle==0 ? Boundary : Internal;
	  }
	  children[1]->firstParticle = lastParticle+1;
          children[1]->particleCount = 0;
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
	  GravityParticle *splitParticle = std::lower_bound(&part[firstParticle],
			&part[lastParticle+1],
			(GravityParticle)(part[lastParticle].key
					  & ((~ (SFC::Key)0) << (62-level))));
	  children[0]->lastParticle = splitParticle - part - 1;
	  children[1]->firstParticle = splitParticle - part;
	  children[0]->particleCount = children[0]->lastParticle - firstParticle + 1;
	  children[1]->particleCount = lastParticle - children[1]->firstParticle + 1;
	  children[0]->myType = children[0]->particleCount==0 ? Empty : (firstParticle==0 ? Boundary : Internal);
	  children[1]->myType = children[1]->particleCount==0 ? Empty : (lastParticle==totalPart+1 ? Boundary : Internal);
          if (children[0]->myType==Empty) children[0]->makeEmpty();
          if (children[1]->myType==Empty) children[1]->makeEmpty();
	  if (firstParticle == 0) children[0]->particleCount --;
	  if (lastParticle == totalPart+1) children[1]->particleCount --;
	}
      } else {
	CkAbort("Ok, the particles should be ordered so, how the hell did we get here?!?");
      }
      children[0]->particlePointer = &part[children[0]->firstParticle];
      children[1]->particlePointer = &part[children[1]->firstParticle];
    }

    // Constructs 2 children of this node based on ORB decomposition
    void makeOrbChildren(GravityParticle *part, int totalPart, int level, int rootsLevel, bool (*compFnPtr[])(GravityParticle, GravityParticle)) {

      children[0] = new BinaryTreeNode();
      children[1] = new BinaryTreeNode();
      children[0]->parent = this;
      children[1]->parent = this;
      children[0]->boundingBox = boundingBox;
      children[1]->boundingBox = boundingBox;

      children[0]->key = key << 1;
      children[1]->key = (key << 1) + 1;
      children[0]->firstParticle = firstParticle;
      children[1]->lastParticle = lastParticle;
    
      //CkPrintf("children keys:%lld,%lld\n",children[0]->key,children[1]->key);
      if(level<rootsLevel){
        //This branch is taken for levels above the TreePiece root level
        //TreePiece root level is the level from where onwards data is internal
        SFC::Key tmp = 1 << (level+1);
        tmp = children[0]->key - tmp;
        SFC::Key tmp2 = part[0].key >> (rootsLevel-level-1);

        if(level+1 != rootsLevel){
          if(tmp==tmp2){
            children[0]->myType = Boundary;
            children[1]->myType = NonLocal;
	          children[1]->firstParticle = lastParticle+1;
	          children[0]->lastParticle = lastParticle;
	          children[0]->particleCount = particleCount;
          }
          else{
            children[0]->myType = NonLocal;
            children[1]->myType = Boundary;
	          children[0]->lastParticle = firstParticle-1;
	          children[1]->firstParticle = firstParticle;
	          children[1]->particleCount = particleCount;
          }
        }
        else{ //Special case for TreePiece root level
          if(tmp==tmp2){
            children[0]->myType = Internal;
            children[1]->myType = NonLocal;
	          children[1]->firstParticle = lastParticle+1;
	          children[0]->lastParticle = lastParticle;
	          children[0]->particleCount = particleCount;
	          if(firstParticle==0)
              children[0]->firstParticle = firstParticle+1;
	          if(lastParticle==totalPart+1)
              children[0]->lastParticle = lastParticle-1;
          }
          else{
            children[0]->myType = NonLocal;
            children[1]->myType = Internal;
	          children[0]->lastParticle = firstParticle-1;
	          children[1]->firstParticle = firstParticle;
	          children[1]->particleCount = particleCount;
	          if(firstParticle==0)
              children[1]->firstParticle = firstParticle+1;
	          if(lastParticle==totalPart+1)
              children[1]->lastParticle = lastParticle-1;
          }
        }
      }
      else{ //Below the TreePiece root level
        float len=0.0,len2=0.0;
        char dim;
				
        len=boundingBox.greater_corner.x-boundingBox.lesser_corner.x;
        dim=0;
        CkAssert(len>=0.0);
        len2=boundingBox.greater_corner.y-boundingBox.lesser_corner.y;
        CkAssert(len2>=0.0);
        if(len2>len) { len = len2; dim=1; }
        len2=boundingBox.greater_corner.z-boundingBox.lesser_corner.z;
        CkAssert(len2>=0.0);
        if(len2>len) { len = len2; dim=2; }
        
        //Sort the particles in longest dimension
        std::sort(&part[firstParticle],&part[lastParticle+1],compFnPtr[dim]);

        //Middle particle is the splitting point to get 2 children

        //Compress the bounding box of this node
        boundingBox.greater_corner[dim] = part[lastParticle].position[dim];
        boundingBox.lesser_corner[dim] = part[firstParticle].position[dim];

        if((lastParticle-firstParticle+1)%2==0){
          children[0]->lastParticle = firstParticle + (lastParticle-firstParticle-1)/2;
          children[1]->firstParticle = firstParticle + (lastParticle-firstParticle+1)/2;
          children[0]->particleCount = (lastParticle-firstParticle+1)/2;
          children[1]->particleCount = (lastParticle-firstParticle+1)/2;
        }
        else{
          children[0]->lastParticle = firstParticle + (lastParticle-firstParticle)/2;
          children[1]->firstParticle = firstParticle + (lastParticle-firstParticle+2)/2;
          children[0]->particleCount = (lastParticle-firstParticle+2)/2;
          children[1]->particleCount = (lastParticle-firstParticle)/2;
        }
				
        children[0]->myType = Internal;
        children[1]->myType = Internal;
       
        //Compress the bounding boxes of both children too
        children[0]->boundingBox.lesser_corner = boundingBox.lesser_corner;
        children[0]->boundingBox.greater_corner = boundingBox.greater_corner;
        children[0]->boundingBox.greater_corner[dim] = part[children[0]->lastParticle].position[dim];
        
        children[1]->boundingBox.lesser_corner = boundingBox.lesser_corner;
        children[1]->boundingBox.lesser_corner[dim] = part[children[1]->firstParticle].position[dim];
        children[1]->boundingBox.greater_corner = boundingBox.greater_corner;
        
      }
      children[0]->particlePointer = &part[children[0]->firstParticle];
      children[1]->particlePointer = &part[children[1]->firstParticle];
      
    }

    // implemented in the .C
    GenericTreeNode *createNew() const;/* {
      return new BinaryTreeNode();
      }*/

    // implemented in the .C
    GenericTreeNode *clone() const;

    /*
    int getNumChunks(int req) {
      int i = 0;
      int num = req;
      while (num > 1) {
	num >>= 1;
	i++;
      }
      int ret = 1 << i;
      if (ret != req) ret <<= 1;
      return ret;
    }
    */

    void getChunks(int num, NodeKey *&ret) {
      int i = 0;
      int base = num;
      while (base > 1) {
	base >>= 1;
	i++;
      }
      base = 1 << i;
      // base now contains the greatest power of two less or equal to num
      int additional = num - base;
      if (ret==NULL) ret = new NodeKey[num];
      for (int j=additional, k=additional*2; j<base; ++j, ++k) ret[k] = j + base;
      base <<= 1;
      for (int j=0; j<additional*2; ++j) ret[j] = j + base;
    }

    int countDepth(int depth) {
      int count = 1;
      if (depth != 0) {
        if (children[0] != NULL) count += children[0]->countDepth(depth-1);
        if (children[1] != NULL) count += children[1]->countDepth(depth-1);
      }
      return count;
    }
    
    int packNodes(BinaryTreeNode *buffer, int depth) {
      //CkPrintf("Entering packNodes: this=%p, buffer=%p, depth=%d\n",this,buffer,depth);
      *buffer = *this;
      buffer->parent = NULL;
      buffer->particlePointer = NULL;
      int used = 1;
      if (depth != 0) {
        if (children[0] != NULL) {
          buffer->children[0] = (BinaryTreeNode*)(((char*)&buffer[used]) - ((char*)buffer));
          //CkPrintf("Entering child 0: offset %ld\n",buffer->children[0]);
          used += children[0]->packNodes(&buffer[used], depth-1);
        } else {
          //CkPrintf("Excluding child 0\n");
          buffer->children[0] = NULL;
        }
        if (children[1] != NULL) {
          buffer->children[1] = (BinaryTreeNode*)(((char*)&buffer[used]) - ((char*)buffer));
          //CkPrintf("Entering child 1: offset %ld\n",buffer->children[1]);
          used += children[1]->packNodes(&buffer[used], depth-1);
        } else {
          //CkPrintf("Excluding child 1\n");
          buffer->children[1] = NULL;
        }
      } else {
        //CkPrintf("Depth reached\n");
        buffer->children[0] = buffer->children[1] = NULL;
      }
      //CkAssert((long int)buffer->children[0] < 10000);
      //CkAssert((long int)buffer->children[1] < 10000);
      //CkPrintf("Returning used = %d\n",used);
      return used;
    }

    void unpackNodes() {
      if (children[0] != NULL) {
        //CkAssert((long int)children[0] < 10000);
        children[0] = (BinaryTreeNode*)(((long int)children[0]) + ((char*)this));
        children[0]->parent = this;
        children[0]->unpackNodes();
      }
      if (children[1] != NULL) {
        //CkAssert((long int)children[1] < 10000);
        children[1] = (BinaryTreeNode*)(((long int)children[1]) + ((char*)this));
        children[1]->parent = this;
        children[1]->unpackNodes();
      }
    }
    
    void pup(PUP::er &p) { pup(p, -1); }
    void pup(PUP::er &p, int depth);/* {
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
      }*/
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
    
    void fullyDelete() {
      if (children[0] != NULL) children[0]->fullyDelete();
      delete children[0];
      if (children[1] != NULL) children[1]->fullyDelete();
      delete children[1];
      if (children[2] != NULL) children[2]->fullyDelete();
      delete children[2];
      if (children[3] != NULL) children[3]->fullyDelete();
      delete children[3];
      if (children[4] != NULL) children[4]->fullyDelete();
      delete children[4];
      if (children[5] != NULL) children[5]->fullyDelete();
      delete children[5];
      if (children[6] != NULL) children[6]->fullyDelete();
      delete children[6];
      if (children[7] != NULL) children[7]->fullyDelete();
      delete children[7];
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
   
#if INTERLIST_VER > 0
    bool contains(NodeKey node) {
      
      return true;
    }
#endif
    
    void makeOctChildren(GravityParticle *part, int totalPart, int level) {}

    void makeOrbChildren(GravityParticle *part, int totalPart, int level, int rootsLevel, bool (*compFnPtr[])(GravityParticle, GravityParticle)) {}

    GenericTreeNode *createNew() const {
      return new OctTreeNode();
    }

    GenericTreeNode *clone() const {
      OctTreeNode *tmp = new OctTreeNode();
      *tmp = *this;
      return tmp;
    }

    /*
    int getNumChunks(int num) {
      int i = 0;
      while (num > 1) {
	num >>= 3;
	i++;
      }
      return 1 << (3*i);
    }
    */

    void getChunks(int num, NodeKey *&ret) {
      return;
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

  /** Added the weight balancer routine*/

  class compare{ //Defines the comparison operator on the map used in balancer
  public:
    compare(){}
    
    bool operator()(NodeKey key1, NodeKey key2) const {

      NodeKey tmp = NodeKey(1);
      int len1=0, len2=0;
      int cnt=1;
  
      while(tmp<=key1 || tmp<=key2){
        tmp<<=1;
        cnt++;
        if(len1==0 && tmp>key1){
          len1=cnt-1;
        }
        if(len2==0 && tmp>key2){
          len2=cnt-1;
        }
      }
  
      if(len1==len2){
        return key1<key2;
      }
      else if(len1>len2){
        key1>>=(len1-len2);
        return key1<key2;
      }
      else{
        key2>>=(len2-len1);
        return key1<key2;
      }
    }
  };

template <typename T>
  class NodeJoint : public CkPool<NodeJoint<T>, 256> {
  public:
  T weight;
  T left;
  T right;
  bool isLeaf;

  NodeJoint(T w, T l, T r) : weight(w), left(l), right(r), isLeaf(false) { }
  NodeJoint(T w) : weight(w), left(0), right(0), isLeaf(true) { }
};

template <typename T>
  class WeightKey {
  public:
  T weight;
  NodeKey key;

  WeightKey(T w, NodeKey k) : weight(w), key(k) {}

  bool operator<(const WeightKey& k2) const {
    if (weight != k2.weight) return weight < k2.weight;
    else return key < k2.key;
  }
};

class NodeKeyClass {
  NodeKey key;
 public:
  NodeKeyClass(NodeKey k) : key(k) {}
  operator NodeKey() const { return key; }
  
  CkHashCode hash(void) const;
  static CkHashCode staticHash(const void *a,size_t);
  int compare(const NodeKeyClass &ind) const;
  static int staticCompare(const void *a,const void *b,size_t);
};

inline CkHashCode NodeKeyClass::hash() const {
  CkHashCode ret = (key >> 32) + key;
  return ret;
}

inline int NodeKeyClass::compare(const NodeKeyClass &ind) const {
  if (key == ind.key) return 1;
  else return 0;
}

template <typename T>
inline void reorderList(NodeKey *nodeKeys, int num, CkVec<T> *zeros, CkVec<NodeKey> *openedNodes, CkHashtableT<NodeKeyClass,NodeJoint<T>*> &nodes) {
  // reorder the output list
  NodeKey current = NodeKey(1);
  NodeJoint<T> *curNode;
  int count = 0;
  if (zeros != NULL) zeros->removeAll();
  if (openedNodes != NULL) openedNodes->removeAll();
  while (1) {
    curNode = nodes.get(current);
    if (curNode != 0) {
      if (curNode->isLeaf || (curNode->left == 0 && curNode->right == 0)) {
        if (curNode->isLeaf) {
          //if (curNode->weight == 0 && zeros != NULL) zeros->push_back(count);
          if (zeros != NULL) zeros->push_back(curNode->weight);
          nodeKeys[count++] = current;
        }
        else if (curNode->weight != 0) {
          // node is not a leaf and has both children weight=0 (while not being itself zero)
          // this means the nodes has just been opened
          nodeKeys[count++] = current << 1;
          nodeKeys[count++] = (current << 1) + 1;
          if (openedNodes != NULL) openedNodes->push_back(current);
        }
        // terminate recursion and get the next sibling
        while (current & 1) current >>= 1;
        if (current == 0) break;
        current++;
      } else {
        // get the first child
        current <<= 1;
      }

      //delete curNode;
    } else {
      // get the first child
      current <<= 1;
    }
  }
  if (count != num) CkPrintf("count = %d, num = %d\n",count,num);
  CkAssert(count == num);
}

template <typename T>
class WeightBalanceState {
 public:
  int realNum; // the total number of keys dividing the domain (counting those which are zero)
  std::set<WeightKey<T> > heaviest;
  std::set<WeightKey<T> > joints;
  CkHashtableT<NodeKeyClass, NodeJoint<T>*> nodes;

  ~WeightBalanceState() {
    NodeJoint<T>::destroyAll();
  }
};

/**
   Function that given a list of NodeKeys and weights associated to them, it
   tries to create a new set of keys with a better distribution of the weights.
   It is used by both Oct domain decomposition and by the CacheManager for the
   chunks.

   When working with persistence (state!=NULL), the first time the function is
   called it expects a complete cut of leaves of the tree. Each subsequent time,
   it espects only the weights of the children of the nodes that have just been
   opened (and therefore had weight zero).

   @param nodeKeys array of keys coming into the function, this array will be
   modified to get new set
   @param weights the weights associated to each key
   @param num the number of input keys
   @param keysTot the size of the nodeKeys array (>=num)
   @param desired how many keys will be present in the output
   @param zeros output a list of indices of nodeKeys that have weight 0. Can be
   null and ignored. If not null "desired" does not consider the keys that are
   NULL
   @param openedNodes a list of all the nodes that have been opened in the execution
   @param state if not NULL the function is working in multiple subsequent calls
   @return true if the keys have been changed, false otherwise
*/
template <typename T>
  inline bool weightBalance(NodeKey *&nodeKeys, T* weights, int &num, int &keysTot, int desired, CkVec<T> *zeros=NULL, CkVec<NodeKey> *openedNodes=NULL, WeightBalanceState<T> *state=NULL) {
  // T can be signed or unsigned
  // if (verbosity>=3) CkPrintf("starting weightBalance iteration\n");

  bool keep = (state != NULL);
  bool isReturning = (state != NULL && state->nodes.numObjects()>0);
  int currentNum; // the number of keys that are currently non zero

  //std::set<WeightKey<T> > heaviest;
  //std::set<WeightKey<T> > joints;
  //CkHashtableT<NodeKeyClass,NodeJoint<T>*> nodes;
  if (state == NULL) {
    state = new WeightBalanceState<T>();
  }

  /*
  if (zeros != NULL) {
    for (int i=0; i<zeros->size(); ++i) {
      nodes.put((*zeros)[i]) = new NodeJoint<T>(0,0,0);
    }
  }
  */

  if (isReturning) {
    currentNum = desired;
    CkAssert (openedNodes != NULL);
    for (int i=0; i<openedNodes->size(); ++i) {
      NodeKey opened = openedNodes->operator[](i);
      NodeKey left = opened << 1;
      NodeKey right = left + 1;
      state->heaviest.insert(WeightKey<T>(weights[i*2],left));
      NodeJoint<T> *node = state->nodes.get(left);
      if (node) {
        node->weight = weights[i*2];
        node->isLeaf = true;
      }
      else state->nodes.put(left) = new NodeJoint<T>(weights[i*2]);
      state->heaviest.insert(WeightKey<T>(weights[i*2+1],right));
      node = state->nodes.get(right);
      if (node) {
        node->weight = weights[i*2+1];
        node->isLeaf = true;
      }
      else state->nodes.put(right) = new NodeJoint<T>(weights[i*2+1]);
      node = state->nodes.getRef(opened);
      node->right = weights[i*2];
      node->left = weights[i*2+1];
      if (weights[i*2] == 0) currentNum--;
      if (weights[i*2+1] == 0) currentNum--;
    }
  } else {
    currentNum = 0;
    state->realNum = num;
    for (int i=0; i<num; ++i) {
      state->heaviest.insert(WeightKey<T>(weights[i],nodeKeys[i]));
      state->nodes.put(nodeKeys[i]) = new NodeJoint<T>(weights[i]);
      if (weights[i] > 0 || zeros == NULL) currentNum++;

      // check if we can have a joint
      NodeKey siblingKey(nodeKeys[i] ^ 1);
      NodeJoint<T> *sibling = state->nodes.get(siblingKey);
      if (sibling != 0) {
        state->joints.insert(WeightKey<T>(weights[i]+sibling->weight,siblingKey>>1));
        T left = (siblingKey&1) ? weights[i] : sibling->weight;
        T right = (siblingKey&1) ? sibling->weight : weights[i];
        state->nodes.put(siblingKey>>1) = new NodeJoint<T>(weights[i]+sibling->weight,left,right);
      }
    }
  }

  /*
  for (int i=0; i<num; ++i) {
    //heaviest.insert(std::pair<T,NodeKey>(weights[i],nodeKeys[i]));
    heaviest.insert(WeightKey<T>(weights[i],nodeKeys[i]));
    nodes.put(nodeKeys[i]) = new NodeJoint<T>(weights[i]);
    NodeKey siblingKey(nodeKeys[i] ^ 1);
    NodeJoint<T> *sibling = nodes.get(siblingKey);
    NodeKey ownKey(nodeKeys[i]);
    T ownWeight = weights[i];
    while (sibling != 0) {
      //CkPrintf("Found possible junction into %llx (w=%d)\n",sibling>>1,siblingWeight->weight);
      if (zeros != NULL && (sibling->weight == 0 || ownWeight == 0)) {
        if (weights[i] != 0) {
          CkPrintf("Zero found1 %llx (%llx: %d)\n",siblingKey,ownKey,ownWeight);
          zeros->push_back(siblingKey);
          heaviest.erase(WeightKey<T>(ownWeight,ownKey));
          heaviest.erase(WeightKey<T>(sibling->weight,siblingKey));
          sibling->isLeaf = false;
          nodes.getRef(nodeKeys[i])->isLeaf = false;
          joints.insert(WeightKey<T>(ownWeight,ownKey));
        } else if (sibling->weight != 0) {
          CkPrintf("Zero found2 %llx (%llx: %d)\n",ownKey,siblingKey,sibling->weight);
          zeros->push_back(ownKey);
          heaviest.erase(WeightKey<T>(ownWeight,ownKey));
          heaviest.erase(WeightKey<T>(sibling->weight,siblingKey));
          sibling->isLeaf = false;
          nodes.getRef(nodeKeys[i])->isLeaf = false;
          joints.insert(WeightKey<T>(sibling->weight,siblingKey));
        } else { // both weight[i] and siblingWeight->weight are zero
          CkPrintf("Zero found3 %llx %llx\n",ownKey,siblingKey);
          // nothing is added to the zeros vector since the fact that both are
          // zero implies that their parent is zero too, and therefore the
          // parent will eventually added to this list
          heaviest.erase(WeightKey<T>(ownWeight,ownKey));
          heaviest.erase(WeightKey<T>(sibling->weight,siblingKey));
          sibling->isLeaf = false;
          nodes.getRef(nodeKeys[i])->isLeaf = false;
          heaviest.insert(WeightKey<T>(0,siblingKey>>1));
          nodes.put(siblingKey>>1) = new NodeJoint<T>(0);
          // prepare for next while loop
          ownKey = siblingKey>>1;
          siblingKey = ownKey ^ 1;
          sibling = nodes.get(siblingKey);
        }
      } else {
        joints.insert(WeightKey<T>(ownWeight+sibling->weight,siblingKey>>1));
        T left = (siblingKey&1) ? ownWeight : sibling->weight;
        T right = (siblingKey&1) ? sibling->weight : ownWeight;
        nodes.put(siblingKey>>1) = new NodeJoint<T>(ownWeight+sibling->weight,left,right);
      }
    }
  }
  */


  // at this point heaviest is a sorted list of the the leaves of the tree;
  // nodes contains all the leaves too, but in a hash table;
  // joints contains all the possible conjunctions of two leaves of the tree (again sorted).

  bool result = false;
  int oldNum;
  //reorderList<T>(nodeKeys, num, nodes);
  while (state->heaviest.rbegin()->weight > state->joints.begin()->weight ||
         currentNum != desired) {
    if (state->joints.empty()) {
      if (currentNum == desired) break;
      else if (currentNum > desired) CmiAbort("weightBalance: possible joints set emply, while convergence not reached!");
    }
    result = true;
    NodeKey openedKey = state->heaviest.rbegin()->key;
    T openedWeight = state->heaviest.rbegin()->weight;
    NodeKey closedKey = state->joints.begin()->key;
    T closedWeight = state->joints.begin()->weight;

    oldNum = currentNum;
    if (oldNum <= desired) {
      // if (verbosity>=3) CkPrintf("[%d] Opening node %llx (%d)\n",CkMyPe(),state->heaviest.rbegin()->key,state->heaviest.rbegin()->weight);
      NodeJoint<T> *opened = state->nodes.getRef(openedKey);
      state->heaviest.erase(WeightKey<T>(openedWeight,openedKey));
      opened->isLeaf = false;
      state->joints.insert(WeightKey<T>(openedWeight,openedKey));
      NodeJoint<T> *parent = state->nodes.get(openedKey>>1);
      if (parent) state->joints.erase(WeightKey<T>(parent->weight,openedKey>>1));
      NodeJoint<T> *child = state->nodes.get(openedKey<<1);
      if (child) {
        child->isLeaf = true;
        state->heaviest.insert(WeightKey<T>(child->weight,openedKey<<1));
      }
      child = state->nodes.get((openedKey<<1) + 1);
      if (child) {
        child->isLeaf = true;
        state->heaviest.insert(WeightKey<T>(child->weight,(openedKey<<1)+1));
      }
      // if both children are non-empty (or if they just don't exist), we have added a new key
      if (zeros == NULL || (opened->left > 0 && opened->right > 0) || !child) currentNum++;
      state->realNum++;

    }
    if (oldNum >= desired) {
      // if (verbosity>=3) CkPrintf("[%d] Closing node %llx (%d)\n",CkMyPe(),state->joints.begin()->key,state->joints.begin()->weight);
      NodeJoint<T> *closed = state->nodes.getRef(closedKey);
      state->joints.erase(state->joints.begin());
      closed->isLeaf = true;
      state->heaviest.insert(WeightKey<T>(closedWeight,closedKey));
      NodeJoint<T> *child = state->nodes.get(closedKey<<1);
      if (child) {
        child->isLeaf = false;
        state->heaviest.erase(WeightKey<T>(closed->left,closedKey<<1));
      }
      child = state->nodes.get((closedKey<<1)+1);
      if (child) {
        child->isLeaf = false;
        state->heaviest.erase(WeightKey<T>(closed->right,(closedKey<<1)+1));
      }
      if (zeros == NULL || (closed->left > 0 && closed->right > 0) || !child) currentNum--;
      state->realNum--;

      // find if there is a new joint
      NodeJoint<T> *sibling = state->nodes.get(closedKey ^ 1);
      if (sibling != 0 && sibling->isLeaf) {
        NodeJoint<T> *parent = state->nodes.get(closedKey>>1);
        if (parent == 0) {
          T left = (closedKey&1) ? sibling->weight : closedWeight;
          T right = (closedKey&1) ? closedWeight : sibling->weight;
          parent = new NodeJoint<T>(closedWeight+sibling->weight,left,right);
          state->nodes.put(closedKey>>1) = parent;
        }
        state->joints.insert(WeightKey<T>(parent->weight,closedKey>>1));
      }

    }
    /*
    else { // normal case of weights unbalanced

    CkPrintf("[%d] Opening node %llx (%d), closing %llx (%d)\n",CkMyPe(),heaviest.rbegin()->key,heaviest.rbegin()->weight,joints.begin()->key,joints.begin()->weight);
    // we can improve the balancing!
    result = true;

    NodeKey openedKey = heaviest.rbegin()->key;
    T openedWeight = heaviest.rbegin()->weight;
    NodeKey closedKey = joints.begin()->key;
    T closedWeight = joints.begin()->weight;
    NodeJoint<T> *opened = nodes.getRef(openedKey);
    NodeJoint<T> *closed = nodes.getRef(closedKey);

    heaviest.erase(WeightKey<T>(openedWeight,openedKey));
    joints.erase(joints.begin());

    CkAssert(closed != 0);
    closed->isLeaf = true;
    NodeJoint<T> *child = nodes.get(closedKey<<1);
    if (child != 0) child->isLeaf = false;
    child = nodes.get((closedKey<<1)+1);
    if (child != 0) child->isLeaf = false;
    heaviest.erase(WeightKey<T>(closed->left,closedKey<<1));
    heaviest.erase(WeightKey<T>(closed->right,(closedKey<<1)+1));
    heaviest.insert(WeightKey<T>(closedWeight,closedKey));

    NodeJoint<T> *sibling = nodes.get(closedKey ^ 1);
    if (sibling != 0 && (sibling->isLeaf || (zeros != NULL && sibling->weight == 0))) {
      NodeJoint<T> *parent = nodes.get(closedKey>>1);
      if (parent == 0) {
        T left = (closedKey&1) ? sibling->weight : closedWeight;
	T right = (closedKey&1) ? closedWeight : sibling->weight;
	parent = new NodeJoint<T>(closedWeight+sibling->weight,left,right);
	nodes.put(closedKey>>1) = parent;
      }
      joints.insert(WeightKey<T>(parent->weight,closedKey>>1));

      if (zeros != NULL && sibling->weight == 0) {
	// found an empty sibling, joining directly to the parent
	// also, delete the empty node from the zero list
	CkPrintf("[%d] WeightBalance: found a node whose sibling is empty!\n",CkMyPe());
	int ptr = 0;
	NodeKey siblingKey = closedKey ^ 1;
	while (ptr < zeros->size() && zeros->operator[](ptr)!=siblingKey) ptr++;
        if (ptr == zeros->size()) {
          CkPrintf("requesting for zero node %llx\n",siblingKey);
          for (int i=0; i<zeros->size(); ++i) CkPrintf("  element %i: %llx\n",i,zeros->operator[](i));
        }
	CkAssert(ptr < zeros->size());
	zeros->remove(ptr);

	closed->isLeaf = false;
	heaviest.erase(WeightKey<T>(closedWeight,closedKey));
	joints.erase(WeightKey<T>(parent->weight,closedKey>>1));
	parent->isLeaf = true;
	heaviest.insert(WeightKey<T>(parent->weight,closedKey>>1));
	// recursively join to the next parent level.... TODO
      }
    }

    CkAssert(opened != 0);
    opened->isLeaf = false;
    joints.insert(WeightKey<T>(openedWeight,openedKey));
    if (opened->left != 0 || opened->right != 0) {
      // the children already exist and they have full counting information
      NodeJoint<T> *openedLeft = nodes.getRef(openedKey<<1);
      NodeJoint<T> *openedRight = nodes.getRef((openedKey<<1)+1);
      openedLeft->isLeaf = true;
      openedRight->isLeaf = true;
      heaviest.insert(WeightKey<T>(openedKey<<1,opened->left));
      heaviest.insert(WeightKey<T>((openedKey<<1)+1,opened->right));

      if (zeros != NULL && (opened->left == 0 || opened->right == 0)) {
	// found a node which is empty, deal with it
	CkPrintf("[%d] WeightBalance: found a node with an empty child!\n",CkMyPe());
      }
    }/ * else if (nodes.get(openedKey<<1) == 0) {
      nodes.put(openedKey<<1) = NodeJoint<T>(0);
      nodes.put((openedKey<<1)+1) = NodeJoint<T>(0);
    }* /
    //reorderList<T>(nodeKeys, num, nodes);
    }
    */
  }

  num = state->realNum;
  if (num > keysTot) {
    delete[] nodeKeys;
    keysTot = num+1;
    nodeKeys = new NodeKey[keysTot];
  }
  reorderList<T>(nodeKeys, num, zeros, openedNodes, state->nodes);
  if (!keep) delete state;
  // if (verbosity>=2) CkPrintf("weightBalance finished\n");
  return result;
}

template <class T>
  inline bool weightBalance(NodeKey *nodeKeys, T* weights, int num, int handleZero){
  //T can be signed or unsigned
	
  NodeKey curHeaviest;
  typename std::map<NodeKey,T,compare>::iterator curLightest;
  T lightestWt= ~T(0);
  T tmpWt=0;
  NodeKey parent,child1,child2;
  int numBalances=0;

  int zeroHandled=0;
  
    //Need to construct a temporary copy of the input data to operate
    //construct a map indexed by the nodekey
    std::map<NodeKey,T,compare> curNodeWts;
    typename std::map<NodeKey,T,compare>::iterator iter;
    typename std::map<NodeKey,T,compare>::iterator iter2;
    curNodeWts.clear();
    for(int i=0;i<num;i++){
      curNodeWts[nodeKeys[i]]=weights[i];
    }

    //loop here
    while(1){
      tmpWt=0;
      lightestWt=~T(0);
      //find the heaviest Node
      for(iter=curNodeWts.begin();iter!=curNodeWts.end();iter++){
	if((*iter).second>tmpWt && (*iter).second!=~T(0)){
	  tmpWt=(*iter).second;
	  curHeaviest=(*iter).first;
	}
      }
      if(tmpWt==0) //In case, no-one had weight > 0
	break;
    
      //find the lightest parent-- implemented only for a binary tree
      iter=curNodeWts.begin();
      iter2=curNodeWts.begin();
      iter2++;
      for( ;iter2!=curNodeWts.end();iter++,iter2++){
	if((*iter).second==~T(0) || (*iter2).second==~T(0))//Ignore those which have been opened
	  continue;
      
	if(handleZero){
	  if((*iter).second==0 || (*iter2).second==0){
	    if((*iter).second==0){
	      curNodeWts.erase(iter);
	    }
	    else if((*iter2).second==0){
	      curNodeWts.erase(iter2);
	    }
          
	    numBalances++;
	    //Open the current heaviest and continue
	    child1=curHeaviest << 1;
	    child2=child1 | NodeKey(1);
	    //Erase the heaviest and add it's two children
	    curNodeWts.erase(curHeaviest);
	    curNodeWts[child1]=~T(0);
	    curNodeWts[child2]=~T(0);
	    zeroHandled=1;
	    break;
	  }
	}
      
	if((*iter).first==curHeaviest || (*iter2).first==curHeaviest)
	  continue;
	child1=(*iter).first;
	child2=(*(iter2)).first;
	child1 >>= 1;
	child2 >>= 1;
	if(child1==child2){
	  tmpWt=(*iter).second+(*iter2).second;
	  if(tmpWt<lightestWt || lightestWt==~T(0)){
	    lightestWt=tmpWt;
	    curLightest=iter;
	  }
	}
      }

      if(handleZero && zeroHandled){
	zeroHandled=0;
	continue;
      }
    
      //balance only if the heaviest is heavier than the lightest possible parent
      if((curNodeWts[curHeaviest] > lightestWt) && lightestWt!=~T(0)){
	numBalances++;
	parent = (*curLightest).first >> 1;
	iter2=curLightest; iter2++; iter2++;
	//Erase the children and add the lightest parent
	curNodeWts.erase(curLightest,iter2);
	//curNodeWts[parent]=lightestWt;
	curNodeWts.insert(std::pair<NodeKey,T>(parent,lightestWt));
	child1=curHeaviest << 1;
	child2=child1 | NodeKey(1);
	//Erase the heaviest and add it's two children
	curNodeWts.erase(curHeaviest);
	curNodeWts[child1]=~T(0);
	curNodeWts[child2]=~T(0);
      }
      else //We are done here
	break;
    }
    //end loop here

    int i=0;
    //construct new node key array before returning
    for(iter=curNodeWts.begin(),i=0;iter!=curNodeWts.end();i++,iter++){
      nodeKeys[i]=(*iter).first;
    }
    if(i!=num)
      CkPrintf("i:%d,num:%d\n",i,num);
    CkAssert(i==num);
  
    if(numBalances>0){
      return true;
    }
    else {
      return false;
    }

  }

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
