/** \file GenericTreeNode.h
    This file defines the generic tree structures.
    @author Filippo Gioachin (previously Graeme Lufkin and Chao Huang)
    @version 2.0
*/

#ifndef GENERICTREENODE_H
#define GENERICTREENODE_H

#include "pup.h"

#include <map>
#include <vector>
#include <algorithm>
#include <set>
#include <list>

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
#ifdef BIGKEYS
  typedef __uint128_t NodeKey;
#else
  typedef CmiUInt8 NodeKey;
#endif
  static const int NodeKeyBits = 8*sizeof(NodeKey);

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

class NodePool;

  class GenericTreeNode {
#ifdef CHANGA_REFACTOR_WALKCHECK
    public:
    bool touched;
    char by;
#endif
  protected:
    NodeType myType;
    NodeKey key;
    CmiUInt8 usedBy;

    GenericTreeNode() : myType(Invalid), key(0), parent(0), firstParticle(0),
	lastParticle(0), remoteIndex(0), usedBy(0), iParticleTypes(0) {
#if COSMO_STATS > 0
      used = false;
#endif
#if INTERLIST_VER > 0
      numBucketsBeneath=0;
      startBucket=-1;
#ifdef CUDA
      nodeArrayIndex = -1;
      bucketArrayIndex = -1;
      wasNeg = true;
#endif
#endif
#ifdef CHANGA_REFACTOR_WALKCHECK
      touched = false;
      by = 'I';
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
    OrientedBox<cosmoType> boundingBox;
    /// The bounding box including search balls of this node
    OrientedBox<double> bndBoxBall;
    /// Mask of particle types contatained in this node
    unsigned int iParticleTypes;
    /// An index for the first particle contained by this node, 0 means outside the node
    /// An index for the first particle contained by this node, 0 means outside the node
    int firstParticle;
    /// An index to the last particle contained by this node, myNumParticles+1 means outside the node
    int lastParticle;
    /// An index to the *real* location of this node if this node is NonLocal, if
    /// it is Boundary or Internal or Bucket it is equal to thisIndex
    /// During Treebuid, it indicates whether remote moments are
    /// needed to calculate this nodes moment.
    int remoteIndex;
    /// Total number of particles contained (across all chares)
    unsigned int particleCount;
    /// Pointer to the first particle in this node
    GravityParticle *particlePointer;

    /// The greatest rung amoung the particles contained in this node (greater
    /// means faster). This information is limited to the nodes in the current
    /// TreePiece, and do not consider non-local data.
    int rungs;

#if INTERLIST_VER > 0
    //int bucketListIndex;
    int numBucketsBeneath;
    int startBucket;
#ifdef CUDA
    // index in nodeinfo array
    int nodeArrayIndex;
    int bucketArrayIndex;
    bool wasNeg;
#endif
#endif
    /// center of smoothActive particles during smooth operation
    Vector3D<double> centerSm;
    /// Radius of bounding sphere of smoothActive particles
    double sizeSm;
    /// Maximum smoothing radius of smoothActive particles
    double fKeyMax;
    /// SMP rank of node owner
    int iRank;

    GenericTreeNode(NodeKey k, NodeType type, int first, int last, GenericTreeNode *p) : myType(type), key(k), parent(p), firstParticle(first), lastParticle(last), remoteIndex(0), usedBy(0) {
#if INTERLIST_VER > 0
      numBucketsBeneath=0;
      startBucket=-1;
#ifdef CUDA
      nodeArrayIndex = -1;
      bucketArrayIndex = -1;
      wasNeg = true;
#endif
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

    bool isValid(){
      return (myType != Invalid);
    }

    bool isCached(){
      return (myType == Cached ||
              myType == CachedBucket ||
              myType == CachedEmpty);
    }

    bool isBucket(){
      return (myType == Bucket ||
              myType == CachedBucket ||
              myType == NonLocalBucket);
    }

    // these two functions are used to track the communication between objects:
    // a nodes is marked usedBy when a local TreePiece has touched it
    void markUsedBy(int index) { usedBy |= (((CmiUInt8)1) << index); }
    bool isUsedBy(int index) { return (usedBy & (((CmiUInt8)1) << index)); }

    /// \brief construct the children of the "this" node following the
    /// given logical criteria (Oct/Orb)
    virtual void makeOctChildren(GravityParticle *part, int totalPart, int level, NodePool *pool = NULL) = 0;
    virtual void makeOrbChildren(GravityParticle *part, int totalPart, int level, int rootsLevel, bool (*compFnPtr[])(GravityParticle, GravityParticle), bool spatial, NodePool *pool = NULL) = 0;

    // get the number of chunks possible for the given request
    // @return a number greater or iqual to th the request
    //virtual int getNumChunks(int num) = 0;
    /// get the nodes corresponding to a particular number of chunks requested
    virtual void getChunks(int num, NodeKey *&ret) = 0;

    /// transform an internal node into a bucket
    inline void makeBucket(GravityParticle *part) {
      myType = Bucket;
      iRank = CkMyRank();
#if INTERLIST_VER > 0
      numBucketsBeneath = 1;
#endif
      calculateRadiusBox(moments, boundingBox);	/* set initial size */
      boundingBox.reset();
      bndBoxBall.reset();
      iParticleTypes = 0;
      rungs = 0;
      for (int i = firstParticle; i <= lastParticle; ++i) {
        moments += part[i];
        boundingBox.grow(part[i].position);
	if(TYPETest(&part[i], TYPE_GAS)) {
	    double fBallMax = part[i].fBallMax();
	    bndBoxBall.grow(part[i].position
			    + Vector3D<double>(fBallMax, fBallMax, fBallMax));
	    bndBoxBall.grow(part[i].position
			    - Vector3D<double>(fBallMax, fBallMax, fBallMax));
	    }
	iParticleTypes |= part[i].iType;
        if (part[i].rung > rungs) rungs = part[i].rung;
      }
      if(particleCount > 1)
	  calculateRadiusFarthestParticle(moments, &part[firstParticle],
					  &part[lastParticle+1]);
    }

    inline void makeEmpty() {
      myType = Empty;
      particleCount = 0;
#if INTERLIST_VER > 0
      numBucketsBeneath = 0;
#endif
      moments.clear();
      boundingBox.reset();
      bndBoxBall.reset();
      iParticleTypes = 0;
    }

    void getGraphViz(std::ostream &out);

    /// @brief return the NodeKey of the lowest common ancestor.
    virtual NodeKey getLongestCommonPrefix(NodeKey k1, NodeKey k2)
    {
      CkAbort("getLongestCommonPrefix not implemented\n");
      return 0;
    }

    virtual int getLevel(NodeKey k) = 0;
    virtual GenericTreeNode *createNew() const = 0;
    virtual GenericTreeNode *clone() const = 0;

    virtual void pup(PUP::er &p, int depth) = 0;
    virtual void pup(PUP::er &p) {
      int iType;
      if(p.isUnpacking()) {
        p | iType;
  	    myType = (NodeType) iType;
#ifdef CUDA
  	    // so that newly shipped nodes are not mistakenly assumed
  	    // to be present on the GPU
  	    nodeArrayIndex = -1;
      bucketArrayIndex = -1;
      wasNeg = true;
#endif
      } else {
        iType = (int) myType;
        p | iType;
      }
      p | key;
      p | moments;
      p | boundingBox;
      p | bndBoxBall;
      p | iParticleTypes;
      p | firstParticle;
      p | lastParticle;
      p | remoteIndex;
      p | particleCount;
#if INTERLIST_VER > 0
      p | numBucketsBeneath;
      p | startBucket;
#endif
      p | centerSm;
      p | sizeSm;
      p | fKeyMax;
#ifdef CHANGA_REFACTOR_WALKCHECK
      p | touched;
      p | by;
#endif
    }
  };

class BinaryTreeNode;

/// Utility to pool allocations of tree nodes
class NodePool {
    int next;
    int szPool;
    /// list of allocated pool blocks.
    std::list<BinaryTreeNode *> pools;
public:
    NodePool() {
        next = szPool = 128;  // Determines the size of the pools
    }
    ~NodePool();
    /// give out one node from the pool, allocating a new pool block if needed.
    BinaryTreeNode *alloc_one();
    BinaryTreeNode *alloc_one(NodeKey k, NodeType type, int first,
			      int nextlast, BinaryTreeNode *p);

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

    NodeKey getLongestCommonPrefix(NodeKey k1, NodeKey k2){

      int l1 = getLevel(k1)+1;
      int l2 = getLevel(k2)+1;

      NodeKey a1 = k1 << (NodeKeyBits-l1);
      NodeKey a2 = k2 << (NodeKeyBits-l2);

      NodeKey a = a1 ^ a2;
      int i = 0;
      while( a < (NodeKey(1)<<(NodeKeyBits-1)) )
      {
        a = a << 1;
        i++;
      }
      return (k1 >> (l1-i)); // or k2 >> (l2-i)
    };


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

    /// Equally divide space into two child nodes.  The split
    /// direction is determined by level.
    /// For each child:
    ///   1. A key is assigned which encodes the tree branching path
    /// to the child
    ///	  2. indices to first and last particles are assigned.
    ///   3. "Type" of node is assigned (Internal, Boundary, etc.)
    ///        This generally depends on the nature of the sibling.
    ///   4. Pointer to particles is assigned.
    ///
    void makeOctChildren(GravityParticle *part, int totalPart, int level,
			 NodePool *pool=NULL) {
      if(pool) {
	  children[0] = pool->alloc_one();
	  children[1] = pool->alloc_one();
	  }
      else {
	  children[0] = new BinaryTreeNode();
	  children[1] = new BinaryTreeNode();
	  }
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

      SFC::Key mask = SFC::Key(1) << ((SFC::KeyBits-1) - level);
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
					  & ((~ (SFC::Key)0) << ((SFC::KeyBits-1)-level))));
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
#if INTERLIST_VER > 0
      if(children[0]->myType == NonLocal) children[0]->numBucketsBeneath = 0;
      if(children[1]->myType == NonLocal) children[1]->numBucketsBeneath = 0;
      // empty nodes are makeEmpty()'ed, so that the numbucketsbeneath them are 0
#endif
    }

    // Constructs 2 children of this node based on ORB decomposition
    void makeOrbChildren(GravityParticle *part, int totalPart, int level,
		         int rootsLevel,
			 bool (*compFnPtr[])(GravityParticle, GravityParticle),
			 bool spatial, NodePool *pool=NULL) {

      if(pool) {
	  children[0] = pool->alloc_one();
	  children[1] = pool->alloc_one();
	  }
      else {
          children[0] = new BinaryTreeNode();
          children[1] = new BinaryTreeNode();
	  }
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
        int dim;
	
	if(spatial) { // "squeeze" box before division
	    boundingBox.reset();
	    for (int i = firstParticle; i <= lastParticle; ++i)
		boundingBox.grow(part[i].position);
	    }

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

	if(!spatial) {
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
	}
	else {
	    GravityParticle dummy;  // dummy for splitting
	    GravityParticle* divStart = &part[firstParticle];
	    Vector3D<double> divide(0.0,0.0,0.0);
	    divide[dim] = part[firstParticle].position[dim] + 0.5*len;
	    dummy.position = divide;
	    GravityParticle* divEnd
		= std::upper_bound(&part[firstParticle],&part[lastParticle+1],
				   dummy,compFnPtr[dim]);
	    children[0]->lastParticle = firstParticle + (divEnd - divStart) - 1;
	    children[1]->firstParticle = firstParticle + (divEnd - divStart);
	    children[0]->particleCount = divEnd - divStart;
	    children[1]->particleCount = 1 + (lastParticle-firstParticle)
		- (divEnd - divStart);
	    CkAssert(children[0]->particleCount > 0);
	    CkAssert(children[1]->particleCount > 0);
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

    // extraSpace is used by the CacheInterface to store a pointer to
    // the message so it can be freed when we are finished.
    int packNodes(BinaryTreeNode *buffer, int depth, int extraSpace=0) {
      //CkPrintf("Entering packNodes: this=%p, buffer=%p, depth=%d\n",this,buffer,depth);
      //*buffer = *this;
      memcpy(buffer, this, sizeof(*this));
      buffer->parent = NULL;
      buffer->particlePointer = NULL;
#if INTERLIST_VER > 0 && defined CUDA
      buffer->nodeArrayIndex = -1;
      buffer->bucketArrayIndex = -1;
      buffer->wasNeg = true;
#endif
      int used = 1;
      if (depth != 0) {
        if (children[0] != NULL) {
          BinaryTreeNode *nextBuf = (BinaryTreeNode *) (((char*)buffer) + used * ALIGN_DEFAULT(sizeof(BinaryTreeNode)+extraSpace));
          buffer->children[0] = (BinaryTreeNode*)(((char*)nextBuf) - ((char*)buffer));
          //CkPrintf("Entering child 0: offset %ld\n",buffer->children[0]);
          used += children[0]->packNodes(nextBuf, depth-1, extraSpace);
        } else {
          //CkPrintf("Excluding child 0\n");
          buffer->children[0] = NULL;
        }
        if (children[1] != NULL) {
          BinaryTreeNode *nextBuf = (BinaryTreeNode *) (((char*)buffer) + used * ALIGN_DEFAULT(sizeof(BinaryTreeNode)+extraSpace));
          buffer->children[1] = (BinaryTreeNode*)(((char*)nextBuf) - ((char*)buffer));
          //CkPrintf("Entering child 1: offset %ld\n",buffer->children[1]);
          used += children[1]->packNodes(nextBuf, depth-1, extraSpace);
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

inline
NodePool::~NodePool() {
        while(!pools.empty()) {
            delete [] pools.front();
            pools.pop_front();
            }
        }

inline BinaryTreeNode *
NodePool::alloc_one() {
        if(next >= szPool) { // need new pool block
            BinaryTreeNode *nodes = new BinaryTreeNode[szPool];
            next = 0;
            pools.push_front(nodes);
            }
        return &pools.front()[next++];
        }

inline BinaryTreeNode *
NodePool::alloc_one(NodeKey k, NodeType type, int first, int nextlast,
		    BinaryTreeNode *p) {
	BinaryTreeNode *one = alloc_one();
	return new (one) BinaryTreeNode(k, type, first, nextlast, p);
	}

  class OctTreeNode : public GenericTreeNode {
  protected:
  public:
    /// the 8 different children
    OctTreeNode* children[8];

    OctTreeNode() : GenericTreeNode() {
      for (int i = 0; i < 8; i++) {
        children[i] = 0;
      }
    }

  public:
    OctTreeNode(NodeKey k, NodeType type, int first, int last, BinaryTreeNode *p) : GenericTreeNode(k, type, first, last, p) {
      for (int i = 0; i < 8; i++) {
        children[i] = 0;
      }
    }

    void fullyDelete() {
      for (int i = 0; i < numChildren(); i++) {
        if (children[i] != NULL) {
          children[i]->fullyDelete();
          delete children[i];
        }
      }
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

    int getLevel(NodeKey k) {
      int i = 0;
      k >>= 3;
      while (k!=0) {
	k >>= 3;
	++i;
      }
      return i;
    }

#if INTERLIST_VER > 0
    bool contains(NodeKey node) {

      return true;
    }
#endif

    void makeOctChildren(GravityParticle *part, int totalPart, int level,
			 NodePool *pool = NULL) {}

    void makeOrbChildren(GravityParticle *part, int totalPart, int level, int rootsLevel, bool (*compFnPtr[])(GravityParticle, GravityParticle), bool spatial,
			 NodePool *pool = NULL) {}

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
