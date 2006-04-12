/** \file GenericTreeNode.h
    This file defines the generic tree structures.
    @author Filippo Gioachin (previously Graeme Lufkin and Chao Huang)
    @version 2.0
*/

#ifndef GENERICTREENODE_H
#define GENERICTREENODE_H

#include "pup.h"
#include "ckpool.h"

#include <algorithm>
#include <map>

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
	
    GenericTreeNode() : myType(Invalid), key(0), parent(0), firstParticle(0), lastParticle(0), remoteIndex(0) {
      #if INTERLIST_VER > 0
        visitedR=false;
        visitedL=false;
			  bucketListIndex=-1;
        startBucket=-1;
      #endif
    }

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
    /// Pointer to the first particle in this node
    GravityParticle *particlePointer;

#if INTERLIST_VER > 0
    bool visitedR;
    bool visitedL;
    int bucketListIndex;
    int startBucket;
#endif
    //GenericTreeNode() : myType(Invalid), key(0), parent(0), beginParticle(0), endParticle(0), remoteIndex(0) { }

    GenericTreeNode(NodeKey k, NodeType type, int first, int last, GenericTreeNode *p) : myType(type), key(k), parent(p), firstParticle(first), lastParticle(last), remoteIndex(0) { 
    #if INTERLIST_VER > 0
      visitedR=false;
      visitedL=false;
			bucketListIndex=-1;
      startBucket=-1;
    #endif
    }

    virtual ~GenericTreeNode() { }
    virtual void fullyDelete() = 0;
    
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
#if INTERLIST_VER > 0
    virtual bool contains(NodeKey nodekey) = 0;
#endif
    
    /// construct the children of the "this" node following the given logical
    /// criteria (Oct/Orb)
    virtual void makeOctChildren(GravityParticle *part, int totalPart, int level) = 0;
    virtual void makeOrbChildren(GravityParticle *part, int totalPart, int level) = 0;

    // get the number of chunks possible for the given request
    // @return a number greater or iqual to th the request
    //virtual int getNumChunks(int num) = 0;
    /// get the nodes corresponding to a particular number of chunks requested
    virtual void getChunks(int num, NodeKey *&ret) = 0;

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
      children[0]->particlePointer = &part[children[0]->firstParticle];
      children[1]->particlePointer = &part[children[1]->firstParticle];
    }

    void makeOrbChildren(GravityParticle *part, int totalPart, int level) {}

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

    void makeOrbChildren(GravityParticle *part, int totalPart, int level) {}

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

template <class T>
inline bool weightBalance(NodeKey *nodeKeys, T* weights, int num){
  //T can be signed or unsigned
	
  NodeKey curHeaviest;
  typename std::map<NodeKey,T,compare>::iterator curLightest;
	T lightestWt= ~T(0);
	T tmpWt=0;
	NodeKey parent,child1,child2;
  int numBalances=0;
  
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
  if(numBalances>0){
    //construct new node key array before returning
	  for(iter=curNodeWts.begin(),i=0;iter!=curNodeWts.end();i++,iter++){
      nodeKeys[i]=(*iter).first;
    }
    CkAssert(i==num);
    return true;
  }
  else { return false; }

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
