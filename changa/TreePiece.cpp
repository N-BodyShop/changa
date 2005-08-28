/** @file TreePiece.cpp
 */

#include <cstdio>
#include <algorithm>
#include <fstream>
#include <assert.h>

//#include "ComlibManager.h"

#include "ParallelGravity.h"
#include "CacheManager.h"

#include "Space.h"
#include "gravity.h"

using namespace std;
using namespace SFC;
using namespace TreeStuff;
using namespace TypeHandling;

int TreeStuff::maxBucketSize;

//forward declaration
string getColor(GenericTreeNode*);

void TreePiece::load(const std::string& fn, const CkCallback& cb) {
  basefilename = fn;
	
  //read in particles
  XDR xdrs;
  FILE* infile = fopen((basefilename + ".mass").c_str(), "rb");
  if(!infile) {
    ckerr << "TreePiece " << thisIndex << ": Couldn't open masses file, aborting" << endl;
    CkAbort("Badness");
  }
  xdrstdio_create(&xdrs, infile, XDR_DECODE);
	
  if(!xdr_template(&xdrs, &fh)) {
    ckerr << "TreePiece " << thisIndex << ": Couldn't read header from masses file, aborting" << endl;
    CkAbort("Badness");
  }
	
  if(fh.magic != FieldHeader::MagicNumber || fh.dimensions != 1 || fh.code != float32) {
    ckerr << "TreePiece " << thisIndex << ": Masses file is corrupt or of incorrect type, aborting" << endl;
    CkAbort("Badness");
  }

  unsigned int *startParticles;
  unsigned int *numParticlesChunk;
  unsigned int excess;
  switch (domainDecomposition) {
  case SFC_dec:
    numParticlesChunk = new unsigned int[2];
    startParticles = new unsigned int[1];
    numParticlesChunk[0] = fh.numParticles / numTreePieces;
    numParticlesChunk[1] = 0; // sentinel for end chunks
    
    excess = fh.numParticles % numTreePieces;
    startParticles[0] = numParticlesChunk[0] * thisIndex;
    if(thisIndex < (int) excess) {
      numParticlesChunk[0]++;
      startParticles[0] += thisIndex;
    } else {
      startParticles[0] += excess;
    }
    myNumParticles = numParticlesChunk[0];
    break;
  case Oct_dec:
    CkAbort("Oct domain decomposition not yet implemented");
    break;
  case ORB_dec:
    CkAbort("ORB domain decomposition not yet implemented");
    break;
  default:
    CkAbort("Invalid domain decomposition requested");
  }

  /* At this point myNumParticles contain the number of particles to be loaded
     into this processro, startParticles and numParticlesChunk are newly
     allocated array containing the first particle and the count of each
     contiguous chunk of particles that has to be loaded. */

  // allocate an array for myParticles
  myParticles = new GravityParticle[myNumParticles + 2];

  if(thisIndex == 0)
    ckerr << " (" << fh.numParticles << ")";

  if(verbosity > 3)
    ckerr << "TreePiece " << thisIndex << ": Of " << fh.numParticles << " particles, taking " << startParticles[0] << " through " << (startParticles[0] + numParticlesChunk[0] - 1) << endl;

  float mass;
  float maxMass;
  if(!xdr_template(&xdrs, &mass) || !xdr_template(&xdrs, &maxMass)) {
    ckerr << "TreePiece " << thisIndex << ": Problem reading beginning of the mass file, aborting" << endl;
    CkAbort("Badness");
  }
  
  if(mass == maxMass) { //all the same mass
    for(u_int64_t i = 0; i < myNumParticles; ++i){
      myParticles[i + 1].mass = mass;
#if COSMO_STATS > 1
      myParticles[i + 1].intcellmass = 0;
      myParticles[i + 1].intpartmass = 0;
      myParticles[i + 1].extcellmass = 0;
      myParticles[i + 1].extpartmass = 0;
#endif
    }
#if COSMO_STATS > 0
    piecemass = myNumParticles*mass;
    //ckerr << "In a tree piece....mass of tree piece particles: " << piecemass << ", single particle; " << mass;
#endif
  } else {

    unsigned int myPart = 0;
    for (int chunkNum = 0; numParticlesChunk[chunkNum] > 0; ++chunkNum) {
      if(!seekField(fh, &xdrs, startParticles[chunkNum])) {
	ckerr << "TreePiece " << thisIndex << ": Could not seek to my part of the mass file, aborting" << endl;
	CkAbort("Badness");
      }

      for(unsigned int i = 0; i < numParticlesChunk[chunkNum]; ++i) {
	if(!xdr_template(&xdrs, &mass)) {
	  ckerr << "TreePiece " << thisIndex << ": Problem reading my part of the mass file, aborting" << endl;
	  CkAbort("Badness");
	}
	myParticles[++myPart].mass = mass;
#if COSMO_STATS > 1
	myParticles[myPart].intcellmass = 0;
	myParticles[myPart].intpartmass = 0;
	myParticles[myPart].extcellmass = 0;
	myParticles[myPart].extpartmass = 0;
#endif
#if COSMO_STATS > 0
	piecemass += mass;
#endif
      }
    }
    CkAssert(myPart == myNumParticles);
  }
  
  xdr_destroy(&xdrs);
  fclose(infile);
  
  for(u_int64_t i = 0; i < myNumParticles; ++i)
    myParticles[i + 1].soft = 0.0;

  infile = fopen((basefilename + ".pos").c_str(), "rb");
  if(!infile) {
    ckerr << "TreePiece " << thisIndex << ": Couldn't open positions file, aborting" << endl;
    CkAbort("Badness");
  }
  xdrstdio_create(&xdrs, infile, XDR_DECODE);
  
  FieldHeader posHeader;
  if(!xdr_template(&xdrs, &posHeader)) {
    ckerr << "TreePiece " << thisIndex << ": Couldn't read header from positions file, aborting" << endl;
    CkAbort("Badness");
  }
  
  if(posHeader.magic != FieldHeader::MagicNumber || posHeader.dimensions != 3 || posHeader.code != float32) {
    ckerr << "TreePiece " << thisIndex << ": Positions file is corrupt or of incorrect type, aborting" << endl;
    CkAbort("Badness");
  }
  
  if(posHeader.time != fh.time || posHeader.numParticles != fh.numParticles) {
    ckerr << "TreePiece " << thisIndex << ": Positions file doesn't match masses file, aborting" << endl;
    CkAbort("Badness");
  }
  
  Vector3D<float> pos;
  Vector3D<float> maxPos;
  if(!xdr_template(&xdrs, &pos) || !xdr_template(&xdrs, &maxPos)) {
    ckerr << "TreePiece " << thisIndex << ": Problem reading beginning of the positions file, aborting" << endl;
    CkAbort("Badness");
  }
  
  boundingBox.lesser_corner = pos;
  boundingBox.greater_corner = maxPos;
  
  if(pos == maxPos) { //all the same position
    //XXX This would be bad!
    Key k = generateKey(pos, boundingBox);
    for(u_int64_t i = 0; i < myNumParticles; ++i) {
      myParticles[i + 1].position = pos;
      myParticles[i + 1].key = k;
    }
  } else {

    unsigned int myPart = 0;
    for (int chunkNum = 0; numParticlesChunk[chunkNum] > 0; ++chunkNum) {
      if(!seekField(posHeader, &xdrs, startParticles[chunkNum])) {
	ckerr << "TreePiece " << thisIndex << ": Could not seek to my part of the positions file, aborting" << endl;
	CkAbort("Badness");
      }

      Key previous = 0;  // hold the last key generated
      Key current;
      //read all my particles' positions and make keys
      for(int i = 0; i < numParticlesChunk[chunkNum]; ++i) {
	if(!xdr_template(&xdrs, &pos)) {
	  ckerr << "TreePiece " << thisIndex << ": Problem reading my part of the positions file, aborting" << endl;
	  CkAbort("Badness");
	}
	myParticles[++myPart].position = pos;
	current = generateKey(pos, boundingBox);
	myParticles[myPart].key = current;
	//CkPrintf("Adding key: %d = %16llx\n",myPart,current);
	if (current < previous) {
	  CkPrintf("TreePiece %d: Key not ordered! (%016llx)\n",thisIndex,current);
	}
	previous = current;
      }
    }
    CkAssert(myPart == myNumParticles);
  }
	
  xdr_destroy(&xdrs);
  fclose(infile);
	
  if(verbosity > 3)
    ckerr << "TreePiece " << thisIndex << ": Read in masses and positions" << endl;
	
  contribute(0, 0, CkReduction::concat, cb);
}

void TreePiece::buildTree(int bucketSize, const CkCallback& cb) {
  maxBucketSize = bucketSize;
  callback = cb;
  // decide which logic are we using to divide the particles: Oct or ORB
  switch (useTree) {
  case Binary_Oct:
  case Oct_Oct:
    Key bounds[2];
    sort(myParticles+1, myParticles+myNumParticles+1);
#ifdef COSMO_PRINT
    CkPrintf("[%d] Keys: %016llx %016llx\n",thisIndex,myParticles[1].key,myParticles[myNumParticles].key);
#endif
    bounds[0] = myParticles[1].key;
    bounds[1] = myParticles[myNumParticles].key;
    contribute(2 * sizeof(Key), bounds, CkReduction::concat, CkCallback(CkIndex_TreePiece::collectSplitters(0), thisArrayID));
    break;
  case Binary_ORB:
    CkAbort("ORB logic for tree-build not yet implemented");
    break;
  }
}

class KeyDouble {
  Key first;
  Key second;
public:
  inline bool operator<(const KeyDouble& k) const {
    return first < k.first;
  }
};

void TreePiece::collectSplitters(CkReductionMsg* m) {
  numSplitters = 2 * numTreePieces;
  splitters = new Key[numSplitters];
  Key* splits = static_cast<Key *>(m->getData());
  copy(splits, splits + numSplitters, splitters);
  KeyDouble* splitters2 = (KeyDouble *)splitters;
  //sort(splitters, splitters + numSplitters);
  sort(splitters2, splitters2 + numTreePieces);
  for (unsigned int i=1; i<numSplitters; ++i) {
    if (splitters[i] < splitters[i-1]) {
      //for (int j=0; j<numSplitters; ++j) CkPrintf("%d: Key %d = %016llx\n",thisIndex,j,splitters[j]);
      CkAbort("Keys not ordered");
    }
  }
  splitters[0] = firstPossibleKey;
  contribute(0, 0, CkReduction::concat, CkCallback(CkIndex_TreePiece::startOctTreeBuild(0), thisArrayID));
  delete m;
  if(verbosity > 3)
    ckerr << "TreePiece " << thisIndex << ": Collected splitters" << endl;
}

void TreePiece::startOctTreeBuild(CkReductionMsg* m) {
  delete m;
	
  if(thisIndex == 0)
    myParticles[0].key = firstPossibleKey;
  else
    myParticles[0].key = splitters[2 * thisIndex - 1];
	
  if(thisIndex == (int) numTreePieces - 1)
    myParticles[myNumParticles + 1].key = lastPossibleKey;
  else
    myParticles[myNumParticles + 1].key = splitters[2 * thisIndex + 2];
	
  //leftBoundary = myParticles;
  //rightBoundary = myParticles + myNumParticles + 1;

  //ckerr << "Piece " << (myParticles + 1)->key << " : " << (myParticles + myNumParticles)->key << " has leftBoundary: " << leftBoundary->key << " rightBoundary: " << rightBoundary->key << endl;

  // create the root of the global tree
  switch (useTree) {
  case Binary_Oct:
    root = new BinaryTreeNode(1, numTreePieces>1?Tree::Boundary:Tree::Internal, 0, myNumParticles+1, 0);
    break;
  case Oct_Oct:
    //root = new OctTreeNode(1, Tree::Boundary, 0, myNumParticles+1, 0);
    break;
  default:
    CkAbort("We should have never reached here!");
  }

  if (thisIndex == 0) root->firstParticle ++;
  if (thisIndex == (int)numTreePieces-1) root->lastParticle --;
  root->particleCount = myNumParticles;
  nodeLookupTable[(Tree::NodeKey)1] = root;

  //root->key = firstPossibleKey;
  root->boundingBox = boundingBox;
  //nodeLookup[root->lookupKey()] = root;
  numBuckets = 0;
  bucketList.clear();
	
  boundaryNodesPending = 0;
	
  if(verbosity > 3)
    ckerr << "TreePiece " << thisIndex << ": Starting tree build" << endl;

  // recursively build the tree
  buildOctTree(root, 0);

  // check all the pending requests in for RemoteMoments
  for (MomentRequestType::iterator iter = momentRequests.begin(); iter != momentRequests.end(); iter++) {
    GenericTreeNode *node = keyToNode(iter->first);
    CkAssert(node != NULL);
    if (node->getType() == Empty || node->moments.totalMass > 0) {
      CkVec<int> *l = iter->second;
      for (int i=0; i<l->length(); ++i) {
	streamingProxy[(*l)[i]].receiveRemoteMoments(iter->first, node->getType(), node->firstParticle, node->particleCount, node->moments);
	//CkPrintf("[%d] sending moments of %s to %d upon treebuild finished\n",thisIndex,keyBits(node->getKey(),63).c_str(),(*l)[i]);
      }
      delete l;
      momentRequests.erase(node->getKey());
    }
  }

  /*
    char fout[100];
    sprintf(fout,"tree.%d.%d",thisIndex,iterationNo);
    ofstream ofs(fout);
    printTree(root,ofs);
    ofs.close();
  */

  //if(boundaryNodesPending == 0)
  //  contribute(0, 0, CkReduction::concat, callback);
  
  if(verbosity > 3)
    ckerr << "TreePiece " << thisIndex << ": Number of buckets: " << numBuckets << endl;
  if(verbosity > 3)
    ckerr << "TreePiece " << thisIndex << ": Finished tree build, resolving boundary nodes" << endl;
}

/*
/// Find what chare this node's left child resides on, and create it
inline SFCTreeNode* TreePiece::lookupLeftChild(SFCTreeNode* node) {
  SFCTreeNode* child = node->createLeftChild();
  nodeLookup[child->lookupKey()] = child;
  child->setType(NonLocal);
  tempNode.key = node->key;
  tempNode.level = node->level + 1;
  if(!nodeOwnership(&tempNode, &child->remoteIndex, &child->numOwners)){
    ckerr << "This is surprising, but may get taken care of." << endl;
    node->leftChild=NULL;
    nodeLookup.erase(child->lookupKey());
    delete child;
    return NULL;
  }	
  return child;
}

inline SFCTreeNode* TreePiece::lookupRightChild(SFCTreeNode* node) {
  SFCTreeNode* child = node->createRightChild();
  nodeLookup[child->lookupKey()] = child;
  child->setType(NonLocal);
  tempNode.key = node->rightChildKey();
  tempNode.level = node->level + 1;
  if(!nodeOwnership(&tempNode, &child->remoteIndex, &child->numOwners)){
    ckerr << "This is surprising, but may get taken care of." << endl;
    node->rightChild=NULL;
    nodeLookup.erase(child->lookupKey());
    delete child;
    return NULL;
  }	
  return child;
}
*/

/// Determine who are all the owners of this node
/// @return true if the caller is part of the owners, false otherwise
inline bool TreePiece::nodeOwnership(const GenericTreeNode *const node, int &firstOwner, int &lastOwner) {
  Key firstKey = Key(node->getKey());
  Key lastKey = Key(node->getKey() + 1);
  const Key mask = Key(1) << 63;
  while (! (firstKey & mask)) {
    firstKey <<= 1;
    lastKey <<= 1;
  }
  firstKey &= ~mask;
  lastKey &= ~mask;
  lastKey -= 1;
  Key *locLeft = lower_bound(splitters, splitters + numSplitters, firstKey);
  Key *locRight = lower_bound(locLeft, splitters + numSplitters, lastKey);
  firstOwner = (locLeft - splitters) >> 1;
  lastOwner = (locRight - splitters - 1) >> 1;
  std::string str = keyBits(node->getKey(),63);
#if COSMO_PRINT > 1
  CkPrintf("[%d] NO: key=%s, first=%d, last=%d\n",thisIndex,str.c_str(),locLeft-splitters,locRight-splitters);
#endif
  return (thisIndex >= firstOwner && thisIndex <= lastOwner);
}

/*
/// Determine if the node is owned, by how many and whom, and designate a "head" owner
inline bool TreePiece::nodeOwnership(SFCTreeNode* node, unsigned int* designatedOwner, unsigned int* numOwners, unsigned int* firstOwner, unsigned int* lastOwner) {
  //find the first place in the array splitters the left boundary of the node can go
  Key* locLeft = upper_bound(splitters, splitters + numSplitters, node->leftBoundary());
  //find the last place the right boundary can go
  Key* locRight = lower_bound(locLeft, splitters + numSplitters, node->rightBoundary());
  if(locLeft == locRight) { //node fits between two splitters
    if((locLeft - splitters) % 2) { //it fits inside a TreePiece
      if(numOwners)
	*numOwners = 1;
      unsigned int owner = (locLeft - splitters) / 2;
      if(firstOwner)
	*firstOwner = owner;
      if(lastOwner)
	*lastOwner = owner;
      if(designatedOwner)
	*designatedOwner = owner;
      return true;
    } else { //it falls between TreePieces
      ckerr << "Wow, I didn't think this could happen.  Live and learn." << endl;
      return false;
    }
  } else {
    //the index of first co-owner of the node
    unsigned int first = (locLeft - splitters) / 2;
    //the index of the last co-owner of the node
    unsigned int last = (locRight - splitters - 1) / 2;
    if(numOwners)
      *numOwners = last - first + 1;
    if(firstOwner)
      *firstOwner = first;
    if(lastOwner)
      *lastOwner = last;
    if(designatedOwner) //the designated owner is the one in the middle
      *designatedOwner = (first + last) / 2;
    return true;
  }
}
*/

/** A recursive algorithm for building my tree.
    Examines successive bits in the particles' keys, looking for splits.
    Each bit is a level of nodes in the tree.  We keep going down until
    we can bucket the particles.
*/
void TreePiece::buildOctTree(GenericTreeNode * node, int level) {

  if (level == 63) {
    ckerr << thisIndex << ": TreePiece: This piece of tree has exhausted all the bits in the keys.  Super double-plus ungood!" << endl;
    ckerr << "Left particle: " << (node->firstParticle) << " Right particle: " << (node->lastParticle) << endl;
    ckerr << "Left key : " << keyBits((myParticles[node->firstParticle]).key, 63).c_str() << endl;
    ckerr << "Right key: " << keyBits((myParticles[node->lastParticle]).key, 63).c_str() << endl;
    return;
  }

  CkAssert(node->getType() == Boundary || node->getType() == Internal);
  
  node->makeOctChildren(myParticles, myNumParticles, level);

  GenericTreeNode *child;
  for (unsigned int i=0; i<node->numChildren(); ++i) {
    child = node->getChildren(i);
    CkAssert(child != NULL);
    nodeLookupTable[child->getKey()] = child;
    if (child->getType() == NonLocal) {
      // find a remote index for the node
      int first, last;
      bool isShared = nodeOwnership(child, first, last);
      CkAssert(!isShared);
      child->remoteIndex = first + (thisIndex & (last-first));
      // if we have a remote child, the node is a Boundary. Thus count that we
      // have to receive one more message for the NonLocal node
      node->remoteIndex --;
      // request the remote chare to fill this node with the Moments
      streamingProxy[child->remoteIndex].requestRemoteMoments(child->getKey(), thisIndex);
      //CkPrintf("[%d] asking for moments of %s to %d\n",thisIndex,keyBits(child->getKey(),63).c_str(),child->remoteIndex);
    } else if (child->getType() == Internal && child->lastParticle - child->firstParticle < maxBucketSize) {
      CkAssert(child->firstParticle != 0 && child->lastParticle != myNumParticles+1);
      child->remoteIndex = thisIndex;
      child->makeBucket(myParticles);
      bucketList.push_back(child);
      numBuckets++;
      if (node->getType() != Boundary) node->moments += child->moments;
    } else if (child->getType() == Empty) {
      child->remoteIndex = thisIndex;
    } else {
      if (child->getType() == Internal) child->remoteIndex = thisIndex;
      // else the index is already 0
      buildOctTree(child, level+1);
      // if we have a Boundary child, we will have to compute it's multipole
      // before we can compute the multipole of the current node (and we'll do
      // it in receiveRemoteMoments)
      if (child->getType() == Boundary) node->remoteIndex --;
      if (node->getType() != Boundary) node->moments += child->moments;
    }
  }

  /* Old version with Boundary node collection. New version collect NonLocal
  if (node->getType() == Boundary) {
    int first, last;
    boundaryNodesPending++;
    bool isShared = nodeOwnership(node, first, last);
    CkAssert(isShared);
    node->numOwners = last - first; // keep track of how many messages we should receive
    int designed = (last+first)>>1;
    //CkPrintf("%016llx [%d] boundary node to %d for %d-%d\n",node->getKey(),thisIndex,designed,first,last);
    if (designed != thisIndex)
      pieces[designed].acceptBoundaryNodeContribution(node->getKey(), node->particleCount, node->moments);
      } else */
  if (node->getType() == Internal) {
    calculateRadiusFarthestCorner(node->moments, node->boundingBox);
  }
}

/** A recursive algorithm for building my tree.
    Examines successive bits in the particles' keys, looking for splits.
    Each bit is a level of nodes in the tree.  We keep going down until
    we can bucket the particles.  The left and right boundaries of this
    piece of tree will point to other pieces on other chares in the array.
*/
/*
void TreePiece::buildTree(GenericTreeNode* node, GravityParticle* leftParticle, GravityParticle* rightParticle) {
  node->beginParticle = leftParticle - myParticles;
  node->endParticle = (rightParticle - myParticles) + 1;
  if(leftParticle == leftBoundary)
    node->beginParticle++;
  if(rightParticle == rightBoundary)
    node->endParticle--;
	
  //check if we should bucket these particles
  if(rightParticle - leftParticle < maxBucketSize) {
    //can't bucket until we've cut at the boundary
    if((leftParticle != leftBoundary) && (rightParticle != rightBoundary)) {
      node->setType(Bucket);
      node->numOwners = 1;
      for(GravityParticle* iter = leftParticle; iter != rightParticle + 1; ++iter)
	node->moments += *iter;
      // calculateRadiusFarthestCorner(node->moments, node->boundingBox);
      calculateRadiusFarthestParticle(node->moments, leftParticle, rightParticle + 1);
      bucketList.push_back(node);
      numBuckets++;
      return;
    }
  } else if(node->level == 63) {
    ckerr << thisIndex << ": TreePiece: This piece of tree has exhausted all the bits in the keys.  Super double-plus ungood!" << endl;
    ckerr << "Left particle: " << (leftParticle - myParticles) << " Right particle: " << (rightParticle - myParticles) << endl;
    ckerr << "Left key : " << keyBits(leftParticle->key, 63) << endl;
    ckerr << "Right key: " << keyBits(rightParticle->key, 63) << endl;
    return;
  }
	
  //this is the bit we are looking at
  Key currentBitMask = static_cast<Key>(1) << (62 - node->level);
  //we need to know the bit values at the left and right
  Key leftBit = leftParticle->key & currentBitMask;
  Key rightBit = rightParticle->key & currentBitMask;
  SFCTreeNode* child;
	
  if(leftBit < rightBit) { //a split at this level
    //find the split by looking for where the key with the bit switched on could go
    GravityParticle* splitParticle = lower_bound(leftParticle, rightParticle + 1, node->key | currentBitMask);
    if(splitParticle == leftBoundary + 1) {
      //we need to make the left child point to a remote chare
      if(thisIndex != 0) //the left-most chare can't point any further left
	lookupLeftChild(node);
      child = node->createRightChild();
      nodeLookup[child->lookupKey()] = child;
      buildTree(child, splitParticle, rightParticle);
    } else if(splitParticle == rightBoundary) {
      //we need to make the right child point to a remote chare
      child = node->createLeftChild();
      nodeLookup[child->lookupKey()] = child;
      buildTree(child, leftParticle, splitParticle - 1);
      if(thisIndex != (int) numTreePieces - 1) //the right-most chare can't point any further right
	lookupRightChild(node);
    } else {
      //neither child is remote, keep going with them
      child = node->createLeftChild();
      nodeLookup[child->lookupKey()] = child;
      buildTree(child, leftParticle, splitParticle - 1);
      child = node->createRightChild();
      nodeLookup[child->lookupKey()] = child;
      buildTree(child, splitParticle, rightParticle);
    }
  } else if(leftBit & rightBit) { //both ones, make a right child
    //should the left child be remote?
    if(leftParticle == leftBoundary && thisIndex != 0)
      lookupLeftChild(node);
    child = node->createRightChild();
    nodeLookup[child->lookupKey()] = child;
    buildTree(child, leftParticle, rightParticle);
  } else if(leftBit > rightBit) {
    ckerr << "Bits not right: " << leftBit << " vs " << rightBit << endl;
    ckerr << "Left particle: " << (leftParticle - myParticles) << " Right particle: " << (rightParticle - myParticles) << endl;
    ckerr << "Left key : " << keyBits(leftParticle->key, 63) << endl;
    ckerr << "Right key: " << keyBits(rightParticle->key, 63) << endl;
    return;
  } else { //both zeros, make a left child
    child = node->createLeftChild();
    nodeLookup[child->lookupKey()] = child;
    buildTree(child, leftParticle, rightParticle);
    //should the right child be remote?
    if(rightParticle == rightBoundary && thisIndex != (int) numTreePieces - 1)
      lookupRightChild(node);
  }
	
  //children have been formed, do bottom-up collection
  if(node->leftChild) {
    node->moments += dynamic_cast<SFCTreeNode *>(node->leftChild)->moments;
  }
  if(node->rightChild) {
    node->moments += dynamic_cast<SFCTreeNode *>(node->rightChild)->moments;
  }
	
  // figure out if this is a boundary node
  if((leftParticle == leftBoundary && thisIndex != 0) || (rightParticle == rightBoundary && thisIndex != (int) numTreePieces - 1)) {
    //send my information about this node to the designated owner
    unsigned int designatedOwner;
    nodeOwnership(node, &designatedOwner, &node->numOwners);
    boundaryNodesPending++;
    //in boundary nodes, remoteIndex contains the total number of particles in the node, from all the co-owners
    //endParticle - beginParticle tells you how many particles are actually on this processor
    node->remoteIndex = node->endParticle - node->beginParticle;
    if((int) designatedOwner != thisIndex) //don't send yourself your contribution
      pieces[designatedOwner].acceptBoundaryNodeContribution(node->lookupKey(), node->remoteIndex, node->moments);
    node->setType(Boundary);
  } else {
    node->numOwners = 1;
    node->setType(Internal);
    calculateRadiusFarthestCorner(node->moments, node->boundingBox);
  }
}
*/

void TreePiece::acceptBoundaryNodeContribution(const Tree::NodeKey key, const int numParticles, const MultipoleMoments &moments) {

  GenericTreeNode *node = keyToNode(key);

  if(node == NULL) {
    //	ckerr << "Key: " << keyBits(lookupKey, 63) << endl;
    pieces[thisIndex].acceptBoundaryNodeContribution(key,numParticles,moments); 
    //	ckerr << "Zeroth particle: " << myParticles << endl;
    //	ckerr << "leftBoundary: " << leftBoundary << endl;
    /*ckerr << "My Left bound : " << keyBits(myParticles[0].key, 63) << endl;
      ckerr << "My Left bound : " << keyBits(leftBoundary->key, 63) << endl;
      ckerr << "My Right bound: " << keyBits(rightBoundary->key, 63) << endl;*/
    return;
  }

  //CkPrintf("%016llx [%d] received, owners=%d\n",key,thisIndex,node->numOwners);
	
  //merge new information
  node->particleCount += numParticles;
  node->moments += moments;
  //decrement number of contributions and, if done, send final information to all co-owners
  if(--node->remoteIndex == 0) {
    calculateRadiusFarthestCorner(node->moments, node->boundingBox);
    int firstOwner, lastOwner;
    //recalculate number of owners, get co-owners
    nodeOwnership(node, firstOwner, lastOwner);
    //send moments, numParticles to all co-owners
    //CkPrintf("%016llx [%d] sending data to %d-%d\n",key,thisIndex,firstOwner,lastOwner);
    for(int i = firstOwner; i <= lastOwner; ++i)
      pieces[i].acceptBoundaryNode(key, node->particleCount, node->moments);
  }
}

void TreePiece::acceptBoundaryNode(const Tree::NodeKey key, const int numParticles, const MultipoleMoments& moments) {

  GenericTreeNode *node = keyToNode(key);

  if(node == NULL) {
    ckerr << "Well crap, how the hell did this happen, especially now? " << key << endl;
    return;
  }
	
  if(node->getType() != Boundary)
    ckerr << "How does this work? " << getColor(node).c_str() << endl;
  //merge final information
  node->particleCount = numParticles;
  node->moments = moments;	
  boundaryNodesPending--;
  if(boundaryNodesPending == 0) {
    CkAbort("Deprecated in favor of requestRemoteMoments");
    calculateRemoteMoments(root);
    contribute(0, 0, CkReduction::concat, callback);
  }
}

void TreePiece::calculateRemoteMoments(GenericTreeNode* node) {
  BinaryTreeNode *bnode = (BinaryTreeNode*)node;
  if(node->getType() == NonLocal) {
    GenericTreeNode* sibling = bnode->getSibling();
    GenericTreeNode* parent = node->parent;
    node->firstParticle = 0;
    if(sibling->getType() == Boundary)
      node->lastParticle = parent->particleCount - sibling->particleCount;
    else
      node->lastParticle = parent->particleCount - (sibling->lastParticle - sibling->firstParticle+1);
    if(node->lastParticle != 0) {
      node->moments = parent->moments - sibling->moments;
      calculateRadiusFarthestCorner(node->moments, node->boundingBox);
    } else {
      node->makeEmpty();
    }
  } else if(node->getType() == Boundary) {
    calculateRemoteMoments(bnode->children[0]);
    calculateRemoteMoments(bnode->children[1]);
  }
}

void TreePiece::requestRemoteMoments(const Tree::NodeKey key, int sender) {
  GenericTreeNode *node = keyToNode(key);
  if (node != NULL && (node->getType() == Empty || node->moments.totalMass > 0)) {
    streamingProxy[sender].receiveRemoteMoments(key, node->getType(), node->firstParticle, node->particleCount, node->moments);
    //CkPrintf("[%d] sending moments of %s to %d directly\n",thisIndex,keyBits(node->getKey(),63).c_str(),sender);
  } else {
    CkVec<int> *l = momentRequests[key];
    if (l == NULL) {
      l = new CkVec<int>();
      momentRequests[key] = l;
      //CkPrintf("[%d] Inserting new CkVec\n",thisIndex);
    }
    l->push_back(sender);
    //CkPrintf("[%d] queued request from %d for %s\n",thisIndex,sender,keyBits(key,63).c_str());
  }
}

void TreePiece::receiveRemoteMoments(const Tree::NodeKey key, Tree::NodeType type, int firstParticle, int numParticles, const MultipoleMoments& moments) {
  GenericTreeNode *node = keyToNode(key);
  CkAssert(node != NULL);
  //CkPrintf("[%d] received moments for %s\n",thisIndex,keyBits(key,63).c_str());
  // assign the incoming moments to the node
  if (type == Empty) node->makeEmpty();
  else {
    if (type == Bucket) {
      node->setType(NonLocalBucket);
      node->firstParticle = firstParticle;
      node->lastParticle = firstParticle + numParticles - 1;
    }
    node->particleCount = numParticles;
    node->moments = moments;
  }
  // look if we can compute the moments of some ancestors, and eventually send
  // them to a requester
  GenericTreeNode *parent = node->parent;
  while (parent != NULL && ++parent->remoteIndex == 0) {
    // compute the multipole for the parent
    //CkPrintf("[%d] computed multipole of %s\n",thisIndex,keyBits(parent->getKey(),63).c_str());
    parent->particleCount = 0;
    parent->remoteIndex = thisIndex; // reset the reference index to ourself
    GenericTreeNode *child;
    for (unsigned int i=0; i<parent->numChildren(); ++i) {
      child = parent->getChildren(i);
      parent->particleCount += child->particleCount;
      parent->moments += child->moments;
    }
    calculateRadiusFarthestCorner(parent->moments, parent->boundingBox);
    // check if someone has requested this node
    MomentRequestType::iterator iter;
    if ((iter = momentRequests.find(parent->getKey())) != momentRequests.end()) {
      CkVec<int> *l = iter->second;
      for (int i=0; i<l->length(); ++i) {
	streamingProxy[(*l)[i]].receiveRemoteMoments(parent->getKey(), parent->getType(), parent->firstParticle, parent->particleCount, parent->moments);
	//CkPrintf("[%d] sending moments of %s to %d\n",thisIndex,keyBits(parent->getKey(),63).c_str(),(*l)[i]);
      }
      delete l;
      momentRequests.erase(parent->getKey());
    }
    // go to the next ancestor
    node = parent;
    parent = node->parent;
  }
  if (parent == NULL) {
    // if we are here then we are at the root, and thus we have finished to get
    // all moments
    //CkPrintf("[%d] contributing after building the tree\n",thisIndex);
    contribute(0, 0, CkReduction::concat, callback);
  }// else CkPrintf("[%d] still missing one child of %s\n",thisIndex,keyBits(parent->getKey(),63).c_str());
}

void TreePiece::calculateGravityDirect(const CkCallback& cb) {
  callback = cb;
	
  GravityRequest req;
  req.requestingPieceIndex = thisIndex;
	
  myNumParticlesPending = myNumParticles * numTreePieces;
	
  for(u_int64_t i = 1; i <= myNumParticles; ++i) {
    req.identifier = i;
    req.soft = myParticles[i].soft;
    req.position = myParticles[i].position;
    myParticles[i].acceleration = 0;
    myParticles[i].potential = 0;
    pieces.fillRequestDirect(req);
  }
	
  started = true;
}

inline void
partForce(GravityParticle *part, GravityRequest& req)
{
  Vector3D<double> r;
  double rsq;
  double twoh, a, b;

  r = part->position - req.position;
  rsq = r.lengthSquared();
  twoh = part->soft + req.soft;
  if(rsq != 0) {
    SPLINE(rsq, twoh, a, b);
    req.acceleration += part->mass * r * b;
    req.potential -= part->mass * a;
  }
}

inline void
partBucketForce(GravityParticle *part, BucketGravityRequest& req)
{
  Vector3D<double> r;
  double rsq;
  double twoh, a, b;

	//SFCTreeNode* reqnode = bucketList[req.identifier];
	
  for(unsigned int j = 0; j < req.numParticlesInBucket; ++j) {
    r = part->position - req.positions[j];
    rsq = r.lengthSquared();
    twoh = part->soft + req.softs[j];
    if(rsq != 0) {
      SPLINE(rsq, twoh, a, b);
      req.accelerations[j] += part->mass * r * b;
			//myParticles[reqnode->beginParticle + j].acceleration += part->mass * r * b;
      req.potentials[j] -= part->mass * a;
    }
  }
}


inline void
nodeForce(GenericTreeNode *node, GravityRequest& req)
{
  Vector3D<double> r;
  double rsq;
  double twoh, a, b, c, d;
  MultipoleMoments m = node->moments;
    
  Vector3D<double> cm(m.cm);

  r = req.position - cm;
  rsq = r.lengthSquared();
  twoh = m.soft + req.soft;
  if(rsq != 0) {
    double dir = 1.0/sqrt(rsq);
    SPLINEQ(dir, rsq, twoh, a, b, c, d);
    double qirx = m.xx*r[0] + m.xy*r[1] + m.xz*r[2];
    double qiry = m.xy*r[0] + m.yy*r[1] + m.yz*r[2];
    double qirz = m.xz*r[0] + m.yz*r[1] + m.zz*r[2];
    double qir = 0.5*(qirx*r[0] + qiry*r[1] + qirz*r[2]);
    double tr = 0.5*(m.xx + m.yy + m.zz);
    double qir3 = b*m.totalMass + d*qir - c*tr;
    req.potential -= m.totalMass * a + c*qir - b*tr;
    req.acceleration[0] -= qir3*r[0] - c*qirx;
    req.acceleration[1] -= qir3*r[1] - c*qiry;
    req.acceleration[2] -= qir3*r[2] - c*qirz;
  }
}

/** @TOOD change the loop so that for one particle in req all the forces are
    compute (i.e outside loop over req, inside loop over nodes (the loop is in
    other functions) */
inline void
nodeBucketForce(GenericTreeNode *node, BucketGravityRequest& req)
{
  Vector3D<double> r;
  double rsq;
  double twoh, a, b, c, d;
  MultipoleMoments m = node->moments;
    
  Vector3D<double> cm(m.cm);

	//SFCTreeNode* reqnode = bucketList[req.identifier];

  for(unsigned int j = 0; j < req.numParticlesInBucket; ++j) {
    r = req.positions[j] - cm;
    rsq = r.lengthSquared();
    twoh = m.soft + req.softs[j];
    if(rsq != 0) {
      double dir = 1.0/sqrt(rsq);
      SPLINEQ(dir, rsq, twoh, a, b, c, d);
      double qirx = m.xx*r[0] + m.xy*r[1] + m.xz*r[2];
      double qiry = m.xy*r[0] + m.yy*r[1] + m.yz*r[2];
      double qirz = m.xz*r[0] + m.yz*r[1] + m.zz*r[2];
      double qir = 0.5*(qirx*r[0] + qiry*r[1] + qirz*r[2]);
      double tr = 0.5*(m.xx + m.yy + m.zz);
      double qir3 = b*m.totalMass + d*qir - c*tr;
      req.potentials[j] -= m.totalMass * a + c*qir - b*tr;
      req.accelerations[j][0] -= qir3*r[0] - c*qirx;
      req.accelerations[j][1] -= qir3*r[1] - c*qiry;
      req.accelerations[j][2] -= qir3*r[2] - c*qirz;
    
			/******************ADDED**********************/
			//SFCTreeNode* reqnode = bucketList[req.identifier];

			//for(unsigned int i = reqnode->beginParticle; i < reqnode->endParticle; ++i){
			//myParticles[reqnode->beginParticle + j].acceleration[0] -= qir3*r[0] - c*qirx;
			//myParticles[reqnode->beginParticle + j].acceleration[1] -= qir3*r[1] - c*qiry;
			//myParticles[reqnode->beginParticle + j].acceleration[2] -= qir3*r[2] - c*qirz;
			//}

		}
  }
}

inline
bool TreePiece::openCriterion(GenericTreeNode *node,
			      GravityRequest& req)
{
  // Note that some of this could be pre-calculated into an "opening radius"
  Sphere<double> s(node->moments.cm,
		   opening_geometry_factor * node->moments.radius / theta);
  return Space::contains(s, req.position);
}

inline
bool TreePiece::openCriterionBucket(GenericTreeNode *node,
				    BucketGravityRequest& req)
{
  // Note that some of this could be pre-calculated into an "opening radius"
  Sphere<double> s(node->moments.cm,
		   opening_geometry_factor * node->moments.radius / theta);
  return Space::intersect(req.boundingBox, s);
}

void TreePiece::fillRequestDirect(GravityRequest req) {
  for(u_int64_t i = 1; i <= myNumParticles; ++i) {
    partForce(&myParticles[i], req);
  }
  streamingProxy[req.requestingPieceIndex].receiveGravityDirect(req);
}

void TreePiece::receiveGravityDirect(const GravityRequest& req) {
  myParticles[req.identifier].acceleration += req.acceleration;
  myParticles[req.identifier].potential += req.potential;
  if(started && --myNumParticlesPending == 0) {
    started = false;
    contribute(0, 0, CkReduction::concat, callback);
  }
}

void TreePiece::startNextParticle() {
  if(nextParticle > myNumParticles)
    return;
  GravityRequest req;
  req.startingNode = root->getKey();
  req.requestingPieceIndex = thisIndex;
  req.identifier = nextParticle;
  req.soft = myParticles[nextParticle].soft;
  req.position = myParticles[nextParticle].position;
  myParticles[nextParticle].treeAcceleration = 0;
  streamingProxy[thisIndex].fillRequestTree(req);
  nextParticle++;
}

void TreePiece::calculateGravityTree(double t, const CkCallback& cb) {
  callback = cb;
  theta = t;
  mySerialNumber = 0;
#if COSMO_STATS > 0
  myNumCellInteractions = 0;
  myNumParticleInteractions = 0;
  myNumMACChecks = 0;
  myNumProxyCalls = 0;
  myNumProxyCallsBack = 0;
#endif

  nextParticle = 1;
  startNextParticle();
	
  myNumParticlesPending = myNumParticles;
  started = true;
	
}

void TreePiece::fillRequestTree(GravityRequest req) {
  //double startTime = CkWallTimer();
  //lookup starting node using startingNode key
  GenericTreeNode *node = keyToNode(req.startingNode);
  if(node == NULL) {
    ckerr << "Well crap, how the hell did this happen here?" << endl;
    return;
  }
	
  //make the request ready to go in the queue
  req.numAdditionalRequests = 1;
  req.acceleration = 0;
  req.potential = 0;
  req.numCellInteractions = 0;
  req.numParticleInteractions = 0;
  req.numMACChecks = 0;
  req.numEntryCalls = 0;
	
  //enter the requests into the queue
  unfilledRequests[mySerialNumber] = req;
	
  // make this request ready to get sent to other pieces
  req.requestingPieceIndex = thisIndex;
  req.identifier = mySerialNumber;
		
  walkTree(node, req);
	
  receiveGravityTree(req);
	
  mySerialNumber++;

  startNextParticle();
		
  //CkPrintf("%g\t%d\t%d\t%d\n", CkWallTimer() - startTime, req.numMACChecks, req.numCellInteractions, req.numParticleInteractions);
  //cout << (CkWallTimer() - startTime) << '\t' << req.numMACChecks << '\t' << req.numCellInteractions << '\t' << req.numParticleInteractions << '\n';
}

void TreePiece::walkTree(GenericTreeNode* node, GravityRequest& req) {
  req.numMACChecks++;
#if COSMO_STATS > 0
  myNumMACChecks++;
#endif
  if(!openCriterion(node, req)) {
    req.numCellInteractions++;
#if COSMO_STATS > 0
    myNumCellInteractions++;
#endif
    nodeForce(node, req);
  } else if(node->getType() == Bucket) {
    req.numParticleInteractions += node->lastParticle - node->firstParticle + 1;
#if COSMO_STATS > 0
    myNumParticleInteractions += node->lastParticle - node->firstParticle + 1;
#endif
    for(int i = node->firstParticle; i <= node->lastParticle; ++i) {
      partForce(&myParticles[i], req);
    }
  } else if(node->getType() == NonLocal || node->getType() == NonLocalBucket) {
    unfilledRequests[mySerialNumber].numAdditionalRequests++;
    req.numEntryCalls++;
    req.startingNode = node->getKey();
    streamingProxy[node->remoteIndex].fillRequestTree(req);
#if COSMO_STATS > 0
    myNumProxyCalls++;
#endif
  } else {
    GenericTreeNode* childIterator;
    for(unsigned int i = 0; i < node->numChildren(); ++i) {
      childIterator = node->getChildren(i);
      if(childIterator)
	walkTree(childIterator, req);
    }
  }
}

void TreePiece::receiveGravityTree(const GravityRequest& req) {
  //lookup request
  UnfilledRequestsType::iterator requestIter = unfilledRequests.find(req.identifier);
  if(requestIter == unfilledRequests.end()) {
    ckerr << "Well crap, how the hell did this happen here and now?" << endl;
    ckerr << "TreePiece " << thisIndex << ": Got request from " << req.requestingPieceIndex << " with id " << req.identifier << endl;
    return;
  }
  GravityRequest& request = requestIter->second;
  request.merge(req);
  if(--request.numAdditionalRequests == 0) {
    if((int) request.requestingPieceIndex == thisIndex) {
      myNumParticlesPending--;
      //this request originated here, it's for one of my particles
      myParticles[request.identifier].update(request);
    } else {
      streamingProxy[request.requestingPieceIndex].receiveGravityTree(request);
#if COSMO_STATS > 0
      myNumProxyCallsBack++;
#endif
    }
		
    unfilledRequests.erase(requestIter);
    if(started && myNumParticlesPending == 0) {
      started = false;
      contribute(0, 0, CkReduction::concat, callback);
#if COSMO_STATS > 0
      cout << "TreePiece " << thisIndex << ": Made " << myNumProxyCalls << " proxy calls forward, " << myNumProxyCallsBack << " to respond" << endl;
#endif
      if(verbosity > 4)
	ckerr << "TreePiece " << thisIndex << ": My particles are done" << endl;
    }
  }
}

#ifdef SEND_VERSION
void TreePiece::startNextBucket() {
  if(currentBucket >= numBuckets)
    return;
	
  SFCTreeNode* node = bucketList[currentBucket++];
  unsigned int numParticlesInBucket = node->endParticle - node->beginParticle;
  BucketGravityRequest req(numParticlesInBucket);
  req.startingNode = root->key;
  req.identifier = node->beginParticle;
  req.requestingPieceIndex = thisIndex;
  for(unsigned int i = node->beginParticle; i < node->endParticle; ++i) {
    req.softs[i - node->beginParticle] = myParticles[i].soft;
    req.positions[i - node->beginParticle] = myParticles[i].position;
    req.boundingBox.grow(myParticles[i].position);
    myParticles[i].treeAcceleration = 0;
  }
  streamingProxy[thisIndex].fillRequestBucketTree(req);
#if COSMO_STATS > 0
  myNumProxyCalls++;
#endif
}
#else
void TreePiece::startNextBucket() {
  if(currentBucket >= numBuckets)
    return;
	
  GenericTreeNode* node = bucketList[currentBucket];
  int numParticlesInBucket = node->particleCount;

  CkAssert(numParticlesInBucket <= maxBucketSize);
  BucketGravityRequest req(numParticlesInBucket);
  req.startingNode = root->getKey();
  req.identifier = currentBucket;
  req.requestingPieceIndex = thisIndex;
  for(int i = node->firstParticle; i <= node->lastParticle; ++i) {
    req.softs[i - node->firstParticle] = myParticles[i].soft;
    req.positions[i - node->firstParticle] = myParticles[i].position;
    req.potentials[i - node->firstParticle] = 0;
    req.boundingBox.grow(myParticles[i].position);
    myParticles[i].treeAcceleration = 0;
  }
  req.finished = 0;
  bucketReqs[currentBucket] = req;
  walkBucketTree(root, bucketReqs[currentBucket]);
  bucketReqs[currentBucket].finished = 1;
  finishBucket(currentBucket);
  //	currentBucket++;
  //	startNextBucket();
}

void TreePiece::finishBucket(int iBucket) {
  BucketGravityRequest *req = &bucketReqs[iBucket];
#ifdef COSMO_PRINT
  CkPrintf("[%d] Is finished %d? %d still missing!\n",thisIndex,iBucket,req->numAdditionalRequests);
#endif

  if(req->finished && req->numAdditionalRequests == 0) {
    myNumParticlesPending -= req->numParticlesInBucket;
#ifdef COSMO_PRINT
    CkPrintf("[%d] Finished bucket %d, %d particles remaining\n",thisIndex,iBucket,myNumParticlesPending);
#endif
    int iStart = bucketList[iBucket]->firstParticle;
    for(unsigned int i = 0; i < req->numParticlesInBucket; ++i) {
      myParticles[iStart + i].treeAcceleration
	+= req->accelerations[i];
      myParticles[iStart + i].potential
	+= req->potentials[i];
    }
    if(started && myNumParticlesPending == 0) {
      started = false;
      contribute(0, 0, CkReduction::concat, callback);
      /*   cout << "TreePiece " << thisIndex << ": Made " << myNumProxyCalls
	   << " proxy calls forward, " << myNumProxyCallsBack
	   << " to respond in finishBucket" << endl;*/
      if(verbosity)
	CkPrintf("[%d] TreePiece %d finished with bucket %d \n",CkMyPe(),thisIndex,iBucket);
      if(verbosity > 4)
	ckerr << "TreePiece " << thisIndex << ": My particles are done"
	     << endl;
    }
  }
}

#endif

void TreePiece::doAllBuckets(){
#if COSMO_DEBUG > 0
    char fout[100];
    sprintf(fout,"tree.%d.%d",thisIndex,iterationNo);
    ofstream ofs(fout);
    printTree(root,ofs);
    ofs.close();
#endif

  /*for(;currentBucket <numBuckets;currentBucket++){
    startNextBucket();
    if(currentBucket%YIELDPERIOD == YIELDPERIOD -1 ){
      CthYield();
      }	
  }*/
  dummyMsg *msg = new (32) dummyMsg;
  *((int *)CkPriorityPtr(msg))=10*(1+thisIndex);
  CkSetQueueing(msg,CK_QUEUEING_IFIFO);
  msg->val=0;
  thisProxy[thisIndex].nextBucket(msg);
}

void TreePiece::nextBucket(dummyMsg *msg){
  unsigned int i=0;
  while(i<_yieldPeriod && currentBucket<numBuckets){
    startNextBucket();
    currentBucket++;
    i++;
  }
  if(currentBucket<numBuckets){
    thisProxy[thisIndex].nextBucket(msg);
  }
}

void TreePiece::calculateGravityBucketTree(double t, const CkCallback& cb) {
  callback = cb;
  theta = t;
  mySerialNumber = 0;
#if COSMO_STATS > 0
  myNumProxyCalls = 0;
  myNumProxyCallsBack = 0;
  myNumCellInteractions=myNumParticleInteractions=myNumMACChecks=0;
  cachecellcount=0;
#endif
  iterationNo++;
  if(localCache == NULL){
    localCache = cacheManagerProxy.ckLocalBranch();
  }
  localCache->cacheSync(iterationNo, root);
  if(verbosity)
    CkPrintf("TreePiece %d: I have %d buckets\n",thisIndex,numBuckets);
	
  bucketReqs = new BucketGravityRequest[numBuckets];
	
  currentBucket = 0;
  myNumParticlesPending = myNumParticles;
  started = true;
#if COSMO_STATS > 0
  countIntersects=0;
  countHits=0;
#endif
  //	startNextBucket();
  doAllBuckets();
}

void TreePiece::fillRequestBucketTree(BucketGravityRequest req) {
  GenericTreeNode *node = keyToNode(req.startingNode);
  if(node == NULL) {
    ckerr << "Well crap, how the hell did this happen here?" << endl;
    return;
  }
	
  //make the request ready to go in the queue
  req.numAdditionalRequests = 1;
  for(unsigned int i = 0; i < req.numParticlesInBucket; ++i) {
    req.accelerations[i] = 0;
    req.potentials[i] = 0;
  }
	
  //enter the requests into the queue
  unfilledBucketRequests[mySerialNumber] = req;
	
  // make this request ready to get sent to other pieces
  req.requestingPieceIndex = thisIndex;
  req.identifier = mySerialNumber;
		
  walkBucketTree(node, req);
	
  receiveGravityBucketTree(req);
	
  startNextBucket();
	
  mySerialNumber++;
}

void TreePiece::startlb(CkCallback &cb){
  callback = cb;
	if(verbosity > 1)
	CkPrintf("[%d] TreePiece %d calling AtSync()\n",CkMyPe(),thisIndex);
  AtSync();
}

void TreePiece::ResumeFromSync(){
	if(verbosity > 1)
	CkPrintf("[%d] TreePiece %d in ResumefromSync\n",CkMyPe(),thisIndex);
  contribute(0, 0, CkReduction::concat, callback);
}

GenericTreeNode *TreePiece::keyToNode(const Tree::NodeKey k) {
  NodeLookupType::iterator iter = nodeLookupTable.find(k);
  if (iter != nodeLookupTable.end()) return iter->second;
  else return NULL;
}

const GenericTreeNode *TreePiece::lookupNode(Tree::NodeKey key){
  return keyToNode(key);
  /*
  SFCTreeNode* node = nodeLookup[lookupKey];
  if(node != NULL){
    copySFCTreeNode(*res,node);
  }else{
    res->setType(Empty);
  }
  */
};

const GravityParticle *TreePiece::lookupParticles(int begin) {
  return &myParticles[begin];
}

#ifdef SEND_VERSION
/*
 * "Send" version of Treewalk.
 * When off non-local data is need the tree walk is continued on the
 * appropriate remote node.
 */
void TreePiece::walkBucketTree(GravityTreeNode* node, BucketGravityRequest& req) {
#if COSMO_STATS > 0
  myNumMACChecks++;
#endif
  if(!openCriterionBucket(node, req)) {
#if COSMO_STATS > 0
    myNumCellInteractions += req.numParticlesInBucket;
#endif
    nodeBucketForce(node, req);
  } else if(node->getType() == Bucket) {
#if COSMO_STATS > 0
    myNumParticleInteractions += req.numParticlesInBucket * (node->endParticle - node->beginParticle);
#endif
    for(unsigned int i = node->beginParticle; i < node->endParticle; ++i) {
      partBucketForce(&myParticles[i], req);
    }
  } else if(node->getType() == NonLocal || node->getType() == NonLocalBucket) {
    unfilledBucketRequests[mySerialNumber].numAdditionalRequests++;
    req.startingNode = dynamic_cast<SFCTreeNode *>(node)->lookupKey();
    streamingProxy[node->remoteIndex].fillRequestBucketTree(req);
#if COSMO_STATS > 0
    myNumProxyCalls++;
#endif
  } else {
    GenericTreeNode** childrenIterator = node->getChildren();
    for(unsigned int i = 0; i < node->numChildren(); ++i) {
      if(childrenIterator[i])
	walkBucketTree(dynamic_cast<GravityTreeNode *>(childrenIterator[i]), req);
    }
  }
}
#else

/*
 * For cached version we have 2 walks: one for on processor and one
 * that hits the cache.
 * When remote data is needed we go to the second version.
 */
void TreePiece::walkBucketTree(GenericTreeNode* node, BucketGravityRequest& req) {
#if COSMO_STATS > 0
  myNumMACChecks++;
#endif
  Vector3D<double> cm(node->moments.cm);
  Vector3D<double> r;
  Sphere<double> s(cm, opening_geometry_factor * node->moments.radius / theta);
  if(!openCriterionBucket(node, req)) {
#if COSMO_STATS > 0
    countIntersects++;
    myNumCellInteractions += req.numParticlesInBucket;
#endif
#if COSMO_STATS > 1
    MultipoleMoments m = node->moments;	
    GenericTreeNode* reqnode = bucketList[req.identifier];
    for(int i = reqnode->firstParticle; i <= reqnode->lastParticle; ++i)
      myParticles[i].intcellmass += m.totalMass;
#endif
#if COSMO_PRINT > 1
  CkPrintf("[%d] walk bucket %s -> node %s\n",thisIndex,keyBits(reqnode->getKey(),63).c_str(),keyBits(node->getKey(),63).c_str());
#endif
    nodeBucketForce(node, req);
  } else if(node->getType() == Bucket) {
#if COSMO_STATS > 0
    myNumParticleInteractions += req.numParticlesInBucket * (node->lastParticle - node->firstParticle + 1);
#endif
    GenericTreeNode* reqnode = bucketList[req.identifier];
    for(int i = node->firstParticle; i <= node->lastParticle; ++i) {
#if COSMO_STATS > 1
      for(int j = reqnode->firstParticle; j <= reqnode->lastParticle; ++j) {
	myParticles[j].intpartmass += myParticles[i].mass;
      }
#endif
  GenericTreeNode *reqnode = bucketList[req.identifier];
#if COSMO_PRINT > 1
  CkPrintf("[%d] walk bucket %s -> part %016llx\n",thisIndex,keyBits(reqnode->getKey(),63).c_str(),myParticles[i].key);
#endif
      partBucketForce(&myParticles[i], req);
    }
  } else if(node->getType() == NonLocal || node->getType() == NonLocalBucket) {
    // Use cachedWalkBucketTree() as callback
    GenericTreeNode *pnode = requestNode(node->remoteIndex, node->getKey(), req);
    if(pnode) {
#if COSMO_STATS > 0
      countHits++;
#endif
      cachedWalkBucketTree(pnode, req);
    }
  } else if (node->getType() != Empty) {
    // here the node can be Internal or Boundary
    GenericTreeNode* childIterator;
    for(unsigned int i = 0; i < node->numChildren(); ++i) {
      childIterator = node->getChildren(i);
      if(childIterator)
	walkBucketTree(childIterator, req);
    }
  }
}

/*
 * Cached version of Tree walk:
 */
void TreePiece::cachedWalkBucketTree(GenericTreeNode* node, BucketGravityRequest& req) {
#if COSMO_STATS > 0
  myNumMACChecks++;
#endif
#if COSMO_PRINT > 1
  CkPrintf("[%d] cachedWalkBucketTree called with node %s of type %s\n",thisIndex,keyBits(node->getKey(),63).c_str(),getColor(node).c_str());
#endif
		
  CkAssert(node->getType() != Invalid);
	
  if(!openCriterionBucket(node, req)) {
#if COSMO_STATS > 0
    myNumCellInteractions += req.numParticlesInBucket;
    cachecellcount+=req.numParticlesInBucket;
#endif
#if COSMO_STATS > 1
    MultipoleMoments m = node->moments;
    GenericTreeNode* reqnode = bucketList[req.identifier];
    for(int i = reqnode->firstParticle; i <= reqnode->lastParticle; ++i)
      myParticles[i].extcellmass += m.totalMass;
#endif
#if COSMO_PRINT > 1
  CkPrintf("[%d] cachedwalk bucket %s -> node %s\n",thisIndex,keyBits(reqnode->getKey(),63).c_str(),keyBits(node->getKey(),63).c_str());
#endif
    nodeBucketForce(node, req);
  } else if(node->getType() == CachedBucket || node->getType() == Bucket || node->getType() == NonLocalBucket) {
    /*
     * Sending the request for all the particles at one go, instead of one by one
     */
    GravityParticle *part = requestParticles(node->getKey(),node->remoteIndex,node->firstParticle,node->lastParticle,req);
    if(part != NULL){
#if COSMO_STATS > 0
      myNumParticleInteractions += req.numParticlesInBucket * (node->lastParticle - node->firstParticle + 1);
#endif
      GenericTreeNode* reqnode = bucketList[req.identifier];

      for(int i = node->firstParticle; i <= node->lastParticle; ++i) {
#if COSMO_STATS > 1
	for(int j = reqnode->firstParticle; j <= reqnode->lastParticle; ++j) {
	  myParticles[j].extpartmass += myParticles[i].mass;
	}
#endif
#if COSMO_PRINT > 1
	CkPrintf("[%d] cachedwalk bucket %s -> part %016llx\n",thisIndex,keyBits(reqnode->getKey(),63).c_str(),part[i-node->firstParticle].key);
#endif
	partBucketForce(&part[i-node->firstParticle], req);
      }
    }	
    /*
  } else if(node->getType() == NonLocal) {
    // Use cachedWalkBucketTree() as callback
    GenericTreeNode *pnode = requestNode(node->remoteIndex, node->getKey(), req);
    if(pnode) {
#if COSMO_PRINT > 1
      CkPrintf("[%d] Requested node %s to %d, got %s of type %s\n",thisIndex,keyBits(node->getKey(),63).c_str(),node->remoteIndex,keyBits(pnode->getKey(),63).c_str(),getColor(pnode).c_str());
#endif
      cachedWalkBucketTree(pnode, req);
    }
    */
  } else if (node->getType() != CachedEmpty && node->getType() != Empty) {
    // Here the type is Cached, Boundary, Internal, NonLocal, which means the
    // node in the global tree has children (it is not a leaf), so we iterate
    // over them. If we get a NULL node, then we missed the cache and we request
    // it

    // Warning, since the cache returns nodes with pointers to other chare
    // elements trees, we could be actually traversing the tree of another chare
    // in this processor.

    // Use cachedWalkBucketTree() as callback
    GenericTreeNode *child;
    for (unsigned int i=0; i<node->numChildren(); ++i) {
      child = node->getChildren(i); //requestNode(node->remoteIndex, node->getChildKey(i), req);
      if (child) {
	cachedWalkBucketTree(child, req);
      } else { //missed the cache
	child = requestNode(node->remoteIndex, node->getChildKey(i), req);
	if (child) {
	  cachedWalkBucketTree(child, req);
	}
      }
    }
  }
}

GenericTreeNode* TreePiece::requestNode(int remoteIndex, Tree::NodeKey key,
				    BucketGravityRequest& req)
{
  // Call proxy on remote node
  CkAssert(remoteIndex < (int) numTreePieces);
  //in the current form it is possible   
  //   assert(remoteIndex != thisIndex);
  if(_cache){	
    if(localCache == NULL){
      localCache = cacheManagerProxy.ckLocalBranch();
    }
#if COSMO_PRINT > 1
    CkPrintf("[%d] requesting node %s to %d for %s\n",thisIndex,keyBits(key,63).c_str(),remoteIndex,keyBits(bucketList[req.identifier]->getKey(),63).c_str());
#endif
    GenericTreeNode *res=localCache->requestNode(thisIndex,remoteIndex,key,&req);
    if(!res){
      req.numAdditionalRequests++;
#if COSMO_STATS > 0
      myNumProxyCalls++;
#endif
    }
    return res;
  }
  else{	
    CkAbort("Non cached version not anymore supported, feel free to fix it!");
    return NULL;
    /*
    req.numAdditionalRequests++;
    streamingProxy[remoteIndex].fillRequestNode(thisIndex, key, req.identifier);
    myNumProxyCalls++;
    return NULL;
    */
  }
}


/*
  copy the data from node to tmp
*/
/*
void TreePiece::copySFCTreeNode(SFCTreeNode &tmp,SFCTreeNode *node){
  if(node == NULL){
    tmp.setType(Empty);
    return;
  }
  tmp.setType(node->getType());
  tmp.moments = node->moments;
  tmp.beginParticle = node->beginParticle;
  tmp.endParticle = node->endParticle;
  tmp.remoteIndex = node->remoteIndex;
  tmp.key = node->key;
  tmp.level = node->level;

  assert(tmp.getType() != Invalid);
	
  if(tmp.getType() == Boundary || tmp.getType() == Internal
     || tmp.getType() == Bucket)
    tmp.remoteIndex = thisIndex;
	
}
*/
/*
  do a prefix traversal starting at node and copy the keys and nodes into 
  the passed arrays
*/
/*
void TreePiece::prefixCopyNode(SFCTreeNode *node,Key lookupKey,Key *cacheKeys,SFCTreeNode *cacheNodes,int *count,int depth){
  if(depth >= _cacheLineDepth){
    return;
  }
  copySFCTreeNode(cacheNodes[*count],node);
  if(node != NULL){
    assert(lookupKey == node->lookupKey());
    cacheKeys[*count] = node->lookupKey();
  }else{
    cacheKeys[*count] = lookupKey;
  }
  (*count)++;
  if(node == NULL){
    return;
  }
  prefixCopyNode(nodeLookup[node->leftChildLookupKey()],node->leftChildLookupKey(),cacheKeys,cacheNodes,count,depth+1);
  prefixCopyNode(nodeLookup[node->rightChildLookupKey()],node->rightChildLookupKey(),cacheKeys,cacheNodes,count,depth+1);
}
*/
//void TreePiece::fillRequestNode(int retIndex, Tree::NodeKey key,
//				BucketGravityRequest& req)
//{

void TreePiece::fillRequestNode(RequestNodeMsg *msg) {
  GenericTreeNode* node = keyToNode(msg->key);
  //GenericTreeNode tmp;
  if(node != NULL) {
    if(_cache) {
      PUP::sizer p1;
      node->pup(p1, msg->depth);
      //CkPrintf("Requested size: %d\n",p1.size());
      FillNodeMsg *reply = new (p1.size(), 0) FillNodeMsg(thisIndex);
      //Key *cacheKeys = new Key[number];
      //SFCTreeNode *cacheNodes = new SFCTreeNode[number];
      //int count=0;
      //prefixCopyNode(node,lookupKey,cacheKeys,cacheNodes,&count,0);
      //	cacheManagerProxy[retIndex].recvNodes(lookupKey,thisIndex,tmp);
      //cacheManagerProxy[msg->retIndex].recvNodes(count,cacheKeys,cacheNodes,thisIndex);

      /// @TODO: check that at destination of "remoteIndex" are correct
      PUP::toMem p2((void*)reply->nodes);
      node->pup(p2, msg->depth);
      //int count = node->copyTo(reply->nodes, msg->depth);
      cacheManagerProxy[msg->retIndex].recvNodes(reply);
		
      //delete [] cacheKeys;
      //delete [] cacheNodes;
    }else{
      CkAbort("Non cached version not anymore supported, feel free to fix it!");
      //copySFCTreeNode(tmp,node);
      //streamingProxy[retIndex].receiveNode(tmp,msg->reqID);
    }
  }
  else {	// Handle NULL nodes
    CkAbort("Ok, before it handled this, but why do we have a null pointer in the tree?!?");
    /*
    tmp.setType(Empty);
    if(_cache){
      cacheManagerProxy[msg->retIndex].recvNodes(msg->key,thisIndex,tmp);
    }else{
      CkAbort("Non cached version not anymore supported, feel free to fix it!");
      //streamingProxy[msg->retIndex].receiveNode(tmp, msg->reqID);
    }
    */
  }
  delete msg;
}

void TreePiece::receiveNode(GenericTreeNode &node, unsigned int reqID)
{
  bucketReqs[reqID].numAdditionalRequests--;
  assert(node.getType() != Invalid);
  if(node.getType() != Empty)	{ // Node could be NULL
    assert((int) node.remoteIndex != thisIndex);
    cachedWalkBucketTree(&node, bucketReqs[reqID]);
  }else{
  }
    
  finishBucket(reqID);
}

void TreePiece::receiveNode_inline(GenericTreeNode &node, unsigned int reqID){
        receiveNode(node,reqID);
}
/*
  This function is not used anymore. It is extremely inefficient
  to request each particle in a node individually instead of requesting all
  the particles in a node as done in requestParticles defined below.
*/
/*
GravityParticle* TreePiece::requestParticle(int remoteIndex, int iPart,
					    BucketGravityRequest& req)
{
  assert(remoteIndex < (int) numTreePieces);
  // Call proxy on remote node
  req.numAdditionalRequests++;
  myNumProxyCalls++;

  streamingProxy[remoteIndex].fillRequestParticle(thisIndex, iPart, req);
    
  return NULL; // If we actually had a cache, this might return something
}
*/

GravityParticle *TreePiece::requestParticles(const Tree::NodeKey &key,int remoteIndex,int begin,int end,BucketGravityRequest& req){
  if (_cache) {
    if(localCache == NULL){
      localCache = cacheManagerProxy.ckLocalBranch();
    }
    GravityParticle *p = localCache->requestParticles(thisIndex,key,remoteIndex,begin,end,&req);
    if (!p) {
      req.numAdditionalRequests += end-begin+1;
    }
    return p;
  } else {
    CkAbort("Non cached version not anymore supported, feel free to fix it!");
    return NULL;
    /*
    req.numAdditionalRequests += end-begin;
    myNumProxyCalls++;
	
    streamingProxy[remoteIndex].fillRequestParticles(key,thisIndex,begin,end,req);
    return NULL;
    */
  }
};

/*
void TreePiece::fillRequestParticle(int retIndex, int iPart,
				    BucketGravityRequest& req)
{
  assert(retIndex < (int) numTreePieces);
  streamingProxy[retIndex].receiveParticle(myParticles[iPart], req);
}
*/

void TreePiece::fillRequestParticles(Key key,int retIndex, int begin,int end,
				     unsigned int reqID)
{
  if(_cache){
    cacheManagerProxy[retIndex].recvParticles(key,&myParticles[begin],end-begin+1,thisIndex);
  }else{
    streamingProxy[retIndex].receiveParticles(&myParticles[begin], end-begin+1,reqID);
  }	
}

/*
void TreePiece::receiveParticle(GravityParticle part,
				BucketGravityRequest& req)
{
  bucketReqs[req.identifier].numAdditionalRequests--;
  myNumParticleInteractions += bucketReqs[req.identifier].numParticlesInBucket;
  partBucketForce(&part, bucketReqs[req.identifier]);
  finishBucket(req.identifier);
}
*/

void TreePiece::receiveParticles(GravityParticle *part,int num,
				 unsigned int reqID)
{
  CkAssert(num > 0);
  bucketReqs[reqID].numAdditionalRequests -= num;
#if COSMO_STATS > 0
  myNumParticleInteractions += bucketReqs[reqID].numParticlesInBucket * num;
#endif

  GenericTreeNode* reqnode = bucketList[reqID];

  for(int i=0;i<num;i++){
#if COSMO_STATS > 1
    for(int j = reqnode->firstParticle; j <= reqnode->lastParticle; ++j) {
      myParticles[j].extpartmass += part[i].mass;
    }
#endif
#if COSMO_PRINT > 1
    CkPrintf("[%d] recvPart bucket %s -> part %016llx\n",thisIndex,keyBits(reqnode->getKey(),63).c_str(),part->key);
#endif
    partBucketForce(&part[i], bucketReqs[reqID]);
  }		
  finishBucket(reqID);
}

void TreePiece::receiveParticles_inline(GravityParticle *part,int num,
					unsigned int reqID){
        receiveParticles(part,num,reqID);
}


#endif

void TreePiece::receiveGravityBucketTree(const BucketGravityRequest& req) {
  //lookup request
  UnfilledBucketRequestsType::iterator requestIter = unfilledBucketRequests.find(req.identifier);
  if(requestIter == unfilledBucketRequests.end()) {
    ckerr << "Well crap, how the hell did this happen here and now?" << endl;
    ckerr << "TreePiece " << thisIndex << ": Got request from " << req.requestingPieceIndex << " with id " << req.identifier << endl;
    return;
  }
  BucketGravityRequest& request = requestIter->second;
  if(request.numParticlesInBucket != req.numParticlesInBucket)
    ckerr << "How could this be?" << endl;
  request.merge(req);
  if(--request.numAdditionalRequests == 0) {
    if((int) request.requestingPieceIndex == thisIndex) {
      myNumParticlesPending -= request.numParticlesInBucket;
      //this request originated here, it's for one of my particles
      for(unsigned int i = 0; i < request.numParticlesInBucket; ++i) {
	myParticles[request.identifier + i].treeAcceleration += request.accelerations[i];
			    
	myParticles[request.identifier + i].potential += request.potentials[i];
      }
    } else {
      streamingProxy[request.requestingPieceIndex].receiveGravityBucketTree(request);
#if COSMO_STATS > 0
      myNumProxyCallsBack++;
#endif
    }
		
    unfilledBucketRequests.erase(requestIter);
    if(started && myNumParticlesPending == 0) {
      started = false;
      contribute(0, 0, CkReduction::concat, callback);
#if COSMO_STATS > 0
      cout << "TreePiece " << thisIndex << ": Made " << myNumProxyCalls << " proxy calls forward, " << myNumProxyCallsBack << " to respond in receiveGravityBucketTree" << endl;
#endif
      if(verbosity > 4)
	ckerr << "TreePiece " << thisIndex << ": My particles are done" << endl;
    }
  }
}

void TreePiece::outputAccelerations(OrientedBox<double> accelerationBox, const string& suffix, const CkCallback& cb) {
  if(thisIndex == 0) {
    if(verbosity > 2)
      ckerr << "TreePiece " << thisIndex << ": Writing header for accelerations file" << endl;
    FILE* outfile = fopen((basefilename + "." + suffix).c_str(), "wb");
    XDR xdrs;
    xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
    fh.code = float64;
    fh.dimensions = 3;
    if(!xdr_template(&xdrs, &fh) || !xdr_template(&xdrs, &accelerationBox.lesser_corner) || !xdr_template(&xdrs, &accelerationBox.greater_corner)) {
      ckerr << "TreePiece " << thisIndex << ": Could not write header to accelerations file, aborting" << endl;
      CkAbort("Badness");
    }
    xdr_destroy(&xdrs);
    fclose(outfile);
  }
	
  if(verbosity > 3)
    ckerr << "TreePiece " << thisIndex << ": Writing my accelerations to disk" << endl;
	
  FILE* outfile = fopen((basefilename + "." + suffix).c_str(), "r+b");
  fseek(outfile, 0, SEEK_END);
  XDR xdrs;
  xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
	
  for(unsigned int i = 1; i <= myNumParticles; ++i) {
    accelerationBox.grow(myParticles[i].acceleration);
    if(!xdr_template(&xdrs, &(myParticles[i].acceleration))) {
      ckerr << "TreePiece " << thisIndex << ": Error writing accelerations to disk, aborting" << endl;
      CkAbort("Badness");
    }
  }
	
  if(thisIndex == (int) numTreePieces - 1) {
    if(!xdr_setpos(&xdrs, FieldHeader::sizeBytes) || !xdr_template(&xdrs, &accelerationBox.lesser_corner) || !xdr_template(&xdrs, &accelerationBox.greater_corner)) {
      ckerr << "TreePiece " << thisIndex << ": Error going back to write the acceleration bounds, aborting" << endl;
      CkAbort("Badness");
    }
    if(verbosity > 2)
      ckerr << "TreePiece " << thisIndex << ": Wrote the acceleration bounds" << endl;
    cb.send();
  }
	
  xdr_destroy(&xdrs);
  fclose(outfile);
	
  if(thisIndex != (int) numTreePieces - 1)
    pieces[thisIndex + 1].outputAccelerations(accelerationBox, suffix, cb);
}

/****************************ADDED***********************************/

void TreePiece::outputAccASCII(OrientedBox<double> accelerationBox, const string& suffix, const CkCallback& cb) {
  if((thisIndex==0 && packed) || (thisIndex==0 && !packed && cnt==0)) {
    if(verbosity > 2)
      ckerr << "TreePiece " << thisIndex << ": Writing header for accelerations file" << endl;
    FILE* outfile = fopen((basefilename + "." + suffix).c_str(), "wb");
		fprintf(outfile,"%d\n",fh.numParticles);
    fclose(outfile);
  }
	
	/*if(thisIndex==0) {
    if(verbosity > 2)
      ckerr << "TreePiece " << thisIndex << ": Writing header for accelerations file" << endl;
    FILE* outfile = fopen((basefilename + "." + suffix).c_str(), "wb");
		fprintf(outfile,"%d\n",fh.numParticles);
    fclose(outfile);
  }*/

  if(verbosity > 3)
    ckerr << "TreePiece " << thisIndex << ": Writing my accelerations to disk" << endl;
	
  FILE* outfile = fopen((basefilename + "." + suffix).c_str(), "r+b");
  fseek(outfile, 0, SEEK_END);
	
  for(unsigned int i = 1; i <= myNumParticles; ++i) {
    accelerationBox.grow(myParticles[i].treeAcceleration);
		Vector3D<double> acc = (myParticles[i].treeAcceleration);
		double val=0.0;
		if(!packed){
			if(cnt==0)
				val=acc.x;
			if(cnt==1)
				val=acc.y;
			if(cnt==2)
				val=acc.z;
		}
		switch(packed){
		case 1:
    if(fprintf(outfile,"%.14g\n",acc.x) < 0) {
      ckerr << "TreePiece " << thisIndex << ": Error writing accelerations to disk, aborting" << endl;
      CkAbort("Badness");
    }
  	if(fprintf(outfile,"%.14g\n",acc.y) < 0) {
      ckerr << "TreePiece " << thisIndex << ": Error writing accelerations to disk, aborting" << endl;
      CkAbort("Badness");
    }
		if(fprintf(outfile,"%.14g\n",acc.z) < 0) {
      ckerr << "TreePiece " << thisIndex << ": Error writing accelerations to disk, aborting" << endl;
      CkAbort("Badness");
    }
		break;
		case 0:
		if(fprintf(outfile,"%.14g\n",val) < 0) {
      ckerr << "TreePiece " << thisIndex << ": Error writing accelerations to disk, aborting" << endl;
      CkAbort("Badness");
    }
		break;
		}
	}
	cnt++;
  /*if(thisIndex==(int)numTreePieces-1) {
    cb.send();
  }*/

	if((thisIndex==(int)numTreePieces-1 && packed) || (thisIndex==(int)numTreePieces-1 && !packed && cnt==3)) {
    cb.send();
  }
  fclose(outfile);
	
	if(thisIndex==(int)numTreePieces-1 && !packed && cnt<3)
		pieces[0].outputAccASCII(accelerationBox, suffix, cb);
		
  if(thisIndex!=(int)numTreePieces-1)
    pieces[thisIndex + 1].outputAccASCII(accelerationBox, suffix, cb);
	
}
/********************************************************************/

void TreePiece::outputStatistics(Interval<unsigned int> macInterval, Interval<unsigned int> cellInterval, Interval<unsigned int> particleInterval, Interval<unsigned int> callsInterval, double totalmass, const CkCallback& cb) {

#if COSMO_STATS > 0
  if(verbosity > 1) {
    ckerr << "TreePiece ";
    ckerr << thisIndex;
    ckerr << ": Statistics\nMy number of MAC checks: ";
    ckerr << myNumMACChecks << endl;
    ckerr << "My number of particle-cell interactions: "
	 << myNumCellInteractions << " Per particle: "
	 << myNumCellInteractions/(double) myNumParticles
	 << "\nCache cell interactions count: " << cachecellcount << endl;
    ckerr << "My number of particle-particle interactions: "
	 << myNumParticleInteractions << " Per Particle: "
	 << myNumParticleInteractions/(double) myNumParticles
	 << endl;
  }
#endif	

#if COSMO_STATS > 1
  /*
	double calmass,prevmass;

	for(int i=1;i<=myNumParticles;i++){
		calmass = (myParticles[i].intcellmass + myParticles[i].intpartmass + myParticles[i].extcellmass + myParticles[i].extpartmass);
		if(i>1)
			prevmass = (myParticles[i-1].intcellmass + myParticles[i-1].intpartmass + myParticles[i-1].extcellmass + myParticles[i-1].extpartmass);
		//CkPrintf("treepiece:%d ,mass:%lf, totalmass:%lf\n",thisIndex,calmass,totalmass);
		if(i>1)
			if(calmass != prevmass)
				CkPrintf("Tree piece:%d -- particles %d and %d differ in calculated total mass\n",thisIndex,i-1,i);
		if(calmass != totalmass)
				CkPrintf("Tree piece:%d -- particle %d differs from total mass\n",thisIndex,i);
	}

	CkPrintf("TreePiece:%d everything seems ok..\n",thisIndex);
  */

  if(thisIndex == 0) {
    macInterval.max = 0;
    macInterval.min = macInterval.max - 1;
    cellInterval = macInterval;
    particleInterval = macInterval;
    callsInterval = macInterval;
		
    if(verbosity > 2)
      ckerr << "TreePiece " << thisIndex << ": Writing headers for statistics files" << endl;
    fh.dimensions = 1;
    fh.code = TypeHandling::uint32;
    FILE* outfile = fopen((basefilename + ".MACs").c_str(), "wb");
    XDR xdrs;
    xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
		
    unsigned int dummy;
    if(!xdr_template(&xdrs, &fh) || !xdr_template(&xdrs, &dummy) || !xdr_template(&xdrs, &dummy)) {
      ckerr << "TreePiece " << thisIndex << ": Could not write header to MAC file, aborting" << endl;
      CkAbort("Badness");
    }
    xdr_destroy(&xdrs);
    fclose(outfile);
		
    outfile = fopen((basefilename + ".cellints").c_str(), "wb");
    xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
    if(!xdr_template(&xdrs, &fh) || !xdr_template(&xdrs, &dummy) || !xdr_template(&xdrs, &dummy)) {
      ckerr << "TreePiece " << thisIndex << ": Could not write header to cell-interactions file, aborting" << endl;
      CkAbort("Badness");
    }
    xdr_destroy(&xdrs);
    fclose(outfile);

    outfile = fopen((basefilename + ".partints").c_str(), "wb");
    xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
    if(!xdr_template(&xdrs, &fh) || !xdr_template(&xdrs, &dummy) || !xdr_template(&xdrs, &dummy)) {
      ckerr << "TreePiece " << thisIndex << ": Could not write header to particle-interactions file, aborting" << endl;
      CkAbort("Badness");
    }
    xdr_destroy(&xdrs);
    fclose(outfile);

    outfile = fopen((basefilename + ".calls").c_str(), "wb");
    xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
    if(!xdr_template(&xdrs, &fh) || !xdr_template(&xdrs, &dummy) || !xdr_template(&xdrs, &dummy)) {
      ckerr << "TreePiece " << thisIndex << ": Could not write header to entry-point calls file, aborting" << endl;
      CkAbort("Badness");
    }
    xdr_destroy(&xdrs);
    fclose(outfile);
  }
	
  if(verbosity > 3)
    ckerr << "TreePiece " << thisIndex << ": Writing my statistics to disk" << endl;
	
  FILE* outfile = fopen((basefilename + ".MACs").c_str(), "r+b");
  fseek(outfile, 0, SEEK_END);
  XDR xdrs;
  xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
	
  for(unsigned int i = 1; i <= myNumParticles; ++i) {
    macInterval.grow(myParticles[i].numMACChecks);
    if(!xdr_template(&xdrs, &(myParticles[i].numMACChecks))) {
      ckerr << "TreePiece " << thisIndex << ": Error writing MAC checks to disk, aborting" << endl;
      CkAbort("Badness");
    }
  }
	
  if(thisIndex == (int) numTreePieces - 1) {
    if(verbosity > 3)
      ckerr << "MAC interval: " << macInterval << endl;
    if(!xdr_setpos(&xdrs, FieldHeader::sizeBytes) || !xdr_template(&xdrs, &macInterval.min) || !xdr_template(&xdrs, &macInterval.max)) {
      ckerr << "TreePiece " << thisIndex << ": Error going back to write the MAC bounds, aborting" << endl;
      CkAbort("Badness");
    }
    if(verbosity > 2)
      ckerr << "TreePiece " << thisIndex << ": Wrote the MAC bounds" << endl;
  }
	
  xdr_destroy(&xdrs);
  fclose(outfile);
	
  outfile = fopen((basefilename + ".cellints").c_str(), "r+b");
  fseek(outfile, 0, SEEK_END);
  xdrstdio_create(&xdrs, outfile, XDR_ENCODE);	
  for(unsigned int i = 1; i <= myNumParticles; ++i) {
    cellInterval.grow(myParticles[i].numCellInteractions);
    if(!xdr_template(&xdrs, &(myParticles[i].numCellInteractions))) {
      ckerr << "TreePiece " << thisIndex << ": Error writing cell interactions to disk, aborting" << endl;
      CkAbort("Badness");
    }
  }
  if(thisIndex == (int) numTreePieces - 1) {
    if(verbosity > 3)
      ckerr << "Cell interactions interval: " << cellInterval << endl;
    if(!xdr_setpos(&xdrs, FieldHeader::sizeBytes) || !xdr_template(&xdrs, &cellInterval.min) || !xdr_template(&xdrs, &cellInterval.max)) {
      ckerr << "TreePiece " << thisIndex << ": Error going back to write the cell interaction bounds, aborting" << endl;
      CkAbort("Badness");
    }
    if(verbosity > 2)
      ckerr << "TreePiece " << thisIndex << ": Wrote the cell interaction bounds" << endl;
  }
  xdr_destroy(&xdrs);
  fclose(outfile);

  outfile = fopen((basefilename + ".calls").c_str(), "r+b");
  fseek(outfile, 0, SEEK_END);
  xdrstdio_create(&xdrs, outfile, XDR_ENCODE);	
  for(unsigned int i = 1; i <= myNumParticles; ++i) {
    callsInterval.grow(myParticles[i].numEntryCalls);
    if(!xdr_template(&xdrs, &(myParticles[i].numEntryCalls))) {
      ckerr << "TreePiece " << thisIndex << ": Error writing entry calls to disk, aborting" << endl;
      CkAbort("Badness");
    }
  }
  if(thisIndex == (int) numTreePieces - 1) {
    if(verbosity > 3)
      ckerr << "Entry call interval: " << callsInterval << endl;
    if(!xdr_setpos(&xdrs, FieldHeader::sizeBytes) || !xdr_template(&xdrs, &callsInterval.min) || !xdr_template(&xdrs, &callsInterval.max)) {
      ckerr << "TreePiece " << thisIndex << ": Error going back to write the entry call bounds, aborting" << endl;
      CkAbort("Badness");
    }
    if(verbosity > 2)
      ckerr << "TreePiece " << thisIndex << ": Wrote the entry call bounds" << endl;
  }
  xdr_destroy(&xdrs);
  fclose(outfile);

  outfile = fopen((basefilename + ".partints").c_str(), "r+b");
  fseek(outfile, 0, SEEK_END);
  xdrstdio_create(&xdrs, outfile, XDR_ENCODE);	
  for(unsigned int i = 1; i <= myNumParticles; ++i) {
    particleInterval.grow(myParticles[i].numParticleInteractions);
    if(!xdr_template(&xdrs, &(myParticles[i].numParticleInteractions))) {
      ckerr << "TreePiece " << thisIndex << ": Error writing particle interactions to disk, aborting" << endl;
      CkAbort("Badness");
    }
  }
  if(thisIndex == (int) numTreePieces - 1) {
    if(verbosity > 3)
      ckerr << "Particle interactions interval: " << particleInterval << endl;
    if(!xdr_setpos(&xdrs, FieldHeader::sizeBytes) || !xdr_template(&xdrs, &particleInterval.min) || !xdr_template(&xdrs, &particleInterval.max)) {
      ckerr << "TreePiece " << thisIndex << ": Error going back to write the particle interaction bounds, aborting" << endl;
      CkAbort("Badness");
    }
    if(verbosity > 2)
      ckerr << "TreePiece " << thisIndex << ": Wrote the particle interaction bounds" << endl;
  }		
  xdr_destroy(&xdrs);
  fclose(outfile);
#endif
	
  if(thisIndex != (int) numTreePieces - 1)
    pieces[thisIndex + 1].outputStatistics(macInterval, cellInterval, particleInterval, callsInterval, totalmass, cb);
  if(thisIndex == (int) numTreePieces - 1) cb.send();
}

void TreePiece::outputRelativeErrors(Interval<double> errorInterval, const CkCallback& cb) {
  if(thisIndex == 0) {
    if(verbosity > 2)
      ckerr << "TreePiece " << thisIndex << ": Writing header for errors file" << endl;
    FILE* outfile = fopen((basefilename + ".error").c_str(), "wb");
    XDR xdrs;
    xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
    fh.code = float64;
    fh.dimensions = 1;
    if(!xdr_template(&xdrs, &fh) || !xdr_template(&xdrs, &errorInterval.min) || !xdr_template(&xdrs, &errorInterval.max)) {
      ckerr << "TreePiece " << thisIndex << ": Could not write header to errors file, aborting" << endl;
      CkAbort("Badness");
    }
    xdr_destroy(&xdrs);
    fclose(outfile);
  }
	
  if(verbosity > 3)
    ckerr << "TreePiece " << thisIndex << ": Writing my errors to disk" << endl;
	
  FILE* outfile = fopen((basefilename + ".error").c_str(), "r+b");
  fseek(outfile, 0, SEEK_END);
  XDR xdrs;
  xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
	
  double error;
	
  for(unsigned int i = 1; i <= myNumParticles; ++i) {
    error = (myParticles[i].treeAcceleration - myParticles[i].acceleration).length() / myParticles[i].acceleration.length();
    errorInterval.grow(error);
    if(!xdr_template(&xdrs, &error)) {
      ckerr << "TreePiece " << thisIndex << ": Error writing errors to disk, aborting" << endl;
      CkAbort("Badness");
    }
  }
	
  if(thisIndex == (int) numTreePieces - 1) {
    if(!xdr_setpos(&xdrs, FieldHeader::sizeBytes) || !xdr_template(&xdrs, &errorInterval.min) || !xdr_template(&xdrs, &errorInterval.max)) {
      ckerr << "TreePiece " << thisIndex << ": Error going back to write the error bounds, aborting" << endl;
      CkAbort("Badness");
    }
    if(verbosity > 2)
      ckerr << "TreePiece " << thisIndex << ": Wrote the error bounds" << endl;
    ckerr << "Error Bounds:" << errorInterval.min << ", "
	 << errorInterval.max << endl;
    cb.send();
  }
	
  xdr_destroy(&xdrs);
  fclose(outfile);
	
  if(thisIndex != (int) numTreePieces - 1)
    pieces[thisIndex + 1].outputRelativeErrors(errorInterval, cb);
}

/// @TODO Fix pup routine to handle correctly the tree
void TreePiece::pup(PUP::er& p) {
  ckout << "TreePiece " << thisIndex << ": Getting PUP'd!" << endl;
  CBase_TreePiece::pup(p);
  p | numTreePieces;
  p | callback;
  p | myNumParticles;
  if(p.isUnpacking()) {
    myParticles = new GravityParticle[myNumParticles + 2];
    //leftBoundary = myParticles;
    //rightBoundary = myParticles + myNumParticles + 1;
  }
  for(int i=0;i<myNumParticles+2;i++){
    p |myParticles[i];
  }
  p | numSplitters;
  if(p.isUnpacking())
    splitters = new Key[numSplitters];
  p(splitters, numSplitters);
  p | pieces;
  //p | streamingProxy;
  p | basefilename;
  p | boundingBox;
  p | fh;
  p | started;
  p | iterationNo;
  if(p.isUnpacking()){
    switch (useTree) {
    case Binary_Oct:
      root = new BinaryTreeNode(1, Tree::Boundary, 0, myNumParticles+1, 0);
      break;
    case Oct_Oct:
      //root = new OctTreeNode(1, Tree::Boundary, 0, myNumParticles+1, 0);
      break;
    default:
      CkAbort("We should have never reached here!");
    }
  }
  p | (*root);
  if(p.isUnpacking()){
    nodeLookupTable[root->getKey()]=root;
  }
  p | boundaryNodesPending;
  //p | tempNode;
  p | theta;
  p | mySerialNumber;
  p | myNumParticlesPending;
  p | numBuckets;
  p | currentBucket;
#if COSMO_STATS > 0
  p | myNumParticleInteractions;
  p | myNumCellInteractions;
  p | myNumMACChecks;
  p | piecemass;
#endif
  if(p.isUnpacking()){
    localCache=cacheManagerProxy.ckLocalBranch();
  }
  if(!(p.isUnpacking())) {
	
    //Pack nodeLookup here
    int num=0;
    for (NodeLookupType::iterator iter=nodeLookupTable.begin();iter!=nodeLookupTable.end();iter++){
      if(iter->second != root && iter->second != NULL){
	num++;
      }	
    }
    p(num);
    for (NodeLookupType::iterator iter=nodeLookupTable.begin();iter!=nodeLookupTable.end();iter++){
      if(iter->second != root && iter->second != NULL){
	Key k = iter->first;
	p | k;
	p | (*(iter->second));
      }	
    }
  }else{
    int num;
    p(num);
    for(int i=0;i<num;i++){
      Key k;
      GenericTreeNode *n = root->createNew();
      p | k;
      p | *n;
      nodeLookupTable[k] = n;
      if(n->getType() == Bucket){
	bucketList.push_back(n);
      }
    }
    int count=0;
    rebuildSFCTree(root,NULL,&count);
    sort(bucketList.begin(),bucketList.end(),compBucket);
    if(verbosity)
			CkPrintf("[%d] TreePiece %d bucketList size %d numBuckets %d nodelookupsize %d count %d\n",CkMyPe(),thisIndex,bucketList.size(),numBuckets,num,count);
  }
  //p | unfilledRequests;
}

void TreePiece::rebuildSFCTree(GenericTreeNode *node,GenericTreeNode *parent,int *count){
  if(node == NULL){
    return;
  }
  (*count)++;
  node->parent = (GenericTreeNode *)parent;
  for (unsigned int i=0; i<node->numChildren(); ++i) {
    GenericTreeNode *child = nodeLookupTable[node->getChildKey(i)];
    switch (useTree) {
    case Binary_Oct:
      ((BinaryTreeNode*)node)->children[i] = (BinaryTreeNode*)child;
      break;
    case Oct_Oct:
      ((OctTreeNode*)node)->children[i] = (OctTreeNode*)child;
      break;
    default:
      CkAbort("We should have never reached here!");
    }
    rebuildSFCTree(child,node,count);
  }
}
bool compBucket(GenericTreeNode *ln,GenericTreeNode *rn){
  return (ln->firstParticle < rn->firstParticle);
}


/** Check that all the particles in the tree are really in their boxes.
    Because the keys are made of only the first 21 out of 23 bits of the
    floating point representation, there can be particles that are outside
    their box by tiny amounts.  Whether this is bad is not yet known. */
void TreePiece::checkTree(GenericTreeNode* node) {
  if(node->getType() == Empty) return;
  if(node->getType() == Bucket) {
    for(unsigned int iter = node->firstParticle; iter <= node->lastParticle; ++iter) {
      if(!node->boundingBox.contains(myParticles[iter].position)) {
	ckerr << "Not in the box: Box: " << node->boundingBox << " Position: " << myParticles[iter].position << "\nNode key: " << keyBits(node->getKey(), 63).c_str() << "\nParticle key: " << keyBits(myParticles[iter].key, 63).c_str() << endl;
      }
    }
  } else if(node->getType() != NonLocal && node->getType() != NonLocalBucket) {
    GenericTreeNode* childIterator;
    for(unsigned int i = 0; i < node->numChildren(); ++i) {
      childIterator = node->getChildren(i);
      if(childIterator)
	checkTree(childIterator);
    }
  }
}

/// Color a node
string getColor(GenericTreeNode* node) {
  ostringstream oss;
  switch(node->getType()) {
  case Bucket:
  case Internal:
    oss << "black";
    break;
  case NonLocal:
  case NonLocalBucket:
    oss << "red";
    break;
  case Boundary:
    oss << "purple";
    break;
  default:
    oss << "yellow";
  }
  return oss.str();
}

/// Make a label for a node
string makeLabel(GenericTreeNode* node) {
  ostringstream oss;
  oss << keyBits(node->getKey(), 63) << "\\n";
  switch(node->getType()) {
  case Invalid:
    oss << "Invalid";
    break;
  case Bucket:
    //oss << "Bucket: " << (node->endParticle - node->beginParticle) << " particles";
    oss << "Bucket";
    break;
  case Internal:
    oss << "Internal";
    break;
  case NonLocal:
    oss << "NonLocal: Chare " << node->remoteIndex;
    break;
  case NonLocalBucket:
    oss << "NonLocalBucket: Chare " << node->remoteIndex;
    break;
  case Empty:
    oss << "Empty";
    break;
  case Boundary:
    oss << "Boundary: Total N " << node->remoteIndex;
    break;
  case Top:
    oss << "Top";
    break;
  default:
    oss << "Unknown NodeType!";
  }
  return oss.str();
}

/// Print a graphviz version of a tree
void TreePiece::printTree(GenericTreeNode* node, ostream& os) {
  if(node == 0)
    return;
	
  string nodeID = keyBits(node->getKey(), 63);
  os << nodeID << " ";
  //os << "\tnode [color=\"" << getColor(node) << "\"]\n";
  //os << "\t\"" << nodeID << "\" [label=\"" << makeLabel(node) << "\\nCM: " << (node->moments.cm) << "\\nM: " << node->moments.totalMass << "\\nN_p: " << (node->endParticle - node->beginParticle) << "\\nOwners: " << node->numOwners << "\"]\n";
  //os << "\t\"" << nodeID << "\" [label=\"" << makeLabel(node) << "\\nLocal N: " << (node->endParticle - node->beginParticle) << "\\nOwners: " << node->numOwners << "\"]\n";
  //os << "\t\"" << nodeID << "\" [label=\"" << keyBits(node->getKey(), 63) << "\\n";
  int first, last;
  switch(node->getType()) {
  case Bucket:
    os << "Bucket: Size=" << (node->lastParticle - node->firstParticle + 1) << "(" << node->firstParticle << "-" << node->lastParticle << ")";
    break;
  case Internal:
    os << "Internal: Size=" << (node->lastParticle - node->firstParticle + 1) << "(" << node->firstParticle << "-" << node->lastParticle << ")";
    break;
  case NonLocal:
    //os << "NonLocal: Chare=" << node->remoteIndex << "\\nRemote N under: " << (node->lastParticle - node->firstParticle + 1) << "\\nOwners: " << node->numOwners;
    nodeOwnership(node, first, last);
    os << "NonLocal: Chare=" << node->remoteIndex << ", Owners=" << first << "-" << last;
    break;
  case NonLocalBucket:
    //os << "NonLocal: Chare=" << node->remoteIndex << "\\nRemote N under: " << (node->lastParticle - node->firstParticle + 1) << "\\nOwners: " << node->numOwners;
    nodeOwnership(node, first, last);
    CkAssert(first == last);
    os << "NonLocalBucket: Chare=" << node->remoteIndex << ", Owner=" << first << ", Size=" << node->particleCount << "(" << node->firstParticle << "-" << node->lastParticle << ")";
    break;
  case Boundary:
    nodeOwnership(node, first, last);
    os << "Boundary: Totalsize=" << node->particleCount << ", Localsize=" << (node->lastParticle - node->firstParticle) << "(" << node->firstParticle + (node->firstParticle==0?1:0) << "-" << node->lastParticle - (node->lastParticle==myNumParticles+1?1:0) << "), Owners=" << first << "-" << last;
    break;
  case Empty:
    os << "Empty "<<node->remoteIndex;
    break;
  }
  if (node->getType() == Bucket || node->getType() == Internal || node->getType() == Boundary || node->getType() == NonLocal || node->getType() == NonLocalBucket) 
    os << " V "<<node->moments.radius<<" "<<node->moments.soft<<" "<<node->moments.cm.x<<" "<<node->moments.cm.y<<" "<<node->moments.cm.z<<" "<<node->moments.xx<<" "<<node->moments.xy<<" "<<node->moments.xz<<" "<<node->moments.yy<<" "<<node->moments.yz<<" "<<node->moments.zz;

  os << "\n";
	
  //if(node->parent)
  //  os << "\t\"" << keyBits(node->parent->getKey(), 63) << "\" -> \"" << nodeID << "\";\n";
	
  if(node->getType() == NonLocal || node->getType() == NonLocalBucket || node->getType() == Bucket || node->getType() == Empty)
    return;

  GenericTreeNode* childIterator;
  for(unsigned int i = 0; i < node->numChildren(); ++i) {
    childIterator = node->getChildren(i);
    if(childIterator)
      printTree(childIterator, os);
    else {
      os << "\tnode [color=\"green\"]\n";
      os << "\t\"" << nodeID << i << "\" [label=\"None\"]\n";
      os << "\t\"" << nodeID << "\" -> \"" << nodeID << i << "\";\n";
    }
  }
}

/// Write a file containing a graphviz dot graph of my tree
void TreePiece::report(const CkCallback& cb) {
  ostringstream outfilename;
  outfilename << "tree_" << thisIndex << ".dot";
  ofstream os(outfilename.str().c_str());

  os << "digraph G" << thisIndex << " {\n";
  os << "\tcenter = \"true\"\n";
  os << "\tsize = \"7.5,10\"\n";
  //os << "\tratio = \"fill\"\n";
  //os << "\tfontname = \"Courier\"\n";
  os << "\tnode [style=\"bold\"]\n";
  os << "\tlabel = \"Piece: " << thisIndex << "\\nParticles: " 
     << myNumParticles << "\"\n";
  /*	os << "\tlabel = \"Piece: " << thisIndex << "\\nParticles: " 
	<< myNumParticles << "\\nLeft Splitter: " << keyBits(myParticles[0].key, 63)
	<< "\\nLeftmost Key: " << keyBits(myParticles[1].key, 63) 
	<< "\\nRightmost Key: " << keyBits(myParticles[myNumParticles].key, 63) 
	<< "\\nRight Splitter: " << keyBits(myParticles[myNumParticles + 1].key, 63) << "\";\n";
  */
  os << "\tfontname = \"Helvetica\"\n";
  printTree(root, os);
  os << "}" << endl;
	
  os.close();
	
  //checkTree(root);
		
  contribute(0, 0, CkReduction::concat, cb);
}

void TreePiece::getPieceValues(piecedata *totaldata){
#if COSMO_STATS > 0
  totaldata->modifypiecedata(myNumCellInteractions,myNumParticleInteractions,myNumMACChecks,piecemass);
  if(thisIndex != (int) numTreePieces - 1)
  	pieces[thisIndex + 1].getPieceValues(totaldata);
  else {
    CkCallback& cb= totaldata->getcallback();
    cb.send(totaldata);
  }
#endif
}
