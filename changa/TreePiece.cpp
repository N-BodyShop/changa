/** @file TreePiece.cpp
 */

#include <cstdio>
#include <algorithm>
#include <fstream>

//#include "ComlibManager.h"

#include "ParallelGravity.h"

#include "Space.h"

using namespace std;
using namespace SFC;
using namespace TreeStuff;
using namespace TypeHandling;

int TreeStuff::maxBucketSize;

void TreePiece::load(const std::string& fn, const CkCallback& cb) {
	basefilename = fn;
	
	//read in particles
	XDR xdrs;
	FILE* infile = fopen((basefilename + ".mass").c_str(), "rb");
	if(!infile) {
		cerr << "TreePiece " << thisIndex << ": Couldn't open masses file, aborting" << endl;
		CkAbort("Badness");
	}
	xdrstdio_create(&xdrs, infile, XDR_DECODE);
	
	if(!xdr_template(&xdrs, &fh)) {
		cerr << "TreePiece " << thisIndex << ": Couldn't read header from masses file, aborting" << endl;
		CkAbort("Badness");
	}
	
	if(fh.magic != FieldHeader::MagicNumber || fh.dimensions != 1 || fh.code != float32) {
		cerr << "TreePiece " << thisIndex << ": Masses file is corrupt or of incorrect type, aborting" << endl;
		CkAbort("Badness");
	}
		
	myNumParticles = fh.numParticles / numTreePieces;
	unsigned int excess = fh.numParticles % numTreePieces;
	unsigned int startParticle = myNumParticles * thisIndex;
	if(thisIndex < excess) {
		myNumParticles++;
		startParticle += thisIndex;
	} else
		startParticle += excess;
	
	if(verbosity > 2)
		cerr << "TreePiece " << thisIndex << ": Of " << fh.numParticles << " particles, taking " << startParticle << " through " << (startParticle + myNumParticles - 1) << endl;

	myParticles = new GravityParticle[myNumParticles + 2];
	
	float mass;
	float maxMass;
	if(!xdr_template(&xdrs, &mass) || !xdr_template(&xdrs, &maxMass)) {
		cerr << "TreePiece " << thisIndex << ": Problem reading beginning of the mass file, aborting" << endl;
		CkAbort("Badness");
	}
	
	if(mass == maxMass) { //all the same mass
		for(u_int64_t i = 0; i < myNumParticles; ++i)
			myParticles[i + 1].mass = mass;
	} else {

		if(!seekField(fh, &xdrs, startParticle)) {
			cerr << "TreePiece " << thisIndex << ": Could not seek to my part of the mass file, aborting" << endl;
			CkAbort("Badness");
		}

		for(u_int64_t i = 0; i < myNumParticles; ++i) {
			if(!xdr_template(&xdrs, &mass)) {
				cerr << "TreePiece " << thisIndex << ": Problem reading my part of the mass file, aborting" << endl;
				CkAbort("Badness");
			}
			myParticles[i + 1].mass = mass;
		}
	}
	
	xdr_destroy(&xdrs);
	fclose(infile);
	
	infile = fopen((basefilename + ".pos").c_str(), "rb");
	if(!infile) {
		cerr << "TreePiece " << thisIndex << ": Couldn't open positions file, aborting" << endl;
		CkAbort("Badness");
	}
	xdrstdio_create(&xdrs, infile, XDR_DECODE);
	
	FieldHeader posHeader;
	if(!xdr_template(&xdrs, &posHeader)) {
		cerr << "TreePiece " << thisIndex << ": Couldn't read header from positions file, aborting" << endl;
		CkAbort("Badness");
	}
	
	if(posHeader.magic != FieldHeader::MagicNumber || posHeader.dimensions != 3 || posHeader.code != float32) {
		cerr << "TreePiece " << thisIndex << ": Positions file is corrupt or of incorrect type, aborting" << endl;
		CkAbort("Badness");
	}
	
	if(posHeader.time != fh.time || posHeader.numParticles != fh.numParticles) {
		cerr << "TreePiece " << thisIndex << ": Positions file doesn't match masses file, aborting" << endl;
		CkAbort("Badness");
	}
	
	Vector3D<float> pos;
	Vector3D<float> maxPos;
	if(!xdr_template(&xdrs, &pos) || !xdr_template(&xdrs, &maxPos)) {
		cerr << "TreePiece " << thisIndex << ": Problem reading beginning of the positions file, aborting" << endl;
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
		if(!seekField(posHeader, &xdrs, startParticle)) {
			cerr << "TreePiece " << thisIndex << ": Could not seek to my part of the positions file, aborting" << endl;
			CkAbort("Badness");
		}

		//read all my particles' positions and make keys
		for(u_int64_t i = 0; i < myNumParticles; ++i) {
			if(!xdr_template(&xdrs, &pos)) {
				cerr << "TreePiece " << thisIndex << ": Problem reading my part of the positions file, aborting" << endl;
				CkAbort("Badness");
			}
			myParticles[i + 1].position = pos;
			myParticles[i + 1].key = generateKey(pos, boundingBox);
		}
	}
	
	xdr_destroy(&xdrs);
	fclose(infile);
	
	if(verbosity > 2)
		cerr << "TreePiece " << thisIndex << ": Read in masses and positions" << endl;
	
	contribute(0, 0, CkReduction::concat, cb);
}

void TreePiece::buildTree(int bucketSize, const CkCallback& cb) {
	maxBucketSize = bucketSize;
	callback = cb;
	Key bounds[2];
	bounds[0] = myParticles[1].key;
	bounds[1] = myParticles[myNumParticles].key;
	contribute(2 * sizeof(Key), bounds, CkReduction::concat, CkCallback(CkIndex_TreePiece::collectSplitters(0), thisArrayID));
}

void TreePiece::collectSplitters(CkReductionMsg* m) {
	numSplitters = 2 * numTreePieces;
	splitters = new Key[numSplitters];
	Key* splits = static_cast<Key *>(m->getData());
	copy(splits, splits + numSplitters, splitters);
	sort(splitters, splitters + numSplitters);
	contribute(0, 0, CkReduction::concat, CkCallback(CkIndex_TreePiece::startTreeBuild(0), thisArrayID));
	delete m;
	if(verbosity > 3)
		cerr << "TreePiece " << thisIndex << ": Collected splitters" << endl;
}

void TreePiece::startTreeBuild(CkReductionMsg* m) {
	delete m;
	
	if(thisIndex == 0)
		myParticles[0].key = firstPossibleKey;
	else
		myParticles[0].key = splitters[2 * thisIndex - 1];
	
	if(thisIndex == numTreePieces - 1)
		myParticles[myNumParticles + 1].key = lastPossibleKey;
	else
		myParticles[myNumParticles + 1].key = splitters[2 * thisIndex + 2];
	
	leftBoundary = myParticles;
	rightBoundary = myParticles + myNumParticles + 1;
	
	root = new SFCTreeNode;
	root->key = firstPossibleKey;
	root->boundingBox = boundingBox;
	nodeLookup[root->lookupKey()] = root;
	numBuckets = 0;
	bucketList.clear();
	
	boundaryNodesPending = 0;
	
	if(verbosity > 3)
		cerr << "TreePiece " << thisIndex << ": Starting tree build" << endl;
	
	buildTree(root, leftBoundary, rightBoundary);
	
	if(boundaryNodesPending == 0)
		contribute(0, 0, CkReduction::concat, callback);
	
	if(verbosity > 3)
		cerr << "TreePiece " << thisIndex << ": Number of buckets: " << numBuckets << endl;
	if(verbosity > 3)
		cerr << "TreePiece " << thisIndex << ": Finished tree build, resolving boundary nodes" << endl;
}

/// Find what chare this node's left child resides on, and create it
inline SFCTreeNode* TreePiece::lookupLeftChild(SFCTreeNode* node) {
	SFCTreeNode* child = node->createLeftChild();
	nodeLookup[child->lookupKey()] = child;
	child->setType(NonLocal);
	tempNode.key = node->key;
	tempNode.level = node->level + 1;
	if(!nodeOwnership(&tempNode, &child->remoteIndex, &child->numOwners))
		cerr << "This is surprising, but may get taken care of." << endl;
	return child;
}

inline SFCTreeNode* TreePiece::lookupRightChild(SFCTreeNode* node) {
	SFCTreeNode* child = node->createRightChild();
	nodeLookup[child->lookupKey()] = child;
	child->setType(NonLocal);
	tempNode.key = node->rightChildKey();
	tempNode.level = node->level + 1;
	if(!nodeOwnership(&tempNode, &child->remoteIndex, &child->numOwners))
		cerr << "This is surprising, but may get taken care of." << endl;
	return child;
}

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
			cerr << "Wow, I didn't think this could happen.  Live and learn." << endl;
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

/** A recursive algorithm for building my tree.
 Examines successive bits in the particles' keys, looking for splits.
 Each bit is a level of nodes in the tree.  We keep going down until
 we can bucket the particles.  The left and right boundaries of this
 piece of tree will point to other pieces on other chares in the array.
 */
void TreePiece::buildTree(SFCTreeNode* node, GravityParticle* leftParticle, GravityParticle* rightParticle) {
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
			calculateRadiusFarthestCorner(node->moments, node->boundingBox);
			bucketList.push_back(node);
			numBuckets++;
			return;
		}
	} else if(node->level == 63) {
		cerr << thisIndex << ": TreePiece: This piece of tree has exhausted all the bits in the keys.  Super double-plus ungood!" << endl;
		cerr << "Left particle: " << (leftParticle - myParticles) << " Right particle: " << (rightParticle - myParticles) << endl;
		cerr << "Left key : " << keyBits(leftParticle->key, 63) << endl;
		cerr << "Right key: " << keyBits(rightParticle->key, 63) << endl;
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
			if(thisIndex != numTreePieces - 1) //the right-most chare can't point any further right
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
		cerr << "Bits not right: " << leftBit << " vs " << rightBit << endl;
		cerr << "Left particle: " << (leftParticle - myParticles) << " Right particle: " << (rightParticle - myParticles) << endl;
		cerr << "Left key : " << keyBits(leftParticle->key, 63) << endl;
		cerr << "Right key: " << keyBits(rightParticle->key, 63) << endl;
		return;
	} else { //both zeros, make a left child
		child = node->createLeftChild();
		nodeLookup[child->lookupKey()] = child;
		buildTree(child, leftParticle, rightParticle);
		//should the right child be remote?
		if(rightParticle == rightBoundary && thisIndex != numTreePieces - 1)
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
	if((leftParticle == leftBoundary && thisIndex != 0) || (rightParticle == rightBoundary && thisIndex != numTreePieces - 1)) {
		//send my information about this node to the designated owner
		unsigned int designatedOwner;
		nodeOwnership(node, &designatedOwner, &node->numOwners);
		boundaryNodesPending++;
		//in boundary nodes, remoteIndex contains the total number of particles in the node, from all the co-owners
		//endParticle - beginParticle tells you how many particles are actually on this processor
		node->remoteIndex = node->endParticle - node->beginParticle;
		if(designatedOwner != thisIndex) //don't send yourself your contribution
			pieces[designatedOwner].acceptBoundaryNodeContribution(node->lookupKey(), node->remoteIndex, node->moments);
		node->setType(Boundary);
	} else {
		node->numOwners = 1;
		node->setType(Internal);
		calculateRadiusFarthestCorner(node->moments, node->boundingBox);
	}
}

void TreePiece::acceptBoundaryNodeContribution(const Key lookupKey, const u_int64_t numParticles, const MultipoleMoments& moments) {
	NodeLookupType::iterator nodeIter = nodeLookup.find(lookupKey);
	
	if(nodeIter == nodeLookup.end()) {
		cerr << "TreePiece " << thisIndex << ": Well crap, how the hell did this happen?" << endl;
		cerr << "Key: " << keyBits(lookupKey, 63) << endl;
		return;
	}
	
	SFCTreeNode* node = nodeIter->second;
	//merge new information
	node->remoteIndex += numParticles;
	node->moments += moments;
	//decrement number of contributions
	node->numOwners--;
	//if done, send final information to all co-owners
	if(node->numOwners == 1) {
		calculateRadiusFarthestCorner(node->moments, node->boundingBox);
		unsigned int firstOwner, lastOwner;
		//recalculate number of owners, get co-owners
		nodeOwnership(node, 0, &node->numOwners, &firstOwner, &lastOwner);
		//send moments, numParticles to all co-owners
		for(unsigned int i = firstOwner; i <= lastOwner; ++i)
			pieces[i].acceptBoundaryNode(lookupKey, node->remoteIndex, node->moments);
	}
}

void TreePiece::acceptBoundaryNode(const Key lookupKey, const u_int64_t numParticles, const MultipoleMoments& moments) {
	NodeLookupType::iterator nodeIter = nodeLookup.find(lookupKey);
	
	if(nodeIter == nodeLookup.end()) {
		cerr << "Well crap, how the hell did this happen, especially now?" << endl;
		return;
	}
	
	SFCTreeNode* node = nodeIter->second;
	if(node->getType() != Boundary)
		cerr << "How does this work?" << endl;
	//merge final information
	node->remoteIndex = numParticles;
	node->moments = moments;	
	boundaryNodesPending--;
	if(boundaryNodesPending == 0) {
		calculateRemoteMoments(root);
		contribute(0, 0, CkReduction::concat, callback);
	}
}

void TreePiece::calculateRemoteMoments(SFCTreeNode* node) {
	if(node->getType() == NonLocal) {
		SFCTreeNode* sibling = dynamic_cast<SFCTreeNode *>(node->getSibling());
		SFCTreeNode* parent = dynamic_cast<SFCTreeNode *>(node->parent);
		node->beginParticle = 0;
		if(sibling->getType() == Boundary)
			node->endParticle = parent->remoteIndex - sibling->remoteIndex;
		else
			node->endParticle = parent->remoteIndex - (sibling->endParticle - sibling->beginParticle);
		if(node->endParticle != 0) {
			node->moments = parent->moments - sibling->moments;
			calculateRadiusFarthestCorner(node->moments, node->boundingBox);
		} else {
			nodeLookup.erase(node->lookupKey());
			if(node->isLeftChild())
				parent->leftChild = 0;
			else
				parent->rightChild = 0;
			delete node;
		}
	} else if(node->getType() == Boundary) {
		if(node->leftChild)
			calculateRemoteMoments(dynamic_cast<SFCTreeNode *>(node->leftChild));
		if(node->rightChild)
			calculateRemoteMoments(dynamic_cast<SFCTreeNode *>(node->rightChild));
	}
}

void TreePiece::calculateGravityDirect(const CkCallback& cb) {
	callback = cb;
	
	GravityRequest req;
	req.requestingPieceIndex = thisIndex;
	
	myNumParticlesPending = myNumParticles * numTreePieces;
	
	for(u_int64_t i = 1; i <= myNumParticles; ++i) {
		req.identifier = i;
		req.position = myParticles[i].position;
		myParticles[i].acceleration = 0;
		pieces.fillRequestDirect(req);
	}
	
	started = true;
}

void TreePiece::fillRequestDirect(GravityRequest req) {
	Vector3D<double> r;
	double rsq;
	for(u_int64_t i = 1; i <= myNumParticles; ++i) {
		r = myParticles[i].position - req.position;
		rsq = r.lengthSquared();
		if(rsq != 0)
			req.acceleration += myParticles[i].mass * r / rsq / sqrt(rsq);
	}
	streamingProxy[req.requestingPieceIndex].receiveGravityDirect(req);
}

void TreePiece::receiveGravityDirect(const GravityRequest& req) {
	myParticles[req.identifier].acceleration += req.acceleration;
	if(started && --myNumParticlesPending == 0) {
		started = false;
		contribute(0, 0, CkReduction::concat, callback);
	}
}

void TreePiece::startNextParticle() {
	if(nextParticle > myNumParticles)
		return;
	GravityRequest req;
	req.startingNode = root->lookupKey();
	req.requestingPieceIndex = thisIndex;
	req.identifier = nextParticle;
	req.position = myParticles[nextParticle].position;
	myParticles[nextParticle].treeAcceleration = 0;
	streamingProxy[thisIndex].fillRequestTree(req);
	nextParticle++;
}

void TreePiece::calculateGravityTree(double t, const CkCallback& cb) {
	callback = cb;
	theta = t;
	mySerialNumber = 0;
	myNumCellInteractions = 0;
	myNumParticleInteractions = 0;
	myNumMACChecks = 0;
	myNumProxyCalls = 0;
	myNumProxyCallsBack = 0;
	
	nextParticle = 1;
	startNextParticle();
	
	myNumParticlesPending = myNumParticles;
	started = true;
	
}

void TreePiece::fillRequestTree(GravityRequest req) {
	//double startTime = CkWallTimer();
	//lookup starting node using startingNode key
	NodeLookupType::iterator nodeIter = nodeLookup.find(req.startingNode);
	if(nodeIter == nodeLookup.end()) {
		cerr << "Well crap, how the hell did this happen here?" << endl;
		return;
	}
	SFCTreeNode* node = nodeIter->second;
	
	//make the request ready to go in the queue
	req.numAdditionalRequests = 1;
	req.acceleration = 0;
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

void TreePiece::walkTree(GravityTreeNode* node, GravityRequest& req) {
	req.numMACChecks++;
	myNumMACChecks++;
	Vector3D<double> r = node->moments.cm / node->moments.totalMass - req.position;
	double rsq = r.lengthSquared();
	double r_open = opening_geometry_factor * node->moments.radius / theta;
	if(rsq > r_open * r_open) {
		req.numCellInteractions++;
		myNumCellInteractions++;
		req.acceleration += node->moments.totalMass * r / rsq / sqrt(rsq);
		rsq = req.acceleration.length();
	} else if(node->getType() == Bucket) {
		req.numParticleInteractions += node->endParticle - node->beginParticle;
		myNumParticleInteractions += node->endParticle - node->beginParticle;
		for(unsigned int i = node->beginParticle; i < node->endParticle; ++i) {
			r = myParticles[i].position - req.position;
			rsq = r.lengthSquared();
			if(rsq != 0)
				req.acceleration += myParticles[i].mass * r / rsq / sqrt(rsq);
		}
	} else if(node->getType() == NonLocal) {
		unfilledRequests[mySerialNumber].numAdditionalRequests++;
		req.numEntryCalls++;
		req.startingNode = dynamic_cast<SFCTreeNode *>(node)->lookupKey();
		streamingProxy[node->remoteIndex].fillRequestTree(req);
		myNumProxyCalls++;
	} else {
		GenericTreeNode** childrenIterator = node->getChildren();
		for(unsigned int i = 0; i < node->numChildren(); ++i) {
			if(childrenIterator[i])
				walkTree(dynamic_cast<GravityTreeNode *>(childrenIterator[i]), req);
		}
	}
}

void TreePiece::receiveGravityTree(const GravityRequest& req) {
	//lookup request
	UnfilledRequestsType::iterator requestIter = unfilledRequests.find(req.identifier);
	if(requestIter == unfilledRequests.end()) {
		cerr << "Well crap, how the hell did this happen here and now?" << endl;
		cerr << "TreePiece " << thisIndex << ": Got request from " << req.requestingPieceIndex << " with id " << req.identifier << endl;
		return;
	}
	GravityRequest& request = requestIter->second;
	request.merge(req);
	if(--request.numAdditionalRequests == 0) {
		if(request.requestingPieceIndex == thisIndex) {
			myNumParticlesPending--;
			//this request originated here, it's for one of my particles
			myParticles[request.identifier].update(request);
		} else {
			streamingProxy[request.requestingPieceIndex].receiveGravityTree(request);
			myNumProxyCallsBack++;
		}
		
		unfilledRequests.erase(requestIter);
		if(started && myNumParticlesPending == 0) {
			started = false;
			contribute(0, 0, CkReduction::concat, callback);
			cout << "TreePiece " << thisIndex << ": Made " << myNumProxyCalls << " proxy calls forward, " << myNumProxyCallsBack << " to respond" << endl;
			if(verbosity > 4)
				cerr << "TreePiece " << thisIndex << ": My particles are done" << endl;
		}
	}
}

void TreePiece::startNextBucket() {
	if(currentBucket >= numBuckets)
		return;
	
	SFCTreeNode* node = bucketList[currentBucket++];
	unsigned int numParticlesInBucket = node->endParticle - node->beginParticle;
	BucketGravityRequest req(numParticlesInBucket);
	req.startingNode = root->lookupKey();
	req.identifier = node->beginParticle;
	req.requestingPieceIndex = thisIndex;
	for(unsigned int i = node->beginParticle; i < node->endParticle; ++i) {
		req.positions[i - node->beginParticle] = myParticles[i].position;
		req.boundingBox.grow(myParticles[i].position);
		myParticles[i].treeAcceleration = 0;
	}
	streamingProxy[thisIndex].fillRequestBucketTree(req);
	myNumProxyCalls++;
}

void TreePiece::calculateGravityBucketTree(double t, const CkCallback& cb) {
	callback = cb;
	theta = t;
	mySerialNumber = 0;
	myNumProxyCalls = 0;
	myNumProxyCallsBack = 0;
	
	cout << "TreePiece " << thisIndex << ": I have " << numBuckets << " buckets" << endl;
	
	currentBucket = 0;
	startNextBucket();
	
	myNumParticlesPending = myNumParticles;
	started = true;
}

void TreePiece::fillRequestBucketTree(BucketGravityRequest req) {
	NodeLookupType::iterator nodeIter = nodeLookup.find(req.startingNode);
	if(nodeIter == nodeLookup.end()) {
		cerr << "Well crap, how the hell did this happen here?" << endl;
		return;
	}
	SFCTreeNode* node = nodeIter->second;
	
	//make the request ready to go in the queue
	req.numAdditionalRequests = 1;
	for(unsigned int i = 0; i < req.numParticlesInBucket; ++i)
		req.accelerations[i] = 0;
	
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

void TreePiece::walkBucketTree(GravityTreeNode* node, BucketGravityRequest& req) {
	myNumMACChecks++;
	Vector3D<double> cm(node->moments.cm / node->moments.totalMass);
	Vector3D<double> r;
	double rsq;
	Sphere<double> s(cm, opening_geometry_factor * node->moments.radius / theta);
	if(!Space::intersect(req.boundingBox, s)) {
		myNumCellInteractions += req.numParticlesInBucket;
		for(unsigned int i = 0; i < req.numParticlesInBucket; ++i) {
			r = cm - req.positions[i];
			rsq = r.lengthSquared();
			req.accelerations[i] += node->moments.totalMass * r / rsq / sqrt(rsq);
		}
	} else if(node->getType() == Bucket) {
		myNumParticleInteractions += req.numParticlesInBucket * (node->beginParticle - node->endParticle);
		for(unsigned int i = node->beginParticle; i < node->endParticle; ++i) {
			for(unsigned int j = 0; j < req.numParticlesInBucket; ++j) {
				r = myParticles[i].position - req.positions[j];
				rsq = r.lengthSquared();
				if(rsq != 0)
					req.accelerations[j] += myParticles[i].mass * r / rsq / sqrt(rsq);
			}
		}
	} else if(node->getType() == NonLocal) {
		unfilledBucketRequests[mySerialNumber].numAdditionalRequests++;
		req.startingNode = dynamic_cast<SFCTreeNode *>(node)->lookupKey();
		streamingProxy[node->remoteIndex].fillRequestBucketTree(req);
		myNumProxyCalls++;
	} else {
		GenericTreeNode** childrenIterator = node->getChildren();
		for(unsigned int i = 0; i < node->numChildren(); ++i) {
			if(childrenIterator[i])
				walkBucketTree(dynamic_cast<GravityTreeNode *>(childrenIterator[i]), req);
		}
	}
}
/*
void TreePiece::cachedWalkBucketTree(GravityTreeNode* node, BucketGravityRequest& req) {
	myNumMACChecks++;
	Vector3D<double> cm(node->moments.cm / node->moments.totalMass);
	Vector3D<double> r;
	double rsq;
	Sphere<double> s(cm, opening_geometry_factor * node->moments.radius / theta);
	if(!Space::intersect(req.boundingBox, s)) {
		myNumCellInteractions += req.numParticlesInBucket;
		for(unsigned int i = 0; i < req.numParticlesInBucket; ++i) {
			r = cm - req.positions[i];
			rsq = r.lengthSquared();
			req.accelerations[i] += node->moments.totalMass * r / rsq / sqrt(rsq);
		}
	} else if(node->getType() == Bucket) {
		myNumParticleInteractions += req.numParticlesInBucket * (node->beginParticle - node->endParticle);
		for(unsigned int i = node->beginParticle; i < node->endParticle; ++i) {
			for(unsigned int j = 0; j < req.numParticlesInBucket; ++j) {
				r = myParticles[i].position - req.positions[j];
				rsq = r.lengthSquared();
				if(rsq != 0)
					req.accelerations[j] += myParticles[i].mass * r / rsq / sqrt(rsq);
			}
		}
	} else if(node->getType() == NonLocal) {
		Key lookupKey = dynamic_cast<SFCTreeNode *>(node)->lookupKey();
		GravityTreeNode* cachedNode = localCache->requestData(thisIndex, lookupKey, node->remoteIndex);
		if(cachedNode) //hit the cache, keep going
			cachedWalkBucketTree(cachedNode, req);
		else { //missed the cache, remember this request
			unfilledBucketRequests[mySerialNumber].numAdditionalRequests++;
			missedRequests.insert(make_pair(lookupKey, mySerialNumber));
		}
	} else {
		GenericTreeNode** childrenIterator = node->getChildren();
		for(unsigned int i = 0; i < node->numChildren(); ++i) {
			if(childrenIterator[i])
				cachedWalkBucketTree(dynamic_cast<GravityTreeNode *>(childrenIterator[i]), req);
		}
	}
}

void TreePiece::receiveRequestedData(Key requestID, GravityTreeNode* node) {
	typedef multimap<unsigned int, BucketGravityRequest>::iterator MI;
	pair<MI, MI> requests = missedRequests.equal_range(requestID);
	for(MI iter = requests.first; iter != requests.second; ++iter)
		cachedWalkBucketTree(node, *iter);
	missedRequests.erase(requestID);
}
*/
void TreePiece::receiveGravityBucketTree(const BucketGravityRequest& req) {
	//lookup request
	UnfilledBucketRequestsType::iterator requestIter = unfilledBucketRequests.find(req.identifier);
	if(requestIter == unfilledBucketRequests.end()) {
		cerr << "Well crap, how the hell did this happen here and now?" << endl;
		cerr << "TreePiece " << thisIndex << ": Got request from " << req.requestingPieceIndex << " with id " << req.identifier << endl;
		return;
	}
	BucketGravityRequest& request = requestIter->second;
	if(request.numParticlesInBucket != req.numParticlesInBucket)
		cerr << "How could this be?" << endl;
	request.merge(req);
	if(--request.numAdditionalRequests == 0) {
		if(request.requestingPieceIndex == thisIndex) {
			myNumParticlesPending -= request.numParticlesInBucket;
			//this request originated here, it's for one of my particles
			for(unsigned int i = 0; i < request.numParticlesInBucket; ++i)
				myParticles[request.identifier + i].treeAcceleration += request.accelerations[i];
		} else {
			streamingProxy[request.requestingPieceIndex].receiveGravityBucketTree(request);
			myNumProxyCallsBack++;
		}
		
		unfilledBucketRequests.erase(requestIter);
		if(started && myNumParticlesPending == 0) {
			started = false;
			contribute(0, 0, CkReduction::concat, callback);
			cout << "TreePiece " << thisIndex << ": Made " << myNumProxyCalls << " proxy calls forward, " << myNumProxyCallsBack << " to respond" << endl;
			if(verbosity > 4)
				cerr << "TreePiece " << thisIndex << ": My particles are done" << endl;
		}
	}
}

void TreePiece::outputAccelerations(OrientedBox<double> accelerationBox, const string& suffix, const CkCallback& cb) {
	if(thisIndex == 0) {
		if(verbosity > 2)
			cerr << "TreePiece " << thisIndex << ": Writing header for accelerations file" << endl;
		FILE* outfile = fopen((basefilename + "." + suffix).c_str(), "wb");
		XDR xdrs;
		xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
		fh.code = float64;
		fh.dimensions = 3;
		if(!xdr_template(&xdrs, &fh) || !xdr_template(&xdrs, &accelerationBox.lesser_corner) || !xdr_template(&xdrs, &accelerationBox.greater_corner)) {
			cerr << "TreePiece " << thisIndex << ": Could not write header to accelerations file, aborting" << endl;
			CkAbort("Badness");
		}
		xdr_destroy(&xdrs);
		fclose(outfile);
	}
	
	if(verbosity > 3)
		cerr << "TreePiece " << thisIndex << ": Writing my accelerations to disk" << endl;
	
	FILE* outfile = fopen((basefilename + "." + suffix).c_str(), "r+b");
	fseek(outfile, 0, SEEK_END);
	XDR xdrs;
	xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
	
	for(unsigned int i = 1; i <= myNumParticles; ++i) {
		accelerationBox.grow(myParticles[i].acceleration);
		if(!xdr_template(&xdrs, &(myParticles[i].acceleration))) {
			cerr << "TreePiece " << thisIndex << ": Error writing accelerations to disk, aborting" << endl;
			CkAbort("Badness");
		}
	}
	
	if(thisIndex == numTreePieces - 1) {
		if(!xdr_setpos(&xdrs, FieldHeader::sizeBytes) || !xdr_template(&xdrs, &accelerationBox.lesser_corner) || !xdr_template(&xdrs, &accelerationBox.greater_corner)) {
			cerr << "TreePiece " << thisIndex << ": Error going back to write the acceleration bounds, aborting" << endl;
			CkAbort("Badness");
		}
		if(verbosity > 2)
			cerr << "TreePiece " << thisIndex << ": Wrote the acceleration bounds" << endl;
		cb.send();
	}
	
	xdr_destroy(&xdrs);
	fclose(outfile);
	
	if(thisIndex != numTreePieces - 1)
		pieces[thisIndex + 1].outputAccelerations(accelerationBox, suffix, cb);
}

void TreePiece::outputStatistics(Interval<unsigned int> macInterval, Interval<unsigned int> cellInterval, Interval<unsigned int> particleInterval, Interval<unsigned int> callsInterval, const CkCallback& cb) {
	if(verbosity > 1) {
		cerr << "TreePiece " << thisIndex << ": Statistics\nMy number of MAC checks: " << myNumMACChecks << endl;
		cerr << "My number of particle-cell interactions: " << myNumCellInteractions << endl;
		cerr << "My number of particle-particle interactions: " << myNumParticleInteractions << endl;
	}
	
	if(thisIndex == 0) {
		macInterval.max = 0;
		macInterval.min = macInterval.max - 1;
		cellInterval = macInterval;
		particleInterval = macInterval;
		callsInterval = macInterval;
		
		if(verbosity > 2)
			cerr << "TreePiece " << thisIndex << ": Writing headers for statistics files" << endl;
		fh.dimensions = 1;
		fh.code = uint32;
		FILE* outfile = fopen((basefilename + ".MACs").c_str(), "wb");
		XDR xdrs;
		xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
		
		unsigned int dummy;
		if(!xdr_template(&xdrs, &fh) || !xdr_template(&xdrs, &dummy) || !xdr_template(&xdrs, &dummy)) {
			cerr << "TreePiece " << thisIndex << ": Could not write header to MAC file, aborting" << endl;
			CkAbort("Badness");
		}
		xdr_destroy(&xdrs);
		fclose(outfile);
		
		outfile = fopen((basefilename + ".cellints").c_str(), "wb");
		xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
		if(!xdr_template(&xdrs, &fh) || !xdr_template(&xdrs, &dummy) || !xdr_template(&xdrs, &dummy)) {
			cerr << "TreePiece " << thisIndex << ": Could not write header to cell-interactions file, aborting" << endl;
			CkAbort("Badness");
		}
		xdr_destroy(&xdrs);
		fclose(outfile);

		outfile = fopen((basefilename + ".partints").c_str(), "wb");
		xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
		if(!xdr_template(&xdrs, &fh) || !xdr_template(&xdrs, &dummy) || !xdr_template(&xdrs, &dummy)) {
			cerr << "TreePiece " << thisIndex << ": Could not write header to particle-interactions file, aborting" << endl;
			CkAbort("Badness");
		}
		xdr_destroy(&xdrs);
		fclose(outfile);

		outfile = fopen((basefilename + ".calls").c_str(), "wb");
		xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
		if(!xdr_template(&xdrs, &fh) || !xdr_template(&xdrs, &dummy) || !xdr_template(&xdrs, &dummy)) {
			cerr << "TreePiece " << thisIndex << ": Could not write header to entry-point calls file, aborting" << endl;
			CkAbort("Badness");
		}
		xdr_destroy(&xdrs);
		fclose(outfile);
	}
	
	if(verbosity > 3)
		cerr << "TreePiece " << thisIndex << ": Writing my statistics to disk" << endl;
	
	FILE* outfile = fopen((basefilename + ".MACs").c_str(), "r+b");
	fseek(outfile, 0, SEEK_END);
	XDR xdrs;
	xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
	
	for(unsigned int i = 1; i <= myNumParticles; ++i) {
		macInterval.grow(myParticles[i].numMACChecks);
		if(!xdr_template(&xdrs, &(myParticles[i].numMACChecks))) {
			cerr << "TreePiece " << thisIndex << ": Error writing MAC checks to disk, aborting" << endl;
			CkAbort("Badness");
		}
	}
	
	if(thisIndex == numTreePieces - 1) {
		if(verbosity > 3)
			cerr << "MAC interval: " << macInterval << endl;
		if(!xdr_setpos(&xdrs, FieldHeader::sizeBytes) || !xdr_template(&xdrs, &macInterval.min) || !xdr_template(&xdrs, &macInterval.max)) {
			cerr << "TreePiece " << thisIndex << ": Error going back to write the MAC bounds, aborting" << endl;
			CkAbort("Badness");
		}
		if(verbosity > 2)
			cerr << "TreePiece " << thisIndex << ": Wrote the MAC bounds" << endl;
	}
	
	xdr_destroy(&xdrs);
	fclose(outfile);
	
	outfile = fopen((basefilename + ".cellints").c_str(), "r+b");
	fseek(outfile, 0, SEEK_END);
	xdrstdio_create(&xdrs, outfile, XDR_ENCODE);	
	for(unsigned int i = 1; i <= myNumParticles; ++i) {
		cellInterval.grow(myParticles[i].numCellInteractions);
		if(!xdr_template(&xdrs, &(myParticles[i].numCellInteractions))) {
			cerr << "TreePiece " << thisIndex << ": Error writing cell interactions to disk, aborting" << endl;
			CkAbort("Badness");
		}
	}
	if(thisIndex == numTreePieces - 1) {
		if(verbosity > 3)
			cerr << "Cell interactions interval: " << cellInterval << endl;
		if(!xdr_setpos(&xdrs, FieldHeader::sizeBytes) || !xdr_template(&xdrs, &cellInterval.min) || !xdr_template(&xdrs, &cellInterval.max)) {
			cerr << "TreePiece " << thisIndex << ": Error going back to write the cell interaction bounds, aborting" << endl;
			CkAbort("Badness");
		}
		if(verbosity > 2)
			cerr << "TreePiece " << thisIndex << ": Wrote the cell interaction bounds" << endl;
	}
	xdr_destroy(&xdrs);
	fclose(outfile);

	outfile = fopen((basefilename + ".calls").c_str(), "r+b");
	fseek(outfile, 0, SEEK_END);
	xdrstdio_create(&xdrs, outfile, XDR_ENCODE);	
	for(unsigned int i = 1; i <= myNumParticles; ++i) {
		callsInterval.grow(myParticles[i].numEntryCalls);
		if(!xdr_template(&xdrs, &(myParticles[i].numEntryCalls))) {
			cerr << "TreePiece " << thisIndex << ": Error writing entry calls to disk, aborting" << endl;
			CkAbort("Badness");
		}
	}
	if(thisIndex == numTreePieces - 1) {
		if(verbosity > 3)
			cerr << "Entry call interval: " << callsInterval << endl;
		if(!xdr_setpos(&xdrs, FieldHeader::sizeBytes) || !xdr_template(&xdrs, &callsInterval.min) || !xdr_template(&xdrs, &callsInterval.max)) {
			cerr << "TreePiece " << thisIndex << ": Error going back to write the entry call bounds, aborting" << endl;
			CkAbort("Badness");
		}
		if(verbosity > 2)
			cerr << "TreePiece " << thisIndex << ": Wrote the entry call bounds" << endl;
	}
	xdr_destroy(&xdrs);
	fclose(outfile);

	outfile = fopen((basefilename + ".partints").c_str(), "r+b");
	fseek(outfile, 0, SEEK_END);
	xdrstdio_create(&xdrs, outfile, XDR_ENCODE);	
	for(unsigned int i = 1; i <= myNumParticles; ++i) {
		particleInterval.grow(myParticles[i].numParticleInteractions);
		if(!xdr_template(&xdrs, &(myParticles[i].numParticleInteractions))) {
			cerr << "TreePiece " << thisIndex << ": Error writing particle interactions to disk, aborting" << endl;
			CkAbort("Badness");
		}
	}
	if(thisIndex == numTreePieces - 1) {
		if(verbosity > 3)
			cerr << "Particle interactions interval: " << particleInterval << endl;
		if(!xdr_setpos(&xdrs, FieldHeader::sizeBytes) || !xdr_template(&xdrs, &particleInterval.min) || !xdr_template(&xdrs, &particleInterval.max)) {
			cerr << "TreePiece " << thisIndex << ": Error going back to write the particle interaction bounds, aborting" << endl;
			CkAbort("Badness");
		}
		if(verbosity > 2)
			cerr << "TreePiece " << thisIndex << ": Wrote the particle interaction bounds" << endl;
		cb.send();
	}		
	xdr_destroy(&xdrs);
	fclose(outfile);
	
	if(thisIndex != numTreePieces - 1)
		pieces[thisIndex + 1].outputStatistics(macInterval, cellInterval, particleInterval, callsInterval, cb);
}

void TreePiece::outputRelativeErrors(Interval<double> errorInterval, const CkCallback& cb) {
	if(thisIndex == 0) {
		if(verbosity > 2)
			cerr << "TreePiece " << thisIndex << ": Writing header for errors file" << endl;
		FILE* outfile = fopen((basefilename + ".error").c_str(), "wb");
		XDR xdrs;
		xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
		fh.code = float64;
		fh.dimensions = 1;
		if(!xdr_template(&xdrs, &fh) || !xdr_template(&xdrs, &errorInterval.min) || !xdr_template(&xdrs, &errorInterval.max)) {
			cerr << "TreePiece " << thisIndex << ": Could not write header to errors file, aborting" << endl;
			CkAbort("Badness");
		}
		xdr_destroy(&xdrs);
		fclose(outfile);
	}
	
	if(verbosity > 3)
		cerr << "TreePiece " << thisIndex << ": Writing my errors to disk" << endl;
	
	FILE* outfile = fopen((basefilename + ".error").c_str(), "r+b");
	fseek(outfile, 0, SEEK_END);
	XDR xdrs;
	xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
	
	double error;
	
	for(unsigned int i = 1; i <= myNumParticles; ++i) {
		error = (myParticles[i].treeAcceleration - myParticles[i].acceleration).length() / myParticles[i].acceleration.length();
		errorInterval.grow(error);
		if(!xdr_template(&xdrs, &error)) {
			cerr << "TreePiece " << thisIndex << ": Error writing errors to disk, aborting" << endl;
			CkAbort("Badness");
		}
	}
	
	if(thisIndex == numTreePieces - 1) {
		if(!xdr_setpos(&xdrs, FieldHeader::sizeBytes) || !xdr_template(&xdrs, &errorInterval.min) || !xdr_template(&xdrs, &errorInterval.max)) {
			cerr << "TreePiece " << thisIndex << ": Error going back to write the error bounds, aborting" << endl;
			CkAbort("Badness");
		}
		if(verbosity > 2)
			cerr << "TreePiece " << thisIndex << ": Wrote the error bounds" << endl;
		cb.send();
	}
	
	xdr_destroy(&xdrs);
	fclose(outfile);
	
	if(thisIndex != numTreePieces - 1)
		pieces[thisIndex + 1].outputRelativeErrors(errorInterval, cb);
}

void TreePiece::pup(PUP::er& p) {
	cerr << "TreePiece " << thisIndex << ": Getting PUP'd!" << endl;
	p | numTreePieces;
	p | callback;
	p | myNumParticles;
	if(p.isUnpacking()) {
		myParticles = new GravityParticle[myNumParticles + 2];
		leftBoundary = myParticles;
		rightBoundary = myParticles + myNumParticles + 1;
	}
	p(myParticles, myNumParticles + 2);
	p | numSplitters;
	if(p.isUnpacking())
		splitters = new Key[numSplitters];
	p(splitters, numSplitters);
	p | pieces;
	p | basefilename;
	p | boundingBox;
	p | fh;
	p | started;
	//p | root;
	p | boundaryNodesPending;
	//p | tempNode;
	p | theta;
	p | mySerialNumber;
	p | myNumParticlesPending;
	
	if(p.isUnpacking()) {
		//XXX Rebuild nodeLookup here
		/*
		stack<TreeNode *> nodeStack;
		nodeStack.push(root);
		TreeNode* node;
		while(!nodeStack.empty()) {
			node = nodeStack.top();
			nodeStack.pop();
			nodeLookup[node->lookupKey()] = node;
			if(!node->isBucket() && !node->isNonLocal()) {
				if(node->rightChild)
					nodeStack.push(node->rightChild);
				if(node->leftChild)
					nodeStack.push(node->leftChild);
			}
		}
		*/
	}
	
	//p | unfilledRequests;
}

/** Check that all the particles in the tree are really in their boxes.
 Because the keys are made of only the first 21 out of 23 bits of the
 floating point representation, there can be particles that are outside
 their box by tiny amounts.  Whether this is bad is not yet known. */
void TreePiece::checkTree(SFCTreeNode* node) {
	if(node->getType() == Bucket) {
		for(unsigned int iter = node->beginParticle; iter != node->endParticle; ++iter) {
			if(!node->boundingBox.contains(myParticles[iter].position))
				cerr << "Not in the box: Box: " << node->boundingBox << " Position: " << myParticles[iter].position << "\nNode key: " << keyBits(node->key, node->level) << "\nParticle key: " << keyBits(myParticles[iter].key, 63) << endl;
		}
	} else if(node->getType() != NonLocal) {
		GenericTreeNode** childrenIterator = node->getChildren();
		for(unsigned int i = 0; i < node->numChildren(); ++i) {
			if(childrenIterator[i])
				checkTree(dynamic_cast<SFCTreeNode *>(childrenIterator[i]));
		}
	}
}

/// Color a node
string getColor(SFCTreeNode* node) {
	ostringstream oss;
	switch(node->getType()) {
		case Bucket:
		case Internal:
			oss << "black";
			break;
		case NonLocal:
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
string makeLabel(SFCTreeNode* node) {
	ostringstream oss;
	oss << keyBits(node->key, node->level) << "\\n";
	switch(node->getType()) {
		case Invalid:
			oss << "Invalid";
			break;
		case Bucket:
			oss << "Bucket: " << (node->endParticle - node->beginParticle) << " particles";
			break;
		case Internal:
			oss << "Internal";
			break;
		case NonLocal:
			oss << "NonLocal: Chare " << node->remoteIndex;
			break;
		case Empty:
			oss << "Empty";
			break;
		case Boundary:
			oss << "Boundary: Total " << node->remoteIndex;
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
void printTree(SFCTreeNode* node, ostream& os) {
	if(node == 0)
		return;
	
	string nodeID = keyBits(node->key, node->level);
	os << "\tnode [color=\"" << getColor(node) << "\"]\n";
	os << "\t\"" << nodeID << "\" [label=\"" << makeLabel(node) << "\\nCM: " << (node->moments.cm / node->moments.totalMass) << "\\nM: " << node->moments.totalMass << "\\nN_p: " << (node->endParticle - node->beginParticle) << "\\nOwners: " << node->numOwners << "\"]\n";
	if(node->parent)
		os << "\t\"" << keyBits(node->key, node->level - 1) << "\" -> \"" << nodeID << "\";\n";
	
	if(node->getType() == NonLocal || node->getType() == Bucket)
		return;

	GenericTreeNode** childrenIterator = node->getChildren();
	for(unsigned int i = 0; i < node->numChildren(); ++i) {
		if(childrenIterator[i])
			printTree(dynamic_cast<SFCTreeNode *>(childrenIterator[i]), os);
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
	os << "\tratio = \"fill\"\n";
	os << "\tfontname = \"Courier\"\n";
	os << "\tlabel = \"Piece: " << thisIndex << "\\nParticles: " 
			<< myNumParticles << "\\nLeft Splitter: " << keyBits(myParticles[0].key, 63)
			<< "\\nLeftmost Key: " << keyBits(myParticles[1].key, 63) 
			<< "\\nRightmost Key: " << keyBits(myParticles[myNumParticles].key, 63) 
			<< "\\nRight Splitter: " << keyBits(myParticles[myNumParticles + 1].key, 63) << "\";\n";
	os << "\tfontname = \"Helvetica\"\n";
	printTree(root, os);
	os << "}" << endl;
	
	os.close();
	
	checkTree(root);
		
	contribute(0, 0, CkReduction::concat, cb);
}
