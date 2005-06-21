/** @file TreePiece.cpp
 @author Graeme Lufkin (gwl@u.washington.edu)
 */

#include <sstream>
#include <algorithm>
#include <fstream>

#include "pup_stl.h"

#include "SFC.h"

#include "TreePiece.h"

#include "Reductions.h"
#include "Tree.h"
#include "DataManager.h"

using std::ostringstream;
using std::vector;
using std::map;
using std::ofstream;

using namespace std;

TreePiece::TreePiece() {
	dm = 0;
	numBoundaries = 0;
	treePieces = CProxy_TreePiece(thisArrayID);
	root = 0;
	myPlace = -1;
}

TreePiece::~TreePiece() {
	if(root)
		delete root;
}
/*
void TreePiece::pup(PUP::er& p) {
	p | numBoundaries;
	p | treePieces;
	p | dataManager;
	if(p.isUnpacking()) {
		dm = dataManager.ckLocalBranch();
		if(dm == 0)
			cerr << "Ack, unpacked and couldn't get the local DataManager" << endl;
	} else if(p.isPacking()) {
		//remove this TreePiece from DataManager's list
		dm->myTreePieces.erase(find(dm->myTreePieces.begin(), dm->myTreePieces.end(), thisIndex));
	}
	p | boundingBox;
	p | myPlace;
	p | callback;
	p | numTreePieces;
	p | myParticles;
	p | mySortedParticles;
	p | myNumParticles;
	p | root;
	p | leftSplitter;
	p | rightSplitter;
	if(p.isUnpacking()) {
		leftBoundary = myParticles.begin();
		rightBoundary = myParticles.end() - 1;
	}
	p | nodeLookup;
}
*/
void TreePiece::registerWithDataManager(const CkGroupID& dataManagerID, const CkCallback& cb) {
	dataManager = CProxy_DataManager(dataManagerID);
	dm = dataManager.ckLocalBranch();
	if(dm == 0) {
		cerr << thisIndex << ": TreePiece: Fatal: Couldn't register with my DataManger" << endl;
		cb.send(0);
		return;
	}
	
	dm->myTreePieces.push_back(thisIndex);

	contribute(0, 0, CkReduction::concat, cb);
}

void TreePiece::loadParticles(const std::string& filename, const int numPieces, const CkCallback& cb) {
	callback = cb;
	
	TipsyReader r(filename);
	if(!r.status()) {
		cerr << thisIndex << ": TreePiece: Fatal: Couldn't open tipsy file!" << endl;
		cb.send(0);
		return;
	}
	
	numTreePieces = numPieces;
	header h = r.getHeader();
	totalNumParticles = h.nbodies;
	myNumParticles = totalNumParticles / numTreePieces;
	unsigned int excess = totalNumParticles % numTreePieces;
	unsigned int startParticle = myNumParticles * thisIndex;
	if(thisIndex < excess) {
		myNumParticles++;
		startParticle += thisIndex;
	} else
		startParticle += excess;
	
	if(verbosity > 2)
		cerr << thisIndex << ": TreePiece: Taking " << myNumParticles << " of " << totalNumParticles << " particles, starting at " << startParticle << endl;

	myParticles.reserve(myNumParticles);
	
	if(!r.seekParticleNum(startParticle)) {
		cerr << thisIndex << ": TreePiece: Fatal: Couldn't seek to my particles!" << endl;
		cb.send(0);
		return;
	}
	
	gas_particle gp;
	dark_particle dp;
	star_particle sp;
	for(int i = 0; i < myNumParticles; ++i) {
		if(i + startParticle < h.nsph) {
			r.getNextGasParticle(gp);
			myParticles.push_back(gp);
		} else if(i + startParticle < h.nsph + h.ndark) {
			r.getNextDarkParticle(dp);
			myParticles.push_back(dp);
		} else {
			r.getNextStarParticle(sp);
			myParticles.push_back(sp);
		}
		boundingBox.grow(myParticles.back().position);
	}
	
	contribute(sizeof(OrientedBox<double>), &boundingBox, growOrientedBox_double, CkCallback(CkIndex_TreePiece::assignKeys(0), treePieces));
}
/*
/// Accept particles from the DataManager
void TreePiece::receiveParticles(const FullParticle* particles, const int n, const CkCallback& cb) {
	myParticles.resize(n);
	copy(particles, particles + n, myParticles.begin());
	mySortedParticles.clear();
	
	if(myParticles.size() == 0) {
		cerr << thisIndex << ": TreePiece: Fatal: I don't have any particles!" << endl;
		cb.send(0);
		return;
	}
	
	//now get my bounding box
	boundingBox.lesser_corner = boundingBox.greater_corner = myParticles.front().position;
	for(vector<FullParticle>::iterator iter = myParticles.begin() + 1; iter != myParticles.end(); ++iter)
		boundingBox.grow(iter->position);
	
	callback = cb;
	//reduce it to find whole simulation bounding box
	contribute(sizeof(boundingBox), &boundingBox, boxGrowthReduction, CkCallback(CkIndex_TreePiece::assignKeys(0), thisArrayID));
	
	if(verbosity >= 5)
		cout << thisIndex << ": TreePiece: Received " << n << " particles from the DataManager" << endl;
}
*/

/// After the bounding box has been found, we can assign keys to the particles
void TreePiece::assignKeys(CkReductionMsg* m) {
	if(m->getSize() != sizeof(OrientedBox<double>)) {
		cerr << thisIndex << ": TreePiece: Fatal: Wrong size reduction message received!" << endl;
		callback.send(0);
		delete m;
		return;
	}
	
	boundingBox = *static_cast<OrientedBox<double> *>(m->getData());
	delete m;
	if(thisIndex == 0 && verbosity)
		cerr << "TreePiece: Bounding box originally: " << boundingBox << endl;
	OrientedBox<double> unitCube(Vector3D<double>(-0.5, -0.5, -0.5), Vector3D<double>(0.5, 0.5, 0.5));
	double diagonal = boundingBox.size().length();
	double epsilon = max(1E-2, 1.0 / totalNumParticles);
	if((boundingBox.lesser_corner - unitCube.lesser_corner).length() / diagonal < epsilon
			&& (boundingBox.greater_corner - unitCube.greater_corner).length() / diagonal < epsilon) {
		boundingBox = unitCube;
	} else {
		OrientedBox<float> fbox(boundingBox);
		cubize(fbox);
		bumpBox(fbox, HUGE_VAL);
		boundingBox = fbox;
	}
	if(thisIndex == 0 && verbosity)
		cerr << "TreePiece: Bounding box resized to: " << boundingBox << endl;
	//give particles keys, using bounding box to scale
	for(vector<FullParticle>::iterator iter = myParticles.begin(); iter != myParticles.end(); ++iter)
		iter->key = generateKey(iter->position, boundingBox);
	sort(myParticles.begin(), myParticles.end());
	
	contribute(0, 0, CkReduction::concat, callback);
	
	if(verbosity >= 5)
		cout << thisIndex << ": TreePiece: Assigned keys to all my particles" << endl;
}

/// Determine my part of the sorting histograms by counting the number of my particles in each bin
void TreePiece::evaluateBoundaries(const CkCallback& cb) {
	int numBins = dm->boundaryKeys.size() - 1;
	//this array will contain the number of particles I own in each bin
	myBinCounts.assign(numBins, 0);
	vector<Key>::const_iterator endKeys = dm->boundaryKeys.end();
	vector<FullParticle>::iterator binBegin = myParticles.begin();
	vector<FullParticle>::iterator binEnd;
	FullParticle dummy;
	vector<int>::iterator binIter = myBinCounts.begin();
	vector<Key>::iterator keyIter = dm->boundaryKeys.begin();
	for(++keyIter; keyIter != endKeys; ++keyIter, ++binIter) {
		dummy.key = *keyIter;
		//find the last place I could put this splitter key in my array of particles
		binEnd = upper_bound(binBegin, myParticles.end(), dummy);
		//this tells me the number of particles between the last two splitter keys
		*binIter = (binEnd - binBegin);
		if(myParticles.end() <= binEnd)
			break;
		binBegin = binEnd;
	}
	
	//send my bin counts back in a reduction
	contribute(numBins * sizeof(int), &(*myBinCounts.begin()), CkReduction::sum_int, cb);
}

/// Once final splitter keys have been decided, I need to give my particles out to the TreePiece responsible for them
void TreePiece::unshuffleParticles(CkReductionMsg* m) {
	callback = *static_cast<CkCallback *>(m->getData());
	
	//find my responsibility
	myPlace = find(dm->responsibleIndex.begin(), dm->responsibleIndex.end(), thisIndex) - dm->responsibleIndex.begin();
	//assign my bounding keys
	leftSplitter = dm->boundaryKeys[myPlace];
	rightSplitter = dm->boundaryKeys[myPlace + 1];
	mySortedParticles.clear();
	mySortedParticles.reserve(dm->particleCounts[myPlace]);
	
	vector<Key>::iterator iter = dm->boundaryKeys.begin();
	vector<Key>::const_iterator endKeys = dm->boundaryKeys.end();
	vector<int>::iterator responsibleIter = dm->responsibleIndex.begin();
	FullParticle dummy;
	vector<FullParticle>::iterator binBegin = myParticles.begin();
	vector<FullParticle>::iterator binEnd;
	for(++iter; iter != endKeys; ++iter, ++responsibleIter) {
		dummy.key = *iter;
		//find particles between this and the last key
		binEnd = upper_bound(binBegin, myParticles.end(), dummy);
		//if I have any particles in this bin, send them to the responsible TreePiece
		if((binEnd - binBegin) > 0) {
			if(*responsibleIter == thisIndex)
				acceptSortedParticles(&(*binBegin), binEnd - binBegin);
			else
				treePieces[*responsibleIter].acceptSortedParticles(&(*binBegin), binEnd - binBegin);
		}
		if(myParticles.end() <= binEnd)
			break;
		binBegin = binEnd;
	}
	//resize myParticles so we can fit the sorted particles and the boundary particles
	myParticles.resize(dm->particleCounts[myPlace] + 2);
}

/*		
void TreePiece::unshuffleParticles(CkReductionMsg* m) {		
	callback = *static_cast<CkCallback *>(m->getData());
	
	//find my responsibility
	myPlace = find(dm->responsibleIndex.begin(), dm->responsibleIndex.end(), thisIndex) - dm->responsibleIndex.begin();
	
	//assign my bounding keys
	leftSplitter = dm->boundaryKeys[myPlace];
	rightSplitter = dm->boundaryKeys[myPlace + 1];
	
	if(numTreePieces > 1) {
		//rotate particles until we're up front
		FullParticle* myFirstParticle = upper_bound(myParticles.begin(), myParticles.end(), FullParticle(leftSplitter));
		rotate(myParticles.begin(), myFirstParticle, myParticles.end());

		//figure out how much extra space to make
		int extraSpace = 2 * myNumParticles / numTreePieces;

		//make space by moving everybody else over
		FullParticle* pastLastParticle = upper_bound(myParticles.begin(), myParticles.end(), FullParticle(rightSplitter));
		myParticles.insert(pastLastParticle, extraSpace, FullParticle());

		beginFreeParticles = myParticles.begin() + myBinCounts[thisIndex];
		beginOthersParticles = beginFreeParticles + extraSpace;
		swapsMade = 0;
		
		giveParticles((thisIndex + 1) % numTreePieces);
	} else
		contribute(0, 0, CkReduction::concat, CkCallback(CkIndex_TreePiece::shareBoundaries(0), thisArrayID));
}

void TreePiece::giveParticles(const int toWhom) {
	int howMany = myBinCounts[toWhom];
	treePieces[toWhom].getParticles(beginOthersParticles, howMany);
	beginOthersParticles += howMany;
}

void TreePiece::getParticles(const FullParticle* particles, int n) {
	if(n > (beginOthersParticles - beginFreeParticles))
		cerr << "Ack!  Fatal" << endl;
	copy(particles, particles + n, beginFreeParticles);
	beginFreeParticles += n;
	
	if(++swapsMade < numTreePieces - 1)
		giveParticles((swapsMade + thisIndex) % numTreePieces);
	else {
		myParticles.resize(beginFreeParticles - myParticles.begin());
		sort(myParticles.begin(), myParticles.end());
		contribute(0, 0, CkReduction::concat, CkCallback(CkIndex_TreePiece::shareBoundaries(0), thisArrayID));
	}
}

void TreePiece::canYouAcceptParticles(const int numParticles, const int piece) {
	if(slopSpace - alreadyPromisedSpace < numParticles) {
		alreadyPromisedSpace += numParticles;
		treePieces[piece].acknowledgeParticles(true);
	} else
		treePieces[piece].acknowledgeParticles(false);
}

void TreePiece::acknowledgeParticles(bool possible) {
	if(possible) {
		treePieces[pushPiece].acceptParticles(beginSendParticle, numParticlesToSend);
		
		//move all unsent particles to the end
		
		//get beginSendParticle and numParticlesToSend off the queue
	} else {
		//put pushPiece and numParticlesToSend back on queue
		//get new pushPiece
	}
}
*/
/// Accept particles from other TreePieces once the sorting has finished
void TreePiece::acceptSortedParticles(const FullParticle* particles, const int n) {
	copy(particles, particles + n, back_inserter(mySortedParticles));
	if(myPlace == -1)
		myPlace = find(dm->responsibleIndex.begin(), dm->responsibleIndex.end(), thisIndex) - dm->responsibleIndex.begin();
	if(dm->particleCounts[myPlace] == mySortedParticles.size()) { //I've got all my particles
		sort(mySortedParticles.begin(), mySortedParticles.end());
		//signify completion with a reduction
		contribute(0, 0, CkReduction::concat, CkCallback(CkIndex_TreePiece::shareBoundaries(0), thisArrayID));
	}
}

/// Give the TreePieces to the left and right of me my boundary keys so we can build correct trees.
void TreePiece::shareBoundaries(CkReductionMsg* m) {
	if(numTreePieces == 1) {
		acceptBoundaryKey(firstPossibleKey);
		acceptBoundaryKey(lastPossibleKey);	
	} else if(myPlace == 0) {
		treePieces[dm->responsibleIndex[myPlace + 1]].acceptBoundaryKey(mySortedParticles.back().key);
		acceptBoundaryKey(firstPossibleKey);
	} else if(myPlace == numTreePieces - 1) {
		treePieces[dm->responsibleIndex[myPlace - 1]].acceptBoundaryKey(mySortedParticles.front().key);
		acceptBoundaryKey(lastPossibleKey);
	} else {
		treePieces[dm->responsibleIndex[myPlace - 1]].acceptBoundaryKey(mySortedParticles.front().key);
		treePieces[dm->responsibleIndex[myPlace + 1]].acceptBoundaryKey(mySortedParticles.back().key);
	}
	delete m;
}

/// Receive my neighbor's boundary key
void TreePiece::acceptBoundaryKey(const Key k) {
	FullParticle dummy;
	dummy.key = k;
	if(k < mySortedParticles.front().key)
		myParticles.front() = dummy;
	else
		myParticles.back() = dummy;
	
	if(++numBoundaries == 2) {
		numBoundaries = 0;
		copy(mySortedParticles.begin(), mySortedParticles.end(), myParticles.begin() + 1);
		leftBoundary = &(*myParticles.begin());
		rightBoundary = &(*myParticles.end()) - 1;
		//I don't own my left and right boundaries, but I need to know them
		myNumParticles = myParticles.size() - 2;
		if(myNumParticles != dm->particleCounts[myPlace])
			cerr << thisIndex << ": TreePiece: Screwed up particle counts!" << endl;
		contribute(0, 0, CkReduction::concat, callback);
	}
}

/// I've got my sorted particles, now start building my tree
void TreePiece::startTreeBuild(const CkCallback& cb) {
	//I don't need this copy of my particles anymore
	mySortedParticles.clear();
	
	//create a root TreeNode
	root = new TreeNode;
	root->key = firstPossibleKey;
	root->box = boundingBox;
	//store this node in my lookup table
	nodeLookup[root->lookupKey()] = root;
	//build it
	buildTree(root, leftBoundary, rightBoundary);
	
	//signify completion with a reduction
	contribute(0, 0, CkReduction::concat, cb);
}

/// Find what chare this node's left child resides on
inline TreeNode* TreePiece::lookupLeftChild(TreeNode* node) {
	//find left child's boundaries
	Key leftNodeBoundary = node->key;
	Key rightNodeBoundary = node->key | ((static_cast<Key>(1) << (62 - node->level)) - 1);
	//we have some extra info, because we have the immediate left key, not just the split
	if(leftBoundary->key <= rightNodeBoundary)
		rightNodeBoundary = leftBoundary->key;
	//find the last place the left boundary of the node could go in the splitters array
	vector<Key>::iterator locLeft = upper_bound(dm->boundaryKeys.begin(), dm->boundaryKeys.end(), leftNodeBoundary);
	//find the first place the right boundary of the node could go in the splitters array
	vector<Key>::iterator locRight = lower_bound(locLeft, dm->boundaryKeys.end(), rightNodeBoundary);
	//pick the index in the middle of the region of chares that hold the node
	// XXX : Perhaps this should be biased to close chares for deep nodes?
	int pickedIndex = (locLeft - 1 + (locRight - locLeft) / 2) - dm->boundaryKeys.begin();
	TreeNode* child = node->createLeftChild();
	nodeLookup[child->lookupKey()] = child;
	child->setType(NonLocal);
	//get the chareID for the chosen region
	child->chareID = dm->responsibleIndex[pickedIndex];
	return child;
}

inline TreeNode* TreePiece::lookupRightChild(TreeNode* node) {
	Key leftNodeBoundary = node->key | (static_cast<Key>(1) << (62 - node->level));
	Key rightNodeBoundary = leftNodeBoundary | ((static_cast<Key>(1) << (62 - node->level)) - 1);
	if(leftNodeBoundary <= rightBoundary->key)
		leftNodeBoundary = rightBoundary->key;
	vector<Key>::iterator locLeft = upper_bound(dm->boundaryKeys.begin(), dm->boundaryKeys.end(), leftNodeBoundary);
	vector<Key>::iterator locRight = lower_bound(locLeft, dm->boundaryKeys.end(), rightNodeBoundary);
	int pickedIndex = (locLeft - 1 + (locRight - locLeft) / 2) - dm->boundaryKeys.begin();
	TreeNode* child = node->createRightChild();
	nodeLookup[child->lookupKey()] = child;
	child->setType(NonLocal);
	child->chareID = dm->responsibleIndex[pickedIndex];
	return child;
}

/** A recursive algorithm for building my tree.
 Examines successive bits in the particles' keys, looking for splits.
 Each bit is a level of nodes in the tree.  We keep going down until
 we can bucket the particles.  The left and right boundaries of this
 piece of tree will point to other pieces on other chares in the array.
 */
void TreePiece::buildTree(TreeNode* node, FullParticle* leftParticle, FullParticle* rightParticle) {
	
	//check if we should bucket these particles
	if(rightParticle - leftParticle <= TreeNode::maxBucketSize) {
		//can't bucket until we've cut at the boundary
		if((leftParticle != leftBoundary) && (rightParticle != rightBoundary)) {
			node->setType(Bucket);
			node->beginBucket = leftParticle;
			node->endBucket = rightParticle + 1;
			return;
		}
	} else if(node->level == 63) {
		cerr << thisIndex << ": TreePiece: This piece of tree has exhausted all the bits in the keys.  Super double-plus ungood!" << endl;
		return;
	}
	
	//this is the bit we are looking at
	Key currentBitMask = static_cast<Key>(1) << (62 - node->level);
	//we need to know the bit values at the left and right
	Key leftBit = leftParticle->key & currentBitMask;
	Key rightBit = rightParticle->key & currentBitMask;
	TreeNode* child;
	
	if((leftParticle == leftBoundary && myPlace != 0) || (rightParticle == rightBoundary && myPlace != numTreePieces - 1))
		node->setType(Boundary);
	else
		node->setType(Internal);
	
	if(leftBit ^ rightBit) { //a split at this level
		//find the split by looking for where the key with the bit switched on could go
		FullParticle* splitParticle = lower_bound(leftParticle, rightParticle + 1, node->key | currentBitMask);
		if(splitParticle == leftBoundary + 1) {
			//we need to make the left child point to a remote chare
			if(myPlace != 0) //the left-most chare can't point any further left
				lookupLeftChild(node);
			child = node->createRightChild();
			nodeLookup[child->lookupKey()] = child;
			buildTree(child, splitParticle, rightParticle);
		} else if(splitParticle == rightBoundary) {
			//we need to make the right child point to a remote chare
			child = node->createLeftChild();
			nodeLookup[child->lookupKey()] = child;
			buildTree(child, leftParticle, splitParticle - 1);
			if(myPlace != numTreePieces - 1) //the right-most chare can't point any further right
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
		if(leftParticle == leftBoundary && myPlace != 0)
			lookupLeftChild(node);
		child = node->createRightChild();
		nodeLookup[child->lookupKey()] = child;
		buildTree(child, leftParticle, rightParticle);
	} else { //both zeros, make a left child
		child = node->createLeftChild();
		nodeLookup[child->lookupKey()] = child;
		buildTree(child, leftParticle, rightParticle);
		//should the right child be remote?
		if(rightParticle == rightBoundary && myPlace != numTreePieces - 1)
			lookupRightChild(node);
	}
	
}

/* //The old tree-building function
void TreePiece::buildTree(TreeNode* node, FullParticle* leftParticle, FullParticle* rightParticle) {
	
	//check if we should bucket these particles
	if(rightParticle - leftParticle <= TreeNode::maxBucketSize) {
		//can't bucket until we've cut at the boundary
		if((leftParticle != leftBoundary) && (rightParticle != rightBoundary)) {
			node->setType(Bucket);
			node->beginBucket = leftParticle;
			node->endBucket = rightParticle + 1;
			return;
		}
	}
	
	//this is the bit we are looking at
	Key currentBitMask = static_cast<Key>(1) << (62 - node->level);
	//we need to know the bit values at the left and right
	Key leftBit = leftParticle->key & currentBitMask;
	Key rightBit = rightParticle->key & currentBitMask;
	
	node->setType(Boundary);
	if(leftBit ^ rightBit) { //a split at this level
		if(leftBit > rightBit)
			cerr << thisIndex << ": How the hell did this happen?" << endl;
		//find the index where the bit changes
		FullParticle* splitParticle = lower_bound(leftParticle, rightParticle + 1, node->key | currentBitMask);
		if(splitParticle - leftBoundary > 1) {
			node->setType(Internal);
			TreeNode* child = node->createLeftChild();
			nodeLookup[child->lookupKey()] = child;
			buildTree(child, leftParticle, splitParticle - 1);
		} else { //left child is non-local
			//find appropriate chare Id
			Key tail = (static_cast<Key>(1) << (62 - node->level)) - static_cast<Key>(1);
			vector<Key>::iterator lowLocation = upper_bound(dm->boundaryKeys.begin(), dm->boundaryKeys.end(), node->key & ~tail);
			vector<Key>::iterator hiLocation = upper_bound(dm->boundaryKeys.begin(), dm->boundaryKeys.end(), node->key | tail);
			int bin = ((lowLocation - dm->boundaryKeys.begin() - 1) + (hiLocation - dm->boundaryKeys.begin() - 1)) / 2;
			if(dm->responsibleIndex[bin] != thisIndex) {
				TreeNode* child = node->createLeftChild();
				child->setType(NonLocal);
				nodeLookup[child->lookupKey()] = child;
				child->chareID = dm->responsibleIndex[bin];
			} else
				node->setType(Internal);
		}
		if(rightBoundary > splitParticle) {
			node->setType(Internal);
			TreeNode* child = node->createRightChild();
			buildTree(child, splitParticle, rightParticle);
		} else { //right child is non-local
			//find appropriate chare Id
			Key tail = (static_cast<Key>(1) << (62 - node->level)) - static_cast<Key>(1);
			Key bit = (static_cast<Key>(1) << (62 - node->level));
			vector<Key>::iterator lowLocation = upper_bound(dm->boundaryKeys.begin(), dm->boundaryKeys.end(), (node->key | bit) & ~tail);
			vector<Key>::iterator hiLocation = upper_bound(dm->boundaryKeys.begin(), dm->boundaryKeys.end(), node->key | bit | tail);
			int bin = ((lowLocation - dm->boundaryKeys.begin() - 1) + (hiLocation - dm->boundaryKeys.begin() - 1)) / 2;
			if(dm->responsibleIndex[bin] == thisIndex) {
				if(((node->key | bit | tail) < rightSplitter) || (myPlace == numTreePieces - 1))
					node->setType(Internal);
			else {
					TreeNode* child = node->createRightChild();
					child->setType(NonLocal);
					nodeLookup[child->lookupKey()] = child;
					child->chareID = dm->responsibleIndex[bin + 1];
				}
			} else {
				TreeNode* child = node->createRightChild();
				child->setType(NonLocal);
				nodeLookup[child->lookupKey()] = child;
				child->chareID = dm->responsibleIndex[bin];
			}
		}
	
	} else if(leftBit & rightBit) { //both ones, make a right child
		TreeNode *child = node->createRightChild();
		nodeLookup[child->lookupKey()] = child;
		buildTree(child, leftParticle, rightParticle);
		if(node->key <= leftBoundary->key) {
			//find appropriate chare ID
			Key tail = (static_cast<Key>(1) << (62 - node->level)) - static_cast<Key>(1);
			vector<Key>::iterator lowLocation = upper_bound(dm->boundaryKeys.begin(), dm->boundaryKeys.end(), node->key & ~tail);
			vector<Key>::iterator hiLocation = upper_bound(dm->boundaryKeys.begin(), dm->boundaryKeys.end(), node->key | tail);
			int bin = ((lowLocation - dm->boundaryKeys.begin() - 1) + (hiLocation - dm->boundaryKeys.begin() - 1)) / 2;
			if(dm->responsibleIndex[bin] != thisIndex) {
				TreeNode *child = node->createLeftChild();
				child->setType(NonLocal);
				nodeLookup[child->lookupKey()] = child;
				child->chareID = dm->responsibleIndex[bin];
			} else
				node->setType(Internal);
		} else
			node->setType(Internal);
	} else { //both zeros, make a left child
		TreeNode* child = node->createLeftChild();
		nodeLookup[child->lookupKey()] = child;
		buildTree(child, leftParticle, rightParticle);
		currentBitMask >>= 1;
		if(((node->key | currentBitMask) >> (62 - node->level)) >= (rightBoundary->key >> (62 - node->level))) {
			//find appropriate chare Id
			Key tail = (static_cast<Key>(1) << (62 - node->level)) - static_cast<Key>(1);
			Key bit = (static_cast<Key>(1) << (62 - node->level));
			vector<Key>::iterator lowLocation = upper_bound(dm->boundaryKeys.begin(), dm->boundaryKeys.end(), (node->key | bit) & ~tail);
			vector<Key>::iterator hiLocation = upper_bound(dm->boundaryKeys.begin(), dm->boundaryKeys.end(), node->key | bit | tail);
			int bin = ((lowLocation - dm->boundaryKeys.begin() - 1) + (hiLocation - dm->boundaryKeys.begin() - 1)) / 2;
			if(dm->responsibleIndex[bin] == thisIndex) {
				if(((node->key | bit | tail) < rightSplitter) || (myPlace == numTreePieces - 1))
					node->setType(Internal);
				else {
					TreeNode* child = node->createRightChild();
					child->setType(NonLocal);
					nodeLookup[child->lookupKey()] = child;
					child->chareID = dm->responsibleIndex[bin + 1];
				}
			} else {
				TreeNode* child = node->createRightChild();
				child->setType(NonLocal);
				nodeLookup[child->lookupKey()] = child;
				child->chareID = dm->responsibleIndex[bin];
			}
		} else
			node->setType(Internal);
	}
	
}
*/

std::string keyBits(const Key k, const int numBits) {
	std::ostringstream oss;
	oss << "N";
	for(int i = 0; i < numBits; i++)
		oss << (k & (static_cast<Key>(1) << (62 - i)) ? 1 : 0);
	return oss.str();
}

/** Check that all the particles in the tree are really in their boxes.
 Because the keys are made of only the first 21 out of 23 bits of the
 floating point representation, there can be particles that are outside
 their box by tiny amounts.  Whether this is bad is not yet known. */
void checkTree(const TreeNode* node) {
	if(node->isBucket()) {
		for(FullParticle* iter = node->beginBucket; iter != node->endBucket; ++iter) {
			if(!node->box.contains(iter->position))
				cerr << "Not in the box: Box: " << node->box << " Position: " << iter->position << "\nNode key: " << keyBits(node->key, node->level) << "\nParticle key: " << keyBits(iter->key, 63) << endl;
		}
	} else if(!node->isNonLocal()) {
		if(node->leftChild)
			checkTree(node->leftChild);
		if(node->rightChild)
			checkTree(node->rightChild);
	}
}

/// Make a label for a node
string makeLabel(TreeNode* node) {
	ostringstream oss;
	oss << keyBits(node->key, node->level) << "\\n";
	switch(node->getType()) {
		case Invalid:
			oss << "Invalid";
			break;
		case Bucket:
			oss << "Bucket: " << (node->endBucket - node->beginBucket) << " particles";
			break;
		case Internal:
			oss << "Internal";
			break;
		case NonLocal:
			oss << "NonLocal: Chare " << node->chareID;
			break;
		case Empty:
			oss << "Empty";
			break;
		case Boundary:
			oss << "Boundary";
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
void printTree(TreeNode* node, ostream& os) {
	if(node == 0)
		return;
	
	string nodeID = keyBits(node->key, node->level);
	os << "\t\"" << nodeID << "\" [label=\"" << makeLabel(node) << "\"]\n";
	
	if(node->isNonLocal() || node->isBucket())
		return;
	
	os << "\t\"" << nodeID << "\" -> \"" << nodeID << "0\";\n";
	os << "\t\"" << nodeID << "\" -> \"" << nodeID << "1\";\n";
	if(node->leftChild == 0)
		os << "\t\"" << nodeID << "0\" [label=\"None\"]\n";
	else
		printTree(node->leftChild, os);
	if(node->rightChild == 0)
		os << "\t\"" << nodeID << "1\" [label=\"None\"]\n";
	else
		printTree(node->rightChild, os);
}
/*
void printBoxes(const TreeNode* root, const int level, ostream& os) {
	static int boxNum = 1;
	
	if(root == NULL)
		return;
	
	if(root->whollyOwned)
		os << root->box.tipsyCommand(boxNum++) << endl;
	else {
		if(root->chareID >= 0) {
			if(root->leftChild == 0) {
				TreeNode tn(const_cast<TreeNode *>(root), root->key);
				os << tn.box.tipsyCommand(boxNum++) << endl;
			} else
				printBoxes(root->leftChild, level + 1, os);
		}
		if(root->chareID <= 0) {
			if(root->rightChild == 0) {
				TreeNode tn(const_cast<TreeNode *>(root), root->key | (static_cast<Key>(1) << (62 - level)));
				os << tn.box.tipsyCommand(boxNum++) << endl;
			} else
				printBoxes(root->rightChild, level + 1, os);

		}
	}
}
*/

/// Write a file containing a graphviz dot graph of my tree
void TreePiece::report(const CkCallback& cb) {
	ostringstream outfilename;
	outfilename << "tree_" << thisIndex << ".dot";
	ofstream os(outfilename.str().c_str());

	os << "digraph G" << thisIndex << " {\n";
	os << "\tcenter = \"true\"\n";
	//os << "\tsize = \"8.5,11\"\n";
	os << "\tlabel = \"Piece: " << thisIndex << "\\nParticles: " 
			<< myNumParticles << "\\nLeft Splitter: " << keyBits(leftSplitter, 63)
			<< "\\nLeftmost Key: " << *(leftBoundary + 1) 
			<< "\\nRightmost Key: " << *(rightBoundary - 1) 
			<< "\\nRight Splitter: " << keyBits(rightSplitter, 63) << "\";\n";
	printTree(root, os);
	os << "}" << endl;
	
	os.close();
	
	checkTree(root);
	
	//also write a tipsy macro file showing the boxes I own
	/*
	outfilename.str("");
	outfilename << "tree_" << thisIndex << ".macro";
	os.open(outfilename.str().c_str());
	os << "tree" << thisIndex<< "\n";
	printBoxes(root, 0, os);
	os << "end" << endl;
	os.close();
	*/
	
	contribute(0, 0, CkReduction::concat, cb);
}
/*
class ImageRequest {
public:
	int width, height;
	Vector3D<double> x, y, corner;
	double lowValue, highValue;
	
	ImageRequest(const liveVizRequest3d& r) {
		width = r.wid;
		height = r.ht;
		x.x = r.x.x;
		x.y = r.x.y;
		x.z = r.x.z;
		y.x = r.y.x;
		y.y = r.y.y;
		y.z = r.y.z;
		//z is ignored currently
		//z.x = r.z.x;
		//z.y = r.z.y;
		//z.z = r.z.z;
		corner.x = r.o.x;
		corner.y = r.o.y;
		corner.z = r.o.z;
		lowValue = r.minZ;
		highValue = r.maxZ;
		
		corner += y - x;
		x *= 0.5 * width / x.lengthSquared();
		y *= -0.5 * height / y.lengthSquared();
	}
};

void TreePiece::generateImage(liveVizRequestMsg* m) {
	ImageRequest r(m->req);
		
	vector<byte> image(r.width * r.height, 0);
	int x, y;
	byte color;
	for(int i = 1; i < myNumParticles + 1; ++i) {
		x = static_cast<int>(dot(myParticles[i].position - r.corner, r.x));
		if(x < 0 || x >= r.width)
			continue;
		y = static_cast<int>(dot(myParticles[i].position - r.corner, r.y));
		if(y < 0 || y >= r.height)
			continue;
		if(myParticles[i].colorValue <= r.lowValue)
			color = 0;
		else if(myParticles[i].colorValue >= r.highValue)
			color = 254;
		else
			color = static_cast<byte>(255 * (myParticles[i].colorValue - r.lowValue) / (r.highValue - r.lowValue));
		if(image[x + y * r.width] < color)
			image[x + y * r.width] = color;
	}
	liveVizDeposit(m, 0, 0, r.width, r.height, image.begin(), this);
}

void TreePiece::applyPerParticleFunction(const CkCallback& cb) {
	perParticleFunction f = dm->userPerParticleFunction;
	for(int i = 1; i < myNumParticles + 1; ++i)
		f(myParticles[i]);
	contribute(0, 0, CkReduction::concat, cb);
}
*/
