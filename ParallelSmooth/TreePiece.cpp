/** \file TreePiece.cpp
 */

#include <sstream>
#include <algorithm>

#include "pup_toNetwork4.h"

#include "TreePiece.h"

#include "Tree.h"
#include "DataManager.h"

using std::ostringstream;
using std::vector;
using std::map;
using std::ofstream;

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

/// After the bounding box has been found, we can assign keys to the particles
void TreePiece::assignKeys(CkReductionMsg* m) {
	if(m->getSize() != sizeof(OrientedBox<float>)) {
		cerr << thisIndex << ": TreePiece: Fatal: Wrong size reduction message received!" << endl;
		callback.send(0);
		delete m;
		return;
	}
	
	// XXX: We may want to round off the bounding box here!
	//boundingBox = *static_cast<OrientedBox<float> *>(m->getData());
	boundingBox = OrientedBox<float>(Vector(-0.5, -0.5, -0.5), Vector(0.5, 0.5, 0.5));
	delete m;
	//give particles keys, using bounding box to scale
	for(vector<FullParticle>::iterator iter = myParticles.begin(); iter != myParticles.end(); ++iter)
		iter->generateKey(boundingBox);
	sort(myParticles.begin(), myParticles.end());
	
	contribute(0, 0, CkReduction::concat, callback);
	
	if(verbosity >= 5)
		cout << thisIndex << ": TreePiece: Assigned keys to all my particles" << endl;
}

/// Determine my part of the sorting histograms by counting the number of my particles in each bin
void TreePiece::evaluateBoundaries(const CkCallback& cb) {
	int numBins = dm->boundaryKeys.size() - 1;
	//this array will contain the number of particles I own in each bin
	vector<int> myBinCounts(numBins, 0);
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
	contribute(numBins * sizeof(int), myBinCounts.begin(), CkReduction::sum_int, cb);
}

/// Once final splitter keys have been decided, I need to give my particles out to the TreePiece responsible for them
void TreePiece::unshuffleParticles(const CkCallback& cb) {
	callback = cb;
	
	//find my responsibility
	myPlace = find(dm->responsibleIndex.begin(), dm->responsibleIndex.end(), thisIndex) - dm->responsibleIndex.begin();
	numTreePieces = dm->responsibleIndex.size();
	//assign my bounding keys
	leftSplitter = dm->boundaryKeys[myPlace];
	rightSplitter = dm->boundaryKeys[myPlace + 1];
	
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
		if((binEnd - binBegin) > 0)
			treePieces[*responsibleIter].acceptSortedParticles(binBegin, binEnd - binBegin);
		if(myParticles.end() <= binEnd)
			break;
		binBegin = binEnd;
	}
	//resize myParticles so we can fit the sorted particles and the boundary particles
	myParticles.resize(dm->particleCounts[myPlace] + 2);
}

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
	if(myPlace == 0) {
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
		leftBoundary = myParticles.begin();
		rightBoundary = myParticles.end() - 1;
		//I don't own my left and right boundaries, but I need to know them
		numParticles = rightBoundary - leftBoundary - 2 + 1;
		if(numParticles != dm->particleCounts[myPlace])
			cerr << thisIndex << ": TreePiece: Screwed up particle counts!" << endl;
		contribute(0, 0, CkReduction::concat, callback);
	}
}

/// I've got my sorted particles, now start building my tree
void TreePiece::startTreeBuild(const CkCallback& cb) {
	//I don't need this copy of my particles anymore
	mySortedParticles.clear();
	
	//create a root TreeNode
	root = new TreeNode(NULL, firstPossibleKey);
	//store this node in my lookup table
	nodeLookup[root->lookupKey()] = root;
	//build it
	buildTree(root, leftBoundary, rightBoundary, 0);
	
	//signify completion with a reduction
	contribute(0, 0, CkReduction::concat, cb);
}

/** A recursive algorithm for building my tree.
 \note This function is way too complicated to explain in words alone.
 */
void TreePiece::buildTree(TreeNode* parent, FullParticle* leftParticle, FullParticle* rightParticle, int level) {
	
	//check if we should bucket these particles
	if(rightParticle - leftParticle < TreeNode::maxBucketSize) {
		if((leftParticle != leftBoundary) && (rightParticle != rightBoundary)) {
			//cerr << "Bucketing: " << leftParticle << " : " << rightParticle << endl;
			//cerr << "Bucket: " << ((rightParticle - leftParticle) ) << endl;
			parent->key |= TreeNode::bucketFlag;
			parent->beginBucket = leftParticle;
			parent->endBucket = rightParticle + 1;
			return;
		}
	}
	
	//this is the bit we are looking at
	Key currentBitMask = static_cast<Key>(1) << (62 - level);
	//we need to know the values at the boundary
	Key leftBit = leftParticle->key & currentBitMask;
	Key rightBit = rightParticle->key & currentBitMask;
	
	if(leftBit ^ rightBit) { //a split at this level
		if(leftBit > rightBit)
			cerr << thisIndex << ": How the hell did this happen?" << endl;
		//find the index where the bit changes (currentKey is a node, we search through particles, so need to use particle bit)
		FullParticle* splitParticle = lower_bound(leftParticle, rightParticle + 1, parent->key | currentBitMask);// | TreeNode::bucketFlag);
		if(splitParticle - leftBoundary > 1) {
			parent->leftChild = new TreeNode(parent, parent->key);
			nodeLookup[parent->leftChild->lookupKey()] = parent->leftChild;
			parent->leftChild->parent = parent;
			buildTree(parent->leftChild, leftParticle, splitParticle - 1, level + 1);
		} else { //left child is non-local
			//find appropriate chare Id
			Key tail = (static_cast<Key>(1) << (62 - level)) - static_cast<Key>(1);
			vector<Key>::iterator lowLocation = upper_bound(dm->boundaryKeys.begin(), dm->boundaryKeys.end(), parent->key & ~tail);
			vector<Key>::iterator hiLocation = upper_bound(dm->boundaryKeys.begin(), dm->boundaryKeys.end(), parent->key | tail);
			int bin = ((lowLocation - dm->boundaryKeys.begin() - 1) + (hiLocation - dm->boundaryKeys.begin() - 1)) / 2;
			if(dm->responsibleIndex[bin] == thisIndex)
				parent->leftChild = 0;
			else
				parent->chareID = -1 - dm->responsibleIndex[bin];
		}
		if(rightBoundary > splitParticle) {
			parent->rightChild = new TreeNode(parent, parent->key | currentBitMask);
			nodeLookup[parent->rightChild->lookupKey()] = parent->rightChild;
			parent->rightChild->parent = parent;
			buildTree(parent->rightChild, splitParticle, rightParticle, level + 1);
		} else { //right child is non-local
			//find appropriate chare Id
			Key tail = (static_cast<Key>(1) << (62 - level)) - static_cast<Key>(1);
			Key bit = (static_cast<Key>(1) << (62 - level));
			vector<Key>::iterator lowLocation = upper_bound(dm->boundaryKeys.begin(), dm->boundaryKeys.end(), (parent->key | bit) & ~tail);
			vector<Key>::iterator hiLocation = upper_bound(dm->boundaryKeys.begin(), dm->boundaryKeys.end(), parent->key | bit | tail);
			int bin = ((lowLocation - dm->boundaryKeys.begin() - 1) + (hiLocation - dm->boundaryKeys.begin() - 1)) / 2;
			if(dm->responsibleIndex[bin] == thisIndex) {
				if(((parent->key | bit | tail) < rightSplitter) || (myPlace == numTreePieces - 1))
					parent->rightChild = 0;
				else
					parent->chareID = 1 + dm->responsibleIndex[bin + 1];
			} else
				parent->chareID = 1 + dm->responsibleIndex[bin];
		}
	
	} else if(leftBit & rightBit) { //both ones, make a right child
		parent->rightChild = new TreeNode(parent, parent->key | currentBitMask);
		nodeLookup[parent->rightChild->lookupKey()] = parent->rightChild;
		parent->rightChild->parent = parent;
		buildTree(parent->rightChild, leftParticle, rightParticle, level + 1);
		if((parent->key) <= leftBoundary->key) {
			//find appropriate chare Id
			Key tail = (static_cast<Key>(1) << (62 - level)) - static_cast<Key>(1);
			vector<Key>::iterator lowLocation = upper_bound(dm->boundaryKeys.begin(), dm->boundaryKeys.end(), parent->key & ~tail);
			vector<Key>::iterator hiLocation = upper_bound(dm->boundaryKeys.begin(), dm->boundaryKeys.end(), parent->key | tail);
			int bin = ((lowLocation - dm->boundaryKeys.begin() - 1) + (hiLocation - dm->boundaryKeys.begin() - 1)) / 2;
			if(dm->responsibleIndex[bin] == thisIndex)
				parent->leftChild = 0;
			else
				parent->chareID = -1 - dm->responsibleIndex[bin];
		}
	} else { //both zeros, make a left child
		parent->leftChild = new TreeNode(parent, parent->key);
		nodeLookup[parent->leftChild->lookupKey()] = parent->leftChild;
		parent->leftChild->parent = parent;
		buildTree(parent->leftChild, leftParticle, rightParticle, level + 1);
		currentBitMask >>= 1;
		if(((parent->key | currentBitMask) >> (62 - level)) >= (rightBoundary->key >> (62 - level))) {
			//find appropriate chare Id
			Key tail = (static_cast<Key>(1) << (62 - level)) - static_cast<Key>(1);
			Key bit = (static_cast<Key>(1) << (62 - level));
			vector<Key>::iterator lowLocation = upper_bound(dm->boundaryKeys.begin(), dm->boundaryKeys.end(), (parent->key | bit) & ~tail);
			vector<Key>::iterator hiLocation = upper_bound(dm->boundaryKeys.begin(), dm->boundaryKeys.end(), parent->key | bit | tail);
			int bin = ((lowLocation - dm->boundaryKeys.begin() - 1) + (hiLocation - dm->boundaryKeys.begin() - 1)) / 2;
			if(dm->responsibleIndex[bin] == thisIndex) {
				if(((parent->key | bit | tail) < rightSplitter) || (myPlace == numTreePieces - 1))
					parent->rightChild = 0;
				else
					parent->chareID = 1 + dm->responsibleIndex[bin + 1];
			} else
				parent->chareID = 1 + dm->responsibleIndex[bin];
		}
	}
	
}

string keyBits(const Key k, const int numBits) {
	ostringstream oss;
	for(int i = 0; i < numBits; i++)
		oss << (k & (static_cast<Key>(1) << (62 - i)) ? 1 : 0);
	return oss.str();
}

void printTree(TreeNode* root, const int level, ostream& os, int index) {
	if(root == 0)
		return;
	
	ostringstream oss;
	oss << "\"N" << keyBits(root->key, level);
	string idPrefix = oss.str();
	oss << "-" << index << "\"";
	string idRoot = oss.str();
	oss.str("");
	oss << "0-" << index << "\"";
	string idLeft = idPrefix + oss.str();
	oss.str("");
	oss << "1-" << index << "\"";
	string idRight = idPrefix + oss.str();
	
	if(root->isBucket()) {
		os << "\t" << idRoot << " [label=\"Bucket: " << (root->endBucket - root->beginBucket) << "\"]\n";
		root->whollyOwned = true;
	} else if(root->chareID != 0) {
		if(root->chareID > 0) {
			os << "\t" << idRoot << " -> " << idLeft << ";\n";
			printTree(root->leftChild, level + 1, os, index);
			os << "\t" << idRoot << " -> " << idPrefix << "1-" << (root->chareID - 1) << "\" [color = \"red\"];\n";
		} else {
			os << "\t" << idRoot << " -> " << idPrefix << "0-" << (-root->chareID - 1) << "\" [color = \"red\"];\n";
			os << "\t" << idRoot << " -> " << idRight << ";\n";
			printTree(root->rightChild, level + 1, os, index);
		}
	} else {
		os << "\t" << idRoot << " -> " << idLeft << ";\n";
		if(root->leftChild != 0)
			printTree(root->leftChild, level + 1, os, index);
		else
			os << "\t" << idLeft << " [label=\"None\"];\n";

		os << "\t" << idRoot << " -> " << idRight << ";\n";
		if(root->rightChild != 0)
			printTree(root->rightChild, level + 1, os, index);
		else
			os << "\t" << idRight << " [label=\"None\"];\n";
		
		if(((root->leftChild == 0) || root->leftChild->whollyOwned) && ((root->rightChild == 0) || root->rightChild->whollyOwned))
			root->whollyOwned = true;
	}
}

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

/// Write a file containing a graphviz dot graph of my tree
void TreePiece::report() {
	ostringstream outfilename;
	outfilename << "tree_" << thisIndex << ".dot";
	ofstream os(outfilename.str().c_str());

	os << "digraph G" << thisIndex << " {\n";
	os << "\tcenter = \"true\"\n\tsize = \"8.5,11\"\n";
	os << "\tlabel = \"Piece: " << thisIndex << "\\nParticles: " 
			<< numParticles << "\\nLeft Splitter: " << keyBits(leftSplitter, 63)
			<< "\\nLeftmost Key: " << *(leftBoundary + 1) 
			<< "\\nRightmost Key: " << *(rightBoundary - 1) 
			<< "\\nRight Splitter: " << keyBits(rightSplitter, 63) << "\";\n";
	printTree(root, 0, os, thisIndex);
	os << "}" << endl;
	
	os.close();
	
	//also write a tipsy macro file showing the boxes I own
	outfilename.str("");
	outfilename << "tree_" << thisIndex << ".macro";
	os.open(outfilename.str().c_str());
	os << "tree" << thisIndex<< "\n";
	printBoxes(root, 0, os);
	os << "end" << endl;
	os.close();
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
	for(int i = 1; i < numParticles + 1; ++i) {
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
	for(int i = 1; i < numParticles + 1; ++i)
		f(myParticles[i]);
	contribute(0, 0, CkReduction::concat, cb);
}
*/
