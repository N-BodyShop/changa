/** @file TreePiece.cpp
 */

#include <cstdio>
#include <algorithm>
#include <fstream>
#include <assert.h>

//#include "ComlibManager.h"

#include "ParallelGravity.h"
#include "CacheManager.h"
#include "DataManager.h"
#include "Reductions.h"
#include "TipsyFile.h"

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
  // mark presence in the cache if we are in the first iteration (indicated by localCache==NULL)
  if (_cache && localCache==NULL) {
    localCache = cacheManagerProxy.ckLocalBranch();
    localCache->markPresence(thisIndex, root);
  }

  basefilename = fn;
  bLoaded = 0;

  //read in particles
  XDR xdrs;
  FILE* infile = fopen((basefilename + ".mass").c_str(), "rb");
  if(!infile) {
    if (thisIndex == 0) ckerr << "TreePiece " << thisIndex << ": Couldn't open masses file, reverting to tipsy..." << endl;
    contribute(0, 0, CkReduction::concat, cb);
    return;
  }
  xdrstdio_create(&xdrs, infile, XDR_DECODE);
	
  if(!xdr_template(&xdrs, &fh)) {
    ckerr << "TreePiece " << thisIndex << ": Couldn't read header from masses file, aborting" << endl;
    contribute(0, 0, CkReduction::concat, cb);
    return;
  }
	
  if(fh.magic != FieldHeader::MagicNumber || fh.dimensions != 1 || fh.code != float32) {
    ckerr << "TreePiece " << thisIndex << ": Masses file is corrupt or of incorrect type, aborting" << endl;
    contribute(0, 0, CkReduction::concat, cb);
    return;
  }

  unsigned int *startParticles;
  unsigned int *numParticlesChunk;
  unsigned int excess;
  switch (domainDecomposition) {
  case SFC_dec:
    numPrefetchReq = 2;
    prefetchReq = new OrientedBox<double>[2];
  case Oct_dec:
  case ORB_dec:
    if (numPrefetchReq == 0) {
      numPrefetchReq = 1;
      prefetchReq = new OrientedBox<double>[1];
    }
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
    //CkAbort("ORB domain decomposition not yet implemented");
    //break;
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
 
  //curBoundingBox = boundingBox;

  if(pos == maxPos) { //all the same position
    //XXX This would be bad!
    Key k;
    if(domainDecomposition!=ORB_dec){
      k = generateKey(pos, boundingBox);
    }
    for(u_int64_t i = 0; i < myNumParticles; ++i) {
      myParticles[i + 1].position = pos;
      if(domainDecomposition!=ORB_dec){
        myParticles[i + 1].key = k;
      }
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
      for(unsigned int i = 0; i < numParticlesChunk[chunkNum]; ++i) {
	      if(!xdr_template(&xdrs, &pos)) {
	        ckerr << "TreePiece " << thisIndex << ": Problem reading my part of the positions file, aborting" << endl;
	        CkAbort("Badness");
	      }
	      myParticles[++myPart].position = pos;
        if(domainDecomposition!=ORB_dec){
	        current = generateKey(pos, boundingBox);
	        myParticles[myPart].key = current;
        }
	      //CkPrintf("Adding key: %d = %16llx\n",myPart,current);
	      /*if (current < previous) {
	        CkPrintf("TreePiece %d: Key not ordered! (%016llx)\n",thisIndex,current);
	      }*/
	      previous = current;
      }
    }
    CkAssert(myPart == myNumParticles);
  }
	
  xdr_destroy(&xdrs);
  fclose(infile);
	
  if(verbosity > 3)
    ckerr << "TreePiece " << thisIndex << ": Read in masses and positions" << endl;
	
  bLoaded = 1;

  //printing all particles
  /*for(int i=0;i<myNumParticles;i++)
    CkPrintf("%lf,%lf,%lf\n",myParticles[i].position.x,myParticles[i].position.y,myParticles[i].position.z);*/
  if(domainDecomposition!=ORB_dec){
    sort(myParticles+1, myParticles+myNumParticles+1);
  }
  contribute(0, 0, CkReduction::concat, cb);
}

void TreePiece::loadTipsy(const std::string& filename, const CkCallback& cb) {
	callback = cb;
	
	bLoaded = 0;
	
	Tipsy::TipsyReader r(filename);
	if(!r.status()) {
		cerr << thisIndex << ": TreePiece: Fatal: Couldn't open tipsy file!" << endl;
		cb.send(0);
		return;
	}
	
	numTreePieces = numTreePieces;
	Tipsy::header h = r.getHeader();
	int totalNumParticles = h.nbodies;
	fh.numParticles = totalNumParticles;
	myNumParticles = totalNumParticles / numTreePieces;
	int excess = totalNumParticles % numTreePieces;
	unsigned int startParticle = myNumParticles * thisIndex;
	if(thisIndex < excess) {
		myNumParticles++;
		startParticle += thisIndex;
	} else
		startParticle += excess;
	
	if(verbosity > 2)
		cerr << thisIndex << ": TreePiece: Taking " << myNumParticles << " of " << totalNumParticles << " particles, starting at " << startParticle << endl;

	// allocate an array for myParticles
	myParticles = new GravityParticle[myNumParticles + 2];
	
	if(!r.seekParticleNum(startParticle)) {
		cerr << thisIndex << ": TreePiece: Fatal: Couldn't seek to my particles!" << endl;
		cb.send(0);
		return;
	}
	
	Tipsy::gas_particle gp;
	Tipsy::dark_particle dp;
	Tipsy::star_particle sp;
	for(unsigned int i = 0; i < myNumParticles; ++i) {
		if(i + startParticle < (unsigned int) h.nsph) {
			r.getNextGasParticle(gp);
			myParticles[i+1].mass = gp.mass;
			myParticles[i+1].position = gp.pos;
			myParticles[i+1].velocity = gp.vel;
			myParticles[i+1].soft = gp.hsmooth;
		} else if(i + startParticle < (unsigned int) h.nsph + h.ndark) {
			r.getNextDarkParticle(dp);
			myParticles[i+1].mass = dp.mass;
			myParticles[i+1].position = dp.pos;
			myParticles[i+1].velocity = dp.vel;
			myParticles[i+1].soft = dp.eps;
		} else {
			r.getNextStarParticle(sp);
			myParticles[i+1].mass = sp.mass;
			myParticles[i+1].position = sp.pos;
			myParticles[i+1].velocity = sp.vel;
			myParticles[i+1].soft = sp.eps;
		}
		myParticles[i+1].iOrder = i + startParticle;
#if COSMO_STATS > 1
		myParticles[i+1].intcellmass = 0;
		myParticles[i+1].intpartmass = 0;
		myParticles[i+1].extcellmass = 0;
		myParticles[i+1].extpartmass = 0;
#endif
#if COSMO_STATS > 0
		piecemass += myParticles[i+1].mass;
#endif
		boundingBox.grow(myParticles[i+1].position);
	}
	
	bLoaded = 1;
  contribute(sizeof(OrientedBox<float>), &boundingBox,
		   growOrientedBox_float,
		   CkCallback(CkIndex_TreePiece::assignKeys(0), pieces));
}

/// After the bounding box has been found, we can assign keys to the particles
void TreePiece::assignKeys(CkReductionMsg* m) {
	if(m->getSize() != sizeof(OrientedBox<float>)) {
		cerr << thisIndex << ": TreePiece: Fatal: Wrong size reduction message received!" << endl;
		callback.send(0);
		delete m;
		return;
	}
	
	boundingBox = *static_cast<OrientedBox<float> *>(m->getData());
	delete m;
  //curBoundingBox = boundingBox;
	if(thisIndex == 0 && verbosity)
		cerr << "TreePiece: Bounding box originally: "
		     << boundingBox << endl;
	//give particles keys, using bounding box to scale
  if(domainDecomposition!=ORB_dec){
	  for(unsigned int i = 0; i < myNumParticles; ++i) {
	    myParticles[i+1].key = generateKey(myParticles[i+1].position,
					       boundingBox);
	  }
  }
	
  if(domainDecomposition!=ORB_dec){
	  sort(&myParticles[1], &myParticles[myNumParticles+1]);
  }
  
	contribute(0, 0, CkReduction::concat, callback);
	
	if(verbosity >= 5)
		cout << thisIndex << ": TreePiece: Assigned keys to all my particles" << endl;
}

/**************ORB Decomposition***************/

/*void TreePiece::getBoundingBox(OrientedBox<float>& box){
  box = boundingBox;
}*/

bool comp_dim0(GravityParticle p1, GravityParticle p2) {
    return p1.position[0] < p2.position[0];
}
bool comp_dim1(GravityParticle p1, GravityParticle p2) {
    return p1.position[1] < p2.position[1];
}
bool comp_dim2(GravityParticle p1, GravityParticle p2) {
    return p1.position[2] < p2.position[2];
}

void TreePiece::initORBPieces(const CkCallback& cb){

  OrientedBox<float> box = boundingBox;
  orbBoundaries.clear();
  orbBoundaries.push_back(myParticles+1);
  orbBoundaries.push_back(myParticles+myNumParticles+1);
  firstTime=true;

  phase=0;

  compFuncPtr[0]= &comp_dim0;
  compFuncPtr[1]= &comp_dim1;
  compFuncPtr[2]= &comp_dim2;

	myBinCounts.clear();
	myBinCounts.push_back(myNumParticles);

  //Find out how many levels will be there before tree goes
  //into a treepiece
  chunkRootLevel=0;
  unsigned int tmp = numTreePieces;
  while(tmp){
    tmp >>= 1;
    chunkRootLevel++;
  }
  chunkRootLevel--;
  
  boxes = new OrientedBox<float>[chunkRootLevel+1];
  splitDims = new char[chunkRootLevel+1];

  boxes[0] = boundingBox;
  
  contribute(sizeof(OrientedBox<float>), &box, boxReduction, cb);
}

/*class Compare{ //Defines the comparison operator on the map used in balancer
  int dim;
public:
  Compare() {}
  Compare(int i) : dim(i) {}
  
  void setDim(int i){ dim = i; }
  
  bool operator()(GravityParticle& p1, GravityParticle& p2) const {
    return p1.position[dim] < p2.position[dim];
  }
};*/

void TreePiece::initBeforeORBSend(unsigned int myCount, const CkCallback& cb, const CkCallback& cback){

  callback = cb;
  //sorterCallBack = sorterCb;
  CkCallback nextCallback = cback;
  myExpectedCount = myCount;
  
  mySortedParticles.clear();
  mySortedParticles.reserve(myExpectedCount);
  
  /*if(myExpectedCount > myNumParticles){
    delete [] myParticles;
    myParticles = new GravityParticle[myExpectedCount + 2];
  }
  myNumParticles = myExpectedCount;*/

  contribute(0, 0, CkReduction::concat, nextCallback);
}

//void TreePiece::sendORBParticles(unsigned int myCount, const CkCallback& cb, const CkCallback& sorterCb){
void TreePiece::sendORBParticles(){

  /*callback = cb;
  sorterCallBack = sorterCb;
  myExpectedCount = myCount;

  std::list<GravityParticle *>::iterator iter;
  std::list<GravityParticle *>::iterator iter2;

  mySortedParticles.clear();
  mySortedParticles.reserve(myExpectedCount);*/
  
  std::list<GravityParticle *>::iterator iter;
  std::list<GravityParticle *>::iterator iter2;
  
  int i=0;
  for(iter=orbBoundaries.begin();iter!=orbBoundaries.end();iter++,i++){
		iter2=iter;
		iter2++;
		if(iter2==orbBoundaries.end())
			break;
    if(i==thisIndex){
			//CkPrintf("[%d] send %d particles to %d\n",thisIndex,myBinCounts[i],i);
      if(myBinCounts[i]>0)
        acceptORBParticles(*iter,myBinCounts[i]);
		}
    else{
			//CkPrintf("[%d] send %d particles to %d\n",thisIndex,myBinCounts[i],i);
      if(myBinCounts[i]>0)
        pieces[i].acceptORBParticles(*iter,myBinCounts[i]);
		}
  }

  if(myExpectedCount > myNumParticles){
    delete [] myParticles;
    myParticles = new GravityParticle[myExpectedCount + 2];
  }
  myNumParticles = myExpectedCount;
}

/// Accept particles from other TreePieces once the sorting has finished
void TreePiece::acceptORBParticles(const GravityParticle* particles, const int n) {
     
  copy(particles, particles + n, back_inserter(mySortedParticles));
  
  //CkPrintf("[%d] accepted %d particles:myexpected:%d,got:%d\n",thisIndex,n,myExpectedCount,mySortedParticles.size());
  if(myExpectedCount == mySortedParticles.size()) {
	  //I've got all my particles
    //Assigning keys to particles
    //Key k = 1 << 63;
    //Key k = thisIndex;
    for(int i=0;i<myExpectedCount;i++){
      mySortedParticles[i].key = thisIndex;
    }
	  //sort(mySortedParticles.begin(), mySortedParticles.end());
	  copy(mySortedParticles.begin(), mySortedParticles.end(), &myParticles[1]);
	  //signify completion with a reduction
	  if(verbosity>1)
      ckout << thisIndex <<" contributing to accept particles"<<endl;
	  if (root != NULL) {
	    root->fullyDelete();
	    delete root;
	    root = NULL;
      nodeLookupTable.clear();
	  }
	  contribute(0, 0, CkReduction::concat, callback);
  }
}

void TreePiece::finalizeBoundaries(ORBSplittersMsg *splittersMsg){
 
  CkCallback& cback = splittersMsg->cb;
  
  std::list<GravityParticle *>::iterator iter;
  std::list<GravityParticle *>::iterator iter2;

  iter = orbBoundaries.begin();
  iter2 = orbBoundaries.begin();
  iter2++;

  phase++;

  int index = thisIndex >> (chunkRootLevel-phase+1);

  Key lastBit;
  lastBit = thisIndex >> (chunkRootLevel-phase);
  lastBit = lastBit & 0x1;
  
  boxes[phase] = boxes[phase-1];
  if(lastBit){
    boxes[phase].lesser_corner[splittersMsg->dim[index]] = splittersMsg->pos[index];
  }
  else{
    boxes[phase].greater_corner[splittersMsg->dim[index]] = splittersMsg->pos[index];
  }

  splitDims[phase-1]=splittersMsg->dim[index];
  
  for(int i=0;i<splittersMsg->length;i++){
    
    int dimen=(int)splittersMsg->dim[i];
	  //Current location of the division is stored in a variable
    //Evaluate the number of particles in each division

    GravityParticle dummy;
    GravityParticle* divStart = *iter;
    Vector3D<double> divide(0.0,0.0,0.0);
    divide[dimen] = splittersMsg->pos[i];
    dummy.position = divide;
    GravityParticle* divEnd = upper_bound(*iter,*iter2,dummy,compFuncPtr[dimen]);
    
    orbBoundaries.insert(iter2,divEnd);
    iter = iter2;
    iter2++;
	}

  firstTime = true;
  
  myBinCounts.assign(2*splittersMsg->length,0);
  copy(tempBinCounts.begin(),tempBinCounts.end(),myBinCounts.begin());

  contribute(0,0,CkReduction::concat,cback);

}

void TreePiece::evaluateParticleCounts(ORBSplittersMsg *splittersMsg){

  //myBinCounts.assign(2*splittersMsg->length,0);
  
  //if(firstTime){
    //myBinCounts.assign(splittersMsg->length,0);
    //copy(tempBinCounts.begin(),tempBinCounts.end(),myBinCounts.begin());
  //}
  CkCallback& cback = splittersMsg->cb;
  
  tempBinCounts.assign(2*splittersMsg->length,0);

  std::list<GravityParticle *>::iterator iter;
  std::list<GravityParticle *>::iterator iter2;

  iter = orbBoundaries.begin();
  iter2 = orbBoundaries.begin();
  iter2++;
  
  for(int i=0;i<splittersMsg->length;i++){
    
    int dimen = (int)splittersMsg->dim[i];
    if(firstTime){
      sort(*iter,*iter2,compFuncPtr[dimen]);
    }
    //curDivision = pos;
    /*if(firstTime){
      curDivision = pos;
      phaseLeader = leader;
      Compare comp(dim);
      sort(myParticles+1, myParticles+myNumParticles+1,comp);
    }*/
	  //Current location of the division is stored in a variable
    //Evaluate the number of particles in each division

		GravityParticle dummy;
    GravityParticle* divStart = *iter;
    Vector3D<double> divide(0.0,0.0,0.0);
    divide[dimen] = splittersMsg->pos[i];
    dummy.position = divide;
    GravityParticle* divEnd = upper_bound(*iter,*iter2,dummy,compFuncPtr[dimen]);
    tempBinCounts[2*i] = divEnd - divStart;
    tempBinCounts[2*i + 1] = myBinCounts[i] - (divEnd - divStart);

    iter++; iter2++;
	}
  
  if(firstTime)
    firstTime=false;
  //thisProxy[phaseLeader].collectORBCounts(firstCnt,secondCnt);
  contribute(2*splittersMsg->length*sizeof(int), &(*tempBinCounts.begin()), CkReduction::sum_int, cback);
}
/**********************************************/

/*
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
*/

/// Determine my part of the sorting histograms by counting the number
/// of my particles in each bin
void TreePiece::evaluateBoundaries(const CkCallback& cb) {
  if (dm == NULL) {
    dm = (CProxy_DataManager(dataManagerID)).ckLocalBranch();
  }
	int numBins = dm->boundaryKeys.size() - 1;
	//this array will contain the number of particles I own in each bin
	myBinCounts.assign(numBins, 0);
	vector<Key>::const_iterator endKeys = dm->boundaryKeys.end();
	GravityParticle *binBegin = &myParticles[1];
	GravityParticle *binEnd;
	GravityParticle dummy;
	vector<int>::iterator binIter = myBinCounts.begin();
	vector<Key>::iterator keyIter = dm->boundaryKeys.begin();
	for(++keyIter; keyIter != endKeys; ++keyIter, ++binIter) {
		dummy.key = *keyIter;
		/// find the last place I could put this splitter key in
		/// my array of particles
		binEnd = upper_bound(binBegin, &myParticles[myNumParticles+1],
				     dummy);
		/// this tells me the number of particles between the
		/// last two splitter keys
		*binIter = (binEnd - binBegin);
		if(&myParticles[myNumParticles+1] <= binEnd)
			break;
		binBegin = binEnd;
	}
	
	//send my bin counts back in a reduction
	contribute(numBins * sizeof(int), &(*myBinCounts.begin()), CkReduction::sum_int, cb);
}

/// Once final splitter keys have been decided, I need to give my
/// particles out to the TreePiece responsible for them

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
	GravityParticle *binBegin = &myParticles[1];
	GravityParticle *binEnd;
	GravityParticle dummy;
	for(++iter; iter != endKeys; ++iter, ++responsibleIter) {
		dummy.key = *iter;
		//find particles between this and the last key
		binEnd = upper_bound(binBegin, &myParticles[myNumParticles+1],
				     dummy);
		//if I have any particles in this bin, send them to the responsible TreePiece
		if((binEnd - binBegin) > 0) {
        //CkPrintf("me:%d to:%d how many:%d\n",thisIndex,*responsibleIter,(binEnd-binBegin));
		    if(*responsibleIter == thisIndex)
			acceptSortedParticles(&(*binBegin), binEnd - binBegin);
		    else
			pieces[*responsibleIter].acceptSortedParticles(&(*binBegin), binEnd - binBegin);
		}
		if(&myParticles[myNumParticles + 1] <= binEnd)
			break;
		binBegin = binEnd;
	}
	//resize myParticles so we can fit the sorted particles and
	//the boundary particles
	if(dm->particleCounts[myPlace] > myNumParticles) {
	    delete[] myParticles;
	    myParticles = new GravityParticle[dm->particleCounts[myPlace] + 2];
	    }
	myNumParticles = dm->particleCounts[myPlace];
}

/// Accept particles from other TreePieces once the sorting has finished
void TreePiece::acceptSortedParticles(const GravityParticle* particles, const int n
) {
        copy(particles, particles + n, back_inserter(mySortedParticles));
        if(myPlace == -1)
	    myPlace = find(dm->responsibleIndex.begin(),
			   dm->responsibleIndex.end(), thisIndex)
		- dm->responsibleIndex.begin();
        if(dm->particleCounts[myPlace] == mySortedParticles.size()) {
	    //I've got all my particles
	    sort(mySortedParticles.begin(), mySortedParticles.end());
	    copy(mySortedParticles.begin(), mySortedParticles.end(),
		 &myParticles[1]);
	    //signify completion with a reduction
	    if(verbosity>1)
        ckout << thisIndex <<" contributing to accept particles"<<endl;
	    if (root != NULL) {
	      root->fullyDelete();
	      delete root;
	      root = NULL;
        nodeLookupTable.clear();
	    }
	    contribute(0, 0, CkReduction::concat, callback);
        }
}

// Sum energies for diagnostics
void TreePiece::calcEnergy(const CkCallback& cb) {
    double dEnergy[6]; // 0 -> kinetic; 1 -> virial ; 2 -> potential
    Vector3D<double> L;
        
    dEnergy[0] = 0.0;
    dEnergy[1] = 0.0;
    dEnergy[2] = 0.0;
    for(unsigned int i = 0; i < myNumParticles; ++i) {
	GravityParticle *p = &myParticles[i+1];
	
	dEnergy[0] += p->mass*p->velocity.lengthSquared();
	dEnergy[1] += p->mass*dot(p->treeAcceleration, p->position);
	dEnergy[2] += p->mass*p->potential;
	L += p->mass*cross(p->position, p->velocity);
	}
    dEnergy[0] *= 0.5;
    dEnergy[2] *= 0.5;
    dEnergy[3] = L.x;
    dEnergy[4] = L.y;
    dEnergy[5] = L.z;
    
    contribute(6*sizeof(double), dEnergy, CkReduction::sum_double, cb);
}

void TreePiece::kick(double dDelta, const CkCallback& cb) {
    for(unsigned int i = 0; i < myNumParticles; ++i)
	myParticles[i+1].velocity += dDelta*myParticles[i+1].treeAcceleration;
    
    contribute(0, 0, CkReduction::concat, cb);
}

void TreePiece::drift(double dDelta, const CkCallback& cb) {
  if (root != NULL) {
    // Delete the tree since is no longer useful
    root->fullyDelete();
    delete root;
    root = NULL;
    nodeLookupTable.clear();
  }

    for(unsigned int i = 0; i < myNumParticles; ++i)
	myParticles[i+1].position += dDelta*myParticles[i+1].velocity;
    
    contribute(0, 0, CkReduction::concat, cb);
}

void TreePiece::setSoft(const double dSoft) {
  for(unsigned int i = 1; i <= myNumParticles; ++i) {
      myParticles[i].soft = dSoft;
  }
}

void TreePiece::buildTree(int bucketSize, const CkCallback& cb) {
  maxBucketSize = bucketSize;
  callback = cb;

#if INTERLIST_VER > 0
  myTreeLevels=-1;
#endif
  //printing all particles
  //CkPrintf("\n\n\nbuilding tree, useTree:%d\n\n\n\n",useTree);
  //for(int i=0;i<myNumParticles+2;i++)
    //CkPrintf("[%d] %016llx  %lf,%lf,%lf\n",thisIndex,myParticles[i].key,myParticles[i].position.x,myParticles[i].position.y,myParticles[i].position.z);
  // decide which logic are we using to divide the particles: Oct or ORB
  switch (useTree) {
  case Binary_Oct:
  case Oct_Oct:
    Key bounds[2];
    //sort(myParticles+1, myParticles+myNumParticles+1);
#ifdef COSMO_PRINT
    CkPrintf("[%d] Keys: %016llx %016llx\n",thisIndex,myParticles[1].key,myParticles[myNumParticles].key);
#endif
    bounds[0] = myParticles[1].key;
    bounds[1] = myParticles[myNumParticles].key;
    contribute(2 * sizeof(Key), bounds, CkReduction::concat, CkCallback(CkIndex_TreePiece::collectSplitters(0), thisArrayID));
    break;
  case Binary_ORB:
    //CkAbort("ORB logic for tree-build not yet implemented");
    //contribute(0,0,CkReduction::concat,sorterCallBack);
    contribute(0, 0, CkReduction::concat, CkCallback(CkIndex_TreePiece::startORBTreeBuild(0), thisArrayID));
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

void TreePiece::quiescence() {
  char fout[100];
  sprintf(fout,"tree.%d.%d",thisIndex,iterationNo);
  ofstream ofs(fout);
  printTree(root,ofs);
  ofs.close();
  /*
  CkPrintf("[%d] quiescence, %d left\n",thisIndex,momentRequests.size());
  for (MomentRequestType::iterator iter = momentRequests.begin(); iter != momentRequests.end(); iter++) {
    CkVec<int> *l = iter->second;
    for (int i=0; i<l->length(); ++i) {
      CkPrintf("[%d] quiescence: %s to %d\n",thisIndex,keyBits(iter->first,63).c_str(),(*l)[i]);
    }
  }
  */
  CkPrintf("[%d] quiescence detected, pending %d\n",thisIndex,myNumParticlesPending);
  for (int i=0; i<numBuckets; ++i) {
    if (bucketReqs[i].numAdditionalRequests != 0)
      CkPrintf("[%d] requests for %d remaining %d\n",thisIndex,i,bucketReqs[i].numAdditionalRequests);
  }
  CkPrintf("quiescence detected!\n");
  mainChare.niceExit();
}

void TreePiece::collectSplitters(CkReductionMsg* m) {
  numSplitters = 2 * numTreePieces;
  delete[] splitters;
  splitters = new Key[numSplitters];
  Key* splits = static_cast<Key *>(m->getData());
  copy(splits, splits + numSplitters, splitters);
  KeyDouble* splitters2 = (KeyDouble *)splitters;
  //sort(splitters, splitters + numSplitters);
  sort(splitters2, splitters2 + numTreePieces);
  for (unsigned int i=1; i<numSplitters; ++i) {
    if (splitters[i] < splitters[i-1]) {
      for (int j=0; j<numSplitters; ++j) CkPrintf("%d: Key %d = %016llx\n",thisIndex,j,splitters[j]);
      if(thisIndex==0)
        CkAbort("Keys not ordered");
    }
  }
  splitters[0] = firstPossibleKey;
  contribute(0, 0, CkReduction::concat, CkCallback(CkIndex_TreePiece::startOctTreeBuild(0), thisArrayID));
  delete m;
  if(verbosity > 3)
    ckerr << "TreePiece " << thisIndex << ": Collected splitters" << endl;
}

/*****************ORB**********************/
/*void TreePiece::receiveBoundingBoxes(BoundingBoxes *msg){

  //OrientedBox<float> *boxesMsg = static_cast<OrientedBox<float>* >(msg->getData())
  //boxes = new OrientedBox<float>[numTreePieces];
  
  OrientedBox<float> *boxesMsg = msg->boxes;
  copy(boxesMsg,boxesMsg+numTreePieces,boxes);
 
  contribute(0, 0, CkReduction::concat, CkCallback(CkIndex_TreePiece::startORBTreeBuild(0), thisArrayID));
}*/

void TreePiece::startORBTreeBuild(CkReductionMsg* m){
  delete m;
 
  myParticles[0].key = thisIndex;
  myParticles[myNumParticles+1].key = thisIndex;

  /*//Find out how many levels will be there before tree goes
  //into a treepiece
  chunkRootLevel=0;
  unsigned int tmp = numTreePieces;
  while(tmp){
    tmp >>= 1;
    chunkRootLevel++;
  }
  chunkRootLevel--;*/
  compFuncPtr[0]= &comp_dim0;
  compFuncPtr[1]= &comp_dim1;
  compFuncPtr[2]= &comp_dim2;
  
  root = new BinaryTreeNode(1, numTreePieces>1?Tree::Boundary:Tree::Internal, 0, myNumParticles+1, 0);
  
  if (thisIndex == 0) root->firstParticle ++;
  if (thisIndex == (int)numTreePieces-1) root->lastParticle --;
  root->particleCount = myNumParticles;
  nodeLookupTable[(Tree::NodeKey)1] = root;

  //root->key = firstPossibleKey;
  root->boundingBox = boundingBox;
  //nodeLookup[root->lookupKey()] = root;
  numBuckets = 0;
  bucketList.clear();
  
  // Fix the root presence in the cache
  localCache->markPresence(thisIndex, root);

  if(verbosity > 3)
    ckerr << "TreePiece " << thisIndex << ": Starting tree build" << endl;

#if INTERLIST_VER > 0
  root->startBucket=0;
#endif
  // recursively build the tree
  buildORBTree(root, 0);

  delete [] boxes;
  delete [] splitDims;
  //Keys to all the particles have been assigned inside buildORBTree
  
  //CkPrintf("[%d] finished building local tree\n",thisIndex);
  
  // check all the pending requests in for RemoteMoments
  for (MomentRequestType::iterator iter = momentRequests.begin(); iter != momentRequests.end(); ) {
    NodeKey nodeKey = iter->first;
    GenericTreeNode *node = keyToNode(nodeKey);
    CkVec<int> *l = iter->second;
    //CkPrintf("[%d] checking moments requests for %s (%s) upon treebuild finished\n",thisIndex,keyBits(iter->first,63).c_str(),keyBits(node->getKey(),63).c_str());
    CkAssert(node != NULL);
    // we actually need to increment the iterator before deleting the element,
    // otherwise the iterator lose its validity!
    iter++;
    if (node->getType() == Empty || node->moments.totalMass > 0) {
      for (int i=0; i<l->length(); ++i) {
	streamingProxy[(*l)[i]].receiveRemoteMoments(nodeKey, node->getType(), node->firstParticle, node->particleCount, node->moments);
	//CkPrintf("[%d] sending moments of %s to %d upon treebuild finished\n",thisIndex,keyBits(node->getKey(),63).c_str(),(*l)[i]);
      }
      delete l;
      momentRequests.erase(node->getKey());
    }
  }

  if(verbosity > 3)
    ckerr << "TreePiece " << thisIndex << ": Number of buckets: " << numBuckets << endl;
  if(verbosity > 3)
    ckerr << "TreePiece " << thisIndex << ": Finished tree build, resolving boundary nodes" << endl;

  if (numTreePieces == 1) contribute(0, 0, CkReduction::concat, callback);

}

OrientedBox<float> TreePiece::constructBoundingBox(GenericTreeNode* node,int level, int numChild){

  OrientedBox<float> tmpBox;
  if(node->getType()==NonLocal){
    if(numChild==0){
      tmpBox = boxes[level];
      tmpBox.greater_corner[splitDims[level]] = boxes[level+1].lesser_corner[splitDims[level]];
    }
    else{
      tmpBox = boxes[level];
      tmpBox.lesser_corner[splitDims[level]] = boxes[level+1].greater_corner[splitDims[level]];
    }
    return tmpBox;
  }
  else{
    return boxes[level+1];
  }
  
  /*SFC::Key tmp = 1 << (level+1);
  tmp = key - tmp;
  
  for(int i=0;i<numTreePieces;i++){
    if(tmp==(i>>(chunkRootLevel-level-1))){
      tmpBox.grow(boxes[i].lesser_corner);
      tmpBox.grow(boxes[i].greater_corner);
    }
  }
  return tmpBox;*/
}

void TreePiece::buildORBTree(GenericTreeNode * node, int level){
  
#if INTERLIST_VER > 0
  if(level>myTreeLevels)
    myTreeLevels=level;
#endif
  //CkPrintf("[%d] in build ORB Tree, level:%d\n",thisIndex,level);
  if (level == 63) {
    ckerr << thisIndex << ": TreePiece(ORB): This piece of tree has exhausted all the bits in the keys.  Super double-plus ungood!" << endl;
    ckerr << "Left particle: " << (node->firstParticle) << " Right particle: " << (node->lastParticle) << endl;
    ckerr << "Left key : " << keyBits((myParticles[node->firstParticle]).key, 63).c_str() << endl;
    ckerr << "Right key: " << keyBits((myParticles[node->lastParticle]).key, 63).c_str() << endl;
    return;
  }

  CkAssert(node->getType() == Boundary || node->getType() == Internal);
  
  node->makeOrbChildren(myParticles, myNumParticles, level, chunkRootLevel, compFuncPtr);

  GenericTreeNode *child;
  for (unsigned int i=0; i<node->numChildren(); ++i) {
    child = node->getChildren(i);
    CkAssert(child != NULL);

    if(level<chunkRootLevel){
      child->boundingBox = constructBoundingBox(child,level,i);
    }

#if INTERLIST_VER > 0
    child->startBucket=numBuckets;
#endif
    nodeLookupTable[child->getKey()] = child;
    if (child->getType() == NonLocal) {
      // find a remote index for the node
      int first, last;
      bool isShared = nodeOwnership(child->getKey(), first, last);
      //CkPrintf("[%d] child Key:%lld, firstOwner:%d, lastOwner:%d\n",thisIndex,child->getKey(),first,last);
      CkAssert(!isShared);
      if (last < first) {
	      // the node is really empty because falling between two TreePieces
	      child->setType(Empty);
	      child->remoteIndex = thisIndex;
      } else {
	      child->remoteIndex = first + (thisIndex & (last-first));
	      // if we have a remote child, the node is a Boundary. Thus count that we
	      // have to receive one more message for the NonLocal node
	      node->remoteIndex --;
	      // request the remote chare to fill this node with the Moments
	      streamingProxy[child->remoteIndex].requestRemoteMoments(child->getKey(), thisIndex);
	      //CkPrintf("[%d] asking for moments of %s to %d\n",thisIndex,keyBits(child->getKey(),63).c_str(),child->remoteIndex);
      }
    } else if (child->getType() == Internal && child->lastParticle - child->firstParticle < maxBucketSize) {
      CkAssert(child->firstParticle != 0 && child->lastParticle != myNumParticles+1);
      child->remoteIndex = thisIndex;
      child->makeBucket(myParticles);
      bucketList.push_back(child);

      //Assign keys to all the particles inside the bucket
      int num = child->lastParticle - child->firstParticle + 1;
      int bits = 0;

      while(num > (1<<bits)){ bits++; }

      Key mask = 1 << (level+1);
      mask = ~mask;
      Key tmpKey = child->getKey() & mask;
      tmpKey = tmpKey << bits;

      for(int i=child->firstParticle;i<=child->lastParticle;i++){
        myParticles[i].key = tmpKey;
        tmpKey++;
      }

#if INTERLIST_VER > 0
      child->bucketListIndex=numBuckets;
      child->startBucket=numBuckets;
#endif
      numBuckets++;
      if (node->getType() != Boundary) node->moments += child->moments;
    } else if (child->getType() == Empty) {
      child->remoteIndex = thisIndex;
    } else {
      if (child->getType() == Internal) child->remoteIndex = thisIndex;
      // else the index is already 0
      buildORBTree(child, level+1);
      // if we have a Boundary child, we will have to compute it's multipole
      // before we can compute the multipole of the current node (and we'll do
      // it in receiveRemoteMoments)
      if (child->getType() == Boundary) node->remoteIndex --;
      if (node->getType() != Boundary) node->moments += child->moments;
    }
  }

  /* The old version collected Boundary nodes, the new version collects NonLocal nodes */

  if (node->getType() == Internal) {
    calculateRadiusFarthestCorner(node->moments, node->boundingBox);
  }

}
/******************************************/

void TreePiece::startOctTreeBuild(CkReductionMsg* m) {
  delete m;
	
  //printing all particles
  /*CkPrintf("\n\n\nbuilding tree\n\n\n\n");
  for(int i=0;i<myNumParticles;i++)
    CkPrintf("%lf,%lf,%lf\n",myParticles[i].position.x,myParticles[i].position.y,myParticles[i].position.z);*/

  if(thisIndex == 0)
    myParticles[0].key = firstPossibleKey;
  else
    myParticles[0].key = splitters[2 * thisIndex - 1];
	
  if(thisIndex == (int) numTreePieces - 1)
    myParticles[myNumParticles + 1].key = lastPossibleKey;
  else
    myParticles[myNumParticles + 1].key = splitters[2 * thisIndex + 2];
	
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
  
  // Fix the root presence in the cache
  localCache->markPresence(thisIndex, root);

  if(verbosity > 3)
    ckerr << "TreePiece " << thisIndex << ": Starting tree build" << endl;

  // set the number of chunks in which we will split the tree for remote computation
  /* OLD, moved to the CacheManager and startIteration
  numChunks = root->getNumChunks(_numChunks);
  assert(numChunks > 0);
  remainingChunk = new int[numChunks];
  root->getChunks(_numChunks, prefetchRoots;
  */

#if INTERLIST_VER > 0
  root->startBucket=0;
#endif
  // recursively build the tree
  buildOctTree(root, 0);

  //CkPrintf("[%d] finished building local tree\n",thisIndex);
  
  // check all the pending requests in for RemoteMoments
  for (MomentRequestType::iterator iter = momentRequests.begin(); iter != momentRequests.end(); ) {
    NodeKey nodeKey = iter->first;
    GenericTreeNode *node = keyToNode(nodeKey);
    CkVec<int> *l = iter->second;
    //CkPrintf("[%d] checking moments requests for %s (%s) upon treebuild finished\n",thisIndex,keyBits(iter->first,63).c_str(),keyBits(node->getKey(),63).c_str());
    CkAssert(node != NULL);
    // we actually need to increment the iterator before deleting the element,
    // otherwise the iterator lose its validity!
    iter++;
    if (node->getType() == Empty || node->moments.totalMass > 0) {
      for (int i=0; i<l->length(); ++i) {
	streamingProxy[(*l)[i]].receiveRemoteMoments(nodeKey, node->getType(), node->firstParticle, node->particleCount, node->moments);
	//CkPrintf("[%d] sending moments of %s to %d upon treebuild finished\n",thisIndex,keyBits(node->getKey(),63).c_str(),(*l)[i]);
      }
      delete l;
      momentRequests.erase(node->getKey());
    }
  }

  if(verbosity > 3)
    ckerr << "TreePiece " << thisIndex << ": Number of buckets: " << numBuckets << endl;
  if(verbosity > 3)
    ckerr << "TreePiece " << thisIndex << ": Finished tree build, resolving boundary nodes" << endl;

  if (numTreePieces == 1) contribute(0, 0, CkReduction::concat, callback);
}

/// Determine who are all the owners of this node
/// @return true if the caller is part of the owners, false otherwise
inline bool TreePiece::nodeOwnership(const Tree::NodeKey nkey, int &firstOwner, int &lastOwner) {
  
  if(useTree == Binary_ORB){ // Added for ORB Trees
    int keyLevel=0;
    Key tmpKey = Key(nkey);
    while(tmpKey > 1){
      tmpKey >>= 1;
      keyLevel++;
    }
    if(keyLevel >= chunkRootLevel){
      tmpKey = nkey >> (keyLevel-chunkRootLevel);
      tmpKey = tmpKey - (1 << chunkRootLevel);
      firstOwner = tmpKey;
      lastOwner = tmpKey;
    }
    else{
      tmpKey = nkey << (chunkRootLevel - keyLevel);
      tmpKey = tmpKey - (1 << chunkRootLevel);
      firstOwner = tmpKey;

      Key mask = (1 << (chunkRootLevel - keyLevel)) - 1;
      tmpKey = nkey << (chunkRootLevel - keyLevel);
      tmpKey = tmpKey - (1 << chunkRootLevel);
      tmpKey = tmpKey + mask;
      lastOwner = tmpKey;
    }
  }
  else{
    Key firstKey = Key(nkey);
    Key lastKey = Key(nkey + 1);
    const Key mask = Key(1) << 63;
    while (! (firstKey & mask)) {
      firstKey <<= 1;
      lastKey <<= 1;
    }
    firstKey &= ~mask;
    lastKey &= ~mask;
    lastKey -= 1;
    Key *locLeft = lower_bound(splitters, splitters + numSplitters, firstKey);
    Key *locRight = upper_bound(locLeft, splitters + numSplitters, lastKey);
    firstOwner = (locLeft - splitters) >> 1;
    lastOwner = (locRight - splitters - 1) >> 1;
#if COSMO_PRINT > 1
    std::string str = keyBits(nkey,63);
    CkPrintf("[%d] NO: key=%s, first=%d, last=%d\n",thisIndex,str.c_str(),locLeft-splitters,locRight-splitters);
#endif
  }
  return (thisIndex >= firstOwner && thisIndex <= lastOwner);
}

/** A recursive algorithm for building my tree.
    Examines successive bits in the particles' keys, looking for splits.
    Each bit is a level of nodes in the tree.  We keep going down until
    we can bucket the particles.
*/
void TreePiece::buildOctTree(GenericTreeNode * node, int level) {

#if INTERLIST_VER > 0
  if(level>myTreeLevels)
    myTreeLevels=level;
#endif
  
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
#if INTERLIST_VER > 0
    child->startBucket=numBuckets;
#endif
    nodeLookupTable[child->getKey()] = child;
    if (child->getType() == NonLocal) {
      // find a remote index for the node
      int first, last;
      bool isShared = nodeOwnership(child->getKey(), first, last);
      CkAssert(!isShared);
      if (last < first) {
	// the node is really empty because falling between two TreePieces
	child->setType(Empty);
	child->remoteIndex = thisIndex;
      } else {
	child->remoteIndex = first + (thisIndex & (last-first));
	// if we have a remote child, the node is a Boundary. Thus count that we
	// have to receive one more message for the NonLocal node
	node->remoteIndex --;
	// request the remote chare to fill this node with the Moments
	streamingProxy[child->remoteIndex].requestRemoteMoments(child->getKey(), thisIndex);
	//CkPrintf("[%d] asking for moments of %s to %d\n",thisIndex,keyBits(child->getKey(),63).c_str(),child->remoteIndex);
      }
    } else if (child->getType() == Internal && child->lastParticle - child->firstParticle < maxBucketSize) {
      CkAssert(child->firstParticle != 0 && child->lastParticle != myNumParticles+1);
      child->remoteIndex = thisIndex;
      child->makeBucket(myParticles);
      bucketList.push_back(child);
#if INTERLIST_VER > 0
      child->bucketListIndex=numBuckets;
      child->startBucket=numBuckets;
#endif
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

  /* The old version collected Boundary nodes, the new version collects NonLocal nodes */

  if (node->getType() == Internal) {
    calculateRadiusFarthestCorner(node->moments, node->boundingBox);
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

inline void partBucketForce(ExternalGravityParticle *part, GenericTreeNode *req, GravityParticle *particles) {
  Vector3D<double> r;
  double rsq;
  double twoh, a, b;

  for(int j = req->firstParticle; j <= req->lastParticle; ++j) {
    r = part->position - particles[j].position;
    rsq = r.lengthSquared();
    twoh = part->soft + particles[j].soft;
    if(rsq != 0) {
      SPLINE(rsq, twoh, a, b);
      particles[j].treeAcceleration += r * (b * part->mass);
      particles[j].potential -= part->mass * a;
    }
  }
}

/** @TOOD change the loop so that for one particle in req all the forces are
    compute (i.e outside loop over req, inside loop over nodes (the loop is in
    other functions) */
/*inline void nodeBucketForce(GenericTreeNode *node, BucketGravityRequest& req) {
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
      double qirx = m.xx*r.x + m.xy*r.y + m.xz*r.z;
      double qiry = m.xy*r.x + m.yy*r.y + m.yz*r.z;
      double qirz = m.xz*r.x + m.yz*r.y + m.zz*r.z;
      double qir = 0.5*(qirx*r.x + qiry*r.y + qirz*r.z);
      double tr = 0.5*(m.xx + m.yy + m.zz);
      double qir3 = b*m.totalMass + d*qir - c*tr;
      req.potentials[j] -= m.totalMass * a + c*qir - b*tr;
      req.accelerations[j][0] -= qir3*r[0] - c*qirx;
      req.accelerations[j][1] -= qir3*r[1] - c*qiry;
      req.accelerations[j][2] -= qir3*r[2] - c*qirz;
  */  
			/******************ADDED**********************/
		/*	//SFCTreeNode* reqnode = bucketList[req.identifier];

			//for(unsigned int i = reqnode->beginParticle; i < reqnode->endParticle; ++i){
			//myParticles[reqnode->beginParticle + j].acceleration[0] -= qir3*r[0] - c*qirx;
			//myParticles[reqnode->beginParticle + j].acceleration[1] -= qir3*r[1] - c*qiry;
			//myParticles[reqnode->beginParticle + j].acceleration[2] -= qir3*r[2] - c*qirz;
			//}

		}
  }
}*/

inline void nodeBucketForce(GenericTreeNode *node, GenericTreeNode *req, GravityParticle *particles) {
  Vector3D<double> r;
  double rsq;
  double twoh, a, b, c, d;
  MultipoleMoments m = node->moments;
    
  Vector3D<double> cm(m.cm);

	//SFCTreeNode* reqnode = bucketList[req.identifier];

  for(int j = req->firstParticle; j <= req->lastParticle; ++j) {
    r = particles[j].position - cm;
    rsq = r.lengthSquared();
    twoh = m.soft + particles[j].soft;
    if(rsq != 0) {
      double dir = 1.0/sqrt(rsq);
      SPLINEQ(dir, rsq, twoh, a, b, c, d);
      double qirx = m.xx*r.x + m.xy*r.y + m.xz*r.z;
      double qiry = m.xy*r.x + m.yy*r.y + m.yz*r.z;
      double qirz = m.xz*r.x + m.yz*r.y + m.zz*r.z;
      double qir = 0.5*(qirx*r.x + qiry*r.y + qirz*r.z);
      double tr = 0.5*(m.xx + m.yy + m.zz);
      double qir3 = b*m.totalMass + d*qir - c*tr;
      particles[j].potential -= m.totalMass * a + c*qir - b*tr;
      particles[j].treeAcceleration.x -= qir3*r.x - c*qirx;
      particles[j].treeAcceleration.y -= qir3*r.y - c*qiry;
      particles[j].treeAcceleration.z -= qir3*r.z - c*qirz;
    }
  }
}

/*******************
#if INTERLIST_VER > 2
// This means always disabled!

//The force calculation routine which calculates force between a particle in my current bucket and all the contents of the interaction list (both node and particle lists)
inline void listPartForce(CkVec<CkVec<GenericTreeNode *> >& cellList, CkVec<CkVec<TreePiece::PartInfo> >& particleList, GravityParticle particle, int levels) {
  Vector3D<double> r;
  double rsq;
  double twoh, a, b, c, d;
  MultipoleMoments m;
    
  //double r0,r1,r2;
  double acc0=0,acc1=0,acc2=0,pot=0;
  
  GenericTreeNode *node;
  
  int len=0;
  
  for(int k=0;k<=levels;k++){
    len=cellList[k].length();
    for(unsigned int i=0;i<len;i++){
      node = cellList[k][i];
    
      m = node->moments;
      Vector3D<double> cm(m.cm);
    
      r = particle.position - cm;
      //r0=r[0]; r1=r[1]; r2=r[2];
      rsq = r.lengthSquared();
      twoh = m.soft + particle.soft;
      if(rsq != 0) {
        double dir = 1.0/sqrt(rsq);
        SPLINEQ(dir, rsq, twoh, a, b, c, d);
        double qirx = m.xx*r.x + m.xy*r.y + m.xz*r.z;
        double qiry = m.xy*r.x + m.yy*r.y + m.yz*r.z;
        double qirz = m.xz*r.x + m.yz*r.y + m.zz*r.z;
        double qir = 0.5*(qirx*r.x + qiry*r.y + qirz*r.z);
        double tr = 0.5*(m.xx + m.yy + m.zz);
        double qir3 = b*m.totalMass + d*qir - c*tr;
        pot -= m.totalMass * a + c*qir - b*tr;
        acc0 -= qir3*r.x - c*qirx;
        acc1 -= qir3*r.y - c*qiry;
        acc2 -= qir3*r.z - c*qirz;
    
		  }
    }
  //}

  double e, f;

  
  //for(int k=0;k<=levels;k++){
    len=particleList[k].length();
    for(unsigned int i=0;i<len;i++){
    
      TreePiece::PartInfo pinfo = particleList[k][i];
      GravityParticle *part;
  	  for(int j = 0; j < pinfo.numParticles; ++j){
        part=&pinfo.particles[j];
        r = part->position - particle.position;
        rsq = r.lengthSquared();
        twoh = part->soft + particle.soft;
        if(rsq != 0) {
          SPLINE(rsq, twoh, e, f);
          acc0 += part->mass * r.x * f;
          acc1 += part->mass * r.y * f;
          acc2 += part->mass * r.z * f;
        
          //particle.acceleration += part->mass * r * f;
          pot -= part->mass * e;
        }
      }
    }
  }
  
  particle.treeAcceleration.x+=acc0;
  particle.treeAcceleration.y+=acc1;
  particle.treeAcceleration.z+=acc2;
  particle.potential+=pot;
}

#endif
******************/
template <typename T>
inline bool TreePiece::openCriterionBucket(GenericTreeNode *node,
					   OrientedBox<T> &boundingBox) {
#if COSMO_STATS > 0
  node->used = true;
#endif
  // Note that some of this could be pre-calculated into an "opening radius"
  Sphere<double> s(node->moments.cm,
		   opening_geometry_factor * node->moments.radius / theta);
  
  /*bool ret = Space::intersect(boundingBox, s);

  double rad = 0.866025403784*pow(node->boundingBox.volume(),0.333333);

  Sphere<double> s1(node->moments.cm,
		   opening_geometry_factor * rad / theta);

  bool ret1 = Space::intersect(boundingBox, s1);

  if(ret1!=ret){
    openingDiffCount++;
  }
  return ret;*/
  
  return Space::intersect(boundingBox, s);
}

inline int TreePiece::openCriterionNode(GenericTreeNode *node,
                                        GenericTreeNode *myNode) {
#if COSMO_STATS > 0
  node->used = true;
#endif
  // Note that some of this could be pre-calculated into an "opening radius"
  Sphere<double> s(node->moments.cm,
		   opening_geometry_factor * node->moments.radius / theta);

	if(myNode->getType()==Bucket || myNode->getType()==CachedBucket || myNode->getType()==NonLocalBucket){
  	//if(Space::intersect(myNode->boundingBox, s))
  	if(Space::intersect(myNode->boundingBox, s))
			return 1;
		else
			return 0;
	}
	else{
		if(Space::intersect(myNode->boundingBox, s)){
			if(Space::contained(myNode->boundingBox,s))
				return 1;
			else
				return -1;
		}
		else
			return 0;
	}
}

void TreePiece::initBuckets() {
  for (int j=0; j<numBuckets; ++j) {
    GenericTreeNode* node = bucketList[j];
    int numParticlesInBucket = node->particleCount;

    CkAssert(numParticlesInBucket <= maxBucketSize);
    //BucketGravityRequest req(numParticlesInBucket);
    //req.startingNode = root->getKey();
    //req.identifier = j;
    //req.numAdditionalRequests = numChunks;
    //req.requestingPieceIndex = thisIndex;
    for(int i = node->firstParticle; i <= node->lastParticle; ++i) {
      /*
      req.softs[i - node->firstParticle] = myParticles[i].soft;
      req.positions[i - node->firstParticle] = myParticles[i].position;
      req.potentials[i - node->firstParticle] = 0;
      req.boundingBox.grow(myParticles[i].position);
      */
      myParticles[i].treeAcceleration = 0;
      myParticles[i].potential = 0;
    }
    bucketReqs[j].finished = 0;
    bucketReqs[j].numAdditionalRequests = numChunks;

/*#if COSMO_DEBUG > 1
    if(iterationNo==1 || listMigrated==true){
      std::set<Tree::NodeKey> *list = new std::set<Tree::NodeKey>();
      bucketcheckList.push_back(*list);
    }
#endif*/
  }
#if COSMO_DEBUG > 1
  bucketcheckList.resize(numBuckets);
#endif
}

void TreePiece::startNextBucket() {
  if(currentBucket >= numBuckets)
    return;

  /*	now done in initBuckets
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
  */

  // start the tree walk from the tree built in the cache
#ifdef CACHE_TREE
  walkBucketTree(localCache->getRoot(), currentBucket);
#else
  walkBucketTree(root, currentBucket);
#endif
  bucketReqs[currentBucket].finished = 1;
  finishBucket(currentBucket);
  //	currentBucket++;
  //	startNextBucket();
}

void TreePiece::finishBucket(int iBucket) {
  BucketGravityRequest *req = &bucketReqs[iBucket];
  GenericTreeNode *node = bucketList[iBucket];
#ifdef COSMO_PRINT
  CkPrintf("[%d] Is finished %d? finished=%d, %d still missing!\n",thisIndex,iBucket,req->finished,req->numAdditionalRequests);
#endif

  if(req->finished && req->numAdditionalRequests == 0) {
    myNumParticlesPending -= node->particleCount;
#ifdef COSMO_PRINT
    CkPrintf("[%d] Finished bucket %d, %d particles remaining\n",thisIndex,iBucket,myNumParticlesPending);
#endif
    /*
    int iStart = bucketList[iBucket]->firstParticle;
    for(unsigned int i = 0; i < req->numParticlesInBucket; ++i) {
      myParticles[iStart + i].treeAcceleration
	+= req->accelerations[i];
      myParticles[iStart + i].potential
	+= req->potentials[i];
    }
    */
    if(started && myNumParticlesPending == 0) {
      started = false;
      contribute(0, 0, CkReduction::concat, callback);
      /*   cout << "TreePiece " << thisIndex << ": Made " << myNumProxyCalls
	   << " proxy calls forward, " << myNumProxyCallsBack
	   << " to respond in finishBucket" << endl;*/
#if COSMO_STATS > 0
      if(verbosity)
	CkPrintf("[%d] TreePiece %d finished with bucket %d , openCriterions:%lld, openingDiffCount:%lld\n",CkMyPe(),thisIndex,iBucket,numOpenCriterionCalls,openingDiffCount);
#else
      if(verbosity)
	CkPrintf("[%d] TreePiece %d finished with bucket %d\n",CkMyPe(),thisIndex,iBucket);
#endif
      if(verbosity > 4)
	ckerr << "TreePiece " << thisIndex << ": My particles are done"
	     << endl;
    }
  }
}

void TreePiece::doAllBuckets(){
#if COSMO_DEBUG > 0
  char fout[100];
  sprintf(fout,"tree.%d.%d",thisIndex,iterationNo);
  ofstream ofs(fout);
  printTree(root,ofs);
  ofs.close();
#endif

  dummyMsg *msg = new (8*sizeof(int)) dummyMsg;
  *((int *)CkPriorityPtr(msg)) = numTreePieces * numChunks + thisIndex + 1;
  CkSetQueueing(msg,CK_QUEUEING_IFIFO);
  msg->val=0;

#if INTERLIST_VER > 0
  checkListLocal[0].length()=0;
  for(int i=0;i<=myTreeLevels;i++){
    cellListLocal[i].length()=0;
    particleListLocal[i].length()=0;
  }
#endif
  
  thisProxy[thisIndex].nextBucket(msg);
}

void TreePiece::nextBucket(dummyMsg *msg){
  unsigned int i=0;

#if INTERLIST_VER > 0
  
  NodeType childType;
 
  while (i<_yieldPeriod && currentBucket < numBuckets) {

    if(currentBucket==0){
      curNodeLocal=root;
      curLevelLocal=0;
    }
    
    //Computes interaction list for this bucket
    preWalkInterTree();
    
    //CkAssert(curNodeLocal->getType()==Bucket);
    CkAssert(checkListLocal[curLevelLocal].length()==0);
    //Go over both the lists to calculate forces with one bucket
    
    if(myLocalCheckListEmpty && curNodeLocal->getType()!=Bucket){
      
      GenericTreeNode *tmpNode;
  
      int startBucket=curNodeLocal->startBucket;
      int tmpBucket;
      int lastBucket;
      int k;
  
      for(k=startBucket+1;k<numBuckets;k++){
        tmpNode = bucketList[k];
        if(tmpNode->lastParticle>curNodeLocal->lastParticle)
          break;
      }
      lastBucket=k-1;
    
      for(k=startBucket;k<=lastBucket;k++){
        calculateForceLocalBucket(k);
        currentBucket++;
        i++;
      }
      myLocalCheckListEmpty=false;
    }
    else{
      CkAssert(curNodeLocal->getType()==Bucket);
      calculateForceLocalBucket(curNodeLocal->bucketListIndex);
      currentBucket++;
      i++;
    }

    /*GenericTreeNode* node = bucketList[curNodeLocal->bucketListIndex];
  
    for(unsigned int j = node->firstParticle; j <= node->lastParticle; ++j) {
      //for(int k=0;k<=curLevelLocal;k++){
        //if(cellListLocal[k].length()>0 || particleListLocal[k].length()>0)
          listPartForce(cellListLocal,particleListLocal,myParticles[j],curLevelLocal);
      //}
    }*/
   
    //Right now, because of -O3 BucketForce routines perform better than partForce routine
    //So, using the following
    if(currentBucket>=numBuckets)
      break;
    
    //Calculate starting Node for next tree walk
    GenericTreeNode *tmpNode=curNodeLocal;
    int flag=0;
    tmpNode->visitedL=true;
    while(tmpNode->parent != NULL){
      
      cellListLocal[curLevelLocal].length()=0;
      particleListLocal[curLevelLocal].length()=0;
      
      tmpNode=tmpNode->parent;
      curLevelLocal--;
      GenericTreeNode* childIterator;
      for(unsigned int j = 0; j < tmpNode->numChildren(); ++j) {
        childIterator = tmpNode->getChildren(j);
        CkAssert (childIterator != NULL);
        childType = childIterator->getType();
        if(childIterator->visitedL == false){
          if(childType == NonLocal || childType == Cached || childType == NonLocalBucket || childType == CachedBucket || childType==Empty || childType==CachedEmpty){
            childIterator->visitedL=true;
          }
          else{
            curNodeLocal=childIterator;
            curLevelLocal++;
            flag=1;
            break;
          }
        }
      }
      if(flag==0){
        tmpNode->visitedL=true;
      }
      if(flag==1){
        flag=0;
        break;
      }
    }
      
	}

#else 
  while(i<_yieldPeriod && currentBucket<numBuckets){
    startNextBucket();
    currentBucket++;
    i++;
  }
#endif

  if (currentBucket<numBuckets) {
    thisProxy[thisIndex].nextBucket(msg);
  } else {
    delete msg;
  }
}

void TreePiece::calculateGravityLocal() {
  doAllBuckets();
}

#if INTERLIST_VER > 0
void TreePiece::initNodeStatus(GenericTreeNode *node){
      
  GenericTreeNode *child;
  node->visitedR=false;
 
  NodeType childType;
 
  if(node->getType()==Bucket)
    return;
  
  for(unsigned int j = 0; j < node->numChildren(); ++j) {
    child = node->getChildren(j);
    CkAssert (child != NULL);
    childType = child->getType();
    if(!(childType == NonLocal || childType == NonLocalBucket || childType == Cached || childType == CachedBucket || childType==Empty || childType==CachedEmpty)){
      initNodeStatus(child);
    }
  }
}

void TreePiece::calculateForceLocalBucket(int bucketIndex){
    int cellListIter;
    int partListIter;
    
    //BucketGravityRequest& req = bucketReqs[bucketIndex];
    for(int k=0;k<=curLevelLocal;k++){
      nodeInterLocal += cellListLocal[k].length()*bucketList[bucketIndex]->particleCount;
      for(cellListIter=0;cellListIter<cellListLocal[k].length();cellListIter++){
        GenericTreeNode *tmp = cellListLocal[k][cellListIter];
#if COSMO_DEBUG > 1
        bucketcheckList[bucketIndex].insert(tmp->getKey());
        combineKeys(tmp->getKey(),bucketIndex);
#endif
        nodeBucketForce(tmp, bucketList[bucketIndex], myParticles);
      }

      LocalPartInfo pinfo;
      for(partListIter=0;partListIter<particleListLocal[k].length();partListIter++){
        pinfo = particleListLocal[k][partListIter];
#if COSMO_DEBUG > 1
        bucketcheckList[bucketIndex].insert((pinfo.nd)->getKey());
        combineKeys((pinfo.nd)->getKey(),bucketIndex);
#endif
        particleInterLocal += bucketList[bucketIndex]->particleCount * (pinfo.numParticles);
        for(int j = 0; j < pinfo.numParticles; ++j){
          partBucketForce(&pinfo.particles[j], bucketList[bucketIndex], myParticles);
        }
      }
    }

    bucketReqs[bucketIndex].finished = 1;
    finishBucket(bucketIndex);

}

void TreePiece::calculateForceRemoteBucket(int bucketIndex, int chunk){
    
  int cellListIter;
  int partListIter;
  
  //BucketGravityRequest& req = bucketReqs[bucketIndex];
  for(int k=0;k<=curLevelRemote;k++){
    nodeInterRemote[chunk] += cellList[k].length()*bucketList[bucketIndex]->particleCount;
    for(cellListIter=0;cellListIter<cellList[k].length();cellListIter++){
      GenericTreeNode *tmp= cellList[k][cellListIter];
#if COSMO_DEBUG > 1
      bucketcheckList[bucketIndex].insert(tmp->getKey());
      combineKeys(tmp->getKey(),bucketIndex);
#endif
      nodeBucketForce(tmp, bucketList[bucketIndex], myParticles);
    }
    for(partListIter=0;partListIter<particleList[k].length();partListIter++){
      RemotePartInfo pinfo = particleList[k][partListIter];
#if COSMO_DEBUG > 1
      bucketcheckList[bucketIndex].insert((pinfo.nd)->getKey());
      combineKeys((pinfo.nd)->getKey(),bucketIndex);
#endif
      particleInterRemote[chunk] += bucketList[bucketIndex]->particleCount * (pinfo.numParticles);
      for(int j = 0; j < pinfo.numParticles; ++j){
        partBucketForce(&pinfo.particles[j], bucketList[bucketIndex], myParticles);
      }
    }
  }

  bucketReqs[bucketIndex].numAdditionalRequests--;
  finishBucket(bucketIndex);
  //remainingChunk[nChunk] -= bucketList[bucketIndex]->particleCount;

}
/*
void TreePiece::calculateForceRemoteBucket(int bucketIndex){
    
  GenericTreeNode* node = bucketList[bucketIndex];
  
  for(unsigned int j = node->firstParticle; j <= node->lastParticle; ++j) {
    //for(int k=0;k<=curLevelRemote;k++){
      //if(cellList[k].length()>0 || particleList[k].length()>0)
        listPartForce(cellList,particleList,myParticles[j],curLevelRemote);
    //}
  }
  bucketReqs[bucketIndex].numAdditionalRequests--;
  finishBucket(bucketIndex);

}
*/
#endif

void TreePiece::calculateGravityRemote(ComputeChunkMsg *msg) {
  unsigned int i=0;
  // cache internal tree: start directly asking the CacheManager
#ifdef CACHE_TREE
  GenericTreeNode *chunkRoot = localCache->chunkRootToNode(prefetchRoots[msg->chunkNum]);
#else
  GenericTreeNode *chunkRoot = keyToNode(prefetchRoots[msg->chunkNum]);
#endif
  if (chunkRoot == NULL) {
    chunkRoot = requestNode(thisIndex, prefetchRoots[msg->chunkNum], msg->chunkNum, -1, true);
  }
  CkAssert(chunkRoot != NULL);
#if COSMO_PRINT > 0
  CkPrintf("[%d] Computing gravity remote for chunk %d with node %016llx\n",thisIndex,msg->chunkNum,chunkRoot->getKey());
#endif

#if INTERLIST_VER > 0

  nChunk = msg->chunkNum;

  NodeType childType;
  //Init for each chunk
  //Possible race condition here....bfore calling calGravRemote...setting some stuff
  if(currentRemoteBucket==0){
    intCheckList[0].length()=0;
    extCheckList[0].length()=0;
    for(i=0;i<=myTreeLevels;i++){
      cellList[i].length()=0;
      particleList[i].length()=0;
    }
    i=0;

    if(nChunk>0)
      initNodeStatus(root);
  }
  
  while (i<_yieldPeriod && currentRemoteBucket < numBuckets) {

    if(currentRemoteBucket==0){
      curNodeRemote=root;
      curLevelRemote=0;
    }
    preWalkRemoteInterTree(chunkRoot,true);
    
    //Everything calculates forces with buckets
    
    CkAssert(intCheckList[curLevelRemote].length()==0 && extCheckList[curLevelRemote].length()==0);
    
    if(myCheckListEmpty && curNodeRemote->getType()!=Bucket){
	
      GenericTreeNode *tmpNode;
  
      int startBucket=curNodeRemote->startBucket;
      int tmpBucket;
      int lastBucket;
      int k;
  
      for(k=startBucket+1;k<numBuckets;k++){
        tmpNode = bucketList[k];
        if(tmpNode->lastParticle>curNodeRemote->lastParticle)
          break;
      }
      lastBucket=k-1;
    
      for(k=startBucket;k<=lastBucket;k++){
        calculateForceRemoteBucket(k, nChunk);
        remainingChunk[nChunk] -= bucketList[k]->particleCount;
        currentRemoteBucket++;
        i++;
      }
      myCheckListEmpty=false;
    }
    else{
      CkAssert(curNodeRemote->getType()==Bucket);
      calculateForceRemoteBucket(curNodeRemote->bucketListIndex, nChunk);
      
      remainingChunk[nChunk] -= bucketList[curNodeRemote->bucketListIndex]->particleCount;
      currentRemoteBucket++;
      i++;
    }

    //Go over both the lists to calculate forces with one bucket
    
    if(currentRemoteBucket>=numBuckets)
      break;
    
    //Calculate starting Node for next tree walk
    GenericTreeNode *tmpNode=curNodeRemote;
    int flag=0;
    tmpNode->visitedR=true;
    while(tmpNode->parent != NULL){
      
      cellList[curLevelRemote].length()=0;
      particleList[curLevelRemote].length()=0;
      
      tmpNode=tmpNode->parent;
      curLevelRemote--;
      GenericTreeNode* childIterator;
      for(unsigned int j = 0; j < tmpNode->numChildren(); ++j) {
        childIterator = tmpNode->getChildren(j);
        CkAssert (childIterator != NULL);
        childType=childIterator->getType();
        if(childIterator->visitedR == false){
          if(childType == NonLocal || childType == Cached || childType == NonLocalBucket || childType == CachedBucket || childType==Empty || childType==CachedEmpty){
            childIterator->visitedR=true;
          }
          else{
            curNodeRemote=childIterator;
            curLevelRemote++;
            flag=1;
            break;
          }
        }
      }
      if(flag==0){
        tmpNode->visitedR=true;
      }
      if(flag==1){
        flag=0;
        break;
      }
    }
      
  }

#else

  while (i<_yieldPeriod && currentRemoteBucket < numBuckets) {
    bucketReqs[currentRemoteBucket].numAdditionalRequests--;
    walkBucketRemoteTree(chunkRoot, msg->chunkNum, currentRemoteBucket, true);

    finishBucket(currentRemoteBucket);
    remainingChunk[msg->chunkNum] -= bucketList[currentRemoteBucket]->particleCount;
    currentRemoteBucket++;
    i++;
  }

#endif

  if (currentRemoteBucket < numBuckets) {
    thisProxy[thisIndex].calculateGravityRemote(msg);
#if COSMO_PRINT > 0
    CkPrintf("{%d} sending self-message chunk %d, prio %d\n",thisIndex,msg->chunkNum,*(int*)CkPriorityPtr(msg));
#endif
  } else {
    currentRemoteBucket = 0;
    
    CkAssert(remainingChunk[msg->chunkNum] >= 0);
    if (remainingChunk[msg->chunkNum] == 0) {
      // we finished completely using this chunk, so we acknowledge the cache
      // if this is not true it means we had some hard misses
#ifdef COSMO_PRINT
      CkPrintf("[%d] Finished chunk %d\n",thisIndex,msg->chunkNum);
#endif
      cacheManagerProxy[CkMyPe()].finishedChunk(msg->chunkNum, nodeInterRemote[msg->chunkNum]+particleInterRemote[msg->chunkNum]);
    }
#if COSMO_PRINT > 0
    CkPrintf("{%d} resetting message chunk %d, prio %d\n",thisIndex,msg->chunkNum,*(int*)CkPriorityPtr(msg));
#endif
    delete msg;
  }
}

#if INTERLIST_VER > 0

void TreePiece::preWalkRemoteInterTree(GenericTreeNode *chunkRoot, bool isRoot){

    //Start copying the checkList of previous level to the next level
    int level;
    GenericTreeNode *child;
    GenericTreeNode *node;
    int flag=0;
    NodeType childType;

    while(1){
      level=curLevelRemote-1;

#if INTERLIST_VER==1
      intPrevListIter=0;
      extPrevListIter=0;
      intCheckList[curLevelRemote].length()=0;
      extCheckList[curLevelRemote].length()=0;

      if(curNodeRemote!=root){
        if(intCheckList[level].length()!=0){
          node=intCheckList[level][0];
          intPrevListIter=1;
        }
        else
          node=NULL;
        if(extCheckList[level].length()!=0){
          //extPrevListIter=0;
          undecidedExtList.enq(extCheckList[level][0]);
          extPrevListIter=1;
        }
			}
      else{
        node=chunkRoot;
      }

      if(node!=NULL || undecidedExtList.isEmpty()==0){
        if(node!=NULL){
          if(node==chunkRoot)
            walkRemoteInterTreeVerI(node,true);
          else
            walkRemoteInterTreeVerI(node,false);
        }
        else
          walkRemoteInterTreeVerI(node,false);
      }
      else{
        myCheckListEmpty=true;
        curNodeRemote = curNodeRemote->parent;
        curLevelRemote--;
        break;
      }
#endif
#if INTERLIST_VER==2

      extPrevListIter=0;
      extCheckList[curLevelRemote].length()=0;

      if(curNodeRemote!=root){
        if(extCheckList[level].length()!=0){
          node=extCheckList[level][0];
          extPrevListIter=1;
        }
        else
          node=NULL;
      }
      else{
        node=chunkRoot;
      }

      if(node!=NULL){
        if(node==chunkRoot)
          walkRemoteInterTreeVerII(node,true);
        else
          walkRemoteInterTreeVerII(node,false);
      }
      else{
        myCheckListEmpty=true;
        curNodeRemote = curNodeRemote->parent;
        curLevelRemote--;
        break;
      }
#endif
      //Loop breaking condition
      if(curNodeRemote->getType()==Bucket)
        break;
      
      for(int i=0;i<curNodeRemote->numChildren();i++){
        child = curNodeRemote->getChildren(i);
      	CkAssert (child != NULL);
        childType = child->getType();
        if(child->visitedR==false){
          if(childType == NonLocal || childType == Cached || childType == NonLocalBucket || childType == CachedBucket || childType==Empty || childType==CachedEmpty){
            child->visitedR=true;
          }
          else{
            flag=1;
            break;
          }
        }
      }
      if(flag==1){
        curNodeRemote=child;
        curLevelRemote++;
        flag=0;
      }
      else{
        CkPrintf("Exceptional case\n");
        CkAssert(false);
        curNodeRemote->visitedR=true;
        break;
      }
    }

}

#endif

#if INTERLIST_VER == 1

void TreePiece::walkRemoteInterTreeVerI(GenericTreeNode *node, bool isRoot) {

  int chunk = nChunk;
  GenericTreeNode* myNode = curNodeRemote;
  int level = curLevelRemote;

  int flag=0;
 
  if(node==NULL){
    if(undecidedExtList.isEmpty())
      CkAssert(false);
    cachedWalkInterTreeVerI(undecidedExtList.deq());
    return;
  }
  
	NodeType nodeType = node->getType();
  //}
  // The order in which the following if are checked is very important for correctness
  if(nodeType == Bucket || nodeType == Internal || nodeType == Empty) {
  } else if(openCriterionNode(node, myNode)==0) {
#if COSMO_STATS > 0
  numOpenCriterionCalls++;
#endif
    if (isRoot && (nodeType == Cached || nodeType == CachedBucket || nodeType == CachedEmpty)) {
      undecidedExtList.enq(node);
    } else {
    }
  } else if (nodeType != Boundary) {
#if COSMO_STATS > 0
  numOpenCriterionCalls++;
#endif
#if COSMO_STATS > 0
  numOpenCriterionCalls++;
#endif
    // means it is some kind of NonLocal or Cached
      if(openCriterionNode(node,myNode)==1 || myNode->getType()==Bucket)
        undecidedExtList.enq(node);
      else
        intCheckList[level].push_back(node);
  } else if(openCriterionNode(node, myNode)==1 || myNode->getType()==Bucket){
    // here the node is Boundary
#if COSMO_STATS > 0
  numOpenCriterionCalls++;
#endif
#if COSMO_STATS > 0
  numOpenCriterionCalls++;
#endif
#if COSMO_STATS > 0
    nodesOpenedRemote++;
#endif
		//Boundary node removed from the interaction list...we dont calculate force with boundary
    GenericTreeNode* childIterator;
    for(unsigned int i = 0; i < node->numChildren(); ++i) {
      childIterator = node->getChildren(i);
      CkAssert (childIterator != NULL);
      undecidedIntList.enq(childIterator);
    }
  }
  else{
#if COSMO_STATS > 0
  numOpenCriterionCalls++;
#endif
#if COSMO_STATS > 0
  numOpenCriterionCalls++;
#endif
    intCheckList[level].push_back(node);
  }
	
  if(myNode!=root){
    if(intPrevListIter<intCheckList[level-1].length()){
      intPrevListIter++;
      walkRemoteInterTreeVerI(intCheckList[level-1][intPrevListIter-1],false);
    }
    else if(!undecidedIntList.isEmpty())
      walkRemoteInterTreeVerI(undecidedIntList.deq(),false);
    else{
      if(!undecidedExtList.isEmpty()){
        cachedWalkInterTreeVerI(undecidedExtList.deq());
      }
    }
  }
  else{
    if(!undecidedIntList.isEmpty())
      walkRemoteInterTreeVerI(undecidedIntList.deq(),false);
    else{
      if(!undecidedExtList.isEmpty()){
        cachedWalkInterTreeVerI(undecidedExtList.deq());
      }
    }
  }
}

#endif

#if INTERLIST_VER==2

void TreePiece::walkRemoteInterTreeVerII(GenericTreeNode *node, bool isRoot) {

  cachedWalkInterTreeVerII(node);
  
}

#endif

void TreePiece::walkBucketRemoteTree(GenericTreeNode *node, int chunk, int reqID, bool isRoot) {
  Vector3D<double> cm(node->moments.cm);
  Vector3D<double> r;
  Sphere<double> s(cm, opening_geometry_factor * node->moments.radius / theta);

  // The order in which the following if are checked is very important for correctness
  if(node->getType() == Bucket || node->getType() == Internal || node->getType() == Empty) {
#if COSMO_PRINT > 0
    CkPrintf("[%d] bucket %d: internal %llx\n",thisIndex,reqID,node->getKey());
#endif
    return;
  } else if(!openCriterionBucket(node, bucketList[reqID]->boundingBox)) {
#if COSMO_STATS > 0
    numOpenCriterionCalls++;
#endif
    if (isRoot && (node->getType() == Cached || node->getType() == CachedBucket || node->getType() == CachedEmpty)) {
#ifdef CACHE_TREE
      GenericTreeNode *nd = localCache->getRoot();
#else
      GenericTreeNode *nd = root;
#endif
      while (nd != NULL) {
	// if one of the ancestors of this node is not opened, then we don't
	// want to duplicate the computation
#if COSMO_STATS > 0
	numOpenCriterionCalls++;
#endif
	if (!openCriterionBucket(nd, bucketList[reqID]->boundingBox)) {
#if COSMO_PRINT > 0
	  CkPrintf("[%d] bucket %d: not opened %llx, found node %llx\n",thisIndex,reqID,node->getKey(),nd->getKey());
#endif
	  return;
	}
	int which = nd->whichChild(node->getKey());
	GenericTreeNode *ndp = nd;
	nd = nd->getChildren(which);
#if COSMO_PRINT > 0
	if (nd!=NULL) CkPrintf("[%d] search: got %llx from %llx, %d's child\n",thisIndex,nd->getKey(),ndp->getKey(),which);
#endif
      }
#if COSMO_PRINT > 0
      CkPrintf("[%d] bucket %d: not opened %llx, called cache\n",thisIndex,reqID,node->getKey());
#endif
      cachedWalkBucketTree(node, chunk, reqID);
    } else {
#if COSMO_PRINT > 0
      CkPrintf("[%d] bucket %d: not opened %llx\n",thisIndex,reqID,node->getKey());
#endif
      return;
    }
  } else if (node->getType() != Boundary) {
    // means it is some kind of NonLocal or Cached
#if COSMO_STATS > 0
    numOpenCriterionCalls++;
#endif
#if COSMO_PRINT > 0
      CkPrintf("[%d] bucket %d: calling cached walk with %llx\n",thisIndex,reqID,node->getKey());
#endif
    if (isRoot) {
#ifdef CACHE_TREE
      GenericTreeNode *nd = localCache->getRoot();
#else
      GenericTreeNode *nd = root;
#endif
      while (nd != NULL && nd->getKey() != node->getKey()) {
	// if one of the ancestors of this node is not opened, then we don't
	// want to duplicate the computation
#if COSMO_STATS > 0
	numOpenCriterionCalls++;
#endif
	if (!openCriterionBucket(nd, bucketList[reqID]->boundingBox)) {
#if COSMO_PRINT > 0
	  CkPrintf("[%d] bucket %d: not opened %llx, found node %llx\n",thisIndex,reqID,node->getKey(),nd->getKey());
#endif
	  return;
	}
	int which = nd->whichChild(node->getKey());
	GenericTreeNode *ndp = nd;
	nd = nd->getChildren(which);
#if COSMO_PRINT > 0
	if (nd!=NULL) CkPrintf("[%d] search: got %llx from %llx, %d's child\n",thisIndex,nd->getKey(),ndp->getKey(),which);
#endif
      }
    }
    cachedWalkBucketTree(node, chunk, reqID);
  } else {
    // here the node is Boundary
#if COSMO_STATS > 0
    nodesOpenedRemote++;
#endif
#if COSMO_STATS > 0
    numOpenCriterionCalls++;
#endif
    GenericTreeNode* childIterator;
    for(unsigned int i = 0; i < node->numChildren(); ++i) {
      childIterator = node->getChildren(i);
      CkAssert (childIterator != NULL);
      walkBucketRemoteTree(childIterator, chunk, reqID, false);
    }
  }
}


void TreePiece::startIteration(double t, int n, Tree::NodeKey *k, const CkCallback& cb) {
  callback = cb;
  theta = t;
  if (n != numChunks && remainingChunk != NULL) {
    // reallocate remaining chunk to the new size
    delete[] remainingChunk;
    remainingChunk = NULL;
    delete[] nodeInterRemote;
    delete[] particleInterRemote;
  }
  numChunks = n;
  prefetchRoots = k;
  if (remainingChunk == NULL) {
    remainingChunk = new int[numChunks];
    nodeInterRemote = new u_int64_t[numChunks];
    particleInterRemote = new u_int64_t[numChunks];
  }
#if COSMO_STATS > 0
  //myNumProxyCalls = 0;
  //myNumProxyCallsBack = 0;
  //myNumCellInteractions=myNumParticleInteractions=myNumMACChecks=0;
  //cachecellcount=0;
  nodesOpenedLocal = 0;
  nodesOpenedRemote = 0;
  numOpenCriterionCalls=0;
#endif
  nodeInterLocal = 0;
  for (int i=0; i<numChunks; ++i) {
    nodeInterRemote[i] = 0;
    particleInterRemote[i] = 0;
  }
  particleInterLocal = 0;
  iterationNo++;
  CkAssert(localCache != NULL);
  if(verbosity)
    CkPrintf("TreePiece %d: I have %d buckets\n",thisIndex,numBuckets);

  if (bucketReqs==NULL) bucketReqs = new BucketGravityRequest[numBuckets];
  
  currentBucket = 0;
  currentRemoteBucket = 0;
  myNumParticlesPending = myNumParticles;
  started = true;

#if INTERLIST_VER > 0
    myCheckListEmpty=false;
    curLevelLocal=0;
    curNodeLocal=NULL;
    curLevelRemote=0;
    curNodeRemote=NULL;
    nChunk=-1;
#endif

  initBuckets();

#if INTERLIST_VER > 0
  //Initialize all the interaction and check lists with empty lists
  myTreeLevels++;
  //myTreeLevels++;
  CkAssert(myTreeLevels>0);
  //if(listMigrated==true || iterationNo==1){
    /*for(int i=0;i<=myTreeLevels;i++){
      CkVec<GenericTreeNode *> *clist = new CkVec<GenericTreeNode *>();
      CkVec<PartInfo> *plist = new CkVec<PartInfo>();
      CkVec<GenericTreeNode *> *ichlist = new CkVec<GenericTreeNode *>();
      CkVec<GenericTreeNode *> *echlist = new CkVec<GenericTreeNode *>();
      CkVec<GenericTreeNode *> *clistl = new CkVec<GenericTreeNode *>();
      CkVec<PartInfo> *plistl = new CkVec<PartInfo>();
      CkVec<GenericTreeNode *> *chlistl = new CkVec<GenericTreeNode *>();
  
      cellList.push_back(*clist);
      //delete clist; 
      particleList.push_back(*plist);
      intCheckList.push_back(*ichlist);
      extCheckList.push_back(*echlist);
      cellListLocal.push_back(*clistl);
      particleListLocal.push_back(*plistl);
      checkListLocal.push_back(*chlistl);
    }*/
    //listMigrated=false;
  //}
  cellList.resize(myTreeLevels+1);
  particleList.resize(myTreeLevels+1);
  intCheckList.resize(myTreeLevels+1);
  extCheckList.resize(myTreeLevels+1);
  cellListLocal.resize(myTreeLevels+1);
  particleListLocal.resize(myTreeLevels+1);
  checkListLocal.resize(myTreeLevels+1);
#endif
  
  for (int i=0; i<numChunks; ++i) remainingChunk[i] = myNumParticles;

  //BucketGravityRequest req0(1);
  //BucketGravityRequest req1(1);
  switch(domainDecomposition){
    case Oct_dec:
    case ORB_dec:
      //Prefetch Roots for Oct
      prefetchReq[0].reset();
      for (int i=1; i<myNumParticles; ++i) {
	prefetchReq[0].grow(myParticles[i].position);
      }
      /*
      req0.positions[0] = myParticles[1].position;
      for(int i=1;i<=myNumParticles;i++)
        req0.boundingBox.grow(myParticles[i].position);
      prefetchReq[0] = req0;
      req1.positions[0] = myParticles[myNumParticles].position;
      req1.boundingBox.grow(myParticles[myNumParticles].position);
      prefetchReq[1] = req1;
      */
      break;
    default:
      //Prefetch Roots for SFC
      prefetchReq[0].reset();
      prefetchReq[0].grow(myParticles[1].position);
      prefetchReq[1].reset();
      prefetchReq[1].grow(myParticles[myNumParticles].position);
      /*
      req0.positions[0] = myParticles[1].position;
      req0.boundingBox.grow(myParticles[1].position);
      prefetchReq[0] = req0;
      req1.positions[0] = myParticles[myNumParticles].position;
      req1.boundingBox.grow(myParticles[myNumParticles].position);
      prefetchReq[1] = req1;
      */
      break;
  }

  prefetchWaiting = 1;
  currentPrefetch = 0;
  int first, last;
#if CACHE_TREE
  GenericTreeNode *child = localCache->chunkRootToNode(prefetchRoots[0]);
#else
  GenericTreeNode *child = keyToNode(prefetchRoots[0]);
#endif
  if (child == NULL) {
    nodeOwnership(prefetchRoots[0], first, last);
    child = requestNode((first+last)>>1, prefetchRoots[0], 0, -1, true);
  }
  if (child != NULL) prefetch(child);

  thisProxy[thisIndex].calculateGravityLocal();
}

void TreePiece::prefetch(GenericTreeNode *node) {
  ///@TODO: all the code that does the prefetching and the chunking
  CkAssert(node->getType() != Invalid);
  //printf("{%d-%d} prefetch %016llx in chunk %d\n",CkMyPe(),thisIndex,node->getKey(),currentPrefetch);

  if (_prefetch) {
    bool needOpened = false;
    for (unsigned int i=0; i<numPrefetchReq; ++i) {
      if (openCriterionBucket(node, prefetchReq[i])) {
	needOpened = true;
	break;
      }
    }
    if (node->getType() != Internal && node->getType() != Bucket && needOpened) {
      if(node->getType() == CachedBucket || node->getType() == NonLocalBucket) {
	// Sending the request for all the particles at one go, instead of one by one
	if (requestParticles(node->getKey(),currentPrefetch,node->remoteIndex,node->firstParticle,node->lastParticle,-1,true) == NULL) {
	  prefetchWaiting ++;
	}
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
	  prefetchWaiting ++;
	  
	  if (child) {
	    prefetch(child);
	  } else { //missed the cache
	    child = requestNode(node->remoteIndex, node->getChildKey(i), currentPrefetch, -1, true);
	    if (child) { // means that node was on a local TreePiece
	      prefetch(child);
	    }
	  }
	}
      }
    }
  }

  prefetchWaiting --;
  //if (prefetchWaiting==0) ckout <<"Waiting for "<<prefetchWaiting<<" more prefetches"<<endl;

  // this means we don't have any more nodes waiting for prefetching
  if (prefetchWaiting == 0) {
    startRemoteChunk();
  }
}

void TreePiece::prefetch(ExternalGravityParticle *node) {
  prefetchWaiting --;
  //if (prefetchWaiting==0) ckout <<"Waiting for "<<prefetchWaiting<<" more prefetches"<<endl;

  // this means we don't have any more nodes waiting for prefetching
  if (prefetchWaiting == 0) {
    startRemoteChunk();
  }
}

void TreePiece::startRemoteChunk() {
  ComputeChunkMsg *msg = new (8*sizeof(int)) ComputeChunkMsg(currentPrefetch);
  *(int*)CkPriorityPtr(msg) = numTreePieces * currentPrefetch + thisIndex + 1;
  CkSetQueueing(msg, CK_QUEUEING_IFIFO);


  thisProxy[thisIndex].calculateGravityRemote(msg);
  
  // start prefetching next chunk
  if (++currentPrefetch < numChunks) {
    int first, last;
    prefetchWaiting = 1;
#if CACHE_TREE
    GenericTreeNode *child = localCache->chunkRootToNode(prefetchRoots[currentPrefetch]);
#else
    GenericTreeNode *child = keyToNode(prefetchRoots[currentPrefetch]);
#endif
    if (child == NULL) {
      nodeOwnership(prefetchRoots[currentPrefetch], first, last);
      child = requestNode((first+last)>>1, prefetchRoots[currentPrefetch], currentPrefetch, -1, true);
    }
    if (child != NULL) prefetch(child);
  }
}

void TreePiece::startlb(CkCallback &cb){
  callback = cb;
  if(verbosity > 1)
    CkPrintf("[%d] TreePiece %d calling AtSync()\n",CkMyPe(),thisIndex);
  localCache->revokePresence(thisIndex);
  AtSync();
}

void TreePiece::ResumeFromSync(){
  if(verbosity > 1)
    CkPrintf("[%d] TreePiece %d in ResumefromSync\n",CkMyPe(),thisIndex);
  // Register to the cache manager
  if (_cache) {
    localCache->markPresence(thisIndex, root);
  }
  contribute(0, 0, CkReduction::concat, callback);
}

GenericTreeNode *TreePiece::keyToNode(const Tree::NodeKey k) {
  NodeLookupType::iterator iter = nodeLookupTable.find(k);
  if (iter != nodeLookupTable.end()) return iter->second;
  else return NULL;
}

const GenericTreeNode *TreePiece::lookupNode(Tree::NodeKey key){
  return keyToNode(key);
};

const GravityParticle *TreePiece::lookupParticles(int begin) {
  return &myParticles[begin];
}

/*
 * For cached version we have 2 walks: one for on processor and one
 * that hits the cache. This does the local computation
 * When remote data is needed we go to the second version.
 */

#if INTERLIST_VER > 0

//Does interaction list computation for each bucket by calling 'walkInterTree'
void TreePiece::preWalkInterTree(){

    //Start copying the checkList of previous level to the next level
    int level;
    GenericTreeNode *child;
    GenericTreeNode *node;
    int flag=0;
    NodeType childType;
    
    while(1){
      level=curLevelLocal-1;
      checkListLocal[curLevelLocal].length()=0;
      prevListIterLocal=0;
      
      if(curNodeLocal!=root){
	if(checkListLocal[level].length()!=0){
          node=checkListLocal[level][0];
          prevListIterLocal=1;
        }
        else node=NULL;
      }
      else{
        GenericTreeNode *nd;
        for(int i=0;i<numChunks;i++){
#ifdef CACHE_TREE
          nd = localCache->chunkRootToNode(prefetchRoots[i]);
#else
          nd = keyToNode(prefetchRoots[i]);
#endif
          if(nd!=NULL)
            undecidedListLocal.enq(nd);
        }
        CkAssert(!undecidedListLocal.isEmpty());
        node=undecidedListLocal.deq();
      }
    	
      //Walks the local tree for my current node
      if(node!=NULL){
#if INTERLIST_VER==1
        walkInterTreeVerI(node);
#else
	walkInterTreeVerII(node);
#endif
      }
      else{
        myLocalCheckListEmpty=true;
        curNodeLocal=curNodeLocal->parent;
        curLevelLocal--;
        break;
      }
      
      CkAssert(undecidedListLocal.isEmpty());
    
      //Loop breaking condition
      if(curNodeLocal->getType()==Bucket)
        break;
      
      //Finds my node on the next level which is not yet visited 
      for(int i=0;i<curNodeLocal->numChildren();i++){
	child = curNodeLocal->getChildren(i);
      	CkAssert (child != NULL);
        childType = child->getType();
	if(child->visitedL==false){
          if(childType == NonLocal || childType == Cached || childType == NonLocalBucket || childType == CachedBucket || childType==Empty || childType==CachedEmpty){
            child->visitedL=true;
          }
          else{
            flag=1;
	    break;
          }
        }
      }
      if(flag==1){
	curNodeLocal=child;
	curLevelLocal++;
        flag=0;
      }
      else{
        CkPrintf("Exceptional case\n");
        CkAssert(false);
        curNodeLocal->visitedL=true;
        break;
      }
    }
}
#endif

#if INTERLIST_VER==1

void TreePiece::walkInterTreeVerI(GenericTreeNode *node) {
  
  //Copy an element from previous level checkList only if everything hasn't been copied till now

  GenericTreeNode *myNode = curNodeLocal;
  int level=curLevelLocal;
  NodeType nodeType = node->getType();
  
#if COSMO_STATS > 0
  numOpenCriterionCalls++;
#endif
  if(openCriterionNode(node, myNode)==0) {
    cellListLocal[level].push_back(node);
  } else if(nodeType == Bucket) {
#if COSMO_STATS > 0
    if(myNode->getType()!=Bucket)
      numOpenCriterionCalls++;
#endif
    if(myNode->getType()==Bucket || openCriterionNode(node,myNode)==1){
      TreePiece::LocalPartInfo pinfo;
      pinfo.particles = &myParticles[node->firstParticle];
#if COSMO_DEBUG>1
      pinfo.nd=node;
#endif
      pinfo.numParticles = node->lastParticle - node->firstParticle + 1;
      particleListLocal[level].push_back(pinfo);
    }
    else{
      checkListLocal[level].push_back(node);
    }
  } else if(nodeType == NonLocal || nodeType == NonLocalBucket) {
#if COSMO_STATS > 0
    if(myNode->getType()!=Bucket)
      numOpenCriterionCalls++;
#endif
    if(myNode->getType()==Bucket || openCriterionNode(node,myNode)==1){
    }
    else{
      checkListLocal[level].push_back(node);
    }
    
    /* DISABLED: this part of the walk is triggered directly by the CacheManager and prefetching
     */
  } else if (nodeType != Empty) {
    // here the node can be Internal or Boundary
#if COSMO_STATS > 0
    if(myNode->getType()!=Bucket)
      numOpenCriterionCalls++;
#endif
    if(myNode->getType()==Bucket || openCriterionNode(node,myNode)==1){
      GenericTreeNode* childIterator;
      for(unsigned int i = 0; i < node->numChildren(); ++i) {
        childIterator = node->getChildren(i);
        if(childIterator)
          undecidedListLocal.enq(childIterator);
      }
    }
    else{
      checkListLocal[level].push_back(node);
    }
  }
  else if (nodeType == Empty) {
    /*#if COSMO_DEBUG > 1
      bucketcheckList[req.identifier].insert(tmp->getKey());
      combineKeys(tmp->getKey(),req.identifier);
    #endif*/
  	//cellListLocal[level].push_back(node);
  }
  
  //if(undecidedListLocal.isEmpty())
    //return;

  if(myNode!=root){
    if(prevListIterLocal>0 && prevListIterLocal<checkListLocal[level-1].length()){
      prevListIterLocal++;
      walkInterTreeVerI(checkListLocal[level-1][prevListIterLocal-1]);
    }
    else if(!undecidedListLocal.isEmpty())
      walkInterTreeVerI(undecidedListLocal.deq());
    else
      return;
  }
  else{
    if(!undecidedListLocal.isEmpty())
      walkInterTreeVerI(undecidedListLocal.deq());
    else
      return;
  }
}
    
#endif

#if INTERLIST_VER==2

void TreePiece::walkInterTreeVerII(GenericTreeNode *node) {
  
  GenericTreeNode *myNode = curNodeLocal;
  int level=curLevelLocal;
  NodeType nodeType = node->getType();
  
  int openValue=-2;

  if(nodeType == NonLocal || nodeType == NonLocalBucket) {
   /* DISABLED: this part of the walk is triggered directly by the CacheManager and prefetching
   */
  } else if(nodeType==Empty){
#ifdef CACHE_TREE 
    if (thisProxy[node->remoteIndex].ckLocal()!=NULL) {
#else
    if (node->remoteIndex==thisIndex) {
#endif
#if COSMO_STATS > 0
      numOpenCriterionCalls++;
#endif
#if COSMO_DEBUG > 1
      cellListLocal[level].push_back(node);
#endif
    }
    else{}
  } else if((openValue=openCriterionNode(node, myNode))==0) {
#if COSMO_STATS > 0
    numOpenCriterionCalls++;
#endif
    if(nodeType!=Boundary)
      cellListLocal[level].push_back(node);
  }  else if (nodeType != Empty) {
#if COSMO_STATS > 0
    numOpenCriterionCalls++;
#endif
    /*#if COSMO_STATS > 0
      if(myNode->getType()!=Bucket)
      numOpenCriterionCalls++;
      #endif*/
    CkAssert(openValue!=-2);
    // here the node can be Internal or Boundary or Bucket
    if(myNode->getType()==Bucket || openValue==1){
      if(nodeType==Bucket){
        TreePiece::LocalPartInfo pinfo;
        //pinfo.particles = &myParticles[node->firstParticle];
        pinfo.particles = node->particlePointer;
#if COSMO_DEBUG>1
        pinfo.nd=node;
#endif
        pinfo.numParticles = node->lastParticle - node->firstParticle + 1;
        particleListLocal[level].push_back(pinfo);
      }
      else{
        GenericTreeNode* childIterator;
        for(unsigned int i = 0; i < node->numChildren(); ++i) {
          childIterator = node->getChildren(i);
          if(childIterator)
            undecidedListLocal.enq(childIterator);
        }
      }
    }
    else{
      checkListLocal[level].push_back(node);
    }
  }
  /*else if (nodeType == Empty) {
#if COSMO_STATS > 0
    numOpenCriterionCalls++;
#endif
#if COSMO_DEBUG > 1
  	cellListLocal[level].push_back(node);
#endif
  }*/
 
  //Call myself if there are still nodes in previous level checklist
  //or my undecided list
  if(myNode!=root){
    if(prevListIterLocal>0 && prevListIterLocal<checkListLocal[level-1].length()){
      prevListIterLocal++;
      walkInterTreeVerII(checkListLocal[level-1][prevListIterLocal-1]);
    }
    else if(!undecidedListLocal.isEmpty())
      walkInterTreeVerII(undecidedListLocal.deq());
    else
      return;
  }
  else{
    if(!undecidedListLocal.isEmpty())
      walkInterTreeVerII(undecidedListLocal.deq());
    else
      return;
  }
}

#endif


void TreePiece::walkBucketTree(GenericTreeNode* node, int reqID) {
#if COSMO_STATS > 0
  myNumMACChecks++;
#endif
  //Vector3D<double> cm(node->moments.cm);
  //Vector3D<double> r;
  //Sphere<double> s(cm, opening_geometry_factor * node->moments.radius / theta);
#if COSMO_PRINT > 0
  if (node->getKey() < numChunks) CkPrintf("[%d] walk: checking %016llx\n",thisIndex,node->getKey());
#endif
#if COSMO_STATS > 0
  numOpenCriterionCalls++;
#endif
  GenericTreeNode* reqnode = bucketList[reqID];
  if(!openCriterionBucket(node, reqnode->boundingBox)) {
    nodeInterLocal += bucketList[reqID]->particleCount;
#if COSMO_STATS > 1
    MultipoleMoments m = node->moments;	
    for(int i = reqnode->firstParticle; i <= reqnode->lastParticle; ++i)
      myParticles[i].intcellmass += m.totalMass;
#endif
#if COSMO_PRINT > 1
  CkPrintf("[%d] walk bucket %s -> node %s\n",thisIndex,keyBits(bucketList[reqID]->getKey(),63).c_str(),keyBits(node->getKey(),63).c_str());
#endif
#if COSMO_DEBUG > 1
  bucketcheckList[reqID].insert(node->getKey());
  combineKeys(node->getKey(),reqID);
#endif
#if COSMO_PRINT > 0
  if (node->getKey() < numChunks) CkPrintf("[%d] walk: computing %016llx\n",thisIndex,node->getKey());
#endif
    nodeBucketForce(node, reqnode, myParticles);
  } else if(node->getType() == Bucket) {
    particleInterLocal += reqnode->particleCount * (node->lastParticle - node->firstParticle + 1);
    for(int i = node->firstParticle; i <= node->lastParticle; ++i) {
#if COSMO_STATS > 1
      for(int j = reqnode->firstParticle; j <= reqnode->lastParticle; ++j) {
	//myParticles[j].intpartmass += myParticles[i].mass;
	myParticles[j].intpartmass += node->particlePointer[i-node->firstParticle].mass;
      }
#endif
#if COSMO_PRINT > 1
      //CkPrintf("[%d] walk bucket %s -> part %016llx\n",thisIndex,keyBits(reqnode->getKey(),63).c_str(),myParticles[i].key);
      CkPrintf("[%d] walk bucket %s -> part %016llx\n",thisIndex,keyBits(reqnode->getKey(),63).c_str(),node->particlePointer[i-node->firstParticle].key);
#endif
      //partBucketForce(&myParticles[i], req);
      partBucketForce(&node->particlePointer[i-node->firstParticle], reqnode, myParticles);
    }
#if COSMO_DEBUG > 1
  bucketcheckList[reqID].insert(node->getKey());
  combineKeys(node->getKey(),reqID);
#endif
  } else if (node->getType() == NonLocal || node->getType() == NonLocalBucket) {
    /* DISABLED: this part of the walk is triggered directly by the CacheManager and prefetching

    // Use cachedWalkBucketTree() as callback
    GenericTreeNode *pnode = requestNode(node->remoteIndex, node->getKey(), req);
    if(pnode) {
      cachedWalkBucketTree(pnode, req);
    }
    */
  } else if (node->getType() != Empty) {
    // here the node can be Internal or Boundary
#if COSMO_STATS > 0
    nodesOpenedLocal++;
#endif
#if COSMO_PRINT > 0
    if (node->getKey() < numChunks) CkPrintf("[%d] walk: opening %016llx\n",thisIndex,node->getKey());
#endif
    GenericTreeNode* childIterator;
    for(unsigned int i = 0; i < node->numChildren(); ++i) {
      childIterator = node->getChildren(i);
      if(childIterator)
	walkBucketTree(childIterator, reqID);
    }
  }
}

 /*
 * Cached version of Tree walk. One characteristic of the tree used is that once
 * we go into cached data, we cannot come back to internal data anymore. Thus we
 * can safely distinguish between local computation done by walkBucketTree and
 * remote computation done by cachedWalkBucketTree.
 */

#if INTERLIST_VER==1

void TreePiece::cachedWalkInterTreeVerI(GenericTreeNode* node) {
  
	int chunk = nChunk;
  GenericTreeNode* myNode = curNodeRemote;
  int level = curLevelRemote;

	NodeType nodeType = node->getType();

  CkAssert(nodeType != Invalid);

#if COSMO_STATS > 0
  numOpenCriterionCalls++;
#endif

  if(openCriterionNode(node, myNode)==0) {
    if(!(nodeType==NonLocal && node->parent->remoteIndex==thisIndex))
      cellList[level].push_back(node);
  } else if(nodeType == CachedBucket || nodeType == Bucket || nodeType == NonLocalBucket) {
    /*
     * Sending the request for all the particles at one go, instead of one by one
     */
    if(myNode->getType()==Bucket){
      ExternalGravityParticle *part = requestParticles(node->getKey(),chunk,node->remoteIndex,node->firstParticle,node->lastParticle,bucketReqs[myNode->bucketListIndex]);
      //PROBLEM: how to form req in this case
      if(part != NULL){
      	TreePiece::RemotePartInfo pinfo;
      	pinfo.particles = part;
#if COSMO_DEBUG>1
        pinfo.nd=node;
#endif
      	pinfo.numParticles = node->lastParticle - node->firstParticle + 1;
      	particleList[level].push_back(pinfo);
    	} else {
      	remainingChunk[chunk] += node->lastParticle - node->firstParticle + 1;
        /*#if COSMO_DEBUG > 1
        bucketReqs[myNode->bucketListIndex].requestedNodes.push_back(node->getKey());
        #endif*/
    	}
    }
    else if(openCriterionNode(node,myNode)==1){
      calculateForces(node,myNode,level,chunk);
#if COSMO_STATS > 0
      numOpenCriterionCalls++;
#endif
    }
    else{
      extCheckList[level].push_back(node);
#if COSMO_STATS > 0
      numOpenCriterionCalls++;
#endif
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
  } else if (nodeType != CachedEmpty && nodeType != Empty) {
    // Here the type is Cached, Boundary, Internal, NonLocal, which means the
    // node in the global tree has children (it is not a leaf), so we iterate
    // over them. If we get a NULL node, then we missed the cache and we request
    // it

    // Warning, since the cache returns nodes with pointers to other chare
    // elements trees, we could be actually traversing the tree of another chare
    // in this processor.

    // Use cachedWalkBucketTree() as callback
		if(myNode->getType()==Bucket){
    	GenericTreeNode *child;
    	//BucketGravityRequest& req = bucketReqs[myNode->bucketListIndex];
    	for (unsigned int i=0; i<node->numChildren(); ++i) {
      	child = node->getChildren(i); //requestNode(node->remoteIndex, node->getChildKey(i), req);
      	if (child) {
          undecidedExtList.enq(child);
      	} else { //missed the cache
					//PROBLEM: how to construct req
	      	child = requestNode(node->remoteIndex, node->getChildKey(i), chunk, bucketReqs[myNode->bucketListIndex]);
	      	if (child) { // means that node was on a local TreePiece
            undecidedExtList.enq(child);
	      	} else { // we completely missed the cache, we will be called back
	        	remainingChunk[chunk] ++;
	      	}
      	}
    	}
		}
		else if(openCriterionNode(node,myNode)==1){
			calculateForcesNode(node,myNode,level,chunk);
      #if COSMO_STATS > 0
        numOpenCriterionCalls++;
      #endif
		}
		else{
      extCheckList[level].push_back(node);
    #if COSMO_STATS > 0
      numOpenCriterionCalls++;
    #endif
    }
  }
  else if(nodeType == CachedEmpty || nodeType == Empty) {
    //Currently, Empty node is being pushed back just to match Filippos nodeBucketCalls
    //I'll remove it later
    //cellList[level].push_back(node);
  }
  
  if(myNode!=root){
    if(extPrevListIter>0 && extPrevListIter<extCheckList[level-1].length()){
      extPrevListIter++;
      cachedWalkInterTreeVerI(extCheckList[level-1][extPrevListIter-1]);
    }
    else if(!undecidedExtList.isEmpty())
      cachedWalkInterTreeVerI(undecidedExtList.deq());
    else
      return;
  }
  else{
    if(!undecidedExtList.isEmpty())
      cachedWalkInterTreeVerI(undecidedExtList.deq());
    else
      return;
  }
}
#endif

#if INTERLIST_VER==2

void TreePiece::cachedWalkInterTreeVerII(GenericTreeNode* node) {
  
  int chunk = nChunk;
  GenericTreeNode* myNode = curNodeRemote;
  int level = curLevelRemote;

  NodeType nodeType = node->getType();
  CkAssert(nodeType != Invalid);

  int openValue=-2;
  
    // Changed for CACHE TREE
#ifdef CACHE_TREE
  if(nodeType == Bucket || nodeType == Internal || (nodeType == Empty && thisProxy[node->remoteIndex].ckLocal()!=NULL)) {
  }
#else
  if(node->remoteIndex==thisIndex && (nodeType == Bucket || nodeType == Internal || nodeType == Empty)) {
  }
#endif
  else if(nodeType == CachedEmpty || nodeType == Empty) {
    //Currently, Empty node is being pushed back just to match Filippos nodeBucketCalls
    //I'll remove it later
#if COSMO_DEBUG > 1
    cellList[level].push_back(node);
#endif
#if COSMO_STATS > 0
    numOpenCriterionCalls++;
#endif
  }
  else if((openValue=openCriterionNode(node, myNode))==0) {
#if COSMO_STATS > 0
    numOpenCriterionCalls++;
#endif
    cellList[level].push_back(node);
  } else if(nodeType == CachedBucket || nodeType == Bucket || nodeType == NonLocalBucket) {
    /*
     * Sending the request for all the particles at one go, instead of one by one
     */
#ifdef CACHE_TREE
    CkAssert(nodeType!=Bucket);
#endif
#if COSMO_STATS > 0
    numOpenCriterionCalls++;
#endif
    if(myNode->getType()==Bucket){
      ExternalGravityParticle *part = requestParticles(node->getKey(),chunk,node->remoteIndex,node->firstParticle,node->lastParticle,myNode->bucketListIndex);
      //PROBLEM: how to form req in this case
      if(part != NULL){
      	TreePiece::RemotePartInfo pinfo;
      	pinfo.particles = part;
#if COSMO_DEBUG>1
        pinfo.nd=node;
#endif
      	pinfo.numParticles = node->lastParticle - node->firstParticle + 1;
      	particleList[level].push_back(pinfo);
      } else {
      	remainingChunk[chunk] += node->lastParticle - node->firstParticle + 1;
        /*#if COSMO_DEBUG > 1
        bucketReqs[myNode->bucketListIndex].requestedNodes.push_back(node->getKey());
        #endif*/
      }
    }
    else if(openValue==1){
      calculateForces(node,myNode,level,chunk);
    }
    else{
      extCheckList[level].push_back(node);
    }
/*#ifdef CACHE_TREE
  } else if(nodeType==Boundary) {
#else
  } else if(node->remoteIndex==thisIndex && nodeType==Boundary) {
#endif
#if COSMO_STATS > 0
  numOpenCriterionCalls++;
#endif
    if(openValue==1 || myNode->getType()==Bucket){
      // here the node is Boundary
#if COSMO_STATS > 0
      nodesOpenedRemote++;
#endif
		  //Boundary node removed from the interaction list...we dont calculate force with boundary
      GenericTreeNode* childIterator;
      for(unsigned int i = 0; i < node->numChildren(); ++i) {
        childIterator = node->getChildren(i);
        CkAssert (childIterator != NULL);
        undecidedExtList.enq(childIterator);
      }
    }
    else{
      //intListIter++;
      extCheckList[level].push_back(node);
    }*/
  } else if (nodeType != CachedEmpty && nodeType != Empty) {
    // Here the type is Cached, Boundary, Internal, NonLocal, which means the
    // node in the global tree has children (it is not a leaf), so we iterate
    // over them. If we get a NULL node, then we missed the cache and we request
    // it

    // Warning, since the cache returns nodes with pointers to other chare
    // elements trees, we could be actually traversing the tree of another chare
    // in this processor.
#if COSMO_STATS > 0
    numOpenCriterionCalls++;
#endif

    // Use cachedWalkBucketTree() as callback
    if(myNode->getType()==Bucket){
      GenericTreeNode *child;
      //BucketGravityRequest& req = bucketReqs[myNode->bucketListIndex];
      for (unsigned int i=0; i<node->numChildren(); ++i) {
	child = node->getChildren(i); //requestNode(node->remoteIndex, node->getChildKey(i), req);
      	if (child) {
          undecidedExtList.enq(child);
      	} else { //missed the cache
	  //PROBLEM: how to construct req
	  child = requestNode(node->remoteIndex, node->getChildKey(i), chunk, myNode->bucketListIndex);
	  if (child) { // means that node was on a local TreePiece
            undecidedExtList.enq(child);
	  } else { // we completely missed the cache, we will be called back
	    remainingChunk[chunk] ++;
	  }
      	}
      }
    }
    else if(openValue==1){
      calculateForcesNode(node,myNode,level,chunk);
    }
    else{
      extCheckList[level].push_back(node);
    }
  }

  if(myNode!=root){
    if(extPrevListIter>0 && extPrevListIter<extCheckList[level-1].length()){
      extPrevListIter++;
      cachedWalkInterTreeVerII(extCheckList[level-1][extPrevListIter-1]);
    }
    else if(!undecidedExtList.isEmpty())
      cachedWalkInterTreeVerII(undecidedExtList.deq());
    else
      return;
  }
  else{
    if(!undecidedExtList.isEmpty())
      cachedWalkInterTreeVerII(undecidedExtList.deq());
    else
      return;
  }
}
#endif

#if INTERLIST_VER > 0
inline void TreePiece::calculateForces(GenericTreeNode *node, GenericTreeNode *myNode,int level,int chunk){

  GenericTreeNode *tmpNode;
  int startBucket=myNode->startBucket;
  int tmpBucket;
  int lastBucket;
  int i;
  
  for(i=startBucket+1;i<numBuckets;i++){
    tmpNode = bucketList[i];
    if(tmpNode->lastParticle>myNode->lastParticle)
      break;
  }
  lastBucket=i-1;
  
  int test=0;
  
  for(i=startBucket;i<=lastBucket;i++){
    ExternalGravityParticle *part = requestParticles(node->getKey(),chunk,node->remoteIndex,node->firstParticle,node->lastParticle,i);
    if(part != NULL){
      CkAssert(test==0);
      TreePiece::RemotePartInfo pinfo;
      pinfo.particles = part;
#if COSMO_DEBUG>1
      pinfo.nd=node;
#endif
      pinfo.numParticles = node->lastParticle - node->firstParticle + 1;
      particleList[level].push_back(pinfo);
      break;
    } else {
      remainingChunk[chunk] += node->lastParticle - node->firstParticle + 1;
      /*#if COSMO_DEBUG > 1
      bucketReqs[i].requestedNodes.push_back(node->getKey());
      #endif*/
    }
    test++;
  }
}

inline void TreePiece::calculateForcesNode(GenericTreeNode *node, GenericTreeNode *myNode,int level,int chunk){

  GenericTreeNode *tmpNode;
  int startBucket=myNode->startBucket;
  int k;
  
	/*for(k=startBucket+1;k<numBuckets;k++){
		tmpNode = bucketList[k];
		if(tmpNode->lastParticle>myNode->lastParticle)
			break;
  }
  lastBucket=k-1;*/
  
  GenericTreeNode *child;
  for (unsigned int i=0; i<node->numChildren(); ++i) {
    child = node->getChildren(i); //requestNode(node->remoteIndex, node->getChildKey(i), req);
    k = startBucket;
    
    if (child) {
      undecidedExtList.enq(child);
    } else { //missed the cache
      //PROBLEM: how to construct req
      child = requestNode(node->remoteIndex, node->getChildKey(i), chunk, k);
      if (child) { // means that node was on a local TreePiece
        undecidedExtList.enq(child);
      } else { // we completely missed the cache, we will be called back
        //We have to queue the requests for all the buckets, all buckets will be called back later
        
        remainingChunk[chunk] ++;
        for(k=startBucket+1;k<numBuckets;k++){
          tmpNode = bucketList[k];
          if(tmpNode->lastParticle>myNode->lastParticle) break;
          child = requestNode(node->remoteIndex, node->getChildKey(i), chunk, k);
          CkAssert(child==NULL);
          remainingChunk[chunk] ++;
        }
      }
    }
  }
}

#endif

void TreePiece::cachedWalkBucketTree(GenericTreeNode* node, int chunk, int reqID) {
  GenericTreeNode *reqnode = bucketList[reqID];
#if COSMO_STATS > 0
  myNumMACChecks++;
#endif
#if COSMO_PRINT > 1
  CkPrintf("[%d] b=%d cachedWalkBucketTree called for %s with node %s of type %s (additional=%d)\n",thisIndex,reqID,keyBits(bucketList[reqID]->getKey(),63).c_str(),keyBits(node->getKey(),63).c_str(),getColor(node).c_str(),bucketReqs[reqID].numAdditionalRequests);
#endif
		
  CkAssert(node->getType() != Invalid);
	
#if COSMO_STATS > 0
  numOpenCriterionCalls++;
#endif
  if(!openCriterionBucket(node, reqnode->boundingBox)) {
    nodeInterRemote[chunk] += reqnode->particleCount;
#if COSMO_STATS > 1
    MultipoleMoments m = node->moments;
    for(int i = reqnode->firstParticle; i <= reqnode->lastParticle; ++i)
      myParticles[i].extcellmass += m.totalMass;
#endif
#if COSMO_PRINT > 1
  CkPrintf("[%d] cachedwalk bucket %s -> node %s\n",thisIndex,keyBits(reqnode->getKey(),63).c_str(),keyBits(node->getKey(),63).c_str());
#endif
#if COSMO_DEBUG > 1
  bucketcheckList[reqID].insert(node->getKey());
  combineKeys(node->getKey(),reqID);
#endif
    nodeBucketForce(node, reqnode, myParticles);
  } else if(node->getType() == CachedBucket || node->getType() == Bucket || node->getType() == NonLocalBucket) {
    /*
     * Sending the request for all the particles at one go, instead of one by one
     */
    //printf("{%d-%d} cachewalk requests for %016llx in chunk %d\n",CkMyPe(),thisIndex,node->getKey(),chunk);
    ExternalGravityParticle *part = requestParticles(node->getKey(),chunk,node->remoteIndex,node->firstParticle,node->lastParticle,reqID);
    if(part != NULL){
      particleInterRemote[chunk] += reqnode->particleCount * (node->lastParticle - node->firstParticle + 1);

      for(int i = node->firstParticle; i <= node->lastParticle; ++i) {
#if COSMO_STATS > 1
	for(int j = reqnode->firstParticle; j <= reqnode->lastParticle; ++j) {
	  myParticles[j].extpartmass += myParticles[i].mass;
	}
#endif
#if COSMO_PRINT > 1
	CkPrintf("[%d] cachedwalk bucket %s -> part %016llx\n",thisIndex,keyBits(reqnode->getKey(),63).c_str(),part[i-node->firstParticle].key);
#endif
	partBucketForce(&part[i-node->firstParticle], reqnode, myParticles);
      }
#if COSMO_DEBUG > 1
      bucketcheckList[reqID].insert(node->getKey());
      combineKeys(node->getKey(),reqID);
#endif
    } else {
      remainingChunk[chunk] += node->lastParticle - node->firstParticle + 1;
      /*#if COSMO_DEBUG > 1
      bucketReqs[req.identifier].requestedNodes.push_back(node->getKey());
      #endif*/
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

#if COSMO_STATS > 0
    nodesOpenedRemote++;
#endif
    // Use cachedWalkBucketTree() as callback
    GenericTreeNode *child;
    for (unsigned int i=0; i<node->numChildren(); ++i) {
      child = node->getChildren(i); //requestNode(node->remoteIndex, node->getChildKey(i), req);
      if (child) {
	cachedWalkBucketTree(child, chunk, reqID);
      } else { //missed the cache
	child = requestNode(node->remoteIndex, node->getChildKey(i), chunk, reqID);
	if (child) { // means that node was on a local TreePiece
	  cachedWalkBucketTree(child, chunk, reqID);
	} else { // we completely missed the cache, we will be called back
	  remainingChunk[chunk] ++;
	}
      }
    }
  }
}


GenericTreeNode* TreePiece::requestNode(int remoteIndex, Tree::NodeKey key, int chunk, int reqID, bool isPrefetch) {

  // Call proxy on remote node
  CkAssert(remoteIndex < (int) numTreePieces);
  CkAssert(chunk < numChunks);
  //in the current form it is possible   
  //   assert(remoteIndex != thisIndex);
  if(_cache){
    CkAssert(localCache != NULL);
    /*
    if(localCache == NULL){
      localCache = cacheManagerProxy.ckLocalBranch();
    }
    */
#if COSMO_PRINT > 1
    CkPrintf("[%d] b=%d requesting node %s to %d for %s (additional=%d)\n",thisIndex,reqID,keyBits(key,63).c_str(),remoteIndex,keyBits(bucketList[reqID]->getKey(),63).c_str(),bucketReqs[reqID].numAdditionalRequests);
#endif
    GenericTreeNode *res=localCache->requestNode(thisIndex,remoteIndex,chunk,key,reqID,isPrefetch);
    if(!res){
      bucketReqs[reqID].numAdditionalRequests++;
      //#if COSMO_STATS > 0
      //myNumProxyCalls++;
      //#endif
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

void TreePiece::fillRequestNode(RequestNodeMsg *msg) {
  GenericTreeNode* node = keyToNode(msg->key);
  //GenericTreeNode tmp;
  if(node != NULL) {
    if(_cache) {
      PUP::sizer p1;
      node->pup(p1, msg->depth);
      FillNodeMsg *reply = new (p1.size(), 0) FillNodeMsg(thisIndex);

      /// @TODO: check that at destination of "remoteIndex" are correct
      PUP::toMem p2((void*)reply->nodes);
      node->pup(p2, msg->depth);
      //int count = node->copyTo(reply->nodes, msg->depth);
      cacheManagerProxy[msg->retIndex].recvNodes(reply);
    }else{
      CkAbort("Non cached version not anymore supported, feel free to fix it!");
      //copySFCTreeNode(tmp,node);
      //streamingProxy[retIndex].receiveNode(tmp,msg->reqID);
    }
  }
  else {	// Handle NULL nodes
    CkAbort("Ok, before it handled this, but why do we have a null pointer in the tree?!?");
  }
  delete msg;
}

void TreePiece::receiveNode(GenericTreeNode &node, int chunk, unsigned int reqID)
{
#if COSMO_PRINT > 1
  CkPrintf("[%d] b=%d, receiveNode, additional=%d\n",thisIndex,reqID,bucketReqs[reqID].numAdditionalRequests);
#endif
  bucketReqs[reqID].numAdditionalRequests--;
  remainingChunk[chunk] --;
  assert(node.getType() != Invalid);
  if(node.getType() != Empty)	{ // Node could be NULL
    assert((int) node.remoteIndex != thisIndex);
    cachedWalkBucketTree(&node, chunk, reqID);
  }else{
#if COSMO_DEBUG > 1
    bucketcheckList[reqID].insert(node.getKey());
    combineKeys(node.getKey(),reqID);
#endif
  }
    
  finishBucket(reqID);
  CkAssert(remainingChunk[chunk] >= 0);
  if (remainingChunk[chunk] == 0) {
#ifdef COSMO_PRINT
    CkPrintf("[%d] Finished chunk %d with a node\n",thisIndex,chunk);
#endif
    cacheManagerProxy[CkMyPe()].finishedChunk(chunk, nodeInterRemote[chunk]+particleInterRemote[chunk]);
  }
}

void TreePiece::receiveNode_inline(GenericTreeNode &node, int chunk, unsigned int reqID){
        receiveNode(node,chunk,reqID);
}

ExternalGravityParticle *TreePiece::requestParticles(const Tree::NodeKey &key,int chunk,int remoteIndex,int begin,int end,int reqID, bool isPrefetch) {
  if (_cache) {
    CkAssert(localCache != NULL);
    /*
    if(localCache == NULL){
      localCache = cacheManagerProxy.ckLocalBranch();
    }
    */
    ExternalGravityParticle *p = localCache->requestParticles(thisIndex,chunk,key,remoteIndex,begin,end,reqID,isPrefetch);
    if (!p) {
#if COSMO_PRINT > 1
      CkPrintf("[%d] b=%d requestParticles: additional=%d\n",thisIndex,reqID,bucketReqs[reqID].numAdditionalRequests);
#endif
      bucketReqs[reqID].numAdditionalRequests += end-begin+1;
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

void TreePiece::fillRequestParticles(RequestParticleMsg *msg) {
  if (_cache) {
    PUP::sizer p1;
    for (int i=msg->begin; i<=msg->end; ++i) myParticles[i].ExternalGravityParticle::pup(p1);
    FillParticleMsg *reply = new (p1.size()) FillParticleMsg(thisIndex, msg->key, msg->end - msg->begin + 1);

    PUP::toMem p2((void*)reply->particles);
    for (int i=msg->begin; i<=msg->end; ++i) myParticles[i].ExternalGravityParticle::pup(p2);
    cacheManagerProxy[msg->retIndex].recvParticles(reply);
    
    //cacheManagerProxy[msg->retIndex].recvParticles(msg->key, &myParticles[msg->begin], msg->end - msg->begin + 1, thisIndex);
  } else {
    CkAbort("Non cached version not anymore supported!");
    //streamingProxy[msg->retIndex].receiveParticles(&myParticles[msg->begin], msg->end - msg->begin + 1, 0, msg->reqID);
  }	
  delete msg;
}

void TreePiece::receiveParticles(ExternalGravityParticle *part,int num,int chunk,
				 unsigned int reqID, Tree::NodeKey remoteBucketID)
{
  CkAssert(num > 0);
#if COSMO_PRINT > 1
  CkPrintf("[%d] b=%d recvPart (additional=%d-%d)\n",thisIndex,reqID,bucketReqs[reqID].numAdditionalRequests,num);
#endif
  bucketReqs[reqID].numAdditionalRequests -= num;
  remainingChunk[chunk] -= num;
  particleInterRemote[chunk] += bucketList[reqID]->particleCount * num;

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
    partBucketForce(&part[i], reqnode, myParticles);
  }
#if COSMO_DEBUG > 1
  /*
  Key mask = Key(~0);
  Tree::NodeKey reqNodeKey;
  std::vector<Tree::NodeKey>::iterator iter;
  //std::vector<Tree::NodeKey> reqNodes = bucketReqs[reqID].requestedNodes;

  if(bucketReqs[reqID].requestedNodes.empty())
    CkPrintf("Error: [%d] bucket:%d has it's list as empty",thisIndex,reqID);
  for(iter = bucketReqs[reqID].requestedNodes.begin(); iter != bucketReqs[reqID].requestedNodes.end(); iter++){
    mask = Key(~0);
    reqNodeKey = (*iter);
    const Key subMask = Key(1) << 63;
    while(!(reqNodeKey & subMask)){
      reqNodeKey <<= 1;
      mask <<= 1;
    }
    reqNodeKey &= ~subMask;
    mask &= ~subMask;
    Key k = part[0].key & mask;
    
    if(k == reqNodeKey){
      break;
    }
  }
  CkAssert(iter!=bucketReqs[reqID].requestedNodes.end());
  bucketcheckList[reqID].insert(*iter);
  combineKeys(*iter,reqID);
  bucketReqs[reqID].requestedNodes.erase(iter);
  */
  bucketcheckList[reqID].insert(remoteBucketID);
  combineKeys(remoteBucketID,reqID);
#endif
  finishBucket(reqID);
  CkAssert(remainingChunk[chunk] >= 0);
  if (remainingChunk[chunk] == 0) {
#ifdef COSMO_PRINT
    CkPrintf("[%d] Finished chunk %d with particle\n",thisIndex,chunk);
#endif
    cacheManagerProxy[CkMyPe()].finishedChunk(chunk, nodeInterRemote[chunk]+particleInterRemote[chunk]);
  }
}

void TreePiece::receiveParticles_inline(ExternalGravityParticle *part,int num,int chunk,
					unsigned int reqID, Tree::NodeKey remoteBucketID){
        receiveParticles(part,num,chunk,reqID,remoteBucketID);
}

#if COSMO_DEBUG > 1

//Recursive routine to combine keys -- Written only for Binary Trees
void TreePiece::combineKeys(Tree::NodeKey key,int bucket){

  Tree::NodeKey mask = Key(1);
  Tree::NodeKey lastBit = key & mask;
  Tree::NodeKey sibKey;
  
  if(lastBit==mask){
    sibKey = key >> 1;
    sibKey <<= 1;
  }
  else{
    sibKey = key | mask;
  }

  std::set<Tree::NodeKey>::iterator iter = (bucketcheckList[bucket]).find(sibKey);

  if(iter==bucketcheckList[bucket].end())
    return;
  else{//Sibling key has been found in the Binary tree
    bucketcheckList[bucket].erase(key);
    bucketcheckList[bucket].erase(sibKey);
    key >>= 1;
    bucketcheckList[bucket].insert(key);
    combineKeys(key,bucket);
  }
}

void TreePiece::checkWalkCorrectness(){

  Tree::NodeKey endKey = Key(1);
  CkPrintf("[%d]checking walk correctness...\n",thisIndex);
  for(int i=0;i<numBuckets;i++){
    if(bucketcheckList[i].size()!=1 || bucketcheckList[i].find(endKey)==bucketcheckList[i].end()){
      CkPrintf("Error: [%d] All the nodes not traversed by bucket no. %d\n",thisIndex,i);
      for (std::set<Tree::NodeKey>::iterator iter=bucketcheckList[i].begin(); iter != bucketcheckList[i].end(); iter++) {
	CkPrintf("       [%d] key %016llx\n",thisIndex,*iter);
      }
      break;
    }
    else { bucketcheckList[i].clear(); }
  }
}
#endif

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
    accelerationBox.grow(myParticles[i].treeAcceleration);
    if(!xdr_template(&xdrs, &(myParticles[i].treeAcceleration))) {
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

void TreePiece::outputAccASCII(const string& suffix, const CkCallback& cb) {
  if((thisIndex==0 && packed) || (thisIndex==0 && !packed && cnt==0)) {
    if(verbosity > 2)
      ckerr << "TreePiece " << thisIndex << ": Writing header for accelerations file" << endl;
    FILE* outfile = fopen((basefilename + "." + suffix).c_str(), "w");
		fprintf(outfile,"%d\n",(int) fh.numParticles);
    fclose(outfile);
  }
	
  if(verbosity > 3)
    ckerr << "TreePiece " << thisIndex << ": Writing my accelerations to disk" << endl;
	
  FILE* outfile = fopen((basefilename + "." + suffix).c_str(), "r+");
  fseek(outfile, 0, SEEK_END);
	
  for(unsigned int i = 1; i <= myNumParticles; ++i) {
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

  fclose(outfile);
  if((thisIndex==(int)numTreePieces-1 && packed)
     || (thisIndex==(int)numTreePieces-1 && !packed && cnt==3)) {
      cb.send();
      }
	
  if(thisIndex==(int)numTreePieces-1 && !packed && cnt<3)
      pieces[0].outputAccASCII(suffix, cb);
		
  if(thisIndex!=(int)numTreePieces-1)
    pieces[thisIndex + 1].outputAccASCII(suffix, cb);
	
}

void TreePiece::outputIOrderASCII(const string& suffix, const CkCallback& cb) {
  if(thisIndex==0) {
    if(verbosity > 2)
      ckerr << "TreePiece " << thisIndex << ": Writing header for iOrder file"
	    << endl;
    FILE* outfile = fopen((basefilename + "." + suffix).c_str(), "w");
    fprintf(outfile,"%d\n",(int) fh.numParticles);
    fclose(outfile);
  }
	
  if(verbosity > 3)
    ckerr << "TreePiece " << thisIndex << ": Writing iOrder to disk" << endl;
	
  FILE* outfile = fopen((basefilename + "." + suffix).c_str(), "r+");
  fseek(outfile, 0, SEEK_END);
  
  for(unsigned int i = 1; i <= myNumParticles; ++i) {
      if(fprintf(outfile,"%d\n", myParticles[i].iOrder) < 0) {
	  ckerr << "TreePiece " << thisIndex
		<< ": Error writing iOrder to disk, aborting" << endl;
	  CkAbort("IO Badness");
	  }
      }
  
  fclose(outfile);
  if(thisIndex==(int)numTreePieces-1) {
      cb.send();
      }
  else
      pieces[thisIndex + 1].outputIOrderASCII(suffix, cb);
  }

/********************************************************************/

void TreePiece::outputStatistics(Interval<unsigned int> macInterval, Interval<unsigned int> cellInterval, Interval<unsigned int> particleInterval, Interval<unsigned int> callsInterval, double totalmass, const CkCallback& cb) {

#if COSMO_STATS > 0
  if(verbosity > 1) {
    u_int64_t nodeInterRemoteTotal = 0;
    u_int64_t particleInterRemoteTotal = 0;
    for (int i=0; i<numChunks; ++i) {
      nodeInterRemoteTotal += nodeInterRemote[i];
      particleInterRemoteTotal += particleInterRemote[i];
    }
    ckerr << "TreePiece ";
    ckerr << thisIndex;
    ckerr << ": Statistics\nMy number of MAC checks: ";
    ckerr << myNumMACChecks << endl;
    ckerr << "My number of opened node: "
	  << nodesOpenedLocal << " local, " << nodesOpenedRemote << " remote." << endl;
    ckerr << "My number of particle-node interactions: "
	  << nodeInterLocal << " local, " << nodeInterRemoteTotal << " remote. Per particle: "
	  << (nodeInterLocal+nodeInterRemoteTotal)/(double) myNumParticles << endl;
    //	 << "\nCache cell interactions count: " << cachecellcount << endl;
    ckerr << "My number of particle-particle interactions: "
	  << particleInterLocal << " local, " << particleInterRemoteTotal
	  << " remote. Per Particle: "
	  << (particleInterLocal+particleInterRemoteTotal)/(double) myNumParticles << endl;
  }
#endif	

#if COSMO_DEBUG > 1
checkWalkCorrectness();
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

  /*
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
  */
#endif
	
  if(thisIndex != (int) numTreePieces - 1)
    pieces[thisIndex + 1].outputStatistics(macInterval, cellInterval, particleInterval, callsInterval, totalmass, cb);
  if(thisIndex == (int) numTreePieces - 1) cb.send();
}

/*
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
*/

/// @TODO Fix pup routine to handle correctly the tree
void TreePiece::pup(PUP::er& p) {
  CBase_TreePiece::pup(p);
  p | numTreePieces;
  p | callback;
  p | myNumParticles;
  if(p.isUnpacking()) {
    myParticles = new GravityParticle[myNumParticles + 2];
  }
  for(unsigned int i=0;i<myNumParticles+2;i++){
    p | myParticles[i];
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
    case Binary_ORB:
      root = new BinaryTreeNode(1, Tree::Boundary, 0, myNumParticles+1, 0);
      break;
    case Oct_Oct:
      //root = new OctTreeNode(1, Tree::Boundary, 0, myNumParticles+1, 0);
      break;
    default:
      CkAbort("We should have never reached here!");
    }
  }

  //PUP components for ORB decomposition
  p | chunkRootLevel;
  if(p.isUnpacking()){
    boxes = new OrientedBox<float>[chunkRootLevel+1];
    splitDims = new char[chunkRootLevel+1];
  }
  for(unsigned int i=0;i<chunkRootLevel;i++){
    p | boxes[i];
    p | splitDims[i];
  }

  p | theta;
  //p | myNumParticlesPending;
  p | prefetchWaiting;
  p | currentPrefetch;
  p | numBuckets;
  p | currentBucket;
  p | currentRemoteBucket;
#if COSMO_STATS > 0
  //p | myNumParticleInteractions;
  //p | myNumCellInteractions;
  p | myNumMACChecks;
  p | nodesOpenedLocal;
  p | nodesOpenedRemote;
  p | numOpenCriterionCalls;
  p | piecemass;
#endif
  // the counters do not need to be pupped!
  //p | nodeInterLocal;
  //p | nodeInterRemote;
  //p | particleInterLocal;
  //p | particleInterRemote;
  if (p.isUnpacking()) {
    particleInterRemote = NULL;
    nodeInterRemote = NULL;

    switch(domainDecomposition) {
    case SFC_dec:
      numPrefetchReq = 2;
      prefetchReq = new OrientedBox<double>[2];
      break;
    case Oct_dec:
    case ORB_dec:
      numPrefetchReq = 1;
      prefetchReq = new OrientedBox<double>[1];
      break;
    default:
      CmiAbort("Pupper has wrong domain decomposition type!\n");
    }
  }

  if(p.isUnpacking()){
    localCache = cacheManagerProxy.ckLocalBranch();
    dm = NULL;

    // reconstruct the data for prefetching
    /* OLD, moved to the cache and startIteration
    numChunks = root->getNumChunks(_numChunks);
    remainingChunk = new int[numChunks];
    root->getChunks(_numChunks, prefetchRoots);
    */
  }

  int notNull = (root==NULL)?0:1;
  p | notNull;
  if (notNull == 1) {
    ckout <<thisIndex<<" pupping tree"<<endl;
    p | (*root);
    if(p.isUnpacking()){
      //  nodeLookupTable[root->getKey()]=root;
      //}

      // reconstruct the nodeLookupTable and the bucketList
      reconstructNodeLookup(root);
    }
  }

  if (verbosity) {
    ckout << "TreePiece " << thisIndex << ": Getting PUP'd!";
    if (p.isSizing()) ckout << " size: " << ((PUP::sizer*)&p)->size();
    ckout << endl;
  }

  /*
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
  */
}

void TreePiece::reconstructNodeLookup(GenericTreeNode *node) {
  nodeLookupTable[node->getKey()] = node;
  node->particlePointer = &myParticles[node->firstParticle];
  if (node->getType() == Bucket) bucketList.push_back(node);
  GenericTreeNode *child;
  for (unsigned int i=0; i<node->numChildren(); ++i) {
    child = node->getChildren(i);
    if (child != NULL) reconstructNodeLookup(child);
  }
}

/*
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
*/

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
    nodeOwnership(node->getKey(), first, last);
    os << "NonLocal: Chare=" << node->remoteIndex << ", Owners=" << first << "-" << last;
    break;
  case NonLocalBucket:
    //os << "NonLocal: Chare=" << node->remoteIndex << "\\nRemote N under: " << (node->lastParticle - node->firstParticle + 1) << "\\nOwners: " << node->numOwners;
    nodeOwnership(node->getKey(), first, last);
    CkAssert(first == last);
    os << "NonLocalBucket: Chare=" << node->remoteIndex << ", Owner=" << first << ", Size=" << node->particleCount << "(" << node->firstParticle << "-" << node->lastParticle << ")";
    break;
  case Boundary:
    nodeOwnership(node->getKey(), first, last);
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

/// Print a graphviz version of a tree
void printTree(GenericTreeNode* node, ostream& os) {
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
    //nodeOwnership(node->getKey(), first, last);
    os << "NonLocal: Chare=" << node->remoteIndex; //<< ", Owners=" << first << "-" << last;
    break;
  case NonLocalBucket:
    //os << "NonLocal: Chare=" << node->remoteIndex << "\\nRemote N under: " << (node->lastParticle - node->firstParticle + 1) << "\\nOwners: " << node->numOwners;
    //nodeOwnership(node->getKey(), first, last);
    //CkAssert(first == last);
    os << "NonLocalBucket: Chare=" << node->remoteIndex << ", Owner=" << first << ", Size=" << node->particleCount; //<< "(" << node->firstParticle << "-" << node->lastParticle << ")";
    break;
  case Boundary:
    //nodeOwnership(node->getKey(), first, last);
    os << "Boundary: Totalsize=" << node->particleCount << ", Localsize=" << (node->lastParticle - node->firstParticle) << "(" << node->firstParticle << "-" << node->lastParticle;
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

/*
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
*/

CkReduction::reducerType TreePieceStatistics::sum;

void TreePiece::collectStatistics(CkCallback& cb) {
#if COSMO_STATS > 0
  u_int64_t nodeInterRemoteTotal = 0;
  u_int64_t particleInterRemoteTotal = 0;
  for (int i=0; i<numChunks; ++i) {
    nodeInterRemoteTotal += nodeInterRemote[i];
    particleInterRemoteTotal += particleInterRemote[i];
  }
  TreePieceStatistics tps(nodesOpenedLocal, nodesOpenedRemote, numOpenCriterionCalls,
      nodeInterLocal, nodeInterRemoteTotal, particleInterLocal, particleInterRemoteTotal);
  contribute(sizeof(TreePieceStatistics), &tps, TreePieceStatistics::sum, cb);
#else
  CkAbort("Invalid call, only valid if COSMO_STATS is defined");
#endif
}
