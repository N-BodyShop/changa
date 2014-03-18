/** @file TreePiece.cpp
 */

#include <cstdio>
#include <algorithm>
#include <fstream>
#include <assert.h>
// jetley
#include "limits.h"

//#include "ComlibManager.h"

#include "ParallelGravity.h"
#include "DataManager.h"
#include "Reductions.h"
// jetley
#include "MultistepLB.h"
#include "MultistepLB_notopo.h"
#include "Orb3dLB.h"
#include "Orb3dLB_notopo.h"
// jetley - refactoring
//#include "codes.h"
#include "Opt.h"
#include "Compute.h"
#include "TreeWalk.h"
//#include "State.h"

#include "Space.h"
#include "gravity.h"
#include "smooth.h"

#include "PETreeMerger.h"

#ifdef CELL
#include "spert_ppu.h"
#include "cell_typedef.h"
#endif

#ifdef CUDA
#ifdef CUDA_INSTRUMENT_WRS
#define GPU_INSTRUMENT_WRS
#include "wr.h"
#endif
#endif

/*
// uncomment when using cuda version of charm but running without GPUs
struct workRequest; 
void kernelSelect(workRequest *wr) {
}
*/

#ifdef PUSH_GRAVITY
#include "ckmulticast.h"
#endif

using namespace std;
using namespace SFC;
using namespace TreeStuff;
using namespace TypeHandling;

int TreeStuff::maxBucketSize;

#ifdef PUSH_GRAVITY
extern CkGroupID ckMulticastGrpId;
#endif

#ifdef CELL
int workRequestOut = 0;
CkVec<CellComputation> ewaldMessages;
#endif

//forward declaration
string getColor(GenericTreeNode*);

const char *typeString(NodeType type);

/*
 * set periodic information in all the TreePieces
 */
void TreePiece::setPeriodic(int nRepsPar, // Number of replicas in
					  // each direction
			    Vector3D<double> fPeriodPar, // Size of periodic box
			    int bEwaldPar,     // Use Ewald summation
			    double fEwCutPar,  // Cutoff on real summation
			    double dEwhCutPar, // Cutoff on Fourier summation
			    int bPeriodPar     // Periodic boundaries
			    )
{
    nReplicas = nRepsPar;
    fPeriod = fPeriodPar;
    bEwald = bEwaldPar;
    fEwCut  = fEwCutPar;
    dEwhCut = dEwhCutPar;
    bPeriodic = bPeriodPar;
    if(ewt == NULL) {
	ewt = new EWT[nMaxEwhLoop];
    }
}

// Scale velocities (needed to convert to canonical momenta for
// comoving coordinates.
void TreePiece::velScale(double dScale, const CkCallback& cb)
{
    for(unsigned int i = 0; i < myNumParticles; ++i) 
	{
	    myParticles[i+1].velocity *= dScale;
	    if(TYPETest(&myParticles[i+1], TYPE_GAS))
		myParticles[i+1].vPred() *= dScale;
	    }
    contribute(cb);
    }

/// After the bounding box has been found, we can assign keys to the particles
void TreePiece::assignKeys(CkReductionMsg* m) {
	if(m->getSize() != sizeof(OrientedBox<float>)) {
		ckerr << thisIndex << ": TreePiece: Fatal: Wrong size reduction message received!" << endl;
		CkAssert(0);
		callback.send(0);
		delete m;
		return;
	}

	boundingBox = *static_cast<OrientedBox<float> *>(m->getData());
	delete m;
	if(thisIndex == 0 && verbosity > 1)
		ckout << "TreePiece: Bounding box originally: "
		     << boundingBox << endl;
	//give particles keys, using bounding box to scale
	if((domainDecomposition!=ORB_dec)
	   && (domainDecomposition!=ORB_space_dec)){
	      // get longest axis
	      Vector3D<float> bsize = boundingBox.size();
	      float max = (bsize.x > bsize.y) ? bsize.x : bsize.y;
	      max = (max > bsize.z) ? max : bsize.z;
	      //
	      // Make the bounding box cubical.
	      //
	      Vector3D<float> bcenter = boundingBox.center();
	      const float fEps = 1.0 + 3.0e-7;  // slop to ensure keys fall
						// between 0 and 1.
	      bsize = Vector3D<float>(fEps*0.5*max);
	      boundingBox = OrientedBox<float>(bcenter-bsize, bcenter+bsize);
	      if(thisIndex == 0 && verbosity > 1)
		      ckout << "TreePiece: Bounding box now: "
			   << boundingBox << endl;

	      for(unsigned int i = 0; i < myNumParticles; ++i) {
		myParticles[i+1].key = generateKey(myParticles[i+1].position,
						   boundingBox);
	      }
	      sort(&myParticles[1], &myParticles[myNumParticles+1]);
	}

#if COSMO_DEBUG > 1
  char fout[100];
  sprintf(fout,"tree.%d.%d.before",thisIndex,iterationNo);
  ofstream ofs(fout);
  for (int i=1; i<=myNumParticles; ++i)
    ofs << keyBits(myParticles[i].key,KeyBits) << " " << myParticles[i].position[0] << " "
        << myParticles[i].position[1] << " " << myParticles[i].position[2] << endl;
  ofs.close();
#endif

	if(verbosity >= 5)
		cout << thisIndex << ": TreePiece: Assigned keys to all my particles" << endl;

/*
  if(myNumParticles > 0) 
    CkPrintf("[%d] %d (%llx) - %d (%llx)\n", thisIndex, 0, myParticles[0].key, myNumParticles+1, myParticles[myNumParticles+1].key);
    */
  contribute(callback);

}



/**************ORB Decomposition***************/

/// Three comparison routines used in sort and upper_bound
/// to order particles in each of three dimensions, respectively
bool comp_dim0(GravityParticle p1, GravityParticle p2) {
    return p1.position[0] < p2.position[0];
}
bool comp_dim1(GravityParticle p1, GravityParticle p2) {
    return p1.position[1] < p2.position[1];
}
bool comp_dim2(GravityParticle p1, GravityParticle p2) {
    return p1.position[2] < p2.position[2];
}

///Initialize stuff before doing ORB decomposition
void TreePiece::initORBPieces(const CkCallback& cb){

  OrientedBox<float> box = boundingBox;
  orbBoundaries.clear();
  orbBoundaries.push_back(myParticles+1);
  orbBoundaries.push_back(myParticles+myNumParticles+1);
  firstTime=true;

  phase=0;

  //Initialize function pointers
  compFuncPtr[0]= &comp_dim0;
  compFuncPtr[1]= &comp_dim1;
  compFuncPtr[2]= &comp_dim2;

	myBinCountsORB.clear();
	myBinCountsORB.push_back(myNumParticles);

  //Find out how many levels will be there before we go into the tree owned
  //completely by the TreePiece
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

/// Allocate memory for sorted particles.
/// @param cb callback after everything is sorted.
/// @param cback callback to perform now.
void TreePiece::initBeforeORBSend(unsigned int myCount,
				  unsigned int myCountGas,
				  unsigned int myCountStar,
				  const CkCallback& cb,
				  const CkCallback& cback){

  callback = cb;
  CkCallback nextCallback = cback;
  if(numTreePieces == 1) {
      myCount = myNumParticles;
      myCountGas = myNumSPH;
      myCountStar = myNumStar;
      }
  myExpectedCount = myCount;
  myExpectedCountSPH = myCountGas;
  myExpectedCountStar = myCountStar;

  mySortedParticles.clear();
  mySortedParticles.reserve(myExpectedCount);
  mySortedParticlesSPH.clear();
  mySortedParticlesSPH.reserve(myExpectedCountSPH);
  mySortedParticlesStar.clear();
  mySortedParticlesStar.reserve(myExpectedCountStar);

  contribute(nextCallback);
}

void TreePiece::sendORBParticles(){

  std::list<GravityParticle *>::iterator iter;
  std::list<GravityParticle *>::iterator iter2;

  int i=0;
  int nCounts = myBinCountsORB.size()/3; // Gas counts are in second
					 // half and star counts are
					 // in third half
  
  for(iter=orbBoundaries.begin();iter!=orbBoundaries.end();iter++,i++){
	iter2=iter;
	iter2++;
	if(iter2==orbBoundaries.end())
		break;
	if(myBinCountsORB[i]>0) {
	    extraSPHData *pGasOut = NULL;
	    if(myBinCountsORB[nCounts+i] > 0) {
		pGasOut = new extraSPHData[myBinCountsORB[nCounts+i]];
		if (verbosity>=3)
		  CkPrintf("me:%d to:%d nPart :%d, nGas:%d\n", thisIndex, i,
			   *iter2 - *iter, myBinCountsORB[nCounts+i]);
		int iGasOut = 0;
		for(GravityParticle *pPart = *iter; pPart < *iter2; pPart++) {
		    if(TYPETest(pPart, TYPE_GAS)) {
			pGasOut[iGasOut] = *(extraSPHData *)pPart->extraData;
			iGasOut++;
			}
		    }
		}
	    extraStarData *pStarOut = NULL;
	    if(myBinCountsORB[2*nCounts+i] > 0) {
		pStarOut = new extraStarData[myBinCountsORB[2*nCounts+i]];
		int iStarOut = 0;
		for(GravityParticle *pPart = *iter; pPart < *iter2; pPart++) {
		    if(pPart->isStar()) {
			pStarOut[iStarOut] = *(extraStarData *)pPart->extraData;
			iStarOut++;
			}
		    }
		}
	
	    if(i==thisIndex){
		acceptORBParticles(*iter,myBinCountsORB[i], pGasOut,
				   myBinCountsORB[nCounts+i], pStarOut,
				   myBinCountsORB[2*nCounts+i]);
		}
	    else{
		pieces[i].acceptORBParticles(*iter,myBinCountsORB[i], pGasOut,
					     myBinCountsORB[nCounts+i],
					     pStarOut,
					     myBinCountsORB[2*nCounts+i]);
		}
	    if(pGasOut != NULL)
		delete[] pGasOut;
	    if(pStarOut != NULL)
		delete[] pStarOut;
	    }
      }

  if(myExpectedCount > (int) myNumParticles){
    delete [] myParticles;
    nStore = (int)((myExpectedCount + 2)*(1.0 + dExtraStore));
    myParticles = new GravityParticle[nStore];
  }
  myNumParticles = myExpectedCount;

  if(myExpectedCountSPH > (int) myNumSPH){
    if(nStoreSPH > 0) delete [] mySPHParticles;
    nStoreSPH = (int)(myExpectedCountSPH*(1.0 + dExtraStore));
    mySPHParticles = new extraSPHData[nStoreSPH];
  }
  myNumSPH = myExpectedCountSPH;

  if(myExpectedCountStar > (int) myNumStar){
    delete [] myStarParticles;
    allocateStars();
  }
  myNumStar = myExpectedCountStar;

  if(myExpectedCount == 0) // No particles.  Make sure transfer is
			   // complete
      acceptORBParticles(NULL, 0, NULL, 0, NULL, 0);
}

/// Accept particles from other TreePieces once the sorting has finished
void TreePiece::acceptORBParticles(const GravityParticle* particles,
				   const int n,
				   const extraSPHData *pGas,
				   const int nGasIn,
				   const extraStarData *pStar,
				   const int nStarIn) {

  copy(particles, particles + n, back_inserter(mySortedParticles));
  copy(pGas, pGas + nGasIn, back_inserter(mySortedParticlesSPH));
  copy(pStar, pStar + nStarIn, back_inserter(mySortedParticlesStar));

  if(myExpectedCount == mySortedParticles.size()) {
      //I've got all my particles
    //Assigning keys to particles
    for(int i=0;i<myExpectedCount;i++){
      mySortedParticles[i].key = thisIndex;
    }
    copy(mySortedParticles.begin(), mySortedParticles.end(), &myParticles[1]);
    copy(mySortedParticlesSPH.begin(), mySortedParticlesSPH.end(),
	 &mySPHParticles[0]);
    copy(mySortedParticlesStar.begin(), mySortedParticlesStar.end(),
	 &myStarParticles[0]);
    // assign gas and star data pointers
    int iGas = 0;
    int iStar = 0;
    for(int i=0;i<myExpectedCount;i++){
	if(myParticles[i+1].isGas()) {
	    myParticles[i+1].extraData
		= (extraSPHData *)&mySPHParticles[iGas];
	    iGas++;
	    }
	else if(myParticles[i+1].isStar()) {
	    myParticles[i+1].extraData
		= (extraStarData *)&myStarParticles[iStar];
	    iStar++;
	    }
	else
	    myParticles[i+1].extraData = NULL;
	}
	  //signify completion with a reduction
    if(verbosity>1)
      ckout << thisIndex <<" contributing to accept particles"<<endl;
    deleteTree();
    contribute(callback);
    }
}

/// @brief Determine my boundaries at the end of ORB decomposition
void TreePiece::finalizeBoundaries(ORBSplittersMsg *splittersMsg){

  CkCallback cback = splittersMsg->cb;

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
    Vector3D<double> divide(0.0,0.0,0.0);
    divide[dimen] = splittersMsg->pos[i];
    dummy.position = divide;
    GravityParticle* divEnd = upper_bound(*iter,*iter2,dummy,compFuncPtr[dimen]);

    orbBoundaries.insert(iter2,divEnd);
    iter = iter2;
    iter2++;
	}

  firstTime = true;

  // First part is total particles, second part is gas counts, third
  // part is star counts
  myBinCountsORB.assign(6*splittersMsg->length,0);
  copy(tempBinCounts.begin(),tempBinCounts.end(),myBinCountsORB.begin());

  delete splittersMsg;
  contribute(cback);
}

/// @brief Evaluate particle counts for ORB decomposition
/// @param m A message containing splitting dimensions and splits, and
///          the callback to contribute
/// Counts the particles of this treepiece on each side of the
/// splits.  These counts are summed in a contribution to the
/// specified callback.
///
void TreePiece::evaluateParticleCounts(ORBSplittersMsg *splittersMsg)
{

  CkCallback& cback = splittersMsg->cb;

  // For each split, BinCounts has total lower, total higher.
  // The second half of the array has the counts for gas particles.
  // The third half of the array has the counts for star particles.
  tempBinCounts.assign(6*splittersMsg->length,0);

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
    //Evaluate the number of particles in each division

    GravityParticle dummy;
    GravityParticle* divStart = *iter;
    Vector3D<double> divide(0.0,0.0,0.0);
    divide[dimen] = splittersMsg->pos[i];
    dummy.position = divide;
    GravityParticle* divEnd = upper_bound(*iter,*iter2,dummy,compFuncPtr[dimen]);
    tempBinCounts[2*i] = divEnd - divStart;
    tempBinCounts[2*i + 1] = myBinCountsORB[i] - (divEnd - divStart);
    int nGasLow = 0;
    int nGasHigh = 0;
    int nStarLow = 0;
    int nStarHigh = 0;
    for(GravityParticle *pPart = divStart; pPart < divEnd; pPart++) {
	// Count gas
	if(TYPETest(pPart, TYPE_GAS))
	    nGasLow++;
	// Count stars
	if(TYPETest(pPart, TYPE_STAR))
	    nStarLow++;
	}
    for(GravityParticle *pPart = divEnd; pPart < *iter2; pPart++) {
	// Count gas
	if(TYPETest(pPart, TYPE_GAS))
	    nGasHigh++;
	// Count stars
	if(TYPETest(pPart, TYPE_STAR))
	    nStarHigh++;
	}
    tempBinCounts[2*splittersMsg->length + 2*i] = nGasLow;
    tempBinCounts[2*splittersMsg->length + 2*i + 1] = nGasHigh;
    tempBinCounts[4*splittersMsg->length + 2*i] = nStarLow;
    tempBinCounts[4*splittersMsg->length + 2*i + 1] = nStarHigh;
    iter++; iter2++;
    }

  if(firstTime)
    firstTime=false;
  contribute(6*splittersMsg->length*sizeof(int), &(*tempBinCounts.begin()), CkReduction::sum_int, cback);
  delete splittersMsg;
}

#ifdef REDUCTION_HELPER
void ReductionHelper::evaluateBoundaries(bool convertToLoad, int skipEvery, const CkCallback& cb) {
  DataManager *localDmObj = dMProxy.ckLocalBranch();
  evaluateBoundaries(convertToLoad, &(localDmObj->boundaryKeys.front()), localDmObj->boundaryKeys.size(), skipEvery, cb);

}

void ReductionHelper::evaluateBoundaries(bool convertToLoad, SFC::Key* keys, const int n, int skipEvery, const CkCallback& cb){

  splitters.assign(keys, keys + n);
  if(localTreePieces.presentTreePieces.size() == 0){
    int numBins = skipEvery ? n - (n-1)/(skipEvery+1) - 1 : n - 1;
    double *dummy = new double[numBins];
    for(int i = 0; i < numBins; i++) dummy[i] = 0;
    contribute(sizeof(double)*numBins, dummy, CkReduction::sum_double, cb);
    delete [] dummy;
    return;
  }

  for(int i = 0; i < localTreePieces.presentTreePieces.size(); i++){
    localTreePieces.presentTreePieces[i]->evaluateBoundaries(convertToLoad, keys, n, skipEvery, cb);
  }
}

void ReductionHelper::evaluateBoundaries(const CkBitVector &binsToSplit, const CkCallback& cb) {
  std::vector<SFC::Key> newSplitters;
  SFC::Key leftBound, rightBound;

  newSplitters.reserve(splitters.size() * 4);
  newSplitters.push_back(SFC::firstPossibleKey);
  for (int i = 0; i < binsToSplit.Length(); i++) {

    if (binsToSplit.Test(i) == true) {
      leftBound = splitters[i];
      rightBound = splitters[i + 1];
      if (newSplitters.back() != (rightBound | 7L) ) {
        if (newSplitters.back() != (leftBound | 7L)) {
          newSplitters.push_back(leftBound | 7L);
        }
        newSplitters.push_back((leftBound / 4 * 3 + rightBound / 4) | 7L);
        newSplitters.push_back((leftBound / 2 + rightBound / 2) | 7L);
        newSplitters.push_back((leftBound / 4 + rightBound / 4 * 3) | 7L);
        newSplitters.push_back(rightBound | 7L);
      }

    }
  }
  if (newSplitters.back() != lastPossibleKey) {
    newSplitters.push_back(lastPossibleKey);
  }
  evaluateBoundaries(true, &newSplitters[0], newSplitters.size(), 0, cb);
}

#endif

/// Determine my part of the sorting histograms by counting the number
/// of my particles in each bin.
/// This routine assumes the particles in key order.
/// The parameter skipEvery means that every "skipEvery" bins counted, one
/// must be skipped.  When skipEvery is set, the keys are in groups of
/// "skipEvery" size, and only splits within each group need to be
/// evaluated.  Hence the counts between the end of one group, and the
/// start of the next group are not evaluated.  This feature is used
/// by the Oct decomposition.
void TreePiece::evaluateBoundaries(bool convertToLoad, SFC::Key* keys, const int n, int skipEvery, const CkCallback& cb){

  if (!warmupFinished) {
    treePieceLoadTmp = (double) myNumParticles;
  }
  else if (treePieceLoadTmp == 0) {
    treePieceLoadTmp = getObjTime();
  }
#ifdef COSMO_EVENT
  double startTimer = CmiWallTimer();
#endif

  /*
  for(int i = 1; i <= myNumParticles; i++){
    CkPrintf("(%d) part %d key %llx\n", CkMyPe(), i, myParticles[i].key);
  }
  */

  int numBins = skipEvery ? n - (n-1)/(skipEvery+1) - 1 : n - 1;

  //this array will contain the number of particles I own in each bin
  
  double *myLoads; 
  int64_t *myCounts;

#ifdef REDUCTION_HELPER
  myCounts = new int64_t[numBins];
#else
  //myBinCounts.assign(numBins, 0);
  myBinCounts.resize(numBins);
  myCounts = myBinCounts.getVec();
#endif
  if (convertToLoad) {
    myLoads = new double[numBins];
    memset(myLoads, 0, numBins * sizeof(double));
  }

  memset(myCounts, 0, numBins*sizeof(int64_t));

  if (myNumParticles > 0) {    
    Key* endKeys = keys+n;
    GravityParticle *binBegin = &myParticles[1];
    GravityParticle *binEnd;
    GravityParticle dummy;
    //int binIter = 0;
    //vector<int>::iterator binIter = myBinCounts.begin();
    //vector<Key>::iterator keyIter = dm->boundaryKeys.begin();
    Key* keyIter = lower_bound(keys, keys+n, binBegin->key);
    int binIter = skipEvery ? (keyIter-keys) - (keyIter-keys-1) / (skipEvery+1) - 1: keyIter - keys - 1;
    int skip = skipEvery ? skipEvery - (keyIter-keys-1) % (skipEvery+1) : -1;
    if (binIter == -1) {
      dummy.key = keys[0];
      binBegin = upper_bound(binBegin, &myParticles[myNumParticles+1], dummy);
      keyIter++;
      binIter++;
      skip = skipEvery ? skipEvery : -1;
    }
    for( ; keyIter != endKeys; ++keyIter) {
      dummy.key = *keyIter;
      /// find the last place I could put this splitter key in
      /// my array of particles
      binEnd = upper_bound(binBegin, &myParticles[myNumParticles+1], dummy);
      /// this tells me the number of particles between the
      /// last two splitter keys
      if (skip != 0) {
        myCounts[binIter] = ((int64_t)(binEnd - binBegin));
        if (convertToLoad) {
          myLoads[binIter] = treePieceLoadTmp / myNumParticles * myCounts[binIter];

          if (havePhaseData(activeRung) && savedPhaseParticle[activeRung] != 0) {
            int activerungbincount = 0;
            for(GravityParticle *pPart = binBegin; pPart < binEnd; pPart++) {
              if (pPart->rung >= activeRung) {
                activerungbincount++;
              }
            }
            myLoads[binIter] = savedPhaseLoad[activeRung] * 
              (activerungbincount / (float) savedPhaseParticle[activeRung]);
          } else if (havePhaseData(0) && myNumParticles != 0) {
            myLoads[binIter] = savedPhaseLoad[0] * 
              (myCounts[binIter] / (float) myNumParticles);
          }
        }
        ++binIter;
        --skip;
      } else {
        skip = skipEvery;
      }
      if(&myParticles[myNumParticles+1] <= binEnd) break;
      binBegin = binEnd;
    }

#ifdef COSMO_EVENTS
    traceUserBracketEvent(boundaryEvaluationUE, startTimer, CmiWallTimer());
#endif
  }
  
  /*
  for(int i = 0; i < numBins; i++){
    CkPrintf("[%d] evaluate bin %d key %llx cnt %d\n", CkMyPe(), i, keys[i+1], myCounts[i]);
  }
  */

  //send my bin counts back in a reduction
#ifdef REDUCTION_HELPER
  /*
  for(int i = 0; i < numBins; i++){
    CkPrintf("[%d] TP::evaluateBoundaries bin %d key %llu count %d\n", thisIndex, i, keys[i], myCounts[i]);
  }
  */
  if (convertToLoad) {
    reductionHelperProxy.ckLocalBranch()->reduceBinLoads(numBins, myLoads, cb);
    delete [] myLoads;
  }
  else {
    treePieceLoadTmp = 0.0;
    reductionHelperProxy.ckLocalBranch()->reduceBinCounts(numBins, myCounts, cb);
  }
  delete [] myCounts;
#else
  if (!convertToLoad) {
    contribute(numBins * sizeof(int), myCounts, CkReduction::sum_int, cb);
  }
  else {
    contribute(numBins * sizeof(double), myLoads, CkReduction::sum_double, cb); 
    delete [] myLoads;
  }
#endif
}

/// Once final splitter keys have been decided, I need to give my
/// particles out to the TreePiece responsible for them

void TreePiece::unshuffleParticles(CkReductionMsg* m){
  double tpLoad;

  warmupFinished = true;

  if (dm == NULL) {
    dm = (DataManager*)CkLocalNodeBranch(dataManagerID);
  }

  tpLoad = getObjTime();
  populateSavedPhaseData(prevLARung, tpLoad, treePieceActivePartsTmp);
  callback = *static_cast<CkCallback *>(m->getData());

  // CkPrintf("[%d] myplace %d\n", thisIndex, myPlace);

  /*
  if(myNumParticles > 0){
    for(int i = 0; i <= myNumParticles+1; i++){
      CkPrintf("[%d] part %d key %llx\n", thisIndex, i, myParticles[i].key);
    }
  }
  */

  if (myNumParticles == 0) {
    incomingParticlesSelf = true;
    acceptSortedParticles(NULL);
    delete m;
    return;
  }

  GravityParticle *binBegin = &myParticles[1];
  vector<Key>::iterator iter =
    lower_bound(dm->boundaryKeys.begin(), dm->boundaryKeys.end(),
                binBegin->key);
  vector<Key>::const_iterator endKeys = dm->boundaryKeys.end();
  int offset = iter - dm->boundaryKeys.begin() - 1;
  vector<int>::iterator responsibleIter = dm->responsibleIndex.begin() + offset;

  GravityParticle *binEnd;
  GravityParticle dummy;
  for( ; iter != endKeys; ++iter, ++responsibleIter) {
    dummy.key = *iter;
    //find particles between this and the last key
    binEnd = upper_bound(binBegin, &myParticles[myNumParticles+1],
        dummy);
    // If I have any particles in this bin, send them to
    // the responsible TreePiece
    int nPartOut = binEnd - binBegin;
    int saved_phase_len = savedPhaseLoad.size();

    if(nPartOut > 0) {
      int nGasOut = 0;
      int nStarOut = 0;
      for(GravityParticle *pPart = binBegin; pPart < binEnd;
          pPart++) {
        if(pPart->isGas())
          nGasOut++;
        if(pPart->isStar())
          nStarOut++;
      }

      ParticleShuffleMsg *shuffleMsg
        = new (saved_phase_len, saved_phase_len, nPartOut, nGasOut, nStarOut)
        ParticleShuffleMsg(saved_phase_len, nPartOut, nGasOut, nStarOut, 0.0);
      memset(shuffleMsg->parts_per_phase, 0, saved_phase_len*sizeof(unsigned int));

      // Calculate the number of particles leaving the treepiece per phase
      for(GravityParticle *pPart = binBegin; pPart < binEnd; pPart++) {
        for(int i = 0; i < saved_phase_len; i++) {
          if (pPart->rung >= i) {
            shuffleMsg->parts_per_phase[i] = shuffleMsg->parts_per_phase[i] + 1;
          }
        }
      }

      shuffleMsg->load = tpLoad * nPartOut / myNumParticles;
      memset(shuffleMsg->loads, shuffleMsg->load, saved_phase_len*sizeof(double));

      // Calculate the partial load per phase
      for (int i = 0; i < saved_phase_len; i++) {
        if (havePhaseData(i) && savedPhaseParticle[i] != 0) {
          shuffleMsg->loads[i] = savedPhaseLoad[i] *
            (shuffleMsg->parts_per_phase[i] / (float) savedPhaseParticle[i]);
        } else if (havePhaseData(0) && myNumParticles != 0) {
          shuffleMsg->loads[i] = savedPhaseLoad[0] *
            (shuffleMsg->parts_per_phase[i] / (float) myNumParticles);
        }
      }

      if (verbosity>=3)
        CkPrintf("me:%d to:%d nPart :%d, nGas:%d, nStar: %d\n",
            thisIndex, *responsibleIter,nPartOut, nGasOut,
            nStarOut);

      int iGasOut = 0;
      int iStarOut = 0;
      GravityParticle *pPartOut = shuffleMsg->particles;
      for(GravityParticle *pPart = binBegin; pPart < binEnd;
          pPart++, pPartOut++) {
        *pPartOut = *pPart;
        if(pPart->isGas()) {
          shuffleMsg->pGas[iGasOut]
            = *(extraSPHData *)pPart->extraData;
          iGasOut++;
        }
        if(pPart->isStar()) {
          shuffleMsg->pStar[iStarOut]
            = *(extraStarData *)pPart->extraData;
          iStarOut++;
        }
      }

      if(*responsibleIter == thisIndex) {
        if (verbosity > 1)
          CkPrintf("TreePiece %d: keeping %d / %d particles: %d\n",
              thisIndex, nPartOut, myNumParticles,
              nPartOut*10000/myNumParticles);
        acceptSortedParticles(shuffleMsg);
      }
      else {
        treeProxy[*responsibleIter].acceptSortedParticles(shuffleMsg);
      }
    }
    if(&myParticles[myNumParticles + 1] <= binEnd)
      break;
    binBegin = binEnd;
  }

  incomingParticlesSelf = true;
  acceptSortedParticles(NULL);

  delete m;
}

/// Accept particles from other TreePieces once the sorting has finished
void TreePiece::acceptSortedParticles(ParticleShuffleMsg *shuffleMsg) {
  //Need to get the place here again.  Getting the place in
  //unshuffleParticles and using it here results in a race condition.
  if (dm == NULL)
    dm = (DataManager*)CkLocalNodeBranch(dataManagerID);
  myPlace = find(dm->responsibleIndex.begin(), dm->responsibleIndex.end(),
      thisIndex) - dm->responsibleIndex.begin();
  if (myPlace == dm->responsibleIndex.size()) myPlace = -2;

  // The following assert does not work anymore when TreePieces can
  //have 0 particles assigned
  //assert(myPlace >= 0 && myPlace < dm->particleCounts.size());
  if (myPlace == -2 || dm->particleCounts[myPlace] == 0) {
    // Special case where no particle is assigned to this TreePiece
    if (myNumParticles > 0){
      delete[] myParticles;
      myParticles = NULL;
    }
    myNumParticles = 0;
    nStore = 0;
    if (nStoreSPH > 0){
      delete[] mySPHParticles;
      mySPHParticles = NULL;
    }
    myNumSPH = 0;
    nStoreSPH = 0;
    if (nStoreStar > 0){
      delete[] myStarParticles;
      myStarParticles = NULL;
    }
    myNumStar = 0;
    nStoreStar = 0;
    incomingParticlesSelf = false;
    incomingParticlesMsg.clear();
    if(verbosity>1) ckout << thisIndex <<" no particles assigned"<<endl;

    deleteTree();
    // We better not have a message with particles for us
    CkAssert(shuffleMsg == NULL);
    contribute(callback);
    return;
  }

  if(shuffleMsg != NULL) {
    incomingParticlesMsg.push_back(shuffleMsg);
    incomingParticlesArrived += shuffleMsg->n;
    treePieceLoadTmp += shuffleMsg->load; 
    savePhaseData(savedPhaseLoadTmp, savedPhaseParticleTmp, shuffleMsg->loads,
      shuffleMsg->parts_per_phase, shuffleMsg->nloads);
  }

  if (verbosity>=3)
    ckout << thisIndex <<" waiting for "
      << dm->particleCounts[myPlace]-incomingParticlesArrived
      << " particles ("<<dm->particleCounts[myPlace]<<"-"
      << incomingParticlesArrived<<")"
      << (incomingParticlesSelf?" self":"")<<endl;


  if(dm->particleCounts[myPlace] == incomingParticlesArrived && incomingParticlesSelf) {
    //I've got all my particles
    if (myNumParticles > 0) delete[] myParticles;

    nStore = (int)((dm->particleCounts[myPlace] + 2)*(1.0 + dExtraStore));
    myParticles = new GravityParticle[nStore];
    myNumParticles = dm->particleCounts[myPlace];
    incomingParticlesArrived = 0;
    incomingParticlesSelf = false;
    treePieceLoad = treePieceLoadTmp;
    treePieceLoadTmp = 0.0;

    savedPhaseLoad.swap(savedPhaseLoadTmp);
    savedPhaseParticle.swap(savedPhaseParticleTmp);
    savedPhaseLoadTmp.clear();
    savedPhaseParticleTmp.clear();

    int nSPH = 0;
    int nStar = 0;
    int iMsg;
    for(iMsg = 0; iMsg < incomingParticlesMsg.size(); iMsg++) {
      nSPH += incomingParticlesMsg[iMsg]->nSPH;
      nStar += incomingParticlesMsg[iMsg]->nStar;
    }
    if (nStoreSPH > 0) delete[] mySPHParticles;
    myNumSPH = nSPH;
    nStoreSPH = (int)(myNumSPH*(1.0 + dExtraStore));
    if(nStoreSPH > 0) mySPHParticles = new extraSPHData[nStoreSPH];
    else mySPHParticles = NULL;

    if (nStoreStar > 0) delete[] myStarParticles;
    myNumStar = nStar;
    allocateStars();

    int nPart = 0;
    nSPH = 0;
    nStar = 0;
    for(iMsg = 0; iMsg < incomingParticlesMsg.size(); iMsg++) {
      memcpy(&myParticles[nPart+1], incomingParticlesMsg[iMsg]->particles,
          incomingParticlesMsg[iMsg]->n*sizeof(GravityParticle));
      nPart += incomingParticlesMsg[iMsg]->n;
      memcpy(&mySPHParticles[nSPH], incomingParticlesMsg[iMsg]->pGas,
          incomingParticlesMsg[iMsg]->nSPH*sizeof(extraSPHData));
      nSPH += incomingParticlesMsg[iMsg]->nSPH;
      memcpy(&myStarParticles[nStar], incomingParticlesMsg[iMsg]->pStar,
          incomingParticlesMsg[iMsg]->nStar*sizeof(extraStarData));
      nStar += incomingParticlesMsg[iMsg]->nStar;
      delete incomingParticlesMsg[iMsg];
    }

    incomingParticlesMsg.clear();
    // assign gas data pointers and determine centroid
    int iGas = 0;
    int iStar = 0;
    Vector3D<double> vCenter(0.0, 0.0, 0.0);
    for(int iPart = 0; iPart < myNumParticles; iPart++) {
      vCenter += myParticles[iPart+1].position;
      if(myParticles[iPart+1].isGas()) {
        myParticles[iPart+1].extraData
          = (extraSPHData *)&mySPHParticles[iGas];
        iGas++;
      }
      else if(myParticles[iPart+1].isStar()) {
        myParticles[iPart+1].extraData
          = (extraStarData *)&myStarParticles[iStar];
        iStar++;
      }
      else
        myParticles[iPart+1].extraData = NULL;
    }

    sort(myParticles+1, myParticles+myNumParticles+1);
    savedCentroid = vCenter/(double)myNumParticles;
    //signify completion with a reduction
    if(verbosity>1) ckout << thisIndex <<" contributing to accept particles"
      <<endl;

    deleteTree();
    //CkPrintf("[%d] accepted %d particles\n", thisIndex, myNumParticles);
    contribute(callback);
  }
}

void TreePiece::savePhaseData(std::vector<double> &loads, std::vector<unsigned int>
  &parts_per_phase, double* shuffleloads, unsigned int *shuffleparts, int
  shufflelen) {
  int len = 0;
  int num_additional;

  for (int i = 0; i < shufflelen; i++) {
    len = loads.size();
    if (i >= len) {
      num_additional = i - len + 1;
      while (num_additional > 0) {
        loads.push_back(-1.0);
        parts_per_phase.push_back(0);
        num_additional--;
      }
    }
    if (loads[i] == -1) {
      loads[i] = 0.0;
    }
    loads[i] += shuffleloads[i];
    parts_per_phase[i] += shuffleparts[i];
  }
}

void TreePiece::populateSavedPhaseData(int phase, double tp_load,
    unsigned int activeparts) {
  int len = savedPhaseLoad.size();
  int num_additional;
  if (phase == -1) {
    phase = 0;
    activeparts = myNumParticles;
    //return;
  }

  if (phase > len-1) {
    num_additional = phase - len + 1;
    while (num_additional > 0) {
      savedPhaseLoad.push_back(-1.0);
      savedPhaseParticle.push_back(0);
      num_additional--;
    }
    len = savedPhaseLoad.size();
  }
  savedPhaseLoad[phase] = tp_load;
  savedPhaseParticle[phase] = activeparts;
}

bool TreePiece::havePhaseData(int phase) {
  return (savedPhaseLoad.size() > phase && savedPhaseLoad[phase] > -0.5);
}

#if 0
void TreePiece::checkin(){
  if(myDecomposer == NULL){
    myDecomposer = decomposerProxy.ckLocalBranch();
  }
  CkPrintf("[%d] checkin\n", thisIndex);
  myDecomposer->checkin();
}

void Decomposer::checkin(){
  numTreePiecesCheckedIn++;
  // +1 for self checkin(), since Decomposer::unshuffleParticles
  // (in which myNumTreePieces is set) may be called
  // after all local TreePieces have called Decomposer::checkin()
  // through TreePiece::acceptSortedParticles();
  CkPrintf("decomposer %d checked in %d/%d\n", CkMyPe(),
  numTreePiecesCheckedIn, myNumTreePieces);
  if(numTreePiecesCheckedIn == myNumTreePieces){
    numTreePiecesCheckedIn = 0;
    myNumTreePieces = -1;

    // return control to mainchare
    //contribute(callback);
  }
}
#endif

// Sum energies for diagnostics
void TreePiece::calcEnergy(const CkCallback& cb) {
    double dEnergy[7]; // 0 -> kinetic; 1 -> virial ; 2 -> potential;
		       // 3-5 -> L; 6 -> thermal
    Vector3D<double> L;

    dEnergy[0] = 0.0;
    dEnergy[1] = 0.0;
    dEnergy[2] = 0.0;
    dEnergy[6] = 0.0;
    for(unsigned int i = 0; i < myNumParticles; ++i) {
	GravityParticle *p = &myParticles[i+1];

	dEnergy[0] += p->mass*p->velocity.lengthSquared();
	dEnergy[1] += p->mass*dot(p->treeAcceleration, p->position);
	dEnergy[2] += p->mass*p->potential;
	L += p->mass*cross(p->position, p->velocity);
	if (TYPETest(p, TYPE_GAS))
	    dEnergy[6] += p->mass*p->u();
	}
    dEnergy[0] *= 0.5;
    dEnergy[2] *= 0.5;
    dEnergy[3] = L.x;
    dEnergy[4] = L.y;
    dEnergy[5] = L.z;

    contribute(7*sizeof(double), dEnergy, CkReduction::sum_double, cb);
}

void TreePiece::kick(int iKickRung, double dDelta[MAXRUNG+1],
		     int bClosing, // Are we at the end of a timestep
		     int bNeedVPred, // do we need to update vpred
		     int bGasIsothermal, // Isothermal EOS
		     double duDelta[MAXRUNG+1], // dts for energy
		     const CkCallback& cb) {
  // LBTurnInstrumentOff();
  for(unsigned int i = 1; i <= myNumParticles; ++i) {
      GravityParticle *p = &myParticles[i];
      if(p->rung >= iKickRung) {
	  if(bNeedVPred && TYPETest(p, TYPE_GAS)) {
	      if(bClosing) { // update predicted quantities to end of step
		  p->vPred() = p->velocity
		      + dDelta[p->rung]*p->treeAcceleration;
		  if(!bGasIsothermal) {
#ifndef COOLING_NONE
		      p->u() = p->u() + p->uDot()*duDelta[p->rung];
		      if (p->u() < 0) {
			  double uold = p->u() - p->uDot()*duDelta[p->rung];
			  p->u() = uold*exp(p->uDot()*duDelta[p->rung]/uold);
			  }
#else /* COOLING_NONE */
		      p->u() += p->PdV()*duDelta[p->rung];
#endif /* COOLING_NONE */
		      p->uPred() = p->u();
		      }
		  }
	      else {	// predicted quantities are at the beginning
			// of step
		  p->vPred() = p->velocity;
		  if(!bGasIsothermal) {
		      p->uPred() = p->u();
#ifndef COOLING_NONE
		      p->u() += p->uDot()*duDelta[p->rung];
		      if (p->u() < 0) {
			  double uold = p->u() - p->uDot()*duDelta[p->rung];
			  p->u() = uold*exp(p->uDot()*duDelta[p->rung]/uold);
			  }
#else /* COOLING_NONE */
		      p->u() += p->PdV()*duDelta[p->rung];
#endif /* COOLING_NONE */
		      }
		  }
	      CkAssert(p->u() >= 0.0);
	      CkAssert(p->uPred() >= 0.0);
	      }
	  p->velocity += dDelta[p->rung]*p->treeAcceleration;
	  }
      }
  contribute(cb);
}

void TreePiece::initAccel(int iKickRung, const CkCallback& cb) 
{
    for(unsigned int i = 1; i <= myNumParticles; ++i) {
	if(myParticles[i].rung >= iKickRung) {
	    myParticles[i].treeAcceleration = 0;
	    myParticles[i].potential = 0;
	    myParticles[i].dtGrav = 0;
	    }
	}

    contribute(cb);
    }

/**
 * Adjust timesteps of active particles.
 * @param iKickRung The rung we are on.
 * @param bEpsAccStep Use sqrt(eps/acc) timestepping
 * @param bGravStep Use sqrt(r^3/GM) timestepping
 * @param bSphStep Use Courant condition
 * @param bViscosityLimitdt Use viscosity in Courant condition
 * @param dEta Factor to use in determing timestep
 * @param dEtaCourant Courant factor to use in determing timestep
 * @param dEtauDot Factor to use in uDot based timestep
 * @param dDelta Base timestep
 * @param dAccFac Acceleration scaling for cosmology
 * @param dCosmoFac Cosmo scaling for Courant
 * @param dhMinOverSoft minimum smoothing parameter.
 * @param bDoGas We are calculating gas forces.
 * @param cb Callback function reduces currrent maximum rung
 */
void TreePiece::adjust(int iKickRung, int bEpsAccStep, int bGravStep,
		       int bSphStep, int bViscosityLimitdt,
		       double dEta, double dEtaCourant, double dEtauDot,
		       double dDelta, double dAccFac,
		       double dCosmoFac, double dhMinOverSoft,
		       int bDoGas,
		       const CkCallback& cb) {
  int iCurrMaxRung = 0;
  int nMaxRung = 0;  // number of particles in maximum rung
  
  for(unsigned int i = 1; i <= myNumParticles; ++i) {
    GravityParticle *p = &myParticles[i];
    if(p->rung >= iKickRung) {
      CkAssert(p->soft > 0.0);
      double dTIdeal = dDelta;
      if(bEpsAccStep) {
	  double acc = dAccFac*p->treeAcceleration.length();
	  double dt;
#ifdef EPSACCH
	  // Set gas timestep based on SPH smoothing.
	  if(bDoGas && p->isGas()) {
	      dt = dEta*sqrt(0.5*p->fBall/acc);
	      }
#else
	  // Set gas timestep based on gravity softening and SPH smoothing.
	  if(bDoGas && p->isGas() && dhMinOverSoft < 1
	     && p->fBall < 2.0*p->soft) {
	      if(p->fBall > 2.0*dhMinOverSoft*p->soft)
		  dt = dEta*sqrt(0.5*p->fBall/acc);
	      else
		  dt = dEta*sqrt(dhMinOverSoft*p->soft/acc);
	      }
#endif
	  else
	      dt = dEta*sqrt(p->soft/acc);
	  if(dt < dTIdeal)
	      dTIdeal = dt;
	  }
      if(bGravStep) {
	  double dt = dEta/sqrt(dAccFac*p->dtGrav);
	  if(dt < dTIdeal)
	      dTIdeal = dt;
	  }
      if(bSphStep && TYPETest(p, TYPE_GAS)) {
	  double dt;
	  double ph = sqrt(0.25*p->fBall*p->fBall);
	  if (p->mumax() > 0.0) {
	      if (bViscosityLimitdt) 
		  dt = dEtaCourant*dCosmoFac*(ph /(p->c() + 0.6*(p->c() + 2*p->BalsaraSwitch()*p->mumax())));
	      else
		  dt = dEtaCourant*dCosmoFac*(ph/(p->c() + 0.6*(p->c() + 2*p->mumax())));
	      }
	  else
	      dt = dEtaCourant*dCosmoFac*(ph/(1.6*p->c()));
	  if(dt < dTIdeal)
	      dTIdeal = dt;

	  if (dEtauDot > 0.0 && p->PdV() < 0.0) { /* Prevent rapid adiabatic cooling */
	      assert(p->u() > 0.0);
	      // Use P/rho as internal energy estimate since "u" may
	      // be driven to 0 with cooling.
	      dt = dEtauDot*p->fDensity*p->PoverRho2()/fabs(p->PdV());
	      if (dt < dTIdeal) 
		  dTIdeal = dt;
	      }
	  }

      int iNewRung = DtToRung(dDelta, dTIdeal);
      if(iNewRung > MAXRUNG) {
	CkError("dt: %g, soft: %g, accel: %g\n", dTIdeal, p->soft,
		p->treeAcceleration.length());
	CkAbort("Timestep too small");
	}
      if(iNewRung < iKickRung) iNewRung = iKickRung;
      if(iNewRung > iCurrMaxRung) {
	  iCurrMaxRung = iNewRung;
	  nMaxRung = 1;
	  }
      else if(iNewRung == iCurrMaxRung)
	  nMaxRung++;
      myParticles[i].rung = iNewRung;
#ifdef NEED_DT
      myParticles[i].dt = dTIdeal;
#endif
    }
  }
  // Pack into array for reduction
  int newcount[2];
  newcount[0] = iCurrMaxRung;
  newcount[1] = nMaxRung;
  contribute(2*sizeof(int), newcount, max_count, cb);
}

void TreePiece::truncateRung(int iCurrMaxRung, const CkCallback& cb) {
    for(unsigned int i = 1; i <= myNumParticles; ++i) {
	GravityParticle *p = &myParticles[i];
	if(p->rung > iCurrMaxRung)
	    p->rung--;
	CkAssert(p->rung <= iCurrMaxRung);
	}
    contribute(cb);
    }

void TreePiece::rungStats(const CkCallback& cb) {
  int64_t nInRung[MAXRUNG+1];

  for(int iRung = 0; iRung <= MAXRUNG; iRung++) nInRung[iRung] = 0;
  for(unsigned int i = 1; i <= myNumParticles; ++i) {
    nInRung[myParticles[i].rung]++;
  }
  contribute((MAXRUNG+1)*sizeof(int64_t), nInRung, CkReduction::sum_long, cb);
}

void TreePiece::countActive(int activeRung, const CkCallback& cb) {
  int64_t nActive[2];

  nActive[0] = nActive[1] = 0;
  for(unsigned int i = 1; i <= myNumParticles; ++i) {
      if(myParticles[i].rung >= activeRung) {
	  nActive[0]++;
	  if(TYPETest(&myParticles[i], TYPE_GAS)) {
	      nActive[1]++;
	      }
	  }
      }
  contribute(2*sizeof(int64_t), nActive, CkReduction::sum_long, cb);
}

/// @brief assign domain number to each particle for diagnostic
void TreePiece::assignDomain(const CkCallback &cb)
{
    for(unsigned int i = 1; i <= myNumParticles; ++i) {
        myParticles[i].interMass = thisIndex;
	}
    contribute(cb);
    }

void TreePiece::drift(double dDelta,  // time step in x containing
				      // cosmo scaling
		      int bNeedVpred, // Update predicted velocities
		      int bGasIsothermal, // Isothermal EOS
		      double dvDelta, // time step in v containing
				      // cosmo scaling
		      double duDelta, // time step for internal energy
		      int nGrowMass,  // GrowMass particles are locked
				      // in place
		      bool buildTree, // is a treebuild happening before the
				      // next drift?
		      const CkCallback& cb) {
  callback = cb;		// called by assignKeys()
  deleteTree();

  if(bucketReqs != NULL) {
    delete[] bucketReqs;
    bucketReqs = NULL;
  }

  boundingBox.reset();
  int bInBox = 1;

  for(unsigned int i = 1; i <= myNumParticles; ++i) {
      GravityParticle *p = &myParticles[i];
      if (p->iOrder >= nGrowMass)
	  p->position += dDelta*p->velocity;
      if(bPeriodic) {
        for(int j = 0; j < 3; j++) {
          if(p->position[j] >= 0.5*fPeriod[j]){
            p->position[j] -= fPeriod[j];
          }
          if(p->position[j] < -0.5*fPeriod[j]){
            p->position[j] += fPeriod[j];
          }

          bool a = (p->position[j] >= -0.5*fPeriod[j]);
          bool b = (p->position[j] < 0.5*fPeriod[j]);

          //CkPrintf("[%d] particle %d test %d: %d,%d\n", thisIndex, i, j, a, b);
          // Sanity Checks
          bInBox = bInBox && a; 
          bInBox = bInBox && b;
        }
        if(!bInBox){
          CkPrintf("[%d] p: %f,%f,%f\n", thisIndex, p->position[0], p->position[1], p->position[2]);
          CkAbort("binbox failed\n");
        }
      }
      boundingBox.grow(p->position);
      if(bNeedVpred && TYPETest(p, TYPE_GAS)) {
	  p->vPred() += dvDelta*p->treeAcceleration;
	  if(!bGasIsothermal) {
#ifndef COOLING_NONE
	      p->uPred() += p->uDot()*duDelta;
	      if (p->uPred() < 0) {
		  // Backout the update to upred
		  double uold = p->uPred() - p->uDot()*duDelta;
		  // uold could be negative because of round-off
		  // error.  If this is the case then uDot*Delta/u is
		  // large, the final uPred will be zero.
		  if(uold <= 0.0) p->uPred() = 0.0;
		  // Cooling rate is large: use an exponential decay
		  // of timescale u/uDot.
		  else p->uPred() = uold*exp(p->uDot()*duDelta/uold);
		  }
#else
	      p->uPred() += p->PdV()*duDelta;
	      if (p->uPred() < 0) {
		  double uold = p->uPred() - p->PdV()*duDelta;
		  p->uPred() = uold*exp(p->PdV()*duDelta/uold);
		  }
#endif
              CkAssert(p->uPred() >= 0.0);
	      }
	  }
      }
  CkAssert(bInBox);
  if(!bInBox){
    CkAbort("binbox2 failed\n");
  }
  if(buildTree)
    contribute(sizeof(OrientedBox<float>), &boundingBox,
      	       growOrientedBox_float,
      	       CkCallback(CkIndex_TreePiece::assignKeys(0), pieces));
  else
    contribute(cb);
}

void TreePiece::newParticle(GravityParticle *p)
{
    if(myNumParticles + 2 >= nStore) {
	CkAbort("No room for new particle: increase dExtraStore");
	}
    // Move Boundary particle
    myParticles[myNumParticles+2] = myParticles[myNumParticles+1];
    myNumParticles++;
    myParticles[myNumParticles] = *p;
    myParticles[myNumParticles].iOrder = -1;
    if(p->isGas()) {
	if(myNumSPH >= nStoreSPH)
	    CkAbort("No room for new SPH particle: increase dExtraStore");
	mySPHParticles[myNumSPH] = *((extraSPHData *) p->extraData);
	myParticles[myNumParticles].extraData = &mySPHParticles[myNumSPH];
	myNumSPH++;
	}
    if(p->isStar()) {
	if(myNumStar >= nStoreStar)
	    CkAbort("No room for new Star particle: increase dExtraStore");
	myStarParticles[myNumStar] = *((extraStarData *) p->extraData);
	myParticles[myNumParticles].extraData = &myStarParticles[myNumStar];
	myNumStar++;
	}
    }

/**
 * Count add/deleted particles, and compact main particle storage.
 * Extradata storage will be reclaimed on the next domain decompose.
 * This contributes to a "set" reduction, which the main chare will
 * iterate through to assign new iOrders.
 */
void TreePiece::colNParts(const CkCallback &cb)
{
    CountSetPart counts;
    counts.index = thisIndex;
    counts.nAddGas = 0;
    counts.nAddDark = 0;
    counts.nAddStar = 0;
    counts.nDelGas = 0;
    counts.nDelDark = 0;
    counts.nDelStar = 0;
    unsigned int i, j;
    int newNPart = myNumParticles; // the number of particles may change.
    
    for(i = 1, j = 1; i <= myNumParticles; ++i) {
	if(j < i)
	    myParticles[j] = myParticles[i];
	GravityParticle *p = &myParticles[i];
	if(p->iOrder == -1) { // A new Particle
	    j++;
	    if (p->isGas())
		counts.nAddGas++;
	    else if(p->isStar())
		counts.nAddStar++;
	    else if(p->isDark())
		counts.nAddDark++;
	    else
		CkAbort("Bad Particle type");
	    }
	else if(TYPETest(p, TYPE_DELETED) ) {
	    newNPart--;
	    if (p->isGas())
		counts.nDelGas++;
	    else if(p->isStar())
		counts.nDelStar++;
	    else if(p->isDark())
		counts.nDelDark++;
	    else
		CkAbort("Bad Particle type");
	    }
	else {
	    j++;
	    }
	}
    if(j < i)  // move boundary particle if needed
	myParticles[j] = myParticles[i];

    myNumParticles = newNPart;
    contribute(sizeof(counts), &counts, CkReduction::concat, cb);
    }

/**
 * Assign iOrders to recently added particles.
 * Also insure keys are OK
 */
void TreePiece::newOrder(int64_t nStartSPH, int64_t nStartDark,
			  int64_t nStartStar, const CkCallback &cb) 
{
    unsigned int i;
    boundingBox.reset();
    for(i = 1; i <= myNumParticles; ++i) {
	GravityParticle *p = &myParticles[i];
	boundingBox.grow(p->position);
	if(p->iOrder == -1) {
	    if (p->isGas()) 
		p->iOrder = nStartSPH++;
	    else if (p->isDark()) 
		p->iOrder = nStartDark++;
	    else {
		/* Also record iOrder in the starLog table. */
		CmiLock(dm->lockStarLog);
		dm->starLog->seTab[dm->starLog->nOrdered].iOrdStar = nStartStar;
		dm->starLog->nOrdered++;
		CmiUnlock(dm->lockStarLog);
		p->iOrder = nStartStar++;
		}
	    }
	}
    callback = cb;		// called by assignKeys()
    // get the new particles into key order
    contribute(sizeof(OrientedBox<float>), &boundingBox,
	       growOrientedBox_float,
	       CkCallback(CkIndex_TreePiece::assignKeys(0), pieces));
    }
    
/**
 * Update total particle numbers.
 */
void TreePiece::setNParts(int64_t _nTotalSPH, int64_t _nTotalDark,
			  int64_t _nTotalStar, const CkCallback &cb) 
{
    nTotalSPH = _nTotalSPH;
    nTotalDark = _nTotalDark;
    nTotalStar = _nTotalStar;
    nTotalParticles = nTotalSPH + nTotalDark + nTotalStar;
    contribute(cb);
    }

/// @param dSoft gravitational softening
/// @param cb callback
void TreePiece::setSoft(const double dSoft, const CkCallback& cb) {
  for(unsigned int i = 1; i <= myNumParticles; ++i) {
#ifdef CHANGESOFT
      myParticles[i].fSoft0 = dSoft;
#endif
      myParticles[i].soft = dSoft;
  }
  contribute(cb);
}

/**
 * \brief Adjust comoving softening to maintain constant physical
 * softening
 */
void TreePiece::physicalSoft(const double dSoftMax, const double dFac,
			     const int bSoftMaxMul, const CkCallback& cb) {
#ifdef CHANGESOFT
    CkAssert(dFac > 0.0);
    if (bSoftMaxMul) {		// dSoftMax is a maximum multiplier
	for(unsigned int i = 1; i <= myNumParticles; ++i) {
	    myParticles[i].soft = myParticles[i].fSoft0*dFac;
	    CkAssert(myParticles[i].soft > 0.0);
	    }
	}
    else {			// dSoftMax is an absolute limit
	CkAssert(dSoftMax > 0.0);
	for(unsigned int i = 1; i <= myNumParticles; ++i) {
	    myParticles[i].soft = myParticles[i].fSoft0*dFac;
	    if(myParticles[i].soft > dSoftMax) myParticles[i].soft = dSoftMax;
	    }
	}
#endif
    contribute(cb);
}

void TreePiece::growMass(int nGrowMass, double dDeltaM, const CkCallback& cb)
{
    for(unsigned int i = 1; i <= myNumParticles; ++i) {
	if(myParticles[i].iOrder < nGrowMass)
	    myParticles[i].mass += dDeltaM;
	}
    contribute(cb);
    }

/*
 * Gathers information for center of mass calculation
 * For each particle type the 0th and first mass moment is summed.
 */
void TreePiece::getCOM(const CkCallback& cb, int bLiveViz) {
    int i;
    double com[12]; // structure is m*position and mass for all
		    // particle types;
    for(i = 0; i < 12; i++)
	com[i] = 0;
    for(i = 1; i <= myNumParticles; ++i) {
	double m = myParticles[i].mass;
	if ( TYPETest(&(myParticles[i]), TYPE_GAS) ) {
	    com[0] += m*myParticles[i].position[0];
	    com[1] += m*myParticles[i].position[1];
	    com[2] += m*myParticles[i].position[2];
	    com[3] += m;
	    }
	else if ( TYPETest(&(myParticles[i]), TYPE_DARK) ) {
	    com[4] += m*myParticles[i].position[0];
	    com[5] += m*myParticles[i].position[1];
	    com[6] += m*myParticles[i].position[2];
	    com[7] += m;
	    }
	else if ( TYPETest(&(myParticles[i]), TYPE_STAR) ) {
	    com[8] += m*myParticles[i].position[0];
	    com[9] += m*myParticles[i].position[1];
	    com[10] += m*myParticles[i].position[2];
	    com[11] += m;
	    }
	}
    if(bLiveViz)		// Use LiveViz array
	lvProxy[thisIndex].ckLocal()->contribute(12*sizeof(double), com,
						 CkReduction::sum_double, cb);
    else
	contribute(12*sizeof(double), com, CkReduction::sum_double, cb);
}

/*
 * Gathers information for center of mass calculation for one type of
 * particle.
 */
void TreePiece::getCOMByType(int iType, const CkCallback& cb, int bLiveViz) {
    int i;
    double com[4]; // structure is m*position and mass

    for(i = 0; i < 4; i++)
	com[i] = 0;
    for(i = 1; i <= myNumParticles; ++i) {
	if ( TYPETest(&(myParticles[i]), iType) ) {
	    double m = myParticles[i].mass;
	    com[0] += m*myParticles[i].position[0];
	    com[1] += m*myParticles[i].position[1];
	    com[2] += m*myParticles[i].position[2];
	    com[3] += m;
	    }
	}
    if(bLiveViz)		// Use LiveViz array
	lvProxy[thisIndex].ckLocal()->contribute(4*sizeof(double), com,
						  CkReduction::sum_double, cb);
    else
	contribute(4*sizeof(double), com, CkReduction::sum_double, cb);
}

/// @brief structure for efficiently reading iOrders.
struct SortStruct {
  int iOrder;
  int iStore;
};

int CompSortStruct(const void * a, const void * b) {
  return ( ( ((struct SortStruct *) a)->iOrder < ((struct SortStruct *) b)->iOrder ? -1 : 1 ) );
}

/// @brief Read file of iOrders to set particle type.
///
/// This is used by the photogenic code of dumpframe.
void TreePiece::SetTypeFromFileSweep(int iSetMask, char *file,
	   struct SortStruct *ss, int nss, int *pniOrder, int *pnSet) {
  int niOrder = 0, nSet = 0;
  int iOrder, iOrderOld, nRet;
  FILE *fp;
  int iss;

  fp = fopen( file, "r" );
  assert( fp != NULL );

  iss = 0;
  iOrderOld = -1;
  while ( (nRet=fscanf( fp, "%d\n", &iOrder )) == 1 ) {
	niOrder++;
	assert( iOrder > iOrderOld );
	iOrderOld = iOrder;
	while (ss[iss].iOrder < iOrder) {
	  iss++;
	  if (iss >= nss) goto DoneSS;
	}
	if (iOrder == ss[iss].iOrder) {
	  TYPESet(&(myParticles[1 + ss[iss].iStore]),iSetMask);
	  nSet++;
	}
  }

DoneSS:
#if 0
  /* The following should only be done for debugging.  It doubles the
     amount of file reading. */
  /* Finish reading file to verify consistency across processors */
  while ( (nRet=fscanf( fp, "%d\n", &iOrder )) == 1 ) {
	niOrder++;
	assert( iOrder > iOrderOld );
	iOrderOld = iOrder;
	}
#endif
  fclose(fp);

  *pniOrder += niOrder;
  *pnSet += nSet;

  return;
}

/*
 * Set particle type by reading in iOrders from file
 */
void
TreePiece::setTypeFromFile(int iSetMask, char *file, const CkCallback& cb)
{
  struct SortStruct *ss;
  int i,nss;

  int niOrder = 0;
  int nSet = 0;

  nss = myNumParticles;
  ss = (struct SortStruct *) malloc(sizeof(*ss)*nss);
  assert( ss != NULL );

  for(i=0;i<nss;++i) {
	ss[i].iOrder = 	myParticles[i+1].iOrder;
	ss[i].iStore = i;
  }

  qsort( ss, nss, sizeof(*ss), CompSortStruct );

  SetTypeFromFileSweep(iSetMask, file, ss, nss, &niOrder, &nSet);

  free( ss );

  int nSetOut[2];
  nSetOut[0] = niOrder;
  nSetOut[1] = nSet;

  contribute(2*sizeof(int), nSetOut, CkReduction::sum_int, cb);
}

#include "DumpFrameData.h"

/*
 * Render this processors portion of the image
 */
void TreePiece::DumpFrame(InDumpFrame in, const CkCallback& cb, int liveVizDump) 
{
    void *Image = dfDataProxy.ckLocalBranch()->Image;
    //
    // XXX The dump frame expects arrays of fixed structures.  This is
    // incompatible with ChaNGa's design of having gas and star data
    // as "decorations".  For now we copy the needed data from the
    // decoration into an unused attribute of the base structure
    // (GravityParticle) so it will be available to dumpframe.
    // Fortunately, we have only one attribute to copy over
    // (fTimeForm), and we have a temporary, dtGrav, that can be
    // overloaded.
    //
    for (int i=1; i<=myNumParticles; ++i)
	if(myParticles[i].isStar())
	    myParticles[i].dtGrav = myParticles[i].fTimeForm();
    
    GravityParticle *p = &(myParticles[1]);

    dfRenderParticlesInit( &in, TYPE_GAS, TYPE_DARK, TYPE_STAR,
			   &p->position[0], &p->mass, &p->soft, &p->fBall,
			   &p->iType,
			   &p->dtGrav, // actually fTimeForm; see above
			   p, sizeof(*p) );
    dfRenderParticles( &in, Image, p, myNumParticles);
    
    if(!liveVizDump) 
	contribute(cb);
    else
      {
	// this is the RGB 3-byte/pixel image created from floating point image
	// data in dfFinishFrame - here we just create the pointer to pass in
	unsigned char *gray;
	

	dfFinishFrame(NULL, 0, 0, &in, Image, true, &gray);

	// final image assembly for liveViz before shipping the image
	// to the client.
	// This calls a reduction which may conflict with a TreePiece
	// reduction.   Hence we use a shadow array "lvProxy" to avoid
	// this conflict.
	liveVizDeposit(savedLiveVizMsg, 0, 0, in.nxPix, in.nyPix, gray,
		       lvProxy[thisIndex].ckLocal(), sum_image_data);

	savedLiveVizMsg = NULL;
	free(gray);
      }
    }

/**
 * Overall start of building Tree.  For Oct trees, this begins by a
 * reduction of each treepieces boundaries to DataManager::collectSplitters.
 *
 * For ORB trees, this continues on to TreePiece::startORBTreeBuild.
 */

#ifdef PUSH_GRAVITY
void TreePiece::buildTree(int bucketSize, const CkCallback& cb, bool _merge) 
#else
void TreePiece::buildTree(int bucketSize, const CkCallback& cb)
#endif
{

#if COSMO_DEBUG > 1
  char fout[100];
  sprintf(fout,"tree.%d.%d.after",thisIndex,iterationNo);
  ofstream ofs(fout);
  for (int i=1; i<=myNumParticles; ++i)
    ofs << keyBits(myParticles[i].key,KeyBits) << " " << myParticles[i].position[0] << " "
        << myParticles[i].position[1] << " " << myParticles[i].position[2] << endl;
  ofs.close();
#endif

  maxBucketSize = bucketSize;
  callback = cb;
  myTreeParticles = myNumParticles;

  deleteTree();
  if(bucketReqs != NULL) {
    delete[] bucketReqs;
    bucketReqs = NULL;
  }
#ifdef PUSH_GRAVITY
  // used to indicate whether trees on SMP node should be
  // merged or not: we do not merge trees when pushing, to
  // increase concurrency
  doMerge = _merge; 
#endif

  // decide which logic are we using to divide the particles: Oct or ORB
  switch (useTree) {
  case Binary_Oct:
  case Oct_Oct:
    Key bounds[2];
    if (myNumParticles > 0) {
#ifdef COSMO_PRINT
      CkPrintf("[%d] Keys: %016llx %016llx\n",thisIndex,myParticles[1].key,myParticles[myNumParticles].key);
#endif
      bounds[0] = myParticles[1].key;
      bounds[1] = myParticles[myNumParticles].key;

      //CkPrintf("[%d] bounds %llx - %llx\n", thisIndex, bounds[0], bounds[1]);

      contribute(2 * sizeof(Key), bounds, CkReduction::concat, CkCallback(CkIndex_DataManager::collectSplitters(0), CProxy_DataManager(dataManagerID)));
    } else {
      // No particles assigned to this TreePiece
      contribute(0, NULL, CkReduction::concat, CkCallback(CkIndex_DataManager::collectSplitters(0), CProxy_DataManager(dataManagerID)));
    }
    break;
  case Binary_ORB:
    // WARNING: ORB trees do not allow TreePieces to have 0 particles!
    contribute(CkCallback(CkIndex_TreePiece::startORBTreeBuild(0), thisArrayID));
    break;
  }
}

/*
class KeyDouble {
  Key first;
  Key second;
public:
  inline bool operator<(const KeyDouble& k) const {
    return first < k.first;
  }
};
*/

void TreePiece::quiescence() {
  /*
  char fout[100];
  sprintf(fout,"tree.%d.%d",thisIndex,iterationNo);
  ofstream ofs(fout);
  printTree(root,ofs);
  ofs.close();
  */
  /*
  CkPrintf("[%d] quiescence, %d left\n",thisIndex,momentRequests.size());
  for (MomentRequestType::iterator iter = momentRequests.begin(); iter != momentRequests.end(); iter++) {
    CkVec<int> *l = iter->second;
    for (int i=0; i<l->length(); ++i) {
      CkPrintf("[%d] quiescence: %s to %d\n",thisIndex,keyBits(iter->first,63).c_str(),(*l)[i]);
    }
  }
  */

  CkPrintf("[%d] quiescence detected, pending %d total %d\n",
                          thisIndex, sLocalGravityState->myNumParticlesPending,
                          numBuckets);

  for (unsigned int i=0; i<numBuckets; ++i) {
    int remaining;
    remaining = sRemoteGravityState->counterArrays[0][i]
                + sLocalGravityState->counterArrays[0][i];
    if (remaining != 0)
      CkPrintf("[%d] requests for %d remaining l %d r %d\n",
                thisIndex,i,
                sRemoteGravityState->counterArrays[0][i],
                sLocalGravityState->counterArrays[0][i]);
  }

  CkPrintf("quiescence detected!\n");
  mainChare.niceExit();
}

GenericTreeNode *TreePiece::get3DIndex() {
  GenericTreeNode *node = root;
  while (node != NULL && node->getType() != Internal && node->getType() != Bucket) {
    int next = -1;
    GenericTreeNode *child;
    for (int i=0; i<node->numChildren(); ++i) {
      child = node->getChildren(i);
      if (child->getType() == Internal || child->getType() == Boundary) {
        if (next != -1) return NULL;
        next = i;
      }
    }
    if (next == -1) return NULL;
    node = node->getChildren(next);
  }
  return node;
}

/*****************ORB**********************/
/*void TreePiece::receiveBoundingBoxes(BoundingBoxes *msg){

  //OrientedBox<float> *boxesMsg = static_cast<OrientedBox<float>* >(msg->getData())
  //boxes = new OrientedBox<float>[numTreePieces];

  OrientedBox<float> *boxesMsg = msg->boxes;
  copy(boxesMsg,boxesMsg+numTreePieces,boxes);

  contribute(CkCallback(CkIndex_TreePiece::startORBTreeBuild(0), thisArrayID));
}*/

void TreePiece::startORBTreeBuild(CkReductionMsg* m){
  delete m;

  if (dm == NULL) {
      dm = (DataManager*)CkLocalNodeBranch(dataManagerID);
      }

  if (myNumParticles == 0) {
    // No particle assigned to this TreePiece
#ifdef PUSH_GRAVITY
    if(doMerge){
#endif
      if (verbosity > 3) ckerr << "TreePiece " << thisIndex << ": No particles, finished tree build" << endl;
      contribute(sizeof(callback), &callback, CkReduction::random, CkCallback(CkIndex_DataManager::combineLocalTrees((CkReductionMsg*)NULL), CProxy_DataManager(dataManagerID)));
#ifdef PUSH_GRAVITY
    }
    else{
      // if not combining trees, return control to mainchare
      contribute(callback);
    }
#endif
    return;
  }
  myParticles[0].key = thisIndex;
  myParticles[myNumParticles+1].key = thisIndex;

  compFuncPtr[0]= &comp_dim0;
  compFuncPtr[1]= &comp_dim1;
  compFuncPtr[2]= &comp_dim2;

  pTreeNodes = new NodePool();
  root = pTreeNodes->alloc_one(1, numTreePieces>1?Tree::Boundary:Tree::Internal,
			       0, myNumParticles+1, 0);

  if (thisIndex == 0) root->firstParticle ++;
  if (thisIndex == (int)numTreePieces-1) root->lastParticle --;
  root->particleCount = myNumParticles;
  nodeLookupTable[(Tree::NodeKey)1] = root;

  //root->key = firstPossibleKey;
  root->boundingBox = boundingBox;
  //nodeLookup[root->lookupKey()] = root;
  numBuckets = 0;
  bucketList.clear();

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
          CkEntryOptions opts;
          opts.setPriority((unsigned int) -100000000);
	  streamingProxy[(*l)[i]].receiveRemoteMoments(nodeKey, node->getType(), node->firstParticle, node->particleCount, node->moments, node->boundingBox, node->bndBoxBall, node->iParticleTypes, &opts);
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

  if (numTreePieces == 1) {
    treeBuildComplete();
  }

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

  //CkPrintf("[%d] in build ORB Tree, level:%d\n",thisIndex,level);
  if (level == NodeKeyBits-1) {
    ckerr << thisIndex << ": TreePiece(ORB): This piece of tree has exhausted all the bits in the keys.  Super double-plus ungood!" << endl;
    ckerr << "Left particle: " << (node->firstParticle) << " Right particle: " << (node->lastParticle) << endl;
    ckerr << "Left key : " << keyBits((myParticles[node->firstParticle]).key, KeyBits).c_str() << endl;
    ckerr << "Right key: " << keyBits((myParticles[node->lastParticle]).key, KeyBits).c_str() << endl;
    return;
  }

  CkAssert(node->getType() == Boundary || node->getType() == Internal);

  node->makeOrbChildren(myParticles, myNumParticles, level, chunkRootLevel,
			compFuncPtr, (domainDecomposition == ORB_space_dec),
			pTreeNodes);
  node->rungs = 0;

  GenericTreeNode *child;
#if INTERLIST_VER > 0
  int bucketsBeneath = 0;
#endif
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
              child->makeEmpty();
	      child->remoteIndex = thisIndex;
      } else {
	      child->remoteIndex = getResponsibleIndex(first, last);
	      // if we have a remote child, the node is a Boundary. Thus count that we
	      // have to receive one more message for the NonLocal node
	      node->remoteIndex --;
	      // request the remote chare to fill this node with the Moments
              CkEntryOptions opts;
              opts.setPriority((unsigned int) -110000000);
	      streamingProxy[child->remoteIndex].requestRemoteMoments(child->getKey(), thisIndex, &opts);
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
      //child->bucketListIndex=numBuckets;
      child->startBucket=numBuckets;
#endif
      numBuckets++;
      if (node->getType() != Boundary) {
	  node->moments += child->moments;
	  node->bndBoxBall.grow(child->bndBoxBall);
	  node->iParticleTypes |= child->iParticleTypes;
	  }
      if (child->rungs > node->rungs) node->rungs = child->rungs;
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
      if (node->getType() != Boundary) {
	  node->moments += child->moments;
	  node->bndBoxBall.grow(child->bndBoxBall);
	  node->iParticleTypes |= child->iParticleTypes;
	  }
      if (child->rungs > node->rungs) node->rungs = child->rungs;
    }
#if INTERLIST_VER > 0
    bucketsBeneath += child->numBucketsBeneath;
#endif
  }

#if INTERLIST_VER > 0
  node->numBucketsBeneath = bucketsBeneath;
#endif

  /* The old version collected Boundary nodes, the new version collects NonLocal nodes */

  if (node->getType() == Internal) {
    calculateRadiusFarthestCorner(node->moments, node->boundingBox);
  }

}
/******************************************/

/**
 * Actual treebuild for each treepiece.  Each treepiece begins its
 * treebuild at the global root.  During the treebuild, the
 * nodeLookupTable, a mapping from keys to nodes, is constructed.
 * Also, the bucket list is constructed.  Constructing moments for the
 * boundary nodes requires requesting the moments of nodes on other
 * pieces.  After all the moments are constructed, the treepiece
 * passes its root to the DataManager via DataManager::notifyPresence.
 *
 * After the local treebuild is finished
 * DataManager::combineLocalTrees is called, and the treebuild phase
 * is finished.
 */

void TreePiece::startOctTreeBuild(CkReductionMsg* m) {
  delete m;

  if (dm == NULL) {
      dm = (DataManager*)CkLocalNodeBranch(dataManagerID);
  }
  
  if (myNumParticles == 0) {
#ifdef PUSH_GRAVITY
    if(doMerge){
#endif
      // No particle assigned to this TreePiece
      if (verbosity > 3) ckerr << "TreePiece " << thisIndex << ": No particles, finished tree build" << endl;
      contribute(sizeof(callback), &callback, CkReduction::random, CkCallback(CkIndex_DataManager::combineLocalTrees((CkReductionMsg*)NULL), CProxy_DataManager(dataManagerID)));
#ifdef PUSH_GRAVITY
    }
    else{
      // if not merging, return control to mainchare
      contribute(callback);
    }
#endif
    return;
  }
  //CmiLock(dm->__nodelock);

  // The boundary particles are set to the splitters of the
  // neighboring pieces, i.e. the nearest particles that are not in
  // this treepiece.
  myPlace = find(dm->responsibleIndex.begin(), dm->responsibleIndex.end(), thisIndex) - dm->responsibleIndex.begin();
  if(myPlace == 0)
    myParticles[0].key = firstPossibleKey;
  else
    myParticles[0].key = dm->splitters[2 * myPlace - 1];

  if(myPlace == dm->responsibleIndex.size() - 1)
    myParticles[myNumParticles + 1].key = lastPossibleKey;
  else
    myParticles[myNumParticles + 1].key = dm->splitters[2 * myPlace + 2];

  CkAssert(myParticles[1].key >= myParticles[0].key);
  CkAssert(myParticles[myNumParticles + 1].key >= myParticles[myNumParticles].key);

  // create the root of the global tree
  switch (useTree) {
  case Binary_Oct:
    pTreeNodes = new NodePool;
    root = pTreeNodes->alloc_one(1, numTreePieces>1?Tree::Boundary:Tree::Internal, 0, myNumParticles+1, 0);
    break;
  case Oct_Oct:
    //root = new OctTreeNode(1, Tree::Boundary, 0, myNumParticles+1, 0);
    break;
  default:
    CkAbort("We should have never reached here!");
  }

  if (myPlace == 0) root->firstParticle ++;
  if (myPlace == dm->responsibleIndex.size()-1) root->lastParticle --;
  root->particleCount = myNumParticles;
  nodeLookupTable[(Tree::NodeKey)1] = root;

  root->boundingBox = boundingBox;
  numBuckets = 0;
  bucketList.clear();

  if(verbosity > 3)
    ckerr << "TreePiece " << thisIndex << ": Starting tree build" << endl;

#if INTERLIST_VER > 0
  root->startBucket=0;
#endif
  // recursively build the tree

  double start;
  try {
#if defined MERGE_REMOTE_REQUESTS
        LocalTreeTraversal traversal;
        // 1. Construct top level of your tree first, being careful
        // not to send out requests for remote moments, and not
        // descending into treepiece-local portions of the tree
        // This is called the Remote tree
        RemoteTreeBuilder remoteTreeBuilder(this,false);
        traversal.dft(root,&remoteTreeBuilder,0);

        // 2. Next, merge these top-level trees, causing requests
        // to be sent out to the owners of the remote moments required
        // on this PE.
        peTreeMergerProxy.ckLocalBranch()->mergeNonLocalRequests(root,this);

        // 3. Then, construct the treepiece-local portions of the tree
        // and respond to remote requests for local moments
        // 4. Finally, do a node-level merge on PE trees

        // Items 3 and 4 are performed after PE-level merge of trees has been done
        return;

#elif defined SPLIT_PHASE_TREE_BUILD
        LocalTreeTraversal traversal; 
        // Try to construct the non-local part of the tree: 
        // this results in early dispatch of remote moment
        // requests
        RemoteTreeBuilder w1(this,true);
        traversal.dft(root,&w1,0);

        // Then construct the local parts of the tree
        LocalTreeBuilder w2(this);
        traversal.dft(root,&w2,0);

#else
        start = CmiWallTimer();
	buildOctTree(root, 0);
        traceUserBracketEvent(tbRecursiveUE,start,CmiWallTimer());
#endif

	}
  catch (std::bad_alloc) {
	CkAbort("Out of memory in treebuild");
	}

  processRemoteRequestsForMoments();

  if(verbosity > 3){
    ckerr << "TreePiece " << thisIndex << ": Number of buckets: " << numBuckets << endl;
    ckerr << "TreePiece " << thisIndex << ": Finished tree build, resolving boundary nodes" << endl;
  }

  //CmiUnlock(dm->__nodelock);
  if (numTreePieces == 1) {
    treeBuildComplete();
  }
}

void TreePiece::recvActiveRung(int arung, const CkCallback& cb) {
  activeRung = arung;
  contribute(cb);
}

void TreePiece::sendRequestForNonLocalMoments(GenericTreeNode *pickedNode){
  int first, last;
  bool isShared = nodeOwnership(pickedNode->getKey(), first, last);
  CkAssert(!isShared);
  if (last >= first) {
    // Choose a piece from among the owners from which to
    // request moments in such a way that if I am a piece with a
    // higher index, I request from a higher indexed treepiece.
    pickedNode->remoteIndex = getResponsibleIndex(first,last);
    // request the remote chare to fill this node with the Moments
    CkEntryOptions opts;
    opts.setPriority((unsigned int) -110000000);
    //CkPrintf("[%d] phase I request moments for %llu from %d\n", tp->thisIndex, node->getKey(), node->remoteIndex);
    treeProxy[pickedNode->remoteIndex].requestRemoteMoments(pickedNode->getKey(), thisIndex, &opts);
  }
}

void TreePiece::processRemoteRequestsForMoments(){
  double start = CmiWallTimer();
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
          CkEntryOptions opts;
          opts.setPriority((unsigned int) -100000000);
	  streamingProxy[(*l)[i]].receiveRemoteMoments(nodeKey, node->getType(), node->firstParticle, node->particleCount, node->moments, node->boundingBox, node->bndBoxBall, node->iParticleTypes, &opts);
	  //CkPrintf("[%d] sending moments of %s to %d upon treebuild finished\n",thisIndex,keyBits(node->getKey(),63).c_str(),(*l)[i]);
      }
      delete l;
      momentRequests.erase(node->getKey());
    }
  }
  traceUserBracketEvent(tbFlushRequestsUE,start,CmiWallTimer());

}

void NonEmptyTreePieceCounter::addLocation(CkLocation &loc){
  const int *indexData = loc.getIndex().data();
  TreePiece *tp = treeProxy[indexData[0]].ckLocal();
  int np = tp->getNumParticles();
  if(np > 0) count++;
}

void NonEmptyTreePieceCounter::reset() {
  count = 0;
}

void TreePiece::mergeNonLocalRequestsDone(){
  // 3. Construct the treepiece-local portions of the tree

  MERGE_REMOTE_REQUESTS_VERBOSE("[%d] mergeNonLocalRequestsDone\n", thisIndex);

  LocalTreeTraversal traversal;
  LocalTreeBuilder localTreeBuilder(this);
  traversal.dft(root,&localTreeBuilder,0);
  localTreeBuildComplete = true;

  // at this point, I might have completed building
  // the entire tree, since:
  // (a) I must have finished building the 
  //     RemoteTree, otherwise the PETreeMerger 
  //     wouldn't have invoked this method on me
  
  // (b) I might have received all the moments for
  //     which I was a client
  
  // 4. Respond to remote requests for local moments
  processRemoteRequestsForMoments();

  // If none of the root's children will be updated
  // with remote moments, we are done
  if(root->remoteIndex == thisIndex){
    treeBuildComplete();
  }
}

/*
void TreePiece::flushNonLocalMomentsClients(){
  // if we are merging remote requests, there might
  // be clients on my PE for some of my internal 
  // nodes
  MERGE_REMOTE_REQUESTS_VERBOSE("[%d] flushNonLocalMomentsClients\n",thisIndex);
  std::map<NodeKey,NonLocalMomentsClientList>::iterator it;
  for(it = nonLocalMomentsClients.begin();
      it != nonLocalMomentsClients.end();
      ++it){
    deliverMomentsToClients(it);
  }
  nonLocalMomentsClients.clear();
}
*/

/// Determine who are all the owners of this node
/// @return true if the caller is part of the owners, false otherwise
bool TreePiece::nodeOwnership(const Tree::NodeKey nkey, int &firstOwner, int &lastOwner) {

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
      // From the nodekey determine the particle keys that bound this node.
    Key firstKey = Key(nkey);
    Key lastKey = Key(nkey + 1);
    const Key mask = Key(1) << KeyBits;
    while (! (firstKey & mask)) {
      firstKey <<= 1;
      lastKey <<= 1;
    }
    firstKey &= ~mask;
    lastKey &= ~mask;
    lastKey -= 1;
    Key *locLeft = lower_bound(dm->splitters, dm->splitters + dm->numSplitters, firstKey);
    Key *locRight = upper_bound(locLeft, dm->splitters + dm->numSplitters, lastKey);
    firstOwner = (locLeft - dm->splitters) >> 1;
    lastOwner = (locRight - dm->splitters - 1) >> 1;
#if COSMO_PRINT > 1
    std::string str = keyBits(nkey,KeyBits);
    CkPrintf("[%d] NO: key=%s, first=%d, last=%d\n",thisIndex,str.c_str(),locLeft-dm->splitters,locRight-dm->splitters);
#endif
  }
  return (myPlace >= firstOwner && myPlace <= lastOwner);
}

/** A recursive algorithm for building my tree.
    Examines successive bits in the particles' keys, looking for splits.
    Each bit is a level of nodes in the tree.  We keep going down until
    we can bucket the particles.
*/
void TreePiece::buildOctTree(GenericTreeNode * node, int level) {

#ifdef TREE_BREADTH_FIRST
  CkQ<GenericTreeNode*> *queue = new CkQ<GenericTreeNode*>(1024);
  CkQ<GenericTreeNode*> *queueNext = new CkQ<GenericTreeNode*>(1024);
  queue->enq(node);

  GenericTreeNode *rootNode = node;
  while (1) {
    node = queue->deq();
    if (node == NULL) {
      node = queueNext->deq();
      CkQ<GenericTreeNode*> *tmp = queue;
      queue = queueNext;
      queueNext = tmp;
      level++;
    }
    if (node == NULL) break;
#endif

  if (level == NodeKeyBits-2) {
    ckerr << thisIndex << ": TreePiece: This piece of tree has exhausted all the bits in the keys.  Super double-plus ungood!" << endl;
    ckerr << "Left particle: " << (node->firstParticle) << " Right particle: " << (node->lastParticle) << endl;
    ckerr << "Left key : " << keyBits((myParticles[node->firstParticle]).key, KeyBits).c_str() << endl;
    ckerr << "Right key: " << keyBits((myParticles[node->lastParticle]).key, KeyBits).c_str() << endl;
    ckerr << "Node type: " << node->getType() << endl;
    ckerr << "myNumParticles: " << myNumParticles << endl;
    CkAbort("Tree is too deep!");
    return;
  }

  CkAssert(node->getType() == Boundary || node->getType() == Internal);

  node->makeOctChildren(myParticles, myNumParticles, level, pTreeNodes);
  // The boundingBox was used above to determine the spacially equal
  // split between the children.  Now reset it so it can be calculated
  // from the particle positions.
  node->boundingBox.reset();
  node->rungs = 0;

  GenericTreeNode *child;
#if INTERLIST_VER > 0
  int bucketsBeneath = 0;
#endif
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
	child->makeEmpty();
	child->remoteIndex = thisIndex;
      } else {
	child->remoteIndex = getResponsibleIndex(first, last);
	// if we have a remote child, the node is a Boundary. Thus count that we
	// have to receive one more message for the NonLocal node
	node->remoteIndex --;
	// request the remote chare to fill this node with the Moments
        CkEntryOptions opts;
        opts.setPriority((unsigned int) -110000000);
	streamingProxy[child->remoteIndex].requestRemoteMoments(child->getKey(), thisIndex, &opts);
	//CkPrintf("[%d] asking for moments of %s to %d\n",thisIndex,keyBits(child->getKey(),63).c_str(),child->remoteIndex);
      }
    } else if (child->getType() == Internal
              && (child->lastParticle - child->firstParticle < maxBucketSize
                  || level >= NodeKeyBits-3)) {
       if(level >= NodeKeyBits-3
	  && child->lastParticle - child->firstParticle >= maxBucketSize)
           ckerr << "Truncated tree with "
                 << child->lastParticle - child->firstParticle
                 << " particle bucket" << endl;
       
      CkAssert(child->firstParticle != 0 && child->lastParticle != myNumParticles+1);
      child->remoteIndex = thisIndex;
      child->makeBucket(myParticles);
      bucketList.push_back(child);
#if INTERLIST_VER > 0
      child->startBucket=numBuckets;
#endif
      numBuckets++;
      if (node->getType() != Boundary) {
        node->moments += child->moments;
        node->boundingBox.grow(child->boundingBox);
        node->bndBoxBall.grow(child->bndBoxBall);
	node->iParticleTypes |= child->iParticleTypes;
      }
      if (child->rungs > node->rungs) node->rungs = child->rungs;
    } else if (child->getType() == Empty) {
      child->remoteIndex = thisIndex;
    } else {
      if (child->getType() == Internal) child->remoteIndex = thisIndex;
      // else the index is already 0
#ifdef TREE_BREADTH_FIRST
      queueNext->enq(child);
#else
      buildOctTree(child, level+1);
#endif
      // if we have a Boundary child, we will have to compute it's multipole
      // before we can compute the multipole of the current node (and we'll do
      // it in receiveRemoteMoments)
      if (child->getType() == Boundary) node->remoteIndex --;
#ifndef TREE_BREADTH_FIRST
      if (node->getType() != Boundary) {
        node->moments += child->moments;
        node->boundingBox.grow(child->boundingBox);
        node->bndBoxBall.grow(child->bndBoxBall);
	node->iParticleTypes |= child->iParticleTypes;
      }
      // for the rung information we can always do now since it is a local property
      if (child->rungs > node->rungs) node->rungs = child->rungs;
#endif
    }

#if INTERLIST_VER > 0
    bucketsBeneath += child->numBucketsBeneath;
#endif
  }

#if INTERLIST_VER > 0
  node->numBucketsBeneath = bucketsBeneath;
#endif

  /* The old version collected Boundary nodes, the new version collects NonLocal nodes */

#ifndef TREE_BREADTH_FIRST
  if (node->getType() == Internal) {
    calculateRadiusFarthestCorner(node->moments, node->boundingBox);
  }
#endif

#ifdef TREE_BREADTH_FIRST
  }
  growBottomUp(rootNode);
#endif
}

#ifdef TREE_BREADTH_FIRST
void TreePiece::growBottomUp(GenericTreeNode *node) {
  GenericTreeNode *child;
  for (unsigned int i=0; i<node->numChildren(); ++i) {
    child = node->getChildren(i);
    if (child->getType() == NonLocal ||
        (child->getType() == Bucket) ||
        child->getType() == Empty) continue;
    growBottomUp(child);
    if (node->getType() != Boundary) {
      node->moments += child->moments;
      node->boundingBox.grow(child->boundingBox);
      node->bndBoxBall.grow(child->bndBoxBall);
      node->iParticleTypes |= child->iParticleTypes;
    }
    if (child->rungs > node->rungs) node->rungs = child->rungs;
  }
  if (node->getType() == Internal) {
    calculateRadiusFarthestCorner(node->moments, node->boundingBox);
  }
}
#endif

/// \brief entry method to obtain the moments of a node
void TreePiece::requestRemoteMoments(const Tree::NodeKey key, int sender) {
  GenericTreeNode *node = keyToNode(key);
  if (node != NULL && (node->getType() == Empty || node->moments.totalMass > 0)) {
      CkEntryOptions opts;
      opts.setPriority((unsigned int) -100000000);
      streamingProxy[sender].receiveRemoteMoments(key, node->getType(), node->firstParticle, node->particleCount, node->moments, node->boundingBox, node->bndBoxBall, node->iParticleTypes, &opts);
    //CkPrintf("[%d] sending moments of %s to %d directly\n",thisIndex,keyBits(node->getKey(),63).c_str(),sender);
  } else {
      // Save request for when we've calculated the moment.
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

/// \brief response from requestRemoteMoments
void TreePiece::receiveRemoteMoments(const Tree::NodeKey key,
				     Tree::NodeType type,
				     int firstParticle,
				     int numParticles,
				     const MultipoleMoments& moments,
				     const OrientedBox<double>& box,
				     const OrientedBox<double>& boxBall,
				     const unsigned int iParticleTypes) {
  GenericTreeNode *node = keyToNode(key);
  CkAssert(node != NULL);
  MERGE_REMOTE_REQUESTS_VERBOSE("[%d] receiveRemoteMoments %llu\n",thisIndex,key);
  // assign the incoming moments to the node
  if (type == Empty) node->makeEmpty();
  else {
    if (type == Bucket || type == NonLocalBucket) {
      node->setType(NonLocalBucket);
      node->firstParticle = firstParticle;
      node->lastParticle = firstParticle + numParticles - 1;
    }
    node->particleCount = numParticles;
    node->moments = moments;
    node->boundingBox = box;
    node->bndBoxBall = boxBall;
    node->iParticleTypes = iParticleTypes;
  }

#ifdef MERGE_REMOTE_REQUESTS
  deliverMomentsToClients(node);
#endif
 
  // look if we can compute the moments of some ancestors, and eventually send
  // them to a requester
  GenericTreeNode *parent = node->parent;
  if(parent != NULL){
    CkAssert(parent->getType() == Boundary);
    parent->remoteIndex++;
  }

#ifdef MERGE_REMOTE_REQUESTS
 // If we are merging remote requests, since we only build the skeleton
  // top-level tree (i.e. RemoteTree) in the first pass (RemoteTreeBuilder)
  // we might not have computed the moments of internal siblings of this node
  // Moreover, we do not decrement the remoteIndex of the parent for these 
  // Internal children whose moments have not yet been computed. Therefore, 
  // we want to avoid the situation wherein a Boundary parent node's non-Internal
  // child's moments are received in this invocation, but the moments of the
  // parent have not yet been computed. We can be sure that the Internal childrens'
  // moments have been computed when the LocalTreeBuilder has finished
  if(!localTreeBuildComplete) return;
#endif

  while (parent != NULL && 
         parent->remoteIndex == 0) {
    GenericTreeNode *parentsParent = boundaryParentReady(parent);
#ifdef MERGE_REMOTE_REQUESTS
    deliverMomentsToClients(parent);
#endif
    parent = parentsParent;
  }

  if (parent == NULL) {
    // if we are here then we are at the root, and thus we have finished to get
    // all moments
    //CkPrintf("[%d] contributing after building the tree\n",thisIndex);
    treeBuildComplete();
  }// else CkPrintf("[%d] still missing one child of %s\n",thisIndex,keyBits(parent->getKey(),63).c_str());
}

GenericTreeNode *TreePiece::boundaryParentReady(GenericTreeNode *parent){
  // compute the multipole for the parent
  MERGE_REMOTE_REQUESTS_VERBOSE("[%d] boundaryParentReady %llu\n",thisIndex,parent->getKey());
  parent->particleCount = 0;
  parent->remoteIndex = thisIndex; // reset the reference index to ourself
  GenericTreeNode *child;
  for (unsigned int i=0; i<parent->numChildren(); ++i) {
    child = parent->getChildren(i);
    parent->particleCount += child->particleCount;
    accumulateMomentsFromChild(parent,child);
  }
  calculateRadiusFarthestCorner(parent->moments, parent->boundingBox);
  // check if someone has requested this node
  MomentRequestType::iterator iter;
  if ((iter = momentRequests.find(parent->getKey())) != momentRequests.end()) {
    CkVec<int> *l = iter->second;
    for (int i=0; i<l->length(); ++i) {
      CkEntryOptions opts;
      opts.setPriority((unsigned int) -100000000);
      streamingProxy[(*l)[i]].receiveRemoteMoments(parent->getKey(), parent->getType(), parent->firstParticle, parent->particleCount, parent->moments, parent->boundingBox, parent->bndBoxBall, parent->iParticleTypes, &opts);
      //CkPrintf("[%d] sending moments of %s to %d\n",thisIndex,keyBits(parent->getKey(),63).c_str(),(*l)[i]);
    }
    delete l;
    momentRequests.erase(parent->getKey());
  }

  // go to the next ancestor
  GenericTreeNode *node = parent;
  parent = node->parent;
  // we just computed the current parent's child's
  // moments. The current parent must be a Boundary
  if(parent != NULL){
    CkAssert(parent->getType() == Boundary);
    parent->remoteIndex++;
  }

  return parent;
}

void TreePiece::accumulateMomentsFromChild(GenericTreeNode *parent, GenericTreeNode *child){
  parent->moments += child->moments;
  parent->boundingBox.grow(child->boundingBox);
  parent->bndBoxBall.grow(child->bndBoxBall);
  parent->iParticleTypes |= child->iParticleTypes;
}

/// @brief determine if moments of node have been requested.
///
/// This is used by MERGE_REMOTE_REQUESTS to distribute the received
/// moments among Treepieces on a core.
void TreePiece::deliverMomentsToClients(GenericTreeNode *node){
  std::map<NodeKey,NonLocalMomentsClientList>::iterator it;
  it = nonLocalMomentsClients.find(node->getKey());
  if(it == nonLocalMomentsClients.end()) return;

  CkAssert(it->second.targetNode == node);
  deliverMomentsToClients(it);
}

/// @brief send moments to all the pieces that have requested them.
///
/// This is used by MERGE_REMOTE_REQUESTS to distribute the received
/// moments among Treepieces on a core.
void TreePiece::deliverMomentsToClients(const std::map<NodeKey,NonLocalMomentsClientList>::iterator &it){
  NonLocalMomentsClientList &entry = it->second;
  CkVec<NonLocalMomentsClient> &clients = entry.clients;
  GenericTreeNode *node = entry.targetNode;

  /*
  if(node->remoteIndex < 0){
    CkPrintf("[%d] sending moments that are not ready: node key %llu\n", thisIndex, node->getKey());
    CkAbort("uninitialized moments being communicated\n");
  }
  */
  CkAssert(node->remoteIndex >= 0);

  for(int i = 0; i < clients.length(); i++){
    MERGE_REMOTE_REQUESTS_VERBOSE("[%d] send %llu (%s) moments to %d\n", thisIndex, node->getKey(), typeString(node->getType()),clients[i].clientTreePiece->getIndex());
    clients[i].clientTreePiece->receiveRemoteMoments(node->getKey(),node->getType(),node->firstParticle,node->particleCount,node->moments,node->boundingBox,node->bndBoxBall,node->iParticleTypes);
  }
  clients.clear();
  nonLocalMomentsClients.erase(it);

}

void TreePiece::treeBuildComplete(){

#ifdef MERGE_REMOTE_REQUESTS
  MERGE_REMOTE_REQUESTS_VERBOSE("[%d] treeBuildComplete\n", thisIndex);
  // reset
  CkAssert(localTreeBuildComplete);
  localTreeBuildComplete = false;
  /*
  if(!nonLocalMomentsClients.empty()){
    std::map<NodeKey,NonLocalMomentsClientList>::iterator it;
    for(it = nonLocalMomentsClients.begin(); it != nonLocalMomentsClients.end(); ++it){
      CkPrintf("[%d] BAD TREE BUILD: outstanding clients for node %llu\n", thisIndex, it->second.targetNode->getKey());
    }
    CkAbort("bad tree build\n");
  }
  */
  CkAssert(nonLocalMomentsClients.empty());
#endif

#ifdef PUSH_GRAVITY
  if(doMerge){
#endif

    dm->notifyPresence(root, this);
    contribute(sizeof(callback), &callback, CkReduction::random, CkCallback(CkIndex_DataManager::combineLocalTrees((CkReductionMsg*)NULL), CProxy_DataManager(dataManagerID)));

#ifdef PUSH_GRAVITY
  }
  else{
    // if not merging, return control to mainchare
    contribute(callback);
  }
#endif
}

bool bIsReplica(int reqID)
{
    int offsetcode = reqID >> 22;
    int x = (offsetcode & 0x7) - 3;
    int y = ((offsetcode >> 3) & 0x7) - 3;
    int z = ((offsetcode >> 6) & 0x7) - 3;

    return x || y || z;
    }

int decodeReqID(int reqID)
{
    const int offsetmask = 0x1ff << 22;

    return reqID & (~offsetmask);
    }

int encodeOffset(int reqID, int x, int y, int z)
{
    // Bit limitations follow
    CkAssert(x > -4 && x < 4);
    CkAssert(y > -4 && y < 4);
    CkAssert(z > -4 && z < 4);
    // Replica in each direction is mapped to 0-7 range
    int offsetcode = (x + 3) | ((y+3) << 3) | ((z+3) << 6);

    // Assume we only have 32 bits to work with (minus sign bit)
    // 4 million buckets is our limitation
    CkAssert(reqID < (1 << 22));
    return reqID | (offsetcode << 22);
    }

// Add reqID to encoded offset
int reEncodeOffset(int reqID, int offsetID)
{
    const int offsetmask = 0x1ff << 22;

    return reqID | (offsetmask & offsetID);
}


/**
 * Initialize all particles for gravity force calculation.
 * This includes zeroing out the acceleration and potential.
 */
void TreePiece::initBuckets() {
  int ewaldCondition = (bEwald ? 0 : 1);
  for (unsigned int j=0; j<numBuckets; ++j) {
    GenericTreeNode* node = bucketList[j];

    // TODO: active bounds may give a performance boost in the
    // multi-timstep regime.
    // node->boundingBox.reset();  // XXXX dangerous should have separate
				// Active bounds
    for(int i = node->firstParticle; i <= node->lastParticle; ++i) {
      if (myParticles[i].rung >= activeRung) {
        myParticles[i].treeAcceleration = 0;
        myParticles[i].potential = 0;
	myParticles[i].dtGrav = 0;
	// node->boundingBox.grow(myParticles[i].position);
      }
    }
    bucketReqs[j].finished = ewaldCondition;

/*#if COSMO_DEBUG > 1
    if(iterationNo==1 || listMigrated==true){
      std::set<Tree::NodeKey> *list = new std::set<Tree::NodeKey>();
      bucketcheckList.push_back(*list);
    }
#endif*/
  }
#if COSMO_DEBUG > 1 || defined CHANGA_REFACTOR_WALKCHECK || defined CHANGA_REFACTOR_WALKCHECK_INTERLIST
  bucketcheckList.resize(numBuckets);
#endif
}

void TreePiece::startNextBucket() {
  int currentBucket = sLocalGravityState->currentBucket;
  if(currentBucket >= numBuckets)
    return;

#if INTERLIST_VER > 0
  // no need to do the following, because interlist walk and interlist compute are only
  // ever associated with each other, unlike the topdown walk, which may be associated
  // with gravity or prefetch objects
  // sInterListWalk->init(sGravity, this);

  GenericTreeNode *lca;		// Least Common Ancestor
  // check whether we have a valid lca. for the first bucket
  // (currentBucket == 0) the lca must be set to the highest point
  // in the local tree that contains bucket 0 (i.e., root).
  lca = getStartAncestor(currentBucket, prevBucket, root);

  int lcaLevel = lca->getLevel(lca->getKey());

  // set opt
  sGravity->init(lca, activeRung, sLocal);
  // must set lowest node here
#ifdef CHANGA_REFACTOR_MEMCHECK
  CkPrintf("startNextBucket memcheck after init\n");
  CmiMemoryCheck();
#endif

#else
  GenericTreeNode *target = bucketList[currentBucket];

  sTopDown->init(sGravity, this);
  sGravity->init(target, activeRung, sLocal);
#endif

#if INTERLIST_VER > 0
    DoubleWalkState *lstate = (DoubleWalkState *)sLocalGravityState;
    if(!lstate->placedRoots[0]){
      lstate->placedRoots[0] = true;
#endif
      for(int cr = 0; cr < numChunks; cr++){
#ifdef DISABLE_NODE_TREE
        GenericTreeNode *chunkRoot = keyToNode(prefetchRoots[cr]);
#else
        GenericTreeNode *chunkRoot = dm->chunkRootToNode(prefetchRoots[cr]);
#endif
        if(chunkRoot != 0){
          for(int x = -nReplicas; x <= nReplicas; x++) {
            for(int y = -nReplicas; y <= nReplicas; y++) {
              for(int z = -nReplicas; z <= nReplicas; z++) {
#if INTERLIST_VER > 0
                // place chunk root on chklist of lcaLevel:
                OffsetNode on;
                on.node = chunkRoot;
                // value of currentBucket doesn't matter
                on.offsetID = encodeOffset(0, x,y,z);
                lstate->chklists[lcaLevel].enq(on);
#else
                // last -1 arg is the activeWalkIndex
                sTopDown->walk(chunkRoot, sLocalGravityState, -1, encodeOffset(currentBucket, x,y,z), -1);
#endif
              }
            }
          }
        }
      }
#if INTERLIST_VER > 0

#ifdef CHANGA_REFACTOR_MEMCHECK
      CkPrintf("startNextBucket memcheck after enq\n");
      CmiMemoryCheck();
#endif
    }
#endif
#if INTERLIST_VER > 0
    // all nodes to walk on have been enqueued, call walk
    // lca is the node the dft starts from.
    // sIntersListStateLocal contains the nodes to be walked on
    // last -1 arg is the activeWalkIndex
    // second last -1 arg is reqID. Because it doesn't make sense in this case,
    // we send the currentBucket number instead (as target)
    // third last -1 arg is chunk
    sInterListWalk->walk(lca, sLocalGravityState, -1, currentBucket, -1);
#ifdef CHANGA_REFACTOR_MEMCHECK
    CkPrintf("startNextBucket memcheck after walk\n");
    CmiMemoryCheck();
#endif

#endif
}

/*inline*/
void TreePiece::finishBucket(int iBucket) {
  BucketGravityRequest *req = &bucketReqs[iBucket];
  int remaining;

  remaining = sRemoteGravityState->counterArrays[0][iBucket]
              + sLocalGravityState->counterArrays[0][iBucket];

  CkAssert(remaining >= 0);
#ifdef COSMO_PRINT
  CkPrintf("[%d] Is finished %d? finished=%d, %d still missing!\n",thisIndex,iBucket,req->finished, remaining);
#endif

  // XXX finished means Ewald is done.
  if(req->finished && remaining == 0) {
    sLocalGravityState->myNumParticlesPending -= 1;

#ifdef COSMO_PRINT_BK
    CkPrintf("[%d] Finished bucket %d, %d particles remaining\n",thisIndex,iBucket, sLocalGravityState->myNumParticlesPending);
#endif

    if(sLocalGravityState->myNumParticlesPending == 0) {
      if(verbosity>1){
#if COSMO_STATS > 0
        CkPrintf("[%d] TreePiece %d finished with bucket %d , openCriterions:%lld\n",CkMyPe(),thisIndex,iBucket,numOpenCriterionCalls);
#else
        CkPrintf("[%d] TreePiece %d finished with bucket %d\n",CkMyPe(),thisIndex,iBucket);
#endif
      }

#if defined CUDA
      // in cuda version, must wait till particle accels.
      // are copied back from gpu; can markwalkdone only
      // after this, otherwise there is a race condition
      // between the start of the next iteration, wherein
      // the treepiece registers itself afresh with the 
      // data manager. if during this time the particles
      // haven't been copied (and updateParticles hasn't
      // been called), the registeredTreePieces list will
      // not be reset, so that the data manager gets 
      // confused.
      dm->transferParticleVarsBack();
      //dm->freeLocalTreeMemory();
#else
      // move on to markwalkdone in non-cuda version
      continueWrapUp();
#endif
    }
  }
}

void TreePiece::continueWrapUp(){
#ifdef CHECK_WALK_COMPLETIONS
  CkPrintf("[%d] markWalkDone TreePiece::continueWrapUp\n", thisIndex);
#endif

#ifdef CACHE_MEM_STATS
  memWithCache = CmiMemoryUsage()/(1024*1024);
#endif
  nNodeCacheEntries = cacheNode.ckLocalBranch()->getCache()->size();
  nPartCacheEntries = cacheGravPart.ckLocalBranch()->getCache()->size();

  markWalkDone();

  if(verbosity > 4){
    ckerr << "TreePiece " << thisIndex << ": My particles are done"
      << endl;
  }
}

void TreePiece::doAllBuckets(){
#if COSMO_DEBUG > 0
  char fout[100];
  sprintf(fout,"tree.%d.%d",thisIndex,iterationNo);
  ofstream ofs(fout);
  printTree(root,ofs);
  ofs.close();
  report();
#endif

  dummyMsg *msg = new (8*sizeof(int)) dummyMsg;
  *((int *)CkPriorityPtr(msg)) = numTreePieces * numChunks + thisIndex + 1;
  CkSetQueueing(msg,CK_QUEUEING_IFIFO);
  msg->val=0;

  thisProxy[thisIndex].nextBucket(msg);
#ifdef CUDA_INSTRUMENT_WRS
  ((DoubleWalkState *)sLocalGravityState)->nodeListConstructionTimeStart();
  ((DoubleWalkState *)sLocalGravityState)->partListConstructionTimeStart();
#endif
}

#ifdef CELL
#define CELLTHREASHOLD 100
#define CELLEWALDTHREASHOLD 30

inline void cellEnableComputation() {
  if (workRequestOut < CELLEWALDTHREASHOLD) {
    for (int i=0; i<ewaldMessages.length(); ++i) {
      CellComputation &comp = ewaldMessages[i];
      comp.owner.calculateEwald(comp.msg);
    }
    ewaldMessages.removeAll();
  }
}
#endif

void TreePiece::nextBucket(dummyMsg *msg){
  unsigned int i=0;

  int currentBucket = sLocalGravityState->currentBucket;
  if(verbosity >= 4)
	CkPrintf("[%d] walking bucket %d\n", thisIndex, currentBucket);	
#if INTERLIST_VER > 0
  sInterListWalk->init(sGravity, this);
#endif
  while(i<_yieldPeriod && currentBucket<numBuckets
#ifdef CELL
    && workRequestOut < CELLTHREASHOLD
#endif
  ){
    GenericTreeNode *target = bucketList[currentBucket];
    if(target->rungs >= activeRung){
#ifdef CHANGA_REFACTOR_INTERLIST_PRINT_BUCKET_START_FIN
      CkPrintf("[%d] local bucket active %d buckRem: %d + %d \n", thisIndex, currentBucket, sRemoteGravityState->counterArrays[0][currentBucket], sLocalGravityState->counterArrays[0][currentBucket]);
#endif
      // construct lists
      startNextBucket();
#if INTERLIST_VER > 0
      // do computation
      GenericTreeNode *lowestNode = ((DoubleWalkState *) sLocalGravityState)->lowestNode;
      int startBucket, end;
      int numActualBuckets = 0;

      getBucketsBeneathBounds(lowestNode, startBucket, end);

      CkAssert(currentBucket >= startBucket);

      sGravity->stateReady(sLocalGravityState, this, -1, currentBucket, end);
#ifdef CHANGA_REFACTOR_MEMCHECK
      CkPrintf("active: nextBucket memcheck after stateReady\n");
      CmiMemoryCheck();
#endif
      // book-keeping
#if COSMO_PRINT_BK > 1 
      CkPrintf("[%d] active local bucket book-keep\n", thisIndex);
#endif
      for(int j = currentBucket; j < end; j++){
#if !defined CELL 
        // for the cuda version, must decrement here
        // since we add to numAdditionalReqs for a bucket
        // each time a work request involving it is sent out.
        // the per bucket, per chunk counts are balanced as follows:
        // initialize in startiteration: +1
        // for request sent to gpu: +1
        // initial walk completed: -1 (this is where we are in the code currently) 
        // for each request completed: -1 (this happens in cudaCallback)
        sLocalGravityState->counterArrays[0][j]--;
#if COSMO_PRINT_BK > 1
        CkPrintf("[%d] bucket %d numAddReq: %d,%d\n", thisIndex, j, sRemoteGravityState->counterArrays[0][j], sLocalGravityState->counterArrays[0][j]);
#endif
        // have to call finishBucket here for inactive buckets.
#ifdef CUDA
        if(bucketList[j]->rungs < activeRung){
#endif
          finishBucket(j);
#ifdef CUDA
        }
#endif

#endif
        if(bucketList[j]->rungs >= activeRung){
          numActualBuckets++;
        }
      }
#ifdef CHANGA_REFACTOR_MEMCHECK
      CkPrintf("active: nextBucket memcheck after book-keeping (prev: %d, curr: %d, i: %d)\n", prevBucket, currentBucket, i);
      CmiMemoryCheck();
#endif
      // current target is active, so it will be
      // ok to mark it as prev active target for next
      // time
      prevBucket = currentBucket;
      currentBucket = end;
      if(end < numBuckets) // finishBucket() could clean up state if
			   // we are at the end
	  sLocalGravityState->currentBucket = end;
      i += numActualBuckets;
#else
      sLocalGravityState->counterArrays[0][currentBucket]--;
      finishBucket(currentBucket);

      currentBucket++;
      if(currentBucket < numBuckets) // state could be deleted in this case.
	  sLocalGravityState->currentBucket++;
      i++;
#endif
    }
    else{
#if INTERLIST_VER > 0
      // target was not active, so tree under lca has
      // not been processed. look for next active bucket
      // keep moving forward until an active bucket is reached
      // all (inactive) buckets encountered meanwhile are
      // finished
#if COSMO_PRINT_BK > 1 
      CkPrintf("[%d] inactive local bucket book-keep\n", thisIndex);
#endif
      while(currentBucket < numBuckets && bucketList[currentBucket]->rungs < activeRung){

        sLocalGravityState->counterArrays[0][currentBucket]--;
#if COSMO_PRINT_BK > 1
        CkPrintf("[%d] bucket %d numAddReq: %d,%d\n", thisIndex, currentBucket, sRemoteGravityState->counterArrays[0][currentBucket], sLocalGravityState->counterArrays[0][currentBucket]);
#endif
        finishBucket(currentBucket);
        currentBucket++;
	if(currentBucket < numBuckets) // state could be deleted in
				       // this case.
	    sLocalGravityState->currentBucket++;
      }
#ifdef CHANGA_REFACTOR_MEMCHECK
      CkPrintf("not active: nextBucket memcheck after book-keeping (prev: %d, curr: %d)\n", prevBucket, currentBucket);
      CmiMemoryCheck();
#endif
#else
      while(currentBucket < numBuckets && bucketList[currentBucket]->rungs < activeRung){

	sLocalGravityState->counterArrays[0][currentBucket]--;
        finishBucket(currentBucket);
        currentBucket++;
	if(currentBucket < numBuckets) // state could be deleted in
				       // this case.
	    sLocalGravityState->currentBucket++;
      }
#endif
      // i isn't incremented because we've skipped inactive buckets
    }// end else (target not active)
  }// end while

  if (currentBucket<numBuckets) {
    thisProxy[thisIndex].nextBucket(msg);
  } else {

#if INTERLIST_VER > 0 && defined CUDA
    // The local walk might have some outstanding computation requests 
    // at this point. Flush these lists to GPU
    DoubleWalkState *ds = (DoubleWalkState *)sLocalGravityState;
    ListCompute *lc = (ListCompute *)sGravity;

    if(lc && ds){
      if(ds->nodeLists.totalNumInteractions > 0){
        lc->sendNodeInteractionsToGpu(ds, this);
        lc->resetCudaNodeState(ds);
      }
      if(ds->particleLists.totalNumInteractions > 0){
        lc->sendPartInteractionsToGpu(ds, this);
        lc->resetCudaPartState(ds);
      }
    }
#endif

    delete msg;
  }
}

/// This function could be replaced by the doAllBuckets() call.
void TreePiece::calculateGravityLocal() {
  doAllBuckets();
}

#ifdef CELL
void cellSPE_ewald(void *data) {
  CellEwaldRequest *cgr = (CellEwaldRequest *)data;
  //CkPrintf("cellSPE_ewald %d\n", cgr->firstBucket);
  int i;
  free_aligned(cgr->roData);
  int offset = (cgr->numActiveData + 3) & ~0x3;
  for (i=0; i<cgr->numActiveData; ++i) {
    // copy the forces calculated to the particle's data
    GravityParticle *dest = cgr->particles[i];
    //CkPrintf("cellSPE_single part %d: %p, %f %f %f\n",i,dest,cr->activeData[i].treeAcceleration.x,cr->activeData[i].treeAcceleration.y,cr->activeData[i].treeAcceleration.z);
    dest->treeAcceleration.x = dest->treeAcceleration.x + cgr->woData[i+offset];
    dest->treeAcceleration.y = dest->treeAcceleration.y + cgr->woData[i+2*offset];
    dest->treeAcceleration.z = dest->treeAcceleration.z + cgr->woData[i+3*offset];
    dest->potential += cgr->woData[i];
  }
  for (i=cgr->firstBucket; i<=cgr->lastBucket; ++i) {
    cgr->tp->bucketReqs[i].finished = 1;
    cgr->tp->finishBucket(i);
  }
  free_aligned(cgr->woData);
  delete cgr->particles;
  workRequestOut --;
  delete cgr;
  cellEnableComputation();
}
#endif

void TreePiece::calculateEwald(dummyMsg *msg) {
#ifdef SPCUDA
  h_idata = (EwaldData*) malloc(sizeof(EwaldData)); 

  CkArrayIndex1D myIndex = CkArrayIndex1D(thisIndex); 
  EwaldHostMemorySetup(h_idata, myNumParticles, nEwhLoop); 
  EwaldGPU();
#else
  unsigned int i=0;
  while (i<_yieldPeriod && ewaldCurrentBucket < numBuckets
#ifdef CELL
	 && workRequestOut < CELLEWALDTHREASHOLD
#endif
	 ) {
#ifdef CELL_EWALD
    int activePart=0, indexActivePart=0;
    for (int k=bucketList[ewaldCurrentBucket]->firstParticle; k<=bucketList[ewaldCurrentBucket]->lastParticle; ++k) {
      if (myParticles[k].rung >= activeRung) activePart++;
    }
    GravityParticle **partList = new GravityParticle*[activePart];
    int outputSize = ROUNDUP_128(4*sizeof(cellSPEtype)*(activePart+3));
    int inputSize = ROUNDUP_128(sizeof(CellEwaldContainer)+nEwhLoop*sizeof(CellEWT)+3*sizeof(cellSPEtype)*(activePart+3));
    cellSPEtype *output = (cellSPEtype*)malloc_aligned(outputSize, 128);
    CellEwaldContainer *input = (CellEwaldContainer*)malloc_aligned(inputSize, 128);
    cellSPEtype *positionX = (cellSPEtype*)(((char*)input)+sizeof(CellEwaldContainer));
    cellSPEtype *positionY = (cellSPEtype*)(((char*)positionX)+((activePart+3)>>2)*(4*sizeof(cellSPEtype)));
    cellSPEtype *positionZ = (cellSPEtype*)(((char*)positionY)+((activePart+3)>>2)*(4*sizeof(cellSPEtype)));
    CellEWT *ewtIn = (CellEWT*)(((char*)positionZ)+((activePart+3)>>2)*(4*sizeof(cellSPEtype)));
    CellEwaldRequest *cr = new CellEwaldRequest(output, activePart, input, partList, this, ewaldCurrentBucket, ewaldCurrentBucket);
    for (int k=bucketList[ewaldCurrentBucket]->firstParticle; k<=bucketList[ewaldCurrentBucket]->lastParticle; ++k) {
      if (myParticles[k].rung >= activeRung) {
	positionX[indexActivePart] = myParticles[k].position.x;
	positionY[indexActivePart] = myParticles[k].position.y;
	positionZ[indexActivePart] = myParticles[k].position.z;
	partList[indexActivePart++] = &myParticles[k];
      }
    }
    input->rootMoments = root->moments;
    input->fEwCut = fEwCut;
    input->fPeriod = fPeriod.x;
    input->numPart = activePart;
    input->nReps = nReplicas;
    input->nEwhLoop = nEwhLoop;
    for (int k=0; k<nEwhLoop; ++k) {
      ewtIn[k] = ewt[k];
    }
    sendWorkRequest (3, NULL, 0, input, inputSize, output, outputSize, (void*)cr, 0, cellSPE_ewald, NULL);
    workRequestOut ++;
#else
    BucketEwald(bucketList[ewaldCurrentBucket], nReplicas, fEwCut);

    bucketReqs[ewaldCurrentBucket].finished = 1;
    finishBucket(ewaldCurrentBucket);
#endif

    ewaldCurrentBucket++;
    i++;
  }

  if (ewaldCurrentBucket<numBuckets) {
#ifdef CELL
    if (workRequestOut < CELLEWALDTHREASHOLD)
#endif
      thisProxy[thisIndex].calculateEwald(msg);
#ifdef CELL
    else
      ewaldMessages.insertAtEnd(CellComputation(thisProxy[thisIndex], msg));
#endif
  } else {
    delete msg;
  }
#endif
}

#if INTERLIST_VER > 0
/*
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
  */

#ifdef CELL
void cellSPE_callback(void *data) {
  //CkPrintf("cellSPE_callback\n");
  CellGroupRequest *cgr = (CellGroupRequest *)data;

  State *state = cgr->state;
  int bucket = cgr->bucket;

  state->counterArrays[0][bucket]--;
  cgr->tp->finishBucket(bucket);

  delete cgr->particles;
  delete cgr;
}

void cellSPE_single(void *data) {
  CellRequest *cr = (CellRequest *)data;
  free_aligned(cr->roData);
  for (int i=0; i<cr->numActiveData; ++i) {
    // copy the forces calculated to the particle's data
    GravityParticle *dest = cr->particles[i];
    //CkPrintf("cellSPE_single part %d: %p, %f %f %f\n",i,dest,cr->activeData[i].treeAcceleration.x,cr->activeData[i].treeAcceleration.y,cr->activeData[i].treeAcceleration.z);
    dest->treeAcceleration = dest->treeAcceleration + cr->activeData[i].treeAcceleration;
    dest->potential += cr->activeData[i].potential;
    if (cr->activeData[i].dtGrav > dest->dtGrav) {
      dest->dtGrav = cr->activeData[i].dtGrav;
    }
  }
  free_aligned(cr->activeData);
  workRequestOut --;
  delete cr;
  cellEnableComputation();
}
#endif

#if 0
#define CELLBUFFERSIZE 16*1024
void TreePiece::calculateForceLocalBucket(int bucketIndex){
    int cellListIter;
    int partListIter;

    /*for (int i=1; i<=myNumParticles; ++i) {
      if (myParticles[i].treeAcceleration.x != 0 ||
	  myParticles[i].treeAcceleration.y != 0 ||
	  myParticles[i].treeAcceleration.z != 0)
	CkPrintf("AAARRRGGGHHH (%d)!!! Particle %d (%d) on TP %d has acc: %f %f %f\n",bucketIndex,i,myParticles[i].iOrder,thisIndex,myParticles[i].treeAcceleration.x,myParticles[i].treeAcceleration.y,myParticles[i].treeAcceleration.z);
    }*/

    if (bucketList[bucketIndex]->rungs >= activeRung) {
#ifdef CELL
      //enableTrace();
      // Combine all the computation in a request to be sent to the cell SPE
      int activePart=0, indexActivePart=0;
      for (int k=bucketList[bucketIndex]->firstParticle; k<=bucketList[bucketIndex]->lastParticle; ++k) {
        if (myParticles[k].rung >= activeRung) activePart++;
      }
      int activePartDataSize = ROUNDUP_128((activePart+3)*sizeof(CellGravityParticle));
      CellGravityParticle *activePartData = (CellGravityParticle*)malloc_aligned(activePartDataSize,128);
      GravityParticle **partList = new GravityParticle*[activePart];
      CellGroupRequest *cgr = new CellGroupRequest(this, bucketIndex, partList);
      WRGroupHandle wrh = createWRGroup(cgr, cellSPE_callback);
      for (int k=bucketList[bucketIndex]->firstParticle; k<=bucketList[bucketIndex]->lastParticle; ++k) {
        if (myParticles[k].rung >= activeRung) {
          activePartData[indexActivePart] = myParticles[k];
	  //CkPrintf("[%d] particle %d: %d on bucket %d, acc: %f %f %f\n",thisIndex,k,myParticles[k].iOrder,bucketIndex,activePartData[indexActivePart].treeAcceleration.x,activePartData[indexActivePart].treeAcceleration.y,activePartData[indexActivePart].treeAcceleration.z);
	  //CkPrintf("[%d] partList %d: particle %d (%p)\n",thisIndex,indexActivePart,k,&myParticles[k]);
          partList[indexActivePart++] = &myParticles[k];
        }
      }
      /*int numPart=0, numNodes=0;
      for(int k=0;k<=curLevelLocal;k++) {
        numNodes += cellListLocal[k].length();
        for(int kk=0; kk<particleListLocal[k].length(); kk++) {
          numPart += particleListLocal[k][kk].numParticles;
        }
      }*/
      int particlesPerRequest = (CELLBUFFERSIZE - 2*sizeof(int)) / sizeof(CellExternalGravityParticle);
      int nodesPerRequest = (CELLBUFFERSIZE - 2*sizeof(int)) / sizeof(CellMultipoleMoments);
      CellContainer *particleContainer = (CellContainer*)malloc_aligned(CELLBUFFERSIZE,128);
      CellExternalGravityParticle *particleData = (CellExternalGravityParticle*)&particleContainer->data;
      CellContainer *nodesContainer = (CellContainer*)malloc_aligned(CELLBUFFERSIZE,128);
      CellMultipoleMoments *nodesData = (CellMultipoleMoments*)&nodesContainer->data;
      int indexNodes=0, indexPart=0;
#endif
      for(int k=0;k<=curLevelLocal;k++){
        int computed=0;
        double startTimer = CmiWallTimer();
#ifdef HPM_COUNTER
        hpmStart(1,"node force");
#endif
        for(cellListIter=0;cellListIter<cellListLocal[k].length();cellListIter++){
          OffsetNode tmp = cellListLocal[k][cellListIter];
#if COSMO_DEBUG > 1
          bucketcheckList[bucketIndex].insert(tmp->getKey());
          combineKeys(tmp->getKey(),bucketIndex);
#endif
#ifdef CELL_NODE
          nodesData[indexNodes] = tmp.node->moments;
          Vector3D<double> tmpoffsetID = decodeOffset(tmp.offsetID);
          nodesData[indexNodes].cm += tmpoffsetID;
          indexNodes++;
          if (indexNodes == nodesPerRequest) {
            // send off request
            void *activeData = malloc_aligned(activePartDataSize,128);
            memcpy(activeData, activePartData, activePart*sizeof(CellGravityParticle));
            CellRequest *userData = new CellRequest((CellGravityParticle*)activeData, activePart, nodesContainer, partList, this);
	    nodesContainer->numInt = activePart;
	    nodesContainer->numExt = indexNodes;
	    //CkPrintf("[%d] sending request 1 %p+%d, %p+%d\n",thisIndex,activeData, activePartDataSize, nodesContainer, CELLBUFFERSIZE);
            sendWorkRequest (1, activeData, activePartDataSize, nodesContainer, CELLBUFFERSIZE, NULL, 0,
                            (void*)userData, WORK_REQUEST_FLAGS_BOTH_CALLBACKS, cellSPE_single, wrh);
	    workRequestOut ++;
            nodeInterLocal += activePart * indexNodes;
	    nodesContainer = (CellContainer*)malloc_aligned(CELLBUFFERSIZE,128);
            nodesData = (CellMultipoleMoments*)&nodesContainer->data;
            indexNodes = 0;
          }
#else
          computed = nodeBucketForce(tmp.node, bucketList[bucketIndex], myParticles,
                                     decodeOffset(tmp.offsetID), activeRung);
#endif
        }
#ifdef CELL
        if (indexNodes > 0 && k==curLevelLocal) {
#ifdef CELL_NODE
          void *activeData = malloc_aligned(activePartDataSize,128);
          memcpy(activeData, activePartData, activePart*sizeof(CellGravityParticle));
          CellRequest *userData = new CellRequest((CellGravityParticle*)activeData, activePart, nodesContainer, partList, this);
	  nodesContainer->numInt = activePart;
	  nodesContainer->numExt = indexNodes;
          //CkPrintf("[%d] sending request 2 %p+%d, %p+%d\n",thisIndex,activeData, activePartDataSize, nodesContainer, ROUNDUP_128(indexNodes*sizeof(CellMultipoleMoments)));
          sendWorkRequest (1, activeData, activePartDataSize, nodesContainer, ROUNDUP_128(indexNodes*sizeof(CellMultipoleMoments)+2*sizeof(int)),
                           NULL, 0, (void*)userData, WORK_REQUEST_FLAGS_BOTH_CALLBACKS, cellSPE_single, wrh);
	  workRequestOut ++;
          nodeInterLocal += activePart * indexNodes;
#endif
        } else if (k==curLevelLocal) {
          free_aligned(nodesContainer);
        }
#endif
#ifdef HPM_COUNTER
        hpmStop(1);
#endif
        nodeInterLocal += cellListLocal[k].length() * computed;

        LocalPartInfo pinfo;
        double newTimer = CmiWallTimer();
        traceUserBracketEvent(nodeForceUE, startTimer, newTimer);
#ifdef HPM_COUNTER
    hpmStart(2,"particle force");
#endif
        for(partListIter=0;partListIter<particleListLocal[k].length();partListIter++){
          pinfo = particleListLocal[k][partListIter];
#if COSMO_DEBUG > 1
          bucketcheckList[bucketIndex].insert((pinfo.nd)->getKey());
          combineKeys((pinfo.nd)->getKey(),bucketIndex);
#endif
          for(int j = 0; j < pinfo.numParticles; ++j){
#ifdef CELL_PART
            particleData[indexPart] = pinfo.particles[j];
            particleData[indexPart].position += pinfo.offset;
            indexPart++;
            if (indexPart == particlesPerRequest) {
              // send off request
              void *activeData = malloc_aligned(activePartDataSize,128);
              memcpy(activeData, activePartData, activePart*sizeof(CellGravityParticle));
              CellRequest *userData = new CellRequest((CellGravityParticle*)activeData, activePart, particleContainer, partList, this);
	      particleContainer->numInt = activePart;
	      particleContainer->numExt = indexPart;
              //CkPrintf("[%d] sending request 3 %p+%d, %p+%d\n",thisIndex,activeData, activePartDataSize, particleContainer, CELLBUFFERSIZE);
              sendWorkRequest (2, activeData, activePartDataSize, particleContainer, CELLBUFFERSIZE, NULL, 0,
                              (void*)userData, WORK_REQUEST_FLAGS_BOTH_CALLBACKS, cellSPE_single, wrh);
	      workRequestOut ++;
	      particleInterLocal += activePart * indexPart;
	      particleContainer = (CellContainer*)malloc_aligned(CELLBUFFERSIZE,128);
              particleData = (CellExternalGravityParticle*)&particleContainer->data;
              indexPart = 0;
            }
#else
            computed = partBucketForce(&pinfo.particles[j], bucketList[bucketIndex], myParticles, pinfo.offset, activeRung);
#endif
          }
          particleInterLocal += pinfo.numParticles * computed;
        }
#ifdef CELL
        if (indexPart > 0 && k==curLevelLocal) {
#ifdef CELL_PART
          //void *activeData = malloc_aligned(activePartDataSize,128);
          //memcpy(activeData, activePartData, activePart*sizeof(GravityParticle));
          CellRequest *userData = new CellRequest(activePartData, activePart, particleContainer, partList, this);
	  particleContainer->numInt = activePart;
	  particleContainer->numExt = indexPart;
          //CkPrintf("[%d] sending request 4 %p+%d, %p+%d (%d int %d ext)\n",thisIndex,activePartData, activePartDataSize, particleContainer, ROUNDUP_128(indexPart*sizeof(CellExternalGravityParticle)),activePart,indexPart);
          sendWorkRequest (2, activePartData, activePartDataSize, particleContainer, ROUNDUP_128(indexPart*sizeof(CellExternalGravityParticle)+2*sizeof(int)),
                           NULL, 0, (void*)userData, WORK_REQUEST_FLAGS_BOTH_CALLBACKS, cellSPE_single, wrh);
	  workRequestOut ++;
          particleInterLocal += activePart * indexPart;
#endif
        } else if (k==curLevelLocal) {
          free_aligned(activePartData);
          free_aligned(particleContainer);
        }
#endif
#ifdef HPM_COUNTER
    hpmStop(2);
#endif
#ifdef COSMO_EVENTS
        traceUserBracketEvent(partForceUE, newTimer, CmiWallTimer());
#endif
      }
#ifdef CELL
      // Now all the requests have been made
      completeWRGroup(wrh);
      OffloadAPIProgress();
#endif
    }

#ifndef CELL
    finishBucket(bucketIndex);
#endif
}
#endif

#if 0
void TreePiece::calculateForceRemoteBucket(int bucketIndex, int chunk){

  int cellListIter;
  int partListIter;

  if (bucketList[bucketIndex]->rungs >= activeRung) {
    for(int k=0;k<=curLevelRemote;k++){
      int computed;
      double startTimer = CmiWallTimer();
#ifdef HPM_COUNTER
      hpmStart(1,"node force");
#endif
      for(cellListIter=0;cellListIter<cellList[k].length();cellListIter++){
        OffsetNode tmp= cellList[k][cellListIter];
#if COSMO_DEBUG > 1
        bucketcheckList[bucketIndex].insert(tmp.node->getKey());
        combineKeys(tmp.node->getKey(),bucketIndex);
#endif
        //computed = nodeBucketForce(tmp.node, bucketList[bucketIndex], myParticles,
        //                           decodeOffset(tmp.offsetID), activeRung);
      }
#ifdef HPM_COUNTER
      hpmStop(1);
#endif
      nodeInterRemote[chunk] += cellList[k].length() * computed;

      double newTimer = CmiWallTimer();
      traceUserBracketEvent(nodeForceUE, startTimer, newTimer);
#ifdef HPM_COUNTER
    hpmStart(2,"particle force");
#endif
      for(partListIter=0;partListIter<particleList[k].length();partListIter++){
        RemotePartInfo pinfo = particleList[k][partListIter];
#if COSMO_DEBUG > 1
        bucketcheckList[bucketIndex].insert((pinfo.nd)->getKey());
        combineKeys((pinfo.nd)->getKey(),bucketIndex);
#endif
        for(int j = 0; j < pinfo.numParticles; ++j){
          //computed = partBucketForce(&pinfo.particles[j], bucketList[bucketIndex],
          //                           myParticles, pinfo.offset, activeRung);
        }
        particleInterRemote[chunk] += pinfo.numParticles * computed;
      }
#ifdef HPM_COUNTER
    hpmStop(2);
#endif
#ifdef COSMO_EVENTS
      traceUserBracketEvent(partForceUE, newTimer, CmiWallTimer());
#endif
    }
  }

  finishBucket(bucketIndex);

}
#endif

#endif

const char *typeString(NodeType type);

void TreePiece::calculateGravityRemote(ComputeChunkMsg *msg) {
  unsigned int i=0;
  // cache internal tree: start directly asking the CacheManager
#ifdef DISABLE_NODE_TREE
  GenericTreeNode *chunkRoot = keyToNode(prefetchRoots[msg->chunkNum]);
#else
  GenericTreeNode *chunkRoot = dm->chunkRootToNode(prefetchRoots[msg->chunkNum]);
#endif

#ifdef CHANGA_REFACTOR_MEMCHECK
  CkPrintf("memcheck right after starting cgr\n");
  CmiMemoryCheck();
#endif

#ifdef CUDA_INSTRUMENT_WRS
  ((DoubleWalkState *)sRemoteGravityState)->nodeListConstructionTimeStart();
  ((DoubleWalkState *)sRemoteGravityState)->partListConstructionTimeStart();
  ((DoubleWalkState *)sInterListStateRemoteResume)->nodeListConstructionTimeStart();
  ((DoubleWalkState *)sInterListStateRemoteResume)->partListConstructionTimeStart();
#endif
  CkAssert(chunkRoot != NULL);

#if INTERLIST_VER > 0
  sInterListWalk->init(sGravity, this);
#else
  sTopDown->init(sGravity, this);
#endif

#ifdef CHANGA_REFACTOR_MEMCHECK
  CkPrintf("memcheck right after init\n");
  CmiMemoryCheck();
#endif


  while (i<_yieldPeriod && sRemoteGravityState->currentBucket < numBuckets
#ifdef CELL
		&& workRequestOut < CELLTHREASHOLD
#endif
  ) {
#ifdef CHANGA_REFACTOR_WALKCHECK
    if(thisIndex == CHECK_INDEX && sRemoteGravityState->currentBucket == CHECK_BUCKET){
      CkPrintf("Starting remote walk\n");
    }
#endif

#ifdef CHANGA_REFACTOR_INTERLIST_PRINT_BUCKET_START_FIN
    CkPrintf("[%d] remote bucket active %d buckRem: %d + %d, remChunk: %d\n", thisIndex, sRemoteGravityState->currentBucket, sRemoteGravityState->counterArrays[0][sRemoteGravityState->currentBucket], sLocalGravityState->counterArrays[0][sRemoteGravityState->currentBucket], sRemoteGravityState->counterArrays[1][msg->chunkNum]);
#endif
    // Interlist and normal versions both have 'target' nodes
    GenericTreeNode *target = bucketList[sRemoteGravityState->currentBucket];

#if INTERLIST_VER > 0
    GenericTreeNode *lca = 0;
#else
    sGravity->init(target, activeRung, sRemote);
#endif

#ifdef CHANGA_REFACTOR_MEMCHECK
    CkPrintf("memcheck in while loop (i: %d, currentRemote Bucket: %d)\n", i, sRemoteGravityState->currentBucket);
    CmiMemoryCheck();
#endif

    if (target->rungs >= activeRung) {
#if INTERLIST_VER > 0
      lca = getStartAncestor(sRemoteGravityState->currentBucket, prevRemoteBucket, root);

      int lcaLevel = lca->getLevel(lca->getKey());
      // set remote opt, etc
      sGravity->init(lca, activeRung, sRemote);

      // need to enqueue the chunkroot replicas only once,
      // before walking the first bucket
      DoubleWalkState *rstate = (DoubleWalkState *)sRemoteGravityState;
      if(!rstate->placedRoots[msg->chunkNum]){
        rstate->placedRoots[msg->chunkNum] = true;
        rstate->level = 0;
        sGravity->initState(rstate);
        
#if COSMO_PRINT_BK > 0
        CkPrintf("[%d] CGR: placing chunk %d root %ld (type %s) replicas\n", thisIndex, msg->chunkNum, chunkRoot->getKey(), typeString(chunkRoot->getType()));
#endif

#endif

        for(int x = -nReplicas; x <= nReplicas; x++) {
          for(int y = -nReplicas; y <= nReplicas; y++) {
    	    for(int z = -nReplicas; z <= nReplicas; z++) {
#if CHANGA_REFACTOR_DEBUG > 1
    	      CkPrintf("[%d]: starting remote walk with chunk=%d, current remote bucket=%d, (%d,%d,%d)\n", thisIndex, msg->chunkNum, sRemoteGravityState->currentBucket, x, y, z);
#endif
#if INTERLIST_VER > 0
              // put chunk root in correct checklist
              // the currentRemote Bucket argument to encodeOffset doesn't really
              // matter; only the offset is of consequence
              OffsetNode on;
              on.node = chunkRoot;
              on.offsetID = encodeOffset(0, x,y,z);
              rstate->chklists[lcaLevel].enq(on);

#else
    	      sTopDown->walk(chunkRoot, sRemoteGravityState, msg->chunkNum,encodeOffset(sRemoteGravityState->currentBucket,x, y, z), remoteGravityAwi);
#endif
    	    }
          }
        }
#if INTERLIST_VER > 0
#ifdef CHANGA_REFACTOR_MEMCHECK
        CkPrintf("active: memcheck after enqueuing (i: %d, currentRemote Bucket: %d, lca: %ld)\n", i, sRemoteGravityState->currentBucket, lca->getKey());
        CmiMemoryCheck();
#endif
      }// if !placedRoots

      // now that all the nodes to walk on have been enqueued in the chklist of the lca, we can call
      // the walk function.
      // lca is the node the dft starts from.
      // sRemoteGravityState contains the nodes to be walked on
      // target for walk has already been set
      // the offset doesn't really matter, so send target bucket number (currentBucketRemote), which is needed to construct the reqID
      // in case nodes are missed and walks need to be resumed.
      sInterListWalk->walk(lca, sRemoteGravityState, msg->chunkNum, sRemoteGravityState->currentBucket, interListAwi);
#ifdef CHANGA_REFACTOR_MEMCHECK
      CkPrintf("active: memcheck after walk\n");
      CmiMemoryCheck();
#endif
      // now that the walk has been completed, do some book-keeping

      // keeping track of the fact that a chunk has been processed by some buckets
      // and that all particles under lowest node have used this chunk.
      GenericTreeNode *lowestNode = rstate->lowestNode;
      int startBucket, end;
      int numActualBuckets = 0;

      getBucketsBeneathBounds(lowestNode, startBucket, end);
      CkAssert(sRemoteGravityState->currentBucket >= startBucket);

      sGravity->stateReady(sRemoteGravityState, this, msg->chunkNum, sRemoteGravityState->currentBucket, end);
#ifdef CHANGA_REFACTOR_MEMCHECK
      CkPrintf("active: memcheck after stateReady\n");
      CmiMemoryCheck();
#endif

#if COSMO_PRINT_BK > 1 
      CkPrintf("[%d] active remote bucket book-keep current: %d end: %d chunk: %d\n", thisIndex, sRemoteGravityState->currentBucket, end, msg->chunkNum);
#endif
      for(int j = sRemoteGravityState->currentBucket; j < end; j++){
#if !defined CELL 
        // see comment in nextBucket for the reason why
        // this decrement is performed for the cuda version
        sRemoteGravityState->counterArrays[0][j]--;
#if COSMO_PRINT_BK > 1
        CkPrintf("[%d] bucket %d numAddReq: %d,%d\n", thisIndex, j, sRemoteGravityState->counterArrays[0][j], sLocalGravityState->counterArrays[0][j]);
#endif

#ifdef CUDA
        if(bucketList[j]->rungs < activeRung){
#endif
          finishBucket(j);
#ifdef CUDA
        }
#endif

#endif
        if(bucketList[j]->rungs >= activeRung){
          numActualBuckets++;
        }
      }
      sRemoteGravityState->counterArrays[1][msg->chunkNum] -= (end-sRemoteGravityState->currentBucket);
#ifdef CHANGA_REFACTOR_MEMCHECK
      CkPrintf("active: memcheck after book-keeping\n");
      CmiMemoryCheck();
#endif

      prevRemoteBucket = sRemoteGravityState->currentBucket;
      sRemoteGravityState->currentBucket = end;
      i += numActualBuckets;

#endif
    }// end if target active
#if INTERLIST_VER > 0
    else{
      // target was not active, so tree under lca has
      // not been processed. look for next active bucket
      // keep moving forward until an active bucket is reached
      // all (inactive) buckets encountered meanwhile are
      // finished
      int bucketsSkipped = 0;
      // prevRemoteBucket keeps track of last *active* bucket
      // if it is -1, the lca turns out to be chunkroot
      // otherwise, the lca is some other valid node
#if COSMO_PRINT_BK > 1 
      CkPrintf("[%d] inactive remote bucket book-keep chunk: %d\n", thisIndex, msg->chunkNum);
#endif
      while(sRemoteGravityState->currentBucket < numBuckets && bucketList[sRemoteGravityState->currentBucket]->rungs < activeRung){
        sRemoteGravityState->counterArrays[0][sRemoteGravityState->currentBucket]--;
#if COSMO_PRINT_BK > 1
        CkPrintf("[%d] bucket %d numAddReq: %d,%d\n", thisIndex, sRemoteGravityState->currentBucket, sRemoteGravityState->counterArrays[0][sRemoteGravityState->currentBucket], sLocalGravityState->counterArrays[0][sRemoteGravityState->currentBucket]);
#endif
        finishBucket(sRemoteGravityState->currentBucket);
        bucketsSkipped++;
        sRemoteGravityState->currentBucket++;
      }
      sRemoteGravityState->counterArrays[1][msg->chunkNum] -= bucketsSkipped;
#ifdef CHANGA_REFACTOR_MEMCHECK
      CkPrintf("inactive: memcheck after book-keeping (prevRemote: %d, current: %d, chunk: %d)\n", prevRemoteBucket, sRemoteGravityState->currentBucket, msg->chunkNum);
      CmiMemoryCheck();
#endif
      // i isn't incremented because we've skipped inactive buckets
    }
#endif


#if INTERLIST_VER > 0
#else
    sRemoteGravityState->counterArrays[0][sRemoteGravityState->currentBucket]--;
    finishBucket(sRemoteGravityState->currentBucket);
    sRemoteGravityState->counterArrays[1][msg->chunkNum] --;
    sRemoteGravityState->currentBucket++;
    i++;
#endif
  }// end while i < yieldPeriod and currentRemote Bucket < numBuckets


  if (sRemoteGravityState->currentBucket < numBuckets) {
    thisProxy[thisIndex].calculateGravityRemote(msg);
#if COSMO_PRINT > 0
    CkPrintf("{%d} sending self-message chunk %d, prio %d\n",thisIndex,msg->chunkNum,*(int*)CkPriorityPtr(msg));
#endif
  } else {
    sRemoteGravityState->currentBucket = 0;

#if INTERLIST_VER > 0 && defined CUDA
    // The remote walk might have some outstanding computation requests 
    // at this point. Flush these lists to GPU
    DoubleWalkState *ds = (DoubleWalkState *)sRemoteGravityState;
    ListCompute *lc = (ListCompute *)sGravity;

    if(lc && ds){
      if(ds->nodeLists.totalNumInteractions > 0){
        lc->sendNodeInteractionsToGpu(ds, this);
        lc->resetCudaNodeState(ds);
      }
      if(ds->particleLists.totalNumInteractions > 0){
        lc->sendPartInteractionsToGpu(ds, this);
        lc->resetCudaPartState(ds);
      }
    }

#endif

    int chunkRemaining;
    chunkRemaining = sRemoteGravityState->counterArrays[1][msg->chunkNum];
#if INTERLIST_VER > 0
    prevRemoteBucket = -1;
#endif
    CkAssert(chunkRemaining >= 0);
#if COSMO_PRINT_BK > 1
    CkPrintf("[%d] cgr chunk: %d remaining Chunk: %d\n", thisIndex, msg->chunkNum, chunkRemaining);
#endif
    if (chunkRemaining == 0) {
      // we finished completely using this chunk, so we acknowledge the cache
      // if this is not true it means we had some hard misses
#ifdef COSMO_PRINT_BK
      CkPrintf("[%d] FINISHED CHUNK %d from calculateGravityRemote\n",thisIndex,msg->chunkNum);
#endif

#ifdef CUDA
      // The remote-resume walk might have some outstanding computation requests 
      // at this point. Flush these lists to GPU
      DoubleWalkState *ds = (DoubleWalkState *)sInterListStateRemoteResume;
      ListCompute *lc = (ListCompute *)sGravity;


      if(lc && ds){
        bool nodeDummy = false;
        bool partDummy = true;

        if(ds->nodeLists.totalNumInteractions > 0){
          nodeDummy = true;
        }
        else if(ds->particleLists.totalNumInteractions > 0){
          partDummy = true;
        }

        if(nodeDummy && partDummy){
          lc->sendNodeInteractionsToGpu(ds, this);
          lc->sendPartInteractionsToGpu(ds, this, true);
          lc->resetCudaNodeState(ds);
          lc->resetCudaPartState(ds);
        }
        else if(nodeDummy){
          lc->sendNodeInteractionsToGpu(ds, this,true);
          lc->resetCudaNodeState(ds);
        }
        else if(partDummy){
          lc->sendPartInteractionsToGpu(ds, this, true);
          lc->resetCudaPartState(ds);
        }
      }
#endif

      cacheGravPart[CkMyPe()].finishedChunk(msg->chunkNum, particleInterRemote[msg->chunkNum]);
#ifdef CHECK_WALK_COMPLETIONS
      CkPrintf("[%d] finishedChunk TreePiece::calculateGravityRemote\n", thisIndex);
#endif
      finishedChunk(msg->chunkNum);
    }
#if COSMO_PRINT > 0
    CkPrintf("{%d} resetting message chunk %d, prio %d\n",thisIndex,msg->chunkNum,*(int*)CkPriorityPtr(msg));
#endif
    delete msg;
  }
}

#ifdef CUDA
void TreePiece::callFreeRemoteChunkMemory(int chunk){
  dm->freeRemoteChunkMemory(chunk);
}
#endif

#if INTERLIST_VER > 0
/// @brief return the largest node which contains current bucket, but
/// which does not contain previous bucket.  If previous is -1 then
/// return the root.
GenericTreeNode *TreePiece::getStartAncestor(int current, int previous, GenericTreeNode *chunkroot){
  GenericTreeNode *lca = 0;

  if(previous == -1){
    return chunkroot;
  }

  GenericTreeNode *target = bucketList[current];
  GenericTreeNode *prevTarget = bucketList[previous];
  Tree::NodeKey targetKey = target->getKey();
  Tree::NodeKey prevTargetKey = prevTarget->getKey();

  // get LCA key
  Tree::NodeKey lcaKey = target->getLongestCommonPrefix(targetKey, prevTargetKey);
  // get LCA node
  lca = keyToNode(lcaKey);
  int whichChild = lca->whichChild(targetKey);
  return lca->getChildren(whichChild);
}
#endif

// We are done with the node Cache

void TreePiece::finishNodeCache(const CkCallback& cb)
{
    int j;
    for (j = 0; j < numChunks; j++) {
	cacheNode.ckLocalBranch()->finishedChunk(j, 0);
	}
    contribute(cb);
    }

#ifdef PUSH_GRAVITY
/*
  This method is intended to calculate forces in the 'small' timesteps,
  i.e. when few particles are active. It causes the TreePiece to broadcast its
  buckets to all other TreePieces. The others compute forces on its buckets
  and contribute forces a reduction for this TreePiece in particular. When all
  work has finished, quiescence is detected and we move on in the small step.
*/

void TreePiece::startPushGravity(int am, double myTheta){
  LBTurnInstrumentOn();
  
  iterationNo++;
  activeRung = am;
  theta = myTheta;
  thetaMono = theta*theta*theta*theta;

  CkAssert(!doMerge);
  if(!createdSpanningTree){
    createdSpanningTree = true;
    allTreePieceSection = CProxySection_TreePiece::ckNew(thisProxy,0,numTreePieces-1,1);
    CkMulticastMgr *mgr = CProxy_CkMulticastMgr(ckMulticastGrpId).ckLocalBranch();
    allTreePieceSection.ckSectionDelegate(mgr);
  }

  BucketMsg *msg = createBucketMsg();
  if(msg != NULL) allTreePieceSection.recvPushBuckets(msg);
}

BucketMsg *TreePiece::createBucketMsg(){
  int saveNumActiveParticles = 0;
  int numActiveParticles = 0;
  int numActiveBuckets = 0;

  // First count the number of active particles and buckets
  for(int i = 0; i < bucketList.size(); i++){
    GenericTreeNode *bucket = bucketList[i];
    int buckStart = bucket->firstParticle; 
    int buckEnd = bucket->lastParticle;
    GravityParticle *buckParticles = bucket->particlePointer;
    saveNumActiveParticles = numActiveParticles;
    for(int j = 0; j <= buckEnd-buckStart; j++){
      if(buckParticles[j].rung >= activeRung){
        numActiveParticles++;
      }
    }
    if(numActiveParticles > saveNumActiveParticles) numActiveBuckets++;
  }

  //CkPrintf("tree %d active particles %d active buckets %d\n", thisIndex, numActiveParticles, numActiveBuckets);

  if(numActiveParticles == 0) return NULL;

  // allocate message
  BucketMsg *msg = new (numActiveBuckets,numActiveParticles) BucketMsg;

  numActiveParticles = 0;
  numActiveBuckets = 0;

  // Copy active particles and buckets into message; change pointer offsets
  // of bucket particle boundaries to integers
  for(int i = 0; i < bucketList.size(); i++){
    // source bucket
    GenericTreeNode *sbucket = bucketList[i];
    int buckStart = sbucket->firstParticle; 
    int buckEnd = sbucket->lastParticle;
    GravityParticle *buckParticles = sbucket->particlePointer;
    saveNumActiveParticles = numActiveParticles;
    for(int j = 0; j <= buckEnd-buckStart; j++){
      if(buckParticles[j].rung >= activeRung){
        // copy active particle to bucket msg
        msg->particles[numActiveParticles] = buckParticles[j]; 
        numActiveParticles++;
      }
    }
    if(numActiveParticles > saveNumActiveParticles){
      // copy active bucket to bucket msg
      msg->buckets[numActiveBuckets] = *sbucket;
      //CkPrintf("tree piece %d bucket %llu %llu active\n", thisIndex, msg->buckets[numActiveBuckets].getKey(), sbucket->getKey());
      GenericTreeNode &tbucket = msg->buckets[numActiveBuckets];
      // set particle bounds for copied bucket (as integers)
      tbucket.particlePointer = NULL;
      tbucket.firstParticle = saveNumActiveParticles;
      tbucket.lastParticle = numActiveParticles-1;
      numActiveBuckets++;
    }
  }

  msg->numBuckets = numActiveBuckets;
  msg->numParticles = numActiveParticles;
  msg->whichTreePiece = thisIndex;

  return msg;
}

void TreePiece::recvPushBuckets(BucketMsg *msg){
  GenericTreeNode *foreignBuckets;
  int numForeignBuckets;


  int numFields = 4;
  // make sure there is enough space for foreignParticles
  foreignParticles.resize(msg->numParticles);
  foreignParticleAccelerations.resize(numFields*msg->numParticles);
  // obtain positions of foreignParticles from message
  unpackBuckets(msg,foreignBuckets,numForeignBuckets);
  if(myNumParticles > 0){
    // If there is a local tree associated with this tree piece,
    // calculate forces on foreignParticles due to it
    calculateForces(foreignBuckets,numForeignBuckets);
  }
  // update cookie
  CkGetSectionInfo(cookieJar[msg->whichTreePiece],msg);

  for(int i = 0; i < msg->numParticles; i++){
    foreignParticleAccelerations[numFields*i] = foreignParticles[i].treeAcceleration.x;
    foreignParticleAccelerations[numFields*i+1] = foreignParticles[i].treeAcceleration.y;
    foreignParticleAccelerations[numFields*i+2] = foreignParticles[i].treeAcceleration.z;
    foreignParticleAccelerations[numFields*i+3] = foreignParticles[i].interMass;
  }

  // contribute accelerations
  CkCallback cb(CkIndex_TreePiece::recvPushAccelerations(NULL),CkArrayIndex1D(msg->whichTreePiece),thisProxy);
  CkMulticastMgr *mgr = CProxy_CkMulticastMgr(ckMulticastGrpId).ckLocalBranch();
  mgr->contribute(foreignParticleAccelerations.length()*sizeof(double),&foreignParticleAccelerations[0],CkReduction::sum_double,cookieJar[msg->whichTreePiece],cb);

  delete msg;
}

void TreePiece::unpackBuckets(BucketMsg *msg, GenericTreeNode *&foreignBuckets, int &numForeignBuckets){
  // Copy foreign particle positions, etc. into local buffer
  for(int i = 0; i < msg->numParticles; i++){
    foreignParticles[i] = msg->particles[i];
    foreignParticles[i].treeAcceleration.x = 0.0;
    foreignParticles[i].treeAcceleration.y = 0.0;
    foreignParticles[i].treeAcceleration.z = 0.0;
    foreignParticles[i].interMass = 0.0;
  }

  // Make buckets point to appropriate positions in local buffer of particles
  foreignBuckets = msg->buckets;
  numForeignBuckets = msg->numBuckets;

  //CkPrintf("tree %d recv %d particles %d buckets from %d\n", thisIndex, msg->numParticles, msg->numBuckets, msg->whichTreePiece);

  GravityParticle *baseParticlePtr = &foreignParticles[0];
  for(int i = 0; i < numForeignBuckets; i++){
    GenericTreeNode &bucket = foreignBuckets[i];
    bucket.particlePointer = baseParticlePtr+bucket.firstParticle;
  }
}

void TreePiece::calculateForces(GenericTreeNode *foreignBuckets, int numForeignBuckets){
  TopDownTreeWalk topdown;
  GravityCompute grav;
  NullState nullState;
  PushGravityOpt pushOpt;

  grav.init(NULL,activeRung,&pushOpt);

  CkAssert(root != NULL);
  for(int i = 0; i < numForeignBuckets; i++){
    GenericTreeNode &target = foreignBuckets[i];
    //CkPrintf("tree piece %d calculate forces on bucket %llu\n", thisIndex, target.getKey());
    grav.setComputeEntity(&target);
    topdown.init(&grav,this);
    // for each replica
    for(int x = -nReplicas; x <= nReplicas; x++){
      for(int y = -nReplicas; y <= nReplicas; y++){
        for(int z = -nReplicas; z <= nReplicas; z++){
          // begin walk at root
          // -1 for chunk and active walk index
          // bucket number 'i' doesn't serve any purpose,
          // since this traversal will not generate any remote requests.
          topdown.walk(root,&nullState,-1,encodeOffset(i,x,y,z),-1);
        }
      }
    }
  }
}

void TreePiece::recvPushAccelerations(CkReductionMsg *msg){
  double *accelerations = (double *) msg->getData();
  int numAccelerations = msg->getSize()/sizeof(double);
  int j = 0;

  int numUpdates = 0;
  int numFields = 4;
  for(int i = 1; i <= myNumParticles; i++){
    if(myParticles[i].rung >= activeRung){ 
      myParticles[i].treeAcceleration.x = accelerations[j];
      myParticles[i].treeAcceleration.y = accelerations[j+1];
      myParticles[i].treeAcceleration.z = accelerations[j+2];

      myParticles[i].interMass = accelerations[j+3]; 
      j += numFields;
      numUpdates++;

      double totalMass = myParticles[i].mass+myParticles[i].interMass;
      if(totalMass != myTotalMass){
        CkPrintf("[%d] particle %d interMass %f should be %f partMass %f\n", thisIndex, i, totalMass, myTotalMass, myParticles[i].mass);
        CkAbort("bad intermass\n");
      }
    }
  }
  CkAssert(numUpdates == numAccelerations/numFields);
}
#endif

void TreePiece::findTotalMass(CkCallback &cb){
  callback = cb;
  myTotalMass = 0;
  for(int i = 1; i <= myNumParticles; i++){
    myTotalMass += myParticles[i].mass;
  }

  contribute(sizeof(double), &myTotalMass, CkReduction::sum_double, CkCallback(CkIndex_TreePiece::recvTotalMass(NULL),thisProxy));
}

void TreePiece::recvTotalMass(CkReductionMsg *msg){
  myTotalMass = *((double *)msg->getData());
  contribute(callback);
}

/// This method starts the tree walk and gravity calculation.  It
/// first registers with the node and particle caches.  It initializes
/// the particle acceleration by calling initBucket().  It initializes
/// treewalk bookkeeping, and starts a prefetch walk which will
/// eventually call a remote walk (a walk on non-local nodes).  It
/// then starts the Ewald initialization, and finally starts the local
/// gravity walk.

void TreePiece::startGravity(int am, // the active mask for multistepping
			       double myTheta, // opening criterion
			       const CkCallback& cb) {
  LBTurnInstrumentOn();
  iterationNo++;

  cbGravity = cb;
  activeRung = am;
  theta = myTheta;
  thetaMono = theta*theta*theta*theta;

  int oldNumChunks = numChunks;
  dm->getChunks(numChunks, prefetchRoots);
  CkArrayIndexMax idxMax = CkArrayIndex1D(thisIndex);
  // The following if is necessary to make nodes containing only TreePieces
  // without particles to get stuck and crash...
  if (numChunks == 0 && myNumParticles == 0) numChunks = 1;
  int dummy;

  cacheNode.ckLocalBranch()->cacheSync(numChunks, idxMax, localIndex);
  cacheGravPart.ckLocalBranch()->cacheSync(numChunks, idxMax, dummy);

  if (myNumParticles == 0) {
    // No particles assigned to this TreePiece
    for (int i=0; i< numChunks; ++i) {
      cacheGravPart.ckLocalBranch()->finishedChunk(i, 0);
    }
    CkCallback cbf = CkCallback(CkIndex_TreePiece::finishWalk(), pieces);
    gravityProxy[thisIndex].ckLocal()->contribute(cbf);
    return;
  }
  
  // allocate and zero out statistics counters
  if (oldNumChunks != numChunks ) {
    delete[] nodeInterRemote;
    delete[] particleInterRemote;
    nodeInterRemote = new u_int64_t[numChunks];
    particleInterRemote = new u_int64_t[numChunks];
  }

  if(nodeInterRemote == NULL)
	nodeInterRemote = new u_int64_t[numChunks];
  if(particleInterRemote == NULL)
	particleInterRemote = new u_int64_t[numChunks];
#if COSMO_STATS > 0
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

  if(verbosity>1)
    CkPrintf("Node: %d, TreePiece %d: I have %d buckets\n", CkMyNode(),
    	     thisIndex,numBuckets);

  if (bucketReqs==NULL) bucketReqs = new BucketGravityRequest[numBuckets];


  ewaldCurrentBucket = 0;

#if INTERLIST_VER > 0
  prevBucket = -1;
  prevRemoteBucket = -1;
#endif

  initBuckets();

  switch(domainDecomposition){
    case Oct_dec:
    case ORB_dec:
    case ORB_space_dec:
      //Prefetch Roots for Oct
      prefetchReq[0].reset();
      for (unsigned int i=1; i<=myNumParticles; ++i) {
        if (myParticles[i].rung >= activeRung) {
          prefetchReq[0].grow(myParticles[i].position);
        }
      }
      break;
    default:
      //Prefetch Roots for SFC
      prefetchReq[0].reset();
      for (unsigned int i=1; i<=myNumParticles; ++i) {
	  // set to first active particle
        if (myParticles[i].rung >= activeRung) {
          prefetchReq[0].grow(myParticles[i].position);
	  break;
	}
      }
      prefetchReq[1].reset();
      for (unsigned int i=myNumParticles; i>=1; --i) {
	  // set to last active particle
        if (myParticles[i].rung >= activeRung) {
	    prefetchReq[1].grow(myParticles[i].position);
	    break;
	}
      }

      break;
  }

#if CHANGA_REFACTOR_DEBUG > 0
  CkPrintf("Beginning prefetch\n");
#endif

  sTopDown = new TopDownTreeWalk;

  sLocal = new LocalOpt;
  sRemote = new RemoteOpt;

  if(_prefetch) sPrefetch = new PrefetchCompute;
  else sPrefetch = new DummyPrefetchCompute;

  sPref = new PrefetchOpt;

  State *remoteWalkState;
  State *localWalkState;
  Compute *compute;
  TreeWalk *walk;

#if INTERLIST_VER > 0
  compute = new ListCompute;
  walk = new LocalTargetWalk;
#else
  compute = new GravityCompute;
#endif

  // much of this part is common to interlist and normal algo.
  compute->init((void *)0, activeRung, sLocal);
  localWalkState = compute->getNewState(numBuckets);

  compute->init((void *)0, activeRung, sRemote);
  remoteWalkState = compute->getNewState(numBuckets, numChunks);



#if INTERLIST_VER > 0
  // interaction list algo needs a separate state for walk resumption
  sInterListStateRemoteResume = compute->getNewState();
  // remote resume and remote share counters
  // but have separate lists
  sInterListStateRemoteResume->counterArrays[0] = remoteWalkState->counterArrays[0];
  sInterListStateRemoteResume->counterArrays[1] = remoteWalkState->counterArrays[1];
#endif

  
  // remaining Chunk[]
  for(int i = 0; i < numChunks; i++) {
    remoteWalkState->counterArrays[1][i] = numBuckets;
    //CkPrintf("[%d] chunk %d init remaining: %d\n", thisIndex, i, remoteWalkState->counterArrays[1][i]);
  }

  for(int i = 0; i < numBuckets; i++){
    remoteWalkState->counterArrays[0][i] = numChunks;
    localWalkState->counterArrays[0][i] = 1;
  }


  sGravity = compute;
  sLocalGravityState = localWalkState;
  sRemoteGravityState = remoteWalkState;
#if INTERLIST_VER > 0
  sInterListWalk = walk;

  // cuda - need number of active buckets this iteration so 
  // that we can flush accumulated interactions to gpu when
  // we realize no more computation requests will be received
  // (i.e. all buckets have finished their RNR/Local walks)
#ifdef CUDA

#ifdef CUDA_INSTRUMENT_WRS
  instrumentId = dm->initInstrumentation();

  localNodeListConstructionTime = 0.0;
  remoteNodeListConstructionTime = 0.0;
  remoteResumeNodeListConstructionTime = 0.0;
  localPartListConstructionTime = 0.0;
  remotePartListConstructionTime = 0.0;
  remoteResumePartListConstructionTime = 0.0;

  nLocalNodeReqs = 0;
  nRemoteNodeReqs = 0;
  nRemoteResumeNodeReqs = 0;
  nLocalPartReqs = 0;
  nRemotePartReqs = 0;
  nRemoteResumePartReqs = 0;
#endif

  numActiveBuckets = 0;
  calculateNumActiveParticles();

  for(int i = 0; i < numBuckets; i++){
    if(bucketList[i]->rungs >= activeRung){
      numActiveBuckets++;
    }
  }
  if(numActiveBuckets > 0 && verbosity > 1){
    CkPrintf("[%d] num active buckets %d avg size: %f\n", thisIndex,
	     numActiveBuckets, 1.0*myNumActiveParticles/numActiveBuckets);
  }


  {
	  DoubleWalkState *state = (DoubleWalkState *)sRemoteGravityState;
	  ((ListCompute *)sGravity)->initCudaState(state, numBuckets, remoteNodesPerReq, remotePartsPerReq, false);
          // no missed nodes/particles
          state->nodes = NULL;
          state->particles = NULL;

	  DoubleWalkState *lstate = (DoubleWalkState *)sLocalGravityState;
	  ((ListCompute *)sGravity)->initCudaState(lstate, numBuckets, localNodesPerReq, localPartsPerReq, false);
          // ditto
          lstate->nodes = NULL;
          // allocate space for local particles 
          // if this is a small phase; we do not transfer
          // all particles owned by this processor in small
          // phases and so refer to particles inside an 
          // auxiliary array shipped with computation requests.
          // this auxiliary array is 'particles', below:
          if(largePhase()){
            lstate->particles = NULL;
          }     
          else{
            // allocate an amount of space that 
            // depends on the rung
            lstate->particles = new CkVec<CompactPartData>(AVG_SOURCE_PARTICLES_PER_ACTIVE*myNumActiveParticles);
            lstate->particles->length() = 0;
            // need to allocate memory for data structure that stores bucket
            // active info (where in the gpu's target particle memory this
            // bucket starts, and its size; strictly speaking, we don't need
            // the size attribute.)
            // XXX - no need to allocate/delete every iteration
            bucketActiveInfo = new BucketActiveInfo[numBuckets];
          }
  }
  {
	  DoubleWalkState *state = (DoubleWalkState *)sInterListStateRemoteResume;
	  ((ListCompute *)sGravity)->initCudaState(state, numBuckets, remoteResumeNodesPerReq, remoteResumePartsPerReq, true);

	  state->nodes = new CkVec<CudaMultipoleMoments>(100000);
          state->nodes->length() = 0;
	  state->particles = new CkVec<CompactPartData>(200000);
          state->particles->length() = 0;
          //state->nodeMap.clear();
          state->nodeMap.resize(remoteResumeNodesPerReq);
          state->nodeMap.length() = 0;
          
          state->partMap.clear();
  }
#endif // CUDA

#if CUDA_STATS
  localNodeInteractions = 0;
  localPartInteractions = 0;
  remoteNodeInteractions = 0;
  remotePartInteractions = 0;
  remoteResumeNodeInteractions = 0;
  remoteResumePartInteractions = 0;
#endif


#endif // INTERLIST_VER

  sLocalGravityState->myNumParticlesPending = numBuckets; 
  sRemoteGravityState->numPendingChunks = numChunks;
  // current remote bucket
  sRemoteGravityState->currentBucket = 0;


  sPrefetch->init((void *)0, activeRung, sPref);
  sTopDown->init(sPrefetch, this);
  sPrefetchState = sPrefetch->getNewState(1);
  // instead of prefetchWaiting, we count through state->counters[0]
  sPrefetchState->counterArrays[0][0] = (2*nReplicas + 1)*(2*nReplicas + 1)*(2*nReplicas + 1);

  activeWalks.reserve(maxAwi);
  addActiveWalk(prefetchAwi, sTopDown,sPrefetch,sPref,sPrefetchState);
#ifdef CHECK_WALK_COMPLETIONS
  CkPrintf("[%d] addActiveWalk prefetch (%d)\n", thisIndex, activeWalks.length());
#endif

#if INTERLIST_VER > 0
  addActiveWalk(interListAwi, sInterListWalk, sGravity, sRemote,
		sInterListStateRemoteResume);
#ifdef CHECK_WALK_COMPLETIONS
  CkPrintf("[%d] addActiveWalk interList (%d)\n", thisIndex, activeWalks.length());
#endif
#else
  addActiveWalk(remoteGravityAwi,sTopDown,sGravity,sRemote,
		sRemoteGravityState);
#ifdef CHECK_WALK_COMPLETIONS
  CkPrintf("[%d] addActiveWalk remoteGravity (%d)\n", thisIndex, activeWalks.length());
#endif
#endif


  // variable currentBucket masquerades as current chunk
  initiatePrefetch(sPrefetchState->currentBucket);

#if CHANGA_REFACTOR_DEBUG > 0
  CkPrintf("[%d]sending message to commence local gravity calculation\n", thisIndex);
#endif

  if (bEwald) thisProxy[thisIndex].EwaldInit();
#if defined CUDA
  // ask datamanager to serialize local trees
  // prefetch can occur concurrently with this, 
  // though calculateGravityLocal can only come
  // afterwards.
  dm->serializeLocalTree();
#else
  thisProxy[thisIndex].commenceCalculateGravityLocal();
#endif


#ifdef CHANGA_PRINT_MEMUSAGE
      int piecesPerPe = numTreePieces/CmiNumPes();
      if(thisIndex % piecesPerPe == 0)
	CkPrintf("[%d]: CmiMaxMemoryUsage: %f M\n", CmiMyPe(),
		 (float)CmiMaxMemoryUsage()/(1 << 20));
#endif
}

// Starts the prefetching of the specified chunk; once the 
// given chunk has been completely prefetched, the prefetch
// compute invokes startRemoteChunk() 
void TreePiece::initiatePrefetch(int chunk){
#ifdef DISABLE_NODE_TREE
  GenericTreeNode *child = keyToNode(prefetchRoots[chunk]);
#else
  GenericTreeNode *child = dm->chunkRootToNode(prefetchRoots[chunk]);
#endif
  CmiAssert(child != NULL);

  int first, last;
  for(int x = -nReplicas; x <= nReplicas; x++) {
    for(int y = -nReplicas; y <= nReplicas; y++) {
      for(int z = -nReplicas; z <= nReplicas; z++) {
        if (child == NULL) {
          nodeOwnership(prefetchRoots[chunk], first, last);
          child = requestNode(dm->responsibleIndex[(first+last)>>1],
              prefetchRoots[chunk], chunk,
              encodeOffset(0, x, y, z),
              prefetchAwi, (void *)0, true);
        }
        if (child != NULL) {
#if CHANGA_REFACTOR_DEBUG > 1
          CkPrintf("[%d] starting prefetch walk with current Prefetch=%d, numPrefetchReq=%d (%d,%d,%d)\n", thisIndex, chunk, numPrefetchReq, x,y,z);
#endif
          sTopDown->walk(child, sPrefetchState, chunk, encodeOffset(0,x,y,z), prefetchAwi);
        }
      }
    }
  }

}

void TreePiece::commenceCalculateGravityLocal(){
#if INTERLIST_VER > 0 
  // must set placedRoots to false before starting local comp.
  DoubleWalkState *lstate = (DoubleWalkState *)sLocalGravityState;
  lstate->placedRoots[0] = false;
#endif
  calculateGravityLocal();
}

void TreePiece::startRemoteChunk() {
    //CkPrintf("[%d] in startRemoteChunk\n", thisIndex);
#if CHANGA_REFACTOR_DEBUG > 0
  CkPrintf("[%d] sending message to commence remote gravity\n", thisIndex);
#endif

  traceUserEvent(prefetchDoneUE);

#ifdef CUDA
  // dm counts until all treepieces have acknowledged prefetch completion
  // it then flattens the tree on the processor, sends it to the device
  // and sends messages to each of the registered treepieces to continueStartRemoteChunk()
  //CkPrintf("[%d] startRemoteChunk done current Prefetch: %d\n", thisIndex, current Prefetch);
  dm->donePrefetch(sPrefetchState->currentBucket);
#else
  continueStartRemoteChunk(sPrefetchState->currentBucket);
#endif
}

/// @brief Main work of StartRemoteChunk()
/// Schedule a TreePiece::calculateGravityRemote() then start
/// prefetching for the next chunk.
void TreePiece::continueStartRemoteChunk(int chunk){
  // FIXME - can value of chunk be different from current Prefetch?
  //CkPrintf("[%d] continueStartRemoteChunk chunk: %d, current Prefetch: %d\n", thisIndex, chunk, current Prefetch);
  ComputeChunkMsg *msg = new (8*sizeof(int)) ComputeChunkMsg(sPrefetchState->currentBucket);
  *(int*)CkPriorityPtr(msg) = numTreePieces * sPrefetchState->currentBucket + thisIndex + 1;
  CkSetQueueing(msg, CK_QUEUEING_IFIFO);

#if INTERLIST_VER > 0
  DoubleWalkState *rstate = (DoubleWalkState *)sRemoteGravityState;
  rstate->placedRoots[msg->chunkNum] = false;
#endif
  thisProxy[thisIndex].calculateGravityRemote(msg);

  // start prefetching next chunk
  if (++sPrefetchState->currentBucket < numChunks) {
    // Nothing needs to be changed for this chunk -
    // the prefetchReqs and their number remains the same
    // We only need to reassociate the tree walk with the
    // prefetch compute object and the prefetch object wiht
    // the prefetch opt object
    sTopDown->reassoc(sPrefetch);
    // prefetch walk isn't associated with any particular bucket
    // but the entire treepiece
    // this method invocation does nothing. indeed, nothing
    // needs to be done because sPrefetch is always associated with
    // sPref
    sPrefetch->reassoc((void *)0,activeRung,sPref);

    // instead of prefetchWaiting, we count through state->counters[0]
    //prefetchWaiting = (2*nReplicas + 1)*(2*nReplicas + 1)*(2*nReplicas + 1);
    sPrefetchState->counterArrays[0][0] = (2*nReplicas + 1)*(2*nReplicas + 1)*(2*nReplicas + 1);

    // variable currentBucket masquerades as current chunk
    initiatePrefetch(sPrefetchState->currentBucket);
  }
}

/*
void TreePiece::startlb(CkCallback &cb){
  callback = cb;
  if(verbosity > 1)
    CkPrintf("[%d] TreePiece %d calling AtSync()\n",CkMyPe(),thisIndex);
  localCache->revokePresence(thisIndex);
  AtSync();
}
*/
  // jetley - contribute your centroid. AtSync is now called by the load balancer (broadcast) when it has
  // all centroids.
void TreePiece::startlb(CkCallback &cb, int activeRung){

  if(verbosity > 1)
     CkPrintf("[%d] load set to: %g, actual: %g\n", thisIndex, treePieceLoad, getObjTime());  

  callback = cb;
  if(verbosity > 1)
    CkPrintf("[%d] TreePiece %d calling AtSync()\n",CkMyPe(),thisIndex);
  
  unsigned int numActiveParticles, i;

  if(activeRung == 0){
    numActiveParticles = myNumParticles;
  }
  else if(activeRung == PHASE_FEEDBACK) {
      numActiveParticles = myNumSPH + myNumStar;
      }
  else{
    for(numActiveParticles = 0, i = 1; i <= myNumParticles; i++)
      if(myParticles[i].rung >= activeRung)
        numActiveParticles++;
  }
  LDObjHandle myHandle = myRec->getLdHandle();
  TaggedVector3D tv(savedCentroid, myHandle, numActiveParticles, myNumParticles, activeRung, prevLARung);
  tv.tp = thisIndex;
  tv.tag = thisIndex;

  treePieceActivePartsTmp = numActiveParticles;
  if (havePhaseData(activeRung)) {
    treePieceLoadExp = savedPhaseLoad[activeRung];
  } else if (havePhaseData(0)) {
    float ratio = 1.0;
    if(myNumParticles != 0){
      ratio = numActiveParticles/(float)myNumParticles;
    }

    treePieceLoadExp  = ratio * savedPhaseLoad[0];
  } else {
    treePieceLoadExp =  treePieceLoad;
  }
  setObjTime(treePieceLoadExp);
  treePieceLoad = 0;

  /*
  CkPrintf("[%d] centroid %f %f %f\n", 
                      thisIndex,
                      tv.vec.x,
                      tv.vec.y,
                      tv.vec.z
                      );
  */

  if(foundLB == Multistep){
    CkCallback cbk(CkIndex_MultistepLB::receiveCentroids(NULL), 0, proxy);
    contribute(sizeof(TaggedVector3D), (char *)&tv, CkReduction::concat, cbk);
  }
  else if(foundLB == Orb3d){
    CkCallback cbk(CkIndex_Orb3dLB::receiveCentroids(NULL), 0, proxy);
    contribute(sizeof(TaggedVector3D), (char *)&tv, CkReduction::concat, cbk);
  }
  else if(foundLB == Multistep_notopo){
    CkCallback lbcb(CkIndex_MultistepLB_notopo::receiveCentroids(NULL), 0, proxy);
    contribute(sizeof(TaggedVector3D), (char *)&tv, CkReduction::concat, lbcb);
  }
  else if(foundLB == Orb3d_notopo){
    CkCallback cbk(CkIndex_Orb3dLB_notopo::receiveCentroids(NULL), 0, proxy);
    contribute(sizeof(TaggedVector3D), (char *)&tv, CkReduction::concat, cbk);
  }
  else if(foundLB == MultistepOrb){
    CkCallback lbcb(CkIndex_MultistepOrbLB::receiveCentroids(NULL), 0, proxy);
    contribute(sizeof(TaggedVector3D), (char *)&tv, CkReduction::concat, lbcb);
  }
  else if(activeRung == 0) {
    doAtSync();
  }
  else {
    contribute(cb);  // Skip the load balancer
  }
  prevLARung = activeRung;
}

void TreePiece::doAtSync(){
  if(verbosity > 1)
      CkPrintf("[%d] TreePiece %d calling AtSync() at %g\n",CkMyPe(),thisIndex, CkWallTimer());
  AtSync();
}

void TreePiece::ResumeFromSync(){
  setObjTime(0.0); 
  if(verbosity > 1)
    CkPrintf("[%d] TreePiece %d in ResumefromSync\n",CkMyPe(),thisIndex);
  contribute(callback);
}

/*
inline GenericTreeNode *TreePiece::keyToNode(const Tree::NodeKey k) {
  NodeLookupType::iterator iter = nodeLookupTable.find(k);
  if (iter != nodeLookupTable.end()) return iter->second;
  else return NULL;
}
*/

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


// Walk a node evalutating its force on a bucket
//
// @param node Node to be walked
// @param reqID request ID which encodes bucket to be worked on and
// a replica offset if we are using periodic boundary conditions
void TreePiece::walkBucketTree(GenericTreeNode* node, int reqID) {
    Vector3D<double> offset = decodeOffset(reqID);
    int reqIDlist = decodeReqID(reqID);
#if COSMO_STATS > 0
  myNumMACChecks++;
#endif
#if COSMO_PRINT > 0
  if (node->getKey() < numChunks) CkPrintf("[%d] walk: checking %016llx\n",thisIndex,node->getKey());
#endif
#if COSMO_STATS > 0
  numOpenCriterionCalls++;
#endif
  GenericTreeNode* reqnode = bucketList[reqIDlist];
  if(!openCriterionBucket(node, reqnode, offset, localIndex)) {
      // Node is well separated; evaluate the force
#if COSMO_STATS > 1
    MultipoleMoments m = node->moments;
    for(int i = reqnode->firstParticle; i <= reqnode->lastParticle; ++i)
      myParticles[i].intcellmass += m.totalMass;
#endif
#if COSMO_PRINT > 1
    CkPrintf("[%d] walk bucket %s -> node %s\n",thisIndex,keyBits(bucketList[reqIDlist]->getKey(),KeyBits).c_str(),keyBits(node->getKey(),KeyBits).c_str());
#endif
#if COSMO_DEBUG > 1
    bucketcheckList[reqIDlist].insert(node->getKey());
    combineKeys(node->getKey(),reqIDlist);
#endif
#if COSMO_PRINT > 0
    if (node->getKey() < numChunks) CkPrintf("[%d] walk: computing %016llx\n",thisIndex,node->getKey());
#endif
#ifdef COSMO_EVENTS
    double startTimer = CmiWallTimer();
#endif
#ifdef HPM_COUNTER
    hpmStart(1,"node force");
#endif
    int computed = nodeBucketForce(node, reqnode, myParticles, offset, activeRung);
#ifdef HPM_COUNTER
    hpmStop(1);
#endif
#ifdef COSMO_EVENTS
    traceUserBracketEvent(nodeForceUE, startTimer, CmiWallTimer());
#endif
    nodeInterLocal += computed;
  } else if(node->getType() == Bucket) {
    int computed = 0;
    for(int i = node->firstParticle; i <= node->lastParticle; ++i) {
#if COSMO_STATS > 1
      for(int j = reqnode->firstParticle; j <= reqnode->lastParticle; ++j) {
        //myParticles[j].intpartmass += myParticles[i].mass;
        myParticles[j].intpartmass += node->particlePointer[i-node->firstParticle].mass;
      }
#endif
#if COSMO_PRINT > 1
      //CkPrintf("[%d] walk bucket %s -> part %016llx\n",thisIndex,keyBits(reqnode->getKey(),63).c_str(),myParticles[i].key);
      CkPrintf("[%d] walk bucket %s -> part %016llx\n",thisIndex,keyBits(reqnode->getKey(),KeyBits).c_str(),node->particlePointer[i-node->firstParticle].key);
#endif
#ifdef COSMO_EVENTS
      double startTimer = CmiWallTimer();
#endif
#ifdef HPM_COUNTER
    hpmStart(2,"particle force");
#endif
      computed = partBucketForce(&node->particlePointer[i-node->firstParticle], reqnode,
                                 myParticles, offset, activeRung);
#ifdef HPM_COUNTER
    hpmStop(2);
#endif
#ifdef COSMO_EVENTS
      traceUserBracketEvent(partForceUE, startTimer, CmiWallTimer());
#endif
    }
    particleInterLocal += node->particleCount * computed;
#if COSMO_DEBUG > 1
    bucketcheckList[reqIDlist].insert(node->getKey());
    combineKeys(node->getKey(),reqIDlist);
#endif
  } else if (node->getType() == NonLocal || node->getType() == NonLocalBucket) {
    /* DISABLED: this part of the walk is triggered directly by the CacheManager and prefetching

    // Use cachedWalkBucketTree() as callback
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

GenericTreeNode* TreePiece::requestNode(int remoteIndex, Tree::NodeKey key, int chunk, int reqID, int awi, void *source, bool isPrefetch) {

  CkAssert(remoteIndex < (int) numTreePieces);
  CkAssert(chunk < numChunks);

  if(_cache){
    //CkAssert(localCache != NULL);
#if COSMO_PRINT > 1

    CkPrintf("[%d] b=%d requesting node %s to %d for %s (additional=%d)\n",thisIndex,reqID,keyBits(key,KeyBits).c_str(),remoteIndex,keyBits(bucketList[reqID]->getKey(),KeyBits).c_str(),
            sRemoteGravityState->counterArrays[0][decodeReqID(reqID)]
            + sLocalGravityState->counterArrays[0][decodeReqID(reqID)]);
#endif
    //GenericTreeNode
    CProxyElement_ArrayElement thisElement(thisProxy[thisIndex]);
    CkCacheUserData userData;
    userData.d0 = (((CmiUInt8) awi)<<32)+reqID;
    userData.d1 = (CmiUInt8) source;

    CkCacheRequestorData<KeyType> request(thisElement, &EntryTypeGravityNode::callback, userData);
    CkArrayIndexMax remIdx = CkArrayIndex1D(remoteIndex);
    GenericTreeNode *res = (GenericTreeNode *) cacheNode.ckLocalBranch()->requestData(key,remIdx,chunk,&gravityNodeEntry,request);

#ifdef CHANGA_REFACTOR_INTERLIST_PRINT_BUCKET_START_FIN
    if(source && !res){
      int start, end;
      GenericTreeNode *testSource = (GenericTreeNode *)source;
      getBucketsBeneathBounds(testSource, start, end);
      CkPrintf("[%d] tgt=%d requesting node %ld source = %ld(%d - %d) buckRem = %d + %d, remChunk = %d\n",thisIndex,decodeReqID(reqID), key, testSource->getKey(), start, end, sRemoteGravityState->counterArrays[0][decodeReqID(reqID)], sLocalGravityState->counterArrays[0][decodeReqID(reqID)], sRemoteGravityState->counterArrays[1][chunk]);

      // check whether source contains target
      GenericTreeNode *target = bucketList[decodeReqID(reqID)];
      if(!testSource->contains(target->getKey())){
        CkAbort("source does not contain target\n");
      }
    }
#endif

    return res;
  }
  else{
    CkAbort("Non cached version not anymore supported, feel free to fix it!");
  }
  return NULL;
}

ExternalGravityParticle *TreePiece::requestParticles(Tree::NodeKey key,int chunk,int remoteIndex,int begin,int end,int reqID, int awi, void *source, bool isPrefetch) {
  if (_cache) {
    //CkAssert(localCache != NULL);
    //ExternalGravityParticle *p = localCache->requestParticles(thisIndex,chunk,key,remoteIndex,begin,end,reqID,awi,isPrefetch);
    CProxyElement_ArrayElement thisElement(thisProxy[thisIndex]);
    CkCacheUserData userData;
    userData.d0 = (((CmiUInt8) awi)<<32)+reqID;
    userData.d1 = (CmiUInt8) source;

    CkCacheRequestorData<KeyType> request(thisElement, &EntryTypeGravityParticle::callback, userData);
    CkArrayIndexMax remIdx = CkArrayIndex1D(remoteIndex);
    CkAssert(key < (((NodeKey) 1) << (NodeKeyBits-1)));
    //
    // Key is shifted to distiguish between nodes and particles
    //
    KeyType ckey = key<<1;
    CacheParticle *p = (CacheParticle *) cacheGravPart.ckLocalBranch()->requestData(ckey,remIdx,chunk,&gravityParticleEntry,request);
    if (p == NULL) {
#ifdef CHANGA_REFACTOR_INTERLIST_PRINT_BUCKET_START_FIN
      if(source){
        int start, end;
        GenericTreeNode *testSource = (GenericTreeNode *)source;
        getBucketsBeneathBounds(testSource, start, end);
        CkPrintf("[%d] tgt=%d requesting particles %ld source = %ld(%d - %d) buckRem = %d + %d, remChunk = %d\n",thisIndex,decodeReqID(reqID), key, testSource->getKey(), start, end, sRemoteGravityState->counterArrays[0][decodeReqID(reqID)], sLocalGravityState->counterArrays[0][decodeReqID(reqID)], sRemoteGravityState->counterArrays[1][chunk]);

        // check whether source contains target
        GenericTreeNode *target = bucketList[decodeReqID(reqID)];
        if(!testSource->contains(target->getKey())){
          CkAbort("source does not contain target\n");
        }
      }
#endif

#if COSMO_PRINT > 1

      CkPrintf("[%d] b=%d requestParticles: additional=%d\n",thisIndex,
	       decodeReqID(reqID),
               sRemoteGravityState->counterArrays[0][decodeReqID(reqID)]
               + sLocalGravityState->counterArrays[0][decodeReqID(reqID)]);
#endif
      return NULL;
    }
    return p->part;
  } else {
    CkAbort("Non cached version not anymore supported, feel free to fix it!");
  }
  return NULL;
}

GravityParticle *
TreePiece::requestSmoothParticles(Tree::NodeKey key,int chunk,int remoteIndex,
				  int begin,int end,int reqID, int awi, void *source,
				  bool isPrefetch) {
  if (_cache) {
    CProxyElement_ArrayElement thisElement(thisProxy[thisIndex]);
    CkCacheUserData userData;
    userData.d0 = (((CmiUInt8) awi)<<32)+reqID;
    userData.d1 = (CmiUInt8) source;

    CkCacheRequestorData<KeyType> request(thisElement, &EntryTypeSmoothParticle::callback, userData);
    CkArrayIndexMax remIdx = CkArrayIndex1D(remoteIndex);
    KeyType ckey = key<<1;
    CacheSmoothParticle *p = (CacheSmoothParticle *) cacheSmoothPart.ckLocalBranch()->requestData(ckey,remIdx,chunk,&smoothParticleEntry,request);
    if (p == NULL) {
      return NULL;
    }
    return p->partCached;
  } else {
    CkAbort("Non cached version not anymore supported, feel free to fix it!");
  }
  return NULL;
}

#if COSMO_DEBUG > 1 || defined CHANGA_REFACTOR_WALKCHECK || defined CHANGA_REFACTOR_WALKCHECK_INTERLIST

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

  std::multiset<Tree::NodeKey>::iterator iter = (bucketcheckList[bucket]).find(sibKey);

  if(iter==bucketcheckList[bucket].end())
    return;
  else{//Sibling key has been found in the Binary tree
    bucketcheckList[bucket].erase(bucketcheckList[bucket].find(key));
    bucketcheckList[bucket].erase(iter);
    if(bucket == TEST_BUCKET && thisIndex == TEST_TP){
      CkPrintf("[%d] combine(%ld, %ld)\n", thisIndex, key, sibKey);
      CkPrintf("[%d] add %ld\n", thisIndex, key >> 1, sibKey);
    }
    key >>= 1;
    bucketcheckList[bucket].insert(key);
    combineKeys(key,bucket);
  }
}

void TreePiece::checkWalkCorrectness(){

  Tree::NodeKey endKey = Key(1);
  int count = (2*nReplicas+1) * (2*nReplicas+1) * (2*nReplicas+1);
  CkPrintf("[%d] checking walk correctness...\n",thisIndex);
  bool someWrong = false;

  for(int i=0;i<numBuckets;i++){
    int wrong = 0;
    if(bucketList[i]->rungs < activeRung) continue;
    if(bucketcheckList[i].size()!=count) wrong = 1;
    for (std::multiset<Tree::NodeKey>::iterator iter = bucketcheckList[i].begin(); iter != bucketcheckList[i].end(); iter++) {
      if (*iter != endKey) wrong = 1;
    }
    if (wrong) {
      someWrong = true;
      CkPrintf("Error: [%d] Not all nodes were traversed by bucket %d\n",thisIndex,i);
      for (std::multiset<Tree::NodeKey>::iterator iter=bucketcheckList[i].begin(); iter != bucketcheckList[i].end(); iter++) {
	CkPrintf("       [%d] key %ld\n",thisIndex,*iter);
      }
    }
    else { bucketcheckList[i].clear(); }
  }
  if(someWrong) CkExit();
}
#endif

/********************************************************************/

void TreePiece::outputStatistics(const CkCallback& cb) {

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

  if(thisIndex != (int) numTreePieces - 1)
    pieces[thisIndex + 1].outputStatistics(cb);
  if(thisIndex == (int) numTreePieces - 1) cb.send();
}

/// @TODO Fix pup routine to handle correctly the tree
void TreePiece::pup(PUP::er& p) {
  CBase_TreePiece::pup(p);

  p | warmupFinished;
  p | treePieceLoad; 
  p | treePieceLoadExp;
  p | treePieceActivePartsTmp;
  p | savedPhaseLoad;
  p | savedPhaseParticle;

  // jetley
  p | proxy;
  p | foundLB;
  p | proxyValid;
  p | proxySet;
  p | savedCentroid;
  p | prevLARung;

  p | callback;
  p | nTotalParticles;
  p | myNumParticles;
  p | nTotalSPH;
  p | myNumSPH;
  p | nTotalDark;
  p | nTotalStar;
  p | myNumStar;
  if(p.isUnpacking()) {
      nStore = (int)((myNumParticles + 2)*(1.0 + dExtraStore));
      myParticles = new GravityParticle[nStore];
      nStoreSPH = (int)(myNumSPH*(1.0 + dExtraStore));
      if(nStoreSPH > 0) mySPHParticles = new extraSPHData[nStoreSPH];
      allocateStars();
  }
  for(unsigned int i=1;i<=myNumParticles;i++){
    p | myParticles[i];
    if(myParticles[i].isGas()) {
	int iSPH;
	if(!p.isUnpacking()) {
	    iSPH = (extraSPHData *)myParticles[i].extraData - mySPHParticles;
	    CkAssert(iSPH < myNumSPH);
	    p | iSPH;
	    }
	else {
	    p | iSPH;
	    myParticles[i].extraData = mySPHParticles + iSPH;
	    CkAssert(iSPH < myNumSPH);
	    }
	}
    if(myParticles[i].isStar()) {
	int iStar;
	if(!p.isUnpacking()) {
	    iStar = (extraStarData *)myParticles[i].extraData - myStarParticles;
	    CkAssert(iStar < myNumStar);
	    p | iStar;
	    }
	else {
	    p | iStar;
	    myParticles[i].extraData = myStarParticles + iStar;
	    CkAssert(iStar < myNumStar);
	    }
	}
  }
  for(unsigned int i=0;i<myNumSPH;i++){
    p | mySPHParticles[i];
  }
  for(unsigned int i=0;i<myNumStar;i++){
    p | myStarParticles[i];
  }
  p | pieces;
  p | basefilename;
  p | boundingBox;
  p | iterationNo;

  p | nSetupWriteStage;

  // Periodic variables
  p | nReplicas;
  p | fPeriod;
  p | bEwald;
  p | fEwCut;
  p | dEwhCut;
  p | bPeriodic;
  p | nMaxEwhLoop;
  if (p.isUnpacking() && bEwald) {
    ewt = new EWT[nMaxEwhLoop];
  }

  p | numBuckets;
#if INTERLIST_VER > 0
  p | prevBucket;
  p | prevRemoteBucket;
#endif
  p | ewaldCurrentBucket;
#if COSMO_STATS > 0
  //p | myNumParticleInteractions;
  //p | myNumCellInteractions;
  p | myNumMACChecks;
  p | nodesOpenedLocal;
  p | nodesOpenedRemote;
  p | numOpenCriterionCalls;
  p | piecemass;
#endif
  if (p.isUnpacking()) {
    particleInterRemote = NULL;
    nodeInterRemote = NULL;

    switch(domainDecomposition) {
    case SFC_dec:
    case SFC_peano_dec:
    case SFC_peano_dec_3D:
    case SFC_peano_dec_2D:
      numPrefetchReq = 2;
      break;
    case Oct_dec:
    case ORB_dec:
    case ORB_space_dec:
      numPrefetchReq = 1;
      break;
    default:
      CmiAbort("Pupper has wrong domain decomposition type!\n");
    }
  }

  p | myPlace;

  p | bGasCooling;
  if(p.isUnpacking()){
    dm = NULL;
#ifndef COOLING_NONE
    if(!_inrestart) {           // not restarting from checkpoint
        dm = (DataManager*)CkLocalNodeBranch(dataManagerID);
        if(bGasCooling)
            CoolData = CoolDerivsInit(dm->Cool);
        }
#endif
  }

  if (verbosity > 1) {
    if (p.isSizing()) {
      ckout << "TreePiece " << thisIndex << ": Getting PUP'd!";
      ckout << " size: " << ((PUP::sizer*)&p)->size();
      ckout << endl;
      }
  }
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

/** Check that all the particles in the tree are really in their boxes.
    Because the keys are made of only the first 21 out of 23 bits of the
    floating point representation, there can be particles that are outside
    their box by tiny amounts.  Whether this is bad is not yet known. */
void TreePiece::checkTree(GenericTreeNode* node) {
  if(node->getType() == Empty) return;
  if(node->getType() == Bucket) {
    for(unsigned int iter = node->firstParticle; iter <= node->lastParticle; ++iter) {
      if(!node->boundingBox.contains(myParticles[iter].position)) {
	ckerr << "Not in the box: Box: " << node->boundingBox << " Position: " << myParticles[iter].position << "\nNode key: " << keyBits(node->getKey(), KeyBits).c_str() << "\nParticle key: " << keyBits(myParticles[iter].key, KeyBits).c_str() << endl;
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
  oss << keyBits(node->getKey(), NodeKeyBits) << "\\n";
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

/// Print a text version of a tree
void TreePiece::printTree(GenericTreeNode* node, ostream& os) {
  if(node == 0)
    return;

  string nodeID = keyBits(node->getKey(), NodeKeyBits);
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
#ifndef HEXADECAPOLE
  if (node->getType() == Bucket || node->getType() == Internal || node->getType() == Boundary || node->getType() == NonLocal || node->getType() == NonLocalBucket)
    os << " V "<<node->moments.radius<<" "<<node->moments.soft<<" "<<node->moments.cm.x<<" "<<node->moments.cm.y<<" "<<node->moments.cm.z<<" "<<node->moments.xx<<" "<<node->moments.xy<<" "<<node->moments.xz<<" "<<node->moments.yy<<" "<<node->moments.yz<<" "<<node->moments.zz<<" "<<node->boundingBox;
#endif
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

/// Print a graphviz version of A tree
void TreePiece::printTreeViz(GenericTreeNode* node, ostream& os) {
  if(node == 0)
    return;

  string nodeID = keyBits(node->getKey(), NodeKeyBits);
  os << "\tnode [color=\"" << getColor(node) << "\"]\n";
  //os << "\t\"" << nodeID << "\" [label=\"" << makeLabel(node) << "\\nCM: " << (node->moments.cm) << "\\nM: " << node->moments.totalMass << "\\nN_p: " << (node->endParticle - node->beginParticle) << "\\nOwners: " << node->numOwners << "\"]\n";
  //os << "\t\"" << nodeID << "\" [label=\"" << makeLabel(node) << "\\nLocal N: " << (node->endParticle - node->beginParticle) << "\\nOwners: " << node->numOwners << "\"]\n";
  os << "\t\"" << nodeID << "\" [label=\"" << keyBits(node->getKey(), NodeKeyBits) << "\\n";
  int first, last;
  switch(node->getType()) {
  case Bucket:
    os << "Bucket\\nSize: " << (node->lastParticle - node->firstParticle + 1);
    break;
  case Internal:
    os << "Internal\\nSize: " << (node->lastParticle - node->firstParticle + 1);
    break;
  case NonLocal:
    nodeOwnership(node->getKey(), first, last);
    os << "NonLocal: Chare " << node->remoteIndex << "\\nOwners: " << (last-first+1) << "\\nRemote size: " << (node->lastParticle - node->firstParticle + 1);
    //os << "NonLocal: Chare=" << node->remoteIndex; //<< ", Owners=" << first << "-" << last;
    break;
  case NonLocalBucket:
    //os << "NonLocal: Chare=" << node->remoteIndex << "\\nRemote N under: " << (node->lastParticle - node->firstParticle + 1) << "\\nOwners: " << node->numOwners;
    nodeOwnership(node->getKey(), first, last);
    //CkAssert(first == last);
    os << "NonLocalBucket: Chare " << node->remoteIndex << "\\nOwner: " << first << "\\nSize: " << node->particleCount; //<< "(" << node->firstParticle << "-" << node->lastParticle << ")";
    break;
  case Boundary:
    //nodeOwnership(node->getKey(), first, last);
    os << "Boundary\\nTotalsize: " << node->particleCount << "\\nLocalsize: " << (node->lastParticle - node->firstParticle);
    break;
  case Empty:
    os << "Empty "<<node->remoteIndex;
    break;
  }
  //if (node->getType() == Bucket || node->getType() == Internal || node->getType() == Boundary || node->getType() == NonLocal || node->getType() == NonLocalBucket)
  //  os << " V "<<node->moments.radius<<" "<<node->moments.soft<<" "<<node->moments.cm.x<<" "<<node->moments.cm.y<<" "<<node->moments.cm.z<<" "<<node->moments.xx<<" "<<node->moments.xy<<" "<<node->moments.xz<<" "<<node->moments.yy<<" "<<node->moments.yz<<" "<<node->moments.zz;

  os << "\"]\n";

  if(node->parent)
    os << "\t\"" << keyBits(node->parent->getKey(), NodeKeyBits) << "\" -> \"" << nodeID << "\";\n";

  if(node->getType() == NonLocal || node->getType() == NonLocalBucket || node->getType() == Bucket || node->getType() == Empty)
    return;

  GenericTreeNode* childIterator;
  for(unsigned int i = 0; i < node->numChildren(); ++i) {
    childIterator = node->getChildren(i);
    if(childIterator)
      printTreeViz(childIterator, os);
    else {
      os << "\tnode [color=\"green\"]\n";
      os << "\t\"" << nodeID << i << "\" [label=\"None\"]\n";
      os << "\t\"" << nodeID << "\" -> \"" << nodeID << i << "\";\n";
    }
  }
}

/// Write a file containing a graphviz dot graph of my tree
void TreePiece::report() {
  ostringstream outfilename;
  outfilename << "tree." << thisIndex << "." << iterationNo << ".dot";
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
  printTreeViz(root, os);
  os << "}" << endl;

  os.close();

  //checkTree(root);

  //contribute(cb);
}

/// Print a text version of a tree
void printGenericTree(GenericTreeNode* node, ostream& os) {
  if(node == 0)
    return;

  string nodeID = keyBits(node->getKey(), NodeKeyBits);
  os << nodeID << " ";
  switch(node->getType()) {
  case Bucket:
    os << "Bucket: Size=" << (node->lastParticle - node->firstParticle + 1) << "(" << node->firstParticle << "-" << node->lastParticle << ")";
    break;
  case Internal:
    os << "Internal: Size=" << (node->lastParticle - node->firstParticle + 1) << "(" << node->firstParticle << "-" << node->lastParticle << ")";
    break;
  case NonLocal:
    os << "NonLocal: Chare=" << node->remoteIndex; //<< ", Owners=" << first << "-" << last;
    break;
  case NonLocalBucket:
    os << "NonLocalBucket: Chare=" << node->remoteIndex << ", Size=" << node->particleCount; //<< "(" << node->firstParticle << "-" << node->lastParticle << ")";
    break;
  case Boundary:
    os << "Boundary: Totalsize=" << node->particleCount << ", Localsize=" << (node->lastParticle - node->firstParticle) << "(" << node->firstParticle << "-" << node->lastParticle;
    break;
  case Empty:
    os << "Empty "<<node->remoteIndex;
    break;
  }
#ifndef HEXADECAPOLE
  if (node->getType() == Bucket || node->getType() == Internal || node->getType() == Boundary || node->getType() == NonLocal || node->getType() == NonLocalBucket)
    os << " V "<<node->moments.radius<<" "<<node->moments.soft<<" "<<node->moments.cm.x<<" "<<node->moments.cm.y<<" "<<node->moments.cm.z<<" "<<node->moments.xx<<" "<<node->moments.xy<<" "<<node->moments.xz<<" "<<node->moments.yy<<" "<<node->moments.yz<<" "<<node->moments.zz;
#endif

  os << "\n";

  //if(node->parent)
  //  os << "\t\"" << keyBits(node->parent->getKey(), 63) << "\" -> \"" << nodeID << "\";\n";

  if(node->getType() == NonLocal || node->getType() == NonLocalBucket || node->getType() == Bucket || node->getType() == Empty)
    return;

  GenericTreeNode* childIterator;
  for(unsigned int i = 0; i < node->numChildren(); ++i) {
    childIterator = node->getChildren(i);
    if(childIterator)
      printGenericTree(childIterator, os);
    else {
      os << "\tnode [color=\"green\"]\n";
      os << "\t\"" << nodeID << i << "\" [label=\"None\"]\n";
      os << "\t\"" << nodeID << "\" -> \"" << nodeID << i << "\";\n";
    }
  }
}

CkReduction::reducerType TreePieceStatistics::sum;

/*
 * Collect treewalking statistics across all TreePieces
 */

void TreePiece::collectStatistics(CkCallback& cb) {
  LBTurnInstrumentOff();
#if COSMO_DEBUG > 1 || defined CHANGA_REFACTOR_WALKCHECK || defined CHANGA_REFACTOR_WALKCHECK_INTERLIST

checkWalkCorrectness();
#endif

#if COSMO_STATS > 0
  u_int64_t nodeInterRemoteTotal = 0;
  u_int64_t particleInterRemoteTotal = 0;
  for (int i=0; i<numChunks; ++i) {
    nodeInterRemoteTotal += nodeInterRemote[i];
    particleInterRemoteTotal += particleInterRemote[i];
  }
  TreePieceStatistics tps(nodesOpenedLocal, nodesOpenedRemote, numOpenCriterionCalls,
      nodeInterLocal, nodeInterRemoteTotal, particleInterLocal, particleInterRemoteTotal, nActive);
  contribute(sizeof(TreePieceStatistics), &tps, TreePieceStatistics::sum, cb);
#else
  CkAbort("Invalid call, only valid if COSMO_STATS is defined");
#endif
}

GenericTreeNode *TreePiece::nodeMissed(int reqID, int remoteIndex, Tree::NodeKey &key, int chunk, bool isPrefetch, int awi, void *source){
  GenericTreeNode *gtn = requestNode(remoteIndex, key, chunk, reqID, awi, source, isPrefetch);
  return gtn;
}

ExternalGravityParticle *TreePiece::particlesMissed(Tree::NodeKey &key, int chunk, int remoteIndex, int firstParticle, int lastParticle, int reqID, bool isPrefetch, int awi, void *source){
  return requestParticles(key, chunk, remoteIndex,firstParticle,lastParticle,reqID, awi, source, isPrefetch);
}

// This is invoked when a remote node is received from the CacheManager
// It sets up a tree walk starting at node and initiates it
void TreePiece::receiveNodeCallback(GenericTreeNode *node, int chunk, int reqID, int awi, void *source){
  int targetBucket = decodeReqID(reqID);
  Vector3D<double> offset = decodeOffset(reqID);

  TreeWalk *tw;
  Compute *compute;
  State *state;

#ifdef CHANGA_REFACTOR_INTERLIST_PRINT_BUCKET_START_FIN
  if(source){
    int start, end;
    GenericTreeNode *testSource = (GenericTreeNode *)source;
    getBucketsBeneathBounds(testSource, start, end);
    CkPrintf("[%d] tgt=%d recv node %ld source = %ld(%d - %d) buckRem = %d + %d, remChunk = %d\n",thisIndex,decodeReqID(reqID), node->getKey(), testSource->getKey(), start, end, sRemoteGravityState->counterArrays[0][decodeReqID(reqID)], sLocalGravityState->counterArrays[0][decodeReqID(reqID)], sRemoteGravityState->counterArrays[1][chunk]);

      // check whether source contains target
      GenericTreeNode *target = bucketList[decodeReqID(reqID)];
      if(!testSource->contains(target->getKey())){
        CkAbort("source does not contain target\n");
      }
  }
#endif
  //nodeMisses--;
  //CkPrintf("- nodeMisses: %d\n", nodeMisses);
  // retrieve the activewalk record
  CkAssert(awi < activeWalks.size());

  ActiveWalk &a = activeWalks[awi];
  tw = a.tw;
  compute = a.c;
  state = a.s;

  // reassociate objects with each other
  tw->reassoc(compute);
  compute->reassoc(source, activeRung, a.o);

  // resume walk
  //CkPrintf("[%d] RECVD NODE (%ld), resuming walk chunk %d\n", thisIndex, node->getKey(), chunk);
  tw->resumeWalk(node, state, chunk, reqID, awi);

  // we need source to update the counters in all buckets
  // underneath the source. note that in the interlist walk,
  // the computeEntity of the compute will likely  have changed as the walk continued.
  // however, the resumeWalk function takes care to set it back to 'source'
  // after it is done walking.
  /*
  if(awi == interListAwi){
    GenericTreeNode *src = (GenericTreeNode *)source;
    CkPrintf("[%d]: recvd node (%ld) for (%ld)\n", thisIndex, node->getKey(), src->getKey());
  }
  */
  compute->nodeRecvdEvent(this,chunk,state,targetBucket);
}

void TreePiece::receiveParticlesCallback(ExternalGravityParticle *egp, int num, int chunk, int reqID, Tree::NodeKey &remoteBucket, int awi, void *source){
  //TreeWalk *tw;
  Compute *c;
  State *state;

  CkAssert(awi < maxAwi);
  
  // retrieve the activewalk record
  ActiveWalk &a = activeWalks[awi];
  //tw = a.tw;
  c = a.c;
  state = a.s;

#ifdef CHANGA_REFACTOR_INTERLIST_PRINT_BUCKET_START_FIN
  if(source){
    int start, end;

    GenericTreeNode *testSource = (GenericTreeNode *)source;
    getBucketsBeneathBounds(testSource, start, end);


    //int targetBucket = decodeReqID(reqID);
    //int start1, end1;
    //GenericTreeNode *lowestNode = ((DoubleWalkState *)state)->lowestNode;
    //GenericTreeNode *testSource1 = lowestNode;
    //getBucketsBeneathBounds(testSource1, start1, end1);

    CkPrintf("[%d] tgt=%d recv particles %ld source = %ld (%d - %d) buckRem = %d + %d, remChunk = %d\n",thisIndex,decodeReqID(reqID), remoteBucket, testSource->getKey(), start, end, sRemoteGravityState->counterArrays[0][decodeReqID(reqID)], sLocalGravityState->counterArrays[0][decodeReqID(reqID)], sRemoteGravityState->counterArrays[1][chunk]);

      // check whether source contains target
      GenericTreeNode *target = bucketList[decodeReqID(reqID)];
      if(!testSource->contains(target->getKey())){
        CkAbort("source does not contain target\n");
      }
  }
#endif
  // Some sanity checks
  if(awi == interListAwi) {
#if INTERLIST_VER > 0
      ListCompute *lc = dynamic_cast<ListCompute *>(c);
      CkAssert(lc != NULL);
#else
      CkAbort("Using ListCompute in non-list version\n");
#endif
      }
  else if(awi == remoteGravityAwi) {
      GravityCompute *gc = dynamic_cast<GravityCompute *>(c);
      CkAssert(gc != NULL);
      }
  else if(awi == prefetchAwi) {
      PrefetchCompute *pc = dynamic_cast<PrefetchCompute *>(c);
      CkAssert(pc != NULL);
      }
  else {
      CkAssert(0);
      }
  
  c->reassoc(source, activeRung, a.o);
  //CkPrintf("[%d] RECVD PARTICLES (%ld), chunk %d\n", thisIndex, remoteBucket, chunk);
  c->recvdParticles(egp,num,chunk,reqID,state,this, remoteBucket);
}

void TreePiece::receiveParticlesFullCallback(GravityParticle *gp, int num,
					     int chunk, int reqID,
					     Tree::NodeKey &remoteBucket,
					     int awi, void *source){
  //TreeWalk *tw;
  Compute *c;
  State *state;

  // retrieve the activewalk record
  ActiveWalk &a = activeWalks[awi];
  //tw = a.tw;
  c = a.c;
  state = a.s;

  c->reassoc(source, activeRung, a.o);
  c->recvdParticlesFull(gp,num,chunk,reqID,state,this, remoteBucket);
}

void TreePiece::addActiveWalk(int iAwi, TreeWalk *tw, Compute *c, Opt *o,
			     State *s){
    activeWalks.insert(iAwi, ActiveWalk(tw,c,o,s));
}

void TreePiece::freeWalkObjects(){
#ifdef CHECK_WALK_COMPLETIONS
  CkPrintf("[%d]  TreePiece::freeWalkObjects\n", thisIndex);
#endif

  if(sTopDown){
    delete sTopDown;
    sTopDown = NULL;
  }

#if INTERLIST_VER > 0

  if(sInterListWalk) {
    delete sInterListWalk;
    sInterListWalk = NULL;
  }
#endif

  if(sGravity){
      sGravity->reassoc(0,0,sRemote);
#if INTERLIST_VER > 0 && defined CUDA
      {
        DoubleWalkState *state = (DoubleWalkState *) sRemoteGravityState;
        DoubleWalkState *rstate = (DoubleWalkState *) sInterListStateRemoteResume;
        delete state->particles;
        delete rstate->nodes;
        delete rstate->particles;
        if(!largePhase()){
          DoubleWalkState *lstate = (DoubleWalkState *) sLocalGravityState;
          delete lstate->particles;
          delete [] bucketActiveInfo;
        }
      }
#endif
    // remote-no-resume state
    sGravity->freeState(sRemoteGravityState);
    sRemoteGravityState = NULL;

#if INTERLIST_VER > 0
    // remote-resume state
    // overwrite copies of counters shared with sRemoteGravityState to
    // avoid double deletion.  See startGravity()
    sInterListStateRemoteResume->counterArrays[0] = NULL;
    sInterListStateRemoteResume->counterArrays[1] = NULL;
    sGravity->freeState(sInterListStateRemoteResume);
    sInterListStateRemoteResume = NULL;
#endif

    // local state
    sGravity->reassoc(0,0,sLocal);
    sGravity->freeState(sLocalGravityState);
    sLocalGravityState = NULL;

    delete sGravity;
    sGravity = NULL;
  }

  if(sPrefetch) {
    sPrefetch->freeState(sPrefetchState);
    delete sPrefetch;
    delete sPref;
    delete sLocal;
    delete sRemote;
    sPrefetch = NULL;
  }
}

void TreePiece::finishedChunk(int chunk){
  sRemoteGravityState->numPendingChunks--;
  if(sRemoteGravityState->numPendingChunks == 0){
#ifdef CHECK_WALK_COMPLETIONS
    CkPrintf("[%d] finishedChunk %d, calling markWalkDone\n", thisIndex, chunk);
#endif
    markWalkDone();
  }
}

void TreePiece::markWalkDone() {
    
    // There are always two walks to wait for: one saying that all
    // computation for our buckets are done and one saying that all
    // outstanding cache requests are satisfied.

    if (++completedActiveWalks == 2) {
	// At this point this treepiece has completed its walk.  However,
	// there may be outstanding requests by other pieces.  We need to
	// wait for all walks to complete before freeing data structures.
#ifdef CHECK_WALK_COMPLETIONS
        CkPrintf("[%d] inside markWalkDone, completedActiveWalks: %d, activeWalks: %d, contrib finishWalk\n", thisIndex, completedActiveWalks, activeWalks.size());
#endif
	CkCallback cb = CkCallback(CkIndex_TreePiece::finishWalk(), pieces);
	gravityProxy[thisIndex].ckLocal()->contribute(cb);
	}
    }

void TreePiece::finishWalk()
{
  if(verbosity > 1)
      CkPrintf("[%d] current load: %g current particles: %d\n", thisIndex,
	       getObjTime(), myNumParticles);
  completedActiveWalks = 0;
  freeWalkObjects();
#ifdef CACHE_MEM_STATS
  memPostCache = CmiMemoryUsage()/(1024*1024);
#endif
#ifdef CHECK_WALK_COMPLETIONS
  CkPrintf("[%d] inside finishWalk contrib callback\n", thisIndex);
#endif

#ifdef CUDA_INSTRUMENT_WRS
  RequestTimeInfo *rti1 = hapi_queryInstrument(instrumentId, DM_TRANSFER_LOCAL, activeRung);
  RequestTimeInfo *rti2 = hapi_queryInstrument(instrumentId, DM_TRANSFER_REMOTE_CHUNK, activeRung);
  RequestTimeInfo *rti3 = hapi_queryInstrument(instrumentId, DM_TRANSFER_BACK, activeRung);
  RequestTimeInfo *rti4 = hapi_queryInstrument(instrumentId, DM_TRANSFER_FREE_LOCAL, activeRung);
  RequestTimeInfo *rti5 = hapi_queryInstrument(instrumentId, DM_TRANSFER_FREE_REMOTE_CHUNK, activeRung);
  
  RequestTimeInfo *rti6 = hapi_queryInstrument(instrumentId, TP_GRAVITY_LOCAL, activeRung);
  RequestTimeInfo *rti7 = hapi_queryInstrument(instrumentId, TP_GRAVITY_REMOTE, activeRung);
  RequestTimeInfo *rti8 = hapi_queryInstrument(instrumentId, TP_GRAVITY_REMOTE_RESUME, activeRung);
  
  RequestTimeInfo *rti9 = hapi_queryInstrument(instrumentId,TP_PART_GRAVITY_LOCAL_SMALLPHASE, activeRung);
  RequestTimeInfo *rti10 = hapi_queryInstrument(instrumentId, TP_PART_GRAVITY_LOCAL, activeRung);
  RequestTimeInfo *rti11 = hapi_queryInstrument(instrumentId, TP_PART_GRAVITY_REMOTE, activeRung);
  RequestTimeInfo *rti12 = hapi_queryInstrument(instrumentId, TP_PART_GRAVITY_REMOTE_RESUME, activeRung);

  RequestTimeInfo *rti13 = hapi_queryInstrument(instrumentId, TOP_EWALD_KERNEL, activeRung);
  RequestTimeInfo *rti14 = hapi_queryInstrument(instrumentId, BOTTOM_EWALD_KERNEL, activeRung);

  if(rti6 != NULL){
    CkPrintf("[%d] (%d) CUDA_INSTRUMENT_WRS localnode: (%f,%f,%f) count: %d\n", thisIndex, activeRung, 
        rti6->transferTime/rti6->n,
        rti6->kernelTime/rti6->n,
        rti6->cleanupTime/rti6->n, rti6->n);
  }
  if(rti7 != NULL){
    CkPrintf("[%d] (%d) CUDA_INSTRUMENT_WRS remotenode: (%f,%f,%f) count: %d\n", thisIndex, activeRung, 
        rti7->transferTime/rti7->n,
        rti7->kernelTime/rti7->n,
        rti7->cleanupTime/rti7->n, rti7->n);
  }
  if(rti8 != NULL){
    CkPrintf("[%d] (%d) CUDA_INSTRUMENT_WRS remoteresumenode: (%f,%f,%f) count: %d\n", thisIndex, activeRung, 
        rti8->transferTime/rti8->n,
        rti8->kernelTime/rti8->n,
        rti8->cleanupTime/rti8->n, rti8->n);
  }

  if(rti10 != NULL){
    CkPrintf("[%d] (%d) CUDA_INSTRUMENT_WRS localpart: (%f,%f,%f) count: %d\n", thisIndex, activeRung, 
        rti10->transferTime/rti10->n,
        rti10->kernelTime/rti10->n,
        rti10->cleanupTime/rti10->n, rti10->n);
  }
  if(rti11 != NULL){
    CkPrintf("[%d] (%d) CUDA_INSTRUMENT_WRS remotepart: (%f,%f,%f) count: %d\n", thisIndex, activeRung, 
        rti11->transferTime/rti11->n,
        rti11->kernelTime/rti11->n,
        rti11->cleanupTime/rti11->n, rti11->n);
  }
  if(rti12 != NULL){
    CkPrintf("[%d] (%d) CUDA_INSTRUMENT_WRS remoteresumepart: (%f,%f,%f) count: %d\n", thisIndex, activeRung, 
        rti12->transferTime/rti12->n,
        rti12->kernelTime/rti12->n,
        rti12->cleanupTime/rti12->n, rti12->n);
  }

  if(nLocalNodeReqs > 0){
    CkPrintf("[%d] (%d) CUDA_INSTRUMENT_WRS construction local node reqs: %d, avg: %f\n", 
              thisIndex, 
              activeRung,
              nLocalNodeReqs,
              localNodeListConstructionTime/nLocalNodeReqs
              );
  }
  if(nRemoteNodeReqs > 0){
    CkPrintf("[%d] (%d) CUDA_INSTRUMENT_WRS construction remote node reqs: %d, avg: %f\n", 
              thisIndex, 
              activeRung,
              nRemoteNodeReqs,
              remoteNodeListConstructionTime/nRemoteNodeReqs
              );
  }
  if(nRemoteResumeNodeReqs > 0){
    CkPrintf("[%d] (%d) CUDA_INSTRUMENT_WRS construction remote resume node reqs: %d, avg: %f\n", 
              thisIndex, 
              activeRung,
              nRemoteResumeNodeReqs,
              remoteResumeNodeListConstructionTime/nRemoteResumeNodeReqs
              );
  }
  if(nLocalPartReqs > 0){
    CkPrintf("[%d] (%d) CUDA_INSTRUMENT_WRS construction local part reqs: %d, avg: %f\n", 
              thisIndex, 
              activeRung,
              nLocalPartReqs,
              localPartListConstructionTime/nLocalPartReqs
              );
  }
  if(nRemotePartReqs > 0){
    CkPrintf("[%d] (%d) CUDA_INSTRUMENT_WRS construction remote part reqs: %d, avg: %f\n", 
              thisIndex, 
              activeRung,
              nRemotePartReqs,
              remotePartListConstructionTime/nRemotePartReqs
              );
  }
  if(nRemoteResumePartReqs > 0){
    CkPrintf("[%d] (%d) CUDA_INSTRUMENT_WRS construction remote resume part reqs: %d, avg: %f\n", 
              thisIndex, 
              activeRung,
              nRemoteResumePartReqs,
              remoteResumePartListConstructionTime/nRemoteResumePartReqs
              );
  }

#endif
#ifdef CUDA_STATS
  CkPrintf("[%d] (%d) CUDA_STATS localnode: %ld\n", thisIndex, activeRung, localNodeInteractions);
  CkPrintf("[%d] (%d) CUDA_STATS remotenode: %ld\n", thisIndex, activeRung, remoteNodeInteractions);
  CkPrintf("[%d] (%d) CUDA_STATS remoteresumenode: %ld\n", thisIndex, activeRung, remoteResumeNodeInteractions);
  CkPrintf("[%d] (%d) CUDA_STATS localpart: %ld\n", thisIndex, activeRung, localPartInteractions);
  CkPrintf("[%d] (%d) CUDA_STATS remotepart: %ld\n", thisIndex, activeRung, remotePartInteractions);
  CkPrintf("[%d] (%d) CUDA_STATS remoteresumepart: %ld\n", thisIndex, activeRung, remoteResumePartInteractions);
  
#endif

  gravityProxy[thisIndex].ckLocal()->contribute(cbGravity);
}

#if INTERLIST_VER > 0
/// \brief get range of bucket numbers beneath a given TreeNode.
/// \param source Given TreeNode
/// \param start Index of first bucket (returned)
/// \param end Index of last bucket (returned)
void TreePiece::getBucketsBeneathBounds(GenericTreeNode *&source, int &start, int &end){
  start = source->startBucket;
  end = start+(source->numBucketsBeneath);
}

// called when missed data is received
void TreePiece::updateBucketState(int start, int end, int n, int chunk, State *state){
#if COSMO_PRINT_BK > 1 
  CkPrintf("[%d] data received book-keep\n", thisIndex);
#endif
  for(int i = start; i < end; i++){
    if(bucketList[i]->rungs >= activeRung){
#if !defined CELL 
       state->counterArrays[0][i] -= n;
#if COSMO_PRINT_BK > 1
       CkPrintf("[%d] bucket %d numAddReq: %d,%d\n", thisIndex, i, sRemoteGravityState->counterArrays[0][i], sLocalGravityState->counterArrays[0][i]);
#endif
       //CkPrintf("[%d] bucket %d numAddReq: %d\n", thisIndex, i, state->counterArrays[0][i]);
#if !defined CUDA
       // updatebucketstart is only called after we have finished
       // with received nodes/particles (i.e. after stateReady in
       // either case.) therefore, some work requests must already
       // have been issued before it was invoked, so there will
       // be a finishBucket call afterwards to ensure progress.
       finishBucket(i);
#endif
#endif
    }
  }
  state->counterArrays[1][chunk] -= n;
  //CkPrintf("[%d] updateBucketState (recvd) remaining: %d\n", thisIndex, state->counterArrays[1][chunk]);
  //addMisses(-n);
  //CkPrintf("- misses: %d\n", getMisses());
}

// called on a miss
void TreePiece::updateUnfinishedBucketState(int start, int end, int n, int chunk, State *state){
#if COSMO_PRINT_BK > 1 
  CkPrintf("[%d] data missed book-keep (%d)\n", thisIndex, chunk);
#endif
  for(int i = start; i < end; i++){
    if(bucketList[i]->rungs >= activeRung){
       state->counterArrays[0][i] += n;
#if COSMO_PRINT_BK > 1
       CkPrintf("[%d] bucket %d numAddReq: %d,%d\n", thisIndex, i, sRemoteGravityState->counterArrays[0][i], sLocalGravityState->counterArrays[0][i]);
#endif
       //CkPrintf("[%d] bucket %d numAddReq: %d\n", thisIndex, i, state->counterArrays[0][i]);
    }
  }
  state->counterArrays[1][chunk] += n;
  //CkPrintf("[%d] updateBucketState (missed) remaining: %d\n", thisIndex, state->counterArrays[1][chunk]);
  //addMisses(n);
  //CkPrintf("+ misses: %d\n", getMisses());
}
#endif // INTERLIST_VER


/* create this chare's portion of the image
   We're lazy, so use the existing DumpFrame framework
*/
 void TreePiece::liveVizDumpFrameInit(liveVizRequestMsg *msg) 
 {
   
   savedLiveVizMsg = msg;
   
   if(thisIndex == 0) {
        ckerr << "Calling liveViz setup on main chare" << endl;
	mainChare.liveVizImagePrep(msg);
   }
 }

/*
 * Utility to turn projections on or off
 * bOn == True => turn on.
 */
void TreePiece::setProjections(int bOn)
{
    if(bOn)
        traceBegin();
    else
        traceEnd();
}

/*
 * Gather memory use statistics recorded during the tree walk
 */
void TreePiece::memCacheStats(const CkCallback &cb) 
{
    int memOut[4];
    memOut[0] = memWithCache;
    memOut[1] = memPostCache;
    memOut[2] = nNodeCacheEntries;
    memOut[3] = nPartCacheEntries;
    contribute(4*sizeof(int), memOut, CkReduction::max_int, cb);
}

void TreePiece::balanceBeforeInitialForces(CkCallback &cb){
  LDObjHandle handle = myRec->getLdHandle();
  LBDatabase *lbdb = LBDatabaseObj();
  int nlbs = lbdb->getNLoadBalancers(); 

  if(nlbs == 0) { // no load balancers.  Skip this
      contribute(cb);
      return;
      }


  Vector3D<float> centroid;
  centroid.x = 0.0;
  centroid.y = 0.0;
  centroid.z = 0.0;
  for(int i = 1; i <= myNumParticles; i++){
    centroid += myParticles[i].position; 
  }
  if(myNumParticles > 0){
    centroid /= myNumParticles;
  }

  TaggedVector3D tv(centroid, handle, myNumParticles, myNumParticles, 0, 0);
  tv.tp = thisIndex;

  /*
  CkPrintf("[%d] centroid %f %f %f\n", 
                      thisIndex,
                      tv.vec.x,
                      tv.vec.y,
                      tv.vec.z
                      );
  */

  string msname("MultistepLB");
  string orb3dname("Orb3dLB");
  string ms_notoponame("MultistepLB_notopo");
  string orb3d_notoponame("Orb3dLB_notopo");
  string msorb_name("MultistepOrbLB");

  BaseLB **lbs = lbdb->getLoadBalancers();
  int i;
  if(foundLB == Null){
    for(i = 0; i < nlbs; i++){
      if(msname == string(lbs[i]->lbName())){ 
      	proxy = lbs[i]->getGroupID();
        foundLB = Multistep;
	proxy = lbs[i]->getGroupID();
        break;
      }
      else if(orb3dname == string(lbs[i]->lbName())){ 
      	proxy = lbs[i]->getGroupID();
        foundLB = Orb3d;
	proxy = lbs[i]->getGroupID();
        break;
      }
     else if(ms_notoponame == string(lbs[i]->lbName())){ 
        proxy = lbs[i]->getGroupID();
        foundLB = Multistep_notopo;
        break;
      }
     else if(msorb_name == string(lbs[i]->lbName())){ 
        proxy = lbs[i]->getGroupID();
        foundLB = MultistepOrb;
        break;
      }
      else if(orb3d_notoponame == string(lbs[i]->lbName())){
      	proxy = lbs[i]->getGroupID();
        foundLB = Orb3d_notopo;
	proxy = lbs[i]->getGroupID();
        break;
      }
    }
  }

  if(foundLB == Multistep){
    CkCallback lbcb(CkIndex_MultistepLB::receiveCentroids(NULL), 0, proxy);
    contribute(sizeof(TaggedVector3D), (char *)&tv, CkReduction::concat, lbcb);
  }
  else if(foundLB == Orb3d){
    CkCallback lbcb(CkIndex_Orb3dLB::receiveCentroids(NULL), 0, proxy);
    contribute(sizeof(TaggedVector3D), (char *)&tv, CkReduction::concat, lbcb);
  }
  else if(foundLB == Multistep_notopo){
    CkCallback lbcb(CkIndex_MultistepLB_notopo::receiveCentroids(NULL), 0, proxy);
    contribute(sizeof(TaggedVector3D), (char *)&tv, CkReduction::concat, lbcb);
  }
  else if(foundLB == Orb3d_notopo){
    CkCallback lbcb(CkIndex_Orb3dLB_notopo::receiveCentroids(NULL), 0, proxy);
    contribute(sizeof(TaggedVector3D), (char *)&tv, CkReduction::concat, lbcb);
  }
  else if(foundLB == MultistepOrb){
    CkCallback lbcb(CkIndex_MultistepOrbLB::receiveCentroids(NULL), 0, proxy);
    contribute(sizeof(TaggedVector3D), (char *)&tv, CkReduction::concat, lbcb);
  }
  else if(foundLB == Null){ 
    // none of the balancers requiring centroids found; go
    // straight to AtSync()
    // At the moment there is no load data: we would have to supply
    // some load data (like particle count) if we want this to work well.
    //  doAtSync();
      contribute(cb);
      return;
  }
  // this will be called in resumeFromSync()
  callback = cb;
}

// Choose a piece from among the owners from which to
// request moments in such a way that they are evenly distributed over
// the owners.
int TreePiece::getResponsibleIndex(int first, int last){
    int which = first + (thisIndex % (1 + last - first));
    if(verbosity > 3) 
	CkPrintf("[tp %d] choosing responsible index %d from %d to %d\n",
		 thisIndex,
		 dm->responsibleIndex[which], first, last);
  return dm->responsibleIndex[which];
}

std::map<NodeKey,NonLocalMomentsClientList>::iterator TreePiece::createTreeBuildMomentsEntry(GenericTreeNode *pickedNode){
  std::pair<std::map<NodeKey,NonLocalMomentsClientList>::iterator,bool> ret;
  ret = nonLocalMomentsClients.insert(make_pair(pickedNode->getKey(),NonLocalMomentsClientList(pickedNode)));
  CkAssert(ret.second);
  return ret.first;
}

/*
void TreePiece::addTreeBuildMomentsClient(GenericTreeNode *targetNode, TreePiece *client, GenericTreeNode *clientNode){
  nonLocalMomentsClients[targetNode->getKey()].addClient(NonLocalMomentsClient(client,clientNode));
}
*/

#ifdef CUDA
void TreePiece::clearMarkedBuckets(CkVec<GenericTreeNode *> &markedBuckets){
  int len = markedBuckets.length();
  for(int i = 0; i < len; i++){
    markedBuckets[i]->bucketArrayIndex = -1;
  }
}

void TreePiece::clearMarkedBucketsAll(){
  for(int i = 0; i < numBuckets; i++){
    bucketList[i]->bucketArrayIndex = -1;
  }
}

#endif

#ifdef REDUCTION_HELPER
ReductionHelper::ReductionHelper(){
}

ReductionHelper::ReductionHelper(CkMigrateMessage *m) : CBase_ReductionHelper(m){
}

void ReductionHelper::pup(PUP::er &p){
    CBase_ReductionHelper::pup(p);
}

void ReductionHelper::countTreePieces(const CkCallback &cb){
  // count the number of tree pieces on this PE
  senseLocalTreePieces();
  // no tree pieces can have called "reduce()" at this point
  numTreePiecesCheckedIn = 0;

  //CkPrintf("ReductionHelper %d counted %d pieces on PE\n", CkMyPe(), myNumTreePieces);

  contribute(cb);
}

void ReductionHelper::reduceBinCounts(int nBins, int64_t *binCounts, const CkCallback &cb){
   numTreePiecesCheckedIn++;
  
  //CkPrintf("ReductionHelper %d recvd %d/%d contributions\n", CkMyPe(), numTreePiecesCheckedIn, myNumTreePieces);

  if(numTreePiecesCheckedIn == 1){
    // resize bin counts vector
    myBinCounts.resize(nBins);
    // initialize counts to contribution to reduction from first tree piece
    memcpy(&myBinCounts[0], binCounts, nBins*sizeof(int64_t));
  }
  else{
    CkAssert(nBins == myBinCounts.size());
    for(int i = 0; i < nBins; i++){
      myBinCounts[i] += binCounts[i];
    }
  }

  // is it time to contribute to PE-wide reduction yet?
  if(numTreePiecesCheckedIn == localTreePieces.presentTreePieces.size()){
    //CkPrintf("ReductionHelper %d contributing to PE-level reduction\n", CkMyPe());
    numTreePiecesCheckedIn = 0;
    contribute(sizeof(int64_t)*myBinCounts.size(), &myBinCounts[0], CkReduction::sum_long, cb);
  }

}

void ReductionHelper::reduceBinLoads(int nBins, double *binLoads, const CkCallback &cb){

  numTreePiecesCheckedIn++;

  if (numTreePiecesCheckedIn == 1) {

    myBinLoads.resize(nBins);
    // initialize to contribution to reduction from first tree piece
    memcpy(&myBinLoads[0], binLoads, nBins*sizeof(double));
  }
  else{
    CkAssert(nBins == myBinLoads.size());
    for(int i = 0; i < nBins; i++){
      myBinLoads[i] += binLoads[i];
    }
  }

  // is it time to contribute to PE-wide reduction yet?
  if(numTreePiecesCheckedIn == localTreePieces.presentTreePieces.size()){
    //CkPrintf("ReductionHelper %d contributing to PE-level reduction\n", CkMyPe());
    numTreePiecesCheckedIn = 0;
    contribute(sizeof(double)*myBinLoads.size(), &myBinLoads[0], CkReduction::sum_double, cb);
  }
}

void TreePieceCounter::addLocation(CkLocation &loc){
  const int *indexData = loc.getIndex().data();
  TreePiece *tp = treeProxy[indexData[0]].ckLocal();
  presentTreePieces.push_back(tp);
}

void TreePieceCounter::reset() {
  presentTreePieces.resize(0);
}

void ReductionHelper::senseLocalTreePieces(){
  localTreePieces.reset();                          
  CkLocMgr *mgr = treeProxy.ckLocMgr();        
  mgr->iterate(localTreePieces);              
}

#endif


