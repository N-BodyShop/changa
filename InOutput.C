#include "ParallelGravity.h"
#include "DataManager.h"
#include "TipsyFile.h"
#include "OrientedBox.h"
#include "Reductions.h"
#include "InOutput.h"
#include <errno.h>

using namespace TypeHandling;
using namespace SFC;
using namespace std;

// Parse NChilada description file
int TreePiece::parseNC(const std::string& fn)
{
    return 0;
    }

void TreePiece::load(const std::string& fn, const CkCallback& cb) {
  // mark presence in the cache if we are in the first iteration (indicated by localCache==NULL)
  LBTurnInstrumentOff();
  
  basefilename = fn;
  bLoaded = 0;

  if(!parseNC(fn)) {	
      contribute(0, 0, CkReduction::concat, cb);
      return;
      }
     
  //read in particles
  XDR xdrs;
  FILE* infile = fopen((basefilename + ".mass").c_str(), "rb");
  if(!infile) {
    contribute(0, 0, CkReduction::concat, cb);
    return;
  }
  xdrstdio_create(&xdrs, infile, XDR_DECODE);
	
  FieldHeader fh;
  
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

  // XXX two problems here: 1. Why do we need to talk about chunks at
  // this point?
  // 2. There is a small memory leak with startParticles and numParticlesChunk
  // 3. Some of this functionality should be moved into a separate
  // function so it can also be used in loadTipsy().

  switch (domainDecomposition) {
  case SFC_dec:
  case SFC_peano_dec:
  case SFC_peano_dec_3D:
  case SFC_peano_dec_2D:
    numPrefetchReq = 2;
  case Oct_dec:
  case ORB_dec:
  case ORB_space_dec:
    numPrefetchReq = 1;
    break;
  default:
      CkAbort("Invalid domain decomposition requested");
  }

  unsigned int numLoadingTreePieces;
  if (domainDecomposition == Oct_dec) {
    numLoadingTreePieces = CkNumPes();
    if (thisIndex >= CkNumPes()) {
      myNumParticles = 0;
      contribute(0, 0, CkReduction::concat, cb);
      return;
    }
  }
  else numLoadingTreePieces = numTreePieces;
  
  unsigned int *startParticles;
  unsigned int *numParticlesChunk;
  unsigned int excess;
    numParticlesChunk = new unsigned int[2];
    startParticles = new unsigned int[1];
    numParticlesChunk[0] = fh.numParticles / numLoadingTreePieces;
    numParticlesChunk[1] = 0; // sentinel for end chunks
    
    excess = fh.numParticles % numLoadingTreePieces;
    startParticles[0] = numParticlesChunk[0] * thisIndex;
    if(thisIndex < (int) excess) {
      numParticlesChunk[0]++;
      startParticles[0] += thisIndex;
    } else {
      startParticles[0] += excess;
    }
    myNumParticles = numParticlesChunk[0];

  /* At this point myNumParticles contain the number of particles to be loaded
     into this processor, startParticles and numParticlesChunk are newly
     allocated array containing the first particle and the count of each
     contiguous chunk of particles that has to be loaded. */

  // allocate an array for myParticles
  nStore = (int)((myNumParticles + 2)*(1.0 + dExtraStore));
  myParticles = new GravityParticle[nStore];
  assert(myParticles != NULL);

  if(thisIndex == 0)
    ckout << " (" << fh.numParticles << ")";

  if(verbosity > 3)
    ckout << "TreePiece " << thisIndex << ": Of " << fh.numParticles << " particles, taking " << startParticles[0] << " through " << (startParticles[0] + numParticlesChunk[0] - 1) << endl;

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
    if((domainDecomposition!=ORB_dec) && (domainDecomposition!=ORB_space_dec)){
      k = generateKey(pos, boundingBox);
    }
    for(u_int64_t i = 0; i < myNumParticles; ++i) {
      myParticles[i + 1].position = pos;
      if((domainDecomposition!=ORB_dec)
	 && (domainDecomposition!=ORB_space_dec)){
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

      Key current;
      //read all my particles' positions and make keys
      for(unsigned int i = 0; i < numParticlesChunk[chunkNum]; ++i) {
	      if(!xdr_template(&xdrs, &pos)) {
	        ckerr << "TreePiece " << thisIndex << ": Problem reading my part of the positions file, aborting" << endl;
	        CkAbort("Badness");
	      }
	      myParticles[++myPart].position = pos;
	      if((domainDecomposition!=ORB_dec)
		 && (domainDecomposition!=ORB_space_dec)){
	        current = generateKey(pos, boundingBox);
	        myParticles[myPart].key = current;
        }
      }
    }
    CkAssert(myPart == myNumParticles);
  }
	
  xdr_destroy(&xdrs);
  fclose(infile);
	
  if(verbosity > 3)
    ckerr << "TreePiece " << thisIndex << ": Read in masses and positions" << endl;
	
  bLoaded = 1;

  if((domainDecomposition!=ORB_dec) && (domainDecomposition!=ORB_space_dec)){
    sort(myParticles+1, myParticles+myNumParticles+1);
  }

  contribute(0, 0, CkReduction::concat, cb);
}

void TreePiece::loadTipsy(const std::string& filename,
			  const double dTuFac, // Convert Temperature
			  const CkCallback& cb) {
	callback = cb;
        CkCallback replyCB(CkIndex_TreePiece::assignKeys(0), pieces);
	
	bLoaded = 0;

	Tipsy::TipsyReader r(filename);
	if(!r.status()) {
		cerr << thisIndex << ": TreePiece: Fatal: Couldn't open tipsy file!" << endl;
		cb.send(0);	// Fire off callback
		return;
	}
	
	Tipsy::header tipsyHeader = r.getHeader();
	nTotalParticles = tipsyHeader.nbodies;
	nTotalSPH = tipsyHeader.nsph;
	nTotalDark = tipsyHeader.ndark;
	nTotalStar = tipsyHeader.nstar;
	dStartTime = tipsyHeader.time;
	int excess;
	unsigned int startParticle;

	switch (domainDecomposition) {
	case SFC_dec:
        case SFC_peano_dec:
	case SFC_peano_dec_3D:
	case SFC_peano_dec_2D:
	    numPrefetchReq = 2;
	case Oct_dec:
	case ORB_dec:
	case ORB_space_dec:
	    numPrefetchReq = 1;
	    break;
	default:
	    CkAbort("Invalid domain decomposition requested");
	    }

        bool skipLoad = false;
        int numLoadingPEs = CkNumPes();

        // check whether this tree piece should load from file
#ifdef ROUND_ROBIN_WITH_OCT_DECOMP 
        if(thisIndex >= numLoadingPEs) skipLoad = true;
#else
        // this is not the best way to divide objects among PEs, 
        // but it is how charm++ does it.
        int numTreePiecesPerPE = numTreePieces/numLoadingPEs;
        int rem = numTreePieces-numTreePiecesPerPE*numLoadingPEs;

        if(rem > 0){
          numTreePiecesPerPE++;
          numLoadingPEs = numTreePieces/numTreePiecesPerPE;
          if(numTreePieces % numTreePiecesPerPE > 0) numLoadingPEs++;
        }

        if(thisIndex % numTreePiecesPerPE > 0) skipLoad = true;
#endif

        /*
        if(thisIndex == 0){
          CkPrintf("[%d] numTreePieces %d numTreePiecesPerPE %d numLoadingPEs %d\n", thisIndex, numTreePieces, numTreePiecesPerPE, numLoadingPEs);
        }
        */


        if(skipLoad){
          myNumParticles = 0;
          contribute(sizeof(OrientedBox<float>), &boundingBox,
                     growOrientedBox_float,
                     replyCB);
          return;
        }

        // find your load offset into input file
        int myIndex = CkMyPe();
	myNumParticles = nTotalParticles / numLoadingPEs;
	excess = nTotalParticles % numLoadingPEs;
	startParticle = myNumParticles * myIndex;
	if(myIndex < (int) excess) {
	    myNumParticles++;
	    startParticle += myIndex;
	    }
	else {
	    startParticle += excess;
	    }

	
	if(verbosity > 2)
		cerr << "TreePiece " << thisIndex << " PE " << CkMyPe() << " Taking " << myNumParticles
		     << " of " << nTotalParticles
		     << " particles: [" << startParticle << "," << startParticle+myNumParticles << ")" << endl;

	// allocate an array for myParticles
	nStore = (int)((myNumParticles + 2)*(1.0 + dExtraStore));
	myParticles = new GravityParticle[nStore];
	// Are we loading SPH?
	if(startParticle < nTotalSPH) {
	    if(startParticle + myNumParticles <= nTotalSPH)
		myNumSPH = myNumParticles;
	    else
		myNumSPH = nTotalSPH - startParticle;
	    }
	else {
	    myNumSPH = 0;
	    }
	nStoreSPH = (int)(myNumSPH*(1.0 + dExtraStore));
	if(nStoreSPH > 0)
	    mySPHParticles = new extraSPHData[nStoreSPH];
	// Are we loading stars?
	if(startParticle + myNumParticles > nTotalSPH + nTotalDark) {
	    if(startParticle <= nTotalSPH + nTotalDark)
		myNumStar = startParticle + myNumParticles
		    - (nTotalSPH + nTotalDark);
	    else
		myNumStar = myNumParticles;
	    }
	else {
	    myNumStar = 0;
	    }
	nStoreStar = (int)(myNumStar*(1.0 + dExtraStore));
	nStoreStar += 12;  // in case we start with 0
	myStarParticles = new extraStarData[nStoreStar];
	
	if(!r.seekParticleNum(startParticle)) {
		CkAbort("Couldn't seek to my particles!");
		return;
		}
	
	Tipsy::gas_particle gp;
	Tipsy::dark_particle dp;
	Tipsy::star_particle sp;

	int iSPH = 0;
	int iStar = 0;
	for(unsigned int i = 0; i < myNumParticles; ++i) {
		if(i + startParticle < (unsigned int) tipsyHeader.nsph) {
			if(!r.getNextGasParticle(gp)) {
			    CkAbort("failed to read gas particle!");
			    }
			myParticles[i+1].mass = gp.mass;
			myParticles[i+1].position = gp.pos;
			myParticles[i+1].velocity = gp.vel;
			myParticles[i+1].soft = gp.hsmooth;
#ifdef CHANGESOFT
			myParticles[i+1].fSoft0 = gp.hsmooth;
#endif
			myParticles[i+1].iType = TYPE_GAS;
			myParticles[i+1].fDensity = gp.rho;
			myParticles[i+1].extraData = &mySPHParticles[iSPH];
			mySPHParticles[iSPH].fMetals() = gp.metals;
			mySPHParticles[iSPH].u() = dTuFac*gp.temp;
			mySPHParticles[iSPH].uPred() = dTuFac*gp.temp;
			mySPHParticles[iSPH].vPred() = gp.vel;
			mySPHParticles[iSPH].fBallMax() = HUGE;
			iSPH++;
		} else if(i + startParticle < (unsigned int) tipsyHeader.nsph
			  + tipsyHeader.ndark) {
			if(!r.getNextDarkParticle(dp)) {
			    CkAbort("failed to read dark particle!");
			    }
			myParticles[i+1].mass = dp.mass;
			myParticles[i+1].position = dp.pos;
			myParticles[i+1].velocity = dp.vel;
			myParticles[i+1].soft = dp.eps;
#ifdef CHANGESOFT
			myParticles[i+1].fSoft0 = dp.eps;
#endif
			myParticles[i+1].iType = TYPE_DARK;
		} else {
			if(!r.getNextStarParticle(sp)) {
			    CkAbort("failed to read star particle!");
			    }
			myParticles[i+1].mass = sp.mass;
			myParticles[i+1].position = sp.pos;
			myParticles[i+1].velocity = sp.vel;
			myParticles[i+1].soft = sp.eps;
#ifdef CHANGESOFT
			myParticles[i+1].fSoft0 = sp.eps;
#endif
			myParticles[i+1].extraData = &myStarParticles[iStar];
			myParticles[i+1].fStarMetals() = sp.metals;
			myParticles[i+1].fTimeForm() = sp.tform;
			myParticles[i+1].iType = TYPE_STAR;
			iStar++;
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
		   replyCB);
}

/// @brief Find starting offsets and begin parallel write.
///
/// Perform Parallel Scan (a log(p) algorithm) to establish start of
/// parallel writes, then do the writing.  Calls writeTipsy() to do
/// the actual writing.
void TreePiece::setupWrite(int iStage, // stage of scan
			   u_int64_t iPrevOffset,
			   const std::string& filename,
			   const double dTime,
			   const double dvFac,
			   const double duTFac,
			   const int bCool,
			   const CkCallback& cb)
{
    if(iStage > nSetupWriteStage + 1) {
	// requeue message
	pieces[thisIndex].setupWrite(iStage, iPrevOffset, filename,
				     dTime, dvFac, duTFac, bCool, cb);
	return;
	}
    nSetupWriteStage++;
    if(iStage == 0)
	nStartWrite = 0;
    
    int iOffset = 1 << nSetupWriteStage;
    nStartWrite += iPrevOffset;
    if(thisIndex+iOffset < (int) numTreePieces) { // Scan on
	if(verbosity > 1) {
	    ckerr << thisIndex << ": stage " << iStage << " sending "
		  << nStartWrite+myNumParticles << " to " << thisIndex+iOffset
		  << endl;
	    }
	pieces[thisIndex+iOffset].setupWrite(iStage+1,
					     nStartWrite+myNumParticles,
					     filename, dTime, dvFac, duTFac,
					     bCool, cb);
	}
    if(thisIndex < iOffset) { // No more messages are coming my way
	// send out all the messages
	for(iStage = iStage+1; (1 << iStage) + thisIndex < (int)numTreePieces;
	    iStage++) {
	    iOffset = 1 << iStage;
	    if(verbosity > 1) {
		ckerr << thisIndex << ": stage " << iStage << " sending "
		      << nStartWrite+myNumParticles << " to "
		      << thisIndex+iOffset << endl;
		}
	    pieces[thisIndex+iOffset].setupWrite(iStage+1,
						 nStartWrite+myNumParticles,
						 filename, dTime, dvFac,
						 duTFac, bCool, cb);
	    }
	if(thisIndex == (int) numTreePieces-1)
	    assert(nStartWrite+myNumParticles == nTotalParticles);
	nSetupWriteStage = -1;	// reset for next time.
	parallelWrite(0, cb, filename, dTime, dvFac, duTFac, bCool);
	}
    }

/// @brief Control the parallelism in the tipsy output by breaking it
/// up into nIOProcessor pieces.
/// @param iPass What pass we are on in the parallel write.  The
/// initial call should be with "0".

void TreePiece::parallelWrite(int iPass, const CkCallback& cb,
			      const std::string& filename, const double dTime,
			      const double dvFac, // scale velocities
			      const double duTFac, // convert temperature
			      const int bCool)
{
    if(nIOProcessor == 0) {	// use them all
	writeTipsy(filename, dTime, dvFac, duTFac, bCool);
	contribute(cb);
	return;
	}
    int nSkip = numTreePieces/nIOProcessor;
    if(nSkip == 0)
	nSkip = 1;
    if(thisIndex%nSkip != iPass) { // N.B. this will only be the case
				   // for the first pass.
	return;
	}
    writeTipsy(filename, dTime, dvFac, duTFac, bCool);
    if(iPass < (nSkip - 1) && thisIndex < (numTreePieces - 1))
	treeProxy[thisIndex+1].parallelWrite(iPass + 1, cb, filename, dTime,
					     dvFac, duTFac, bCool);
    contribute(cb);
    }


// Serial output of tipsy file.
void TreePiece::serialWrite(const u_int64_t iPrevOffset, // previously written
							 //particles
			    const std::string& filename,  // output file
			    const double dTime,	 // time or expansion
			    const double dvFac,  // velocity conversion
			    const double duTFac, // temperature conversion
			    const int bCool,
			    const CkCallback& cb)
{
    int *piSph = NULL;
    if(myNumSPH > 0)
	piSph = new int[myNumSPH];
    int *piStar = NULL;
    if(myNumStar > 0)
	piStar = new int[myNumStar];
    /*
     * Calculate offsets to send across processors
     */
    int iGasOut = 0;
    int iStarOut = 0;
    for(int iPart = 1; iPart <= myNumParticles ; iPart++) {
	if(myParticles[iPart].isGas()) {
	    piSph[iGasOut] = (extraSPHData *)myParticles[iPart].extraData
		- mySPHParticles;
	    iGasOut++;
	    }
	if(myParticles[iPart].isStar()) {
	    piStar[iStarOut] = (extraStarData *)myParticles[iPart].extraData
		- myStarParticles;
	    iStarOut++;
	    }
	}
    pieces[0].oneNodeWrite(thisIndex, myNumParticles, myNumSPH, myNumStar,
			   myParticles, mySPHParticles, myStarParticles,
			   piSph, piStar, iPrevOffset, filename, dTime,
			   dvFac, duTFac, bCool, cb);
    if(myNumSPH > 0)
	delete [] piSph;
    if(myNumStar > 0)
	delete [] piStar;
    }

void
TreePiece::oneNodeWrite(int iIndex, // Index of Treepiece
			int iOutParticles, // number of particles
			int iOutSPH, // number of SPH particles
			int iOutStar, // number of Star particles
			GravityParticle* particles, // particles to write
			extraSPHData *pGas, // SPH data
			extraStarData *pStar, // Star data
			int *piSPH, // SPH data offsets
			int *piStar, // Star data offsets
			const u_int64_t iPrevOffset, // previously written
						    //particles
			const std::string& filename,  // output file
			const double dTime,	 // time or expansion
			const double dvFac,  // velocity conversion
			const double duTFac, // temperature conversion
			const int bCool,
			const CkCallback& cb)
{
    /*
     * Calculate offsets from across processors
     */
    int iGasOut = 0;
    int iStarOut = 0;
    for(int iPart = 1; iPart <= iOutParticles; iPart++) {
	if(particles[iPart].isGas()) {
	    particles[iPart].extraData = pGas+ piSPH[iGasOut];
	    iGasOut++;
	    }
	if(particles[iPart].isStar()) {
	    particles[iPart].extraData = pStar+ piStar[iStarOut];
	    iStarOut++;
	    }
	}
    /*
     * setup pointers/data for writeTipsy()
     * XXX yes, this is gross.
     */
    int saveNumParticles = myNumParticles;
    GravityParticle *saveMyParticles = myParticles;
    extraSPHData *savemySPHParticles = mySPHParticles;
    extraStarData *savemyStarParticles = myStarParticles;
    myNumParticles = iOutParticles;
    myParticles = particles;
    mySPHParticles = pGas;
    myStarParticles = pStar;
    nStartWrite = iPrevOffset;
    writeTipsy(filename, dTime, dvFac, duTFac, bCool);
    /*
     * Restore pointers/data
     */
    myNumParticles = saveNumParticles;
    myParticles = saveMyParticles;
    mySPHParticles = savemySPHParticles;
    myStarParticles = savemyStarParticles;
    
    if(iIndex < (numTreePieces - 1))
	treeProxy[iIndex+1].serialWrite(iPrevOffset + iOutParticles, filename,
					dTime, dvFac, duTFac, bCool, cb);
    else
	cb.send();  // we are done.
    }
    
void TreePiece::writeTipsy(const std::string& filename, const double dTime,
			   const double dvFac, // scale velocities
			   const double duTFac, // convert temperature
			   const int bCool)
{
    Tipsy::header tipsyHeader;

    tipsyHeader.time = dTime;
    tipsyHeader.nbodies = nTotalParticles;
    tipsyHeader.nsph = nTotalSPH;
    tipsyHeader.nstar = nTotalStar;
    tipsyHeader.ndark = nTotalParticles - (nTotalSPH + nTotalStar);
    
    Tipsy::TipsyWriter w(filename, tipsyHeader);
    
    if(thisIndex == 0)
	w.writeHeader();
    if(!w.seekParticleNum(nStartWrite))
	CkAbort("bad seek");
    for(unsigned int i = 0; i < myNumParticles; i++) {
	if(myParticles[i+1].isGas()) {
	    Tipsy::gas_particle gp;
	    gp.mass = myParticles[i+1].mass;
	    gp.pos = myParticles[i+1].position;
	    gp.vel = myParticles[i+1].velocity*dvFac;
#ifdef CHANGESOFT
	    gp.hsmooth = myParticles[i+1].fSoft0;
#else
	    gp.hsmooth = myParticles[i+1].soft;
#endif
	    gp.phi = myParticles[i+1].potential;
	    gp.rho = myParticles[i+1].fDensity;
	    gp.metals = myParticles[i+1].fMetals();
	    if(bCool) {
#ifndef COOLING_NONE
		gp.temp = CoolCodeEnergyToTemperature(dm->Cool,
						      &myParticles[i+1].CoolParticle(),
						      myParticles[i+1].u(),
						      myParticles[i+1].fMetals());
#else
		CkAbort("cooling output without cooling code");
#endif
		}
	    else 
		gp.temp = duTFac*myParticles[i+1].u();

	    if(!w.putNextGasParticle(gp))
		CkAbort("Bad Write");
	    }
	else if(myParticles[i+1].isDark()) {
	    Tipsy::dark_particle dp;
	    dp.mass = myParticles[i+1].mass;
	    dp.pos = myParticles[i+1].position;
	    dp.vel = myParticles[i+1].velocity*dvFac;
#ifdef CHANGESOFT
	    dp.eps = myParticles[i+1].fSoft0;
#else
	    dp.eps = myParticles[i+1].soft;
#endif
	    dp.phi = myParticles[i+1].potential;

	    if(!w.putNextDarkParticle(dp))
		CkAbort("Bad Write");
	    }
	else if(myParticles[i+1].isStar()) {
	    Tipsy::star_particle sp;
	    sp.mass = myParticles[i+1].mass;
	    sp.pos = myParticles[i+1].position;
	    sp.vel = myParticles[i+1].velocity*dvFac;
#ifdef CHANGESOFT
	    sp.eps = myParticles[i+1].fSoft0;
#else
	    sp.eps = myParticles[i+1].soft;
#endif
	    sp.phi = myParticles[i+1].potential;
	    sp.metals = myParticles[i+1].fStarMetals();
	    sp.tform = myParticles[i+1].fTimeForm();

	    if(!w.putNextStarParticle(sp))
		CkAbort("Bad Write");
	    }
	else {
	    CkAbort("Bad particle type in tipsyWrite");
	    }
	}
    }

static bool compIOrder(const GravityParticle& first,
		       const GravityParticle & second)
{
    if(first.iOrder < second.iOrder)
	return true;
    return false;
    }

///
/// @brief Calculate boundaries based on iOrder
///
inline int64_t *iOrderBoundaries(int nPieces, int64_t nMaxOrder) 
{
    int64_t *startParticle = new int64_t[nPieces+1];
    int iPiece;

    // Note that with particle creation/destruction this will not be
    // perfectly balanced.
    //
    for(iPiece = 0; iPiece < nPieces; iPiece++) {
	int nOutParticles = nMaxOrder/ nPieces;
	startParticle[iPiece] = nOutParticles * iPiece;
	}
    startParticle[nPieces] = nMaxOrder+1;
    return startParticle;
    }

///
/// @brief Reorder particles for output
///
void TreePiece::reOrder(int64_t _nMaxOrder, CkCallback& cb)
{
    callback = cb; // Save callback for after shuffle
    int *counts = new int[numTreePieces];
    int iPiece;
    
    nMaxOrder = _nMaxOrder;
    for(iPiece = 0; iPiece < numTreePieces; iPiece++) {
	counts[iPiece] = 0;
	}

    int64_t *startParticle = iOrderBoundaries(numTreePieces, nMaxOrder);

    if (myNumParticles > 0) {
	// Sort particles in iOrder
	sort(myParticles+1, myParticles+myNumParticles+1, compIOrder);

	// Tag boundary particle to avoid overruns
	myParticles[myNumParticles+1].iOrder = nMaxOrder+1;
    
	// Loop through to get particle counts.
	GravityParticle *binBegin = &myParticles[1];
	GravityParticle *binEnd;
	for(iPiece = 0; iPiece < numTreePieces; iPiece++) {
	    for(binEnd = binBegin; binEnd->iOrder < startParticle[iPiece+1];
		binEnd++);
	    int nPartOut = binEnd - binBegin;
	    counts[iPiece] = nPartOut;
	    if(&myParticles[myNumParticles + 1] <= binEnd)
		break;
	    binBegin = binEnd;
	    }
	}
    myIOParticles = -1;
    CkCallback cbShuffle = CkCallback(CkIndex_TreePiece::ioShuffle(NULL),
				      pieces);
    contribute(numTreePieces*sizeof(int), counts, CkReduction::sum_int,
	       cbShuffle);
    delete [] startParticle;
    delete [] counts;
    }

///
/// @brief Perform the shuffle for reOrder
///
void TreePiece::ioShuffle(CkReductionMsg *msg) 
{
    int *counts = (int *)msg->getData();
    myIOParticles = counts[thisIndex];
    delete msg;
    
    int iPiece;
    
    //
    // Calculate iOrder boundaries for all processors
    //
    int64_t *startParticle = iOrderBoundaries(numTreePieces, nMaxOrder);

    if (myNumParticles > 0) {
	// Particles have been sorted in reOrder()
    // Loop through sending particles to correct processor.
    GravityParticle *binBegin = &myParticles[1];
    GravityParticle *binEnd;
    for(iPiece = 0; iPiece < numTreePieces; iPiece++) {
	for(binEnd = binBegin; binEnd->iOrder < startParticle[iPiece+1];
	    binEnd++);
	int nPartOut = binEnd - binBegin;
	if(nPartOut > 0) {
	    int nGasOut = 0;
	    int nStarOut = 0;
	    for(GravityParticle *pPart = binBegin; pPart < binEnd; pPart++) {
		if(pPart->isGas())
		    nGasOut++;
		if(pPart->isStar())
		    nStarOut++;
		}
	    ParticleShuffleMsg *shuffleMsg
		= new (nPartOut, nGasOut, nStarOut)
		ParticleShuffleMsg(nPartOut, nGasOut, nStarOut, 0.0);
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
	    if (verbosity>=3)
		CkPrintf("me:%d to:%d how many:%d\n",thisIndex, iPiece,
			 (binEnd-binBegin));
	    if(iPiece == thisIndex) {
		ioAcceptSortedParticles(shuffleMsg);
		}
	    else {
		pieces[iPiece].ioAcceptSortedParticles(shuffleMsg);
		}
	    }
	if(&myParticles[myNumParticles + 1] <= binEnd)
	    break;
	binBegin = binEnd;
	}
    }
	
    delete[] startParticle;

    // signify completion
    incomingParticlesSelf = true;
    ioAcceptSortedParticles(NULL);
    }

/// Accept particles from other TreePieces once the sorting has finished
void TreePiece::ioAcceptSortedParticles(ParticleShuffleMsg *shuffleMsg) {

    if(shuffleMsg != NULL) {
	incomingParticlesMsg.push_back(shuffleMsg);
	incomingParticlesArrived += shuffleMsg->n;
	}

    if(verbosity > 2)
	ckout << thisIndex << ": incoming: " << incomingParticlesArrived
	      << " myIO: " << myIOParticles << endl;
    
  if(myIOParticles == incomingParticlesArrived && incomingParticlesSelf) {
      //I've got all my particles, now count them
    if(verbosity>1) ckout << thisIndex <<" got ioParticles"
			  <<endl;
    int nTotal = 0;
    int nSPH = 0;
    int nStar = 0;
    int iMsg;
    for(iMsg = 0; iMsg < incomingParticlesMsg.size(); iMsg++) {
	nTotal += incomingParticlesMsg[iMsg]->n;
	nSPH += incomingParticlesMsg[iMsg]->nSPH;
	nStar += incomingParticlesMsg[iMsg]->nStar;
	}
      
    if (myNumParticles > 0) delete[] myParticles;
    nStore = (int) ((nTotal + 2)*(1.0 + dExtraStore));
    myParticles = new GravityParticle[nStore];
    myNumParticles = nTotal;
    // reset for next time
    incomingParticlesArrived = 0;
    incomingParticlesSelf = false;

    myNumSPH = nSPH;
    if (nStoreSPH > 0) delete[] mySPHParticles;
    nStoreSPH = (int) (myNumSPH*(1.0 + dExtraStore));
    mySPHParticles = new extraSPHData[nStoreSPH];

    myNumStar = nStar;
    if(nStoreStar > 0) delete[] myStarParticles;
    nStoreStar = (int) (myNumStar*(1.0 + dExtraStore));
    nStoreStar += 12;
    myStarParticles = new extraStarData[nStoreStar];

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

    // assign gas data pointers
    int iGas = 0;
    int iStar = 0;
    for(int iPart = 0; iPart < myNumParticles; iPart++) {
	if(myParticles[iPart+1].isGas()) {
	    myParticles[iPart+1].extraData
		= (extraSPHData *)&mySPHParticles[iGas];
	    iGas++;
	    }
	if(myParticles[iPart+1].isStar()) {
	    myParticles[iPart+1].extraData
		= (extraStarData *)&myStarParticles[iStar];
	    iStar++;
	    }
	}

    sort(myParticles+1, myParticles+myNumParticles+1, compIOrder);
    //signify completion with a reduction
    if(verbosity>1) ckout << thisIndex <<" contributing to ioAccept particles"
			  <<endl;

    if (root != NULL) {
      root->fullyDelete();
      delete root;
      root = NULL;
      nodeLookupTable.clear();
    }
    contribute(0, 0, CkReduction::concat, callback);
  }
}

void TreePiece::outputAccelerations(OrientedBox<double> accelerationBox, const string& suffix, const CkCallback& cb) {
  FieldHeader fh;
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

// Output a Tipsy ASCII file.
void TreePiece::outputASCII(OutputParams& params, // specifies
						  // filename, format,
						  // and quantity to
						  // be output
			    int bParaWrite,	  // Every processor
						  // can write.  If
						  // false, all output
						  // gets sent to
						  // treepiece "0" for writing.
			    const CkCallback& cb) {
  FILE* outfile;
  double *adOut;	// array for oneNode I/O
  Vector3D<double> *avOut;	// array for one node I/O
  params.dm = dm; // pass cooling information
  
  if((thisIndex==0 && packed) || (thisIndex==0 && !packed && cnt==0)) {
    if(verbosity > 2)
      ckout << "TreePiece " << thisIndex << ": Writing header for output file" << endl;
    outfile = fopen(params.fileName.c_str(), "w");
    CkAssert(outfile != NULL);
    fprintf(outfile,"%d\n",(int) nTotalParticles);
    fclose(outfile);
  }
	
  if(verbosity > 3)
    ckout << "TreePiece " << thisIndex << ": Writing output to disk" << endl;
	
  if(bParaWrite) {
      outfile = fopen(params.fileName.c_str(), "r+");
      if(outfile == NULL)
	    ckerr << "Treepiece " << thisIndex << " failed to open "
		  << params.fileName.c_str() << " : " << errno << endl;
      CkAssert(outfile != NULL);
      int result = fseek(outfile, 0L, SEEK_END);
      CkAssert(result == 0);
      }
  else {
      if(params.bVector && packed)
	  avOut = new Vector3D<double>[myNumParticles];
      else
	  adOut = new double[myNumParticles];
      }
    for(unsigned int i = 1; i <= myNumParticles; ++i) {
      Vector3D<double> vOut;
      double dOut;
      if(params.bVector) {
	  vOut = params.vValue(&myParticles[i]);
	  if(!packed){
	      if(cnt==0)
		  dOut = vOut.x;
	      if(cnt==1)
		  dOut = vOut.y;
	      if(cnt==2)
		  dOut = vOut.z;
	      }
	  }
      else
	  dOut = params.dValue(&myParticles[i]);
      
      if(bParaWrite) {
	  if(params.bVector && packed){
	      if(fprintf(outfile,"%.14g\n",vOut.x) < 0) {
		  ckerr << "TreePiece " << thisIndex << ": Error writing array to disk, aborting" << endl;
		    CkAbort("Badness");
		    }
	      if(fprintf(outfile,"%.14g\n",vOut.y) < 0) {
		  ckerr << "TreePiece " << thisIndex << ": Error writing array to disk, aborting" << endl;
		  CkAbort("Badness");
		  }
	      if(fprintf(outfile,"%.14g\n",vOut.z) < 0) {
		  ckerr << "TreePiece " << thisIndex << ": Error writing array to disk, aborting" << endl;
		  CkAbort("Badness");
		  }
	      }
	  else {
	      if(fprintf(outfile,"%.14g\n",dOut) < 0) {
		  ckerr << "TreePiece " << thisIndex << ": Error writing array to disk, aborting" << endl;
		  CkAbort("Badness");
		  }
	      }
	  }
      else {
	  if(params.bVector && packed)
	      avOut[i-1] = vOut;
	  else
	      adOut[i-1] = dOut;
	  }
      }
  cnt++;
  if(cnt == 3 || !params.bVector)
      cnt = 0;
  
  if(bParaWrite) {
      int result = fclose(outfile);
      if(result != 0)
	    ckerr << "Bad close: " << strerror(errno) << endl;
      CkAssert(result == 0);

      if(thisIndex!=(int)numTreePieces-1) {
	  pieces[thisIndex + 1].outputASCII(params, bParaWrite, cb);
	  return;
	  }

      if(packed || !params.bVector || (!packed && cnt==0)) {
	  cb.send(); // We are done.
	  return;
	  }
      // go through pieces again for unpacked vector.
      pieces[0].outputASCII(params, bParaWrite, cb);
      }
  else {
      int bDone = packed || !params.bVector || (!packed && cnt==0); // flag for last time
      if(params.bVector && packed) {
	  pieces[0].oneNodeOutVec(params, avOut, myNumParticles, thisIndex,
				  bDone, cb);
	  delete [] avOut;
	  }
      else {
	  pieces[0].oneNodeOutArr(params, adOut, myNumParticles, thisIndex,
				  bDone, cb);
	  delete [] adOut;
	  }
      }
}

// Receives an array of vectors to write out in ASCII format
// Assumed to be called from outputASCII() and will continue with the
// next tree piece unless "bDone".
void TreePiece::oneNodeOutVec(OutputParams& params,
			      Vector3D<double>* avOut, // array to be output
			      int nPart, // number of elements in avOut
			      int iIndex, // treepiece which called me
			      int bDone, // Last call
			      CkCallback& cb) 
{
    FILE* outfile = fopen(params.fileName.c_str(), "r+");
    if(outfile == NULL)
	ckerr << "Treepiece " << thisIndex << " failed to open "
	      << params.fileName.c_str() << " : " << errno << endl;
    CkAssert(outfile != NULL);
    int result = fseek(outfile, 0L, SEEK_END);
    CkAssert(result == 0);
    for(int i = 0; i < nPart; ++i) {
	if(fprintf(outfile,"%.14g\n",avOut[i].x) < 0) {
	  ckerr << "TreePiece " << thisIndex << ": Error writing array to disk, aborting" << endl;
	    CkAbort("Badness");
	    }
	if(fprintf(outfile,"%.14g\n",avOut[i].y) < 0) {
	  ckerr << "TreePiece " << thisIndex << ": Error writing array to disk, aborting" << endl;
	  CkAbort("Badness");
	  }
	if(fprintf(outfile,"%.14g\n",avOut[i].z) < 0) {
	  ckerr << "TreePiece " << thisIndex << ": Error writing array to disk, aborting" << endl;
	  CkAbort("Badness");
	  }
	}
    result = fclose(outfile);
    if(result != 0)
	ckerr << "Bad close: " << strerror(errno) << endl;
    CkAssert(result == 0);

    if(iIndex!=(int)numTreePieces-1) {
	  pieces[iIndex + 1].outputASCII(params, 0, cb);
	  return;
	  }

    if(bDone) {
	  cb.send(); // We are done.
	  return;
	  }
    CkAbort("packed array Done logic wrong");
    }

// Receives an array of doubles to write out in ASCII format
// Assumed to be called from outputASCII() and will continue with the
// next tree piece unless "bDone"

void TreePiece::oneNodeOutArr(OutputParams& params,
			      double *adOut, // array to be output
			      int nPart, // length of adOut
			      int iIndex, // treepiece which called me
			      int bDone, // Last call
			      CkCallback& cb) 
{
    FILE* outfile = fopen(params.fileName.c_str(), "r+");
    if(outfile == NULL)
	ckerr << "Treepiece " << thisIndex << " failed to open "
	      << params.fileName.c_str() << " : " << errno << endl;
    CkAssert(outfile != NULL);
    int result = fseek(outfile, 0L, SEEK_END);
    CkAssert(result == 0);
    for(int i = 0; i < nPart; ++i) {
	if(fprintf(outfile,"%.14g\n",adOut[i]) < 0) {
	  ckerr << "TreePiece " << thisIndex << ": Error writing array to disk, aborting" << endl;
	    CkAbort("Badness");
	    }
	}
    result = fclose(outfile);
    if(result != 0)
	ckerr << "Bad close: " << strerror(errno) << endl;
    CkAssert(result == 0);

    if(iIndex!=(int)numTreePieces-1) {
	  pieces[iIndex + 1].outputASCII(params, 0, cb);
	  return;
	  }

    if(bDone) {
	  cb.send(); // We are done.
	  return;
	  }
    // go through pieces again for unpacked vector.
    pieces[0].outputASCII(params, 0, cb);
    }

void TreePiece::outputIntASCII(OutputIntParams& params, // specifies
						  // filename, format,
						  // and quantity to
						  // be output
			    int bParaWrite,	  // Every processor
						  // can write.  If
						  // false, all output
						  // gets sent to
						  // treepiece "0" for writing.
			    const CkCallback& cb) {
  FILE* outfile;
  int *aiOut;	// array for oneNode I/O
  
  if(thisIndex==0) {
    if(verbosity > 2)
      ckout << "TreePiece " << thisIndex << ": Writing header for output file" << endl;
    outfile = fopen(params.fileName.c_str(), "w");
    CkAssert(outfile != NULL);
    fprintf(outfile,"%d\n",(int) nTotalParticles);
    fclose(outfile);
  }
	
  if(verbosity > 3)
    ckout << "TreePiece " << thisIndex << ": Writing output to disk" << endl;
	
  if(bParaWrite) {
      outfile = fopen(params.fileName.c_str(), "r+");
      if(outfile == NULL)
	    ckerr << "Treepiece " << thisIndex << " failed to open "
		  << params.fileName.c_str() << " : " << errno << endl;
      CkAssert(outfile != NULL);
      int result = fseek(outfile, 0L, SEEK_END);
      CkAssert(result == 0);
      }
  else {
      aiOut = new int[myNumParticles];
      }
  for(unsigned int i = 1; i <= myNumParticles; ++i) {
      int iOut;
      iOut = params.iValue(&myParticles[i]);
      
      if(bParaWrite) {
	  if(fprintf(outfile,"%d\n", iOut) < 0) {
	      ckerr << "TreePiece " << thisIndex << ": Error writing array to disk, aborting" << endl;
	      CkAbort("Badness");
	      }
	  }
      else {
	  aiOut[i-1] = iOut;
	  }
      }
     
  if(bParaWrite) {
      int result = fclose(outfile);
      if(result != 0)
	    ckerr << "Bad close: " << strerror(errno) << endl;
      CkAssert(result == 0);

      if(thisIndex!=(int)numTreePieces-1) {
	  pieces[thisIndex + 1].outputIntASCII(params, bParaWrite, cb);
	  return;
	  }

      cb.send(); // We are done.
      return;
      }
  else {
      pieces[0].oneNodeOutIntArr(params, aiOut, myNumParticles, thisIndex, cb);
      delete [] aiOut;
      }
}

// Receives an array of ints to write out in ASCII format
// Assumed to be called from outputIntASCII() and will continue with the
// next tree piece.

void TreePiece::oneNodeOutIntArr(OutputIntParams& params,
			      int *aiOut, // array to be output
			      int nPart, // length of adOut
			      int iIndex, // treepiece which called me
			      CkCallback& cb) 
{
    FILE* outfile = fopen(params.fileName.c_str(), "r+");
    if(outfile == NULL)
	ckerr << "Treepiece " << thisIndex << " failed to open "
	      << params.fileName.c_str() << " : " << errno << endl;
    CkAssert(outfile != NULL);
    int result = fseek(outfile, 0L, SEEK_END);
    CkAssert(result == 0);
    for(int i = 0; i < nPart; ++i) {
	if(fprintf(outfile,"%d\n",aiOut[i]) < 0) {
	  ckerr << "TreePiece " << thisIndex << ": Error writing array to disk, aborting" << endl;
	    CkAbort("Badness");
	    }
	}
    result = fclose(outfile);
    if(result != 0)
	ckerr << "Bad close: " << strerror(errno) << endl;
    CkAssert(result == 0);

    if(iIndex!=(int)numTreePieces-1) {
	  pieces[iIndex + 1].outputIntASCII(params, 0, cb);
	  return;
	  }

    cb.send(); // We are done.
    }
