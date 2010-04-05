#include "ParallelGravity.h"
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
    prefetchReq = new OrientedBox<double>[2];
  case Oct_dec:
  case ORB_dec:
    if (numPrefetchReq == 0) {
      numPrefetchReq = 1;
      prefetchReq = new OrientedBox<double>[1];
    }
    break;
  default:
      CkAbort("Invalid domain decomposition requested");
  }

  unsigned int *startParticles;
  unsigned int *numParticlesChunk;
  unsigned int excess;
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

  /* At this point myNumParticles contain the number of particles to be loaded
     into this processro, startParticles and numParticlesChunk are newly
     allocated array containing the first particle and the count of each
     contiguous chunk of particles that has to be loaded. */

  // allocate an array for myParticles
  myParticles = new GravityParticle[myNumParticles + 2];
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
      }
    }
    CkAssert(myPart == myNumParticles);
  }
	
  xdr_destroy(&xdrs);
  fclose(infile);
	
  if(verbosity > 3)
    ckerr << "TreePiece " << thisIndex << ": Read in masses and positions" << endl;
	
  bLoaded = 1;

  if(domainDecomposition!=ORB_dec){
    sort(myParticles+1, myParticles+myNumParticles+1);
  }
  contribute(0, 0, CkReduction::concat, cb);
}

void TreePiece::loadTipsy(const std::string& filename,
			  const double dTuFac, // Convert Temperature
			  const CkCallback& cb) {
	callback = cb;
	
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
	nTotalStars = tipsyHeader.nstar;
	dStartTime = tipsyHeader.time;
	int excess;
	unsigned int startParticle;

	switch (domainDecomposition) {
	case SFC_dec:
        case SFC_peano_dec:
	case SFC_peano_dec_3D:
	case SFC_peano_dec_2D:
	    numPrefetchReq = 2;
	    prefetchReq = new OrientedBox<double>[2];
	case Oct_dec:
	case ORB_dec:
	    if (numPrefetchReq == 0) {
		numPrefetchReq = 1;
		prefetchReq = new OrientedBox<double>[1];
		}
	    break;
	default:
	    CkAbort("Invalid domain decomposition requested");
	    }

	myNumParticles = nTotalParticles / numTreePieces;
    
	excess = nTotalParticles % numTreePieces;
	startParticle = myNumParticles * thisIndex;
	if(thisIndex < (int) excess) {
	    myNumParticles++;
	    startParticle += thisIndex;
	    }
	else {
	    startParticle += excess;
	    }
	
	if(verbosity > 2)
		cerr << thisIndex << ": TreePiece: Taking " << myNumParticles
		     << " of " << nTotalParticles
		     << " particles, starting at " << startParticle << endl;

	// allocate an array for myParticles
	myParticles = new GravityParticle[myNumParticles + 2];
	if(startParticle < nTotalSPH) {
	    if(startParticle + myNumParticles <= nTotalSPH)
		myNumSPH = myNumParticles;
	    else
		myNumSPH = nTotalSPH - startParticle;
	    }
	else {
	    myNumSPH = 0;
	    }
	mySPHParticles = new extraSPHData[myNumSPH];
	
	if(!r.seekParticleNum(startParticle)) {
		CkAbort("Couldn't seek to my particles!");
		return;
		}
	
	Tipsy::gas_particle gp;
	Tipsy::dark_particle dp;
	Tipsy::star_particle sp;

	int iSPH = 0;
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
			myParticles[i+1].iType = TYPE_STAR;
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

// Perform Parallel Scan to establish start of parallel writes, then
// do the writing
void TreePiece::setupWrite(int iStage, // stage of scan
			   u_int64_t iPrevOffset,
			   const std::string& filename,
			   const double dTime,
			   const double dvFac,
			   const double duTFac,
			   const CkCallback& cb)
{
    if(iStage > nSetupWriteStage + 1) {
	// requeue message
	pieces[thisIndex].setupWrite(iStage, iPrevOffset, filename,
				     dTime, dvFac, duTFac, cb);
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
					     cb);
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
						 duTFac, cb);
	    }
	if(thisIndex == (int) numTreePieces-1)
	    assert(nStartWrite+myNumParticles == nTotalParticles);
	nSetupWriteStage = -1;	// reset for next time.
	writeTipsy(filename, dTime, dvFac, duTFac);
	contribute(0, 0, CkReduction::concat, cb);
	}
    }

// Serial output of tipsy file.
void TreePiece::serialWrite(const u_int64_t iPrevOffset, // previously written
							 //particles
			    const std::string& filename,  // output file
			    const double dTime,	 // time or expansion
			    const double dvFac,  // velocity conversion
			    const double duTFac, // temperature conversion
			    const CkCallback& cb)
{
    int *piSph = NULL;
    if(myNumSPH > 0)
	piSph = new int[myNumSPH];
    /*
     * Calculate offsets to send across processors
     */
    int iGasOut = 0;
    for(int iPart = 1; iPart <= myNumParticles ; iPart++) {
	if(TYPETest(&myParticles[iPart], TYPE_GAS)) {
	    piSph[iGasOut] = (extraSPHData *)myParticles[iPart].extraData
		- mySPHParticles;
	    iGasOut++;
	    }
	}
    pieces[0].oneNodeWrite(thisIndex, myNumParticles, myNumSPH, myParticles,
			   mySPHParticles, piSph, iPrevOffset, filename, dTime,
			   dvFac, duTFac, cb);
    if(myNumSPH > 0)
	delete [] piSph;
    }

void
TreePiece::oneNodeWrite(int iIndex, // Index of Treepiece
			int iOutParticles, // number of particles
			int iOutSPH, // number of SPH particles
			GravityParticle* particles, // particles to write
			extraSPHData *pGas, // SPH data
			int *piSPH, // SPH data offsets
			const u_int64_t iPrevOffset, // previously written
						    //particles
			const std::string& filename,  // output file
			const double dTime,	 // time or expansion
			const double dvFac,  // velocity conversion
			const double duTFac, // temperature conversion
			const CkCallback& cb)
{
    /*
     * Calculate offsets from across processors
     */
    int iGasOut = 0;
    for(int iPart = 1; iPart <= iOutParticles; iPart++) {
	if(TYPETest(&particles[iPart], TYPE_GAS)) {
	    particles[iPart].extraData = pGas+ piSPH[iGasOut];
	    iGasOut++;
	    }
	}
    /*
     * setup pointers/data for writeTipsy()
     * XXX yes, this is gross.
     */
    int saveNumParticles = myNumParticles;
    GravityParticle *saveMyParticles = myParticles;
    extraSPHData *savemySPHParticles = mySPHParticles;
    myNumParticles = iOutParticles;
    myParticles = particles;
    mySPHParticles = pGas;
    nStartWrite = iPrevOffset;
    writeTipsy(filename, dTime, dvFac, duTFac);
    /*
     * Restore pointers/data
     */
    myNumParticles = saveNumParticles;
    myParticles = saveMyParticles;
    mySPHParticles = savemySPHParticles;
    
    if(iIndex < (numTreePieces - 1))
	treeProxy[iIndex+1].serialWrite(iPrevOffset + iOutParticles, filename,
					dTime, dvFac, duTFac, cb);
    else
	cb.send();  // we are done.
    }
    
void TreePiece::writeTipsy(const std::string& filename, const double dTime,
			   const double dvFac, // scale velocities
			   const double duTFac) // convert temperature
{
    Tipsy::header tipsyHeader;

    tipsyHeader.time = dTime;
    tipsyHeader.nbodies = nTotalParticles;
    tipsyHeader.nsph = nTotalSPH;
    tipsyHeader.nstar = nTotalStars;
    tipsyHeader.ndark = nTotalParticles - (nTotalSPH + nTotalStars);
    
    Tipsy::TipsyWriter w(filename, tipsyHeader);
    
    if(thisIndex == 0)
	w.writeHeader();
    if(!w.seekParticleNum(nStartWrite))
	CkAbort("bad seek");
    for(unsigned int i = 0; i < myNumParticles; i++) {
	if(myParticles[i+1].iOrder < tipsyHeader.nsph) {
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
	    gp.temp = duTFac*myParticles[i+1].u();

	    if(!w.putNextGasParticle(gp))
		CkAbort("Bad Write");
	    }
	else if(myParticles[i+1].iOrder < tipsyHeader.nsph
		      + tipsyHeader.ndark) {
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
	else {
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

	    if(!w.putNextStarParticle(sp))
		CkAbort("Bad Write");
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
/// @brief Reorder particles for output
///
void TreePiece::reOrder(CkCallback& cb)
{
    callback = cb;
    int64_t *startParticle = new int64_t[numTreePieces+1];
    int iPiece;
    
    //
    // Calculate iOrder boundaries for all processors
    // @TODO: assumes no particle creation/destruction
    //
    for(iPiece = 0; iPiece < numTreePieces; iPiece++) {
	int nOutParticles = nTotalParticles/ numTreePieces;
	int excess = nTotalParticles % numTreePieces;
	startParticle[iPiece] = nOutParticles * iPiece;
	if(iPiece < (int) excess) {
	    startParticle[iPiece] += iPiece;
	    }
	else {
	    startParticle[iPiece] += excess;
	    }
	}
    startParticle[numTreePieces] = nTotalParticles; // @TODO: replace
						    // with MaxIOrder
						    // for particle
						    // creation/deletion
    // Sort particles in iOrder
    sort(myParticles+1, myParticles+myNumParticles+1, compIOrder);

    // Tag boundary particle to avoid overruns
    myParticles[myNumParticles+1].iOrder = nTotalParticles;
    
    // Loop through sending particles to correct processor.
    GravityParticle *binBegin = &myParticles[1];
    GravityParticle *binEnd;
    for(iPiece = 0; iPiece < numTreePieces; iPiece++) {
	for(binEnd = binBegin; binEnd->iOrder < startParticle[iPiece+1];
	    binEnd++);
	if((binEnd - binBegin) > 0) {
	    int nGasOut = 0;
	    for(GravityParticle *pPart = binBegin; pPart < binEnd; pPart++) {
		if(TYPETest(pPart, TYPE_GAS))
		    nGasOut++;
		}
	    extraSPHData *pGasOut = NULL;
	    if(nGasOut > 0) {
		pGasOut = new extraSPHData[nGasOut];
		int iGasOut = 0;
		for(GravityParticle *pPart = binBegin; pPart < binEnd; pPart++) {
		    if(TYPETest(pPart, TYPE_GAS)) {
			pGasOut[iGasOut] = *(extraSPHData *)pPart->extraData;
			iGasOut++;
			}
		    }
		}
	    if (verbosity>=3)
		CkPrintf("me:%d to:%d how many:%d\n",thisIndex, iPiece,
			 (binEnd-binBegin));
	    if(iPiece == thisIndex) {
		ioAcceptSortedParticles(binBegin, binEnd - binBegin, pGasOut,
					nGasOut);
		}
	    else {
		pieces[iPiece].ioAcceptSortedParticles(binBegin,
				       binEnd - binBegin, pGasOut, nGasOut);
		}
	    if(nGasOut > 0)
		delete pGasOut;
	    }
	if(&myParticles[myNumParticles + 1] <= binEnd)
	    break;
	binBegin = binEnd;
	}

    delete[] startParticle;

    // signify completion
    incomingParticlesSelf = true;
    ioAcceptSortedParticles(binBegin, 0, NULL, 0);
    }

/// Accept particles from other TreePieces once the sorting has finished
void TreePiece::ioAcceptSortedParticles(const GravityParticle* particles,
					const int n, const extraSPHData *pGas,
					const int nGasIn) {

    int myIOParticles = nTotalParticles / numTreePieces;
    
    int excess = nTotalParticles % numTreePieces;
    if(thisIndex < (int) excess) {
	    myIOParticles++;
	    }
  // allocate new particles array on first call
  if (incomingParticles == NULL) {
    incomingParticles = new GravityParticle[myIOParticles + 2];
    incomingGas = new std::vector<extraSPHData>;
  }

  memcpy(&incomingParticles[incomingParticlesArrived+1], particles,
	 n*sizeof(GravityParticle));
  incomingParticlesArrived += n;
  int nLastGas = incomingGas->size();
  incomingGas->resize(nLastGas + nGasIn);
  memcpy(&((*incomingGas)[nLastGas]), pGas, nGasIn*sizeof(extraSPHData));

  assert(incomingParticlesArrived <= myIOParticles);
  
  if(verbosity >= 3) {
      ckerr << thisIndex << ": received " << incomingParticlesArrived << endl;
      }
  
  if(myIOParticles == incomingParticlesArrived && incomingParticlesSelf) {
    //I've got all my particles
    delete[] myParticles;
    myParticles = incomingParticles;
    incomingParticles = NULL;
    myNumParticles = myIOParticles;
    // reset for next time
    incomingParticlesArrived = 0;
    incomingParticlesSelf = false;

    delete[] mySPHParticles;
    myNumSPH = incomingGas->size();
    mySPHParticles = new extraSPHData[myNumSPH];
    memcpy(mySPHParticles, &((*incomingGas)[0]),
	   incomingGas->size()*sizeof(extraSPHData));
    delete incomingGas;

    // assign gas data pointers
    int iGas = 0;
    for(int iPart = 0; iPart < myNumParticles; iPart++) {
	if(TYPETest(&myParticles[iPart+1], TYPE_GAS)) {
	    myParticles[iPart+1].extraData
		= (extraSPHData *)&mySPHParticles[iGas];
	    iGas++;
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

void TreePiece::outputIOrderASCII(const string& fileName, const CkCallback& cb) {
  if(thisIndex==0) {
    if(verbosity > 2)
      ckerr << "TreePiece " << thisIndex << ": Writing header for iOrder file"
	    << endl;
    FILE* outfile = fopen(fileName.c_str(), "w");
    fprintf(outfile,"%d\n",(int) nTotalParticles);
    fclose(outfile);
  }
	
  if(verbosity > 3)
    ckerr << "TreePiece " << thisIndex << ": Writing iOrder to disk" << endl;
	
  FILE* outfile = fopen(fileName.c_str(), "r+");
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
      pieces[thisIndex + 1].outputIOrderASCII(fileName, cb);
  }
