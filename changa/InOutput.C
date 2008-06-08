#include "ParallelGravity.h"
#include "TipsyFile.h"
#include "OrientedBox.h"
#include "Reductions.h"

using namespace TypeHandling;

// Parse NChilada description file
int TreePiece::parseNC(const std::string& fn)
{
    return 0;
    }

void TreePiece::load(const std::string& fn, const CkCallback& cb) {
  // mark presence in the cache if we are in the first iteration (indicated by localCache==NULL)
  LBTurnInstrumentOff();
  if (_cache && localCache==NULL) {
    localCache = cacheManagerProxy.ckLocalBranch();
    localCache->markPresence(thisIndex, root);
  }

  
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
	tipsyHeader = r.getHeader();
	int totalNumParticles = tipsyHeader.nbodies;
	fh.numParticles = totalNumParticles;
	fh.time = tipsyHeader.time;
	int excess;
	unsigned int startParticle;

	switch (domainDecomposition) {
	case SFC_dec:
        case SFC_peano_dec:
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

	myNumParticles = totalNumParticles / numTreePieces;
    
	excess = totalNumParticles % numTreePieces;
	startParticle = myNumParticles * thisIndex;
	if(thisIndex < (int) excess) {
	    myNumParticles++;
	    startParticle += thisIndex;
	    }
	else {
	    startParticle += excess;
	    }
	
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
		if(i + startParticle < (unsigned int) tipsyHeader.nsph) {
			r.getNextGasParticle(gp);
			myParticles[i+1].mass = gp.mass;
			myParticles[i+1].position = gp.pos;
			myParticles[i+1].velocity = gp.vel;
			myParticles[i+1].soft = gp.hsmooth;
			myParticles[i+1].iType = TYPE_GAS;
		} else if(i + startParticle < (unsigned int) tipsyHeader.nsph
			  + tipsyHeader.ndark) {
			r.getNextDarkParticle(dp);
			myParticles[i+1].mass = dp.mass;
			myParticles[i+1].position = dp.pos;
			myParticles[i+1].velocity = dp.vel;
			myParticles[i+1].soft = dp.eps;
			myParticles[i+1].iType = TYPE_DARK;
		} else {
			r.getNextStarParticle(sp);
			myParticles[i+1].mass = sp.mass;
			myParticles[i+1].position = sp.pos;
			myParticles[i+1].velocity = sp.vel;
			myParticles[i+1].soft = sp.eps;
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
			   int64_t iPrevOffset,
			   const std::string& filename,
			   const double dTime,
			   const double dvFac,
			   const CkCallback& cb)
{
    if(iStage > nSetupWriteStage + 1) {
	// requeue message
	pieces[thisIndex].setupWrite(iStage, iPrevOffset, filename,
					  dTime, dvFac, cb);
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
					     filename, dTime, dvFac, cb);
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
						 filename, dTime, dvFac, cb);
	    }
	if(thisIndex == (int) numTreePieces-1)
	    assert(nStartWrite+myNumParticles == fh.numParticles);
	nSetupWriteStage = -1;	// reset for next time.
	writeTipsy(filename, dTime, dvFac);
	contribute(0, 0, CkReduction::concat, cb);
	}
    }

void TreePiece::writeTipsy(const std::string& filename, const double dTime,
			   const double dvFac)
{
    tipsyHeader.time = dTime;
    
    Tipsy::TipsyWriter w(filename, tipsyHeader);
    
    if(thisIndex == 0)
	w.writeHeader();
    w.seekParticleNum(nStartWrite);
    for(unsigned int i = 0; i < myNumParticles; i++) {
	if(myParticles[i+1].iOrder < (unsigned int) tipsyHeader.nsph) {
	    Tipsy::gas_particle gp;
	    gp.mass = myParticles[i+1].mass;
	    gp.pos = myParticles[i+1].position;
	    gp.vel = myParticles[i+1].velocity*dvFac;
	    gp.hsmooth = myParticles[i+1].soft;
	    gp.phi = myParticles[i+1].potential;

	    w.putNextGasParticle(gp);
	    }
	else if(myParticles[i+1].iOrder < (unsigned int) tipsyHeader.nsph
		      + tipsyHeader.ndark) {
	    Tipsy::dark_particle dp;
	    dp.mass = myParticles[i+1].mass;
	    dp.pos = myParticles[i+1].position;
	    dp.vel = myParticles[i+1].velocity*dvFac;
	    dp.eps = myParticles[i+1].soft;
	    dp.phi = myParticles[i+1].potential;

	    w.putNextDarkParticle(dp);
	    }
	else {
	    Tipsy::star_particle sp;
	    sp.mass = myParticles[i+1].mass;
	    sp.pos = myParticles[i+1].position;
	    sp.vel = myParticles[i+1].velocity*dvFac;
	    sp.eps = myParticles[i+1].soft;
	    sp.phi = myParticles[i+1].potential;

	    w.putNextStarParticle(sp);
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
	int nOutParticles = fh.numParticles/ numTreePieces;
	int excess = fh.numParticles % numTreePieces;
	startParticle[iPiece] = nOutParticles * iPiece;
	if(iPiece < (int) excess) {
	    startParticle[iPiece] += iPiece;
	    }
	else {
	    startParticle[iPiece] += excess;
	    }
	}
    startParticle[numTreePieces] = fh.numParticles; // @TODO: replace
						    // with MaxIOrder
						    // for particle
						    // creation/deletion
    // Sort particles in iOrder
    sort(myParticles+1, myParticles+myNumParticles+1, compIOrder);

    // Tag boundary particle to avoid overruns
    myParticles[myNumParticles+1].iOrder = fh.numParticles;
    
    // Loop through sending particles to correct processor.
    GravityParticle *binBegin = &myParticles[1];
    GravityParticle *binEnd;
    for(iPiece = 0; iPiece < numTreePieces; iPiece++) {
	for(binEnd = binBegin; binEnd->iOrder < startParticle[iPiece+1];
	    binEnd++);
	if((binEnd - binBegin) > 0) {
	    if (verbosity>=3)
		CkPrintf("me:%d to:%d how many:%d\n",thisIndex, iPiece,
			 (binEnd-binBegin));
	    if(iPiece == thisIndex) {
		ioAcceptSortedParticles(binBegin, binEnd - binBegin);
		}
	    else {
		pieces[iPiece].ioAcceptSortedParticles(binBegin,
						       binEnd - binBegin);
		}
	    }
	if(&myParticles[myNumParticles + 1] <= binEnd)
	    break;
	binBegin = binEnd;
	}

    delete[] startParticle;

    // signify completion
    incomingParticlesSelf = true;
    ioAcceptSortedParticles(binBegin, 0);
    }

/// Accept particles from other TreePieces once the sorting has finished
void TreePiece::ioAcceptSortedParticles(const GravityParticle* particles,
					const int n) {

    int myIOParticles = fh.numParticles / numTreePieces;
    
    int excess = fh.numParticles % numTreePieces;
    if(thisIndex < (int) excess) {
	    myIOParticles++;
	    }
  // allocate new particles array on first call
  if (incomingParticles == NULL) {
    incomingParticles = new GravityParticle[myIOParticles + 2];
  }

  memcpy(&incomingParticles[incomingParticlesArrived+1], particles,
	 n*sizeof(GravityParticle));
  incomingParticlesArrived += n;

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

void TreePiece::outputDensityASCII(const string& suffix, const CkCallback& cb) {
  if(thisIndex==0) {
    if(verbosity > 2)
      ckerr << "TreePiece " << thisIndex << ": Writing header for Density file"
	    << endl;
    FILE* outfile = fopen((basefilename + "." + suffix).c_str(), "w");
    fprintf(outfile,"%d\n",(int) fh.numParticles);
    fclose(outfile);
  }
	
  if(verbosity > 3)
    ckerr << "TreePiece " << thisIndex << ": Writing Density to disk" << endl;
	
  FILE* outfile = fopen((basefilename + "." + suffix).c_str(), "r+");
  fseek(outfile, 0, SEEK_END);
  
  for(unsigned int i = 1; i <= myNumParticles; ++i) {
      if(fprintf(outfile,"%g\n", myParticles[i].fDensity) < 0) {
	  ckerr << "TreePiece " << thisIndex
		<< ": Error writing density to disk, aborting" << endl;
	  CkAbort("IO Badness");
	  }
      }
  
  fclose(outfile);
  if(thisIndex==(int)numTreePieces-1) {
      cb.send();
      }
  else
      pieces[thisIndex + 1].outputDensityASCII(suffix, cb);
  }
