/// IO routines
#ifdef HAVE_VALUES_H
#include <values.h>
#else
#include <float.h>
#endif

#include "ParallelGravity.h"
#include "DataManager.h"
#include "TipsyFile.h"
#include "OrientedBox.h"
#include "Reductions.h"
#include "InOutput.h"
#include "ckio.h"
#include <errno.h>
#include <float.h>

using namespace TypeHandling;
using namespace SFC;
using namespace std;

template <typename TPos, typename TVel>
void load_tipsy_gas(Tipsy::TipsyReader &r, GravityParticle &p, double dTuFac) 
{
    Tipsy::gas_particle_t<TPos, TVel> gp;
    
    if(!r.getNextGasParticle_t(gp)) {
        CkAbort("failed to read gas particle!");
        }
    p.mass = gp.mass;
    p.position = gp.pos;
    p.velocity = gp.vel;
    p.soft = gp.hsmooth;
#ifdef CHANGESOFT
    p.fSoft0 = gp.hsmooth;
#endif
    p.iType = TYPE_GAS;
    p.fDensity = gp.rho;
    p.fMetals() = gp.metals;
    // O and Fe ratio based on Asplund et al 2009
    p.fMFracOxygen() = 0.43*gp.metals;
    p.fMFracIron() = 0.098*gp.metals;
#ifdef DIFFUSION
    p.fMetalsPred() = gp.metals;
    p.fMFracOxygenPred() = 0.43*gp.metals;
    p.fMFracIronPred() = 0.098*gp.metals;
#endif
    p.u() = dTuFac*gp.temp;
    p.uPred() = dTuFac*gp.temp;
    // Initial estimate of sound speed.
    double gamma = GAMMA_NONCOOL;
    double gammam1 = gamma - 1.0;
    p.c() = sqrt(gamma*gammam1*p.uPred());
    p.vPred() = gp.vel;
    p.fBallMax() = FLT_MAX;     // N.B. don't use DOUBLE_MAX here:
                                // fBallMax*fBallMax should not overflow.
    p.fESNrate() = 0.0;
    p.fTimeCoolIsOffUntil() = 0.0;
    p.dTimeFB() = 0.0;
#ifdef NEED_DT
    p.dt = FLT_MAX;
#endif
#ifdef DTADJUST
    p.dtNew() = FLT_MAX;
#ifndef COOLING_NONE
    p.uDot() = 0.0;  // Used in initial timestep
#endif
    p.PdV() = 0.0;  // Used in initial timestep
#endif
#ifdef SUPERBUBBLE
    p.cpHotInit() = 0;
    p.uHot() = 0.0;
    p.uHotPred() = 0.0;
    p.uHotDot() = 0.0;
    p.massHot() = 0.0;
    p.fThermalCond() = 0.0;
    p.fThermalLength() = 0.0;
    p.fPromoteSum() = 0.0;
    p.fPromoteSumuPred() = 0.0;
    p.fPromoteuPredInit() = 0.0;
#endif
}

template <typename TPos, typename TVel>
void load_tipsy_dark(Tipsy::TipsyReader &r, GravityParticle &p) 
{
    Tipsy::dark_particle_t<TPos, TVel> dp;
    if(!r.getNextDarkParticle_t(dp)) {
        CkAbort("failed to read dark particle!");
        }
	p.mass = dp.mass;
	p.position = dp.pos;
	p.velocity = dp.vel;
	p.soft = dp.eps;
#ifdef CHANGESOFT
	p.fSoft0 = dp.eps;
#endif
	p.fDensity = 0.0;
	p.iType = TYPE_DARK;
}

template <typename TPos, typename TVel>
void load_tipsy_star(Tipsy::TipsyReader &r, GravityParticle &p) 
{
    Tipsy::star_particle_t<TPos, TVel> sp;

    if(!r.getNextStarParticle_t(sp)) {
        CkAbort("failed to read star particle!");
        }
    p.mass = sp.mass;
    p.position = sp.pos;
    p.velocity = sp.vel;
    p.soft = sp.eps;
#ifdef CHANGESOFT
    p.fSoft0 = sp.eps;
#endif
    p.fDensity = 0.0;
    p.iType = TYPE_STAR;
    p.fStarMetals() = sp.metals;
    // Metals to O and Fe based on Asplund et al 2009
    p.fStarMFracOxygen() = 0.43*sp.metals;
    p.fStarMFracIron() = 0.098*sp.metals;
    p.fMassForm() = sp.mass;
    p.fTimeForm() = sp.tform;
#ifdef COOLING_MOLECULARH 
    p.dStarLymanWerner() = 0.0;
#endif
}


void TreePiece::loadTipsy(const std::string& filename,
			  const double dTuFac, // Convert Temperature
                          const bool bDoublePos,
                          const bool bDoubleVel,
			  const CkCallback& cb) {
        LBTurnInstrumentOff();
        basefilename = filename;

	Tipsy::TipsyReader r(filename, bDoublePos, bDoubleVel);
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

        bool skipLoad = false;
        int numLoadingPEs = CkNumPes();

        // check whether this tree piece should load from file
#ifdef ROUND_ROBIN_WITH_OCT_DECOMP 
        if(thisIndex >= numLoadingPEs) skipLoad = true;
#else
        int numTreePiecesPerPE = numTreePieces/numLoadingPEs;
        int rem = numTreePieces-numTreePiecesPerPE*numLoadingPEs;

#ifdef DEFAULT_ARRAY_MAP
        if (rem > 0) {
          int sizeSmallBlock = numTreePiecesPerPE; 
          int numLargeBlocks = rem; 
          int sizeLargeBlock = numTreePiecesPerPE + 1; 
          int largeBlockBound = numLargeBlocks * sizeLargeBlock; 
          
          if (thisIndex < largeBlockBound) {
            if (thisIndex % sizeLargeBlock > 0) {
              skipLoad = true; 
            }
          }
          else {
            if ((thisIndex - largeBlockBound) % sizeSmallBlock > 0) {
              skipLoad = true; 
            }
          }
        }
        else {
          if ( (thisIndex % numTreePiecesPerPE) > 0) {
            skipLoad = true; 
          }
        }
#else
        // this is not the best way to divide objects among PEs, 
        // but it is how charm++ BlockMap does it.
        if(rem > 0){
          numTreePiecesPerPE++;
          numLoadingPEs = numTreePieces/numTreePiecesPerPE;
          if(numTreePieces % numTreePiecesPerPE > 0) numLoadingPEs++;
        }

        if(thisIndex % numTreePiecesPerPE > 0 || thisIndex >= numLoadingPEs * numTreePiecesPerPE) skipLoad = true;
#endif
#endif

        /*
        if(thisIndex == 0){
          CkPrintf("[%d] numTreePieces %d numTreePiecesPerPE %d numLoadingPEs %d\n", thisIndex, numTreePieces, numTreePiecesPerPE, numLoadingPEs);
        }
        */


        if(skipLoad){
          myNumParticles = 0;
          nStartRead = -1;
          contribute(cb);
          return;
        }

        // find your load offset into input file
        int myIndex = CkMyPe();
	myNumParticles = nTotalParticles / numLoadingPEs;
	int excess = nTotalParticles % numLoadingPEs;
	int64_t startParticle = ((int64_t) myNumParticles) * myIndex;
	if(myIndex < excess) {
	    myNumParticles++;
	    startParticle += myIndex;
	    }
	else {
	    startParticle += excess;
	    }
	if(startParticle >= nTotalParticles) {
	    CkError("Bad startParticle: %ld, nPart: %ld, myIndex: %d, nLoading: %d\n",
		    startParticle, nTotalParticles, myIndex, numLoadingPEs);
	    }
	CkAssert(startParticle < nTotalParticles);
        nStartRead = startParticle;
	
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
	allocateStars();
	
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
                    myParticles[i+1].extraData = &mySPHParticles[iSPH];
                    if(!bDoublePos)
                        load_tipsy_gas<float,float>(r, myParticles[i+1],
                                                    dTuFac) ;
                    else if(!bDoubleVel)
                        load_tipsy_gas<double,float>(r, myParticles[i+1],
                                                     dTuFac) ;
                    else
                        load_tipsy_gas<double,double>(r, myParticles[i+1],
                                                      dTuFac) ;
                    iSPH++;
		} else if(i + startParticle < (unsigned int) tipsyHeader.nsph
			  + tipsyHeader.ndark) {
                    if(!bDoublePos)
                        load_tipsy_dark<float,float>(r, myParticles[i+1]);
                    else if(!bDoubleVel)
                        load_tipsy_dark<double,float>(r, myParticles[i+1]);
                    else
                        load_tipsy_dark<double,double>(r, myParticles[i+1]);
		} else {
                    myParticles[i+1].extraData = &myStarParticles[iStar];
                    if(!bDoublePos)
                        load_tipsy_star<float,float>(r, myParticles[i+1]);
                    else if(!bDoubleVel)
                        load_tipsy_star<double,float>(r, myParticles[i+1]);
                    else
                        load_tipsy_star<double,double>(r, myParticles[i+1]);
                    iStar++;
		}
#ifdef SIDMINTERACT
		myParticles[i+1].iNSIDMInteractions = 0;
#endif
		myParticles[i+1].rung = 0;
		myParticles[i+1].fBall = 0.0;
		myParticles[i+1].iOrder = i + startParticle;
#ifdef SPLITGAS
		if(myParticles[i+1].iOrder < tipsyHeader.nsph) myParticles[i+1].iSplitOrder() = i + startParticle;
		if(myParticles[i+1].iOrder >= tipsyHeader.nsph) myParticles[i+1].iOrder +=  tipsyHeader.nsph;
#endif
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
  myParticles[0].key = firstPossibleKey;
  myParticles[myNumParticles+1].key = lastPossibleKey;
  contribute(cb);
}

/// @brief return maximum iOrders.
/// Like the above, but the file has already been read
void TreePiece::getMaxIOrds(const CkCallback& cb)
{
    CmiInt8 nMaxOrd[3] = {0, 0, 0}; // 0 -> gas, 1 -> dark, 2 -> all

    for(int i = 0; i < myNumParticles; i++) {
        int64_t dummy;
        dummy = myParticles[i+1].iOrder;
        if(dummy > nMaxOrd[0] && myParticles[i+1].isGas())
            nMaxOrd[0] = dummy;
        if(dummy > nMaxOrd[1] && myParticles[i+1].isDark())
            nMaxOrd[1] = dummy;
        if(dummy > nMaxOrd[2])
            nMaxOrd[2] = dummy;
        }
    contribute(3*sizeof(CmiInt8), nMaxOrd, CkReduction::max_long, cb);
    }

void TreePiece::readTipsyArray(OutputParams& params, const CkCallback& cb)
{
    params.dm = dm; // pass cooling information
    FILE *infile = CmiFopen((params.fileName+"." + params.sTipsyExt).c_str(),
                            "r+");
    CkAssert(infile != NULL);
    // Check if its a binary file
    unsigned int iDum;
    XDR xdrs;
    xdrstdio_create(&xdrs, infile, XDR_DECODE);
    xdr_u_int(&xdrs,&iDum);
    xdr_destroy(&xdrs);
    if(iDum == nTotalParticles) { // We've got a binary file; read it
        int64_t seek_pos;
        if(params.bVector)
            seek_pos = sizeof(iDum) + nStartRead*(int64_t)sizeof(float)*3;
        else
            seek_pos = sizeof(iDum) + nStartRead*(int64_t)sizeof(float);
        fseek(infile, seek_pos, SEEK_SET);
        xdrstdio_create(&xdrs, infile, XDR_DECODE);
        for(unsigned int i = 0; i < myNumParticles; ++i) {
            if(params.bFloat) {
                if(params.bVector) {
                    Vector3D<float> vValue;
                    xdr_template(&xdrs, &vValue);
                    params.setVValue(&myParticles[i+1], vValue);
                }
                else {
                    float dValue;
                    xdr_float(&xdrs, &dValue);
                    params.setDValue(&myParticles[i+1], dValue);
                }
            }
            else {
                int iValue;
                xdr_template(&xdrs, &iValue);
                params.setIValue(&myParticles[i+1], iValue);
                }
            }
        xdr_destroy(&xdrs);
    }
    else {                      // assume we've got an ASCII file
        int64_t nTot;
        fseek(infile, 0, SEEK_SET);
        int nread;
        nread = fscanf(infile, "%ld\n", &nTot);
        CkAssert(nread == 1);
        int nDim = 1;           // Dimensions to read
        if(params.bVector) {
            nDim = 3;
            CkAssert(packed == 0);
        }
        for(int iDim = 0; iDim < nDim; iDim++) {
            for(int i = 0; i < nStartRead; i++) {
                double dummy;
                nread = fscanf(infile, "%lf\n", &dummy);
                CkAssert(nread == 1);
            }
            for(int i = 0; i < myNumParticles; i++) {
                if(params.bFloat) {
                    double dDummy;
                    nread = fscanf(infile, "%lf\n", &dDummy);
                    CkAssert(nread == 1);
                    if(params.bVector) {
                        Vector3D<double> vDummy
                            = params.vValue(&myParticles[i+1]);
                        vDummy[iDim] = dDummy;
                        params.setVValue(&myParticles[i+1], vDummy);
                    }
                    else
                        params.setDValue(&myParticles[i+1], dDummy);
                }
                else {
                    int64_t iDummy;
                    nread = fscanf(infile, "%ld\n", &iDummy);
                    CkAssert(nread == 1);
                    params.setIValue(&myParticles[i+1], iDummy);
                }
            }
        }
    }
        
    CmiFclose(infile);
    contribute(cb);
    }

static double fh_time; // gross, but quick way to get time

/// Returns total number of particles in a given file
/// @param filename data file of interest
int64_t ncGetCount(std::string filename)
{
    FILE* infile = CmiFopen(filename.c_str(), "rb");
    if(!infile) {
	return 0;  // Assume there is none of this particle type
	}

    XDR xdrs;
    FieldHeader fh;
    xdrstdio_create(&xdrs, infile, XDR_DECODE);

    if(!xdr_template(&xdrs, &fh)) {
	throw XDRException("Couldn't read header from file!");
	}
    if(fh.magic != FieldHeader::MagicNumber) {
	throw XDRException("This file does not appear to be a field file (magic number doesn't match).");
	}
    if(fh.dimensions != 3 && fh.dimensions != 1) {
	throw XDRException("Wrong dimension.");
	}
    fh_time = fh.time;
    xdr_destroy(&xdrs);
    fclose(infile);
    return fh.numParticles;
    }

static void *readFieldData(const std::string filename, FieldHeader &fh, unsigned int dim,
    int64_t numParticles, int64_t startParticle)
{
    FILE* infile = CmiFopen(filename.c_str(), "rb");
    if(!infile) {
	string smess("Couldn't open field file: ");
	smess += filename;
	throw XDRException(smess);
	}

    XDR xdrs;
    xdrstdio_create(&xdrs, infile, XDR_DECODE);

    if(!xdr_template(&xdrs, &fh)) {
	throw XDRException("Couldn't read header from file!");
	}
    if(fh.magic != FieldHeader::MagicNumber) {
	throw XDRException("This file does not appear to be a field file (magic number doesn't match).");
	}
    if(fh.dimensions != dim) {
	throw XDRException("Wrong dimension of positions.");
	}

    void* data = readField(fh, &xdrs, numParticles, startParticle);
	
    if(data == 0) {
	throw XDRException("Had problems reading in the field");
	}
    xdr_destroy(&xdrs);
    fclose(infile);
    return data;
    }

/// @brief load attributes common to all particles
static void load_NC_base(std::string filename, int64_t startParticle,
                        int myNum, GravityParticle *myParts)
{
    FieldHeader fh;
    // Positions
    if(verbosity && startParticle == 0)
        CkPrintf("loading positions\n");

    void *data = readFieldData(filename + "/pos", fh, 3, myNum,
        startParticle);
    for(int i = 0; i < myNum; ++i) {
	switch(fh.code) {
	case float32:
            myParts[i].position = static_cast<Vector3D<float> *>(data)[i];
	    break;
	case float64:
            myParts[i].position = static_cast<Vector3D<double> *>(data)[i];
	    break;
	default:
	    throw XDRException("I don't recognize the type of this field!");
	    }
        }
    deleteField(fh, data);
    // velocities
    if(verbosity && startParticle == 0)
        CkPrintf("loading velocities\n");
    data = readFieldData(filename + "/vel", fh, 3, myNum,
        startParticle);
    for(int i = 0; i < myNum; ++i) {
	switch(fh.code) {
	case float32:
            myParts[i].velocity = static_cast<Vector3D<float> *>(data)[i];
	    break;
	case float64:
            myParts[i].velocity = static_cast<Vector3D<double> *>(data)[i];
	    break;
	default:
	    throw XDRException("I don't recognize the type of this field!");
	    }
        }
    deleteField(fh, data);
    // masses
    if(verbosity && startParticle == 0)
        CkPrintf("loading masses\n");

    data = readFieldData(filename + "/mass", fh, 1, myNum,
        startParticle);
    for(int i = 0; i < myNum; ++i) {
	switch(fh.code) {
	case float32:
            myParts[i].mass = static_cast<float *>(data)[i];
	    break;
	case float64:
            myParts[i].mass = static_cast<double *>(data)[i];
	    break;
	default:
	    throw XDRException("I don't recognize the type of this field!");
	    }
        }
    deleteField(fh, data);
    // softenings
    if(verbosity && startParticle == 0)
        CkPrintf("loading softenings\n");

    data = readFieldData(filename + "/soft", fh, 1, myNum,
        startParticle);
    for(int i = 0; i < myNum; ++i) {
	switch(fh.code) {
	case float32:
            myParts[i].soft = static_cast<float *>(data)[i];
#ifdef CHANGESOFT
            myParts[i].fSoft0 = static_cast<float *>(data)[i];
#endif
	    break;
	case float64:
            myParts[i].soft = static_cast<double *>(data)[i];
#ifdef CHANGESOFT
            myParts[i].fSoft0 = static_cast<double *>(data)[i];
#endif
	    break;
	default:
	    throw XDRException("I don't recognize the type of this field!");
	    }
        }
    deleteField(fh, data);
}

static void load_NC_gas(std::string filename, int64_t startParticle,
                        int myNumSPH, GravityParticle *myParts,
                        extraSPHData *mySPHParts, double dTuFac)
{
    if(verbosity && startParticle == 0)
        CkPrintf("loading gas\n");

    load_NC_base(filename, startParticle, myNumSPH, myParts);

    for(int i = 0; i < myNumSPH; ++i) {
        myParts[i].iType = TYPE_GAS;
        myParts[i].extraData = &mySPHParts[i];
        myParts[i].fBallMax() = FLT_MAX;  // N.B. don't use DOUBLE_MAX here:
                                          // fBallMax*fBallMax should not overflow.
        myParts[i].fESNrate() = 0.0;
        myParts[i].fTimeCoolIsOffUntil() = 0.0;
        myParts[i].dTimeFB() = 0.0;
#ifdef NEED_DT
        myParts[i].dt = FLT_MAX;
#endif
#ifdef DTADJUST
        myParts[i].dtNew() = FLT_MAX;
#ifndef COOLING_NONE
        myParts[i].uDot() = 0.0;  // Used in initial timestep
#endif
        myParts[i].PdV() = 0.0;  // Used in initial timestep
#endif
#ifdef SUPERBUBBLE
        myParts[i].cpHotInit() = 0;
        myParts[i].uHot() = 0.0;
        myParts[i].uHotPred() = 0.0;
        myParts[i].uHotDot() = 0.0;
        myParts[i].massHot() = 0.0;
        myParts[i].fThermalCond() = 0.0;
        myParts[i].fThermalLength() = 0.0;
        myParts[i].fPromoteSum() = 0.0;
        myParts[i].fPromoteSumuPred() = 0.0;
        myParts[i].fPromoteuPredInit() = 0.0;
#endif
        }

    FieldHeader fh;
    // Density
    if(verbosity && startParticle == 0)
        CkPrintf("loading densities\n");

    void *data = readFieldData(filename + "/GasDensity", fh, 1, myNumSPH,
                               startParticle);
    for(int i = 0; i < myNumSPH; ++i) {
	switch(fh.code) {
	case float32:
            myParts[i].fDensity = static_cast<float *>(data)[i];
	    break;
	case float64:
            myParts[i].fDensity = static_cast<double *>(data)[i];
	    break;
	default:
	    throw XDRException("I don't recognize the type of this field!");
	    }
        }
    deleteField(fh, data);
  if(ncGetCount(filename + "/OxMassFrac") > 0) {
    // Oxygen
    if(verbosity && startParticle == 0)
        CkPrintf("loading Oxygen\n");

    data = readFieldData(filename + "/OxMassFrac", fh, 1, myNumSPH,
                               startParticle);
    for(int i = 0; i < myNumSPH; ++i) {
	switch(fh.code) {
	case float32:
            myParts[i].fMFracOxygen() = static_cast<float *>(data)[i];
	    break;
	case float64:
            myParts[i].fMFracOxygen() = static_cast<double *>(data)[i];
	    break;
	default:
	    throw XDRException("I don't recognize the type of this field!");
	    }
        }
    deleteField(fh, data);
  }
  if(ncGetCount(filename + "/FeMassFrac") > 0) {
    // Iron
    if(verbosity && startParticle == 0)
        CkPrintf("loading Iron\n");

    data = readFieldData(filename + "/FeMassFrac", fh, 1, myNumSPH,
                               startParticle);
    for(int i = 0; i < myNumSPH; ++i) {
	switch(fh.code) {
	case float32:
            myParts[i].fMFracIron() = static_cast<float *>(data)[i];
	    break;
	case float64:
            myParts[i].fMFracIron() = static_cast<double *>(data)[i];
	    break;
	default:
	    throw XDRException("I don't recognize the type of this field!");
	    }
        }
    deleteField(fh, data);
  }
    // Temperature
    if(verbosity && startParticle == 0)
        CkPrintf("loading temperature\n");

    data = readFieldData(filename + "/temperature", fh, 1, myNumSPH,
                               startParticle);
    for(int i = 0; i < myNumSPH; ++i) {
	switch(fh.code) {
	case float32:
            myParts[i].u() = dTuFac*static_cast<float *>(data)[i];
	    break;
	case float64:
            myParts[i].u() = dTuFac*static_cast<double *>(data)[i];
	    break;
	default:
	    throw XDRException("I don't recognize the type of this field!");
	    }
        myParts[i].uPred() = myParts[i].u();
        // Initial estimate of sound speed.
        double gamma = GAMMA_NONCOOL;
        double gammam1 = gamma - 1.0;
        myParts[i].c() = sqrt(gamma*gammam1*myParts[i].uPred());
        }
    deleteField(fh, data);
}

static void load_NC_dark(std::string filename, int64_t startParticle,
                        int myNumDark, GravityParticle *myParts)
{
    if(verbosity && startParticle == 0)
        CkPrintf("loading darks\n");
    
    load_NC_base(filename, startParticle, myNumDark, myParts);

    for(int i = 0; i < myNumDark; ++i) {
        myParts[i].fDensity = 0.0;
        myParts[i].iType = TYPE_DARK;
        }
}

static void load_NC_star(std::string filename, int64_t startParticle,
                         int myNumStar, GravityParticle *myParts,
                         extraStarData *myStarParts)
{
    if(verbosity && startParticle == 0)
        CkPrintf("loading stars\n");

    for(int i = 0; i < myNumStar; ++i) {
        myParts[i].fDensity = 0.0;
        myParts[i].iType = TYPE_STAR;
        myParts[i].extraData = &myStarParts[i];
        }
    load_NC_base(filename, startParticle, myNumStar, myParts);

    FieldHeader fh;
    // Oxygen
    if(verbosity && startParticle == 0)
        CkPrintf("loading Oxygen\n");

    void *data = readFieldData(filename + "/OxMassFrac", fh, 1, myNumStar,
                               startParticle);
    for(int i = 0; i < myNumStar; ++i) {
	switch(fh.code) {
	case float32:
            myParts[i].fStarMFracOxygen() = static_cast<float *>(data)[i];
	    break;
	case float64:
            myParts[i].fStarMFracOxygen() = static_cast<double *>(data)[i];
	    break;
	default:
	    throw XDRException("I don't recognize the type of this field!");
	    }
        }
    deleteField(fh, data);
    // Iron
    if(verbosity && startParticle == 0)
        CkPrintf("loading Iron\n");

    data = readFieldData(filename + "/FeMassFrac", fh, 1, myNumStar,
                               startParticle);
    for(int i = 0; i < myNumStar; ++i) {
	switch(fh.code) {
	case float32:
            myParts[i].fStarMFracIron() = static_cast<float *>(data)[i];
	    break;
	case float64:
            myParts[i].fStarMFracIron() = static_cast<double *>(data)[i];
	    break;
	default:
	    throw XDRException("I don't recognize the type of this field!");
	    }
        }
    deleteField(fh, data);
    // Formation Time
    if(verbosity && startParticle == 0)
        CkPrintf("loading timeform\n");

    data = readFieldData(filename + "/timeform", fh, 1, myNumStar,
                               startParticle);
    for(int i = 0; i < myNumStar; ++i) {
	switch(fh.code) {
	case float32:
            myParts[i].fTimeForm() = static_cast<float *>(data)[i];
	    break;
	case float64:
            myParts[i].fTimeForm() = static_cast<double *>(data)[i];
	    break;
	default:
	    throw XDRException("I don't recognize the type of this field!");
	    }
        }
    deleteField(fh, data);
    for(int i = 0; i < myNumStar; ++i) {
        myParts[i].fMassForm() = myParts[i].mass;
#ifdef COOLING_MOLECULARH 
        myParts[i].dStarLymanWerner() = 0.0;
#endif
        }
}

void TreePiece::loadNChilada(const std::string& filename,
                             const double dTuFac, // Convert Temperature
			     const CkCallback& cb) {
        LBTurnInstrumentOff();
        basefilename = filename;

	nTotalSPH = ncGetCount(filename + "/gas/pos");
	nTotalDark = ncGetCount(filename + "/dark/pos");
	nTotalStar = ncGetCount(filename + "/star/pos");
	nTotalParticles = nTotalSPH + nTotalDark + nTotalStar;
        if(nTotalParticles <= 0)
            CkAbort("No particles can be read.  Check file permissions\n");
	dStartTime = fh_time;

        bool skipLoad = false;
        int numLoadingPEs = CkNumPes();

        // check whether this tree piece should load from file
#ifdef ROUND_ROBIN_WITH_OCT_DECOMP 
        if(thisIndex >= numLoadingPEs) skipLoad = true;
#else
        int numTreePiecesPerPE = numTreePieces/numLoadingPEs;
        int rem = numTreePieces-numTreePiecesPerPE*numLoadingPEs;

#ifdef DEFAULT_ARRAY_MAP
        if (rem > 0) {
          int sizeSmallBlock = numTreePiecesPerPE; 
          int numLargeBlocks = rem; 
          int sizeLargeBlock = numTreePiecesPerPE + 1; 
          int largeBlockBound = numLargeBlocks * sizeLargeBlock; 
          
          if (thisIndex < largeBlockBound) {
            if (thisIndex % sizeLargeBlock > 0) {
              skipLoad = true; 
            }
          }
          else {
            if ((thisIndex - largeBlockBound) % sizeSmallBlock > 0) {
              skipLoad = true; 
            }
          }
        }
        else {
          if ( (thisIndex % numTreePiecesPerPE) > 0) {
            skipLoad = true; 
          }
        }
#else
        // this is not the best way to divide objects among PEs, 
        // but it is how charm++ BlockMap does it.
        if(rem > 0){
          numTreePiecesPerPE++;
          numLoadingPEs = numTreePieces/numTreePiecesPerPE;
          if(numTreePieces % numTreePiecesPerPE > 0) numLoadingPEs++;
        }

        if(thisIndex % numTreePiecesPerPE > 0 || thisIndex >= numLoadingPEs * numTreePiecesPerPE) skipLoad = true;
#endif
#endif

        if(skipLoad){
          myNumParticles = 0;
          nStartRead = -1;
          contribute(cb);
          return;
        }

        // find your load offset into input file
        int myIndex = CkMyPe();
	myNumParticles = nTotalParticles / numLoadingPEs;
	int excess = nTotalParticles % numLoadingPEs;
	int64_t startParticle = ((int64_t)myNumParticles) * myIndex;
	if(myIndex < excess) {
	    myNumParticles++;
	    startParticle += myIndex;
	    }
	else {
	    startParticle += excess;
	    }
	if(startParticle >= nTotalParticles) {
	    CkError("Bad startParticle: %ld, nPart: %ld, myIndex: %d, nLoading: %d\n",
		    startParticle, nTotalParticles, myIndex, numLoadingPEs);
	    }
	CkAssert(startParticle < nTotalParticles);
        nStartRead = startParticle;
	
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
	allocateStars();
	
        if(myNumSPH > 0) {
            load_NC_gas(filename + "/gas", startParticle, myNumSPH,
                        &myParticles[1], mySPHParticles, dTuFac);
            }
        int myNumDark = myNumParticles - myNumSPH - myNumStar;
        startParticle -= nTotalSPH;
        if(startParticle < 0)
            startParticle = 0;
        if(myNumDark > 0) {
            load_NC_dark(filename + "/dark", startParticle, myNumDark,
                         &myParticles[myNumSPH + 1]);
            }
        startParticle = nStartRead - nTotalSPH - nTotalDark;
        if(startParticle < 0)
            startParticle = 0;
        if(myNumStar > 0) {
            load_NC_star(filename + "/star", startParticle, myNumStar,
                         &myParticles[myNumSPH + myNumDark + 1],
                         myStarParticles);
            }
        for(int i = 0; i < myNumParticles; ++i) {
            myParticles[i+1].rung = 0;
            myParticles[i+1].fBall = 0.0;
            myParticles[i+1].iOrder = i + nStartRead;
#ifdef SIDMINTERACT
            myParticles[i+1].iNSIDMInteractions = 0;
#endif
#ifdef SPLITGAS
            if(myParticles[i+1].iOrder < nTotalSPH) myParticles[i+1].iSplitOrder() = i + nStartRead;
            if(myParticles[i+1].iOrder >= nTotalSPH) myParticles[i+1].iOrder  += nTotalSPH;
#endif
            boundingBox.grow(myParticles[i+1].position);
        }
        
  myParticles[0].key = firstPossibleKey;
  myParticles[myNumParticles+1].key = lastPossibleKey;
  contribute(cb);
}

/// Generic read of binary (NChilada) array format into floating point
/// or integer particle attribute  (NOTE MISNOMER)
void TreePiece::readFloatBinary(OutputParams& params, int bParaRead,
    const CkCallback& cb)
{
    FieldHeader fh;
    void *data;
    int64_t startParticle = nStartRead;
    params.dm = dm; // pass cooling information
    
    if((params.iType & TYPE_GAS) && (myNumSPH > 0)) {
        data = readFieldData(params.fileName + "/gas/" + params.sNChilExt, fh,
            1, myNumSPH, nStartRead);
        for(int i = 0; i < myNumSPH; ++i) {
            switch(fh.code) {
            case int32:
                params.setIValue(&myParticles[i+1], static_cast<int *>(data)[i]);
                break;
            case int64:
                params.setIValue(&myParticles[i+1], static_cast<int64_t *>(data)[i]);
                break;
            case float32:
                params.setDValue(&myParticles[i+1], static_cast<float *>(data)[i]);
                break;
            case float64:
                params.setDValue(&myParticles[i+1], static_cast<double *>(data)[i]);
                break;
            default:
                throw XDRException("I don't recognize the type of this field!");
                }
            }
        deleteField(fh, data);
        }
    int myNumDark = myNumParticles - myNumSPH - myNumStar;
    startParticle -= nTotalSPH;
    if(startParticle < 0) startParticle = 0;
    if((params.iType & TYPE_DARK) && (myNumDark > 0)) {
        data = readFieldData(params.fileName + "/dark/" + params.sNChilExt, fh,
            1, myNumDark, startParticle);
        for(int i = 0; i < myNumDark; ++i) {
            switch(fh.code) {
            case int32:
                params.setIValue(&myParticles[myNumSPH+i+1], static_cast<int *>(data)[i]);
                break;
            case int64:
                params.setIValue(&myParticles[myNumSPH+i+1], static_cast<int64_t *>(data)[i]);
                break;
            case float32:
                params.setDValue(&myParticles[myNumSPH+i+1], static_cast<float *>(data)[i]);
                break;
            case float64:
                params.setDValue(&myParticles[myNumSPH+i+1], static_cast<double *>(data)[i]);
                break;
            default:
                throw XDRException("I don't recognize the type of this field!");
                }
            }
        deleteField(fh, data);
        }
    startParticle = nStartRead - nTotalSPH - nTotalDark;
    if(startParticle < 0)
        startParticle = 0;
    if((params.iType & TYPE_STAR) && (myNumStar > 0)) {
        data = readFieldData(params.fileName + "/star/" + params.sNChilExt, fh,
            1, myNumStar, startParticle);
        for(int i = 0; i < myNumStar; ++i) {
            switch(fh.code) {
            case int32:
                params.setIValue(&myParticles[myNumSPH+myNumDark+i+1], static_cast<int *>(data)[i]);
                break;
            case int64:
                params.setIValue(&myParticles[myNumSPH+myNumDark+i+1], static_cast<int64_t *>(data)[i]);
                break;
            case float32:
                params.setDValue(&myParticles[myNumSPH+myNumDark+i+1], static_cast<float *>(data)[i]);
                break;
            case float64:
                params.setDValue(&myParticles[myNumSPH+myNumDark+i+1], static_cast<double *>(data)[i]);
                break;
            default:
                throw XDRException("I don't recognize the type of this field!");
                }
            }
        deleteField(fh, data);
        }
    contribute(cb);
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
                           const bool bDoublePos,
                           const bool bDoubleVel,
			   const int bCool,
			   const CkCallback& cb)
{
    if(iStage > nSetupWriteStage + 1) {
	// requeue message
	pieces[thisIndex].setupWrite(iStage, iPrevOffset, filename,
				     dTime, dvFac, duTFac, bDoublePos,
                                     bDoubleVel, bCool, cb);
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
					     bDoublePos, bDoubleVel, bCool, cb);
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
						 duTFac, bDoublePos,
                                                 bDoubleVel, bCool, cb);
	    }
	if(thisIndex == (int) numTreePieces-1)
	    assert(nStartWrite+myNumParticles == nTotalParticles);
	nSetupWriteStage = -1;	// reset for next time.
	parallelWrite(0, cb, filename, dTime, dvFac, duTFac, bDoublePos,
                      bDoubleVel, bCool);
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
                              const bool bDoublePos,
                              const bool bDoubleVel,
			      const int bCool)
{
    Tipsy::header tipsyHeader;

    tipsyHeader.time = dTime;
    tipsyHeader.nbodies = nTotalParticles;
    tipsyHeader.nsph = nTotalSPH;
    tipsyHeader.nstar = nTotalStar;
    tipsyHeader.ndark = nTotalParticles - (nTotalSPH + nTotalStar);

    if(nIOProcessor == 0) {	// use them all
	Tipsy::TipsyWriter wAll(filename, tipsyHeader, false, bDoublePos,
                                bDoubleVel);
	writeTipsy(wAll, dvFac, duTFac, bDoublePos, bDoubleVel, bCool);
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

    Tipsy::TipsyWriter w(filename, tipsyHeader, false, bDoublePos, bDoubleVel);
    writeTipsy(w, dvFac, duTFac, bDoublePos, bDoubleVel, bCool);
    if(iPass < (nSkip - 1) && thisIndex < (numTreePieces - 1))
	treeProxy[thisIndex+1].parallelWrite(iPass + 1, cb, filename, dTime,
					     dvFac, duTFac, bDoublePos,
                                             bDoubleVel, bCool);
    contribute(cb);
    }


/// @brief Serial output of tipsy file.
///
/// Send my particles to piece 0 for output.
///
void TreePiece::serialWrite(const u_int64_t iPrevOffset, // previously written
							 //particles
			    const std::string& filename,  // output file
			    const double dTime,	 // time or expansion
			    const double dvFac,  // velocity conversion
			    const double duTFac, // temperature conversion
                            const bool bDoublePos,
                            const bool bDoubleVel,
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
			   dvFac, duTFac, bDoublePos, bDoubleVel, bCool, cb);
    if(myNumSPH > 0)
	delete [] piSph;
    if(myNumStar > 0)
	delete [] piStar;
    }

/// @brief a global variable to store the tipsy writing state for
/// serial I/O.
///
/// Since this is for serial I/O there are no issues for multiple
/// cores accessing this.
///
Tipsy::TipsyWriter *globalTipsyWriter;

/// @brief write out the particles I have been sent
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
                        const bool bDoublePos,
                        const bool bDoubleVel,
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
    if(iIndex == 0) {
        // Create and truncate output file.    
        FILE *fp = CmiFopen(filename.c_str(), "w");
        CmiFclose(fp);
        
	Tipsy::header tipsyHeader;

	tipsyHeader.time = dTime;
	tipsyHeader.nbodies = nTotalParticles;
	tipsyHeader.nsph = nTotalSPH;
	tipsyHeader.nstar = nTotalStar;
	tipsyHeader.ndark = nTotalParticles - (nTotalSPH + nTotalStar);
    
	globalTipsyWriter = new Tipsy::TipsyWriter(filename, tipsyHeader,
                                                   false, bDoublePos,
                                                   bDoubleVel);
	}
	    
    writeTipsy(*globalTipsyWriter, dvFac, duTFac, bDoublePos, bDoubleVel,
               bCool);
    /*
     * Restore pointers/data
     */
    myNumParticles = saveNumParticles;
    myParticles = saveMyParticles;
    mySPHParticles = savemySPHParticles;
    myStarParticles = savemyStarParticles;
    
    if(iIndex < (numTreePieces - 1))
	treeProxy[iIndex+1].serialWrite(iPrevOffset + iOutParticles, filename,
					dTime, dvFac, duTFac, 
                                        bDoublePos, bDoubleVel, bCool, cb);
    else {
	delete globalTipsyWriter;
	cb.send();  // we are done.
	}
    }
    
template <typename TPos, typename TVel>
void write_tipsy_gas(Tipsy::TipsyWriter &w, GravityParticle &p,
                     const double dvFac,
                     const double duTFac,
                     const int bCool,
                     COOL *Cool) 
{
    Tipsy::gas_particle_t<TPos, TVel> gp;

    gp.mass = p.mass;
    gp.pos = p.position;
    gp.vel = p.velocity*dvFac;
#ifdef CHANGESOFT
    gp.hsmooth = p.fSoft0;
#else
    gp.hsmooth = p.soft;
#endif
    gp.phi = p.potential;
    gp.rho = p.fDensity;
    gp.metals = p.fMetals();
    if(bCool) {
#ifndef COOLING_NONE
#ifdef COOLING_GRACKLE
        gp.temp = CoolCodeEnergyToTemperature(Cool, &p.CoolParticle(), p.u(),
                                              p.fDensity, p.fMetals());
#else
        gp.temp = CoolCodeEnergyToTemperature(Cool, &p.CoolParticle(), p.u(),
                                              p.fMetals());
#endif
#else
        CkAbort("cooling output without cooling code");
#endif
        }
    else 
        gp.temp = duTFac*p.u();

    if(!w.putNextGasParticle_t(gp)) {
        CkError("[%d] Write gas failed, errno %d\n", CkMyPe(), errno);
        CkAbort("Bad Write");
        }
}

template <typename TPos, typename TVel>
void write_tipsy_dark(Tipsy::TipsyWriter &w, GravityParticle &p,
                     const double dvFac) 
{
    Tipsy::dark_particle_t<TPos, TVel> dp;

    dp.mass = p.mass;
    dp.pos = p.position;
    dp.vel = p.velocity*dvFac;
#ifdef CHANGESOFT
    dp.eps = p.fSoft0;
#else
    dp.eps = p.soft;
#endif
    dp.phi = p.potential;

    if(!w.putNextDarkParticle_t(dp)) {
        CkError("[%d] Write dark failed, errno %d\n", CkMyPe(), errno);
        CkAbort("Bad Write");
    }
}

template <typename TPos, typename TVel>
void write_tipsy_star(Tipsy::TipsyWriter &w, GravityParticle &p,
                     const double dvFac) 
{
    Tipsy::star_particle_t<TPos, TVel> sp;

    sp.mass = p.mass;
    sp.pos = p.position;
    sp.vel = p.velocity*dvFac;
#ifdef CHANGESOFT
    sp.eps = p.fSoft0;
#else
    sp.eps = p.soft;
#endif
    sp.phi = p.potential;
    sp.metals = p.fStarMetals();
    sp.tform = p.fTimeForm();

    if(!w.putNextStarParticle_t(sp)) {
        CkError("[%d] Write star failed, errno %d\n", CkMyPe(), errno);
        CkAbort("Bad Write");
    }
}

void TreePiece::writeTipsy(Tipsy::TipsyWriter& w,
			   const double dvFac, // scale velocities
			   const double duTFac, // convert temperature
                           const bool bDoublePos,
                           const bool bDoubleVel,
			   const int bCool)
{
    
    if(nStartWrite == 0)
	w.writeHeader();
    if(!w.seekParticleNum(nStartWrite))
	CkAbort("bad seek");
    for(unsigned int i = 0; i < myNumParticles; i++) {
	if(myParticles[i+1].isGas()) {
            if(!bDoublePos)
                write_tipsy_gas<float,float>(w, myParticles[i+1], dvFac,
                                             duTFac, bCool, dm->Cool);
            else if(!bDoubleVel)
                write_tipsy_gas<double,float>(w, myParticles[i+1], dvFac,
                                              duTFac, bCool, dm->Cool);
            else
                write_tipsy_gas<double,double>(w, myParticles[i+1], dvFac,
                                               duTFac, bCool, dm->Cool);
                
	    }
	else if(myParticles[i+1].isDark()) {
            if(!bDoublePos)
                write_tipsy_dark<float,float>(w, myParticles[i+1], dvFac);
            else if(!bDoubleVel)
                write_tipsy_dark<double,float>(w, myParticles[i+1], dvFac);
            else
                write_tipsy_dark<double,double>(w, myParticles[i+1], dvFac);
	    }
	else if(myParticles[i+1].isStar()) {
            if(!bDoublePos)
                write_tipsy_star<float,float>(w, myParticles[i+1], dvFac);
            else if(!bDoubleVel)
                write_tipsy_star<double,float>(w, myParticles[i+1], dvFac);
            else
                write_tipsy_star<double,double>(w, myParticles[i+1], dvFac);
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
	int64_t nOutParticles = nMaxOrder/ nPieces;
	startParticle[iPiece] = nOutParticles * iPiece;
	}
    startParticle[nPieces] = nMaxOrder+1;
    return startParticle;
    }

///
/// @brief Reorder particles for output
///
void TreePiece::reOrder(int64_t _nMaxOrder, const CkCallback& cb)
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

  double tpLoad = getObjTime();
  populateSavedPhaseData(prevLARung, tpLoad, treePieceActivePartsTmp);

    if (myNumParticles > 0) {
	// Particles have been sorted in reOrder()
    // Loop through sending particles to correct processor.
    GravityParticle *binBegin = &myParticles[1];
    GravityParticle *binEnd;
    for(iPiece = 0; iPiece < numTreePieces; iPiece++) {
	for(binEnd = binBegin; binEnd->iOrder < startParticle[iPiece+1];
	    binEnd++);
	int nPartOut = binEnd - binBegin;
  int saved_phase_len = savedPhaseLoad.size();
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
		= new (saved_phase_len, saved_phase_len, nPartOut, nGasOut, nStarOut)
		ParticleShuffleMsg(saved_phase_len, nPartOut, nGasOut, nStarOut, 0.0);
    memset(shuffleMsg->parts_per_phase, 0, saved_phase_len*sizeof(unsigned int));
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

        for(int i = 0; i < saved_phase_len; i++) {
          if (pPart->rung >= i) {
            shuffleMsg->parts_per_phase[i] = shuffleMsg->parts_per_phase[i] + 1;
          }
        }
        if(havePhaseData(PHASE_FEEDBACK)
           && (pPart->isGas() || pPart->isStar()))
            shuffleMsg->parts_per_phase[PHASE_FEEDBACK] += 1;
		}
      shuffleMsg->load = tpLoad * nPartOut / myNumParticles;
      memset(shuffleMsg->loads, 0.0, saved_phase_len*sizeof(double));

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
		CkPrintf("me:%d to:%d how many:%ld\n",thisIndex, iPiece,
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
  treePieceLoadTmp += shuffleMsg->load;
  savePhaseData(savedPhaseLoadTmp, savedPhaseParticleTmp, shuffleMsg->loads,
      shuffleMsg->parts_per_phase, shuffleMsg->nloads);
	}

    if(verbosity > 2)
	ckout << thisIndex << ": incoming: " << incomingParticlesArrived
	      << " myIO: " << myIOParticles << endl;
    
  if(myIOParticles == incomingParticlesArrived && incomingParticlesSelf) {
      //I've got all my particles, now count them
    if(verbosity>1) ckout << thisIndex <<" got ioParticles"
			  <<endl;

    treePieceLoad = treePieceLoadTmp;
    treePieceLoadTmp = 0.0;

    savedPhaseLoad.swap(savedPhaseLoadTmp);
    savedPhaseParticle.swap(savedPhaseParticleTmp);
    savedPhaseLoadTmp.clear();
    savedPhaseParticleTmp.clear();

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
    if(nStoreSPH > 0)
        mySPHParticles = new extraSPHData[nStoreSPH];
    else
        mySPHParticles = NULL;

    myNumStar = nStar;
    if(nStoreStar > 0) delete[] myStarParticles;
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

    deleteTree();
    contribute(callback);
  }
}

// Output a Tipsy ASCII array file.
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
  int *aiOut;	// array for oneNode I/O
  double *adOut;	// array for oneNode I/O
  Vector3D<double> *avOut;	// array for one node I/O
  params.dm = dm; // pass cooling information
  
  if((thisIndex==0 && packed) || (thisIndex==0 && !packed && cnt==0)) {
    if(verbosity > 2)
      ckout << "TreePiece " << thisIndex << ": Writing header for output file" << endl;
    outfile = CmiFopen((params.fileName+"."+params.sTipsyExt).c_str(), "w");
    CkAssert(outfile != NULL);
    fprintf(outfile,"%d\n",(int) nTotalParticles);
    CmiFclose(outfile);
  }
	
  if(verbosity > 3)
    ckout << "TreePiece " << thisIndex << ": Writing output to disk" << endl;
	
  if(bParaWrite) {
      outfile = CmiFopen((params.fileName+"."+params.sTipsyExt).c_str(), "a");
      if(outfile == NULL)
	    ckerr << "Treepiece " << thisIndex << " failed to open "
		  << params.fileName.c_str() << " : " << errno << endl;
      CkAssert(outfile != NULL);
      }
  else {
      if(params.bFloat) {
          if(params.bVector && packed)
              avOut = new Vector3D<double>[myNumParticles];
          else
              adOut = new double[myNumParticles];
          }
      else
          aiOut = new int[myNumParticles];
      }
    for(unsigned int i = 1; i <= myNumParticles; ++i) {
      Vector3D<double> vOut;
      double dOut;
      int iOut;
      if(params.bFloat) {
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
          }
      else
          iOut = params.iValue(&myParticles[i]);
      
      if(bParaWrite) {
        if(params.bFloat) {
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
            if(fprintf(outfile,"%d\n", iOut) < 0) {
                ckerr << "TreePiece " << thisIndex << ": Error writing array to disk, aborting" << endl;
                CkAbort("Badness");
                }
            }
	  }
      else {
        if(params.bFloat) {
	  if(params.bVector && packed)
	      avOut[i-1] = vOut;
	  else
	      adOut[i-1] = dOut;
	  }
        else
            aiOut[i-1] = iOut;
        }
      }
  cnt++;
  if(cnt == 3 || !params.bVector)
      cnt = 0;
  
  if(bParaWrite) {
      int result = CmiFclose(outfile);
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
      if(params.bFloat) {
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
      else {
          pieces[0].oneNodeOutIntArr(params, aiOut, myNumParticles, thisIndex, cb);
          delete [] aiOut;
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
			      const CkCallback& cb)
{
    FILE* outfile = CmiFopen((params.fileName+"."+params.sTipsyExt).c_str(), "a");
    if(outfile == NULL)
	ckerr << "Treepiece " << thisIndex << " failed to open "
	      << params.fileName.c_str() << " : " << errno << endl;
    CkAssert(outfile != NULL);
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
    int result = CmiFclose(outfile);
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
			      const CkCallback& cb)
{
    FILE* outfile = CmiFopen((params.fileName+"."+params.sTipsyExt).c_str(), "a");
    if(outfile == NULL)
	ckerr << "Treepiece " << thisIndex << " failed to open "
	      << params.fileName.c_str() << " : " << errno << endl;
    CkAssert(outfile != NULL);
    for(int i = 0; i < nPart; ++i) {
	if(fprintf(outfile,"%.14g\n",adOut[i]) < 0) {
	  ckerr << "TreePiece " << thisIndex << ": Error writing array to disk, aborting" << endl;
	    CkAbort("Badness");
	    }
	}
    int result = CmiFclose(outfile);
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

// Receives an array of ints to write out in ASCII format
// Assumed to be called from outputASCII() and will continue with the
// next tree piece.

void TreePiece::oneNodeOutIntArr(OutputParams& params,
			      int *aiOut, // array to be output
			      int nPart, // length of adOut
			      int iIndex, // treepiece which called me
			      const CkCallback& cb)
{
    FILE* outfile = CmiFopen((params.fileName+"."+params.sTipsyExt).c_str(), "a");
    if(outfile == NULL)
	ckerr << "Treepiece " << thisIndex << " failed to open "
	      << params.fileName.c_str() << " : " << errno << endl;
    CkAssert(outfile != NULL);
    for(int i = 0; i < nPart; ++i) {
	if(fprintf(outfile,"%d\n",aiOut[i]) < 0) {
	  ckerr << "TreePiece " << thisIndex << ": Error writing array to disk, aborting" << endl;
	    CkAbort("Badness");
	    }
	}
    int result = CmiFclose(outfile);
    if(result != 0)
	ckerr << "Bad close: " << strerror(errno) << endl;
    CkAssert(result == 0);
    if(iIndex!=(int)numTreePieces-1) {
	  pieces[iIndex + 1].outputASCII(params, 0, cb);
	  return;
	  }

    cb.send(); // We are done.
    }

/// Output a Tipsy or NChilada XDR binary float array file.
void Main::outputBinary(OutputParams& params, // specifies
                                              // filename, format,
                                              // and quantity to
                                              // be output
                        int bParaWrite,	      // Every processor can write.
                        const CkCallback& cb) {

    // Save params and callback
    pOutput = &params;
    cbIO = cb;
    Ck::IO::Options opts;
        opts.basePE = 0;
    if(bParaWrite) {
        if(param.nIOProcessor == 0) {
            opts.activePEs = CkNumNodes();
            }
        else {
            opts.activePEs = param.nIOProcessor;
            }
        }
    else {
        opts.activePEs = 1;
        }
    
    if(params.iBinaryOut != 6) {
        Ck::IO::open(params.fileName+"."+params.sTipsyExt,
                     CkCallback(CkIndex_Main::cbOpen(0), thishandle),
                     opts);
        }
    else {
        params.iTypeWriting = 0;
        std::string strFile = getNCNextOutput(params);
        if(strFile.empty()) { // No data to write
            cb.send();
            return;
            }
        Ck::IO::open(strFile,
                     CkCallback(CkIndex_Main::cbOpen(0), thishandle),
                     opts);
        }
}

/// Return the file name of the next nchilada output to write (gas,
/// dark or star).  Returns empty if we are done.  Also advances the
/// iTypeWriting in params.
std::string Main::getNCNextOutput(OutputParams& params) 
{
    if(params.iTypeWriting == 0) {
        params.iTypeWriting = TYPE_GAS;
        if((params.iType & TYPE_GAS) && (nTotalSPH > 0)) {
            return params.fileName+"/gas/"+params.sNChilExt;
            }
        }
    if(params.iTypeWriting == TYPE_GAS) {
        params.iTypeWriting = TYPE_DARK;
        if((params.iType & TYPE_DARK) && (nTotalDark > 0)) {
            return params.fileName+"/dark/"+params.sNChilExt;
            }
        }
    if(params.iTypeWriting == TYPE_DARK) {
        params.iTypeWriting = TYPE_STAR;
        if((params.iType & TYPE_STAR) && (nTotalStar > 0)) {
            return params.fileName+"/star/"+params.sNChilExt;
            }
        }
    return "";
}

/// Determine offsets for all pieces, then start the write session.
/// This needs to be a threaded entry method.
void Main::cbOpen(Ck::IO::FileReadyMsg *msg)
{
    size_t nHeader;
    size_t nBytes;
    int64_t nParts;
    if(pOutput->iBinaryOut != 6) {
        nHeader = sizeof(int);
        nParts = nTotalParticles;
    }
    else {
        nHeader = FieldHeader::sizeBytes;
        if(pOutput->bFloat) {
            if(pOutput->bVector) 
                nHeader += 6*sizeof(float);
            else
                nHeader += 2*sizeof(float);
        }
        else
            nHeader += 2*sizeof(int64_t);
        if(pOutput->iTypeWriting == TYPE_GAS)
            nParts = nTotalSPH;
        else if(pOutput->iTypeWriting == TYPE_DARK)
            nParts = nTotalDark;
        else if(pOutput->iTypeWriting == TYPE_STAR)
            nParts = nTotalStar;
        }
    if(pOutput->bFloat)
        nBytes = nParts*sizeof(float);
    else if(pOutput->iBinaryOut != 6)
        nBytes = nParts*sizeof(int); // Tipsy binary does ints
    else
        nBytes = nParts*sizeof(int64_t);
    if(pOutput->bVector) 
        nBytes *= 3;
    // XXX This would be better as a parallel prefix operation       
    treeProxy[0].outputBinaryStart(*pOutput, 0, CkCallbackResumeThread());
    fIOFile = msg->file;
    Ck::IO::startSession(msg->file, nBytes + nHeader, 0,
                         CkCallback(CkIndex_Main::cbIOReady(NULL), thishandle),
                         CkCallback(CkIndex_Main::cbIOComplete(NULL), thishandle));
    delete msg;
}

/// Determine start offsets for this piece based on previous pieces.
void TreePiece::outputBinaryStart(OutputParams& params, int64_t nStart,
                                  const CkCallback& cb)
{
    nStartWrite = nStart;

    if(thisIndex == numTreePieces - 1) {
        cb.send();
        return;
        }

    int64_t nMyParts;
    if(params.iBinaryOut != 6) {
        nMyParts = myNumParticles;
    }
    else {
        if(params.iTypeWriting == TYPE_GAS)
            nMyParts = myNumSPH;
        else if(params.iTypeWriting == TYPE_DARK)
            nMyParts = myNumParticles - myNumSPH - myNumStar;
        else if(params.iTypeWriting == TYPE_STAR)
            nMyParts = myNumStar;
        else
            CkAbort("Bad writing type.");
        }
    pieces[thisIndex+1].outputBinaryStart(params, nStart + nMyParts, cb);
}

/// Session is ready; write my data.
void Main::cbIOReady(Ck::IO::SessionReadyMsg *msg)
{
    treeProxy.outputBinary(msg->session, *pOutput);
    delete msg;
}

/// All IO has completed.  Close the file
void Main::cbIOComplete(CkMessage *msg)
{
    Ck::IO::close(fIOFile, CkCallback(CkIndex_Main::cbIOClosed(NULL),
                                      thishandle));
    delete msg;
}

/// File is closed.  Update header for NChilada.  Resume main program
void Main::cbIOClosed(CkMessage *msg)
{
    FILE *outfile;
    XDR xdrs;
    
    delete msg;
    
    if(pOutput->iBinaryOut != 6) {
        int iDum;
        outfile = CmiFopen((pOutput->fileName+"."+pOutput->sTipsyExt).c_str(),
                           "r+");
        CkAssert(outfile != NULL);
        xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
        iDum = (int)nTotalParticles;
        xdr_int(&xdrs,&iDum);
        xdr_destroy(&xdrs);
        CmiFclose(outfile);
    }
    else {
        CkReductionMsg *msgMinMax;
        treeProxy.minmaxNCOut(*pOutput,
                           CkCallbackResumeThread((void *&)msgMinMax));
        pOutput->iTypeWriting >>= 1;  // Kludge quickly get the
                                         // current filename from
                                         // getNCNextOutput()
        
        FieldHeader fh;
        std::string sOutFile = getNCNextOutput(*pOutput);
        CkAssert(!sOutFile.empty());
        if(pOutput->iTypeWriting == TYPE_GAS)
            fh.numParticles = nTotalSPH;
        else if(pOutput->iTypeWriting == TYPE_DARK)
            fh.numParticles = nTotalDark;
        else if(pOutput->iTypeWriting == TYPE_STAR)
            fh.numParticles = nTotalStar;
        
        fh.time = pOutput->dTime;
        if(pOutput->bFloat)
            fh.code = Type2Code<float>::code;
        else
            fh.code = Type2Code<int64_t>::code;
        if(pOutput->bVector)
            fh.dimensions = 3;
        else
            fh.dimensions = 1;

        outfile = CmiFopen(sOutFile.c_str(), "r+");
        CkAssert(outfile != NULL);
        xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
        xdr_template(&xdrs, &fh);
        if(pOutput->bFloat) {
            if(pOutput->bVector) {
                OrientedBox<float> bbVec = *((OrientedBox<float>*)
                                             msgMinMax->getData());
                xdr_template(&xdrs, &bbVec.lesser_corner);
                xdr_template(&xdrs, &bbVec.greater_corner);
                }
            else {
                float* minmax = (float *) msgMinMax->getData();
                xdr_float(&xdrs, &minmax[0]);
                xdr_float(&xdrs, &minmax[1]);
                }
            }
        else {
            int64_t* minmax = (int64_t *) msgMinMax->getData();
            xdr_template(&xdrs, &minmax[0]);
            xdr_template(&xdrs, &minmax[1]);
            }
        xdr_destroy(&xdrs);
        CmiFclose(outfile);
        delete msgMinMax;
        /// Continue to next particle type
        sOutFile = getNCNextOutput(*pOutput);
        if(!sOutFile.empty()) {
            Ck::IO::Options opts;
            if(param.bParaWrite) {
                if(param.nIOProcessor == 0) {
                    opts.activePEs = CkNumNodes();
                    }
                else {
                    opts.activePEs = param.nIOProcessor;
                    }
                }
            else {
                opts.activePEs = 1;
                }
            Ck::IO::open(sOutFile,
                         CkCallback(CkIndex_Main::cbOpen(0), thishandle),
                         opts);
        
            return;
            }
        }
    // Finally: resume the main program.
    cbIO.send();
}

void TreePiece::minmaxNCOut(OutputParams& params, const CkCallback& cb)
{
    OrientedBox<float> bbVec;
    float minmax[2] = {FLT_MAX, -FLT_MAX};
    int64_t iminmax[2] = {INT_MAX, -INT_MAX};
    params.dm = dm; // pass cooling information

    for(unsigned int i = 1; i <= myNumParticles; ++i) {
        if((TYPETest(&myParticles[i], params.iType)
            && TYPETest(&myParticles[i], params.iTypeWriting))
           || params.iBinaryOut != 6) {
          if(params.bFloat) {
            if(params.bVector)
                bbVec.grow(params.vValue(&myParticles[i]));
            else {
                float dValue = params.dValue(&myParticles[i]);
                if(dValue < minmax[0])
                    minmax[0] = dValue;
                if(dValue > minmax[1])
                    minmax[1] = dValue;
                }
            }
          else {
              int64_t iValue = params.iValue(&myParticles[i]);
              if(iValue < iminmax[0])
                  iminmax[0] = iValue;
              if(iValue > iminmax[1])
                  iminmax[1] = iValue;
              }
          }
        }
  if(params.bFloat) {
    if(params.bVector)
        contribute(sizeof(OrientedBox<float>), &bbVec,
                   growOrientedBox_float, cb);
    else
        contribute(2*sizeof(float), minmax, minmax_float, cb);
    }
  else
      contribute(2*sizeof(int64_t), iminmax, minmax_long, cb);
}
 
/// Output a Tipsy XDR binary float array file.
void TreePiece::outputBinary(Ck::IO::Session session, OutputParams& params)
{
    XDR xdrs;
    params.dm = dm; // pass cooling information

    if(verbosity > 3)
	CkPrintf("TreePiece %d: Writing output to disk\n", thisIndex);
    
    int64_t nMyParts;
    size_t nHeader;
    if(params.iBinaryOut != 6) {
        nMyParts = myNumParticles;
        nHeader = sizeof(int);
    }
    else {
        nHeader = FieldHeader::sizeBytes;
        if(params.bFloat) {
            if(params.bVector) 
                nHeader += 6*sizeof(float);
            else
                nHeader += 2*sizeof(float);
            }
        else
            nHeader += 2*sizeof(int64_t);
        if(params.iTypeWriting == TYPE_GAS)
            nMyParts = myNumSPH;
        else if(params.iTypeWriting == TYPE_DARK)
            nMyParts = myNumParticles - myNumSPH - myNumStar;
        else if(params.iTypeWriting == TYPE_STAR)
            nMyParts = myNumStar;
        else
            CkAbort("Bad writing type.");
        }

    size_t nBytes = nMyParts*sizeof(float);
    size_t iOffset = nStartWrite*sizeof(float);
    if(!params.bFloat) {
        if(params.iBinaryOut == 6) {
            nBytes = nMyParts*sizeof(int64_t);
            iOffset = nStartWrite*sizeof(int64_t);
            }
        else {  // Tipsy only does ints.
            nBytes = nMyParts*sizeof(int);
            iOffset = nStartWrite*sizeof(int);
            }
        }
    if(params.bVector) {
        nBytes *= 3;
        iOffset *= 3;
        }
    iOffset += nHeader;

    // CkIO Offset seems to have problems; hence the following work-around:
    if(thisIndex == 0) {
        iOffset = 0;
        nBytes += nHeader;
        }

    if(nBytes == 0)     // nothing for this piece to do.
        return;

    char *buf = new char[nBytes];
    CkAssert(nBytes < UINT_MAX); // Documentation for xdrmem_create()
                                 // specifies unsigned int.  Could be
                                 // a problem.
    xdrmem_create(&xdrs, buf, nBytes, XDR_ENCODE);

    // CkIO Offset seems to have problems; hence more work-around:
    // write a dummy header
    if(thisIndex == 0) {
        if(params.iBinaryOut != 6) {
            int iDum = 0;
            xdr_int(&xdrs, &iDum);
            }
        else {
            FieldHeader fh;
            fh.time = 0.0;
            fh.numParticles = 0;
            fh.dimensions = 0;
            if(params.bFloat)
                fh.code = Type2Code<float>::code;
            else
                fh.code = Type2Code<int64_t>::code;
            xdr_template(&xdrs, &fh);
            if(params.bFloat) {
                if(params.bVector) {
                    Vector3D<float> vMin(0.0);
                    xdr_template(&xdrs, &vMin);
                    xdr_template(&xdrs, &vMin);
                    }
                else {
                    float fMin = 0.0;
                    xdr_float(&xdrs, &fMin);
                    xdr_float(&xdrs, &fMin);
                    }
                }
            else {
                int64_t iMin = 0;
                xdr_hyper(&xdrs, &iMin);
                xdr_hyper(&xdrs, &iMin);
                }
            }
        }
    // end of CkIO Offset workaround.
    
    int nOut = 0;

    for(unsigned int i = 1; i <= myNumParticles; ++i) {
        if((TYPETest(&myParticles[i], params.iType)
            && TYPETest(&myParticles[i], params.iTypeWriting))
           || params.iBinaryOut != 6) {
            if(params.bFloat) {
                if(params.bVector) {
                    Vector3D<float> vec = params.vValue(&myParticles[i]);
                    xdr_template(&xdrs, &vec);
                    }
                else {
                    float dValue = params.dValue(&myParticles[i]);
                    xdr_float(&xdrs, &dValue);
                    }
                }
            else {
                if(params.iBinaryOut == 6) {
                    int64_t iValue = params.iValue(&myParticles[i]);
                    xdr_template(&xdrs, &iValue);
                    }
                else {  // tipsy binary array only does 32 bit ints
                    int iValue = params.iValue(&myParticles[i]);
                    xdr_template(&xdrs, &iValue);
                    }
                }
            nOut++;
            }
        }
    CkAssert(nOut == nMyParts);
    xdr_destroy(&xdrs);
    Ck::IO::write(session, buf, nBytes, iOffset);
    delete [] buf;
}
