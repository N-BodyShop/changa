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
    p.u() = dTuFac*gp.temp;
    p.uPred() = dTuFac*gp.temp;
    p.vPred() = gp.vel;
    p.fBallMax() = HUGE;
    p.fESNrate() = 0.0;
    p.fTimeCoolIsOffUntil() = 0.0;
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
	int64_t startParticle = myNumParticles * myIndex;
	if(myIndex < excess) {
	    myNumParticles++;
	    startParticle += myIndex;
	    }
	else {
	    startParticle += excess;
	    }
	if(startParticle >= nTotalParticles) {
	    CkError("Bad startParticle: %d, nPart: %d, myIndex: %d, nLoading: %d\n",
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
		myParticles[i+1].rung = 0;
		myParticles[i+1].fBall = 0.0;
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
	
  contribute(cb);
}

/// @brief read iOrder file
void TreePiece::readIOrd(const std::string& filename, const CkCallback& cb)
{
    CmiInt8 nMaxOrd[3] = {0, 0, 0}; // 0 -> gas, 1 -> dark, 2 -> all

    if(nStartRead >= 0) {
        FILE *fp = CmiFopen(filename.c_str(), "r");
        CkAssert(fp != NULL);
        int64_t nTot;
        int nread;
        nread = fscanf(fp, "%ld\n", &nTot);
        CkAssert(nread == 1);
        for(int i = 0; i < nStartRead; i++) {
            int64_t dummy;
            nread = fscanf(fp, "%ld\n", &dummy);
            CkAssert(nread == 1);
            }
        for(int i = 0; i < myNumParticles; i++) {
            int64_t dummy;
            nread = fscanf(fp, "%ld\n", &dummy);
            CkAssert(nread == 1);
            myParticles[i+1].iOrder = dummy;
            if(dummy > nMaxOrd[0] && myParticles[i+1].isGas())
                nMaxOrd[0] = dummy;
            if(dummy > nMaxOrd[1] && myParticles[i+1].isDark())
                nMaxOrd[1] = dummy;
            if(dummy > nMaxOrd[2])
                nMaxOrd[2] = dummy;
            
            }
        CmiFclose(fp);
        }
    contribute(3*sizeof(CmiInt8), nMaxOrd, CkReduction::max_long, cb);
    }

/// @brief read iGasOrder file
void TreePiece::readIGasOrd(const std::string& filename, const CkCallback& cb)
{
    if(nStartRead >= 0 && myNumStar > 0) {
        FILE *fp = CmiFopen(filename.c_str(), "r");
        CkAssert(fp != NULL);
        int64_t nTot;
        int nread;
        nread = fscanf(fp, "%ld\n", &nTot);
        CkAssert(nread == 1);
        for(int i = 0; i < nStartRead; i++) {
            int64_t dummy;
            nread = fscanf(fp, "%ld\n", &dummy);
            CkAssert(nread == 1);
            }
        for(int i = 0; i < myNumParticles; i++) {
            int64_t dummy;
            nread = fscanf(fp, "%ld\n", &dummy);
            CkAssert(nread == 1);
            if(myParticles[i+1].isStar())
                myParticles[i+1].iGasOrder() = dummy;
            }
        CmiFclose(fp);
        }
    contribute(cb);
    }

/// @brief read OxMassFrac file
void TreePiece::readOxMassFrac(const std::string& filename, const CkCallback& cb)
{
    if(nStartRead >= 0 && (myNumStar > 0 || myNumSPH > 0)) {
        FILE *fp = CmiFopen(filename.c_str(), "r");
        CkAssert(fp != NULL);
        int64_t nTot;
        int nread;
        nread = fscanf(fp, "%ld\n", &nTot);
        CkAssert(nread == 1);
        for(int i = 0; i < nStartRead; i++) {
            double dummy;
            nread = fscanf(fp, "%lf\n", &dummy);
            CkAssert(nread == 1);
            }
        for(int i = 0; i < myNumParticles; i++) {
            double dummy;
            nread = fscanf(fp, "%lf\n", &dummy);
            CkAssert(nread == 1);
            if(myParticles[i+1].isGas())
                myParticles[i+1].fMFracOxygen() = dummy;
            if(myParticles[i+1].isStar())
                myParticles[i+1].fStarMFracOxygen() = dummy;
            }
        CmiFclose(fp);
        }
    contribute(cb);
    }

/// @brief read FeMassFrac file
void TreePiece::readFeMassFrac(const std::string& filename, const CkCallback& cb)
{
    if(nStartRead >= 0 && (myNumStar > 0 || myNumSPH > 0)) {
        FILE *fp = CmiFopen(filename.c_str(), "r");
        CkAssert(fp != NULL);
        int64_t nTot;
        int nread;
        nread = fscanf(fp, "%ld\n", &nTot);
        CkAssert(nread == 1);
        for(int i = 0; i < nStartRead; i++) {
            double dummy;
            nread = fscanf(fp, "%lf\n", &dummy);
            CkAssert(nread == 1);
            }
        for(int i = 0; i < myNumParticles; i++) {
            double dummy;
            nread = fscanf(fp, "%lf\n", &dummy);
            CkAssert(nread == 1);
            if(myParticles[i+1].isGas())
                myParticles[i+1].fMFracIron() = dummy;
            if(myParticles[i+1].isStar())
                myParticles[i+1].fStarMFracIron() = dummy;
            }
        CmiFclose(fp);
        }
    contribute(cb);
    }

/// @brief read ESNRate file
void TreePiece::readESNrate(const std::string& filename, const CkCallback& cb)
{
    if(nStartRead >= 0 && (myNumStar > 0 || myNumSPH > 0)) {
        FILE *fp = CmiFopen(filename.c_str(), "r");
        CkAssert(fp != NULL);
        int64_t nTot;
        int nread;
        nread = fscanf(fp, "%ld\n", &nTot);
        CkAssert(nread == 1);
        for(int i = 0; i < nStartRead; i++) {
            double dummy;
            nread = fscanf(fp, "%lf\n", &dummy);
            CkAssert(nread == 1);
            }
        for(int i = 0; i < myNumParticles; i++) {
            double dummy;
            nread = fscanf(fp, "%lf\n", &dummy);
            CkAssert(nread == 1);
            if(myParticles[i+1].isGas())
                myParticles[i+1].fESNrate() = dummy;
            if(myParticles[i+1].isStar())
                myParticles[i+1].fStarESNrate() = dummy;
            }
        CmiFclose(fp);
        }
    contribute(cb);
    }

/// @brief read Formation Mass file
void TreePiece::readMassForm(const std::string& filename, const CkCallback& cb)
{
    if(nStartRead >= 0 && myNumStar > 0) {
        FILE *fp = CmiFopen(filename.c_str(), "r");
        CkAssert(fp != NULL);
        int64_t nTot;
        int nread;
        nread = fscanf(fp, "%ld\n", &nTot);
        CkAssert(nread == 1);
        for(int i = 0; i < nStartRead; i++) {
            double dummy;
            nread = fscanf(fp, "%lf\n", &dummy);
            CkAssert(nread == 1);
            }
        for(int i = 0; i < myNumParticles; i++) {
            double dummy;
            nread = fscanf(fp, "%lf\n", &dummy);
            CkAssert(nread == 1);
            if(myParticles[i+1].isStar())
                myParticles[i+1].fMassForm() = dummy;
            }
        CmiFclose(fp);
        }
    contribute(cb);
    }

/// @brief read coolontime file
void TreePiece::readCoolOnTime(const std::string& filename, const CkCallback& cb)
{
    if(nStartRead >= 0 && myNumSPH > 0) {
        FILE *fp = CmiFopen(filename.c_str(), "r");
        CkAssert(fp != NULL);
        int64_t nTot;
        int nread;
        nread = fscanf(fp, "%ld\n", &nTot);
        CkAssert(nread == 1);
        for(int i = 0; i < nStartRead; i++) {
            double dummy;
            nread = fscanf(fp, "%lf\n", &dummy);
            CkAssert(nread == 1);
            }
        for(int i = 0; i < myNumParticles; i++) {
            double dummy;
            nread = fscanf(fp, "%lf\n", &dummy);
            CkAssert(nread == 1);
            if(myParticles[i+1].isGas())
                myParticles[i+1].fTimeCoolIsOffUntil() = dummy;
            }
        CmiFclose(fp);
        }
    contribute(cb);
    }

/// @brief read CoolArray file
void TreePiece::readCoolArray0(const std::string& filename, const CkCallback& cb)
{
#ifndef COOLING_NONE
    if(nStartRead >= 0 && myNumSPH > 0) {
        FILE *fp = CmiFopen(filename.c_str(), "r");
        CkAssert(fp != NULL);
        int64_t nTot;
        int nread;
        nread = fscanf(fp, "%ld\n", &nTot);
        CkAssert(nread == 1);
        for(int i = 0; i < nStartRead; i++) {
            double dummy;
            nread = fscanf(fp, "%lf\n", &dummy);
            CkAssert(nread == 1);
            }
        for(int i = 0; i < myNumParticles; i++) {
            double dummy;
            nread = fscanf(fp, "%lf\n", &dummy);
            CkAssert(nread == 1);
            if(myParticles[i+1].isGas())
                COOL_ARRAY0(unused1,&myParticles[i+1].CoolParticle(),unused2)
                    = dummy;
            }
        CmiFclose(fp);
        }
#endif
    contribute(cb);
    }

/// @brief read CoolArray file
void TreePiece::readCoolArray1(const std::string& filename, const CkCallback& cb)
{
#ifndef COOLING_NONE
    if(nStartRead >= 0 && myNumSPH > 0) {
        FILE *fp = CmiFopen(filename.c_str(), "r");
        CkAssert(fp != NULL);
        int64_t nTot;
        int nread;
        nread = fscanf(fp, "%ld\n", &nTot);
        CkAssert(nread == 1);
        for(int i = 0; i < nStartRead; i++) {
            double dummy;
            nread = fscanf(fp, "%lf\n", &dummy);
            CkAssert(nread == 1);
            }
        for(int i = 0; i < myNumParticles; i++) {
            double dummy;
            nread = fscanf(fp, "%lf\n", &dummy);
            CkAssert(nread == 1);
            if(myParticles[i+1].isGas())
                COOL_ARRAY1(unused1,&myParticles[i+1].CoolParticle(),unused2)
                    = dummy;
            }
        CmiFclose(fp);
        }
#endif
    contribute(cb);
    }
/// @brief read CoolArray file
void TreePiece::readCoolArray2(const std::string& filename, const CkCallback& cb)
{
#ifndef COOLING_NONE
    if(nStartRead >= 0 && myNumSPH > 0) {
        FILE *fp = CmiFopen(filename.c_str(), "r");
        CkAssert(fp != NULL);
        int64_t nTot;
        int nread;
        nread = fscanf(fp, "%ld\n", &nTot);
        CkAssert(nread == 1);
        for(int i = 0; i < nStartRead; i++) {
            double dummy;
            nread = fscanf(fp, "%lf\n", &dummy);
            CkAssert(nread == 1);
            }
        for(int i = 0; i < myNumParticles; i++) {
            double dummy;
            nread = fscanf(fp, "%lf\n", &dummy);
            CkAssert(nread == 1);
            if(myParticles[i+1].isGas())
                COOL_ARRAY2(unused1,&myParticles[i+1].CoolParticle(),unused2)
                    = dummy;
            }
        CmiFclose(fp);
        }
#endif
    contribute(cb);
    }

/// @brief read CoolArray file
void TreePiece::readCoolArray3(const std::string& filename, const CkCallback& cb)
{
#ifndef COOLING_NONE
    if(nStartRead >= 0 && myNumSPH > 0) {
        FILE *fp = CmiFopen(filename.c_str(), "r");
        CkAssert(fp != NULL);
        int64_t nTot;
        int nread;
        nread = fscanf(fp, "%ld\n", &nTot);
        CkAssert(nread == 1);
        for(int i = 0; i < nStartRead; i++) {
            double dummy;
            nread = fscanf(fp, "%lf\n", &dummy);
            CkAssert(nread == 1);
            }
        for(int i = 0; i < myNumParticles; i++) {
            double dummy;
            nread = fscanf(fp, "%lf\n", &dummy);
            CkAssert(nread == 1);
#ifndef COOLING_COSMO
            if(myParticles[i+1].isGas())
                COOL_ARRAY3(unused1,&myParticles[i+1].CoolParticle(),unused2)
                    = dummy;
#endif
            }
        CmiFclose(fp);
        }
#endif
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
    load_NC_base(filename, startParticle, myNumSPH, myParts);

    for(int i = 0; i < myNumSPH; ++i) {
        myParts[i].iType = TYPE_GAS;
        myParts[i].extraData = &mySPHParts[i];
        myParts[i].fBallMax() = HUGE;
        myParts[i].fESNrate() = 0.0;
        myParts[i].fTimeCoolIsOffUntil() = 0.0;
        }

    FieldHeader fh;
    // Density
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
    // Oxygen
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
    // Iron
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
    // Temperature
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
        }
    deleteField(fh, data);
}

static void load_NC_dark(std::string filename, int64_t startParticle,
                        int myNumDark, GravityParticle *myParts)
{
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
    for(int i = 0; i < myNumStar; ++i) {
        myParts[i].fDensity = 0.0;
        myParts[i].iType = TYPE_STAR;
        myParts[i].extraData = &myStarParts[i];
        }
    load_NC_base(filename, startParticle, myNumStar, myParts);

    FieldHeader fh;
    // Oxygen
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
	dStartTime = fh_time;

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
	int64_t startParticle = myNumParticles * myIndex;
	if(myIndex < excess) {
	    myNumParticles++;
	    startParticle += myIndex;
	    }
	else {
	    startParticle += excess;
	    }
	if(startParticle >= nTotalParticles) {
	    CkError("Bad startParticle: %d, nPart: %d, myIndex: %d, nLoading: %d\n",
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
        if(myNumDark > 0) {
            load_NC_dark(filename + "/dark", startParticle-nTotalSPH,
                         myNumDark, &myParticles[myNumSPH + 1]);
            }
        if(myNumStar > 0) {
            load_NC_star(filename + "/star", startParticle-nTotalSPH-nTotalDark,
                         myNumStar, &myParticles[myNumSPH + myNumDark + 1],
                         myStarParticles);
            }
        for(int i = 0; i < myNumParticles; ++i) {
            myParticles[i+1].rung = 0;
            myParticles[i+1].fBall = 0.0;
            myParticles[i+1].iOrder = i + startParticle;
            boundingBox.grow(myParticles[i+1].position);
        }
        
  contribute(cb);
}

/// Generic read of binary (NChilada) array format into integer
/// particle attribute
void TreePiece::readIntBinary(OutputIntParams& params, int bParaRead,
    const CkCallback& cb)
{
    FieldHeader fh;
    void *data;
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
            default:
                throw XDRException("I don't recognize the type of this field!");
                }
            }
        deleteField(fh, data);
        }
    int myNumDark = myNumParticles - myNumSPH - myNumStar;
    if((params.iType & TYPE_DARK) && (myNumDark > 0)) {
        data = readFieldData(params.fileName + "/dark/" + params.sNChilExt, fh,
            1, myNumDark, nStartRead+myNumSPH);
        for(int i = 0; i < myNumDark; ++i) {
            switch(fh.code) {
            case int32:
                params.setIValue(&myParticles[myNumSPH+i+1], static_cast<int *>(data)[i]);
                break;
            case int64:
                params.setIValue(&myParticles[myNumSPH+i+1], static_cast<int64_t *>(data)[i]);
                break;
            default:
                throw XDRException("I don't recognize the type of this field!");
                }
            }
        deleteField(fh, data);
        }
    if((params.iType & TYPE_STAR) && (myNumStar > 0)) {
        data = readFieldData(params.fileName + "/star/" + params.sNChilExt, fh,
            1, myNumStar, nStartRead+myNumSPH+myNumDark);
        for(int i = 0; i < myNumStar; ++i) {
            switch(fh.code) {
            case int32:
                params.setIValue(&myParticles[myNumSPH+myNumDark+i+1], static_cast<int *>(data)[i]);
                break;
            case int64:
                params.setIValue(&myParticles[myNumSPH+myNumDark+i+1], static_cast<int64_t *>(data)[i]);
                break;
            default:
                throw XDRException("I don't recognize the type of this field!");
                }
            }
        deleteField(fh, data);
        }
    contribute(cb);
}

/// Generic read of binary (NChilada) array format into floating point
/// particle attribute
void TreePiece::readFloatBinary(OutputParams& params, int bParaRead,
    const CkCallback& cb)
{
    FieldHeader fh;
    void *data;
    if((params.iType & TYPE_GAS) && (myNumSPH > 0)) {
        data = readFieldData(params.fileName + "/gas/" + params.sNChilExt, fh,
            1, myNumSPH, nStartRead);
        for(int i = 0; i < myNumSPH; ++i) {
            switch(fh.code) {
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
    if((params.iType & TYPE_DARK) && (myNumDark > 0)) {
        data = readFieldData(params.fileName + "/dark/" + params.sNChilExt, fh,
            1, myNumDark, nStartRead+myNumSPH);
        for(int i = 0; i < myNumDark; ++i) {
            switch(fh.code) {
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
    if((params.iType & TYPE_STAR) && (myNumStar > 0)) {
        data = readFieldData(params.fileName + "/star/" + params.sNChilExt, fh,
            1, myNumStar, nStartRead+myNumSPH+myNumDark);
        for(int i = 0; i < myNumStar; ++i) {
            switch(fh.code) {
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
        gp.temp = CoolCodeEnergyToTemperature(Cool, &p.CoolParticle(), p.u(),
                                              p.fMetals());
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
        CkError("[%d] Write dark failed, errno %d\n", CkMyPe(), errno);
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
		= new (0, 0, nPartOut, nGasOut, nStarOut)
		ParticleShuffleMsg(0, nPartOut, nGasOut, nStarOut, 0.0);
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
			      CkCallback& cb) 
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
      int result = CmiFclose(outfile);
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
	  pieces[iIndex + 1].outputIntASCII(params, 0, cb);
	  return;
	  }

    cb.send(); // We are done.
    }

// Output a Tipsy XDR binary float array file.
void TreePiece::outputBinary(OutputParams& params, // specifies
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
    XDR xdrs;
    float *afOut;	// array for oneNode I/O
    Vector3D<float> *avOut;	// array for one node I/O
    params.dm = dm; // pass cooling information
    int iDum;
    
    if(thisIndex==0) {
	if(verbosity > 2)
	    ckout << "TreePiece " << thisIndex << ": Writing header for output file" << endl;
        if(params.iBinaryOut != 6) {
            outfile = CmiFopen((params.fileName+"."+params.sTipsyExt).c_str(), "w");
            CkAssert(outfile != NULL);
            xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
            iDum = (int)nTotalParticles;
            xdr_int(&xdrs,&iDum);
            xdr_destroy(&xdrs);
            CmiFclose(outfile);
            }
	}
    
    if(verbosity > 3)
	ckout << "TreePiece " << thisIndex << ": Writing output to disk" << endl;
    
  if(params.iBinaryOut != 6) {
    if(bParaWrite) {
        outfile = CmiFopen((params.fileName+"."+params.sTipsyExt).c_str(), "a");
	if(outfile == NULL)
	    ckerr << "Treepiece " << thisIndex << " failed to open "
		  << params.fileName.c_str() << " : " << errno << endl;
	CkAssert(outfile != NULL);
	xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
	}
    else {
	if(params.bVector && packed)
	    avOut = new Vector3D<float>[myNumParticles];
	else
	    afOut = new float[myNumParticles];
	}
    for(unsigned int i = 1; i <= myNumParticles; ++i) {
	Vector3D<float> vOut;
	float fOut;
	if(params.bVector) {
	    vOut = params.vValue(&myParticles[i]);
	    if(!packed){
		if(cnt==0)
		    fOut = vOut.x;
		if(cnt==1)
		    fOut = vOut.y;
		if(cnt==2)
		    fOut = vOut.z;
		}
	    }
	else
	    fOut = params.dValue(&myParticles[i]);
	
	if(bParaWrite) {
	    if(params.bVector && packed){
		if(!xdr_float(&xdrs,&vOut.x)) {
		    ckerr << "TreePiece " << thisIndex << ": Error writing array to disk, aborting" << endl;
		    CkAbort("Badness");
		    }
		if(!xdr_float(&xdrs,&vOut.y)) {
		    ckerr << "TreePiece " << thisIndex << ": Error writing array to disk, aborting" << endl;
		    CkAbort("Badness");
		    }
		if(!xdr_float(&xdrs,&vOut.z)) {
		    ckerr << "TreePiece " << thisIndex << ": Error writing array to disk, aborting" << endl;
		    CkAbort("Badness");
		    }
		}
	    else {
		if(!xdr_float(&xdrs,&fOut)) {
		    ckerr << "TreePiece " << thisIndex << ": Error writing array to disk, aborting" << endl;
		    CkAbort("Badness");
		    }
		}
	    }
	else {
	    if(params.bVector && packed)
		avOut[i-1] = vOut;
	    else
		afOut[i-1] = fOut;
	    }
	}
    cnt++;
    if(cnt == 3 || !params.bVector)
	cnt = 0;
    
    if(bParaWrite) {
	xdr_destroy(&xdrs);
	int result = CmiFclose(outfile);
	if(result != 0)
	    ckerr << "Bad close: " << strerror(errno) << endl;
	CkAssert(result == 0);
	
	if(thisIndex!=(int)numTreePieces-1) {
	    pieces[thisIndex + 1].outputBinary(params, bParaWrite, cb);
	    return;
	    }
	
	if(packed || !params.bVector || (!packed && cnt==0)) {
	    cb.send(); // We are done.
	    return;
	    }
	// go through pieces again for unpacked vector.
	pieces[0].outputBinary(params, bParaWrite, cb);
	}
    else {
	int bDone = packed || !params.bVector || (!packed && cnt==0); // flag for last time
	if(params.bVector && packed) {
	    pieces[0].oneNodeOutBinVec(params, avOut, myNumParticles, thisIndex,
				       bDone, cb);
	    delete [] avOut;
	    }
	else {
	    pieces[0].oneNodeOutBinArr(params, afOut, myNumParticles, thisIndex,
				       bDone, cb);
	    delete [] afOut;
	    }
	}
    }
  else { // N-Chilada output
      int nOutDark = 0;
      int nOutGas = 0;
      int nOutStar = 0;
      if(params.bVector) {
          Vector3D<float> *avOutGas;
          Vector3D<float> *avOutDark;
          Vector3D<float> *avOutStar;
          if(params.iType & TYPE_GAS)
              avOutGas = new Vector3D<float>[myNumSPH];
          if(params.iType & TYPE_DARK)
              avOutDark = new Vector3D<float>[myNumParticles - myNumSPH - myNumStar];
          if(params.iType & TYPE_STAR)
              avOutStar = new Vector3D<float>[myNumStar];
          for(unsigned int i = 1; i <= myNumParticles; ++i) {
              if(TYPETest(&myParticles[i], params.iType)) {
                  if(myParticles[i].iType & TYPE_GAS) {
                      avOutGas[nOutGas] = params.vValue(&myParticles[i]);
                      nOutGas++;
                      }
                  if(myParticles[i].iType & TYPE_DARK) {
                      avOutDark[nOutDark] = params.vValue(&myParticles[i]);
                      nOutDark++;
                      }
                  if(myParticles[i].iType & TYPE_STAR) {
                      avOutStar[nOutStar] = params.vValue(&myParticles[i]);
                      nOutStar++;
                      }
                  }
              }
          pieces[0].oneNodeOutNCVec(params, avOutGas, avOutDark, avOutStar,
                                    nOutGas, nOutDark, nOutStar, thisIndex, cb);
          if(params.iType & TYPE_GAS)
              delete [] avOutGas;
          if(params.iType & TYPE_DARK)
              delete [] avOutDark;
          if(params.iType & TYPE_STAR)
              delete[] avOutStar;
          }
      else {
          float *afOutGas;
          float *afOutDark;
          float *afOutStar;
          if(params.iType & TYPE_GAS)
              afOutGas = new float[myNumSPH];
          if(params.iType & TYPE_DARK)
              afOutDark = new float[myNumParticles - myNumSPH - myNumStar];
          if(params.iType & TYPE_STAR)
              afOutStar = new float[myNumStar];
          for(unsigned int i = 1; i <= myNumParticles; ++i) {
              if(TYPETest(&myParticles[i], params.iType)) {
                  if(myParticles[i].iType & TYPE_GAS) {
                      afOutGas[nOutGas] = params.dValue(&myParticles[i]);
                      nOutGas++;
                      }
                  if(myParticles[i].iType & TYPE_DARK) {
                      afOutDark[nOutDark] = params.dValue(&myParticles[i]);
                      nOutDark++;
                      }
                  if(myParticles[i].iType & TYPE_STAR) {
                      afOutStar[nOutStar] = params.dValue(&myParticles[i]);
                      nOutStar++;
                      }
                  }
              }
          pieces[0].oneNodeOutNCArr(params, afOutGas, afOutDark, afOutStar,
                                    nOutGas, nOutDark, nOutStar, thisIndex, cb);
          if(params.iType & TYPE_GAS)
              delete [] afOutGas;
          if(params.iType & TYPE_DARK)
              delete [] afOutDark;
          if(params.iType & TYPE_STAR)
              delete [] afOutStar;
          }
      }
}

/// @brief Receives an array of vectors to write out in NChilada format.
///
/// Assumes we are working on one node and will keep the file open
/// from one array call to the next for efficiency's sake.
void TreePiece::oneNodeOutNCVec(OutputParams& params,
                                Vector3D<float>* avOutGas, // array to be output
                                Vector3D<float>* avOutDark, // array to be output
                                Vector3D<float>* avOutStar, // array to be output
                                int nGas, // number of elements in avOut
                                int nDark, // number of elements in avOut
                                int nStar, // number of elements in avOut
                                int iIndex, // treepiece which called me
                                CkCallback& cb) 
{
    static FILE *outfileGas;
    static XDR xdrsGas;
    static OrientedBox<float> bbGas;
    static FILE *outfileDark;
    static XDR xdrsDark;
    static OrientedBox<float> bbDark;
    static FILE *outfileStar;
    static XDR xdrsStar;
    static OrientedBox<float> bbStar;
    FieldHeader fh;
    if(iIndex == 0) { // first time
        fh.time = params.dTime;
        fh.dimensions = 3;
        fh.code = Type2Code<float>::code;
        if((params.iType & TYPE_GAS) && (nTotalSPH > 0)) {
            outfileGas = CmiFopen((params.fileName+"/gas/"+params.sNChilExt).c_str(), "w");
            CkAssert(outfileGas != NULL);
            xdrstdio_create(&xdrsGas, outfileGas, XDR_ENCODE);
            fh.numParticles = nTotalSPH;
            xdr_template(&xdrsGas, &fh);
            bbGas.reset();
            xdr_template(&xdrsGas, &bbGas.lesser_corner);
            xdr_template(&xdrsGas, &bbGas.greater_corner);
            }
        if((params.iType & TYPE_DARK) && (nTotalDark > 0)) {
            outfileDark = CmiFopen((params.fileName+"/dark/"+params.sNChilExt).c_str(), "w");
            CkAssert(outfileDark != NULL);
            xdrstdio_create(&xdrsDark, outfileDark, XDR_ENCODE);
            fh.numParticles = nTotalDark;
            xdr_template(&xdrsDark, &fh);
            bbDark.reset();
            xdr_template(&xdrsDark, &bbDark.lesser_corner);
            xdr_template(&xdrsDark, &bbDark.greater_corner);
            }
        if((params.iType & TYPE_STAR) && (nTotalStar > 0)) {
            outfileStar = CmiFopen((params.fileName+"/star/"+params.sNChilExt).c_str(), "w");
            CkAssert(outfileStar != NULL);
            xdrstdio_create(&xdrsStar, outfileStar, XDR_ENCODE);
            fh.numParticles = nTotalStar;
            xdr_template(&xdrsStar, &fh);
            bbStar.reset();
            xdr_template(&xdrsStar, &bbStar.lesser_corner);
            xdr_template(&xdrsStar, &bbStar.greater_corner);
        }
    }
    if(params.iType & TYPE_GAS) {
        bool result = writeField(&xdrsGas, nGas, avOutGas);
        for(u_int64_t i = 0; i < nGas; ++i)
            bbGas.grow(avOutGas[i]);
        CkAssert(result);
        }
    if(params.iType & TYPE_DARK) {
        bool result = writeField(&xdrsDark, nDark, avOutDark);
        for(u_int64_t i = 0; i < nDark; ++i)
            bbDark.grow(avOutDark[i]);
        CkAssert(result);
        }
    if(params.iType & TYPE_STAR) {
        bool result = writeField(&xdrsStar, nStar, avOutStar);
        for(u_int64_t i = 0; i < nStar; ++i)
            bbStar.grow(avOutStar[i]);
        CkAssert(result);
        }
    if(iIndex!=(int)numTreePieces-1) {
	  pieces[iIndex + 1].outputBinary(params, 0, cb);
	  return;
	  }
    else {
        if((params.iType & TYPE_GAS) && (nTotalSPH > 0)) {
            xdr_setpos(&xdrsGas, fh.sizeBytes);
            xdr_template(&xdrsGas, &bbGas.lesser_corner);
            xdr_template(&xdrsGas, &bbGas.greater_corner);
            xdr_destroy(&xdrsGas);
            CmiFclose(outfileGas);
            }
        if((params.iType & TYPE_DARK) && (nTotalDark > 0)) {
            xdr_setpos(&xdrsDark, fh.sizeBytes);
            xdr_template(&xdrsDark, &bbDark.lesser_corner);
            xdr_template(&xdrsDark, &bbDark.greater_corner);
            xdr_destroy(&xdrsDark);
            CmiFclose(outfileDark);
            }
        if((params.iType & TYPE_STAR) && (nTotalStar > 0)) {
            xdr_setpos(&xdrsStar, fh.sizeBytes);
            xdr_template(&xdrsStar, &bbStar.lesser_corner);
            xdr_template(&xdrsStar, &bbStar.greater_corner);
            xdr_destroy(&xdrsStar);
            CmiFclose(outfileStar);
            }
        cb.send(); // We are done.
        return;
    }
}

/// @brief Receives an array of floats to write out in NChilada format.
///
/// Assumes we are working on one node and will keep the file open
/// from one array call to the next for efficiency's sake.
void TreePiece::oneNodeOutNCArr(OutputParams& params,
                                float* afOutGas, // array to be output
                                float* afOutDark, // array to be output
                                float* afOutStar, // array to be output
                                int nGas, // number of elements in avOut
                                int nDark, // number of elements in avOut
                                int nStar, // number of elements in avOut
                                int iIndex, // treepiece which called me
                                CkCallback& cb) 
{
    static FILE *outfileGas;
    static XDR xdrsGas;
    static float bbGas[2];
    static FILE *outfileDark;
    static XDR xdrsDark;
    static float bbDark[2];
    static FILE *outfileStar;
    static XDR xdrsStar;
    static float bbStar[2];
    FieldHeader fh;
    if(iIndex == 0) { // first time
        fh.time = params.dTime;
        fh.dimensions = 1;
        fh.code = Type2Code<float>::code;
        if((params.iType & TYPE_GAS) && (nTotalSPH > 0)) {
            outfileGas = CmiFopen((params.fileName+"/gas/"+params.sNChilExt).c_str(), "w");
            CkAssert(outfileGas != NULL);
            xdrstdio_create(&xdrsGas, outfileGas, XDR_ENCODE);
            fh.numParticles = nTotalSPH;
            xdr_template(&xdrsGas, &fh);
            bbGas[0] = HUGE_VAL;
            bbGas[1] = -HUGE_VAL;
            xdr_template(&xdrsGas, &bbGas[0]);
            xdr_template(&xdrsGas, &bbGas[1]);
            }
        if((params.iType & TYPE_DARK) && (nTotalDark > 0)) {
            outfileDark = CmiFopen((params.fileName+"/dark/"+params.sNChilExt).c_str(), "w");
            CkAssert(outfileDark != NULL);
            xdrstdio_create(&xdrsDark, outfileDark, XDR_ENCODE);
            fh.numParticles = nTotalDark;
            xdr_template(&xdrsDark, &fh);
            bbDark[0] = HUGE_VAL;
            bbDark[1] = -HUGE_VAL;
            xdr_template(&xdrsDark, &bbDark[0]);
            xdr_template(&xdrsDark, &bbDark[1]);
            }
        if((params.iType & TYPE_STAR) && (nTotalStar > 0)) {
            outfileStar = CmiFopen((params.fileName+"/star/"+params.sNChilExt).c_str(), "w");
            CkAssert(outfileStar != NULL);
            xdrstdio_create(&xdrsStar, outfileStar, XDR_ENCODE);
            fh.numParticles = nTotalStar;
            xdr_template(&xdrsStar, &fh);
            bbStar[0] = HUGE_VAL;
            bbStar[1] = -HUGE_VAL;
            xdr_template(&xdrsStar, &bbStar[0]);
            xdr_template(&xdrsStar, &bbStar[1]);
        }
    }
    if(params.iType & TYPE_GAS) {
        bool result = writeField(&xdrsGas, nGas, afOutGas);
        CkAssert(result);
        for(u_int64_t i = 0; i < nGas; ++i) {
            if(afOutGas[i] < bbGas[0])
                bbGas[0] = afOutGas[i];
            if(afOutGas[i] > bbGas[1])
                bbGas[1] = afOutGas[i];
            }
        }
    if(params.iType & TYPE_DARK) {
        bool result = writeField(&xdrsDark, nDark, afOutDark);
        CkAssert(result);
        for(u_int64_t i = 0; i < nDark; ++i) {
            if(afOutDark[i] < bbDark[0])
                bbDark[0] = afOutDark[i];
            if(afOutDark[i] > bbDark[1])
                bbDark[1] = afOutDark[i];
            }
        }
    if(params.iType & TYPE_STAR) {
        bool result = writeField(&xdrsStar, nStar, afOutStar);
        CkAssert(result);
        for(u_int64_t i = 0; i < nStar; ++i) {
            if(afOutStar[i] < bbStar[0])
                bbStar[0] = afOutStar[i];
            if(afOutStar[i] > bbStar[1])
                bbStar[1] = afOutStar[i];
            }
        }
    if(iIndex!=(int)numTreePieces-1) {
	  pieces[iIndex + 1].outputBinary(params, 0, cb);
	  return;
	  }
    else {
        if((params.iType & TYPE_GAS) && (nTotalSPH > 0)) {
            xdr_setpos(&xdrsGas, fh.sizeBytes);
            xdr_template(&xdrsGas, &bbGas[0]);
            xdr_template(&xdrsGas, &bbGas[1]);
            xdr_destroy(&xdrsGas);
            CmiFclose(outfileGas);
            }
        if((params.iType & TYPE_DARK) && (nTotalDark > 0)) {
            xdr_setpos(&xdrsDark, fh.sizeBytes);
            xdr_template(&xdrsDark, &bbDark[0]);
            xdr_template(&xdrsDark, &bbDark[1]);
            xdr_destroy(&xdrsDark);
            CmiFclose(outfileDark);
            }
        if((params.iType & TYPE_STAR) && (nTotalStar > 0)) {
            xdr_setpos(&xdrsStar, fh.sizeBytes);
            xdr_template(&xdrsStar, &bbStar[0]);
            xdr_template(&xdrsStar, &bbStar[1]);
            xdr_destroy(&xdrsStar);
            CmiFclose(outfileStar);
            }
        cb.send(); // We are done.
        return;
    }
}

/// @brief Receives an array of ints to write out in NChilada format.
///
/// Assumes we are working on one node and will keep the file open
/// from one array call to the next for efficiency's sake.
void TreePiece::oneNodeOutNCInt(OutputIntParams& params,
                                int* aiOutGas, // array to be output
                                int* aiOutDark, // array to be output
                                int* aiOutStar, // array to be output
                                int nGas, // number of elements in avOut
                                int nDark, // number of elements in avOut
                                int nStar, // number of elements in avOut
                                int iIndex, // treepiece which called me
                                CkCallback& cb) 
{
    static FILE *outfileGas;
    static XDR xdrsGas;
    static int bbGas[2];
    static FILE *outfileDark;
    static XDR xdrsDark;
    static int bbDark[2];
    static FILE *outfileStar;
    static XDR xdrsStar;
    static int bbStar[2];
    FieldHeader fh;
    if(iIndex == 0) { // first time
        fh.time = params.dTime;
        fh.dimensions = 1;
        fh.code = Type2Code<int>::code;
        if((params.iType & TYPE_GAS) && (nTotalSPH > 0)) {
            outfileGas = CmiFopen((params.fileName+"/gas/"+params.sNChilExt).c_str(), "w");
            CkAssert(outfileGas != NULL);
            xdrstdio_create(&xdrsGas, outfileGas, XDR_ENCODE);
            fh.numParticles = nTotalSPH;
            xdr_template(&xdrsGas, &fh);
            bbGas[0] = INT_MAX;
            bbGas[1] = -INT_MAX;
            xdr_template(&xdrsGas, &bbGas[0]);
            xdr_template(&xdrsGas, &bbGas[1]);
            }
        if((params.iType & TYPE_DARK) && (nTotalDark > 0)) {
            outfileDark = CmiFopen((params.fileName+"/dark/"+params.sNChilExt).c_str(), "w");
            CkAssert(outfileDark != NULL);
            xdrstdio_create(&xdrsDark, outfileDark, XDR_ENCODE);
            fh.numParticles = nTotalDark;
            xdr_template(&xdrsDark, &fh);
            bbDark[0] = INT_MAX;
            bbDark[1] = -INT_MAX;
            xdr_template(&xdrsDark, &bbDark[0]);
            xdr_template(&xdrsDark, &bbDark[1]);
            }
        if((params.iType & TYPE_STAR) && (nTotalStar > 0)) {
            outfileStar = CmiFopen((params.fileName+"/star/"+params.sNChilExt).c_str(), "w");
            CkAssert(outfileStar != NULL);
            xdrstdio_create(&xdrsStar, outfileStar, XDR_ENCODE);
            fh.numParticles = nTotalStar;
            xdr_template(&xdrsStar, &fh);
            bbStar[0] = INT_MAX;
            bbStar[1] = -INT_MAX;
            xdr_template(&xdrsStar, &bbStar[0]);
            xdr_template(&xdrsStar, &bbStar[1]);
        }
    }
    if(params.iType & TYPE_GAS) {
        bool result = writeField(&xdrsGas, nGas, aiOutGas);
        CkAssert(result);
        for(u_int64_t i = 0; i < nGas; ++i) {
            if(aiOutGas[i] < bbGas[0])
                bbGas[0] = aiOutGas[i];
            if(aiOutGas[i] > bbGas[1])
                bbGas[1] = aiOutGas[i];
            }
        }
    if(params.iType & TYPE_DARK) {
        bool result = writeField(&xdrsDark, nDark, aiOutDark);
        CkAssert(result);
        for(u_int64_t i = 0; i < nDark; ++i) {
            if(aiOutDark[i] < bbDark[0])
                bbDark[0] = aiOutDark[i];
            if(aiOutDark[i] > bbDark[1])
                bbDark[1] = aiOutDark[i];
            }
        }
    if(params.iType & TYPE_STAR) {
        bool result = writeField(&xdrsStar, nStar, aiOutStar);
        CkAssert(result);
        for(u_int64_t i = 0; i < nStar; ++i) {
            if(aiOutStar[i] < bbStar[0])
                bbStar[0] = aiOutStar[i];
            if(aiOutStar[i] > bbStar[1])
                bbStar[1] = aiOutStar[i];
            }
        }
    if(iIndex!=(int)numTreePieces-1) {
	  pieces[iIndex + 1].outputIntBinary(params, 0, cb);
	  return;
	  }
    else {
        if((params.iType & TYPE_GAS) && (nTotalSPH > 0)) {
            xdr_setpos(&xdrsGas, fh.sizeBytes);
            xdr_template(&xdrsGas, &bbGas[0]);
            xdr_template(&xdrsGas, &bbGas[1]);
            xdr_destroy(&xdrsGas);
            CmiFclose(outfileGas);
            }
        if((params.iType & TYPE_DARK) && (nTotalDark > 0)) {
            xdr_setpos(&xdrsDark, fh.sizeBytes);
            xdr_template(&xdrsDark, &bbDark[0]);
            xdr_template(&xdrsDark, &bbDark[1]);
            xdr_destroy(&xdrsDark);
            CmiFclose(outfileDark);
            }
        if((params.iType & TYPE_STAR) && (nTotalStar > 0)) {
            xdr_setpos(&xdrsStar, fh.sizeBytes);
            xdr_template(&xdrsStar, &bbStar[0]);
            xdr_template(&xdrsStar, &bbStar[1]);
            xdr_destroy(&xdrsStar);
            CmiFclose(outfileStar);
            }
        cb.send(); // We are done.
        return;
    }
}

// Receives an array of vectors to write out in Binary format
// Assumed to be called from outputBinary() and will continue with the
// next tree piece unless "bDone".
void TreePiece::oneNodeOutBinVec(OutputParams& params,
				 Vector3D<float>* avOut, // array to be output
				 int nPart, // number of elements in avOut
				 int iIndex, // treepiece which called me
				 int bDone, // Last call
				 CkCallback& cb) 
{
    FILE* outfile = CmiFopen((params.fileName+"."+params.sTipsyExt).c_str(), "a");
    XDR xdrs;
    if(outfile == NULL)
	ckerr << "Treepiece " << thisIndex << " failed to open "
	      << params.fileName.c_str() << " : " << errno << endl;
    CkAssert(outfile != NULL);
    xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
    for(int i = 0; i < nPart; ++i) {
	if(!xdr_float(&xdrs,&(avOut[i].x))) {
	    ckerr << "TreePiece " << thisIndex << ": Error writing array to disk, aborting" << endl;
	    CkAbort("Badness");
	    }
	if(!xdr_float(&xdrs,&(avOut[i].y))) {
	    ckerr << "TreePiece " << thisIndex << ": Error writing array to disk, aborting" << endl;
	    CkAbort("Badness");
	    }
	if(!xdr_float(&xdrs,&(avOut[i].z))) {
	    ckerr << "TreePiece " << thisIndex << ": Error writing array to disk, aborting" << endl;
	    CkAbort("Badness");
	    }
	}
    xdr_destroy(&xdrs);
    int result = CmiFclose(outfile);
    if(result != 0)
	ckerr << "Bad close: " << strerror(errno) << endl;
    CkAssert(result == 0);

    if(iIndex!=(int)numTreePieces-1) {
	  pieces[iIndex + 1].outputBinary(params, 0, cb);
	  return;
	  }

    if(bDone) {
	  cb.send(); // We are done.
	  return;
	  }
    CkAbort("packed array Done logic wrong");
    }

// Receives an array of floats to write out in Binary format
// Assumed to be called from outputBinary() and will continue with the
// next tree piece unless "bDone"

void TreePiece::oneNodeOutBinArr(OutputParams& params,
			      float *afOut, // array to be output
			      int nPart, // length of afOut
			      int iIndex, // treepiece which called me
			      int bDone, // Last call
			      CkCallback& cb) 
{
    FILE* outfile = CmiFopen((params.fileName+"."+params.sTipsyExt).c_str(), "a");
    XDR xdrs;
    if(outfile == NULL)
	ckerr << "Treepiece " << thisIndex << " failed to open "
	      << params.fileName.c_str() << " : " << errno << endl;
    CkAssert(outfile != NULL);
    xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
    for(int i = 0; i < nPart; ++i) {
	if(!xdr_float(&xdrs,&(afOut[i]))) {
	  ckerr << "TreePiece " << thisIndex << ": Error writing array to disk, aborting" << endl;
	    CkAbort("Badness");
	    }
	}
    xdr_destroy(&xdrs);
    int result = CmiFclose(outfile);
    if(result != 0)
	ckerr << "Bad close: " << strerror(errno) << endl;
    CkAssert(result == 0);

    if(iIndex!=(int)numTreePieces-1) {
	  pieces[iIndex + 1].outputBinary(params, 0, cb);
	  return;
	  }

    if(bDone) {
	  cb.send(); // We are done.
	  return;
	  }
    // go through pieces again for unpacked vector.
    pieces[0].outputBinary(params, 0, cb);
    }

void TreePiece::outputIntBinary(OutputIntParams& params, // specifies
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
    XDR xdrs;
  
    if(thisIndex==0) {
        int iDum;
        
        if(verbosity > 2)
            ckout << "TreePiece " << thisIndex << ": Writing header for output file" << endl;
        if(params.iBinaryOut != 6) {
            outfile = CmiFopen((params.fileName+"."+params.sTipsyExt).c_str(), "wb");
            CkAssert(outfile != NULL);
            xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
            iDum = (int)nTotalParticles;
            xdr_int(&xdrs,&iDum);
            xdr_destroy(&xdrs);
            CmiFclose(outfile);
            }
    }
	
    if(verbosity > 3)
        ckout << "TreePiece " << thisIndex << ": Writing output to disk" << endl;
	
  if(params.iBinaryOut != 6) {
    if(bParaWrite) {
        outfile = CmiFopen((params.fileName+"."+params.sTipsyExt).c_str(), "a");
        if(outfile == NULL)
	    ckerr << "Treepiece " << thisIndex << " failed to open "
		  << params.fileName.c_str() << " : " << errno << endl;
        CkAssert(outfile != NULL);
	xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
        }
    else {
        aiOut = new int[myNumParticles];
        }
    for(unsigned int i = 1; i <= myNumParticles; ++i) {
        int iOut;
        iOut = params.iValue(&myParticles[i]);
      
        if(bParaWrite) {
            if(!xdr_int(&xdrs,&iOut)) {
                ckerr << "TreePiece " << thisIndex << ": Error writing array to disk, aborting" << endl;
                CkAbort("Badness");
                }
            }
        else {
            aiOut[i-1] = iOut;
            }
          }
     
    if(bParaWrite) {
	xdr_destroy(&xdrs);
        int result = CmiFclose(outfile);
        if(result != 0)
	    ckerr << "Bad close: " << strerror(errno) << endl;
        CkAssert(result == 0);

        if(thisIndex!=(int)numTreePieces-1) {
            pieces[thisIndex + 1].outputIntBinary(params, bParaWrite, cb);
            return;
            }

        cb.send(); // We are done.
        return;
        }
    else {
        pieces[0].oneNodeOutBinInt(params, aiOut, myNumParticles, thisIndex,
                                   cb);
        delete [] aiOut;
      }
    }
  else { // N-Chilada output
      int nOutDark = 0;
      int nOutGas = 0;
      int nOutStar = 0;
      int *aiOutGas;
      int *aiOutDark;
      int *aiOutStar;
      if(params.iType & TYPE_GAS)
          aiOutGas = new int[myNumSPH];
      if(params.iType & TYPE_DARK)
          aiOutDark = new int[myNumParticles - myNumSPH - myNumStar];
      if(params.iType & TYPE_STAR)
          aiOutStar = new int[myNumStar];
      for(unsigned int i = 1; i <= myNumParticles; ++i) {
          if(TYPETest(&myParticles[i], params.iType)) {
              if(myParticles[i].iType & TYPE_GAS) {
                  aiOutGas[nOutGas] = params.iValue(&myParticles[i]);
                  nOutGas++;
                  }
              if(myParticles[i].iType & TYPE_DARK) {
                  aiOutDark[nOutDark] = params.iValue(&myParticles[i]);
                  nOutDark++;
                  }
              if(myParticles[i].iType & TYPE_STAR) {
                  aiOutStar[nOutStar] = params.iValue(&myParticles[i]);
                  nOutStar++;
                  }
              }
          }
      pieces[0].oneNodeOutNCInt(params, aiOutGas, aiOutDark, aiOutStar,
                                nOutGas, nOutDark, nOutStar, thisIndex, cb);
      if(params.iType & TYPE_GAS)
          delete [] aiOutGas;
      if(params.iType & TYPE_DARK)
          delete [] aiOutDark;
      if(params.iType & TYPE_STAR)
          delete [] aiOutStar;
      }
}

// Receives an array of ints to write out in binary format
// Assumed to be called from outputIntBinary() and will continue with the
// next tree piece.

void TreePiece::oneNodeOutBinInt(OutputIntParams& params,
			      int *aiOut, // array to be output
			      int nPart, // length of adOut
			      int iIndex, // treepiece which called me
			      CkCallback& cb) 
{
    FILE* outfile = CmiFopen((params.fileName+"."+params.sTipsyExt).c_str(), "a");
    XDR xdrs;
    if(outfile == NULL)
	ckerr << "Treepiece " << thisIndex << " failed to open "
	      << params.fileName.c_str() << " : " << errno << endl;
    CkAssert(outfile != NULL);
    xdrstdio_create(&xdrs, outfile, XDR_ENCODE);
    for(int i = 0; i < nPart; ++i) {
	if(!xdr_int(&xdrs,&(aiOut[i]))) {
            ckerr << "TreePiece " << thisIndex << ": Error writing array to disk, aborting" << endl;
	    CkAbort("Badness");
	    }
        }
    xdr_destroy(&xdrs);
    int result = CmiFclose(outfile);
    if(result != 0)
	ckerr << "Bad close: " << strerror(errno) << endl;
    CkAssert(result == 0);
    if(iIndex!=(int)numTreePieces-1) {
	  pieces[iIndex + 1].outputIntBinary(params, 0, cb);
	  return;
	  }
    cb.send(); // We are done.
    }
