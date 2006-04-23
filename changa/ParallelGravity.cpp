/** @file ParallelGravity.cpp
 */
 
#include <iostream>

//#include <popt.h>
#include <unistd.h>

#include "Sorter.h"
#include "ParallelGravity.h"
#include "DataManager.h"
#include "CacheManager.h"
#include "TipsyFile.h"
#include "param.h"

extern char *optarg;
extern int optind, opterr, optopt;

using namespace std;

CProxy_Main mainChare;
int verbosity;
CProxy_TreePiece treeProxy;
//CkReduction::reducerType callbackReduction;
bool _cache;
int _nocache;
int _cacheLineDepth;
unsigned int _yieldPeriod;
DomainsDec domainDecomposition;
GenericTrees useTree;
CProxy_TreePiece streamingProxy;
int _prefetch;
int _numChunks;
int _randChunks;

CkGroupID dataManagerID;
CkArrayID treePieceID;

void _Leader(void) {
    puts("USAGE: ParallelGravity [SETTINGS | FLAGS] [SIM_FILE]");
    puts("SIM_FILE: Configuration file of a particular simulation, which");
    puts("          includes desired settings and relevant input and");
    puts("          output files. Settings specified in this file override");
    puts("          the default settings.");
    puts("SETTINGS");
    puts("or FLAGS: Command line settings or flags for a simulation which");
    puts("          will override any defaults and any settings or flags");
    puts("          specified in the SIM_FILE.");
}


void _Trailer(void) {
	puts("(see man page (Ha!) for more information)");
}

Main::Main(CkArgMsg* m) {
	_cache = true;
	mainChare = thishandle;
	
	prmInitialize(&prm,_Leader,_Trailer);

	param.dTimeStep = 0.0;
	prmAddParam(prm, "dTimeStep", paramDouble, &param.dTimeStep,
		    sizeof(double),"dT", "Base Timestep for integration");
	param.nSteps = 1;
	prmAddParam(prm, "nSteps", paramInt, &param.nSteps,
		    sizeof(int),"n", "Number of Timesteps");
	param.dSoft = 0.0;
	prmAddParam(prm, "dSoftening", paramDouble, &param.dSoft,
		    sizeof(double),"S", "Gravitational softening");
	theta = 0.7;
	prmAddParam(prm, "dTheta", paramDouble, &theta,
		    sizeof(double),"t", "Opening angle");
	printBinaryAcc=1;
	prmAddParam(prm, "bPrintBinary", paramBool, &printBinaryAcc,
		    sizeof(int),"z", "Print accelerations in Binary");
	param.achInFile[0] = '\0';
	prmAddParam(prm,"achInFile",paramString,param.achInFile,
		    256, "I", "input file name (or base file name)");
	param.bStaticTest = 0;
	prmAddParam(prm, "bStaticTest", paramBool, &param.bStaticTest,
		    sizeof(int),"st", "Static test of performance");
	
	_randChunks = 0;
	prmAddParam(prm, "bRandChunks", paramBool, &_randChunks,
		    sizeof(int),"rand", "Randomize the order of remote chunk computation");
	
	_nocache = 0;
	prmAddParam(prm, "bNoCache", paramBool, &_nocache,
		    sizeof(int),"nc", "Disable the CacheManager caching behaviour");
	
        verbosity = 0;
	prmAddParam(prm, "iVerbosity", paramInt, &verbosity,
		    sizeof(int),"v", "Verbosity");
	numTreePieces = CkNumPes();
	prmAddParam(prm, "nTreePieces", paramInt, &numTreePieces,
		    sizeof(int),"p", "Number of TreePieces");
	bucketSize = 12;
	prmAddParam(prm, "nTreePieces", paramInt, &bucketSize,
		    sizeof(int),"b", "Particles per Bucket");
	_numChunks = 1;
	prmAddParam(prm, "nChunks", paramInt, &_numChunks,
		    sizeof(int),"c", "Chunks per TreePiece");
	_yieldPeriod=100000000;
	prmAddParam(prm, "nYield", paramInt, &_yieldPeriod,
		    sizeof(int),"y", "Yield Period");
	_cacheLineDepth=1;
	prmAddParam(prm, "nCacheDepth", paramInt, &_cacheLineDepth,
		    sizeof(double),"d", "Chunks per TreePiece");
	_prefetch=false;
	prmAddParam(prm, "bPrefetch", paramBool, &_prefetch,
		    sizeof(int),"f", "Enable prefetching in the cache");
	int cacheSize = 100000000;
	prmAddParam(prm, "nCacheSize", paramInt, &cacheSize,
		    sizeof(int),"s", "Size of cache");
	domainDecomposition=SFC_dec;
	prmAddParam(prm, "nDomainDecompose", paramInt, &domainDecomposition,
		    sizeof(int),"D", "Kind of domain decomposition of particles");
  
	if(!prmArgProc(prm,m->argc,m->argv)) {
	    CkExit();
	    }
	
	// hardcoding some parameters, later may be full options
	useTree = Binary_Oct;

	//if (verbosity) 
	  ckerr << "yieldPeriod set to " << _yieldPeriod << endl;
	if(_cacheLineDepth < 0)
		CkAbort("Cache Line depth must be greater than or equal to 0");

        //if(verbosity)
          ckerr << "Prefetching..." << (_prefetch?"ON":"OFF") << endl;
  
          //if (verbosity)
	  ckerr << "Number of chunks for remote tree walk set to " << _numChunks << endl;
	if (_numChunks <= 0)
	  CkAbort("Number of chunks for remote tree walk must be greater than 0");

        //if(verbosity)
          ckerr << "Chunk Randomization..." << (_randChunks?"ON":"OFF") << endl;
  
	if(param.achInFile[0]) {
	    basefilename = param.achInFile;
	}else{
		ckerr<<"Base file name missing\n";
		CkExit();
	}

        if (_nocache) _cacheLineDepth = 0;

	//if(verbosity) {
	  ckerr<<"cache "<<_cache << (_nocache?" (disabled)":"") <<endl;
	  ckerr<<"cacheLineDepth "<<_cacheLineDepth<<endl;
        //}
	if (verbosity) {
	  if(printBinaryAcc==1)
	    ckerr<<"particle accelerations to be printed in binary format..."<<endl;
	  else
	    ckerr<<"particles accelerations to be printed in ASCII format..."<<endl;
        
	  ckerr << "Verbosity level " << verbosity << endl;
        }
        //if(verbosity) {
          switch(domainDecomposition){
          case SFC_dec:
            ckerr << "Domain decomposition...SFC" << endl;
            break;
          case Oct_dec:
            ckerr << "Domain decomposition...Oct" << endl;
            break;
          case ORB_dec:
            ckerr << "Domain decomposition...ORB" << endl;
            break;
          default:
            CkAbort("None of the implemented decompositions specified");
          }
        //}

	cacheManagerProxy = CProxy_CacheManager::ckNew(cacheSize);

	CProxy_BlockMap myMap=CProxy_BlockMap::ckNew(); 
	CkArrayOptions opts(numTreePieces); 
	opts.setMap(myMap);

	pieces = CProxy_TreePiece::ckNew(numTreePieces,opts);
	treeProxy = pieces;
	
	//create the DataManager
	dataManager = CProxy_DataManager::ckNew(pieces);
	dataManagerID = dataManager;

	streamingProxy = pieces;

	//create the Sorter
	sorter = CProxy_Sorter::ckNew();

	//StreamingStrategy* strategy = new StreamingStrategy(10,50);
	//ComlibAssociateProxy(strategy, streamingProxy);

	if(verbosity)
	  ckerr << "Created " << numTreePieces << " pieces of tree" << endl;
	
	CProxy_Main(thishandle).nextStage();
}

void Main::niceExit() {
  static unsigned int count = 0;
  if (++count == numTreePieces) CkExit();
}

void Main::nextStage() {
	double startTime;
	double tolerance = 0.01;	// tolerance for domain decomposition
	
	//piecedata *totaldata = new piecedata;
	
	// DEBUGGING
	//CkStartQD(CkCallback(CkIndex_TreePiece::quiescence(),pieces));

	//pieces.registerWithDataManager(dataManager, CkCallbackResumeThread());

	/******** Particles Loading ********/
	ckerr << "Loading particles ...";
	startTime = CkWallTimer();
	pieces.load(basefilename, CkCallbackResumeThread());
	if(!(pieces[0].ckLocal()->bLoaded)) {
	    // Try loading Tipsy format
	    ckerr << " trying Tipsy ..." << endl;
	    
	    Tipsy::PartialTipsyFile ptf(basefilename, 0, 1);
	    if(!ptf.loadedSuccessfully()) {
	      ckerr << endl << "Couldn't load the tipsy file \""
		    << basefilename.c_str()
		    << "\". Maybe it's not a tipsy file?" << endl;
	      CkExit();
	      return;
	    }
	    pieces.loadTipsy(basefilename, CkCallbackResumeThread());
	}	
	ckerr << " took " << (CkWallTimer() - startTime) << " seconds."
	      << endl;

	if(prmSpecified(prm,"dSoftening")) {
	    ckerr << "Set Softening...\n";
	    pieces.setSoft(param.dSoft);
	}
	
	for(int i=0; i<=param.nSteps; i++){
	  /******** Resorting of particles and Domain Decomposition ********/
	  ckerr << "Domain decomposition ...";
	  startTime = CkWallTimer();
	  sorter.startSorting(dataManager, numTreePieces, tolerance,
			      CkCallbackResumeThread());
	  ckerr << " took " << (CkWallTimer() - startTime) << " seconds."
		<< endl;

	  if (i > 0) {
	    /********* Load balancer ********/
	    ckerr << "Load balancer ...";
	    startTime = CkWallTimer();
	    pieces.startlb(CkCallbackResumeThread());
	    ckerr<< " took "<<(CkWallTimer() - startTime)
		 << " seconds." << endl;
	  }

	  /******** Tree Build *******/
	  ckerr << "Building trees ...";
	  startTime = CkWallTimer();
	  pieces.buildTree(bucketSize, CkCallbackResumeThread());
	  ckerr << " took " << (CkWallTimer() - startTime) << " seconds." << endl;

	  /******** Force Computation ********/
	  ckerr << "Calculating gravity (tree bucket, theta = " << theta << ") ...";
	  startTime = CkWallTimer();
	  //pieces.calculateGravityBucketTree(theta, CkCallbackResumeThread());
	  cacheManagerProxy.cacheSync(theta, CkCallbackResumeThread());
	  ckerr << " took " << (CkWallTimer() - startTime) << " seconds."
		<< endl;

#if COSMO_STATS > 0
	  /********* TreePiece Statistics ********/
	  ckerr << "Total statistics iteration " << i << ":" << endl;
	  CkCallback cb(CkCallback::resumeThread);
	  pieces.collectStatistics(cb);
	  CkReductionMsg *tps = (CkReductionMsg *) cb.thread_delay();
	  ((TreePieceStatistics*)tps->getData())->printTo(ckerr);

	  /********* Cache Statistics ********/
	  CkCallback ccb(CkCallback::resumeThread);
	  cacheManagerProxy.collectStatistics(ccb);
	  CkReductionMsg *cs = (CkReductionMsg *) ccb.thread_delay();
	  ((CacheStatistics*)cs->getData())->printTo(ckerr);
#endif

	  /******** Update of Particle Positions and Energy Estimation *******/
	  if(!param.bStaticTest && i > 0) // Closing Kick
	      pieces.kick(0.5*param.dTimeStep, CkCallbackResumeThread());
	  calcEnergy();
	  /*
	   * Writing of intermediate outputs can be done here.
	   */
	  if(!param.bStaticTest && i < param.nSteps) { // Opening Kick
						       // and Drift
	      pieces.kick(0.5*param.dTimeStep, CkCallbackResumeThread());
	      pieces.drift(param.dTimeStep, CkCallbackResumeThread());
	  }
	  
	} // End of main computation loop

	/******** Shutdown process ********/

	if(param.nSteps == 0) {
	    ckerr << "Outputting accelerations ...";
	    startTime = CkWallTimer();
	    if(printBinaryAcc)
		pieces[0].outputAccelerations(OrientedBox<double>(),
					      "acc2", CkCallbackResumeThread());
	    else
		pieces[0].outputAccASCII("acc2", CkCallbackResumeThread());
	
	    pieces[0].outputIOrderASCII("iord", CkCallbackResumeThread());
	    ckerr << " took " << (CkWallTimer() - startTime) << " seconds."
		  << endl;
	    }
	
#if COSMO_STATS > 0
	ckerr << "Outputting statistics ...";
	startTime = CkWallTimer();
	Interval<unsigned int> dummy;
	
	//pieces[0].outputStatistics(dummy, dummy, dummy, dummy, totaldata->totalmass, CkCallbackResumeThread());
	pieces[0].outputStatistics(dummy, dummy, dummy, dummy, 0, CkCallbackResumeThread());

	ckerr << " took " << (CkWallTimer() - startTime) << " seconds." << endl;
#endif

#if COSMO_DEBUG > 1
	ckerr << "Outputting relative errors ...";
	startTime = CkWallTimer();
	pieces[0].outputRelativeErrors(Interval<double>(), CkCallbackResumeThread());
	ckerr << " took " << (CkWallTimer() - startTime) << " seconds." << endl;
#endif

	ckerr << "Done." << endl;
	
	ckerr << endl << "******************" << endl << endl; 
	CkExit();
}

void Main::calcEnergy() {
    CkCallback ccb(CkCallback::resumeThread);
    pieces.calcEnergy(ccb);
    CkReductionMsg *msg = (CkReductionMsg *) ccb.thread_delay();
    double *dEnergy = (double *) msg->getData();
    CkPrintf("Energy: %g %g %g %g %g %g %g %g\n", dEnergy[0] + dEnergy[1],
	     dEnergy[0]+dEnergy[2], dEnergy[0],
	     dEnergy[1], dEnergy[2], dEnergy[3], dEnergy[4], dEnergy[5]);
    delete msg;
}

void registerStatistics() {
#if COSMO_STATS > 0
  CacheStatistics::sum = CkReduction::addReducer(CacheStatistics::sumFn);
  TreePieceStatistics::sum = CkReduction::addReducer(TreePieceStatistics::sumFn);
#endif
}

#include "ParallelGravity.def.h"
#include "CacheManager.def.h"

/*
#define CK_TEMPLATES_ONLY
#include "CacheManager.def.h"
#undef CK_TEMPLATES_ONLY
#include "CacheManager.def.h"
*/
