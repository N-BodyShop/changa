/** @file ParallelGravity.cpp
 */
 
#include <stdint.h>
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
CProxy_CacheManager streamingCache;
int _prefetch;
int _numChunks;
int _randChunks;

ComlibInstanceHandle cinst1, cinst2;

int boundaryEvaluationUE;
int weightBalanceUE;
int networkProgressUE;
int nodeForceUE;
int partForceUE;

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

int killAt;

Main::Main(CkArgMsg* m) {
	_cache = true;
	mainChare = thishandle;

        boundaryEvaluationUE = traceRegisterUserEvent("Evaluating Boudaries");
        weightBalanceUE = traceRegisterUserEvent("Weight Balancer");
        networkProgressUE = traceRegisterUserEvent("CmiNetworkProgress");
        nodeForceUE = traceRegisterUserEvent("Node interaction");
        partForceUE = traceRegisterUserEvent("Particle interaction");
	
	prmInitialize(&prm,_Leader,_Trailer);
	csmInitialize(&param.csm);

	param.dTimeStep = 0.0;
	prmAddParam(prm, "dTimeStep", paramDouble, &param.dTimeStep,
		    sizeof(double),"dT", "Base Timestep for integration");
	param.nSteps = 1;
	prmAddParam(prm, "nSteps", paramInt, &param.nSteps,
		    sizeof(int),"n", "Number of Timesteps");
	param.iStartStep = 0;
	prmAddParam(prm, "iStartStep", paramInt, &param.iStartStep,
		    sizeof(int),"nstart", "Initial step numbering");

        killAt = 0;
        prmAddParam(prm, "killAt", paramInt, &killAt,
		    sizeof(int),"killat", "Killing after this step");

	param.dEta = 0.03;
	prmAddParam(prm, "dEta", paramDouble, &param.dEta,
		    sizeof(double),"eta", "Time integration accuracy");
	param.iOutInterval = 10;
	prmAddParam(prm, "iOutInterval", paramInt, &param.iOutInterval,
		    sizeof(int),"oi", "Output Interval");
	param.iLogInterval = 1;
	prmAddParam(prm, "iLogInterval", paramInt, &param.iLogInterval,
		    sizeof(int),"log", "Log Interval");
	param.dSoft = 0.0;
	prmAddParam(prm, "dSoftening", paramDouble, &param.dSoft,
		    sizeof(double),"S", "Gravitational softening");
	theta = 0.7;
	prmAddParam(prm, "dTheta", paramDouble, &theta,
		    sizeof(double),"t", "Opening angle");
	param.bPeriodic = 0;
	prmAddParam(prm, "bPeriodic", paramBool, &param.bPeriodic,
		    sizeof(int),"per", "Periodic Boundaries");
	param.nReplicas = 1;
	prmAddParam(prm, "nReplicas", paramInt, &param.nReplicas,
		    sizeof(int),"nrep", "Number of periodic replicas");
	param.fPeriod = 1.0;
	prmAddParam(prm, "fPeriod", paramDouble, &param.fPeriod,
		    sizeof(double),"fper", "Periodic size");
	param.bEwald = 1;
	prmAddParam(prm,"bEwald",paramBool, &param.bEwald, sizeof(int),
		    "ewald", "enable/disable Ewald correction = +ewald");
	param.dEwCut = 2.6;
	prmAddParam(prm,"dEwCut", paramDouble, &param.dEwCut, sizeof(double),
		    "ewc", "<dEwCut> = 2.6");
	param.dEwhCut = 2.8;
	prmAddParam(prm,"dEwhCut", paramDouble, &param.dEwhCut, sizeof(double),
		    "ewh", "<dEwhCut> = 2.8");
	param.csm->bComove = 0;
	prmAddParam(prm, "bComove", paramBool, &param.csm->bComove,
		    sizeof(int),"cm", "Comoving coordinates");
	param.csm->dHubble0 = 0.0;
	prmAddParam(prm,"dHubble0",paramDouble,&param.csm->dHubble0, 
		    sizeof(double),"Hub", "<dHubble0> = 0.0");
	param.csm->dOmega0 = 1.0;
	prmAddParam(prm,"dOmega0",paramDouble,&param.csm->dOmega0,
		    sizeof(double),"Om", "<dOmega0> = 1.0");
	param.csm->dLambda = 0.0;
	prmAddParam(prm,"dLambda",paramDouble,&param.csm->dLambda,
		    sizeof(double),"Lambda", "<dLambda> = 0.0");
	param.csm->dOmegaRad = 0.0;
	prmAddParam(prm,"dOmegaRad",paramDouble,&param.csm->dOmegaRad,
		    sizeof(double),"Omrad", "<dOmegaRad> = 0.0");
	param.csm->dOmegab = 0.0;
	prmAddParam(prm,"dOmegab",paramDouble,&param.csm->dOmegab,
		    sizeof(double),"Omb", "<dOmegab> = 0.0");
	param.csm->dQuintess = 0.0;
	prmAddParam(prm,"dQuintess",paramDouble,&param.csm->dQuintess,
		    sizeof(double),"Quint",
		    "<dQuintessence (constant w = -1/2) > = 0.0");
	param.dRedTo = 0.0;
	prmAddParam(prm,"dRedTo",paramDouble,&param.dRedTo,sizeof(double),
		    "zto", "specifies final redshift for the simulation");
	

	printBinaryAcc=1;
	prmAddParam(prm, "bPrintBinary", paramBool, &printBinaryAcc,
		    sizeof(int),"z", "Print accelerations in Binary");
	param.achInFile[0] = '\0';
	prmAddParam(prm,"achInFile",paramString,param.achInFile,
		    256, "I", "input file name (or base file name)");
	strcpy(param.achOutName,"pargrav");
	prmAddParam(prm,"achOutName",3,param.achOutName,256,"o",
				"output name for snapshots and logfile");
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
	prmAddParam(prm, "nBucket", paramInt, &bucketSize,
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
	
	if(!param.bPeriodic) {
	    param.nReplicas = 0;
	    param.fPeriod = 1.0e38;
	    param.bEwald = 0;
	    }
	    
	// hardcoding some parameters, later may be full options
	if(domainDecomposition==ORB_dec){ useTree = Binary_ORB; }
  else { useTree = Binary_Oct; }

        ckerr << "Running of " << CkNumPes() << " processors with " << numTreePieces << " TreePieces" << endl;

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
        streamingCache = cacheManagerProxy;

	//create the Sorter
	sorter = CProxy_Sorter::ckNew();

	StreamingStrategy* strategy = new StreamingStrategy(10,50);
        cinst1 = ComlibRegister(strategy);
	//ComlibAssociateProxy(strategy, streamingProxy);
        strategy = new StreamingStrategy(10,50);
        cinst2 = ComlibRegister(strategy);
        //ComlibAssociateProxy(strategy, streamingCache);

	if(verbosity)
	  ckerr << "Created " << numTreePieces << " pieces of tree" << endl;
	
	CProxy_Main(thishandle).nextStage();
}

void Main::niceExit() {
  static unsigned int count = 0;
  if (++count == numTreePieces) CkExit();
}

/// Function to determine the start time of the simulation.
/// Modifies dTime member of main.
void Main::getStartTime() 
{
	// Grab starting time.
	if(param.csm->bComove) {
	    if(param.csm->dHubble0 == 0.0) {
		ckerr << "No hubble constant specified" << endl;
		CkExit();
		return;
		}
	    /* input time is expansion factor.  Convert */
	    double dAStart = pieces[0].ckLocal()->fh.time; 
	    double z = 1.0/dAStart - 1.0;
	    double aTo, tTo;
	    dTime = csmExp2Time(param.csm, dAStart);
	    if (verbosity > 0)
		ckout << "Input file, Time:" << dTime << "Redshift:" << z
		      << "Expansion factor:" << dAStart << endl;
	    if (prmSpecified(prm,"dRedTo")) {
		if (!prmArgSpecified(prm,"nSteps") &&
		    prmArgSpecified(prm,"dTimeStep")) {
		    aTo = 1.0/(param.dRedTo + 1.0);
		    tTo = csmExp2Time(param.csm,aTo);
		    if (verbosity > 0)
			ckout << "Simulation to Time:" << tTo << " Redshift:"
			      << 1.0/aTo - 1.0 << " Expansion factor:" << aTo
			      << endl;
		    if (tTo < dTime) {
			ckerr << "Badly specified final redshift, check -zto parameter."
			      << endl;
			CkExit();
			}
		    param.nSteps = (int)ceil((tTo-dTime)/param.dTimeStep);
		    }
		else if (!prmArgSpecified(prm,"dTimeStep") &&
			 prmArgSpecified(prm,"nSteps")) {
		    aTo = 1.0/(param.dRedTo + 1.0);
		    tTo = csmExp2Time(param.csm,aTo);
		    if (verbosity > 0)
			ckout << "Simulation to Time:" << tTo << " Redshift:"
			      << 1.0/aTo - 1.0 << " Expansion factor:" << aTo
			      << endl;
		    if (tTo < dTime) {
			ckerr << "Badly specified final redshift, check -zto parameter."
			      << endl;
			CkExit();
			}
		    if(param.nSteps != 0)
			param.dTimeStep =
				(tTo-dTime)/(param.nSteps - param.iStartStep);
		    else
			param.dTimeStep = 0.0;
		    }
		else if (!prmSpecified(prm,"nSteps") &&
			 prmFileSpecified(prm,"dTimeStep")) {
		    aTo = 1.0/(param.dRedTo + 1.0);
		    tTo = csmExp2Time(param.csm,aTo);
		    if (verbosity > 0)
			ckout << "Simulation to Time:" << tTo << " Redshift:"
			      << 1.0/aTo - 1.0 << " Expansion factor:" << aTo
			      << endl;
		    if (tTo < dTime) {
			ckerr << "Badly specified final redshift, check -zto parameter."
			      << endl;
			CkExit();
			}
		    param.nSteps = (int)ceil((tTo-dTime)/param.dTimeStep);
		    }
		else if (!prmSpecified(prm,"dTimeStep") &&
			 prmFileSpecified(prm,"nSteps")) {
		    aTo = 1.0/(param.dRedTo + 1.0);
		    tTo = csmExp2Time(param.csm,aTo);
		    if (verbosity > 0)
			ckout << "Simulation to Time:" << tTo << " Redshift:"
			      << 1.0/aTo - 1.0 << " Expansion factor:" << aTo
			      << endl;
		    if (tTo < dTime) {
			ckerr << "Badly specified final redshift, check -zto parameter."
			      << endl;
			CkExit();
			}
		    if(param.nSteps != 0)
			param.dTimeStep =	(tTo-dTime)/(param.nSteps
							 - param.iStartStep);
		    else
			param.dTimeStep = 0.0;
		    }
		}
	    else {
		tTo = dTime + param.nSteps*param.dTimeStep;
		aTo = csmTime2Exp(param.csm,tTo);
		if (verbosity > 0)
		    ckout << "Simulation to Time:" << tTo << " Redshift:"
			  << 1.0/aTo - 1.0 << " Expansion factor:"
			  << aTo << endl;
		}
	    // convert to canonical comoving coordinates.
	    pieces.velScale(dAStart*dAStart);
	    }
	else {  // not Comove
	    dTime = pieces[0].ckLocal()->fh.time; 
	    if(verbosity > 0)
		ckout << "Input file, Time:" << dTime << endl;
	    double tTo = dTime + param.nSteps*param.dTimeStep;
	    if (verbosity > 0) {
		ckout << "Simulation to Time:" << tTo << endl;
		}
	    }
    }

void Main::advanceBigStep(int iStep) {
  double tolerance = 0.01;	// tolerance for domain decomposition
  int currentStep = 0; // the current timestep within the big step
  int activeRung = 0;  // the minimum rung that is active
  int nextMaxRung = 0; // the rung that determines the smallest time for advancing

  while (currentStep < MAXSUBSTEPS) {

    if(!param.bStaticTest) {
      // Find new rung for active particles
      nextMaxRung = adjust(activeRung);
      if(currentStep == 0) rungStats();
      if(verbosity) ckerr << "MaxRung: " << nextMaxRung << endl;
      CkAssert(nextMaxRung <= MAXRUNG); // timesteps WAY too short
 
      // Opening Kick
      double dKickFac[MAXRUNG+1];
      for(int iRung = activeRung; iRung <= nextMaxRung; iRung++) {
        double dTimeSub = RungToDt(param.dTimeStep, iRung);
        dKickFac[iRung] = csmComoveKickFac(param.csm, dTime, 0.5*dTimeSub);
      }
      pieces.kick(activeRung, dKickFac, CkCallbackResumeThread());

      // Drift of smallest step
      double dTimeSub = RungToDt(param.dTimeStep, nextMaxRung);
      double dDriftFac = csmComoveDriftFac(param.csm, dTime, dTimeSub);
      pieces.drift(dDriftFac, CkCallbackResumeThread());

      // Advance time to end of smallest step
      dTime += dTimeSub;
    }
    currentStep += RungToSubsteps(nextMaxRung);

    // determine largest timestep that needs a kick
    activeRung = 0;
    int tmpRung = currentStep;
    while (tmpRung &= ~MAXSUBSTEPS) {
      activeRung++;
      tmpRung <<= 1;
    }

    if(verbosity)
      ckerr << "Step: " << iStep << " Substep: " << currentStep
            << " Time: " << dTime
            << " Gravity rung " << activeRung << " to "
            << nextMaxRung << endl;


    /***** Resorting of particles and Domain Decomposition *****/
    ckerr << "Domain decomposition ...";
    double startTime = CkWallTimer();
    sorter.startSorting(dataManager, numTreePieces, tolerance,
                        CkCallbackResumeThread());
    ckerr << " took " << (CkWallTimer() - startTime) << " seconds."
          << endl;

    /********* Load balancer ********/
    ckerr << "Load balancer ...";
    startTime = CkWallTimer();
    pieces.startlb(CkCallbackResumeThread());
    ckerr<< " took "<<(CkWallTimer() - startTime) << " seconds."
         << endl;

    /******** Tree Build *******/
    ckerr << "Building trees ...";
    startTime = CkWallTimer();
    pieces.buildTree(bucketSize, CkCallbackResumeThread());
    ckerr << " took " << (CkWallTimer() - startTime) << " seconds."
          << endl;

    /******** Force Computation ********/
    ckerr << "Calculating gravity (tree bucket, theta = " << theta
          << ") ...";
    startTime = CkWallTimer();
    // Set up Ewald Tables
    if(param.bPeriodic && param.bEwald)
      pieces.EwaldInit(param.dEwhCut);
    //pieces.calculateGravityBucketTree(theta, CkCallbackResumeThread());
    cacheManagerProxy.cacheSync(theta, activeRung, CkCallbackResumeThread());
    ckerr << " took " << (CkWallTimer() - startTime) << " seconds."
          << endl;

#if COSMO_STATS > 0
    /********* TreePiece Statistics ********/
    ckerr << "Total statistics iteration " << iStep << "/";
    int printingStep = currentStep;
    while (printingStep) {
      ckerr << ((printingStep & MAXSUBSTEPS) ? "1" : "0");
      printingStep = (printingStep & ~MAXSUBSTEPS) << 1;
    }
    ckerr << " (rungs " << activeRung << "-" << nextMaxRung << ")" << ":" << endl;
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



    if(!param.bStaticTest) {
      // Closing Kick
      double dKickFac[MAXRUNG+1];
      for(int iRung = activeRung; iRung <= nextMaxRung; iRung++) {
        double dTimeSub = RungToDt(param.dTimeStep, iRung);
        dKickFac[iRung] = csmComoveKickFac(param.csm,
                                           dTime - 0.5*dTimeSub,
                                           0.5*dTimeSub);
      }
      pieces.kick(activeRung, dKickFac, CkCallbackResumeThread());
    }
		
  }
}
    
void Main::nextStage() {
  double startTime;
  double tolerance = 0.01;	// tolerance for domain decomposition

  // DEBUGGING
  //CkStartQD(CkCallback(CkIndex_TreePiece::quiescence(),pieces));

  pieces.setPeriodic(param.nReplicas, param.fPeriod, param.bEwald,
                     param.dEwCut, param.bPeriodic);

  /******** Particles Loading ********/
  ckerr << "Loading particles ...";
  startTime = CkWallTimer();
  pieces.load(basefilename, CkCallbackResumeThread());
  if(!(pieces[0].ckLocal()->bLoaded)) {
    // Try loading Tipsy format
    ckerr << " trying Tipsy ...";
	    
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
	    
  getStartTime();
	
  char achLogFileName[MAXPATHLEN];
  sprintf(achLogFileName, "%s.log", param.achOutName);
  ofstream ofsLog(achLogFileName, ios_base::trunc);
  ofsLog << "# Starting ParallelGravity" << endl;
  ofsLog.close();
	
  prmLogParam(prm, achLogFileName);
	
  if(prmSpecified(prm,"dSoftening")) {
    ckerr << "Set Softening...\n";
    pieces.setSoft(param.dSoft);
  }
	
  if(param.bPeriodic) {	// puts all particles within the boundary
    pieces.drift(0.0, CkCallbackResumeThread());
  }
  
  /***** Initial sorting of particles and Domain Decomposition *****/
  ckerr << "Initial domain decomposition ...";
  startTime = CkWallTimer();
  sorter.startSorting(dataManager, numTreePieces, tolerance,
                      CkCallbackResumeThread());
  ckerr << " took " << (CkWallTimer() - startTime) << " seconds."
        << endl;
  
  /******** Tree Build *******/
  ckerr << "Building trees ...";
  startTime = CkWallTimer();
  pieces.buildTree(bucketSize, CkCallbackResumeThread());
  ckerr << " took " << (CkWallTimer() - startTime) << " seconds."
        << endl;
  
  /******** Force Computation ********/
  ckerr << "Calculating gravity (theta = " << theta
        << ") ...";
  startTime = CkWallTimer();
  // Set up Ewald Tables
  if(param.bPeriodic && param.bEwald)
    pieces.EwaldInit(param.dEwhCut);
  //pieces.calculateGravityBucketTree(theta, CkCallbackResumeThread());
  cacheManagerProxy.cacheSync(theta, 0, CkCallbackResumeThread());
  ckerr << " took " << (CkWallTimer() - startTime) << " seconds."
        << endl;
  
#if COSMO_STATS > 0
  /********* TreePiece Statistics ********/
  ckerr << "Total statistics initial iteration :" << endl;
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
  
  for(int iStep = param.iStartStep+1; iStep <= param.nSteps; iStep++){

    if (verbosity) ckerr << "Starting big step " << iStep << endl;
    startTime = CkWallTimer();
    advanceBigStep(iStep-1);
    ckerr << "Big step " << iStep << " took " << (CkWallTimer() - startTime) << "seconds." <<endl;

    if(iStep%param.iLogInterval == 0) {
      calcEnergy(dTime, achLogFileName);
    }
    /*
     * Writing of intermediate outputs can be done here.
     */
    if(iStep%param.iOutInterval == 0) {
      writeOutput(iStep);
    }
	  
    if (killAt > 0 && killAt == iStep) break;

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

  /*#if COSMO_DEBUG > 1
    ckerr << "Outputting relative errors ...";
    startTime = CkWallTimer();
    pieces[0].outputRelativeErrors(Interval<double>(), CkCallbackResumeThread());
    ckerr << " took " << (CkWallTimer() - startTime) << " seconds." << endl;
    #endif*/

  ckerr << "Done." << endl;
	
  ckerr << endl << "******************" << endl << endl; 
  CkExit();
}

void Main::calcEnergy(double dTime, char *achLogFileName) {
    CkCallback ccb(CkCallback::resumeThread);
    pieces.calcEnergy(ccb);
    CkReductionMsg *msg = (CkReductionMsg *) ccb.thread_delay();
    double *dEnergy = (double *) msg->getData();

    FILE *fpLog = fopen(achLogFileName, "a");
    
    static int first = 1;
    
    double a = csmTime2Exp(param.csm, dTime);
    
    if(first) {
	fprintf(fpLog, "# time redshift TotalEVir TotalE Kinetic Virial Potential TotalECosmo Lx Ly Lz\n");
	first = 0;
	dEcosmo = 0.0;
	}
    else {
	if(param.csm->bComove) {
	    dEcosmo += 0.5*(a - csmTime2Exp(param.csm, dTimeOld))
		*(dEnergy[2] + dUOld);
	    }
	}

    dUOld = dEnergy[2];
    dTimeOld = dTime;
    double z = 1.0/a - 1.0;
    
    fprintf(fpLog, "%g %g %g %g %g %g %g %g %g %g %g\n", dTime, z,
	    dEnergy[0] + dEnergy[1], dEnergy[0]+dEnergy[2], dEnergy[0],
	    dEnergy[1], dEnergy[2], dEcosmo, dEnergy[3], dEnergy[4],
	    dEnergy[5]);
    fclose(fpLog);
    
    delete msg;
}

void Main::writeOutput(int iStep) 
{
    char achFile[256];
    double dOutTime;
    double dvFac;
    
    sprintf(achFile,"%s.%06i",param.achOutName,iStep);
    if(param.csm->bComove) {
	dOutTime = csmTime2Exp(param.csm, dTime);
	dvFac = 1.0/(dOutTime*dOutTime);
	}
    else {
	dOutTime = dTime;
	dvFac = 1.0;
	}
    
    pieces.setupWrite(0, 0, achFile, dOutTime, dvFac,
		      CkCallbackResumeThread());
    }

int Main::adjust(int iKickRung) 
{
    CkCallback ccb(CkCallback::resumeThread);
    double a = csmTime2Exp(param.csm,dTime);
    
    pieces.adjust(iKickRung, param.dEta, param.dTimeStep, 1.0/(a*a*a), ccb);

    CkReductionMsg *msg = (CkReductionMsg *) ccb.thread_delay();
    int iCurrMaxRung = *(int *)msg->getData();
    delete msg;
    return iCurrMaxRung;
    }

void Main::rungStats() 
{
    CkCallback ccb(CkCallback::resumeThread);
    pieces.rungStats(ccb);
    CkReductionMsg *msg = (CkReductionMsg *) ccb.thread_delay();
    int *nInRung = (int *)msg->getData();
    
    ckout << "Rung distribution: (";
    for(int iRung = 0; iRung <= MAXRUNG; iRung++)
	ckout << " " << nInRung[iRung] << ",";
    ckout << ")" << endl;
    
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
