/** @file ParallelGravity.cpp
 */
 
#include <stdint.h>
#include <iostream>

#include <unistd.h>
#include <sys/param.h>

// Debug floating point problems
// #include <fenv.h>

#include "Sorter.h"
#include "ParallelGravity.h"
#include "DataManager.h"
#include "TipsyFile.h"
#include "param.h"

extern char *optarg;
extern int optind, opterr, optopt;

using namespace std;

CProxy_Main mainChare;
int verbosity;
CProxy_TreePiece treeProxy;
CProxy_CkCacheManager cacheManagerProxy;
//CkReduction::reducerType callbackReduction;
bool _cache;
int _nocache;
int _cacheLineDepth;
unsigned int _yieldPeriod;
DomainsDec domainDecomposition;
int peanoKey;
GenericTrees useTree;
CProxy_TreePiece streamingProxy;
CProxy_CkCacheManager streamingCache;
unsigned int numTreePieces;
int _prefetch;
int _numChunks;
int _randChunks;
unsigned int bucketSize;
int lbcomm_cutoff_msgs;

/* readonly */ double theta;
/* readonly */ double thetaMono;
/* readonly */ int nSmooth;

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
	bIsRestarting = 0;
	bChkFirst = 1;
	dSimStartTime = CkWallTimer();

	// Floating point exceptions.
	// feenableexcept(FE_OVERFLOW | FE_DIVBYZERO | FE_INVALID);

        boundaryEvaluationUE = traceRegisterUserEvent("Evaluating Boudaries");
        weightBalanceUE = traceRegisterUserEvent("Weight Balancer");
        networkProgressUE = traceRegisterUserEvent("CmiNetworkProgress");
        nodeForceUE = traceRegisterUserEvent("Node interaction");
        partForceUE = traceRegisterUserEvent("Particle interaction");
	
	prmInitialize(&prm,_Leader,_Trailer);
	csmInitialize(&param.csm);

	param.dDelta = 0.0;
	prmAddParam(prm, "dDelta", paramDouble, &param.dDelta,
		    sizeof(double),"dt", "Base Timestep for integration");
	param.iMaxRung = MAXRUNG;
	prmAddParam(prm, "iMaxRung", paramInt, &param.iMaxRung,
		    sizeof(int),"mrung", "Maximum timestep rung (IGNORED)");
	param.nSteps = 1;
	prmAddParam(prm, "nSteps", paramInt, &param.nSteps,
		    sizeof(int),"n", "Number of Timesteps");
	param.iStartStep = 0;
	prmAddParam(prm, "iStartStep", paramInt, &param.iStartStep,
		    sizeof(int),"nstart", "Initial step numbering");
	param.iWallRunTime = 0;
	prmAddParam(prm,"iWallRunTime",paramInt,&param.iWallRunTime,
		    sizeof(int),"wall",
		    "<Maximum Wallclock time (in minutes) to run> = 0 = infinite");

        killAt = 0;
        prmAddParam(prm, "killAt", paramInt, &killAt,
		    sizeof(int),"killat", "Killing after this step");

	param.bEpsAccStep = 1;
	prmAddParam(prm, "bEpsAccStep", paramBool, &param.bEpsAccStep,
		    sizeof(int),"epsacc", "Use sqrt(eps/a) timestepping");
	param.bGravStep = 0;
	prmAddParam(prm, "bGravStep", paramBool, &param.bGravStep,
		    sizeof(int),"gravstep",
		    "Use gravity interaction timestepping");
	param.dEta = 0.03;
	prmAddParam(prm, "dEta", paramDouble, &param.dEta,
		    sizeof(double),"eta", "Time integration accuracy");
	param.bCannonical = 1;
	prmAddParam(prm, "bCannonical", paramBool, &param.bCannonical,
		    sizeof(int),"can", "Cannonical Comoving eqns (IGNORED)");
	param.bKDK = 1;
	prmAddParam(prm, "bKDK", paramBool, &param.bKDK,
		    sizeof(int),"kdk", "KDK timestepping (IGNORED)");
	param.bBenchmark = 0;
	prmAddParam(prm, "bBenchmark", paramBool, &param.bBenchmark,
		    sizeof(int),"bench", "Benchmark only; no output or checkpoints");
	param.iOutInterval = 10;
	prmAddParam(prm, "iOutInterval", paramInt, &param.iOutInterval,
		    sizeof(int),"oi", "Output Interval");
	param.iLogInterval = 1;
	prmAddParam(prm, "iLogInterval", paramInt, &param.iLogInterval,
		    sizeof(int),"ol", "Log Interval");
	param.iCheckInterval = 10;
	prmAddParam(prm, "iCheckInterval", paramInt, &param.iCheckInterval,
		    sizeof(int),"oc", "Checkpoint Interval");
	param.bDoDensity = 1;
	prmAddParam(prm, "bDoDensity", paramBool, &param.bDoDensity,
		    sizeof(int),"den", "Enable Density outputs");
	nSmooth = 32;
	prmAddParam(prm, "nSmooth", paramInt, &nSmooth,
		    sizeof(int),"nsm", "Number of neighbors for smooth");
	param.bDoGravity = 1;
	prmAddParam(prm, "bDoGravity", paramBool, &param.bDoGravity,
		    sizeof(int),"g", "Enable Gravity");
	param.dSoft = 0.0;
	prmAddParam(prm, "dSoft", paramDouble, &param.dSoft,
		    sizeof(double),"e", "Gravitational softening");
	theta = 0.7;
	prmAddParam(prm, "dTheta", paramDouble, &theta,
		    sizeof(double), "theta", "Opening angle");
	param.dTheta2 = 0.7;
	prmAddParam(prm, "dTheta2", paramDouble, &param.dTheta2,
		    sizeof(double),"theta2",
		    "Opening angle after switchTheta (IGNORED)");
#ifdef HEXADECAPOLE
	param.iOrder = 4;
#else
	param.iOrder = 2;
#endif
	prmAddParam(prm, "iOrder", paramInt, &param.iOrder,
		    sizeof(int), "or", "Multipole expansion order(IGNORED)");
	param.bPeriodic = 0;
	prmAddParam(prm, "bPeriodic", paramBool, &param.bPeriodic,
		    sizeof(int),"per", "Periodic Boundaries");
	param.nReplicas = 1;
	prmAddParam(prm, "nReplicas", paramInt, &param.nReplicas,
		    sizeof(int),"nrep", "Number of periodic replicas");
	param.fPeriod = 1.0;
	prmAddParam(prm, "dPeriod", paramDouble, &param.fPeriod,
		    sizeof(double),"L", "Periodic size");
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
	
	param.dMsolUnit = 1.0;
	prmAddParam(prm,"dMsolUnit",paramDouble,&param.dMsolUnit,
		    sizeof(double),"msu", "<Solar mass/system mass unit>");
	param.dKpcUnit = 1000.0;
	prmAddParam(prm,"dKpcUnit",paramDouble,&param.dKpcUnit,
		    sizeof(double),"kpcu", "<Kiloparsec/system length unit>");

	printBinaryAcc=1;
	prmAddParam(prm, "bPrintBinary", paramBool, &printBinaryAcc,
		    sizeof(int),"z", "Print accelerations in Binary");
	param.bStandard = 1;
	prmAddParam(prm, "bStandard", paramBool, &param.bStandard,sizeof(int),
		    "std", "output in standard TIPSY binary format (IGNORED)");
	param.bOverwrite = 1;
	prmAddParam(prm, "bOverwrite", paramBool, &param.bOverwrite,sizeof(int),
		    "overwrite", "overwrite outputs (IGNORED)");
	param.bParaRead = 1;
	prmAddParam(prm, "bParaRead", paramBool, &param.bParaRead,sizeof(int),
		    "par", "enable/disable parallel reading of files (IGNORED)");
	param.bParaWrite = 0;
	prmAddParam(prm, "bParaWrite", paramBool, &param.bParaWrite,sizeof(int),
		    "paw", "enable/disable parallel writing of files (IGNORED)");
	param.achInFile[0] = '\0';
	prmAddParam(prm,"achInFile",paramString,param.achInFile,
		    256, "I", "input file name (or base file name)");
	strcpy(param.achOutName,"pargrav");
	prmAddParam(prm,"achOutName",paramString,param.achOutName,256,"o",
				"output name for snapshots and logfile");
	param.dExtraStore = 0;
	prmAddParam(prm,"dExtraStore",paramDouble,&param.dExtraStore,
		    sizeof(double), NULL, "IGNORED");
	
	bDumpFrame = 0;
	df = NULL;
	param.dDumpFrameStep = -1.0;
	prmAddParam(prm,"dDumpFrameStep",paramDouble, &param.dDumpFrameStep,
		    sizeof(double), "dfi",
		    "<number of steps between dumped frames> = -1 (disabled)");
	param.dDumpFrameTime = -1.0;
	prmAddParam(prm,"dDumpFrameTime",paramDouble,&param.dDumpFrameTime,
		    sizeof(double), "dft",
		    "<time interval between dumped frames> = -1 (disabled)");
	param.iDirector = 1;
	prmAddParam(prm,"iDirector",paramInt,&param.iDirector,sizeof(int),
		    "idr","<number of director files: 1, 2, 3> = 1");

	
	param.bStaticTest = 0;
	prmAddParam(prm, "bStaticTest", paramBool, &param.bStaticTest,
		    sizeof(int),"st", "Static test of performance");
	
	_randChunks = 1;
	prmAddParam(prm, "bRandChunks", paramBool, &_randChunks,
		    sizeof(int),"rand", "Randomize the order of remote chunk computation (default: ON)");
	
	_nocache = 0;
	prmAddParam(prm, "bNoCache", paramBool, &_nocache,
		    sizeof(int),"nc", "Disable the CacheManager caching behaviour");
	
        verbosity = 0;
	prmAddParam(prm, "iVerbosity", paramInt, &verbosity,
		    sizeof(int),"v", "Verbosity");
	numTreePieces = 8 * CkNumPes();
	prmAddParam(prm, "nTreePieces", paramInt, &numTreePieces,
		    sizeof(int),"p", "Number of TreePieces (default: 8*procs)");
	bucketSize = 12;
	prmAddParam(prm, "nBucket", paramInt, &bucketSize,
		    sizeof(int),"b", "Particles per Bucket (default: 12)");
	_numChunks = 1;
	prmAddParam(prm, "nChunks", paramInt, &_numChunks,
		    sizeof(int),"c", "Chunks per TreePiece (default: 1)");
	_yieldPeriod=5;
	prmAddParam(prm, "nYield", paramInt, &_yieldPeriod,
		    sizeof(int),"y", "Yield Period (default: 5)");
	_cacheLineDepth=4;
	prmAddParam(prm, "nCacheDepth", paramInt, &_cacheLineDepth,
		    sizeof(int),"d", "Cache Line Depth (default: 4)");
	_prefetch=true;
	prmAddParam(prm, "bPrefetch", paramBool, &_prefetch,
		    sizeof(int),"f", "Enable prefetching in the cache (default: ON)");
	int cacheSize = 100000000;
	prmAddParam(prm, "nCacheSize", paramInt, &cacheSize,
		    sizeof(int),"s", "Size of cache");
	domainDecomposition=SFC_dec;
        peanoKey=0;
	prmAddParam(prm, "nDomainDecompose", paramInt, &domainDecomposition,
		    sizeof(int),"D", "Kind of domain decomposition of particles");
        lbcomm_cutoff_msgs = 1;
	prmAddParam(prm, "lbcommCutoffMsgs", paramInt, &lbcomm_cutoff_msgs,
		    sizeof(int),"lbcommcut", "Cutoff for communication recording");
    
	
	if(!prmArgProc(prm,m->argc,m->argv)) {
	    CkExit();
	    }
	
	thetaMono = theta*theta*theta*theta;
	if(prmSpecified(prm, "iMaxRung")) {
	    ckerr << "WARNING: ";
	    ckerr << "iMaxRung parameter ignored. MaxRung is " << MAXRUNG
		  << endl;
	    }
	if(prmSpecified(prm, "bKDK")) {
	    ckerr << "WARNING: ";
	    ckerr << "bKDK parameter ignored; KDK is always used." << endl;
	    }
	if(prmSpecified(prm, "bStandard")) {
	    ckerr << "WARNING: ";
	    ckerr << "bStandard parameter ignored; Output is always standard."
		  << endl;
	    }
	if(prmSpecified(prm, "iOrder")) {
	    ckerr << "WARNING: ";
#ifdef HEXADECAPOLE
	    ckerr << "iOrder parameter ignored; expansion order is 4."
#else
	    ckerr << "iOrder parameter ignored; expansion order is 2."
#endif
		  << endl;
	    }
	if(prmSpecified(prm, "dTheta2")) {
	    ckerr << "WARNING: ";
	    ckerr << "dTheta2 parameter ignored."
		  << endl;
	    }
	if(prmSpecified(prm, "dExtraStore")) {
	    ckerr << "WARNING: ";
	    ckerr << "dExtraStore parameter ignored."
		  << endl;
	    }
	if(prmSpecified(prm, "bCannonical")) {
	    ckerr << "WARNING: ";
	    ckerr << "bCannonical parameter ignored."
		  << endl;
	    }
	if(prmSpecified(prm, "bOverwrite")) {
	    ckerr << "WARNING: ";
	    ckerr << "bOverwrite parameter ignored."
		  << endl;
	    }
	    
	if(!param.bPeriodic) {
	    param.nReplicas = 0;
	    param.fPeriod = 1.0e38;
	    param.bEwald = 0;
	    }
	    
	/* bolzman constant in cgs */
#define KBOLTZ	1.38e-16
	/* mass of hydrogen atom in grams */
#define MHYDR 1.67e-24
	/* solar mass in grams */
#define MSOLG 1.99e33
	/* G in cgs */
#define GCGS 6.67e-8
	/* kiloparsec in centimeters */
#define KPCCM 3.085678e21
	/* Thompson cross-section (cm^2) */
#define SIGMAT 6.6524e-25
	/* Speed of Light cm/s */
#define LIGHTSPEED 2.9979e10

#define SECONDSPERYEAR   31557600.
	/*
	 ** Convert kboltz/mhydrogen to system units, assuming that
	 ** G == 1.
	 */
	if(prmSpecified(prm, "dMsolUnit") &&
	   prmSpecified(prm, "dKpcUnit")) {
		param.dGasConst = param.dKpcUnit*KPCCM*KBOLTZ
			/MHYDR/GCGS/param.dMsolUnit/MSOLG;
		/* code energy per unit mass --> erg per g */
		param.dErgPerGmUnit = GCGS*param.dMsolUnit*MSOLG
		    /(param.dKpcUnit*KPCCM);
		/* code density --> g per cc */
		param.dGmPerCcUnit = (param.dMsolUnit*MSOLG)
		    /pow(param.dKpcUnit*KPCCM,3.0);
		/* code time --> seconds */
		param.dSecUnit = sqrt(1/(param.dGmPerCcUnit*GCGS));
		/* code comoving density --> g per cc = param.dGmPerCcUnit (1+z)^3 */
		param.dComovingGmPerCcUnit = param.dGmPerCcUnit;
		}

        if (domainDecomposition == SFC_peano_dec) peanoKey = 1;

	// hardcoding some parameters, later may be full options
	if(domainDecomposition==ORB_dec){ useTree = Binary_ORB; }
	else { useTree = Binary_Oct; }

        ckerr << "Running on " << CkNumPes() << " processors with " << numTreePieces << " TreePieces" << endl;

	if (verbosity) 
	  ckerr << "yieldPeriod set to " << _yieldPeriod << endl;
	if(_cacheLineDepth < 0)
		CkAbort("Cache Line depth must be greater than or equal to 0");

        if(verbosity)
          ckerr << "Prefetching..." << (_prefetch?"ON":"OFF") << endl;
  
        if (verbosity)
	  ckerr << "Number of chunks for remote tree walk set to " << _numChunks << endl;
	if (_numChunks <= 0)
	  CkAbort("Number of chunks for remote tree walk must be greater than 0");

        if(verbosity)
          ckerr << "Chunk Randomization..." << (_randChunks?"ON":"OFF") << endl;
  
	if(param.achInFile[0]) {
	    basefilename = param.achInFile;
	}else{
		ckerr<<"Base file name missing\n";
		CkExit();
	}

        if (_nocache) _cacheLineDepth = 0;

	if(verbosity) {
	  ckerr<<"cache "<<_cache << (_nocache?" (disabled)":"") <<endl;
	  ckerr<<"cacheLineDepth "<<_cacheLineDepth<<endl;
        }
	if (verbosity) {
	  if(printBinaryAcc==1)
	    ckerr<<"particle accelerations to be printed in binary format..."<<endl;
	  else
	    ckerr<<"particles accelerations to be printed in ASCII format..."<<endl;
        
	  ckerr << "Verbosity level " << verbosity << endl;
        }
        if(verbosity) {
          switch(domainDecomposition){
          case SFC_dec:
            ckerr << "Domain decomposition...SFC Morton" << endl;
            break;
          case Oct_dec:
            ckerr << "Domain decomposition...Oct" << endl;
            break;
          case ORB_dec:
            ckerr << "Domain decomposition...ORB" << endl;
            break;
          case SFC_peano_dec:
            ckerr << "Domain decomposition...SFC Peano-Hilbert" << endl;
            break;
          default:
            CkAbort("None of the implemented decompositions specified");
          }
        }

	CProxy_BlockMap myMap=CProxy_BlockMap::ckNew(); 
	CkArrayOptions opts(numTreePieces); 
	opts.setMap(myMap);

	CProxy_TreePiece pieces = CProxy_TreePiece::ckNew(numTreePieces,opts);
	treeProxy = pieces;

	// create the CacheManager
	cacheManagerProxy = CProxy_CkCacheManager::ckNew(cacheSize, pieces.ckLocMgr()->getGroupID());

	//create the DataManager
	CProxy_DataManager dataManager = CProxy_DataManager::ckNew(pieces);
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
	
	CProxy_Main(thishandle).setupICs();
}

Main::Main(CkMigrateMessage* m) {
    mainChare = thishandle;
    bIsRestarting = 1;
    sorter = CProxy_Sorter::ckNew();
    CkPrintf("Main(CkMigrateMessage) called\n");
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
	    double dAStart = treeProxy[0].ckLocal()->fh.time; 
	    double z = 1.0/dAStart - 1.0;
	    double aTo, tTo;
	    dTime = csmExp2Time(param.csm, dAStart);
	    if (verbosity > 0)
		ckout << "Input file, Time:" << dTime << " Redshift:" << z
		      << " Expansion factor:" << dAStart << endl;
	    if (prmSpecified(prm,"dRedTo")) {
		if (!prmArgSpecified(prm,"nSteps") &&
		    prmArgSpecified(prm,"dDelta")) {
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
		    param.nSteps = (int)ceil((tTo-dTime)/param.dDelta);
		    }
		else if (!prmArgSpecified(prm,"dDelta") &&
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
			param.dDelta =
				(tTo-dTime)/(param.nSteps - param.iStartStep);
		    else
			param.dDelta = 0.0;
		    }
		else if (!prmSpecified(prm,"nSteps") &&
			 prmFileSpecified(prm,"dDelta")) {
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
		    param.nSteps = (int)ceil((tTo-dTime)/param.dDelta);
		    }
		else if (!prmSpecified(prm,"dDelta") &&
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
			param.dDelta =	(tTo-dTime)/(param.nSteps
							 - param.iStartStep);
		    else
			param.dDelta = 0.0;
		    }
		}
	    else {
		tTo = dTime + param.nSteps*param.dDelta;
		aTo = csmTime2Exp(param.csm,tTo);
		if (verbosity > 0)
		    ckout << "Simulation to Time:" << tTo << " Redshift:"
			  << 1.0/aTo - 1.0 << " Expansion factor:"
			  << aTo << endl;
		}
	    // convert to canonical comoving coordinates.
	    treeProxy.velScale(dAStart*dAStart);
	    }
	else {  // not Comove
	    dTime = treeProxy[0].ckLocal()->fh.time; 
	    if(verbosity > 0)
		ckout << "Input file, Time:" << dTime << endl;
	    double tTo = dTime + param.nSteps*param.dDelta;
	    if (verbosity > 0) {
		ckout << "Simulation to Time:" << tTo << endl;
		}
	    }
    }

// determine if we need a smaller step for dumping frames
inline int Main::nextMaxRungIncDF(int nextMaxRung) 
{
    if (bDumpFrame && nextMaxRung < df[0]->iMaxRung)
	nextMaxRung = df[0]->iMaxRung;
    return nextMaxRung;
}

void Main::advanceBigStep(int iStep) {
  double tolerance = 0.01;	// tolerance for domain decomposition
  int currentStep = 0; // the current timestep within the big step
  int activeRung = 0;  // the minimum rung that is active
  int nextMaxRung = 0; // the rung that determines the smallest time for advancing

  while (currentStep < MAXSUBSTEPS) {

    if(!param.bStaticTest) {
      CkAssert(param.dDelta != 0.0);
      // Find new rung for active particles
      nextMaxRung = adjust(activeRung);
      if(currentStep == 0) rungStats();
      if(verbosity) ckerr << "MaxRung: " << nextMaxRung << endl;
      CkAssert(nextMaxRung <= MAXRUNG); // timesteps WAY too short
 
      // Opening Kick
      if(verbosity)
	  ckerr << "Kick Open:" << endl;
      double dKickFac[MAXRUNG+1];
      for(int iRung = activeRung; iRung <= nextMaxRung; iRung++) {
        double dTimeSub = RungToDt(param.dDelta, iRung);
	if(verbosity) {
	    ckerr << " Rung " << iRung << ": " << 0.5*dTimeSub << endl;
	    }
        dKickFac[iRung] = csmComoveKickFac(param.csm, dTime, 0.5*dTimeSub);
      }
      treeProxy.kick(activeRung, dKickFac, CkCallbackResumeThread());

      // Dump frame may require a smaller step
      int driftRung = nextMaxRungIncDF(nextMaxRung);
      int driftSteps = 1;
      while(driftRung > nextMaxRung) {
	  driftRung--;
	  driftSteps >>= 1;
	  }
      driftRung = nextMaxRungIncDF(nextMaxRung);
      
      double dTimeSub = RungToDt(param.dDelta, driftRung);
      // Drift of smallest step
      for(int iSub = 0; iSub < driftSteps; iSub++) 
	  {
	      if(verbosity)
		  ckerr << "Drift: Rung " << driftRung << " Delta "
			<< dTimeSub << endl;
	      double dDriftFac = csmComoveDriftFac(param.csm, dTime, dTimeSub);
	      treeProxy.drift(dDriftFac, CkCallbackResumeThread());

	      // Advance time to end of smallest step
	      dTime += dTimeSub;
	      currentStep += RungToSubsteps(driftRung);
	      /* 
	       ** Dump Frame
	       */
	      double dStep = iStep + ((double) currentStep)/MAXSUBSTEPS;
	      if (param.dDumpFrameTime > 0 && dTime >= df[0]->dTime)
		  DumpFrame(dTime, dStep );
	      else if(param.dDumpFrameStep > 0 && dStep >= df[0]->dStep) 
		  DumpFrame(dTime, dStep );
	      }
    }

    int lastActiveRung = activeRung;

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
    sorter.startSorting(dataManagerID, numTreePieces, tolerance,
                        CkCallbackResumeThread(), true);
    ckerr << " took " << (CkWallTimer() - startTime) << " seconds."
          << endl;

    /********* Load balancer ********/
    // jetley - commenting out lastActiveRung == 0 check, balance load even for fast rungs
    //if(lastActiveRung == 0) {
	ckerr << "Load balancer ...";
	startTime = CkWallTimer();
	treeProxy.startlb(CkCallbackResumeThread(), activeRung);
	ckerr<< " took "<<(CkWallTimer() - startTime) << " seconds."
	     << endl;
    //	}

    /******** Tree Build *******/
    ckerr << "Building trees ...";
    startTime = CkWallTimer();
    treeProxy.buildTree(bucketSize, CkCallbackResumeThread());
    ckerr << " took " << (CkWallTimer() - startTime) << " seconds."
          << endl;

    if(param.bDoGravity) {
	/******** Force Computation ********/
	ckerr << "Calculating gravity (tree bucket, theta = " << theta
	      << ") ...";
	startTime = CkWallTimer();
	treeProxy.startIteration(activeRung, CkCallbackResumeThread());
	ckerr << " took " << (CkWallTimer() - startTime) << " seconds."
	      << endl;
    }

#if COSMO_STATS > 0
    /********* TreePiece Statistics ********/
    ckerr << "Total statistics iteration " << iStep << "/";
    int printingStep = currentStep;
    while (printingStep) {
      ckerr << ((printingStep & MAXSUBSTEPS) ? "1" : "0");
      printingStep = (printingStep & ~MAXSUBSTEPS) << 1;
    }
    ckerr << " (rungs " << activeRung << "-" << nextMaxRung << ")" << ":" << endl;
    CkReductionMsg *tps;
    treeProxy.collectStatistics(CkCallbackResumeThread((void*&)tps));
    ((TreePieceStatistics*)tps->getData())->printTo(ckerr);

    /********* Cache Statistics ********/
    CkReductionMsg *cs;
    streamingCache.collectStatistics(CkCallbackResumeThread((void*&)cs));
    ((CacheStatistics*)cs->getData())->printTo(ckerr);
#endif



    if(!param.bStaticTest) {
      // Closing Kick
      double dKickFac[MAXRUNG+1];
      if(verbosity)
	  ckerr << "Kick Close:" << endl;
      for(int iRung = activeRung; iRung <= nextMaxRung; iRung++) {
        double dTimeSub = RungToDt(param.dDelta, iRung);
	if(verbosity) {
	    ckerr << " Rung " << iRung << ": " << 0.5*dTimeSub << endl;
	    }
        dKickFac[iRung] = csmComoveKickFac(param.csm,
                                           dTime - 0.5*dTimeSub,
                                           0.5*dTimeSub);
      }
      treeProxy.kick(activeRung, dKickFac, CkCallbackResumeThread());
    }
		
  }
}
    
// Load particles into pieces

void Main::setupICs() {
  double startTime;

  // DEBUGGING
  //CkStartQD(CkCallback(CkIndex_TreePiece::quiescence(),treeProxy));

  treeProxy.setPeriodic(param.nReplicas, param.fPeriod, param.bEwald,
			param.dEwCut, param.dEwhCut, param.bPeriodic);

  /******** Particles Loading ********/
  ckerr << "Loading particles ...";
  startTime = CkWallTimer();
  treeProxy.load(basefilename, CkCallbackResumeThread());
  if(!(treeProxy[0].ckLocal()->bLoaded)) {
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
    treeProxy.loadTipsy(basefilename, CkCallbackResumeThread());
    ckerr << treeProxy[0].ckLocal()->tipsyHeader << endl;
  }	
  ckerr << " took " << (CkWallTimer() - startTime) << " seconds."
        << endl;
	    
  getStartTime();
	
  char achLogFileName[MAXPATHLEN];
  sprintf(achLogFileName, "%s.log", param.achOutName);
  ofstream ofsLog(achLogFileName, ios_base::trunc);
  ofsLog << "# Starting ChaNGa" << endl;
  ofsLog.close();
	
  prmLogParam(prm, achLogFileName);
	
  if(prmSpecified(prm,"dSoft")) {
    ckerr << "Set Softening...\n";
    treeProxy.setSoft(param.dSoft);
  }
	
  if(param.bPeriodic) {	// puts all particles within the boundary
    treeProxy.drift(0.0, CkCallbackResumeThread());
  }
  
  initialForces();
}

int CheckForStop()
{
	/*
	 ** Checks for existence of STOPFILE in run directory. If
	 ** found, the file is removed and the return status is set
	 ** to 1, otherwise 0. 
	 */

	FILE *fp = NULL;
	const char *achFile = "STOP";

	if ((fp = fopen(achFile,"r")) != NULL) {
		(void) printf("User interrupt detected.\n");
		(void) fclose(fp);
		(void) unlink(achFile);
		return 1;
		}
	return 0;
	}

// Callback to restart simulation after a checkpoint or after writing
// a checkpoint.
void
Main::restart() 
{
    if(bIsRestarting) {
	dSimStartTime = CkWallTimer();
	ckout << "Restarting at " << param.iStartStep << endl;
	char achLogFileName[MAXPATHLEN];
	sprintf(achLogFileName, "%s.log", param.achOutName);
	ofstream ofsLog(achLogFileName, ios_base::app);
	ofsLog << "# ReStarting ChaNGa" << endl;
	ofsLog.close();
	treeProxy.drift(0.0, CkCallbackResumeThread());
	mainChare.initialForces();
	}
    else {
	ofstream ofsCheck("lastcheckpoint", ios_base::trunc);
	ofsCheck << bChkFirst << endl;
	if(iStop)
	    CkExit();
	else
	    mainChare.doSimulation();
	}
    }

// The restart callback needs to be an array entry, so we have a short
// entry that simply calls the main entry.

void
TreePiece::restart() 
{
    mainChare.restart();
    }

void
Main::initialForces()
{
  double startTime;
  double tolerance = 0.01;	// tolerance for domain decomposition

  /***** Initial sorting of particles and Domain Decomposition *****/
  ckerr << "Initial domain decomposition ...";
  startTime = CkWallTimer();
  sorter.startSorting(dataManagerID, numTreePieces, tolerance,
	 	      CkCallbackResumeThread(), true);
  ckerr << " took " << (CkWallTimer() - startTime) << " seconds."
        << endl;
  /******** Tree Build *******/
  ckerr << "Building trees ...";
  startTime = CkWallTimer();
  treeProxy.buildTree(bucketSize, CkCallbackResumeThread());
  ckerr << " took " << (CkWallTimer() - startTime) << " seconds."
        << endl;
  
  if(param.bDoGravity) {
      /******** Force Computation ********/
      ckerr << "Calculating gravity (theta = " << theta
	    << ") ...";
      startTime = CkWallTimer();
      treeProxy.startIteration(0, CkCallbackResumeThread());
      ckerr << " took " << (CkWallTimer() - startTime) << " seconds."
	    << endl;
      }
  
#if COSMO_STATS > 0
  /********* TreePiece Statistics ********/
  ckerr << "Total statistics initial iteration :" << endl;
  CkReductionMsg *tps;
  treeProxy.collectStatistics(CkCallbackResumeThread((void*&)tps));
  ((TreePieceStatistics*)tps->getData())->printTo(ckerr);
  
  /********* Cache Statistics ********/
  CkReductionMsg *cs;
  streamingCache.collectStatistics(CkCallbackResumeThread((void*&)cs));
  ((CacheStatistics*)cs->getData())->printTo(ckerr);
#endif
  
  /* 
   ** Dump Frame Initialization
   */
  // Currently for restarts, we have to set iStartStep.  Once we have more
  // complete restarts, this may be changed.
  if(DumpFrameInit(dTime, 0.0, param.iStartStep > 0)
     && df[0]->dDumpFrameStep > 0) {
      /* Bring frame count up to correct place for restart. */
      while(df[0]->dStep + df[0]->dDumpFrameStep < param.iStartStep) {
	  df[0]->dStep += df[0]->dDumpFrameStep;
	  df[0]->nFrame++;
	  }
      // initialize the rest of the dumpframes

      if (param.iDirector > 1) {
	  int j;
	  for(j=0; j < param.iDirector; j++) {
	      df[j]->dStep = df[0]->dStep;
	      df[j]->dDumpFrameStep = df[0]->dDumpFrameStep;
	      df[j]->nFrame = df[0]->nFrame;
	      } 
	  }
      }
     
  doSimulation();
}

// Principal method which does all the coordination of the simulation
// over timesteps.

void
Main::doSimulation()
{
  double startTime;
  char achLogFileName[MAXPATHLEN];
  sprintf(achLogFileName, "%s.log", param.achOutName);


  for(int iStep = param.iStartStep+1; iStep <= param.nSteps; iStep++){
    if (killAt > 0 && killAt == iStep) break;
    
    if (verbosity) ckerr << "Starting big step " << iStep << endl;
    startTime = CkWallTimer();
    advanceBigStep(iStep-1);
    double stepTime = CkWallTimer() - startTime;
    ckerr << "Big step " << iStep << " took " << stepTime << " seconds."
	  << endl;

    if(iStep%param.iLogInterval == 0) {
	calcEnergy(dTime, stepTime, achLogFileName);
    }
    iStop = CheckForStop();
    /*
     * Writing of intermediate outputs can be done here.
     */
    if((param.bBenchmark == 0)
       && (iStep == param.nSteps || iStop || iStep%param.iOutInterval == 0)) {
	writeOutput(iStep);
    }
	  
    if(!iStop && param.iWallRunTime > 0) {
	if (param.iWallRunTime*60. - (CkWallTimer()-dSimStartTime)
	    < stepTime*1.5) {
	    ckout << "RunTime limit exceeded.  Writing checkpoint and exiting.\n";
	    ckout << "    iWallRunTime(sec): " << param.iWallRunTime*60
		  << " Time running: " << CkWallTimer() - dSimStartTime
		  << " Last step: " << stepTime << endl;
	    iStop = 1;
	    }
	}
    
    if((param.bBenchmark == 0)
       && (iStop || iStep%param.iCheckInterval == 0)) {
	char achCheckFileName[MAXPATHLEN];
	if(bChkFirst) {
	    sprintf(achCheckFileName, "%s.chk0", param.achOutName);
	    bChkFirst = 0;
	    }
	else {
	    sprintf(achCheckFileName, "%s.chk1", param.achOutName);
	    bChkFirst = 1;
	    }
	param.iStartStep = iStep; // update so that restart continues on
	bIsRestarting = 0;
	CkCallback cb(CkIndex_TreePiece::restart(), treeProxy[0]);
	CkStartCheckpoint(achCheckFileName, cb);
	return;
    }
    if (iStop) break;

  } // End of main computation loop

  /******** Shutdown process ********/

  if(param.nSteps == 0) {
      if(param.bDoDensity) {
	  ckout << "Calculating densities ...";
	  startTime = CkWallTimer();
	  treeProxy.startIterationSmooth(0, CkCallbackResumeThread());
	  ckout << " took " << (CkWallTimer() - startTime) << " seconds."
		<< endl;
	  }
      
    ckout << "Outputting accelerations and densities ...";
    startTime = CkWallTimer();
    if(printBinaryAcc)
      treeProxy[0].outputAccelerations(OrientedBox<double>(),
				    "acc2", CkCallbackResumeThread());
    else {
	treeProxy.reOrder(CkCallbackResumeThread());
	treeProxy[0].outputAccASCII("acc2", CkCallbackResumeThread());
	}

    treeProxy[0].outputDensityASCII("den", CkCallbackResumeThread());
    treeProxy[0].outputIOrderASCII("iord", CkCallbackResumeThread());
    ckout << " took " << (CkWallTimer() - startTime) << " seconds."
	  << endl;
  }
	
#if COSMO_STATS > 0
  ckerr << "Outputting statistics ...";
  startTime = CkWallTimer();
  Interval<unsigned int> dummy;
	
  treeProxy[0].outputStatistics(dummy, dummy, dummy, dummy, 0, CkCallbackResumeThread());

  ckerr << " took " << (CkWallTimer() - startTime) << " seconds." << endl;
#endif

  ckerr << "Done." << endl;
	
#ifdef HPM_COUNTER
  streamingCache.stopHPM(CkCallbackResumeThread());
#endif
  ckerr << endl << "******************" << endl << endl; 
  CkExit();
}

/**
 * Calculate various energy and momentum quantities, and output them
 * to a log file.
 */

void Main::calcEnergy(double dTime, double wallTime, char *achLogFileName) {
    CkReductionMsg *msg;
    treeProxy.calcEnergy(CkCallbackResumeThread((void*&)msg));
    double *dEnergy = (double *) msg->getData();

    FILE *fpLog = fopen(achLogFileName, "a");
    
    static int first = 1;
    
    double a = csmTime2Exp(param.csm, dTime);
    
    if(first) {
	fprintf(fpLog, "# time redshift TotalEVir TotalE Kinetic Virial Potential TotalECosmo Lx Ly Lz Wallclock\n");
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
    
    fprintf(fpLog, "%g %g %g %g %g %g %g %g %g %g %g %g\n", dTime, z,
	    dEnergy[0] + dEnergy[1], dEnergy[0]+dEnergy[2], dEnergy[0],
	    dEnergy[1], dEnergy[2], dEcosmo, dEnergy[3], dEnergy[4],
	    dEnergy[5], wallTime);
    fclose(fpLog);
    
    delete msg;
}

void Main::writeOutput(int iStep) 
{
    char achFile[256];
    double dOutTime;
    double dvFac;
    double startTime;
    
    sprintf(achFile,"%s.%06i",param.achOutName,iStep);
    if(param.csm->bComove) {
	dOutTime = csmTime2Exp(param.csm, dTime);
	dvFac = 1.0/(dOutTime*dOutTime);
	}
    else {
	dOutTime = dTime;
	dvFac = 1.0;
	}
    
    if(verbosity) {
	ckerr << "ReOrder particles ...";
	startTime = CkWallTimer();
	}
    
    treeProxy.reOrder(CkCallbackResumeThread());
    if(verbosity)
	ckerr << " took " << (CkWallTimer() - startTime) << " seconds."
	      << endl;
    
    if(verbosity) {
	ckerr << "Writing output ...";
	startTime = CkWallTimer();
	}
    treeProxy.setupWrite(0, 0, achFile, dOutTime, dvFac,
		      CkCallbackResumeThread());
    if(verbosity)
	ckerr << " took " << (CkWallTimer() - startTime) << " seconds."
	      << endl;
    }

int Main::adjust(int iKickRung) 
{
    CkReductionMsg *msg;
    double a = csmTime2Exp(param.csm,dTime);
    
    treeProxy.adjust(iKickRung, param.bEpsAccStep, param.bGravStep, param.dEta,
		  param.dDelta, 1.0/(a*a*a), CkCallbackResumeThread((void*&)msg));

    int iCurrMaxRung = *(int *)msg->getData();
    delete msg;
    return iCurrMaxRung;
    }

void Main::rungStats() 
{
    CkReductionMsg *msg;
    treeProxy.rungStats(CkCallbackResumeThread((void*&)msg));
    int *nInRung = (int *)msg->getData();
    
    ckout << "Rung distribution: (";
    for(int iRung = 0; iRung <= MAXRUNG; iRung++)
	ckout << " " << nInRung[iRung] << ",";
    ckout << ")" << endl;
    
    delete msg;
    }

/*
 * Taken from PKDGRAV msrDumpFrameInit()
 */
int
Main::DumpFrameInit(double dTime, double dStep, int bRestart) {
	char achFile[256];
	
	if (param.dDumpFrameStep > 0 || param.dDumpFrameTime > 0) {
		bDumpFrame = 1;
		int i;

		df = (DumpFrameContext **) malloc(param.iDirector*sizeof(*df));
		for(i = 0; i < param.iDirector; i++) {
		    achFile[0] = '\0';
		    df[i] = NULL;
		    if(i == 0) {
			sprintf(achFile,"%s.director", param.achOutName);
			}
		    else {
			sprintf(achFile,"%s.director%d", param.achOutName, i+1);
			}
		    dfInitialize( &df[i], param.dSecUnit/SECONDSPERYEAR, 
				  dTime, param.dDumpFrameTime, dStep, 
				  param.dDumpFrameStep, param.dDelta, 
				  param.iMaxRung, verbosity, achFile );
		    }
		    
		/* Read in photogenic particle list */
		if (df[0]->bGetPhotogenic) {
		  achFile[0] = 0;
		  sprintf(achFile,"%s.photogenic", param.achOutName);

          CkReductionMsg *msg;
		  treeProxy.setTypeFromFile(TYPE_PHOTOGENIC, achFile, CkCallbackResumeThread((void*&)msg));
		  int *nSet = (int *)msg->getData();
		  if (verbosity)
		      ckout << nSet[0] << " iOrder numbers read. " << nSet[1]
			    <<" direct iOrder photogenic particles selected."
			    << endl;
		  delete msg;
		  }

		if(!bRestart)
		    DumpFrame(dTime, dStep );
                return 1;
		}
	else { return 0; }
    }

void Main::DumpFrame(double dTime, double dStep)
{
    double dExp;
    int i;

    for(i = 0; i < param.iDirector; i++) {
	if (df[i]->iDimension == DF_3D) {
#ifdef VOXEL
	    assert(0);		// Unimplemented
		/* 3D Voxel Projection */
		struct inDumpVoxel in;
		assert(0);

		dfSetupVoxel( msr->df, dTime, dStep, &in );

		pstDumpVoxel(msr->pst, &in, sizeof(struct inDumpVoxel), NULL, NULL );
		dsec1 = msrTime() - sec;
		
		dfFinishVoxel( msr->df, dTime, dStep, &in );
		
		dsec2 = msrTime() - sec;
		
		printf("DF Dumped Voxel %i at %g (Wallclock: Render %f tot %f secs).\n",
			   df->nFrame-1,dTime,dsec1,dsec2);
#endif
		}
	else {
		/* 2D Projection */
	    struct inDumpFrame in;
		double *com;
		CkReductionMsg *msgCOM;
		CkReductionMsg *msgCOMbyType;
		
		if (df[i]->bGetCentreOfMass) {
		    treeProxy.getCOM(CkCallbackResumeThread((void*&)msgCOM));
		    com = (double *)msgCOM->getData();
		    }

		if (df[i]->bGetPhotogenic) {
		    int iType = TYPE_PHOTOGENIC;
		    treeProxy.getCOMByType(iType, CkCallbackResumeThread((void*&)msgCOMbyType));
		    com = (double *)msgCOMbyType->getData();
		  }

#if (0)
		if (df[i]->bGetOldestStar) {
		  pstOldestStar(msr->pst, NULL, 0, &com[0], NULL);
		  }
#endif

		dExp = csmTime2Exp(param.csm,dTime);
		dfSetupFrame(df[i], dTime, dStep, dExp, com, &in );

        CkReductionMsg *msgDF;
		treeProxy.DumpFrame(in, CkCallbackResumeThread((void*&)msgDF));
		// N.B. Beginning of message contains the DumpFrame
		// parameters needed for proper merging.
		void *Image = ((char *)msgDF->getData())
		    + sizeof(struct inDumpFrame);

		dfFinishFrame(df[i], dTime, dStep, &in, Image );
		
		delete msgDF;

		if (df[i]->bGetCentreOfMass) {
		    delete msgCOM;
		    }
		if (df[i]->bGetPhotogenic) {
		    delete msgCOMbyType;
		    }
		}
	}
    }

void registerStatistics() {
#if COSMO_STATS > 0
  CacheStatistics::sum = CkReduction::addReducer(CacheStatistics::sumFn);
  TreePieceStatistics::sum = CkReduction::addReducer(TreePieceStatistics::sumFn);
#endif
#ifdef HPM_COUNTER
  hpmInit(1,"ChaNGa");
#endif
}

void Main::pup(PUP::er& p) 
{
    CBase_Main::pup(p);
    p | basefilename;
    p | theta;
    p | dTime;
    p | dEcosmo;
    p | dUOld;
    p | dTimeOld;
    p | printBinaryAcc;
    p | param;
    p | bDumpFrame;
    p | bChkFirst;
    }

#include "ParallelGravity.def.h"

/*
#define CK_TEMPLATES_ONLY
#include "CacheManager.def.h"
#undef CK_TEMPLATES_ONLY
#include "CacheManager.def.h"
*/
