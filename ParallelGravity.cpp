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
#include "smooth.h"
#include "Sph.h"

#ifdef CUDA
// for default per-list parameters
#include "cuda_typedef.h"
#endif

extern char *optarg;
extern int optind, opterr, optopt;

using namespace std;

CProxy_Main mainChare;
int verbosity;
int bVDetails;
CProxy_TreePiece treeProxy; // Proxy for the TreePiece chare array
CProxy_LvArray lvProxy;	    // Proxy for the liveViz array
CProxy_LvArray smoothProxy; // Proxy for smooth reductions
CProxy_CkCacheManager cacheGravPart;
CProxy_CkCacheManager cacheSmoothPart;
CProxy_CkCacheManager cacheNode;
CProxy_DataManager dMProxy;
bool _cache;
int _nocache;
int _cacheLineDepth;
unsigned int _yieldPeriod;
DomainsDec domainDecomposition;
int peanoKey;
GenericTrees useTree;
CProxy_TreePiece streamingProxy;
unsigned int numTreePieces;
unsigned int particlesPerChare;
int _prefetch;
int _numChunks;
int _randChunks;
unsigned int bucketSize;
int lbcomm_cutoff_msgs;

double dFracNoDomainDecomp; // Dummy for backward compatibility

//jetley
int localNodesPerReq;
int remoteNodesPerReq;
int remoteResumeNodesPerReq;
int localPartsPerReq;
int remotePartsPerReq;
int remoteResumePartsPerReq;

// multi-stepping particle transfer strategy 
// switch threshold
double largePhaseThreshold;

double theta;
double thetaMono;
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
    puts("USAGE: ChaNGa [SETTINGS | FLAGS] [PARAM_FILE]");
    puts("PARAM_FILE: Configuration file of a particular simulation, which");
    puts("          includes desired settings and relevant input and");
    puts("          output files. Settings specified in this file override");
    puts("          the default settings.");
    puts("SETTINGS");
    puts("or FLAGS: Command line settings or flags for a simulation which");
    puts("          will override any defaults and any settings or flags");
    puts("          specified in the PARAM_FILE.");
}


void _Trailer(void) {
	puts("(see the web page at\nhttp://librarian.phys.washington.edu/astro/index.php/Research:ChaNGa\nfor more information)");
}

int killAt;
int cacheSize;

Main::Main(CkArgMsg* m) {
	args = m;
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
#ifdef CUDA_TRACE 
        traceRegisterUserEvent("Tree Serialization", CUDA_SER_TREE);
        traceRegisterUserEvent("List Serialization", CUDA_SER_LIST);

        traceRegisterUserEvent("Local Node", CUDA_LOCAL_NODE_KERNEL);
        traceRegisterUserEvent("Remote Node", CUDA_REMOTE_NODE_KERNEL);
        traceRegisterUserEvent("Remote Resume Node", CUDA_REMOTE_RESUME_NODE_KERNEL);
        traceRegisterUserEvent("Local Particle", CUDA_LOCAL_PART_KERNEL);
        traceRegisterUserEvent("Remote Particle", CUDA_REMOTE_PART_KERNEL);
        traceRegisterUserEvent("Remote Resume Particle", CUDA_REMOTE_RESUME_PART_KERNEL);
#endif
	
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
		    sizeof(int),"killat", "Stop the simulation after this step");

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
	param.dSoftMax = 0.0;
	prmAddParam(prm,"dSoftMax",paramDouble, &param.dSoftMax,
		    sizeof(double),"eMax", "maximum comoving gravitational softening length (abs or multiplier)");
	param.bPhysicalSoft = 0;
	prmAddParam(prm,"bPhysicalSoft", paramBool, &param.bPhysicalSoft,
		    sizeof(int),"PhysSoft", "Physical gravitational softening length");
	param.bSoftMaxMul = 1;
	prmAddParam(prm,"bSoftMaxMul", paramBool, &param.bSoftMaxMul,
		    sizeof(int),"SMM", "Use maximum comoving gravitational softening length as a multiplier +SMM");
	param.dTheta = 0.7;
	prmAddParam(prm, "dTheta", paramDouble, &param.dTheta,
		    sizeof(double), "theta", "Opening angle");
	param.dTheta2 = 0.7;
	prmAddParam(prm, "dTheta2", paramDouble, &param.dTheta2,
		    sizeof(double),"theta2",
		    "Opening angle after switchTheta");
	param.daSwitchTheta = 1./3.;
	prmAddParam(prm,"daSwitchTheta",paramDouble,&param.daSwitchTheta,
		    sizeof(double),"aSwitchTheta",
		    "<a to switch theta at> = 1./3.");
#ifdef HEXADECAPOLE
	param.iOrder = 4;
#else
	param.iOrder = 2;
#endif
	prmAddParam(prm, "iOrder", paramInt, &param.iOrder,
		    sizeof(int), "or", "Multipole expansion order(IGNORED)");
	//
	// Cosmology parameters
	//
	param.bPeriodic = 0;
	prmAddParam(prm, "bPeriodic", paramBool, &param.bPeriodic,
		    sizeof(int),"per", "Periodic Boundaries");
	param.nReplicas = 0;
	prmAddParam(prm, "nReplicas", paramInt, &param.nReplicas,
		    sizeof(int),"nrep", "Number of periodic replicas");
	param.fPeriod = 1.0;
	prmAddParam(prm, "dPeriod", paramDouble, &param.fPeriod,
		    sizeof(double),"L", "Periodic size");
	param.vPeriod.x = 1.0;
	prmAddParam(prm,"dxPeriod",paramDouble,&param.vPeriod.x,
		    sizeof(double),"Lx",
		    "<periodic box length in x-dimension> = 1.0");
	param.vPeriod.y = 1.0;
	prmAddParam(prm, "dyPeriod",paramDouble,&param.vPeriod.y,
		    sizeof(double),"Ly",
		    "<periodic box length in y-dimension> = 1.0");
	param.vPeriod.z = 1.0;
	prmAddParam(prm,"dzPeriod",paramDouble,&param.vPeriod.z,
		    sizeof(double),"Lz",
		    "<periodic box length in z-dimension> = 1.0");
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
	
	//
	// Parameters for GrowMass: slowly growing mass of particles.
	//
	param.bDynGrowMass = 1;
	prmAddParam(prm,"bDynGrowMass",paramBool,&param.bDynGrowMass,
		    sizeof(int),"gmd","<dynamic growmass particles>");
	param.nGrowMass = 0;
	prmAddParam(prm,"nGrowMass",paramInt,&param.nGrowMass,sizeof(int),
		    "gmn","<number of particles to increase mass> = 0");
	param.dGrowDeltaM = 0.0;
	prmAddParam(prm,"dGrowDeltaM",paramDouble,&param.dGrowDeltaM,
		    sizeof(double),"gmdm",
		    "<Total growth in mass/particle> = 0.0");
	param.dGrowStartT = 0.0;
	prmAddParam(prm,"dGrowStartT",paramDouble,&param.dGrowStartT,
		    sizeof(double),"gmst",
		    "<Start time for growing mass> = 0.0");
	param.dGrowEndT = 1.0;
	prmAddParam(prm,"dGrowEndT",paramDouble,&param.dGrowEndT,
		    sizeof(double),"gmet","<End time for growing mass> = 1.0");

	//
	// Gas parameters
	//
	param.bDoGas = 0;
	prmAddParam(prm, "bDoGas", paramBool, &param.bDoGas,
		    sizeof(int),"gas", "Enable Gas Calculation");
	param.bFastGas = 0;
	prmAddParam(prm, "bFastGas", paramBool, &param.bFastGas,
		    sizeof(int),"Fgas", "Fast Gas Method");
	param.dFracFastGas = 0.0;
	prmAddParam(prm,"dFracFastGas",paramDouble,&param.dFracFastGas,
		    sizeof(double),"ffg",
		    "<Fraction of Active Particles for Fast Gas>");
	param.dhMinOverSoft = 0.0;
	prmAddParam(prm,"dhMinOverSoft",paramDouble,&param.dhMinOverSoft,
		    sizeof(double),"hmin",
		    "<Minimum h as a fraction of Softening> = 0.0");
	param.bViscosityLimiter = 1;
	prmAddParam(prm,"bViscosityLimiter",paramBool,&param.bViscosityLimiter,
		    sizeof(int), "vlim","<Viscosity Limiter> = 1");
	param.iViscosityLimiter = 1;
	prmAddParam(prm,"iViscosityLimiter",paramInt,&param.iViscosityLimiter,
		    sizeof(int), "ivlim","<Viscosity Limiter Type> = 1");
	param.bGeometric = 0;
	prmAddParam(prm,"bGeometric",paramBool,&param.bGeometric,sizeof(int),
		    "geo","geometric/arithmetic mean to calc Grad(P/rho) = +geo");
	param.bGasAdiabatic = 0;
	prmAddParam(prm,"bGasAdiabatic",paramBool,&param.bGasAdiabatic,
		    sizeof(int),"GasAdiabatic",
		    "<Gas is Adiabatic> = +GasAdiabatic");
	param.bGasIsothermal = 0;
	prmAddParam(prm,"bGasIsothermal",paramBool,&param.bGasIsothermal,
		    sizeof(int),"GasIsothermal",
		    "<Gas is Isothermal> = -GasIsothermal");
	param.bGasCooling = 0;
	prmAddParam(prm,"bGasCooling",paramBool,&param.bGasCooling,
		    sizeof(int),"GasCooling",
		    "<Gas is Cooling> = +GasCooling");
	CoolAddParams(&param.CoolParam, prm);
	param.dMsolUnit = 1.0;
	prmAddParam(prm,"dMsolUnit",paramDouble,&param.dMsolUnit,
		    sizeof(double),"msu", "<Solar mass/system mass unit>");
	param.dKpcUnit = 1000.0;
	prmAddParam(prm,"dKpcUnit",paramDouble,&param.dKpcUnit,
		    sizeof(double),"kpcu", "<Kiloparsec/system length unit>");
	param.ddHonHLimit = 0.1;
	prmAddParam(prm,"ddHonHLimit",paramDouble,&param.ddHonHLimit,
		    sizeof(double),"dhonh", "<|dH|/H Limiter> = 0.1");
	param.dConstAlpha = 1.0;
	prmAddParam(prm,"dConstAlpha",paramDouble,&param.dConstAlpha,
		    sizeof(double),"alpha",
		    "<Alpha constant in viscosity> = 1.0");
	param.dConstBeta = 2.0;
	prmAddParam(prm,"dConstBeta",paramDouble,&param.dConstBeta,
		    sizeof(double),"beta",
		    "<Beta constant in viscosity> = 2.0");
	param.dConstGamma = 5.0/3.0;
	prmAddParam(prm,"dConstGamma",paramDouble,&param.dConstGamma,
		    sizeof(double),"gamma", "<Ratio of specific heats> = 5/3");
	param.dMeanMolWeight = 1.0;
	prmAddParam(prm,"dMeanMolWeight",paramDouble,&param.dMeanMolWeight,
		    sizeof(double),"mmw",
		    "<Mean molecular weight in amu> = 1.0");
	param.dGasConst = 1.0;
	prmAddParam(prm,"dGasConst",paramDouble,&param.dGasConst,
		    sizeof(double),"gcnst", "<Gas Constant>");
	param.bBulkViscosity = 0;
	prmAddParam(prm,"bBulkViscosity",paramBool,&param.bBulkViscosity,
		    sizeof(int), "bulk","<Bulk Viscosity> = 0");
	// SPH timestepping
	param.bSphStep = 1;
	prmAddParam(prm,"bSphStep",paramBool,&param.bSphStep,sizeof(int),
		    "ss","<SPH timestepping>");
	param.bViscosityLimitdt = 0;
	prmAddParam(prm,"bViscosityLimitdt",paramBool,
		    &param.bViscosityLimitdt,sizeof(int),
		    "vlim","<Balsara Viscosity Limit dt> = 0");
	param.dEtaCourant = 0.4;
	prmAddParam(prm,"dEtaCourant",paramDouble,&param.dEtaCourant,
		    sizeof(double),"etaC", "<Courant criterion> = 0.4");
	param.dEtauDot = 0.25;
	prmAddParam(prm,"dEtauDot",paramDouble,&param.dEtauDot,
		    sizeof(double),"etau", "<uDot criterion> = 0.25");
	param.nTruncateRung = 0;
        prmAddParam(prm,"nTruncateRung",paramInt,&param.nTruncateRung,
		    sizeof(int),"nTR",
		    "<number of MaxRung particles to delete MaxRung> = 0");

	//
	// Output parameters
	//
	printBinaryAcc=0;
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
	param.bParaWrite = 1;
	prmAddParam(prm, "bParaWrite", paramBool, &param.bParaWrite,sizeof(int),
		    "paw", "enable/disable parallel writing of files");
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

	param.bLiveViz = 0;
	prmAddParam(prm, "bLiveViz", paramInt,&param.bLiveViz, sizeof(int),
		    "liveviz", "enable real-time simulation render support (disabled)");

	
	param.bStaticTest = 0;
	prmAddParam(prm, "bStaticTest", paramBool, &param.bStaticTest,
		    sizeof(int),"st", "Static test of performance");
	
	_randChunks = 1;
	prmAddParam(prm, "bRandChunks", paramBool, &_randChunks,
		    sizeof(int),"rand", "Randomize the order of remote chunk computation (default: ON)");
	
	_nocache = 0;
	// prmAddParam(prm, "bNoCache", paramBool, &_nocache,
	//	    sizeof(int),"nc", "Disable the CacheManager caching behaviour");
	
        verbosity = 0;
	prmAddParam(prm, "iVerbosity", paramInt, &verbosity,
		    sizeof(int),"v", "Verbosity");
	bVDetails = 0;
	prmAddParam(prm, "bVDetails", paramBool, &bVDetails,
		    sizeof(int),"vdetails", "Verbosity");
	numTreePieces = 8 * CkNumPes();
	prmAddParam(prm, "nTreePieces", paramInt, &numTreePieces,
		    sizeof(int),"p", "Number of TreePieces (default: 8*procs)");
	particlesPerChare = 0;
	prmAddParam(prm, "nPartPerChare", paramInt, &particlesPerChare,
            sizeof(int),"ppc", "Average number of particles per TreePiece");
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
	cacheSize = 100000000;
	prmAddParam(prm, "nCacheSize", paramInt, &cacheSize,
		    sizeof(int),"s", "Size of cache (IGNORED)");
	domainDecomposition=SFC_dec;
        peanoKey=0;
	prmAddParam(prm, "nDomainDecompose", paramInt, &domainDecomposition,
		    sizeof(int),"D", "Kind of domain decomposition of particles");
	dFracNoDomainDecomp = 0.0;
	prmAddParam(prm, "dFracNoDomainDecomp", paramDouble, &dFracNoDomainDecomp,
		    sizeof(double),"fndd", "(IGNORED)");
        lbcomm_cutoff_msgs = 1;
	prmAddParam(prm, "lbcommCutoffMsgs", paramInt, &lbcomm_cutoff_msgs,
		    sizeof(int),"lbcommcut", "Cutoff for communication recording (IGNORED)");
	param.bConcurrentSph = 0;
	prmAddParam(prm, "bConcurrentSph", paramBool, &param.bConcurrentSph,
		    sizeof(int),"consph", "Enable SPH running concurrently with Gravity");
    
          // jetley - cuda parameters
#ifdef CUDA

          localNodesPerReqDouble = NODE_INTERACTIONS_PER_REQUEST_L;
	  prmAddParam(prm, "localNodesPerReq", paramDouble, &localNodesPerReqDouble,
                sizeof(double),"localnodes", "Num. local node interactions allowed per CUDA request");

          remoteNodesPerReqDouble = NODE_INTERACTIONS_PER_REQUEST_RNR;
	  prmAddParam(prm, "remoteNodesPerReq", paramDouble, &remoteNodesPerReqDouble,
                sizeof(double),"remotenodes", "Num. remote node interactions allowed per CUDA request");

          remoteResumeNodesPerReqDouble = NODE_INTERACTIONS_PER_REQUEST_RR;
	  prmAddParam(prm, "remoteResumeNodesPerReq", paramDouble, &remoteResumeNodesPerReqDouble,
                sizeof(double),"remoteresumenodes", "Num. remote resume node interactions allowed per CUDA request");

          localPartsPerReqDouble = PART_INTERACTIONS_PER_REQUEST_L;
            prmAddParam(prm, "localPartsPerReq", paramDouble, &localPartsPerReqDouble,
                sizeof(double),"localparts", "Num. local particle interactions allowed per CUDA request");

          remotePartsPerReqDouble = PART_INTERACTIONS_PER_REQUEST_RNR;
            prmAddParam(prm, "remotePartsPerReq", paramDouble, &remotePartsPerReqDouble,
                sizeof(double),"remoteparts", "Num. remote particle interactions allowed per CUDA request");

          remoteResumePartsPerReqDouble = PART_INTERACTIONS_PER_REQUEST_RR;
          prmAddParam(prm, "remoteResumePartsPerReq", paramDouble, &remoteResumePartsPerReqDouble,
              sizeof(double),"remoteresumeparts", "Num. remote resume particle interactions allowed per CUDA request");

          largePhaseThreshold = TP_LARGE_PHASE_THRESHOLD_DEFAULT;
          prmAddParam(prm, "largePhaseThreshold", paramDouble, &largePhaseThreshold,
              sizeof(double),"largephasethresh", "Ratio of active to total particles at which all particles (not just active ones) are sent to gpu in the target buffer (No source particles are sent.)");

#endif

	if(!prmArgProc(prm,m->argc,m->argv)) {
	    CkExit();
        }
	
	if(bVDetails)
	    verbosity = 1;
	
	if(prmSpecified(prm, "iMaxRung")) {
	    ckerr << "WARNING: ";
	    ckerr << "iMaxRung parameter ignored. MaxRung is " << MAXRUNG
		  << endl;
	    }
	if(prmSpecified(prm, "nTruncateRung")) {
	    ckerr << "WARNING: ";
	    ckerr << "nTruncateRung parameter ignored." << endl;
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
	if(!prmSpecified(prm, "dTheta2")) {
	    param.dTheta2 = param.dTheta;
	    }
	theta = param.dTheta;
        thetaMono = theta*theta*theta*theta;
	if(prmSpecified(prm, "dExtraStore")) {
	    ckerr << "WARNING: ";
	    ckerr << "dExtraStore parameter ignored."
		  << endl;
	    }
	if(prmSpecified(prm, "bCannonical")) {
	    ckerr << "WARNING: ";
	    ckerr << "bCannonical parameter ignored; integration is always cannonical"
		  << endl;
	    }
	if(prmSpecified(prm, "bOverwrite")) {
	    ckerr << "WARNING: ";
	    ckerr << "bOverwrite parameter ignored."
		  << endl;
	    }
	if(prmSpecified(prm, "dFracNoDomainDecomp")) {
	    ckerr << "WARNING: ";
	    ckerr << "dFracNoDomainDecomp parameter ignored."
		  << endl;
	    }
	    
	if (param.bPhysicalSoft) {
	    if (!param.csm->bComove) {
		ckerr << "WARNING: bPhysicalSoft reset to 0 for non-comoving (bComove == 0)\n";
		param.bPhysicalSoft = 0;
		}
#ifndef CHANGESOFT
	    ckerr << "ERROR: You must compile with -DCHANGESOFT to use changing softening options\n";
	    CkAbort("CHANGESOFT needed");
#endif
	    }
	
	/*
	 ** Make sure that we have some setting for nReplicas if
	 ** bPeriodic is set.
	 */
	if (param.bPeriodic && !prmSpecified(prm,"nReplicas")) {
	    param.nReplicas = 1;
	    }
	/*
	 ** Warn that we have a setting for nReplicas if bPeriodic NOT set.
	 */
	if (!param.bPeriodic && param.nReplicas != 0) {
	    CkPrintf("WARNING: nReplicas set to non-zero value for non-periodic!\n");
	    }
	/*
	 ** Determine the period of the box that we are using.
	 ** Set the new d[xyz]Period parameters which are now used instead
	 ** of a single dPeriod, but we still want to have compatibility
	 ** with the old method of setting dPeriod.
	 */
	if (prmSpecified(prm,"dPeriod") && !prmSpecified(prm,"dxPeriod")) {
	    param.vPeriod.x = param.fPeriod;
	    }
	if (prmSpecified(prm,"dPeriod") && !prmSpecified(prm,"dyPeriod")) {
	    param.vPeriod.y = param.fPeriod;
	    }
	if (prmSpecified(prm,"dPeriod") && !prmSpecified(prm,"dzPeriod")) {
	    param.vPeriod.z = param.fPeriod;
	    }
	/*
	 ** Periodic boundary conditions can be disabled along any of the
	 ** x,y,z axes by specifying a period of zero for the given axis.
	 ** Internally, the period is set to infinity (Cf. pkdBucketWalk()
	 ** and pkdDrift(); also the INTERSECT() macro in smooth.h).
	 */
	if (param.fPeriod  == 0) param.fPeriod  = 1.0e38;
	if (param.vPeriod.x == 0) param.vPeriod.x = 1.0e38;
	if (param.vPeriod.y == 0) param.vPeriod.y = 1.0e38;
	if (param.vPeriod.z == 0) param.vPeriod.z = 1.0e38;
	if(!param.bPeriodic) {
	    param.nReplicas = 0;
	    param.fPeriod = 1.0e38;
	    param.vPeriod = Vector3D<double>(1.0e38);
	    param.bEwald = 0;
	    }

#ifdef CUDA
          double mil = 1e6;
          localNodesPerReq = (int) (localNodesPerReqDouble * mil);
          remoteNodesPerReq = (int) (remoteNodesPerReqDouble * mil);
          remoteResumeNodesPerReq = (int) (remoteResumeNodesPerReqDouble * mil);
          localPartsPerReq = (int) (localPartsPerReqDouble * mil);
          remotePartsPerReq = (int) (remotePartsPerReqDouble * mil);
          remoteResumePartsPerReq = (int) (remoteResumePartsPerReqDouble * mil);

          ckout << "INFO: ";
          ckout << "localNodes: " << localNodesPerReqDouble << endl;
          ckout << "INFO: ";
          ckout << "remoteNodes: " << remoteNodesPerReqDouble << endl;
          ckout << "INFO: ";
          ckout << "remoteResumeNodes: " << remoteResumeNodesPerReqDouble << endl;
          ckout << "INFO: ";
          ckout << "localParts: " << localPartsPerReqDouble << endl;
          ckout << "INFO: ";
          ckout << "remoteParts: " << remotePartsPerReqDouble << endl;
          ckout << "INFO: ";
          ckout << "remoteResumeParts: " << remoteResumePartsPerReqDouble << endl;

          ckout << "INFO: largePhaseThreshold: " << largePhaseThreshold << endl;

          if(numTreePieces > CkNumPes() || numTreePieces <= 0){
            numTreePieces = CkNumPes();
          }
#endif

	if(prmSpecified(prm, "lbcommCutoffMsgs")) {
	    ckerr << "WARNING: ";
	    ckerr << "lbcommcut parameter ignored." << endl;
	    }
	    
	if(prmSpecified(prm, "bGeometric")) {
	    ckerr << "WARNING: ";
	    ckerr << "bGeometric parameter ignored." << endl;
	    }
	if(prmSpecified(prm, "dhMinOverSoft")) {
	    ckerr << "WARNING: ";
	    ckerr << "dhMinOverSoft parameter ignored." << endl;
	    }
	if (prmSpecified(prm,"bViscosityLimiter")) {
	    if (!param.bViscosityLimiter) param.iViscosityLimiter=0;
	    }

	if(param.bGasCooling && (param.bGasIsothermal || param.bGasAdiabatic)) {
	    ckerr << "WARNING: ";
	    ckerr << "More than one gas model set.  ";
	    ckerr << "Defaulting to Cooling.";
	    ckerr << endl;
	    param.bGasCooling = 1;
	    param.bGasAdiabatic = 0;
	    param.bGasIsothermal = 0;
	    }
	if(param.bGasAdiabatic && param.bGasIsothermal) {
	    ckerr << "WARNING: ";
	    ckerr << "Both adiabatic and isothermal set.";
	    ckerr << "Defaulting to Adiabatic.";
	    ckerr << endl;
	    param.bGasAdiabatic = 1;
	    param.bGasIsothermal = 0;
	    }
	if(prmSpecified(prm, "bBulkViscosity")) {
	    ckerr << "WARNING: ";
	    ckerr << "bBulkViscosity parameter ignored." << endl;
	    }
#ifdef COOLING_NONE
        if(param.bGasCooling)
	    CkAbort("Gas cooling requested but not compiled in");
#endif
	if(param.bGasCooling) {
	    CkAssert(prmSpecified(prm, "dMsolUnit")
		     && prmSpecified(prm, "dKpcUnit"));
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
        if (domainDecomposition == SFC_peano_dec_2D) peanoKey = 2;
        if (domainDecomposition == SFC_peano_dec_3D) peanoKey = 3;

	// hardcoding some parameters, later may be full options
	if(domainDecomposition==ORB_dec){
	    useTree = Binary_ORB;
	    CkAbort("ORB decomposition known to be bad and not implemented");
	    }
	else { useTree = Binary_Oct; }

        ckerr << "Running on " << CkNumPes() << " processors/ "<< CkNumNodes() <<" nodes with " << numTreePieces << " TreePieces" << endl;

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
          case SFC_peano_dec_3D:
            ckerr << "Domain decomposition... new SFC Peano-Hilbert" << endl;
            break;
          case SFC_peano_dec_2D:
            ckerr << "Domain decomposition... 2D SFC Peano-Hilbert" << endl;
            break;
          default:
            CkAbort("None of the implemented decompositions specified");
          }
        }
	
	CkArrayOptions opts(numTreePieces); 
	if (domainDecomposition == Oct_dec) {
	  CProxy_RRMap myMap=CProxy_RRMap::ckNew(); 
	  opts.setMap(myMap);
	} else {
	  CProxy_BlockMap myMap=CProxy_BlockMap::ckNew(); 
	  opts.setMap(myMap);
	}

	CProxy_TreePiece pieces = CProxy_TreePiece::ckNew(opts);
	treeProxy = pieces;

	opts.bindTo(treeProxy);
	lvProxy = CProxy_LvArray::ckNew(opts);
	// Create an array for the smooth reductions
	smoothProxy = CProxy_LvArray::ckNew(opts);
	
	// create CacheManagers
	// Gravity particles
	cacheGravPart = CProxy_CkCacheManager::ckNew(cacheSize, pieces.ckLocMgr()->getGroupID());
	// Smooth particles
	cacheSmoothPart = CProxy_CkCacheManager::ckNew(cacheSize, pieces.ckLocMgr()->getGroupID());
	// Nodes: we need the right number of phases to keep the
	// nodes.
	nPhases = 0;
	if(param.bDoGravity) nPhases++;
	if(param.bDoGas) {
	    nPhases += 2;
	    if(param.bFastGas)
		nPhases += 2;
	    }
	if(nPhases == 0) {
	    ckerr << "Neither bDoGravity or bDoGas are set!" << endl;
	    CkAbort("Nothing to do!");
	    }
	CkGroupID *gids = new CkGroupID[nPhases];
	int i;
	for(i = 0; i < nPhases; i++)
	    gids[i] = pieces.ckLocMgr()->getGroupID();
	cacheNode = CProxy_CkCacheManager::ckNew(cacheSize, nPhases, gids);
	delete gids;

	//create the DataManager
	CProxy_DataManager dataManager = CProxy_DataManager::ckNew(pieces);
	dataManagerID = dataManager;
        dMProxy = dataManager;

	streamingProxy = pieces;

	//create the Sorter
	sorter = CProxy_Sorter::ckNew(0);

	if(verbosity)
	  ckerr << "Created " << numTreePieces << " pieces of tree" << endl;
	
	CProxy_Main(thishandle).setupICs();
}

Main::Main(CkMigrateMessage* m) {
    args = (CkArgMsg *)malloc(sizeof(*args));
    args->argv = CmiCopyArgs(((CkArgMsg *) m)->argv);
    mainChare = thishandle;
    bIsRestarting = 1;
    CkPrintf("Main(CkMigrateMessage) called\n");
    sorter = CProxy_Sorter::ckNew(0);
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
	    double dAStart = treeProxy[0].ckLocal()->dStartTime;
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
	    dTime = treeProxy[0].ckLocal()->dStartTime; 
	    if(verbosity > 0)
		ckout << "Input file, Time:" << dTime << endl;
	    double tTo = dTime + param.nSteps*param.dDelta;
	    if (verbosity > 0) {
		ckout << "Simulation to Time:" << tTo << endl;
		}
	    }
    }

// Return true if we need to write an output
int Main::bOutTime()
{	
    if (iOut < vdOutTime.size()) {
	if (dTime >= vdOutTime[iOut]) {
	    ++iOut;
	    return(1);
	    }
	else return(0);
	}
    else return(0);
    }

// Read in desired output times and reshifts from a file
//
void Main::getOutTimes()
{
    FILE *fp;
    int ret;
    double z,a,n,t;
    char achIn[80];
    char achFileName[MAXPATHLEN];
    sprintf(achFileName, "%s.red", param.achOutName);
	
    fp = fopen(achFileName,"r");
    if (!fp) {
	if (verbosity)
	    cerr << "WARNING: Could not open redshift input file: "
		 << achFileName << endl;
	return;
	}
    while (1) {
	if (!fgets(achIn,80,fp))
	    break;
	
	switch (achIn[0]) {
	case 'z':
	    ret = sscanf(&achIn[1],"%lf",&z);
	    if (ret != 1) break;
	    a = 1.0/(z+1.0);
	    vdOutTime.push_back(csmExp2Time(param.csm,a));
	    break;
	case 'a':
	    ret = sscanf(&achIn[1],"%lf",&a);
	    if (ret != 1) break;
	    vdOutTime.push_back(csmExp2Time(param.csm,a));
	    break;
	case 't':
	    ret = sscanf(&achIn[1],"%lf",&t);
	    if (ret != 1) break;
	    vdOutTime.push_back(csmExp2Time(param.csm,t));
	    break;
	case 'n':
	    ret = sscanf(&achIn[1],"%lf",&n);
	    if (ret != 1) break;
	    vdOutTime.push_back(csmExp2Time(param.csm,
					    dTime + (n-0.5)*param.dDelta));
	    break;
	default:
	    ret = sscanf(achIn,"%lf",&z);
	    if (ret != 1) break;
	    a = 1.0/(z+1.0);
	    vdOutTime.push_back(csmExp2Time(param.csm,a));
	    }
	if (ret != 1)
	    break;
	}
    /*
     ** Now sort the array of output times into ascending order.
     */
    vdOutTime.quickSort();
    fclose(fp);
    }

// determine if we need a smaller step for dumping frames
inline int Main::nextMaxRungIncDF(int nextMaxRung) 
{
    if (bDumpFrame && nextMaxRung < df[0]->iMaxRung)
	nextMaxRung = df[0]->iMaxRung;
    return nextMaxRung;
}

void Main::advanceBigStep(int iStep) {
  int currentStep = 0; // the current timestep within the big step
  int activeRung = 0;  // the minimum rung that is active
  int nextMaxRung = 0; // the rung that determines the smallest time for advancing

  while (currentStep < MAXSUBSTEPS) {

    if(!param.bStaticTest) {
      CkAssert(param.dDelta != 0.0);
      // Find new rung for active particles
      nextMaxRung = adjust(activeRung);
      if(currentStep == 0) rungStats();
      if(verbosity) ckout << "MaxRung: " << nextMaxRung << endl;

      CkAssert(nextMaxRung <= MAXRUNG); // timesteps WAY too short
 
      // Opening Kick
      if(verbosity)
	  ckout << "Kick Open:" << endl;
      double dKickFac[MAXRUNG+1];
      double duKick[MAXRUNG+1];
      for (int iRung=0; iRung<=MAXRUNG; iRung++) {
        duKick[iRung]=dKickFac[iRung]=0;
      }
      for(int iRung = activeRung; iRung <= nextMaxRung; iRung++) {
        double dTimeSub = RungToDt(param.dDelta, iRung);
	if(verbosity) {
	    ckout << " Rung " << iRung << ": " << 0.5*dTimeSub << endl;
	    }
        duKick[iRung] = 0.5*dTimeSub;
        dKickFac[iRung] = csmComoveKickFac(param.csm, dTime, 0.5*dTimeSub);
      }
      if(param.bDoGas) {
	  if(verbosity)
	      ckout << "uDot update:" << endl;
	  double z = 1.0/csmTime2Exp(param.csm,dTime) - 1.0;
	  if(param.bGasCooling)
	      dMProxy.CoolingSetTime(z, dTime, CkCallbackResumeThread());
	  treeProxy.updateuDot(activeRung, duKick, dTime, z, param.bGasCooling,
			       1, CkCallbackResumeThread());
	  }
      treeProxy.kick(activeRung, dKickFac, 0, param.bDoGas,
		     param.bGasIsothermal, duKick, CkCallbackResumeThread());

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
		  ckout << "Drift: Rung " << driftRung << " Delta "
			<< dTimeSub << endl;

	      // Only effective if growmass parameters have been set.
	      growMass(dTime, dTimeSub);
	      // Are the GrowMass particles locked in place?
	      int nGrowMassDrift = param.nGrowMass;
	      if(param.bDynGrowMass) nGrowMassDrift = 0;
	      
	      double dDriftFac = csmComoveDriftFac(param.csm, dTime, dTimeSub);
	      double dKickFac = csmComoveKickFac(param.csm, dTime, dTimeSub);
	      treeProxy.drift(dDriftFac, param.bDoGas, param.bGasIsothermal,
			      dKickFac, dTimeSub, nGrowMassDrift,
			      CkCallbackResumeThread());

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

    // int lastActiveRung = activeRung;

    // determine largest timestep that needs a kick
    activeRung = 0;
    int tmpRung = currentStep;
    while (tmpRung &= ~MAXSUBSTEPS) {
      activeRung++;
      tmpRung <<= 1;
    }

    if(param.bDoGas) {
	double duKick[MAXRUNG+1];
	if(verbosity)
	    ckout << "uDot update:" << endl;
	for(int iRung = activeRung; iRung <= nextMaxRung; iRung++) {
	    double dTimeSub = RungToDt(param.dDelta, iRung);
	    if(verbosity) {
		ckout << " Rung " << iRung << ": " << 0.5*dTimeSub << endl;
		}
	    duKick[iRung] = 0.5*dTimeSub;
	  }
	double z = 1.0/csmTime2Exp(param.csm,dTime) - 1.0;
	if(param.bGasCooling)
	    dMProxy.CoolingSetTime(z, dTime, CkCallbackResumeThread());
	treeProxy.updateuDot(activeRung, duKick, dTime, z, param.bGasCooling,
			     0, CkCallbackResumeThread());
	}

    ckout << "Step: " << (iStep + ((double) currentStep)/MAXSUBSTEPS)
          << " Time: " << dTime
          << " Rungs " << activeRung << " to "
          << nextMaxRung << ". ";

    countActive(activeRung);
    
    /***** Resorting of particles and Domain Decomposition *****/
    ckout << "Domain decomposition ...";
    double startTime = CkWallTimer();
    double tolerance = 0.01;	// tolerance for domain decomposition
    sorter.startSorting(dataManagerID, tolerance,
                        CkCallbackResumeThread(), true);
    ckout << " took " << (CkWallTimer() - startTime) << " seconds."
          << endl;

    /********* Load balancer ********/
    // jetley - commenting out lastActiveRung == 0 check, balance load even for fast rungs
    //if(lastActiveRung == 0) {
	ckout << "Load balancer ...";
	startTime = CkWallTimer();
	treeProxy.startlb(CkCallbackResumeThread(), activeRung);
	ckout << " took "<<(CkWallTimer() - startTime) << " seconds."
	     << endl;
    //	}

    /******** Tree Build *******/
    ckout << "Building trees ...";
    startTime = CkWallTimer();
    treeProxy.buildTree(bucketSize, CkCallbackResumeThread());
    iPhase = 0;
    ckout << " took " << (CkWallTimer() - startTime) << " seconds."
          << endl;

    CkCallback cbGravity(CkCallback::resumeThread);

    if(param.bDoGravity) {
	updateSoft();
	if(param.csm->bComove) {
	    double a = csmTime2Exp(param.csm,dTime);
	    if (a >= param.daSwitchTheta) theta = param.dTheta2; 
	    }
	/******** Force Computation ********/
	ckout << "Calculating gravity (tree bucket, theta = " << theta
	      << ") ...";
	startTime = CkWallTimer();
	if(param.bConcurrentSph) {
	    ckout << endl;
	    treeProxy.startIteration(activeRung, theta, cbGravity);
#ifdef CUDA_INSTRUMENT_WRS
            dMProxy.clearInstrument(cbGravity);
#endif
	    }
	else {
	    treeProxy.startIteration(activeRung, theta,
				     CkCallbackResumeThread());
#ifdef CUDA_INSTRUMENT_WRS
            dMProxy.clearInstrument(CkCallbackResumeThread());
#endif
	    ckout << " took " << (CkWallTimer() - startTime) << " seconds."
		  << endl;
	    }
	iPhase++;
    }
    else {
	treeProxy.initAccel(activeRung, CkCallbackResumeThread());
	}
    if(param.bDoGas) {
	doSph(activeRung);
	}
    
    if(param.bConcurrentSph && param.bDoGravity) {
	CkFreeMsg(cbGravity.thread_delay());
	ckout << "Calculating gravity and SPH took "
	      << (CkWallTimer() - startTime) << " seconds." << endl;
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
    cacheNode.collectStatistics(CkCallbackResumeThread((void*&)cs));
    ((CkCacheStatistics*)cs->getData())->printTo(ckerr);
    cacheGravPart.collectStatistics(CkCallbackResumeThread((void*&)cs));
    ((CkCacheStatistics*)cs->getData())->printTo(ckerr);
    cacheSmoothPart.collectStatistics(CkCallbackResumeThread((void*&)cs));
    ((CkCacheStatistics*)cs->getData())->printTo(ckerr);
#endif

    CkAssert(iPhase <= nPhases);
    
    if(iPhase < nPhases)
	treeProxy.finishNodeCache(nPhases-iPhase, CkCallbackResumeThread());

    if(!param.bStaticTest) {
      // Closing Kick
      double dKickFac[MAXRUNG+1];
      double duKick[MAXRUNG+1];
      for (int iRung=0; iRung<=MAXRUNG; iRung++) {
        duKick[iRung]=dKickFac[iRung]=0;
      }
      if(verbosity)
	  ckout << "Kick Close:" << endl;
      for(int iRung = activeRung; iRung <= nextMaxRung; iRung++) {
        double dTimeSub = RungToDt(param.dDelta, iRung);
	if(verbosity) {
	    ckout << " Rung " << iRung << ": " << 0.5*dTimeSub << endl;
	    }
        duKick[iRung] = 0.5*dTimeSub;
        dKickFac[iRung] = csmComoveKickFac(param.csm,
                                           dTime - 0.5*dTimeSub,
                                           0.5*dTimeSub);
      }
      if(param.bDoGas) {
	  double z = 1.0/csmTime2Exp(param.csm,dTime) - 1.0;
	  if(param.bGasCooling)
	      dMProxy.CoolingSetTime(z, dTime, CkCallbackResumeThread());
	  treeProxy.updateuDot(activeRung, duKick, dTime, z, param.bGasCooling,
			       1, CkCallbackResumeThread());
	  }
      treeProxy.kick(activeRung, dKickFac, 1, param.bDoGas,
		     param.bGasIsothermal, duKick, CkCallbackResumeThread());
    }
		
  }
}
    
// Load particles into pieces

void Main::setupICs() {
  double startTime;

  // DEBUGGING
  //CkStartQD(CkCallback(CkIndex_TreePiece::quiescence(),treeProxy));

  treeProxy.setPeriodic(param.nReplicas, param.vPeriod, param.bEwald,
			param.dEwCut, param.dEwhCut, param.bPeriodic);

  /******** Particles Loading ********/
  ckout << "Loading particles ...";
  startTime = CkWallTimer();
  treeProxy.load(basefilename, CkCallbackResumeThread());
  if(!(treeProxy[0].ckLocal()->bLoaded)) {
    // Try loading Tipsy format; first just grab the header.
    ckout << " trying Tipsy ...";
	    
    Tipsy::PartialTipsyFile ptf(basefilename, 0, 1);
    if(!ptf.loadedSuccessfully()) {
      ckerr << endl << "Couldn't load the tipsy file \""
            << basefilename.c_str()
            << "\". Maybe it's not a tipsy file?" << endl;
      CkExit();
      return;
    }
    double dTuFac = param.dGasConst/(param.dConstGamma-1)/param.dMeanMolWeight;
    treeProxy.loadTipsy(basefilename, dTuFac, CkCallbackResumeThread());
  }	
  ckout << " took " << (CkWallTimer() - startTime) << " seconds."
        << endl;
	    
  nTotalParticles = treeProxy[0].ckLocal()->nTotalParticles;
  nTotalSPH = treeProxy[0].ckLocal()->nTotalSPH;
  nActiveGrav = nTotalParticles;
  nActiveSPH = nTotalSPH;
  ckout << "N: " << nTotalParticles << endl;
  if(particlesPerChare == 0)
      particlesPerChare = nTotalParticles/numTreePieces;
  
  if(nTotalSPH == 0 && param.bDoGas) {
      ckerr << "WARNING: no SPH particles and bDoGas is set\n";
      param.bDoGas = 0;
      }
  getStartTime();
  if(param.nSteps > 0) getOutTimes();
  for(iOut = 0; iOut < vdOutTime.size(); iOut++) {
      if(dTime < vdOutTime[iOut]) break;
      }
	
  if(param.bGasCooling) 
      initCooling();
  
  char achLogFileName[MAXPATHLEN];
  sprintf(achLogFileName, "%s.log", param.achOutName);
  ofstream ofsLog(achLogFileName, ios_base::trunc);
  ofsLog << "# Starting ChaNGa version 2.00" << endl;
  ofsLog << "#";		// Output command line
  for (int i = 0; i < args->argc; i++)
      ofsLog << " " << args->argv[i];
  ofsLog << endl;
  ofsLog << "# Running on " << CkNumPes() << " processors" << endl;
#ifdef __DATE__
#ifdef __TIME__
  ofsLog <<"# Code compiled: " << __DATE__ << " " << __TIME__ << endl;
#endif
#endif
  ofsLog << "# Preprocessor macros:";
#ifdef CHANGESOFT
  ofsLog << " CHANGESOFT";
#endif
#ifdef COOLING_NONE
  ofsLog << " COOLING_NONE";
#endif
#ifdef COOLING_COSMO
  ofsLog << " COOLING_COSMO";
#endif
#ifdef COOLING_PLANET
  ofsLog << " COOLING_PLANET";
#endif
#ifdef HEXADECAPOLE
  ofsLog << " HEXADECAPOLE";
#endif
#ifdef INTERLIST_VER
  ofsLog << " INTERLIST_VER:" << INTERLIST_VER;
#endif
  ofsLog << endl;
  ofsLog.close();
	
  prmLogParam(prm, achLogFileName);
	
  ofsLog.open(achLogFileName, ios_base::app);
  if(param.csm->bComove) {
      ofsLog << "# RedOut:";
      if(vdOutTime.size() == 0) ofsLog << " none";
      for(int i = 0; i < vdOutTime.size(); i++) {
	  double z = 1.0/csmTime2Exp(param.csm, vdOutTime[i]) - 1.0;
	  ofsLog << " " << z;
	  }
    }
  else {
      ofsLog << "# TimeOut:";
      if(vdOutTime.size() == 0) ofsLog << " none";
      for(int i = 0; i < vdOutTime.size(); i++) {
	  ofsLog << " " << vdOutTime[i];
	  }
    }
  ofsLog << endl;
  ofsLog.close();

  if(prmSpecified(prm,"dSoft")) {
    ckout << "Set Softening...\n";
    treeProxy.setSoft(param.dSoft);
  }
	
  if(param.bPeriodic) {	// puts all particles within the boundary
      ckout << "drift particles to reset" << endl;
      treeProxy.drift(0.0, 0, 0, 0.0, 0.0, 0, CkCallbackResumeThread());
      ckout << "end drift particles to reset" << endl;

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
	ofsLog << "#";
	for (int i = 0; i < CmiGetArgc(args->argv); i++)
	    ofsLog << " " << args->argv[i];
	ofsLog << endl;
	ofsLog.close();
	/*
	 * Parse command line parameters
	 */
	prmInitialize(&prm,_Leader,_Trailer);
	/*
	 * Add a subset of the parameters.
	 */
	prmAddParam(prm, "dDelta", paramDouble, &param.dDelta,
		    sizeof(double),"dt", "Base Timestep for integration");
	prmAddParam(prm, "nSteps", paramInt, &param.nSteps,
		    sizeof(int),"n", "Number of Timesteps");
	prmAddParam(prm,"iWallRunTime",paramInt,&param.iWallRunTime,
		    sizeof(int),"wall",
		    "<Maximum Wallclock time (in minutes) to run> = 0 = infinite");
	prmAddParam(prm, "nBucket", paramInt, &bucketSize,
		    sizeof(int),"b", "Particles per Bucket (default: 12)");
	prmAddParam(prm, "iOutInterval", paramInt, &param.iOutInterval,
		    sizeof(int),"oi", "Output Interval");
	prmAddParam(prm, "iCheckInterval", paramInt, &param.iCheckInterval,
		    sizeof(int),"oc", "Checkpoint Interval");
	prmAddParam(prm, "iVerbosity", paramInt, &verbosity,
		    sizeof(int),"v", "Verbosity");
	
	if(!prmArgOnlyProc(prm,CmiGetArgc(args->argv),args->argv)) {
	    CkExit();
	}
	treeProxy.drift(0.0, 0, 0, 0.0, 0.0, 0, CkCallbackResumeThread());
	if(param.bGasCooling) 
	    initCooling();
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
  ckout << "Initial domain decomposition ...";
  startTime = CkWallTimer();
  sorter.startSorting(dataManagerID, tolerance,
	 	      CkCallbackResumeThread(), true);
  ckout << " took " << (CkWallTimer() - startTime) << " seconds."
        << endl;
  /******** Tree Build *******/
  ckout << "Building trees ...";
  startTime = CkWallTimer();
  treeProxy.buildTree(bucketSize, CkCallbackResumeThread());
  ckout << " took " << (CkWallTimer() - startTime) << " seconds."
        << endl;
  iPhase = 0;
  
#ifdef CUDA
  ckout << "Init. Accel. ...";
  treeProxy.initAccel(0, CkCallbackResumeThread());
  ckout << " took " << (CkWallTimer() - startTime) << " seconds."
        << endl;
#endif
      
  CkCallback cbGravity(CkCallback::resumeThread);  // needed below to wait for gravity

  if(param.bDoGravity) {
      updateSoft();
      if(param.csm->bComove) {
	  double a = csmTime2Exp(param.csm,dTime);
	  if (a >= param.daSwitchTheta) theta = param.dTheta2; 
	  }
      /******** Force Computation ********/
      ckout << "Calculating gravity (theta = " << theta
	    << ") ...";
      startTime = CkWallTimer();
      if(param.bConcurrentSph) {
	  treeProxy.startIteration(0, theta, cbGravity);
#ifdef CUDA_INSTRUMENT_WRS
            dMProxy.clearInstrument(cbGravity);
#endif
	  }
      else {
	  treeProxy.startIteration(0, theta, CkCallbackResumeThread());
#ifdef CUDA_INSTRUMENT_WRS
            dMProxy.clearInstrument(CkCallbackResumeThread());
#endif
	  ckout << " took " << (CkWallTimer() - startTime) << " seconds."
		<< endl;
	  }
      iPhase++;
      }
  else {
      treeProxy.initAccel(0, CkCallbackResumeThread());
      }
  if(param.bDoGas) {
      initSph();
      }
  
  if(param.bConcurrentSph && param.bDoGravity) {
      CkFreeMsg(cbGravity.thread_delay());
	ckout << "Calculating gravity and SPH took "
	      << (CkWallTimer() - startTime) << " seconds." << endl;
      }
  
  CkAssert(iPhase <= nPhases);
  if(iPhase < nPhases)
      treeProxy.finishNodeCache(nPhases-iPhase, CkCallbackResumeThread());

  // Initial Log entry
  char achLogFileName[MAXPATHLEN];
  sprintf(achLogFileName, "%s.log", param.achOutName);
  calcEnergy(dTime, CkWallTimer() - startTime, achLogFileName);

#if COSMO_STATS > 0
  /********* TreePiece Statistics ********/
  ckerr << "Total statistics initial iteration :" << endl;
  CkReductionMsg *tps;
  treeProxy.collectStatistics(CkCallbackResumeThread((void*&)tps));
  ((TreePieceStatistics*)tps->getData())->printTo(ckerr);
  
  /********* Cache Statistics ********/
  CkReductionMsg *cs;
  cacheNode.collectStatistics(CkCallbackResumeThread((void*&)cs));
  ((CkCacheStatistics*)cs->getData())->printTo(ckerr);
  cacheGravPart.collectStatistics(CkCallbackResumeThread((void*&)cs));
  ((CkCacheStatistics*)cs->getData())->printTo(ckerr);
  cacheSmoothPart.collectStatistics(CkCallbackResumeThread((void*&)cs));
  ((CkCacheStatistics*)cs->getData())->printTo(ckerr);
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
     

  if (param.bLiveViz > 0) {
    ckout << "Initializing liveViz module..." << endl;

  /* Initialize the liveViz module */

    liveVizConfig lvConfig(true,true);
    
    int funcIndex = CkIndex_TreePiece::liveVizDumpFrameInit(NULL);
    CkCallback* lvcb = new CkCallback(funcIndex, treeProxy);

    liveVizInit(lvConfig, treeProxy, *lvcb);
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
    if (killAt > 0 && killAt == iStep) {
      ckout << "KillAT: Stopping after " << (CkWallTimer()-dSimStartTime) << " seconds\n";
      break;
    }
    
    if (verbosity) ckout << "Starting big step " << iStep << endl;
    startTime = CkWallTimer();
    advanceBigStep(iStep-1);
    double stepTime = CkWallTimer() - startTime;
    ckout << "Big step " << iStep << " took " << stepTime << " seconds."
	  << endl;

    if(iStep%param.iLogInterval == 0) {
	calcEnergy(dTime, stepTime, achLogFileName);
    }
    iStop = CheckForStop();
    /*
     * Writing of intermediate outputs can be done here.
     */
    if((param.bBenchmark == 0)
       && (bOutTime() || iStep == param.nSteps || iStop
	   || iStep%param.iOutInterval == 0)) {
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
      char achFile[256];
    
      sprintf(achFile,"%s.%06i",param.achOutName,0);
      if((!param.bDoGas) && param.bDoDensity) {
	  iPhase = 0;
	  // If gas isn't being calculated, we can do the total
	  // densities before we start the output.
	  ckout << "Calculating total densities ...";
	  DensitySmoothParams pDen(TYPE_GAS|TYPE_DARK|TYPE_STAR, 0);
	  startTime = CkWallTimer();
	  treeProxy.startIterationSmooth(&pDen, CkCallbackResumeThread());
	  iPhase++;
	  ckout << " took " << (CkWallTimer() - startTime) << " seconds."
		<< endl;
#if 0
	  // For testing concurrent Sph, we don't want to do the
	  // resmooth.

          ckout << "Recalculating densities ...";
          startTime = CkWallTimer();
	  treeProxy.startIterationReSmooth(&pDen, CkCallbackResumeThread());
	  iPhase++;
          ckout << " took " << (CkWallTimer() - startTime) << " seconds." << endl;
#endif
	  CkAssert(iPhase <= nPhases);
	  if(iPhase < nPhases)
	      treeProxy.finishNodeCache(nPhases-iPhase, CkCallbackResumeThread());
          ckout << "Reodering ...";
          startTime = CkWallTimer();
	  treeProxy.reOrder(CkCallbackResumeThread());
          ckout << " took " << (CkWallTimer() - startTime) << " seconds." << endl;
	  ckout << "Outputting densities ...";
	  startTime = CkWallTimer();
	  DenOutputParams pDenOut(string(achFile) + ".den");
	  treeProxy[0].outputASCII(pDenOut, param.bParaWrite, CkCallbackResumeThread());
	  ckout << " took " << (CkWallTimer() - startTime) << " seconds."
		<< endl;
	  ckout << "Outputting hsmooth ...";
	  HsmOutputParams pHsmOut(string(achFile) + ".hsmall");
	  treeProxy[0].outputASCII(pHsmOut, param.bParaWrite, CkCallbackResumeThread());
	  }
      else {
	  treeProxy.reOrder(CkCallbackResumeThread());
	  }

      writeOutput(0);
      if(param.bDoGas) {
	  ckout << "Outputting gas properties ...";
	  if(printBinaryAcc)
	      CkAssert(0);
	  else {
	      DenOutputParams pDenOut(string(achFile) + ".gasden");
	      treeProxy[0].outputASCII(pDenOut, param.bParaWrite, CkCallbackResumeThread());
	      PresOutputParams pPresOut(string(achFile) + ".pres");
	      treeProxy[0].outputASCII(pPresOut, param.bParaWrite, CkCallbackResumeThread());
	      HsmOutputParams pSphHOut(string(achFile) + ".SphH");
	      treeProxy[0].outputASCII(pSphHOut, param.bParaWrite, CkCallbackResumeThread());
	      DivVOutputParams pDivVOut(string(achFile) + ".divv");
	      treeProxy[0].outputASCII(pDivVOut, param.bParaWrite, CkCallbackResumeThread());
	      PDVOutputParams pPDVOut(string(achFile) + ".PdV");
	      treeProxy[0].outputASCII(pPDVOut, param.bParaWrite, CkCallbackResumeThread());
	      MuMaxOutputParams pMuMaxOut(string(achFile) + ".mumax");
	      treeProxy[0].outputASCII(pMuMaxOut, param.bParaWrite, CkCallbackResumeThread());
	      BSwOutputParams pBSwOut(string(achFile) + ".BSw");
	      treeProxy[0].outputASCII(pBSwOut, param.bParaWrite, CkCallbackResumeThread());
	      CsOutputParams pCsOut(string(achFile) + ".c");
	      treeProxy[0].outputASCII(pCsOut, param.bParaWrite, CkCallbackResumeThread());
#ifndef COOLING_NONE
	      EDotOutputParams pEDotOut(string(achFile) + ".eDot");
	      treeProxy[0].outputASCII(pEDotOut, param.bParaWrite,
				       CkCallbackResumeThread());
#endif
	      }
	  }
      ckout << "Outputting accelerations  ...";
      if(printBinaryAcc)
	  treeProxy[0].outputAccelerations(OrientedBox<double>(),
					   "acc2", CkCallbackResumeThread());
      else {
	  AccOutputParams pAcc(string(achFile) + ".acc2");
	  treeProxy[0].outputASCII(pAcc, param.bParaWrite, CkCallbackResumeThread());
	  }
#ifdef NEED_DT
      ckout << "Outputting dt ...";
      adjust(0);
      DtOutputParams pDt(string(achFile) + ".dt");
      treeProxy[0].outputASCII(pDt, param.bParaWrite, CkCallbackResumeThread());
#endif
      RungOutputParams pRung(string(achFile) + ".rung");
      treeProxy[0].outputASCII(pRung, param.bParaWrite, CkCallbackResumeThread());
      if(param.bDoGas && param.bDoDensity) {
	  double tolerance = 0.01;	// tolerance for domain decomposition
	  // The following call is to get the particles in key order
	  // before the sort.
	  treeProxy.drift(0.0, 0, 0, 0.0, 0.0, 0, CkCallbackResumeThread());
	  sorter.startSorting(dataManagerID, tolerance,
			      CkCallbackResumeThread(), true);
	  treeProxy.buildTree(bucketSize, CkCallbackResumeThread());
	  iPhase = 0;
	  
	  ckout << "Calculating total densities ...";
	  DensitySmoothParams pDen(TYPE_GAS|TYPE_DARK|TYPE_STAR, 0);
	  startTime = CkWallTimer();
	  treeProxy.startIterationSmooth(&pDen, CkCallbackResumeThread());
	  iPhase++;
	  ckout << " took " << (CkWallTimer() - startTime) << " seconds."
		<< endl;
          ckout << "Recalculating densities ...";
          startTime = CkWallTimer();
	  treeProxy.startIterationReSmooth(&pDen, CkCallbackResumeThread());
	  iPhase++;
          ckout << " took " << (CkWallTimer() - startTime) << " seconds." << endl;
	  CkAssert(iPhase <= nPhases);
	  if(iPhase < nPhases)
	      treeProxy.finishNodeCache(nPhases-iPhase, CkCallbackResumeThread());
          ckout << "Reodering ...";
          startTime = CkWallTimer();
	  treeProxy.reOrder(CkCallbackResumeThread());
          ckout << " took " << (CkWallTimer() - startTime) << " seconds." << endl;
	  ckout << "Outputting densities ...";
	  startTime = CkWallTimer();
	  DenOutputParams pDenOut(string(achFile) + ".den");
	  treeProxy[0].outputASCII(pDenOut, param.bParaWrite, CkCallbackResumeThread());
	  ckout << " took " << (CkWallTimer() - startTime) << " seconds."
		<< endl;
	  ckout << "Outputting hsmooth ...";
	  HsmOutputParams pHsmOut(string(achFile) + ".hsmall");
	  treeProxy[0].outputASCII(pHsmOut, param.bParaWrite, CkCallbackResumeThread());
	  }
      treeProxy[0].outputIOrderASCII(string(achFile) + ".iord",
				     CkCallbackResumeThread());
  }
	
#if COSMO_STATS > 0
  ckerr << "Outputting statistics ...";
  startTime = CkWallTimer();
  //Interval<unsigned int> dummy;
	
  treeProxy[0].outputStatistics(CkCallbackResumeThread());
  //treeProxy[0].outputStatistics(dummy, dummy, dummy, dummy, 0, CkCallbackResumeThread());

  ckerr << " took " << (CkWallTimer() - startTime) << " seconds." << endl;
#endif

  ckout << "Done." << endl;
	
#ifdef HPM_COUNTER
  cacheNode.stopHPM(CkCallbackResumeThread());
  cacheGravPart.stopHPM(CkCallbackResumeThread());
  cacheSmoothPart.stopHPM(CkCallbackResumeThread());
#endif
  ckout << endl << "******************" << endl << endl; 
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
    static int first = 1;

    FILE *fpLog = fopen(achLogFileName, "a");
    
    double a = csmTime2Exp(param.csm, dTime);
    
    if(first && !bIsRestarting) {
	fprintf(fpLog, "# time redshift TotalEVir TotalE Kinetic Virial Potential TotalECosmo Ethermal Lx Ly Lz Wallclock\n");
	dEcosmo = 0.0;
	first = 0;
	}
    else {
	/*
	 * Estimate integral (\dot a*U*dt) over the interval.
	 * Note that this is equal to integral (W*da) and the latter
	 * is more accurate when a is changing rapidly.
	 */
	if(param.csm->bComove) {
	    dEcosmo += 0.5*(a - csmTime2Exp(param.csm, dTimeOld))
		*(dEnergy[2] + dUOld);
	    }
	first = 0;
	}

    dUOld = dEnergy[2];
    dEnergy[2] *= a;
    dEnergy[1] *= a;
    dTimeOld = dTime;
    double z = 1.0/a - 1.0;
    double dEtotal = dEnergy[0] + dEnergy[2] - dEcosmo + a*a*dEnergy[6];
    // Energy using the virial of Clausius
    double dEtotalVir = dEnergy[0] + dEnergy[1] - dEcosmo + a*a*dEnergy[6];
    
    fprintf(fpLog, "%g %g %g %g %g %g %g %g %g %g %g %g %g\n", dTime, z,
	    dEtotalVir, dEtotal, dEnergy[0], dEnergy[1], dEnergy[2], dEcosmo,
	    dEnergy[6], dEnergy[3], dEnergy[4],
	    dEnergy[5], wallTime);
    fclose(fpLog);
    
    delete msg;
}

void Main::writeOutput(int iStep) 
{
    char achFile[256];
    double dOutTime;
    double dvFac;
    double startTime = 0.0;
    
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
	ckout << "ReOrder particles ...";
	startTime = CkWallTimer();
	}
    
    treeProxy.reOrder(CkCallbackResumeThread());
    if(verbosity)
	ckout << " took " << (CkWallTimer() - startTime) << " seconds."
	      << endl;
    
    if(verbosity) {
	ckout << "Writing output ...";
	startTime = CkWallTimer();
	}
    double duTFac = (param.dConstGamma-1)*param.dMeanMolWeight/param.dGasConst;
    if(param.bParaWrite)
    	treeProxy.setupWrite(0, 0, achFile, dOutTime, dvFac, duTFac,
			     param.bGasCooling, CkCallbackResumeThread());
    else
	treeProxy[0].serialWrite(0, achFile, dOutTime, dvFac, duTFac,
				 param.bGasCooling, CkCallbackResumeThread());
    if(verbosity)
	ckout << " took " << (CkWallTimer() - startTime) << " seconds."
	      << endl;
#ifndef COOLING_NONE
      Cool0OutputParams pCool0Out(string(achFile) + "." + COOL_ARRAY0_EXT);
      treeProxy[0].outputASCII(pCool0Out, param.bParaWrite,
			       CkCallbackResumeThread());
      Cool1OutputParams pCool1Out(string(achFile) + "." + COOL_ARRAY1_EXT);
      treeProxy[0].outputASCII(pCool1Out, param.bParaWrite,
			       CkCallbackResumeThread());
      Cool2OutputParams pCool2Out(string(achFile) + "." + COOL_ARRAY2_EXT);
      treeProxy[0].outputASCII(pCool2Out, param.bParaWrite,
			       CkCallbackResumeThread());
#endif
    if(param.nSteps != 0 && param.bDoDensity) {
	double tolerance = 0.01;	// tolerance for domain decomposition
	// The following call is to get the particles in key order
	// before the sort.
	treeProxy.drift(0.0, 0, 0, 0.0, 0.0, 0, CkCallbackResumeThread());
	sorter.startSorting(dataManagerID, tolerance,
			    CkCallbackResumeThread(), true);
	treeProxy.buildTree(bucketSize, CkCallbackResumeThread());

	iPhase = 0;
	ckout << "Calculating total densities ...";
	DensitySmoothParams pDen(TYPE_GAS|TYPE_DARK|TYPE_STAR, 0);
	startTime = CkWallTimer();
	treeProxy.startIterationSmooth(&pDen, CkCallbackResumeThread());
	iPhase++;
	if(iPhase < nPhases)
	    treeProxy.finishNodeCache(nPhases-iPhase, CkCallbackResumeThread());
	ckout << " took " << (CkWallTimer() - startTime) << " seconds."
	      << endl;
	ckout << "Reodering ...";
	startTime = CkWallTimer();
	treeProxy.reOrder(CkCallbackResumeThread());
	ckout << " took " << (CkWallTimer() - startTime) << " seconds."
	      << endl;
	ckout << "Outputting densities ...";
	startTime = CkWallTimer();
	DenOutputParams pDenOut(string(achFile) + ".den");
	treeProxy[0].outputASCII(pDenOut, param.bParaWrite, CkCallbackResumeThread());
	ckout << " took " << (CkWallTimer() - startTime) << " seconds."
	      << endl;
	}
    }

int Main::adjust(int iKickRung) 
{
    CkReductionMsg *msg;
    double a = csmTime2Exp(param.csm,dTime);
    
    treeProxy.adjust(iKickRung, param.bEpsAccStep, param.bGravStep,
		     param.bSphStep, param.bViscosityLimitdt,
		     param.dEta, param.dEtaCourant,
		     param.dEtauDot, param.dDelta, 1.0/(a*a*a), a,
		     CkCallbackResumeThread((void*&)msg));

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

void Main::countActive(int activeRung) 
{
    CkReductionMsg *msg;
    treeProxy.countActive(activeRung, CkCallbackResumeThread((void*&)msg));
    int *nActive = (int *)msg->getData();
    
    nActiveGrav = nActive[0];
    nActiveSPH = nActive[1];
    
    ckout << "Gravity Active: " << nActive[0]
	  << ", Gas Active: " << nActive[1] << endl ;
    
    delete msg;
    }

void Main::updateSoft()
{
#ifdef CHANGESOFT
    if (param.bPhysicalSoft) {
	 double dFac = 1./csmTime2Exp(param.csm,dTime);

	 if (param.bSoftMaxMul && dFac > param.dSoftMax) dFac = param.dSoftMax;

	 treeProxy.physicalSoft(param.dSoftMax, dFac, param.bSoftMaxMul,
				CkCallbackResumeThread());
	}
#endif
    }

//
// Slowly increase mass of a subset of particles.
//
void Main::growMass(double dTime, double dDelta)
{
    if (param.nGrowMass > 0 && dTime > param.dGrowStartT
	&& dTime <= param.dGrowEndT) {
	double dDeltaM = param.dGrowDeltaM*dDelta
			    /(param.dGrowEndT - param.dGrowStartT);
	treeProxy.growMass(param.nGrowMass, dDeltaM, CkCallbackResumeThread());
	}
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
		double *com = NULL;
		CkReductionMsg *msgCOM;
		CkReductionMsg *msgCOMbyType;
		
		if (df[i]->bGetCentreOfMass) {
		    treeProxy.getCOM(CkCallbackResumeThread((void*&)msgCOM), 0);
		    com = (double *)msgCOM->getData();
		    }

		if (df[i]->bGetPhotogenic) {
		    int iType = TYPE_PHOTOGENIC;
		    treeProxy.getCOMByType(iType, CkCallbackResumeThread((void*&)msgCOMbyType), 0);
		    com = (double *)msgCOMbyType->getData();
		  }

#if (0)
		if (df[i]->bGetOldestStar) {
		  pstOldestStar(msr->pst, NULL, 0, &com[0], NULL);
		  }
#endif

		dExp = csmTime2Exp(param.csm,dTime);
		dfSetupFrame(df[i], dTime, dStep, dExp, com, &in, 0, 0 );

        CkReductionMsg *msgDF;
	treeProxy.DumpFrame(in, CkCallbackResumeThread((void*&)msgDF), false);
		// N.B. Beginning of message contains the DumpFrame
		// parameters needed for proper merging.
		void *Image = ((char *)msgDF->getData())
		    + sizeof(struct inDumpFrame);

		unsigned char *gray;

		dfFinishFrame(df[i], dTime, dStep, &in, Image, false, &gray);
		
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
  CkCacheStatistics::sum = CkReduction::addReducer(CkCacheStatistics::sumFn);
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
    p | nTotalParticles;
    p | nTotalSPH;
    p | theta;
    p | dTime;
    p | dEcosmo;
    p | dUOld;
    p | dTimeOld;
    p | printBinaryAcc;
    p | param;
    p | vdOutTime;
    p | iOut;
    p | bDumpFrame;
    p | bChkFirst;
    p | nPhases;
    p | sorter;
    }


void Main::liveVizImagePrep(liveVizRequestMsg *msg) 
{
  double dExp;

  struct inDumpFrame in;
  double *com = NULL;
  CkReductionMsg *msgCOM;
  CkReductionMsg *msgCOMbyType;
  struct DumpFrameContext df_temp;

  // only do the first Director 

  df_temp = *df[0];  

  if (df_temp.bGetCentreOfMass) {
      treeProxy.getCOM(CkCallbackResumeThread((void*&)msgCOM), 1);
      com = (double *)msgCOM->getData();
      }

  if (df_temp.bGetPhotogenic) {
      int iType = TYPE_PHOTOGENIC;
      treeProxy.getCOMByType(iType,
			     CkCallbackResumeThread((void*&)msgCOMbyType), 1);
      com = (double *)msgCOMbyType->getData();
      }

#if (0)
    if (df_temp.bGetOldestStar) {
 	 pstOldestStar(msr->pst, NULL, 0, &com[0], NULL);
    }
#endif

  dExp = csmTime2Exp(param.csm,dTime);
  dfSetupFrame(&df_temp, dTime, 0.0, dExp, com, &in, 
		     msg->req.wid, msg->req.ht );
  CkCallback cb(CkCallback::ignore);
  treeProxy.DumpFrame(in, cb, true);

  if (df_temp.bGetCentreOfMass)
      delete msgCOM;
  if (df_temp.bGetPhotogenic)
      delete msgCOMbyType;
  
}

#include "ParallelGravity.def.h"

/*
#define CK_TEMPLATES_ONLY
#include "CacheManager.def.h"
#undef CK_TEMPLATES_ONLY
#include "CacheManager.def.h"
*/
