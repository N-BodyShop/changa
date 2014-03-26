/**
 * @mainpage ChaNGa: Charm++ N-body Gravity 
 * 
 * For internal documentation of the operation of this code please see the
 * following classes:
 * 
 * Main: Overall control flow of the program.
 * 
 * TreePiece: Structure that holds the particle data and tree
 * structure in each object.
 *
 * Sorter: Organize domain decomposition.
 *
 * DataManager: Organize data across processors.
 *
 * Tree::GenericTreeNode: Interface for the tree data structure.
 *
 * GravityParticle: Per particle data structure.
 */
 
#include <stdint.h>
#include <iostream>

#include <unistd.h>
#include <sys/param.h>
#include <sys/stat.h>

// Debug floating point problems
// #include <fenv.h>

#include "BaseLB.h"

#include "Sorter.h"
#include "ParallelGravity.h"
#include "DataManager.h"
#include "TipsyFile.h"
#include "param.h"
#include "smooth.h"
#include "Sph.h"
#include "starform.h"
#include "feedback.h"

#include "PETreeMerger.h"

#ifdef CUDA
// for default per-list parameters
#include "cuda_typedef.h"
#endif

extern char *optarg;
extern int optind, opterr, optopt;
extern const char * const Cha_CommitID;

using namespace std;

CProxy_Main mainChare;
int verbosity;
int bVDetails;
CProxy_TreePiece treeProxy; // Proxy for the TreePiece chare array
#ifdef REDUCTION_HELPER
CProxy_ReductionHelper reductionHelperProxy;
#endif
CProxy_LvArray lvProxy;	    // Proxy for the liveViz array
CProxy_LvArray smoothProxy; // Proxy for smooth reductions
CProxy_LvArray gravityProxy; // Proxy for gravity reductions
CProxy_CkCacheManager<KeyType> cacheGravPart;
CProxy_CkCacheManager<KeyType> cacheSmoothPart;
CProxy_CkCacheManager<KeyType> cacheNode;
CProxy_DataManager dMProxy;

CProxy_DumpFrameData dfDataProxy;
CProxy_PETreeMerger peTreeMergerProxy;




bool _cache;
int _nocache;
int _cacheLineDepth;
unsigned int _yieldPeriod;
DomainsDec domainDecomposition;
double dExtraStore;		// fraction of extra particle storage
double dMaxBalance;		// Max piece imbalance for load balancing
int iGasModel; 			// For backward compatibility
int peanoKey;
GenericTrees useTree;
CProxy_TreePiece streamingProxy;
unsigned int numTreePieces;
unsigned int particlesPerChare;
int nIOProcessor;		// Number of pieces to be doing I/O at once
int _prefetch;
int _numChunks;
int _randChunks;
unsigned int bucketSize;
int lbcomm_cutoff_msgs;

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

int boundaryEvaluationUE;
int weightBalanceUE;
int networkProgressUE;
int nodeForceUE;
int partForceUE;

int tbRecursiveUE;
int tbFlushRequestsUE;
int prefetchDoneUE;

CkGroupID dataManagerID;
CkArrayID treePieceID;

#ifdef PUSH_GRAVITY
#include "ckmulticast.h"
CkGroupID ckMulticastGrpId;
#endif

CProxy_ProjectionsControl prjgrp;

/* The following is for backward compatibility and deprecated */
#define GASMODEL_UNSET -1
enum GasModel {
	GASMODEL_ADIABATIC, 
	GASMODEL_ISOTHERMAL, 
	GASMODEL_COOLING, 
	GASMODEL_GLASS
	}; 

bool doDumpLB;
int lbDumpIteration;
bool doSimulateLB;

// Number of bins to use for the first iteration
// of every Oct decomposition step
int numInitDecompBins;

// Specifies the number of sub-bins a bin is split into
//  for Oct decomposition
int octRefineLevel;

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

///
/// @brief Main routine to start simulation.
///
/// This routine parses the command line and file specified parameters,
/// and allocates the charm structures.  The charm "readonly" variables
/// are writable in this method, and are broadcast globally once this
/// method exits.  Note that the method finishes with an asynchronous
/// call to setupICs().
///
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

        tbRecursiveUE = traceRegisterUserEvent("TreeBuild::buildOctTree::recursive");
        tbFlushRequestsUE = traceRegisterUserEvent("TreeBuild::buildOctTree::flushRequests");

        prefetchDoneUE = traceRegisterUserEvent("PrefetchDone");
	
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
	//
	// Output flags
	//
	param.iOutInterval = 10;
	prmAddParam(prm, "iOutInterval", paramInt, &param.iOutInterval,
		    sizeof(int),"oi", "Output Interval");
	param.iLogInterval = 1;
	prmAddParam(prm, "iLogInterval", paramInt, &param.iLogInterval,
		    sizeof(int),"ol", "Log Interval");
	param.iCheckInterval = 10;
	prmAddParam(prm, "iCheckInterval", paramInt, &param.iCheckInterval,
		    sizeof(int),"oc", "Checkpoint Interval");
	param.iBinaryOut = 0;
	prmAddParam(prm, "iBinaryOutput", paramInt, &param.iBinaryOut,
		    sizeof(int), "binout",
		    "<array outputs 0 ascii, 1 float, 2 double, 3 FLOAT(internal)> = 0");
	param.bDoDensity = 1;
	prmAddParam(prm, "bDoDensity", paramBool, &param.bDoDensity,
		    sizeof(int),"den", "Enable Density outputs");
	param.bDoIOrderOutput = 0;
	prmAddParam(prm,"bDoIOrderOutput",paramBool,&param.bDoIOrderOutput,
		    sizeof(int), "iordout","enable/disable iOrder outputs = -iordout");
	param.bDoSoftOutput = 0;
	prmAddParam(prm,"bDoSoftOutput",paramBool,&param.bDoSoftOutput,
		    sizeof(int),
		    "softout","enable/disable soft outputs = -softout");
	param.bDohOutput = 0;
	prmAddParam(prm,"bDohOutput",paramBool,&param.bDohOutput,sizeof(int),
		    "hout","enable/disable h outputs = -hout");
	param.bDoCSound = 0;
	prmAddParam(prm,"bDoCSound",paramBool,&param.bDoCSound,sizeof(int),
		    "csound","enable/disable sound speed outputs = -csound");
	param.nSmooth = 32;
	prmAddParam(prm, "nSmooth", paramInt, &param.nSmooth,
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
	param.dFracFastGas = 0.1;
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
	iGasModel = GASMODEL_UNSET;	/* Deprecated in for backwards
					   compatibility */
	prmAddParam(prm,"iGasModel", paramInt, &iGasModel, sizeof(int),
		    "GasModel", "<Gas model employed> = 0 (Adiabatic)");
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

	param.bStarForm = 0;
	prmAddParam(prm,"bStarForm",paramBool,&param.bStarForm,sizeof(int),
		    "stfm","<Star Forming> = 0");

	param.stfm = new Stfm();
	param.stfm->AddParams(prm);

	param.bFeedback = 0;
	prmAddParam(prm,"bFeedBack",paramBool,&param.bFeedback,sizeof(int),
		    "fdbk","<Stars provide feedback> = 0");

	param.feedback = new Fdbk();
	param.feedback->AddParams(prm);

	param.iRandomSeed = 1;
	prmAddParam(prm,"iRandomSeed", paramInt, &param.iRandomSeed,
		    sizeof(int), "iRand", "<Feedback random Seed> = 1");
	

	
	//
	// Output parameters
	//
	printBinaryAcc=0;
	prmAddParam(prm, "bPrintBinary", paramBool, &printBinaryAcc,
		    sizeof(int),"z", "Print accelerations in Binary");
	param.bStandard = 1;
	prmAddParam(prm, "bStandard", paramBool, &param.bStandard,sizeof(int),
		    "std", "output in standard TIPSY binary format (IGNORED)");
        param.bDoublePos = 0;
        prmAddParam(prm, "bDoublePos", paramBool, &param.bDoublePos,
                    sizeof(int), "dp",
                    "input/output double precision positions = -dp");
        param.bDoubleVel = 0;
        prmAddParam(prm,"bDoubleVel", paramBool, &param.bDoubleVel,sizeof(int),
                    "dv", "input/output double precision velocities = -dv");
	param.bOverwrite = 1;
	prmAddParam(prm, "bOverwrite", paramBool, &param.bOverwrite,sizeof(int),
		    "overwrite", "overwrite outputs (IGNORED)");
	param.bParaRead = 1;
	prmAddParam(prm, "bParaRead", paramBool, &param.bParaRead,sizeof(int),
		    "par", "enable/disable parallel reading of files (IGNORED)");
	param.bParaWrite = 1;
	prmAddParam(prm, "bParaWrite", paramBool, &param.bParaWrite,sizeof(int),
		    "paw", "enable/disable parallel writing of files");
	param.nIOProcessor = 0;
	prmAddParam(prm,"nIOProcessor",paramInt,&param.nIOProcessor,
		    sizeof(int), "npio",
		    "number of simultaneous I/O processors = 0 (all)");
	param.achInFile[0] = '\0';
	prmAddParam(prm,"achInFile",paramString,param.achInFile,
		    256, "I", "input file name (or base file name)");
	strcpy(param.achOutName,"pargrav");
	prmAddParam(prm,"achOutName",paramString,param.achOutName,256,"o",
				"output name for snapshots and logfile");
	param.dExtraStore = 0.01;
	prmAddParam(prm,"dExtraStore",paramDouble,&param.dExtraStore,
		    sizeof(double), "estore",
		    "Extra memory for new particles");
	param.dMaxBalance = 1e10;
	prmAddParam(prm,"dMaxBalance",paramDouble,&param.dMaxBalance,
		    sizeof(double), "maxbal",
		    "Maximum piece ratio for load balancing");
	
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
	param.cacheLineDepth=4;
	prmAddParam(prm, "nCacheDepth", paramInt, &param.cacheLineDepth,
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
	param.dFracNoDomainDecomp = 0.0;
	prmAddParam(prm, "dFracNoDomainDecomp", paramDouble,
		    &param.dFracNoDomainDecomp, sizeof(double),"fndd",
		    "Fraction of active particles for no new DD = 0.0");
        lbcomm_cutoff_msgs = 1;
	prmAddParam(prm, "lbcommCutoffMsgs", paramInt, &lbcomm_cutoff_msgs,
		    sizeof(int),"lbcommcut", "Cutoff for communication recording (IGNORED)");
	param.bConcurrentSph = 0;
	prmAddParam(prm, "bConcurrentSph", paramBool, &param.bConcurrentSph,
		    sizeof(int),"consph", "Enable SPH running concurrently with Gravity");

#ifdef PUSH_GRAVITY
        param.dFracPushParticles = 0.0;
	prmAddParam(prm, "dFracPush", paramDouble,
		    &param.dFracPushParticles, sizeof(double),"fPush",
		    "Maximum proportion of active to total particles for push-based force evaluation = 0.0");
#endif


#ifdef SELECTIVE_TRACING
        /* tracing control */ 
        monitorRung = 12;
        prmAddParam(prm, "iTraceRung", paramInt, &monitorRung,
              sizeof(int),"traceRung", "Gravity starting rung to trace selectively");

        monitorStart = 0;
        prmAddParam(prm, "iTraceStart", paramInt, &monitorStart,
              sizeof(int),"traceStart", "When to start selective tracing");

        numTraceIterations = 3;
        prmAddParam(prm, "iTraceFor", paramInt, &numTraceIterations,
              sizeof(int),"traceFor", "Trace this many instances of the selected rungs");

        numSkipIterations = 400;
        prmAddParam(prm, "iTraceSkip", paramInt, &numSkipIterations,
              sizeof(int),"traceSkip", "Skip tracing for these many iterations");

        numMaxTrace = 5;
        prmAddParam(prm, "iTraceMax", paramInt, &numMaxTrace,
              sizeof(int),"traceMax", "Max. num. iterations traced");


        traceIteration = 0;
        traceState = TraceNormal;
        projectionsOn = false;
#endif

        numInitDecompBins = (1<<11);
        prmAddParam(prm, "iInitDecompBins", paramInt, &numInitDecompBins,
              sizeof(int),"initDecompBins", "Number of bins to use for the first iteration of every Oct decomposition step");

        octRefineLevel = 1;
        prmAddParam(prm, "iOctRefineLevel", paramInt, &octRefineLevel,
                    sizeof(int),"octRefineLevel", "Binary logarithm of the number of sub-bins a bin is split into for Oct decomposition (e.g. octRefineLevel 3 splits into 8 sub-bins) (default: 1)");

        doDumpLB = false;
        prmAddParam(prm, "bdoDumpLB", paramBool, &doDumpLB,
              sizeof(bool),"doDumpLB", "Should Orb3dLB dump LB database to text file and stop?");

        lbDumpIteration = 0;
        prmAddParam(prm, "ilbDumpIteration", paramInt, &lbDumpIteration,
              sizeof(int),"lbDumpIteration", "Load balancing iteration for which to dump database");

        doSimulateLB = false;
        prmAddParam(prm, "bDoSimulateLB", paramBool, &doSimulateLB,
              sizeof(bool),"doSimulateLB", "Should Orb3dLB simulate LB decisions from dumped text file and stop?");

        CkAssert(!(doDumpLB && doSimulateLB));

    
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

        int processSimfile = 1;
	if(!prmArgProc(prm,m->argc,m->argv, processSimfile)) {
	    CkExit();
        }

        // ensure that number of bins is a power of 2
        if (numInitDecompBins > numTreePieces) {
          unsigned int numBins = 1;
          while (numBins < numTreePieces) {
            numBins <<= 1;
          }
          numInitDecompBins = numBins;
        }
	
	if(bVDetails && !verbosity)
	    verbosity = 1;
	
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
	if(!prmSpecified(prm, "dTheta2")) {
	    param.dTheta2 = param.dTheta;
	    }
	/* set readonly/global variables */
	theta = param.dTheta;
        thetaMono = theta*theta*theta*theta;
	dExtraStore = param.dExtraStore;
	dMaxBalance = param.dMaxBalance;
	_cacheLineDepth = param.cacheLineDepth;
	nIOProcessor = param.nIOProcessor;
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
	if(param.iOutInterval <= 0) {
	    ckerr << "WARNING: ";
	    ckerr << "iOutInterval is non-positive; setting to 1."
		  << endl;
	    param.iOutInterval = 1;
	    }
	    
	if(param.iLogInterval <= 0) {
	    ckerr << "WARNING: ";
	    ckerr << "iLogInterval is non-positive; setting to 1."
		  << endl;
	    param.iLogInterval = 1;
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
	
        dTimeOld = 0.0; // Reset if this is restarting from an output.

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
	if (prmSpecified(prm,"bViscosityLimiter")) {
	    if (!param.bViscosityLimiter) param.iViscosityLimiter=0;
	    }

	if(prmSpecified(prm, "iGasModel")) {
	    switch(iGasModel) {
	    case GASMODEL_ADIABATIC:
		param.bGasAdiabatic = 1;	
		break;
	    case GASMODEL_ISOTHERMAL:
		param.bGasIsothermal = 1;
		break;
	    case GASMODEL_COOLING:
		param.bGasCooling = 1;
		break;
		}
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
	if((param.bFeedback || param.bStarForm) && !param.bDoGas) {
	    ckerr << "WARNING: star formation set without enabling SPH" << endl;
	    ckerr << "Enabling SPH" << endl;
	    param.bDoGas = 1;
	    }
	if(param.bDoGas && !(param.bGasCooling || param.bGasAdiabatic
			     || param.bGasIsothermal)) {
	    ckerr << "Defaulting to Adiabatic Gas Model." << endl;
	    param.bGasAdiabatic = 1;
	    }
	if(!param.bDoGas)
		param.bSphStep = 0;
#include "physconst.h"
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

        if (domainDecomposition == SFC_peano_dec) peanoKey = 3;
        if (domainDecomposition == SFC_peano_dec_2D) peanoKey = 2;
        if (domainDecomposition == SFC_peano_dec_3D) peanoKey = 3;

	// hardcoding some parameters, later may be full options
	if(domainDecomposition==ORB_dec || domainDecomposition==ORB_space_dec){
	    useTree = Binary_ORB;
	    // CkAbort("ORB decomposition known to be bad and not implemented");
	    }
	else { useTree = Binary_Oct; }

#define xstr(s) str(s)
#define str(s) #s
	ckerr << "ChaNGa version " << xstr(NBODY_PACKAGE_VERSION) << ", commit "
              << Cha_CommitID << endl;
#undef str
#undef xstr
        ckerr << "Running on " << CkNumPes() << " processors/ "<< CkNumNodes() <<" nodes with " << numTreePieces << " TreePieces" << endl;
        ckout << "Running on " << CkNumPes() << " processors/ "<< CkNumNodes() <<" nodes with " << numTreePieces << " TreePieces" << endl;

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
          case ORB_space_dec:
            ckerr << "Domain decomposition...ORB space" << endl;
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
#ifdef ROUND_ROBIN_WITH_OCT_DECOMP
	if (domainDecomposition == Oct_dec) {
	  CProxy_RRMap myMap=CProxy_RRMap::ckNew(); 
	  opts.setMap(myMap);
	} else {
#endif
#ifdef DEFAULT_ARRAY_MAP
          CProxy_DefaultArrayMap myMap = CProxy_DefaultArrayMap::ckNew();
#else
	  CProxy_BlockMap myMap = CProxy_BlockMap::ckNew(); 
#endif
	  opts.setMap(myMap);
#ifdef ROUND_ROBIN_WITH_OCT_DECOMP
	}
#endif

	CProxy_TreePiece pieces = CProxy_TreePiece::ckNew(opts);

        prjgrp = CProxy_ProjectionsControl::ckNew();

	treeProxy = pieces;
#ifdef REDUCTION_HELPER
        reductionHelperProxy = CProxy_ReductionHelper::ckNew();
#endif

	opts.bindTo(treeProxy);
	lvProxy = CProxy_LvArray::ckNew(opts);
	// Create an array for the smooth reductions
	smoothProxy = CProxy_LvArray::ckNew(opts);
	// Create an array for the gravity reductions
	gravityProxy = CProxy_LvArray::ckNew(opts);

#ifdef PUSH_GRAVITY
        ckMulticastGrpId = CProxy_CkMulticastMgr::ckNew();
#endif

        peTreeMergerProxy = CProxy_PETreeMerger::ckNew();
        dfDataProxy = CProxy_DumpFrameData::ckNew();
	
	// create CacheManagers
	// Gravity particles
	cacheGravPart = CProxy_CkCacheManager<KeyType>::ckNew(cacheSize, pieces.ckLocMgr()->getGroupID());
	// Smooth particles
	cacheSmoothPart = CProxy_CkCacheManager<KeyType>::ckNew(cacheSize, pieces.ckLocMgr()->getGroupID());
	// Nodes
	cacheNode = CProxy_CkCacheManager<KeyType>::ckNew(cacheSize, pieces.ckLocMgr()->getGroupID());

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

/// 
/// @brief Restart Main constructor
///
/// This is only called when restarting from a checkpoint.
///
Main::Main(CkMigrateMessage* m) : CBase_Main(m) {
    args = (CkArgMsg *)malloc(sizeof(*args));
    args->argv = CmiCopyArgs(((CkArgMsg *) m)->argv);
    mainChare = thishandle;
    bIsRestarting = 1;
    CkPrintf("Main(CkMigrateMessage) called\n");
    sorter = CProxy_Sorter::ckNew(0);
    }


/// @brief entry method to cleanly shutdown.  Only used for debugging.
void Main::niceExit() {
  static unsigned int count = 0;
  if (++count == numTreePieces) CkExit();
}

/// @brief determine start time of simulation
///
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

		if (!prmArgSpecified(prm,"nSteps") &&
		    prmArgSpecified(prm,"dDelta")) {

		    param.nSteps = (int)ceil((tTo-dTime)/param.dDelta)+param.iStartStep;
		    }
		else if (!prmArgSpecified(prm,"dDelta") &&
			 prmArgSpecified(prm,"nSteps")) {
		    if(param.nSteps != 0)
			param.dDelta =
				(tTo-dTime)/(param.nSteps - param.iStartStep);
		    else
			param.dDelta = 0.0;
		    }
		else if (!prmSpecified(prm,"nSteps") &&
			 prmFileSpecified(prm,"dDelta")) {
		    param.nSteps = (int)ceil((tTo-dTime)/param.dDelta) + param.iStartStep;
		    }
		else if (!prmSpecified(prm,"dDelta") &&
			 prmFileSpecified(prm,"nSteps")) {
		    if(param.nSteps != 0)
			param.dDelta =	(tTo-dTime)/(param.nSteps
							 - param.iStartStep);
		    else
			param.dDelta = 0.0;
		    }
		}
	    else {
		tTo = dTime + (param.nSteps-param.iStartStep)*param.dDelta;
		aTo = csmTime2Exp(param.csm,tTo);
		if (verbosity > 0)
		    ckout << "Simulation to Time:" << tTo << " Redshift:"
			  << 1.0/aTo - 1.0 << " Expansion factor:"
			  << aTo << endl;
		}
	    // convert to canonical comoving coordinates.
	    treeProxy.velScale(dAStart*dAStart, CkCallbackResumeThread());
	    }
	else {  // not Comove
	    dTime = treeProxy[0].ckLocal()->dStartTime; 
	    if(verbosity > 0)
		ckout << "Input file, Time:" << dTime << endl;
	    double tTo = dTime + (param.nSteps-param.iStartStep)*param.dDelta;
	    if (verbosity > 0) {
		ckout << "Simulation to Time:" << tTo << endl;
		}
	    }
    }

/// @brief Return true if we need to write an output
/// 
/// Advances iOut attribute, therefore this can only be called once
/// per timestep.
///
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

///
/// @brief Read in desired output times and reshifts from a file
///
/// Fills in the vdOutTime vector by reading the .red file
///
void Main::getOutTimes()
{
    FILE *fp;
    int ret;
    double z,a,n,t;
    char achIn[80];
    string achFileName = string(param.achOutName) + ".red";
	
    fp = fopen(achFileName.c_str(),"r");
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

/// @brief wait for gravity in the case of concurrent SPH
inline void Main::waitForGravity(const CkCallback &cb, double startTime) 
{
    if(param.bConcurrentSph && param.bDoGravity) {
#ifdef PUSH_GRAVITY
      if(bDoPush){
        CkWaitQD();
      }
      else{
#endif
        CkFreeMsg(cb.thread_delay());
#ifdef PUSH_GRAVITY
      }
#endif
        CkPrintf("Calculating gravity and SPH took %g seconds.\n",
		 CkWallTimer()-startTime);
#ifdef SELECTIVE_TRACING
        turnProjectionsOff();
#endif
	}
    }

///
/// @brief Take one base timestep of the simulation.
/// @param iStep The current step number.
///
/// This method implements the standard "Kick Drift Kick" (Quinn et al
/// 1997) hierarchical timestepping algorithm.  It assumes that the
/// forces for the first opening kick have already been calculated.
///

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
 
      if(verbosity)
	  memoryStats();
      // Opening Kick
      if(verbosity)
	  ckout << "Kick Open:" << endl;
      double dKickFac[MAXRUNG+1];
      double duKick[MAXRUNG+1];
      double dStartTime[MAXRUNG+1];
      for (int iRung=0; iRung<=MAXRUNG; iRung++) {
        duKick[iRung]=dKickFac[iRung]=0;
      }
      for(int iRung = activeRung; iRung <= nextMaxRung; iRung++) {
        double dTimeSub = RungToDt(param.dDelta, iRung);
	if(verbosity) {
	    CkPrintf(" Rung %d: %g\n", iRung, 0.5*dTimeSub);
	    }
        duKick[iRung] = 0.5*dTimeSub;
        dKickFac[iRung] = csmComoveKickFac(param.csm, dTime, 0.5*dTimeSub);
	dStartTime[iRung] = dTime;
        }
      if(verbosity > 1)
	  memoryStats();
      if(param.bDoGas) {
	  double startTime = CkWallTimer();
	  if(verbosity)
	      CkPrintf("uDot update: Rung %d ... ", activeRung);
	  double z = 1.0/csmTime2Exp(param.csm,dTime) - 1.0;
	  if(param.bGasCooling)
	      dMProxy.CoolingSetTime(z, dTime, CkCallbackResumeThread());
	  treeProxy.updateuDot(activeRung, duKick, dStartTime,
			       param.bGasCooling, 1, 1,
			       CkCallbackResumeThread());
	  if(verbosity)
	      CkPrintf("took %g seconds.\n", CkWallTimer() - startTime);
	  }
      treeProxy.kick(activeRung, dKickFac, 0, param.bDoGas,
		     param.bGasIsothermal, duKick, CkCallbackResumeThread());

      if(verbosity > 1)
	  memoryStats();
      // Dump frame may require a smaller step
      int driftRung = nextMaxRungIncDF(nextMaxRung);
      // Get number of drift steps from driftRung
      int driftSteps = 1;
      while(driftRung > nextMaxRung) {
	  driftRung--;
	  driftSteps <<= 1;
	  }
      driftRung = nextMaxRungIncDF(nextMaxRung);
      if(param.bDoGas && (driftRung == nextMaxRung)) {
	  driftRung++;
	  driftSteps <<=1;
	  }
      
      double dTimeSub = RungToDt(param.dDelta, driftRung);
      // Drift of smallest step
      for(int iSub = 0; iSub < driftSteps; iSub++) 
	  {
	      if(verbosity)
		  CkPrintf("Drift: Rung %d Delta %g\n", driftRung, dTimeSub);

	      // Only effective if growmass parameters have been set.
	      growMass(dTime, dTimeSub);
	      // Are the GrowMass particles locked in place?
	      int nGrowMassDrift = param.nGrowMass;
	      if(param.bDynGrowMass) nGrowMassDrift = 0;
	      
	      double dDriftFac = csmComoveDriftFac(param.csm, dTime, dTimeSub);
	      double dKickFac = csmComoveKickFac(param.csm, dTime, dTimeSub);
	      bool buildTree = (iSub + 1 == driftSteps);
	      treeProxy.drift(dDriftFac, param.bDoGas, param.bGasIsothermal,
			      dKickFac, dTimeSub, nGrowMassDrift, buildTree,
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
	      // Update uDot at the 1/2 way point.
	      if(param.bDoGas && (iSub+1 == driftSteps >> 1)) {
		  double startTime = CkWallTimer();
		  if(verbosity)
		      CkPrintf("uDot' update: Rung %d ... ", nextMaxRung);
		  double z = 1.0/csmTime2Exp(param.csm,dTime) - 1.0;
		  if(param.bGasCooling)
		      dMProxy.CoolingSetTime(z, dTime,
					     CkCallbackResumeThread());
		  dStartTime[nextMaxRung] = dTime; // only updating this rung
		  treeProxy.updateuDot(nextMaxRung, duKick, dStartTime,
				       param.bGasCooling, 0, 0,
				       CkCallbackResumeThread());
		  if(verbosity)
		      CkPrintf("took %g seconds.\n", CkWallTimer() - startTime);
		  }
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
    /*
     * Form stars at user defined intervals
     */
    double dTimeSF = RungToDt(param.dDelta, nextMaxRung);
    if((param.bStarForm || param.bFeedback)
       && param.stfm->isStarFormRung(activeRung)) {
        double startTime = CkWallTimer();
        CkPrintf("Domain decomposition for star formation/feedback... ");
        sorter.startSorting(dataManagerID, ddTolerance,
                            CkCallbackResumeThread(), true);
        CkPrintf("took %g seconds.\n", CkWallTimer()-startTime);
        CkPrintf("Load balancer for star formation/feedback... ");
        startTime = CkWallTimer();
        treeProxy.startlb(CkCallbackResumeThread(), PHASE_FEEDBACK);
        CkPrintf("took %g seconds.\n", CkWallTimer()-startTime);
        if(param.bStarForm)
            FormStars(dTime, max(dTimeSF, param.stfm->dDeltaStarForm));
        if(param.bFeedback) 
            StellarFeedback(dTime, max(dTimeSF, param.stfm->dDeltaStarForm));
        }

    ckout << "\nStep: " << (iStep + ((double) currentStep)/MAXSUBSTEPS)
          << " Time: " << dTime
          << " Rungs " << activeRung << " to "
          << nextMaxRung << ". ";

    countActive(activeRung);
    
    if(verbosity > 1)
	memoryStats();

    /***** Resorting of particles and Domain Decomposition *****/
    CkPrintf("Domain decomposition ... \n");
    double startTime;
    bool bDoDD = param.dFracNoDomainDecomp*nTotalParticles < nActiveGrav;

    startTime = CkWallTimer();
    if (iStep % 5 != 0) {
      bDoDD = false;
    }
    if (bDoDD) {
      sorter.startSorting(dataManagerID, ddTolerance,
          CkCallbackResumeThread(), bDoDD);
    } else {
      treeProxy.migrateAndUnshuffle(CkCallbackResumeThread());
      CkPrintf("Doing DD with migrateAndUnshuffle\n");
    }
    /*
    ckout << " took " << (CkWallTimer() - startTime) << " seconds."
          << endl;
          */
    CkPrintf("DD ending at %g took total %g seconds.\n", CkWallTimer(), CkWallTimer()-startTime);

    if(!bDoDD)
	CkPrintf("Skipped DD\n");

    if(verbosity > 1)
	memoryStats();
    /********* Load balancer ********/
    //ckout << "Load balancer ...";
    CkPrintf("Load balancer ... \n");
    startTime = CkWallTimer();
    treeProxy.startlb(CkCallbackResumeThread(), activeRung);
    /*
    ckout << " took "<<(CkWallTimer() - startTime) << " seconds."
	     << endl;
             */
    CkPrintf("LB ending at %g took total %g seconds.\n", CkWallTimer(), CkWallTimer()-startTime);

    if(verbosity > 1)
	memoryStats();


#ifdef PUSH_GRAVITY
    bool bDoPush = param.dFracPushParticles*nTotalParticles > nActiveGrav;
    if(bDoPush) CkPrintf("[main] fracActive %f PUSH_GRAVITY\n", 1.0*nActiveGrav/nTotalParticles);
#endif

    /******** Tree Build *******/
    //ckout << "Building trees ...";
    CkPrintf("Building trees ... \n");
    startTime = CkWallTimer();
#ifdef PUSH_GRAVITY
    treeProxy.buildTree(bucketSize, CkCallbackResumeThread(),!bDoPush);
#else
    //treeProxy.buildTree(bucketSize, CkCallbackResumeThread());
    treeProxy.buildTreeOverlap(bucketSize, CkCallbackResumeThread());
#endif
    CkPrintf("TB ending at %g took total %g seconds.\n", CkWallTimer(), CkWallTimer()-startTime);

    CkCallback cbGravity(CkCallback::resumeThread);

    if(verbosity > 1)
	memoryStats();
    if(param.bDoGravity) {
	updateSoft();
	if(param.csm->bComove) {
	    double a = csmTime2Exp(param.csm,dTime);
	    if (a >= param.daSwitchTheta) theta = param.dTheta2; 
	    }
	/******** Force Computation ********/
	//ckout << "Calculating gravity (tree bucket, theta = " << theta
	//      << ") ...";
#ifdef SELECTIVE_TRACING
        turnProjectionsOn(activeRung);
#endif

        CkPrintf("Calculating gravity (tree bucket, theta = %f) ... \n", theta);
	startTime = CkWallTimer();
	if(param.bConcurrentSph) {
	    ckout << endl;

#ifdef PUSH_GRAVITY
            if(bDoPush){ 
              treeProxy.startPushGravity(activeRung, theta);
            }
	    else{ 
#endif
              treeProxy.startGravity(activeRung, theta, cbGravity);
#ifdef PUSH_GRAVITY
            }
#endif


#ifdef CUDA_INSTRUMENT_WRS
            dMProxy.clearInstrument(cbGravity);
#endif
	    }
	else {
#ifdef PUSH_GRAVITY
	    if(bDoPush){
              treeProxy.startPushGravity(activeRung, theta);
              CkWaitQD();
            }
	    else{
#endif
              treeProxy.startGravity(activeRung, theta, CkCallbackResumeThread());
#ifdef PUSH_GRAVITY
            }
#endif

#ifdef CUDA_INSTRUMENT_WRS
            dMProxy.clearInstrument(CkCallbackResumeThread());
#endif
	    //ckout << " took " << (CkWallTimer() - startTime) << " seconds."
	    //	  << endl;
            CkPrintf("CG ending at %g took total %g seconds\n", CkWallTimer(), CkWallTimer()-startTime);
#ifdef SELECTIVE_TRACING
            turnProjectionsOff();
#endif
            if(verbosity)
                memoryStatsCache();
	    }
        }
    else {
	treeProxy.initAccel(activeRung, CkCallbackResumeThread());
	}
    if(verbosity > 1)
	memoryStats();
    if(param.bDoGas) {
        doSph(activeRung);
        if(verbosity && !param.bConcurrentSph)
            memoryStatsCache();
        }
    
    if(!param.bStaticTest) {
      // Closing Kick
      double dKickFac[MAXRUNG+1];
      double duKick[MAXRUNG+1];
      double dStartTime[MAXRUNG+1];
      for (int iRung=0; iRung<=MAXRUNG; iRung++) {
        duKick[iRung]=dKickFac[iRung]=0;
      }
      if(verbosity > 1)
	  memoryStats();
      if(verbosity)
	  ckout << "Kick Close:" << endl;
      for(int iRung = activeRung; iRung <= nextMaxRung; iRung++) {
        double dTimeSub = RungToDt(param.dDelta, iRung);
	if(verbosity) {
	    CkPrintf(" Rung %d: %g\n", iRung, 0.5*dTimeSub);
	    }
        duKick[iRung] = 0.5*dTimeSub;
	dStartTime[iRung] = dTime - 0.5*dTimeSub;
        dKickFac[iRung] = csmComoveKickFac(param.csm,
                                           dTime - 0.5*dTimeSub,
                                           0.5*dTimeSub);
      }
      if(param.bDoGas) {
	  double startTime = CkWallTimer();
	  if(verbosity)
	      CkPrintf("uDot update: Rung %d ... ", activeRung);
	  double z = 1.0/csmTime2Exp(param.csm,dTime) - 1.0;
	  if(param.bGasCooling)
	      dMProxy.CoolingSetTime(z, dTime, CkCallbackResumeThread());
	  treeProxy.updateuDot(activeRung, duKick, dStartTime,
			       param.bGasCooling, 1, 1,
			       CkCallbackResumeThread());
	  if(verbosity)
	      CkPrintf("took %g seconds.\n", CkWallTimer() - startTime);
	  }
      waitForGravity(cbGravity, startTime);
      treeProxy.kick(activeRung, dKickFac, 1, param.bDoGas,
		     param.bGasIsothermal, duKick, CkCallbackResumeThread());
      // 1/2 step uDot update
      if(activeRung > 0 && param.bDoGas) {
	  double startTime = CkWallTimer();
	  if(verbosity)
	      CkPrintf("uDot' update: Rung %d ... ", activeRung-1);
	  duKick[activeRung-1] = 0.5*RungToDt(param.dDelta, activeRung-1);
	  dStartTime[activeRung-1] = dTime;
	  treeProxy.updateuDot(activeRung-1, duKick, dStartTime,
			       param.bGasCooling, 0, 0,
			       CkCallbackResumeThread());
	  if(verbosity)
	      CkPrintf("took %g seconds.\n", CkWallTimer() - startTime);
	  }
    }
    else
	waitForGravity(cbGravity, startTime);

#if COSMO_STATS > 0
    /********* TreePiece Statistics ********/
    ckerr << "Total statistics iteration " << iStep << "/";
    int printingStep = currentStep;
    while (printingStep) {
      ckerr << ((printingStep & MAXSUBSTEPS) ? "1" : "0");
      printingStep = (printingStep & ~MAXSUBSTEPS) << 1;
    }
    //ckerr << " (rungs " << activeRung << "-" << nextMaxRung << ")" << ":" << endl;
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

    treeProxy.finishNodeCache(CkCallbackResumeThread());

#ifdef CHECK_TIME_WITHIN_BIGSTEP
    if(param.iWallRunTime > 0 && ((CkWallTimer()-wallTimeStart) > param.iWallRunTime*60.)){
      CkPrintf("wall time %g exceeded limit %g within advancestep\n", CkWallTimer()-wallTimeStart, param.iWallRunTime*60.);
      CkExit();
    }
#endif

		
  }
}
    
///
/// @brief Load particles into pieces
///
/// Reads the particle data in from a file.  Since the full
/// information about the run (including starting time) isn't known
/// until the particles are loaded, this routine also completes the
/// specification of the run details and writes out the log file
/// entry.  It concludes by calling initialForces()

void Main::setupICs() {
  double startTime;

  treeProxy.setPeriodic(param.nReplicas, param.vPeriod, param.bEwald,
			param.dEwCut, param.dEwhCut, param.bPeriodic);

  /******** Particles Loading ********/
  CkPrintf("Loading particles ...");
  startTime = CkWallTimer();

  try {
      struct stat s;
      int err = stat(basefilename.c_str(), &s);
      if(err == -1) {
          ckerr << "File error: " << basefilename.c_str() << endl;
          CkExit();
          return;
          }
      double dTuFac = param.dGasConst/(param.dConstGamma-1)/param.dMeanMolWeight;
      if(S_ISDIR(s.st_mode)) {
          // Assume its an NChilada directory
          treeProxy.loadNChilada(basefilename, dTuFac,
              CkCallbackResumeThread());
          }
      else {
          // Assume Tipsy format
          // Try loading Tipsy format; first just grab the header.
          CkPrintf(" trying Tipsy ...");
	    
          bool bInDPos = false;
          bool bInDVel = false;
  
          Tipsy::PartialTipsyFile ptf(basefilename, 0, 0);
          if(!ptf.loadedSuccessfully()) {
              ckerr << endl << "Couldn't load the tipsy file \""
                    << basefilename.c_str()
                    << "\". Maybe it's not a tipsy file?" << endl;
              CkExit();
              return;
              }
          bInDPos = ptf.isDoublePos();
          bInDVel = ptf.isDoubleVel();
          if(bInDPos) CkPrintf("assumed double positions...");
          if(bInDVel) CkPrintf("assumed double velocities...");

          treeProxy.loadTipsy(basefilename, dTuFac, bInDPos, bInDVel,
              CkCallbackResumeThread());
          }
  }
  catch (std::ios_base::failure e) {
    ckerr << "File read: " << basefilename.c_str() << ": " << e.what()
          << endl;
    CkExit();
    return;
  }

  ckout << " took " << (CkWallTimer() - startTime) << " seconds."
        << endl;
	    
  nTotalParticles = treeProxy[0].ckLocal()->nTotalParticles;
  nTotalSPH = treeProxy[0].ckLocal()->nTotalSPH;
  nTotalDark = treeProxy[0].ckLocal()->nTotalDark;
  nTotalStar = treeProxy[0].ckLocal()->nTotalStar;
  nMaxOrderGas = nTotalSPH - 1;
  nMaxOrderDark = nTotalSPH + nTotalDark - 1;
  nMaxOrder = nTotalParticles - 1;
  nActiveGrav = nTotalParticles;
  nActiveSPH = nTotalSPH;
  ckout << "N: " << nTotalParticles << endl;
  if(particlesPerChare == 0)
      particlesPerChare = nTotalParticles/numTreePieces;
  
  if(nTotalSPH == 0 && param.bDoGas) {
      ckerr << "WARNING: no SPH particles and bDoGas is set\n";
      param.bDoGas = 0;
      }
  if(nTotalSPH > 0 && !param.bDoGas) {
      if(prmSpecified(prm, "bDoGas"))
	  ckerr << "WARNING: SPH particles present and bDoGas is set off\n";
      else
	  param.bDoGas = 1;
      }
  getStartTime();
  if(param.nSteps > 0) getOutTimes();
  for(iOut = 0; iOut < vdOutTime.size(); iOut++) {
      if(dTime < vdOutTime[iOut]) break;
      }
	
  if(param.iStartStep) bIsRestarting = true;

  if(param.bGasCooling || param.bStarForm) {
      initCooling();
      if(param.iStartStep) restartGas();
      }
  
  if(param.bStarForm || param.bFeedback)
      param.stfm->CheckParams(prm, param);

  if(param.bStarForm)
      initStarLog();
	
  if(param.bFeedback)
      param.feedback->CheckParams(prm, param);
  else
      param.feedback->NullFeedback();

  string achLogFileName = string(param.achOutName) + ".log";
  ofstream ofsLog;
  if(bIsRestarting)
      ofsLog.open(achLogFileName.c_str(), ios_base::app);
  else
      ofsLog.open(achLogFileName.c_str(), ios_base::trunc);
  if(!ofsLog)
      CkAbort("Error opening log file.");
      
#define xstr(s) str(s)
#define str(s) #s
  ofsLog << "# Starting ChaNGa version " << xstr(NBODY_PACKAGE_VERSION)
         << " commit " << Cha_CommitID << endl;
#undef str
#undef xstr
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
#ifdef EPSACCH
  ofsLog << " EPSACCH";
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
#ifdef BIGKEYS
  ofsLog << " BIGKEYS";
#endif
  ofsLog << endl;
  ofsLog << "# Key sizes: " << sizeof(KeyType) << " bytes particle "
         << sizeof(NodeKey) << " bytes node" << endl;

  // Print out load balance information
  LBDatabase *lbdb = LBDatabaseObj();
  int nlbs = lbdb->getNLoadBalancers(); 
  if(nlbs == 0) {
      ofsLog << "# No load balancer in use" << endl;
      }
  else {
      int ilb;
      BaseLB **lbs = lbdb->getLoadBalancers();
      ofsLog << "# Load balancers:";
      for(ilb = 0; ilb < nlbs; ilb++){
	  ofsLog << " " << lbs[ilb]->lbName();
	  }
      ofsLog << endl;
      }
  
  ofsLog.close();
	
  prmLogParam(prm, achLogFileName.c_str());
	
  ofsLog.open(achLogFileName.c_str(), ios_base::app);
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
  if(!ofsLog)
      CkAbort("Error closing log file");

  if(prmSpecified(prm,"dSoft")) {
    ckout << "Set Softening...\n";
    treeProxy.setSoft(param.dSoft, CkCallbackResumeThread());
  }
	
// for periodic, puts all particles within the boundary
// Also assigns keys and sorts.
  treeProxy.drift(0.0, 0, 0, 0.0, 0.0, 0, true, CkCallbackResumeThread());

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

/// @brief Callback to restart simulation after a checkpoint or after writing
/// a checkpoint.
///
/// A restart looks like we've just finished writing a checkpoint.  We
/// can tell the difference by the bIsRestarting flag set in the Main
/// CkMigrate constructor.  In that case we do some simple parameter
/// parsing and go to InitialForces.  Otherwise we return to the
/// doSimulation() loop.

void
Main::restart() 
{
    if(bIsRestarting) {
	dSimStartTime = CkWallTimer();
#define xstr(s) str(s)
#define str(s) #s
	ckout << "ChaNGa version " << xstr(NBODY_PACKAGE_VERSION)
              << ", commit " << Cha_CommitID << endl;
#undef str
#undef xstr
	ckout << "Restarting at " << param.iStartStep << endl;
	string achLogFileName = string(param.achOutName) + ".log";
	ofstream ofsLog(achLogFileName.c_str(), ios_base::app);
	ofsLog << "# ReStarting ChaNGa commit " << Cha_CommitID << endl;
	ofsLog << "#";
	for (int i = 0; i < CmiGetArgc(args->argv); i++)
	    ofsLog << " " << args->argv[i];
	ofsLog << endl;
        ofsLog << "# Running on " << CkNumPes() << " processors/ "
               << CkNumNodes() << " nodes" << endl;
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
        prmAddParam(prm, "dEta", paramDouble, &param.dEta,
                    sizeof(double),"eta", "Time integration accuracy");
	prmAddParam(prm,"nIOProcessor",paramInt,&param.nIOProcessor,
		    sizeof(int), "npio",
		    "number of simultaneous I/O processors = 0 (all)");
	prmAddParam(prm, "nBucket", paramInt, &bucketSize,
		    sizeof(int),"b", "Particles per Bucket (default: 12)");
	prmAddParam(prm, "nCacheDepth", paramInt, &param.cacheLineDepth,
		    sizeof(int),"d", "Cache Line Depth (default: 4)");
	prmAddParam(prm, "bConcurrentSph", paramBool, &param.bConcurrentSph,
		    sizeof(int),"consph",
		    "Enable SPH running concurrently with Gravity");
	prmAddParam(prm, "bFastGas", paramBool, &param.bFastGas,
		    sizeof(int),"Fgas", "Fast Gas Method");
	prmAddParam(prm,"dFracFastGas",paramDouble,&param.dFracFastGas,
		    sizeof(double),"ffg",
		    "<Fraction of Active Particles for Fast Gas>");
	prmAddParam(prm,"ddHonHLimit",paramDouble,&param.ddHonHLimit,
		    sizeof(double),"dhonh", "<|dH|/H Limiter> = 0.1");
	prmAddParam(prm, "iOutInterval", paramInt, &param.iOutInterval,
		    sizeof(int),"oi", "Output Interval");
	prmAddParam(prm, "iBinaryOutput", paramInt, &param.iBinaryOut,
		    sizeof(int), "binout",
		    "<array outputs 0 ascii, 1 float, 2 double, 3 FLOAT(internal)> = 0");
	prmAddParam(prm, "iCheckInterval", paramInt, &param.iCheckInterval,
		    sizeof(int),"oc", "Checkpoint Interval");
	prmAddParam(prm, "iVerbosity", paramInt, &verbosity,
		    sizeof(int),"v", "Verbosity");
	prmAddParam(prm,"dExtraStore",paramDouble,&param.dExtraStore,
		    sizeof(double), "estore",
		    "Extra memory for new particles");
	prmAddParam(prm,"dMaxBalance",paramDouble,&param.dMaxBalance,
		    sizeof(double), "maxbal",
		    "Maximum piece ratio for load balancing");

        int processSimfile = 0; 
	if(!prmArgProc(prm,CmiGetArgc(args->argv),args->argv,processSimfile)) {
	    CkExit();
	}
	
	dMProxy.resetReadOnly(param, CkCallbackResumeThread());
	treeProxy.drift(0.0, 0, 0, 0.0, 0.0, 0, true, CkCallbackResumeThread());
	if(param.bGasCooling || param.bStarForm) 
	    initCooling();
	if(param.bStarForm)
	    initStarLog();
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

/// For checkpointing, the restart callback needs to be an array entry,
/// so we have a short entry that simply calls the main entry.

void
TreePiece::restart() 
{
    mainChare.restart();
    }

///
/// @brief Initial calculation of forces.
///
/// This is called both when starting or restarting a run.  It
/// concludes by calling doSimulation(), the main simulation loop.
///
void
Main::initialForces()
{
  double startTime;

  // DEBUGGING
  // CkStartQD(CkCallback(CkIndex_TreePiece::quiescence(),treeProxy));

  /***** Initial sorting of particles and Domain Decomposition *****/
  CkPrintf("Initial domain decomposition ... ");

  startTime = CkWallTimer();
  sorter.startSorting(dataManagerID, ddTolerance,
	 	      CkCallbackResumeThread(), true);
  CkPrintf("total %g seconds.\n", CkWallTimer()-startTime);
  /*
  ckout << " took " << (CkWallTimer() - startTime) << " seconds."
        << endl;
        */

#ifdef PUSH_GRAVITY
  treeProxy.findTotalMass(CkCallbackResumeThread());
#endif
  
  // Balance load initially after decomposition
  //ckout << "Initial load balancing ..." << endl;
  CkPrintf("Initial load balancing ... ");
  startTime = CkWallTimer();
  treeProxy.balanceBeforeInitialForces(CkCallbackResumeThread());
  /*
  ckout << " took " << (CkWallTimer() - startTime) << " seconds."
        << endl;
        */
  CkPrintf("took %g seconds.\n", CkWallTimer()-startTime);

  /******** Tree Build *******/
  //ckout << "Building trees ...";
  CkPrintf("Building trees ... ");
  startTime = CkWallTimer();
#ifdef PUSH_GRAVITY
  treeProxy.buildTree(bucketSize, CkCallbackResumeThread(),true);
#else
  treeProxy.buildTree(bucketSize, CkCallbackResumeThread());
#endif
  CkPrintf("took %g seconds.\n", CkWallTimer()-startTime);

  if(verbosity)
      memoryStats();
  
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
      //ckout << "Calculating gravity (theta = " << theta
      //    << ") ...";
#ifdef SELECTIVE_TRACING
      turnProjectionsOn(0);
#endif
      CkPrintf("Calculating gravity (theta = %f) ... \n", theta);
      startTime = CkWallTimer();
      if(param.bConcurrentSph) {

	  treeProxy.startGravity(0, theta, cbGravity);

#ifdef CUDA_INSTRUMENT_WRS
            dMProxy.clearInstrument(cbGravity);
#endif
	  }
      else {
	  treeProxy.startGravity(0, theta, CkCallbackResumeThread());

#ifdef CUDA_INSTRUMENT_WRS
            dMProxy.clearInstrument(CkCallbackResumeThread());
#endif
          CkPrintf("took %g seconds.\n", CkWallTimer()-startTime);
#ifdef SELECTIVE_TRACING
          turnProjectionsOff();
#endif
	  if(verbosity)
	      memoryStatsCache();
	  }
      }
  else {
      treeProxy.initAccel(0, CkCallbackResumeThread());
      }
  if(param.bDoGas) {
      initSph();
      }
  
  if(param.bConcurrentSph && param.bDoGravity) {
      CkFreeMsg(cbGravity.thread_delay());
	//ckout << "Calculating gravity and SPH took "
	//      << (CkWallTimer() - startTime) << " seconds." << endl;
        CkPrintf("Initial gravity and SPH took %f seconds.\n", CkWallTimer()-startTime);
#ifdef SELECTIVE_TRACING
        turnProjectionsOff();
#endif
      }
  
  treeProxy.finishNodeCache(CkCallbackResumeThread());

  // Initial Log entry
  string achLogFileName = string(param.achOutName) + ".log";
  calcEnergy(dTime, CkWallTimer() - startTime, achLogFileName.c_str());

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

///
/// \brief Principal method which does all the coordination of the
/// simulation over timesteps.
///
/// This routine calls advanceBigStep() for each step, logs
/// statistics, determines if output is needed, and halts the
/// simulation when done.
///

void
Main::doSimulation()
{
  double startTime;
  string achLogFileName = string(param.achOutName) + ".log";

#ifdef CHECK_TIME_WITHIN_BIGSTEP
  wallTimeStart = CkWallTimer();
#endif

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
	calcEnergy(dTime, stepTime, achLogFileName.c_str());
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
       && (iStop || (param.iCheckInterval && iStep%param.iCheckInterval == 0))) {
	string achCheckFileName(param.achOutName);
	if(bChkFirst) {
	    achCheckFileName += ".chk0";
	    bChkFirst = 0;
	    }
	else {
	    achCheckFileName += ".chk1";
	    bChkFirst = 1;
	    }
	// The following drift is called because it deletes the tree
	// so it won't be saved on disk.
	treeProxy.drift(0.0, 0, 0, 0.0, 0.0, 0, false, CkCallbackResumeThread());
	treeProxy[0].flushStarLog(CkCallbackResumeThread());
	param.iStartStep = iStep; // update so that restart continues on
	bIsRestarting = 0;
	CkCallback cb(CkIndex_TreePiece::restart(), treeProxy[0]);
	CkStartCheckpoint(achCheckFileName.c_str(), cb);
	return;
    }
    if (iStop) break;

  } // End of main computation loop

  /******** Shutdown process ********/

  if(param.nSteps == 0) {
      string achFile = string(param.achOutName) + ".000000";
      // assign each particle its domain for diagnostic.
      treeProxy.assignDomain(CkCallbackResumeThread());

      if((!param.bDoGas) && param.bDoDensity) {
	  // If gas isn't being calculated, we can do the total
	  // densities before we start the output.
	  ckout << "Calculating total densities ...";
	  DensitySmoothParams pDen(TYPE_GAS|TYPE_DARK|TYPE_STAR, 0);
	  startTime = CkWallTimer();
	  double dfBall2OverSoft2 = 0.0;
	  treeProxy.startSmooth(&pDen, 1, param.nSmooth, dfBall2OverSoft2,
				CkCallbackResumeThread());
	  ckout << " took " << (CkWallTimer() - startTime) << " seconds."
		<< endl;
#if 0
	  // For testing concurrent Sph, we don't want to do the
	  // resmooth.

          ckout << "Recalculating densities ...";
          startTime = CkWallTimer();
	  treeProxy.startIterationReSmooth(&pDen, CkCallbackResumeThread());
          ckout << " took " << (CkWallTimer() - startTime) << " seconds." << endl;
#endif
	  treeProxy.finishNodeCache(CkCallbackResumeThread());
          ckout << "Reordering ...";
          startTime = CkWallTimer();
	  treeProxy.reOrder(nMaxOrder, CkCallbackResumeThread());
          ckout << " took " << (CkWallTimer() - startTime) << " seconds." << endl;
	  ckout << "Outputting densities ...";
	  startTime = CkWallTimer();
	  DenOutputParams pDenOut(achFile, 0, 0.0);
	  treeProxy[0].outputASCII(pDenOut, param.bParaWrite, CkCallbackResumeThread());
	  ckout << " took " << (CkWallTimer() - startTime) << " seconds."
		<< endl;
	  ckout << "Outputting hsmooth ...";
	  HsmOutputParams pHsmOut(achFile, 0, 0.0);
	  treeProxy[0].outputASCII(pHsmOut, param.bParaWrite, CkCallbackResumeThread());
	  }
      else {
	  treeProxy.reOrder(nMaxOrder, CkCallbackResumeThread());
	  }

      writeOutput(0);
      if(param.bDoGas) {
	  ckout << "Outputting gas properties ...";
	  if(printBinaryAcc)
	      CkAssert(0);
	  else {
	      GasDenOutputParams pDenOut(achFile, param.iBinaryOut, 0.0);
	      PresOutputParams pPresOut(achFile, param.iBinaryOut, 0.0);
	      HsmOutputParams pSphHOut(achFile, param.iBinaryOut, 0.0);
	      DivVOutputParams pDivVOut(achFile, param.iBinaryOut, 0.0);
	      PDVOutputParams pPDVOut(achFile, param.iBinaryOut, 0.0);
	      MuMaxOutputParams pMuMaxOut(achFile, param.iBinaryOut, 0.0);
	      BSwOutputParams pBSwOut(achFile, param.iBinaryOut, 0.0);
	      CsOutputParams pCsOut(achFile, param.iBinaryOut, 0.0);
              if (param.iBinaryOut) {
                  treeProxy[0].outputBinary(pPresOut, param.bParaWrite,
                      CkCallbackResumeThread());
                  treeProxy[0].outputBinary(pDivVOut, param.bParaWrite,
                      CkCallbackResumeThread());
                  treeProxy[0].outputBinary(pPDVOut, param.bParaWrite,
                      CkCallbackResumeThread());
                  treeProxy[0].outputBinary(pMuMaxOut, param.bParaWrite,
                      CkCallbackResumeThread());
                  treeProxy[0].outputBinary(pBSwOut, param.bParaWrite,
                      CkCallbackResumeThread());
                  treeProxy[0].outputBinary(pCsOut, param.bParaWrite,
                      CkCallbackResumeThread());
#ifndef COOLING_NONE
                  if(param.bGasCooling) {
                      EDotOutputParams pEDotOut(achFile, param.iBinaryOut, 0.0);
                      treeProxy[0].outputBinary(pEDotOut, param.bParaWrite,
                          CkCallbackResumeThread());
                      }
#endif
                  }
              else {
                  treeProxy[0].outputASCII(pDenOut, param.bParaWrite, CkCallbackResumeThread());
                  treeProxy[0].outputASCII(pPresOut, param.bParaWrite, CkCallbackResumeThread());
                  treeProxy[0].outputASCII(pSphHOut, param.bParaWrite, CkCallbackResumeThread());
                  treeProxy[0].outputASCII(pDivVOut, param.bParaWrite, CkCallbackResumeThread());
                  treeProxy[0].outputASCII(pPDVOut, param.bParaWrite, CkCallbackResumeThread());
                  treeProxy[0].outputASCII(pMuMaxOut, param.bParaWrite, CkCallbackResumeThread());
                  treeProxy[0].outputASCII(pBSwOut, param.bParaWrite, CkCallbackResumeThread());
                  treeProxy[0].outputASCII(pCsOut, param.bParaWrite, CkCallbackResumeThread());
#ifndef COOLING_NONE
                  if(param.bGasCooling) {
                      EDotOutputParams pEDotOut(achFile, 0, 0.0);
                      treeProxy[0].outputASCII(pEDotOut, param.bParaWrite,
					   CkCallbackResumeThread());
                      }
#endif
                  }
	      }
	  }
      ckout << "Outputting accelerations  ...";
      if(printBinaryAcc)
	  treeProxy[0].outputAccelerations(OrientedBox<double>(),
					   "acc2", CkCallbackResumeThread());
      else {
	  AccOutputParams pAcc(achFile, param.iBinaryOut, 0.0);
          if(param.iBinaryOut)
              treeProxy[0].outputBinary(pAcc, param.bParaWrite, CkCallbackResumeThread());
          else
              treeProxy[0].outputASCII(pAcc, param.bParaWrite, CkCallbackResumeThread());
	  }
#ifdef NEED_DT
      ckout << "Outputting dt ...";
      adjust(0);
      DtOutputParams pDt(achFile, param.iBinaryOut, 0.0);
      if(param.iBinaryOut)
          treeProxy[0].outputBinary(pDt, param.bParaWrite, CkCallbackResumeThread());
      else
          treeProxy[0].outputASCII(pDt, param.bParaWrite, CkCallbackResumeThread());
#endif
      RungOutputParams pRung(achFile, param.iBinaryOut, 0.0);
      KeyOutputParams pKey(achFile, param.iBinaryOut, 0.0);
      DomainOutputParams pDomain(achFile, param.iBinaryOut, 0.0);
      if(param.iBinaryOut) {
          treeProxy[0].outputIntBinary(pRung, param.bParaWrite, CkCallbackResumeThread());
          treeProxy[0].outputBinary(pKey, param.bParaWrite, CkCallbackResumeThread());
          treeProxy[0].outputBinary(pDomain, param.bParaWrite, CkCallbackResumeThread());
          }
      else {
          treeProxy[0].outputIntASCII(pRung, param.bParaWrite, CkCallbackResumeThread());
          treeProxy[0].outputASCII(pKey, param.bParaWrite, CkCallbackResumeThread());
          treeProxy[0].outputASCII(pDomain, param.bParaWrite, CkCallbackResumeThread());
          }
      if(param.bDoGas && param.bDoDensity) {
	  // The following call is to get the particles in key order
	  // before the sort.
	  treeProxy.drift(0.0, 0, 0, 0.0, 0.0, 0, true, CkCallbackResumeThread());
	  sorter.startSorting(dataManagerID, ddTolerance,
			      CkCallbackResumeThread(), true);
#ifdef PUSH_GRAVITY
	  treeProxy.buildTree(bucketSize, CkCallbackResumeThread(),true);
#else
	  treeProxy.buildTree(bucketSize, CkCallbackResumeThread());
#endif
	  
	  ckout << "Calculating total densities ...";
	  DensitySmoothParams pDen(TYPE_GAS|TYPE_DARK|TYPE_STAR, 0);
	  startTime = CkWallTimer();
	  double dfBall2OverSoft2 = 0.0;
	  treeProxy.startSmooth(&pDen, 1, param.nSmooth, dfBall2OverSoft2,
					 CkCallbackResumeThread());
	  ckout << " took " << (CkWallTimer() - startTime) << " seconds."
		<< endl;
          ckout << "Recalculating densities ...";
          startTime = CkWallTimer();
	  treeProxy.startReSmooth(&pDen, CkCallbackResumeThread());
          ckout << " took " << (CkWallTimer() - startTime) << " seconds." << endl;
	  treeProxy.finishNodeCache(CkCallbackResumeThread());
          ckout << "Reodering ...";
          startTime = CkWallTimer();
	  treeProxy.reOrder(nMaxOrder, CkCallbackResumeThread());
          ckout << " took " << (CkWallTimer() - startTime) << " seconds." << endl;
	  ckout << "Outputting densities ...";
	  startTime = CkWallTimer();
	  DenOutputParams pDenOut(achFile, param.iBinaryOut, 0.0);
          if(param.iBinaryOut)
              treeProxy[0].outputBinary(pDenOut, param.bParaWrite, CkCallbackResumeThread());
          else
              treeProxy[0].outputASCII(pDenOut, param.bParaWrite, CkCallbackResumeThread());
	  ckout << " took " << (CkWallTimer() - startTime) << " seconds."
		<< endl;
	  ckout << "Outputting hsmooth ...";
	  HsmOutputParams pHsmOut(achFile, param.iBinaryOut, 0.0);
          if(param.iBinaryOut)
              treeProxy[0].outputBinary(pHsmOut, param.bParaWrite, CkCallbackResumeThread());
          else
              treeProxy[0].outputASCII(pHsmOut, param.bParaWrite, CkCallbackResumeThread());
	  }
  }
	
#if COSMO_STATS > 0
  ckerr << "Outputting statistics ...";
  startTime = CkWallTimer();
  //Interval<unsigned int> dummy;
	
  treeProxy[0].outputStatistics(CkCallbackResumeThread());

  ckerr << " took " << (CkWallTimer() - startTime) << " seconds." << endl;
#endif

  ckout << "Done." << endl;
	
#ifdef HPM_COUNTER
  cacheNode.stopHPM(CkCallbackResumeThread());
  cacheGravPart.stopHPM(CkCallbackResumeThread());
  cacheSmoothPart.stopHPM(CkCallbackResumeThread());
#endif
  ckout << endl << "******************" << endl << endl; 
  // Some memory cleanup
  delete param.stfm;
  treeProxy.ckDestroy();
  CkWaitQD();
  CkExit();
}

///
/// \brief Calculate various energy and momentum quantities, and output them
/// to a log file.
///

void
Main::calcEnergy(double dTime, double wallTime, const char *achLogFileName)
{
    CkReductionMsg *msg;
    treeProxy.calcEnergy(CkCallbackResumeThread((void*&)msg));
    double *dEnergy = (double *) msg->getData();
    static int first = 1;

    FILE *fpLog = fopen(achLogFileName, "a");
    
    double a = csmTime2Exp(param.csm, dTime);
    
    if(first && (!bIsRestarting || dTimeOld == 0.0)) {
	fprintf(fpLog, "# time redshift TotalEVir TotalE Kinetic Virial Potential TotalECosmo Ethermal Lx Ly Lz Wallclock\n");
	dEcosmo = 0.0;
	first = 0;
	}
    else {
	/*
	 * Estimate integral (\\dot a*U*dt) over the interval.
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

#include <errno.h>
/// @brief mkdir with error checking
inline int safeMkdir(const char *achFile) {
    int ret = mkdir(achFile, 0755);
    if(ret != 0 && errno == EEXIST) {
        CkError("WARNING: overwriting existing directory or file\n");
        ret = 0;
        }
    return ret;
    }

///
/// @brief Output a snapshot
/// @param iStep Timestep we are outputting, used for file name.
///

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
    
    treeProxy.reOrder(nMaxOrder, CkCallbackResumeThread());
    if(verbosity)
	ckout << " took " << (CkWallTimer() - startTime) << " seconds."
	      << endl;

    double duTFac = (param.dConstGamma-1)*param.dMeanMolWeight/param.dGasConst;
    
    if(param.iBinaryOut != 6)
        {
        if(verbosity) {
            ckout << "Writing Tipsy file ...";
            startTime = CkWallTimer();
            }
        if(param.bParaWrite)
            treeProxy.setupWrite(0, 0, achFile, dOutTime, dvFac, duTFac,
                                 param.bDoublePos, param.bDoubleVel,
                                 param.bGasCooling, CkCallbackResumeThread());
        else
            treeProxy[0].serialWrite(0, achFile, dOutTime, dvFac, duTFac,
                                     param.bDoublePos, param.bDoubleVel,
                                     param.bGasCooling, CkCallbackResumeThread());
        if(verbosity)
            ckout << " took " << (CkWallTimer() - startTime) << " seconds."
                  << endl;
        }
    else { // N-Chilada output
        // Set up N-Chilada directory structure
        CkAssert(safeMkdir(achFile) == 0);
        if(nTotalSPH > 0) {
            string dirname(string(achFile) + "/gas");
            CkAssert(safeMkdir(dirname.c_str()) == 0);
            }
        if(nTotalDark > 0) {
            string dirname(string(achFile) + "/dark");
            CkAssert(safeMkdir(dirname.c_str()) == 0);
            }
        if(nTotalStar > 0) {
            string dirname(string(achFile) + "/star");
            CkAssert(safeMkdir(dirname.c_str()) == 0);
            }
        MassOutputParams pMassOut(achFile, param.iBinaryOut, dOutTime);
        treeProxy[0].outputBinary(pMassOut, param.bParaWrite,
                                  CkCallbackResumeThread());
        PosOutputParams pPosOut(achFile, param.iBinaryOut, dOutTime);
        treeProxy[0].outputBinary(pPosOut, param.bParaWrite,
                                  CkCallbackResumeThread());
        VelOutputParams pVelOut(achFile, param.iBinaryOut, dOutTime, dvFac);
        treeProxy[0].outputBinary(pVelOut, param.bParaWrite,
                                  CkCallbackResumeThread());
        SoftOutputParams pSoftOut(achFile, param.iBinaryOut, dOutTime);
        treeProxy[0].outputBinary(pSoftOut, param.bParaWrite,
                                  CkCallbackResumeThread());
        PotOutputParams pPotOut(achFile, param.iBinaryOut, dOutTime);
        treeProxy[0].outputBinary(pPotOut, param.bParaWrite,
                                  CkCallbackResumeThread());
        if(nTotalSPH > 0) {
            GasDenOutputParams pGasDenOut(achFile, param.iBinaryOut, dOutTime);
            treeProxy[0].outputBinary(pGasDenOut, param.bParaWrite,
                                      CkCallbackResumeThread());
            TempOutputParams pTempOut(achFile, param.iBinaryOut, dOutTime,
                param.bGasCooling, duTFac);
            treeProxy[0].outputBinary(pTempOut, param.bParaWrite,
                                      CkCallbackResumeThread());
            HsmOutputParams pHsmOut(achFile, param.iBinaryOut, dOutTime);
            treeProxy[0].outputBinary(pHsmOut, param.bParaWrite,
                                      CkCallbackResumeThread());
            }
        if(nTotalSPH > 0 || nTotalStar > 0) {
            MetalsOutputParams pMetalsOut(achFile, param.iBinaryOut, dOutTime);
            treeProxy[0].outputBinary(pMetalsOut, param.bParaWrite,
                                      CkCallbackResumeThread());
            }
        if(nTotalStar > 0) {
            TimeFormOutputParams pTimeFormOut(achFile, param.iBinaryOut, dOutTime);
            treeProxy[0].outputBinary(pTimeFormOut, param.bParaWrite,
                                      CkCallbackResumeThread());
            }
        }
    
    if(verbosity) {
	ckout << "Writing arrays ...";
	startTime = CkWallTimer();
	}

    OxOutputParams pOxOut(achFile, param.iBinaryOut, dOutTime);
    FeOutputParams pFeOut(achFile, param.iBinaryOut, dOutTime);
    MFormOutputParams pMFormOut(achFile, param.iBinaryOut, dOutTime);
    coolontimeOutputParams pcoolontimeOut(achFile, param.iBinaryOut, dOutTime);
    ESNRateOutputParams pESNRateOut(achFile, param.iBinaryOut, dOutTime);
#ifndef COOLING_NONE
    Cool0OutputParams pCool0Out(achFile, param.iBinaryOut, dOutTime);
    Cool1OutputParams pCool1Out(achFile, param.iBinaryOut, dOutTime);
    Cool2OutputParams pCool2Out(achFile, param.iBinaryOut, dOutTime);
#endif
    SoftOutputParams pSoftOut(achFile, param.iBinaryOut, dOutTime);
    HsmOutputParams pHsmOut(achFile, param.iBinaryOut, dOutTime);
    CsOutputParams pCSOut(achFile, param.iBinaryOut, dOutTime);

    if (param.iBinaryOut) {
	if (param.bStarForm || param.bFeedback) {
	    treeProxy[0].outputBinary(pOxOut, param.bParaWrite,
				      CkCallbackResumeThread());
	    treeProxy[0].outputBinary(pFeOut, param.bParaWrite,
				      CkCallbackResumeThread());
	    treeProxy[0].outputBinary(pMFormOut, param.bParaWrite,
				      CkCallbackResumeThread());
	    treeProxy[0].outputBinary(pcoolontimeOut, param.bParaWrite,
				      CkCallbackResumeThread());
	    treeProxy[0].outputBinary(pESNRateOut, param.bParaWrite,
				      CkCallbackResumeThread());
	    }
#ifndef COOLING_NONE
	if(param.bGasCooling) {
	    treeProxy[0].outputBinary(pCool0Out, param.bParaWrite,
				     CkCallbackResumeThread());
	    treeProxy[0].outputBinary(pCool1Out, param.bParaWrite,
				     CkCallbackResumeThread());
	    treeProxy[0].outputBinary(pCool2Out, param.bParaWrite,
				      CkCallbackResumeThread());
	    }
#endif
	if(param.bDoSoftOutput && param.iBinaryOut != 6) {
	    treeProxy[0].outputBinary(pSoftOut, param.bParaWrite,
				      CkCallbackResumeThread());
            }
	    
	if(param.bDohOutput && param.iBinaryOut != 6) {
	    treeProxy[0].outputBinary(pHsmOut, param.bParaWrite,
				      CkCallbackResumeThread());
            }
	if(param.bDoCSound)
	    treeProxy[0].outputBinary(pCSOut, param.bParaWrite,
				      CkCallbackResumeThread());
	if(param.bDoIOrderOutput || param.bStarForm || param.bFeedback) {
	    IOrderOutputParams pIOrdOut(achFile, param.iBinaryOut, dOutTime);
	    treeProxy[0].outputIntBinary(pIOrdOut, param.bParaWrite,
                                         CkCallbackResumeThread());
	    if(param.bStarForm) {
		IGasOrderOutputParams pIGasOrdOut(achFile, param.iBinaryOut,
                    dOutTime);
		treeProxy[0].outputIntBinary(pIGasOrdOut, param.bParaWrite,
                                             CkCallbackResumeThread());
                }
            }
	} else {
	if (param.bStarForm || param.bFeedback) {
	    treeProxy[0].outputASCII(pOxOut, param.bParaWrite,
				     CkCallbackResumeThread());
	    treeProxy[0].outputASCII(pFeOut, param.bParaWrite,
				     CkCallbackResumeThread());
	    treeProxy[0].outputASCII(pMFormOut, param.bParaWrite,
				      CkCallbackResumeThread());
	    treeProxy[0].outputASCII(pcoolontimeOut, param.bParaWrite,
				     CkCallbackResumeThread());
	    treeProxy[0].outputASCII(pESNRateOut, param.bParaWrite,
				      CkCallbackResumeThread());
	    }
#ifndef COOLING_NONE
	if(param.bGasCooling) {
	    treeProxy[0].outputASCII(pCool0Out, param.bParaWrite,
				     CkCallbackResumeThread());
	    treeProxy[0].outputASCII(pCool1Out, param.bParaWrite,
				     CkCallbackResumeThread());
	    treeProxy[0].outputASCII(pCool2Out, param.bParaWrite,
				     CkCallbackResumeThread());
	    }
#endif
	if(param.bDoSoftOutput)
	    treeProxy[0].outputASCII(pSoftOut, param.bParaWrite,
				      CkCallbackResumeThread());
	if(param.bDohOutput)
	    treeProxy[0].outputASCII(pHsmOut, param.bParaWrite,
				     CkCallbackResumeThread());
	if(param.bDoCSound)
	    treeProxy[0].outputASCII(pCSOut, param.bParaWrite,
				      CkCallbackResumeThread());
	if(param.bDoIOrderOutput || param.bStarForm || param.bFeedback) {
	    IOrderOutputParams pIOrdOut(achFile, param.iBinaryOut, dOutTime);
	    treeProxy[0].outputIntASCII(pIOrdOut, param.bParaWrite,
					CkCallbackResumeThread());
	    if(param.bStarForm) {
		IGasOrderOutputParams pIGasOrdOut(achFile, param.iBinaryOut,
                    dOutTime);
		treeProxy[0].outputIntASCII(pIGasOrdOut, param.bParaWrite,
					    CkCallbackResumeThread());
	      }
	  }
	}
      
    if(verbosity)
	ckout << " took " << (CkWallTimer() - startTime) << " seconds."
	      << endl;
    if(param.nSteps != 0 && param.bDoDensity) {
	// The following call is to get the particles in key order
	// before the sort.
	treeProxy.drift(0.0, 0, 0, 0.0, 0.0, 0, true, CkCallbackResumeThread());
	sorter.startSorting(dataManagerID, ddTolerance,
			    CkCallbackResumeThread(), true);
#ifdef PUSH_GRAVITY
	treeProxy.buildTree(bucketSize, CkCallbackResumeThread(),true);
#else
	treeProxy.buildTree(bucketSize, CkCallbackResumeThread());
#endif

	if(verbosity)
	    ckout << "Calculating total densities ...";
	DensitySmoothParams pDen(TYPE_GAS|TYPE_DARK|TYPE_STAR, 0);
	startTime = CkWallTimer();
	double dfBall2OverSoft2 = 0.0;
	treeProxy.startSmooth(&pDen, 1, param.nSmooth, dfBall2OverSoft2,
			      CkCallbackResumeThread());
	treeProxy.finishNodeCache(CkCallbackResumeThread());
	if(verbosity) {
	    ckout << " took " << (CkWallTimer() - startTime) << " seconds."
		  << endl;
	    ckout << "Reodering ...";
	    }
	startTime = CkWallTimer();
	treeProxy.reOrder(nMaxOrder, CkCallbackResumeThread());
	if(verbosity) {
	    ckout << " took " << (CkWallTimer() - startTime) << " seconds."
		  << endl;
	    ckout << "Outputting densities ...";
	    }
	startTime = CkWallTimer();
	DenOutputParams pDenOut(string(achFile), param.iBinaryOut, dOutTime);
	if (param.iBinaryOut)
	    treeProxy[0].outputBinary(pDenOut, param.bParaWrite, CkCallbackResumeThread());
	else
	    treeProxy[0].outputASCII(pDenOut, param.bParaWrite, CkCallbackResumeThread());
	ckout << " took " << (CkWallTimer() - startTime) << " seconds."
	      << endl;

	if(param.bDoGas) {	// recalculate gas densities for timestep
	    if(verbosity)
		ckout << "Recalculating gas densities ...";
	    startTime = CkWallTimer();
	    // The following call is to get the particles in key order
	    // before the sort.
	    treeProxy.drift(0.0, 0, 0, 0.0, 0.0, 0, true, CkCallbackResumeThread());
	    sorter.startSorting(dataManagerID, ddTolerance,
				CkCallbackResumeThread(), true);
#ifdef PUSH_GRAVITY
	    treeProxy.buildTree(bucketSize, CkCallbackResumeThread(),true);
#else
	    treeProxy.buildTree(bucketSize, CkCallbackResumeThread());
#endif
	    DensitySmoothParams pDenGas(TYPE_GAS, 0);
	    dfBall2OverSoft2 = 4.0*param.dhMinOverSoft*param.dhMinOverSoft;
	    treeProxy.startSmooth(&pDenGas, 1, param.nSmooth, dfBall2OverSoft2,
				  CkCallbackResumeThread());
	    treeProxy.finishNodeCache(CkCallbackResumeThread());
	    if(verbosity)
		ckout << " took " << (CkWallTimer() - startTime) << " seconds."
		      << endl;
	    }
	}
    }

///
/// @brief Calculate timesteps of particles.
///
/// Particles on the KickRung and shorter have their timesteps
/// adjusted.  Calls TreePiece::adjust().
///
/// @param iKickRung Rung (and above) about to be kicked.
///

int Main::adjust(int iKickRung) 
{
    CkReductionMsg *msg;
    double a = csmTime2Exp(param.csm,dTime);
    
    treeProxy.adjust(iKickRung, param.bEpsAccStep, param.bGravStep,
		     param.bSphStep, param.bViscosityLimitdt,
		     param.dEta, param.dEtaCourant,
		     param.dEtauDot, param.dDelta, 1.0/(a*a*a), a,
		     0.0,  /* set to dhMinOverSoft if we implement
			      Gasoline's LowerSoundSpeed. */
		     param.bDoGas,
		     CkCallbackResumeThread((void*&)msg));

    int iCurrMaxRung = ((int *)msg->getData())[0];
    int nMaxRung = ((int *)msg->getData())[1];
    delete msg;
    if(nMaxRung <= param.nTruncateRung && iCurrMaxRung > iKickRung) {
	if(verbosity)
	    CkPrintf("n_CurrMaxRung = %d, iCurrMaxRung = %d: promoting particles\n", nMaxRung, iCurrMaxRung);
	iCurrMaxRung--;
	treeProxy.truncateRung(iCurrMaxRung, CkCallbackResumeThread());
	}
    return iCurrMaxRung;
    }

///
/// @brief Count and print out the number of particles in each
/// timestep bin.
///
/// This routine is for information only.
///
void Main::rungStats() 
{
    CkReductionMsg *msg;
    treeProxy.rungStats(CkCallbackResumeThread((void*&)msg));
    int64_t *nInRung = (int64_t *)msg->getData();
    
    ckout << "Rung distribution: (";
    for(int iRung = 0; iRung <= MAXRUNG; iRung++)
	ckout << " " << nInRung[iRung] << ",";
    ckout << ")" << endl;
    
    delete msg;
    }

///
/// @brief determine number of active particles in the given rung
///
/// Calls TreePiece::countActive() and sets Main::nActiveGrav and
/// Main::nActiveSPH.
///
void Main::countActive(int activeRung) 
{
    CkReductionMsg *msg;
    treeProxy.countActive(activeRung, CkCallbackResumeThread((void*&)msg));
    int64_t *nActive = (int64_t *)msg->getData();
    
    nActiveGrav = nActive[0];
    nActiveSPH = nActive[1];
    
    ckout << "Gravity Active: " << nActive[0]
	  << ", Gas Active: " << nActive[1] << endl ;
    
    delete msg;
    }

/**
 * \brief Change the softening in comoving coordinates
 *
 * When compiled with -DCHANGESOFT, and bPhysicalSoft is set, the
 * (comoving) softening is changed so that it is constant in physical units.
 */
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

///
/// @brief Slowly increase mass of a subset of particles.
///
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
		  FILE *fp = fopen(achFile, "r" );
		  if(fp == NULL)
		      CkAbort("DumpFrame: photogenic specified, but no photogenic file\n");
		  fclose(fp);

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

#include "DumpFrameData.h"

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
	dfDataProxy.clearFrame(in, CkCallbackResumeThread());
	treeProxy.DumpFrame(in, CkCallbackResumeThread(), false);
	dfDataProxy.combineFrame(in, CkCallbackResumeThread((void*&)msgDF));
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

/**
 * \brief Coalesce all added and deleted particles and update global
 * quantities.
 */

void Main::addDelParticles()
{
    CkReductionMsg *msg;
    
    if(verbosity)
	CkPrintf("Changing Particle number\n");
	
    treeProxy.colNParts(CkCallbackResumeThread((void*&)msg));
    // an array of one element per treepiece produced by CkReduction::concat
    // @TODO: replace this with a proper scan
    CountSetPart *counts = (CountSetPart *) msg->getData();
    CkAssert(msg->getSize() == numTreePieces*sizeof(*counts));

    // Callback for neworder
    CkCallbackResumeThread cb;

    int iPiece = 0;
    for(iPiece = 0; iPiece < numTreePieces; iPiece++) {
	treeProxy[counts[iPiece].index].newOrder(nMaxOrderGas+1,
						 nMaxOrderDark+1,
						 nMaxOrder+1, cb);

	nMaxOrderGas += counts[iPiece].nAddGas;
	nMaxOrderDark += counts[iPiece].nAddDark;
	nMaxOrder += counts[iPiece].nAddStar;
	nTotalSPH += counts[iPiece].nAddGas - counts[iPiece].nDelGas;
	nTotalDark += counts[iPiece].nAddDark - counts[iPiece].nDelDark;
	nTotalStar += counts[iPiece].nAddStar - counts[iPiece].nDelStar;
	}
    nTotalParticles = nTotalSPH + nTotalDark + nTotalStar;
    delete msg;
    if (verbosity)
	CkPrintf("New numbers of particles: %d gas %d dark %d star\n",
		 nTotalSPH, nTotalDark, nTotalStar);
    
    cb.thread_delay();
    treeProxy.setNParts(nTotalSPH, nTotalDark, nTotalStar,
			CkCallbackResumeThread());
    }

/**
 * Diagnostic function to summmarize memory usage across all
 * processors
 */
void Main::memoryStats() 
{
    CkReductionMsg *msg;
    dMProxy.memoryStats(CkCallbackResumeThread((void*&)msg));
    // simple maximum for now.
    int *memMax = (int *)msg->getData();
    CkPrintf("Maximum memory: %d MB\n", *memMax);
    delete msg;
    }

/**
 * Diagnostic function to summmarize cache memory usage across all
 * processors
 */
void Main::memoryStatsCache() 
{
    CkReductionMsg *msg;
    treeProxy.memCacheStats(CkCallbackResumeThread((void*&)msg));
    // simple maximum for now.
    int *memMax = (int *)msg->getData();
    CkPrintf("Maximum memory with cache: %d MB, after: %d MB\n", memMax[0],
	     memMax[1]);
    CkPrintf("Maximum cache entries: node: %d , particle: %d\n", memMax[2],
	     memMax[3]);
    delete msg;
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
    p | nTotalDark;
    p | nTotalStar;
    p | nMaxOrderGas;
    p | nMaxOrderDark;
    p | nMaxOrder;
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

#ifdef SELECTIVE_TRACING
void Main::turnProjectionsOn(int activeRung){
  CkAssert(!projectionsOn);
  if(numMaxTrace <= 0) return;
  if(traceState == TraceNormal){
    if(activeRung != monitorRung){
    }
    else if(traceIteration < numTraceIterations){
      if(monitorStart == 0){
        prjgrp.on(CkCallbackResumeThread());
        projectionsOn = true;
        traceIteration++;
        numMaxTrace--;
      }
      else{
        monitorStart--;
      }
    }
    else{
      traceState = TraceSkip;
      traceIteration = 1;
    }
  }
  else if(traceState == TraceSkip){
    if(activeRung != monitorRung){
    }
    else if(traceIteration < numSkipIterations){
      traceIteration++;
    }
    else{
      traceState = TraceNormal;
      if(monitorStart == 0){
        prjgrp.on(CkCallbackResumeThread());
        projectionsOn = true;
        traceIteration = 1;
        numMaxTrace--;
      }
      else{
        monitorStart--;
      }
    }
  }
}


void Main::turnProjectionsOff(){
  if(projectionsOn){
    prjgrp.off(CkCallbackResumeThread());
    projectionsOn = false;
  }
}
#endif

const char *typeString(NodeType type);

void GenericTreeNode::getGraphViz(std::ostream &out){
#ifdef BIGKEYS
  uint64_t lower = getKey();
  uint64_t upper = getKey() >> 64;
  out << upper << lower
    << "[label=\"" << upper << lower
#else    
  out << getKey()
    << "[label=\"" << getKey()
#endif
    << " " 
    << typeString(getType())
    << "\\n"
    << moments.totalMass << " "
    << "(" << moments.cm.x << ","
    << moments.cm.y << ","
    << moments.cm.z << ")"
    << "\"];"
    << std::endl;

  if(getType() == Boundary){
    for(int i = 0; i < numChildren(); i++){
      GenericTreeNode *child = getChildren(i);
#ifdef BIGKEYS
      uint64_t ch_lower = getKey();
      uint64_t ch_upper = getKey() >> 64;
      out << upper << lower << "->" << ch_upper << ch_lower << ";" << std::endl;
#else
      out << getKey() << "->" << child->getKey() << ";" << std::endl;
#endif
    }
  }
}

void printTreeGraphVizRecursive(GenericTreeNode *node, ostream &out){
  node->getGraphViz(out);
  out << endl;
  if(node->getType() == Boundary){
    for(int i = 0; i < node->numChildren(); i++){
      printTreeGraphVizRecursive(node->getChildren(i),out);
    }
  }
}

void printTreeGraphViz(GenericTreeNode *node, ostream &out, const string &name){
  out << "digraph " << name << " {" << endl;
  printTreeGraphVizRecursive(node,out);
  out << "}" << endl;
}


#include "ParallelGravity.def.h"

/*
#define CK_TEMPLATES_ONLY
#include "CacheManager.def.h"
#undef CK_TEMPLATES_ONLY
#include "CacheManager.def.h"
*/
