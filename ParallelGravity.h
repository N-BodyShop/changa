/** @file ParallelGravity.h
 */

#ifndef PARALLELGRAVITY_H
#define PARALLELGRAVITY_H

#include "config.h"
#include <string>
#include <map>
#include <vector>
#include <algorithm>

#include "pup_stl.h"
#include "ckio.h"

#include "Vector3D.h"
#include "tree_xdr.h"
#include "TipsyFile.h"
#include "SFC.h"
#include "TreeNode.h"
#include "GenericTreeNode.h"
#include "Interval.h"
#include "parameters.h"
#include "param.h"
#include "dumpframe.h"
#include <liveViz.h>

#include "TaggedVector3D.h"

#include "codes.h"
#include "CacheInterface.h"

#ifdef SPCUDA
#include "EwaldCUDA.h"
#endif

#ifdef CUDA
#include "cuda_typedef.h"
#endif

#include "keytype.h"

PUPbytes(InDumpFrame);
PUPbytes(COOL);
PUPbytes(COOLPARAM);

#ifdef HPM_COUNTER
#include <libhpm.h>
#endif

#include <map>

#define MERGE_REMOTE_REQUESTS_VERBOSE /*CkPrintf*/

using namespace std;

using namespace Tree;

/// Load balancers that need the spatial information.
enum LBStrategy{
  Null=0,
  Multistep,
  Orb3d,
  Multistep_notopo,
  MultistepNode_notopo,
  Orb3d_notopo,
  MultistepOrb,
  HierarchOrb
};
PUPbytes(LBStrategy);

#ifdef SELECTIVE_TRACING
enum TraceState {
  TraceNormal = 0,
  TraceSkip
};
#endif

/// Possible domain decomposition methods
enum DomainsDec {
    SFC_dec=0,	// Space Filling Curve with Morton ordering
    Oct_dec=1, 	// Oct tree
    ORB_dec=2,	// Bisect the longest axis, balancing particles
    SFC_peano_dec=3,	// SFC with Peano-Hilbert ordering
    SFC_peano_dec_3D=4, // Joachim Stadel's implementation of P-H ordering
    SFC_peano_dec_2D=5,	// 2D version of Peano-Hilbert ordering
    ORB_space_dec=6		// Bisect space
};

/// Directions for sending boundaries
enum NborDir {
  LEFT = 0,
  RIGHT
};
PUPbytes(NborDir);

/// tolerance for unequal pieces in SFC based decompositions.
const double ddTolerance = 0.1;

inline void operator|(PUP::er &p,DomainsDec &d) {
  int di;
  if (p.isUnpacking()) {
    p | di;
    d = (DomainsDec)di;
  } else {
    di = (int)d;
    p | di;
  }
}

#include "GravityParticle.h"

class SmoothParams;

///  Class for new maxOrder broadcast
class NewMaxOrder
{
 public:
    int64_t nMaxOrderGas;
    int64_t nMaxOrderDark;
    int64_t nMaxOrder;
    void pup(PUP::er& p) {
        p| nMaxOrderGas;
        p| nMaxOrderDark;
        p| nMaxOrder;
        }
    };

#include "InOutput.h"

#include "ParallelGravity.decl.h"

extern CProxy_Main mainChare;
extern int verbosity;
extern bool _cache;
extern int _nocache;
extern int _cacheLineDepth;
extern unsigned int _yieldPeriod;
extern DomainsDec domainDecomposition;
extern double dExtraStore;
extern double dMaxBalance;
extern double dFracLoadBalance;
extern double dGlassDamper;
extern int bUseCkLoopPar;
extern GenericTrees useTree;
extern CProxy_TreePiece treeProxy;
#ifdef REDUCTION_HELPER
extern CProxy_ReductionHelper reductionHelperProxy;
#endif
extern CProxy_LvArray lvProxy;	    // Proxy for the liveViz array
extern CProxy_LvArray smoothProxy;  // Proxy for smooth reduction
extern CProxy_LvArray gravityProxy; // Proxy for gravity reduction
extern CProxy_TreePiece streamingProxy;
extern CProxy_DataManager dMProxy;
extern CProxy_IntraNodeLBManager nodeLBMgrProxy;
extern unsigned int numTreePieces;
extern unsigned int particlesPerChare;
extern int nIOProcessor;

extern CProxy_DumpFrameData dfDataProxy;
extern CProxy_PETreeMerger peTreeMergerProxy;
extern CProxy_CkCacheManager<KeyType> cacheGravPart;
extern CProxy_CkCacheManager<KeyType> cacheSmoothPart;
extern CProxy_CkCacheManager<KeyType> cacheNode;

/// The group ID of your DataManager.  You must set this!
extern CkGroupID dataManagerID;

extern int boundaryEvaluationUE;
extern int weightBalanceUE;
extern int networkProgressUE;
extern int nodeForceUE;
extern int partForceUE;

extern int tbFlushRequestsUE;
extern int prefetchDoneUE;

extern int _prefetch;
extern int _randChunks;
extern int _numChunks;
extern unsigned int bucketSize;

//jetley
extern int localNodesPerReq;
extern int remoteNodesPerReq;
extern int remoteResumeNodesPerReq;
extern int localPartsPerReq;
extern int remotePartsPerReq;
extern int remoteResumePartsPerReq;

extern double largePhaseThreshold;

extern cosmoType theta;
extern cosmoType thetaMono;

extern int numInitDecompBins;
extern int octRefineLevel;

/// @brief Message to efficiently start entry methods with no arguments.
class dummyMsg : public CMessage_dummyMsg{
public:
};

#if COSMO_STATS > 0
class TreePieceStatistics {
  u_int64_t nodesOpenedLocal;
  u_int64_t nodesOpenedRemote;
  u_int64_t nodeInterLocal;
  u_int64_t nodeInterRemote;
  u_int64_t particleInterLocal;
  u_int64_t particleInterRemote;
  u_int64_t openCriterionCalls;
  int nActive;

  TreePieceStatistics() : nodesOpenedLocal(0), nodesOpenedRemote(0),
      nodeInterLocal(0), nodeInterRemote(0), particleInterLocal(0),
      particleInterRemote(0), openCriterionCalls(0), nActive(0) { }

 public:
  TreePieceStatistics(u_int64_t nol, u_int64_t nor, u_int64_t occ, u_int64_t nil, u_int64_t nir,
		      u_int64_t pil, u_int64_t pir, int na) :
      nodesOpenedLocal(nol), nodesOpenedRemote(nor), nodeInterLocal(nil),
      nodeInterRemote(nir), particleInterLocal(pil), particleInterRemote(pir),
      openCriterionCalls(occ), nActive(na) { }

  void printTo(CkOStream &os) {
    os << "  TreePiece: " << nActive << " particles active." << endl;
    os << "  TreePiece: " << nodesOpenedLocal << " local nodes opened, ";
    os << nodesOpenedRemote << " remote" << endl;
    os << "  TreePiece: " << openCriterionCalls << " num of open criterion calls" << endl;
    os << "  TreePiece: " << nodeInterLocal << " local particle-node interactions, ";
    os << nodeInterRemote << " remote" << endl;
    os << "  TreePiece: " << particleInterLocal << " local particle-particle interactions, ";
    os << particleInterRemote << " remote" << endl;
    os << "  TreePiece: "
       << (particleInterLocal + particleInterRemote)/(double) nActive
       << " particles, "
       << (nodeInterLocal + nodeInterRemote)/(double) nActive
       << " nodes per particle" << endl;
  }

  static CkReduction::reducerType sum;

  static CkReductionMsg *sumFn(int nMsg, CkReductionMsg **msgs) {
    TreePieceStatistics ret;
    for (int i=0; i<nMsg; ++i) {
      CkAssert(msgs[i]->getSize() == sizeof(TreePieceStatistics));
      TreePieceStatistics *data = (TreePieceStatistics *)msgs[i]->getData();
      ret.nodesOpenedLocal += data->nodesOpenedLocal;
      ret.nodesOpenedRemote += data->nodesOpenedRemote;
      ret.openCriterionCalls += data->openCriterionCalls;
      ret.nodeInterLocal += data->nodeInterLocal;
      ret.nodeInterRemote += data->nodeInterRemote;
      ret.particleInterLocal += data->particleInterLocal;
      ret.particleInterRemote += data->particleInterRemote;
      ret.nActive += data->nActive;
    }
    return CkReductionMsg::buildNew(sizeof(TreePieceStatistics), &ret);
  }
};
#endif

/// @brief Message to start a remote gravity walk.
class ComputeChunkMsg : public CMessage_ComputeChunkMsg {
  ComputeChunkMsg() {} // not available
 public:
  int chunkNum;

  ComputeChunkMsg(int i) : chunkNum(i) {
  }
};

/// @brief Message for evaluating splits for the ORB domain decomposition
///
/// This message contains the splitting dimensions and values for an
/// ORB tree.  The size arrays pos and dim are set by arguments to
/// the new() operator, but should be equal to length.
class ORBSplittersMsg : public CMessage_ORBSplittersMsg{
public:
  /// Number of splits
  int length;
  /// Positions of splits
  double *pos;
  /// Dimension of splits
  char *dim;
  /// Callback for reduction of particle counts
  CkCallback cb;

  ORBSplittersMsg(int len, CkCallback callback): length (len), cb(callback) {}

};

/// Message for shuffling particles during domain decomposition
class ParticleShuffleMsg : public CMessage_ParticleShuffleMsg{
public:
    int nloads;
    int n;
    int nSPH;
    int nStar;
    double *loads;
    unsigned int *parts_per_phase;
    GravityParticle *particles;
    extraSPHData *pGas;
    extraStarData *pStar;
    ParticleShuffleMsg(int nload, int npart, int nsph, int nstar): 
      nloads(nload), n(npart), nSPH(nsph), nStar(nstar) {}
};

#ifdef PUSH_GRAVITY
#include "ckmulticast.h"

struct BucketMsg : public CkMcastBaseMsg, public CMessage_BucketMsg {
  GenericTreeNode *buckets;
  int numBuckets;
  ExternalGravityParticle *particles;
  int numParticles;
  int whichTreePiece;
};
#endif
    
/// Class to count added and deleted particles
class CountSetPart 
{
 public:
    int index;			/* chare index */
    int nAddGas;
    int nDelGas;
    int nAddDark;
    int nDelDark;
    int nAddStar;
    int nDelStar;

    void pup(PUP::er& p) {
	p | index;
	p | nAddGas;
	p | nDelGas;
	p | nAddDark;
	p | nDelDark;
	p | nAddStar;
	p | nDelStar;
	}
    };

/*
 * Multistepping routines
 *
 * Each major timestep can have MAXSUBSTEPS substeps where MAXSUBSTEPS
 * is a large power of 2 that can fit in an integer: 1 << MAXRUNG,
 * with MAXRUNG something like 30.
 * A given particle is on a "Rung" such that it takes 1 << Rung
 * substeps per major timestep.  That is, it's force is updated every
 * 1 << (MAXRUNG - Rung) smallest substeps.
 *
 * Needed routines:
 * DtToRung(): take an ideal timestep, the major timestep and
 * determine a rung.
 */

const int MAXRUNG = 30;
const int MAXSUBSTEPS = 1 << MAXRUNG;
const double MAXSUBSTEPS_INV = 1 / (double)MAXSUBSTEPS;

/// @brief Given a rung, return the number of substeps in one big step.
inline int RungToSubsteps(int iRung) {
  CkAssert(iRung <= MAXRUNG);
  return 1 << (MAXRUNG - iRung);
}

/// @brief Given the size of the big step, and a desired timestep,
/// return the rung of the largest timestep less than dTideal.
inline int DtToRung(double dDelta, double dTideal) {
  int iSteps = (int) ceil(dDelta/dTideal);

  int iRung = 0;
  iSteps--;
  while(iSteps > 0) {
    iRung++;
    iSteps >>= 1;
  }
  return iRung;
}

/// @brief Given the size of the big step, and a rung, return the
/// corresponding timestep size.
inline double RungToDt(double dDelta, int iRung) {
  return dDelta*RungToSubsteps(iRung)*MAXSUBSTEPS_INV;
}

/// @brief slot in MultistepLB to hold feedback phase load information
const int PHASE_FEEDBACK = MAXRUNG + 1;

/// @brief Pressure floor to force Jeans length to be larger than the
/// spatial resolution.

inline double PoverRhoFloorJeans(double dResolveJeans, GravityParticle *p)
{
    /*
     * Add pressure floor to keep Jeans Mass
     * resolved.  In comparison with Agertz et
     * al. 2009, dResolveJeans should be 3.0:
     * P_min = 3*G*max(h,eps)^2*rho^2
     * Note that G = 1 in our code
     */
#ifdef JEANSSOFTONLY
    double l2 = p->soft*p->soft;
#else
    double l2 = 0.25*p->fBall*p->fBall;
#ifdef JEANSSOFT
    double e2 = p->soft*p->soft; 
    if (l2 < e2) l2 = e2; /* Jeans scale can't be smaller than softening */
#endif
#endif
    return l2*dResolveJeans*p->fDensity;
}

/// @brief Adiabatic index to use with the Jeans pressure floor.
const double GAMMA_JEANS = 2.0;
const double GAMMA_NONCOOL = 5.0/3.0;


/// @brief Overall flow control of the simulation.
///
/// As well as controlling the overall flow of the simulation, the
/// constructors are the main entry points into the program.
/// The sequence of tasks is: read the simulation parameters (Main()),
/// read in the initial conditions (setupICs()), calculate the initial
/// forces (initialForces()), then iterate across timesteps and write
/// the final output (doSimulation()).
///
class Main : public CBase_Main {
	CkArgMsg *args;
	std::string basefilename;
        /// Save parameters for output
        OutputParams *pOutput;
    // NChilada file names used to generate the XML description
    CkVec<std::string> *NCgasNames;
    CkVec<std::string> *NCdarkNames;
    CkVec<std::string> *NCstarNames;
	/// globally finished IO
	CkCallback cbIO;
        /// Save file token for CkIO
        Ck::IO::File fIOFile;
	CProxy_Sorter sorter;
	int64_t nTotalParticles;
	int64_t nTotalSPH;
	int64_t nTotalDark;
	int64_t nTotalStar;
	/// Total Sink Particles
	int64_t nSink;
	int64_t nMaxOrderGas;  /* Maximum iOrders */
	int64_t nMaxOrderDark;
	int64_t nMaxOrder;
	
	double dTime;		/* Simulation time */
	double dTime0;		///< Simulation time at dStep = 0
	double dEcosmo;		/* variables for integrating
				   Lazer-Irvine eq. */
	double dUOld;
	double dTimeOld;
	PRM prm;		/* parameter parsing info */
	Parameters param; /* actual parameters */
	CkVec<double> vdOutTime; // Desired output times
	int iOut;
	/*
	 ** Tracking for frame dumping function
	 */
	int bDumpFrame;
	struct DumpFrameContext **df;
	int bIsRestarting;
        /// SPH Alpha has been read in.
        int bHaveAlpha;
	int bChkFirst;		/* alternate between 0 and 1 for checkpoint */
	double dSimStartTime;   // Start time for entire simulation
	int iStop;		/* indicate we're stopping the
				   simulation early */
	int64_t nActiveGrav;
	int64_t nActiveSPH;

#ifdef CUDA
          double localNodesPerReqDouble;
          double remoteNodesPerReqDouble;
          double remoteResumeNodesPerReqDouble;
          double localPartsPerReqDouble;
          double remotePartsPerReqDouble;
          double remoteResumePartsPerReqDouble;
#endif

#ifdef CHECK_TIME_WITHIN_BIGSTEP
       double wallTimeStart;
#endif

       /// @brief Hold wall clock timings of calculation phases for a
       /// given rung.
       class timing_fields {
       public:
           int count;           ///< number of times on this rung
           double tGrav;        ///< Gravity time
           double tuDot;        ///< Energy integration
           double tDD;          ///< Domain Decomposition
           double tLoadB;       ///< Load Balancing
           double tTBuild;      ///< Tree Building
           double tAdjust;      ///< Timestep adjustment
           double tEmergAdjust; ///< Emergency Timestep adjustment
           double tKick;        ///< Kick time
           double tDrift;       ///< Drift time
           double tCache;       ///< Cache teardown
       public:
           ///@brief Zero out fields
           void clear() {
               count = 0;
               tGrav = tuDot = tDD = tLoadB = tTBuild = tAdjust
                   = tEmergAdjust = tKick = tDrift = tCache = 0.0;
               }
           };
       
       CkVec<timing_fields> timings;  ///< One element for each rung.
       void writeTimings(int iStep);

#ifdef SELECTIVE_TRACING
       int monitorRung;
       int monitorStart;
       int numTraceIterations;
       int numSkipIterations;
       int numMaxTrace;
       int traceIteration;
       int traceState;
       bool projectionsOn;

       void turnProjectionsOn(int activeRung);
       void turnProjectionsOff();
#endif


public:

	Main(CkArgMsg* m);
	Main(CkMigrateMessage *m);

	void niceExit();
	void setupICs();
	void initialForces();
	void doSimulation();
	void restart(CkCheckpointStatusMsg *msg);
	void waitForGravity(const CkCallback &cb, double startTime,
            int activeRung);
        void advanceBigStep(int);
        void domainDecomp(int iPhase);
        void loadBalance(int iPhase);
        void buildTree(int iPhase);
        void startGravity(const CkCallback& cbGravity, int iActiveRung,
            double *startTime) ;
        void externalGravity(int iActiveRung);
        void updateuDot(int iActiveRung, const double duKick[],
            const double dStartTime[], int bUpdateState, int bAll);
        void kick(bool bClosing, int iActiveRung, int nextMaxRung,
            const CkCallback &cbGravity, double gravStartTime);
	int adjust(int iKickRung);
	void rungStats();
	void countActive(int activeRung);
        void emergencyAdjust(int iRung);
	void starCenterOfMass();
	void calcEnergy(double, double, const char *);
	void getStartTime();
	void getOutTimes();
	int bOutTime();
	void writeOutput(int iStep) ;
        void outputBinary(OutputParams& params, int bParaWrite,
                          const CkCallback& cb);
        void cbOpen(Ck::IO::FileReadyMsg *msg);
        void cbIOReady(Ck::IO::SessionReadyMsg *msg);
        void cbIOComplete(CkMessage *msg);
        void cbIOClosed(CkMessage *msg);
        std::string getNCNextOutput(OutputParams& params);
    void writeNCXML(std::string filename);
    void NCXMLattrib(ofstream *desc, CkVec<std::string> *names, std::string family);
	void updateSoft();
	void growMass(double dTime, double dDelta);
	void initSph();
	void initCooling();
	void initStarLog();
	int ReadASCII(char *extension, int nDataPerLine, double *dDataOut);
        void restartGas();
	void doSph(int activeRung, int bNeedDensity = 1);
	void AGORAfeedbackPreCheck(double dTime, double dDelta, double dTimeToSF);
	void FormStars(double dTime, double dDelta);
	void StellarFeedback(double dTime, double dDelta);
	void outputBlackHoles(double dTime);
	void SetSink();
	void FormSinks(double dTime, double dDelta, int iKickRung);
	void doSinks(double dTime, double dDelta, int iKickRung);
	int DumpFrameInit(double dTime, double dStep, int bRestart);
	void DumpFrame(double dTime, double dStep);
	int nextMaxRungIncDF(int nextMaxRung);
	void addDelParticles();
	void memoryStats();
	void memoryStatsCache();
	void pup(PUP::er& p);
	void liveVizImagePrep(liveVizRequestMsg *msg);
        void doSIDM(double dTime,double dDelta, int activeRung); /* SIDM */
};

/* IBM brain damage */
#undef hz
/// @brief Coefficients for the Fourier space part of the Ewald sum.
typedef struct ewaldTable {
  double hx,hy,hz;
  double hCfac,hSfac;
} EWT;

// jetley
class MissRecord;
class State;

///Remote Cell interaction lists for all tree levels
typedef struct OffsetNodeStruct
{
      GenericTreeNode *node;
      int offsetID;
}OffsetNode;

#if INTERLIST_VER > 0
/// @brief Queue of nodes to check for interactions.
typedef CkQ<OffsetNode> CheckList;
/// @brief Vector of nodes that are undecided at this level.
typedef CkVec<OffsetNode> UndecidedList;
/// @brief Vector of undecided lists, one for each level.
typedef CkVec<UndecidedList> UndecidedLists;
#endif

/// @brief Remote particles in an interaction list.
typedef struct particlesInfoR{
    ExternalGravityParticle* particles;
    int numParticles;
    Vector3D<cosmoType> offset;
#if defined CHANGA_REFACTOR_PRINT_INTERACTIONS || defined CHANGA_REFACTOR_WALKCHECK_INTERLIST || defined CUDA
    NodeKey key;
#endif
#if COSMO_DEBUG > 1
    GenericTreeNode *nd;
#endif
} RemotePartInfo;

/// @brief Local particles in an interaction list.
typedef struct particlesInfoL{
    GravityParticle* particles;
    int numParticles;
    Vector3D<cosmoType> offset;
#if defined CHANGA_REFACTOR_PRINT_INTERACTIONS || defined CHANGA_REFACTOR_WALKCHECK_INTERLIST || defined CUDA
    NodeKey key;
#endif
#if COSMO_DEBUG > 1 || defined CUDA
    GenericTreeNode *nd;
#endif
} LocalPartInfo;

/** @brief Data needed for the CkLoop intranode parallelization.
 *
 * This structure holds the data that needs to be passed to a free
 * processor so that it can calculate the gravity interactions.
 * Each attribute is a list so that multiple buckets can be operated on.
 */
typedef struct LoopParDataStruct {
  CkVec<GenericTreeNode*> lowNodes;  ///< Lowest node containing the
                                     ///  buckets to interact
  CkVec<int> bucketids;              ///< startBucket number
  CkVec<int> chunkids;               ///< remote walk chunk number
  CkVec<CkVec<OffsetNode> > clists;  ///< Cell interactions
  CkVec<CkVec<RemotePartInfo> > rpilists;  ///< Remote particle interactions
  CkVec<CkVec<LocalPartInfo> > lpilists;   ///< Local particle interactions
  TreePiece* tp;                           ///< Treepiece that owns this data
} LoopParData;



#ifdef CUDA
struct BucketActiveInfo{
  int start;
  int size;
};
#endif

class SmoothCompute;
#include "Compute.h"

#if INTERLIST_VER > 0 && defined CUDA
template<typename T> class GenericList;
#endif

/// @brief client that has requested a moment.
struct NonLocalMomentsClient {
  TreePiece *clientTreePiece;
  GenericTreeNode *clientNode;

  NonLocalMomentsClient() :
    clientTreePiece(NULL),
    clientNode(NULL)
  {}

  NonLocalMomentsClient(TreePiece *tp, GenericTreeNode *node) : 
    clientTreePiece(tp),
    clientNode(node)
  {}
};

/// @brief List of clients needing a particular moment
struct NonLocalMomentsClientList {
  GenericTreeNode *targetNode;
  CkVec<NonLocalMomentsClient> clients;

  NonLocalMomentsClientList() : 
    targetNode(NULL)
  {}

  NonLocalMomentsClientList(GenericTreeNode *node) :
    targetNode(node)
  {}

  void addClient(const NonLocalMomentsClient &cli){
    clients.push_back(cli);
  }
};

/// Fundamental structure that holds particle and tree data.
class TreePiece : public CBase_TreePiece {
   // jetley
   friend class PrefetchCompute;
   friend class GravityCompute;
   friend class SmoothCompute;
   friend class KNearestSmoothCompute;
   friend class ReSmoothCompute;
   friend class MarkSmoothCompute;
   friend class ListCompute;
   friend class NearNeighborState;
   friend class ReNearNeighborState;
   friend class MarkNeighborState;
   friend class BottomUpTreeWalk;
#if INTERLIST_VER > 0 && defined CUDA
   friend class DataManager;
   template<typename T> friend class GenericList;
#endif

   friend class RemoteTreeBuilder; 
   friend class LocalTreeBuilder; 

   /// @brief Walk for gravity prefetch
   TreeWalk *sTopDown;
   TreeWalk *twSmooth;
#if INTERLIST_VER > 0
   TreeWalk *sInterListWalk;
   // clearable, used for resumed walks
   State *sInterListStateRemoteResume;
#endif
   Compute *sGravity, *sPrefetch;
   SmoothCompute *sSmooth;
   
   Opt *sLocal, *sRemote, *sPref;
   Opt *optSmooth;

   State *sPrefetchState;
   /// Keeps track of the gravity walks over the local tree.
   State *sLocalGravityState, *sRemoteGravityState, *sSmoothState;
   typedef std::map<KeyType, CkVec<int>* > SmPartRequestType;
   // buffer of requests for smoothParticles.
   SmPartRequestType smPartRequests;

   CkVec<ActiveWalk> activeWalks;
   int completedActiveWalks; // XXX this should be part of the gravity
			     // walk state.
   int nCacheAccesses; // keep track of outstanding cache accesses to
		       // know when writebacks complete.  XXX this
		       // should be part of the smooth state
   
   /// number of active particles on the last active rung for load balancing
   unsigned int nPrevActiveParts;
   std::vector<double> savedPhaseLoad;
   std::vector<unsigned int> savedPhaseParticle;
   /// temporary accumulator for phase load information during domain decomposition
   std::vector<double> savedPhaseLoadTmp;
   /// temporary accumulator for phase particle counts during domain decomposition
   std::vector<unsigned int> savedPhaseParticleTmp;

   int memWithCache, memPostCache;  // store memory usage.
   int nNodeCacheEntries, nPartCacheEntries;  // store memory usage.

#ifdef PUSH_GRAVITY
   bool doMerge;
   bool createdSpanningTree;
   CProxySection_TreePiece allTreePieceSection;
   CkVec<GravityParticle> foreignParticles;
   CkVec<double> foreignParticleAccelerations;

   map<int,CkSectionInfo> cookieJar;
  
   BucketMsg *createBucketMsg();
   void unpackBuckets(BucketMsg *, GenericTreeNode *&foreignBuckets, int &numForeignBuckets);
   void calculateForces(GenericTreeNode *foreignBuckets, int numForeignBuckets);

#endif

 public:

#ifdef PUSH_GRAVITY
  void startPushGravity(int am, double myTheta);
  void recvPushBuckets(BucketMsg *);
  void recvPushAccelerations(CkReductionMsg *);
#endif

#if COSMO_PRINT_BK > 1
  State *getSRemoteGravityState(){ return sRemoteGravityState; }
  State *getSLocalGravityState(){ return sLocalGravityState; }
#endif
  void memCacheStats(const CkCallback &cb);
  void addActiveWalk(int iAwi, TreeWalk *tw, Compute *c, Opt *o, State *s);

  /// @brief Called when walk on the current TreePiece is done.
  void markWalkDone();
  /// @brief Called when walk on all TreePieces is done.
  void finishWalk();
  /// @brief Called when smooth walk on the current TreePiece is done.
  void markSmoothWalkDone();
  /// @brief Called when smooth walk on all TreePieces is done.
  void finishSmoothWalk();

  int getIndex() {
    return thisIndex;
  }

  /// @brief accumulate node interaction count for statistics
  void addToNodeInterRemote(int chunk, int howmany){
    nodeInterRemote[chunk] += howmany;
  }

  /// @brief accumulate particle interaction count for statistics
  void addToParticleInterRemote(int chunk, int howmany){
    particleInterRemote[chunk] += howmany;
  }

  /// @brief accumulate node interaction count for statistics
  void addToNodeInterLocal(int howmany){
    nodeInterLocal += howmany;
  }

  /// @brief accumulate particle interaction count for statistics
  void addToParticleInterLocal(int howmany){
    particleInterLocal += howmany;
  }

  /// Start prefetching the specfied chunk; prefetch compute
  /// calls startRemoteChunk() once chunk prefetch is complete
  void initiatePrefetch(int chunk);
  /// Start a new remote computation upon prefetch finished
  void startRemoteChunk();

  /// Return the number of particles on this TreePiece.
  int getNumParticles(){
    return myNumParticles;
  }

  /// Return the pointer to the particles on this TreePiece.
  GravityParticle *getParticles(){return myParticles;}


#ifdef CUDA
        // this variable holds the number of buckets active at
        // the start of an iteration
        // it is used to ascertain how many buckets still need to 
        // be processed via the stateReady function with regard
        // to their local and remote-no-resume walks. 
        // if all numActiveBuckets
        // have been processed but we still have leftover nodes/particles
        // in the list of interations to the sent to the gpu, we flush
        // the list
        int numActiveBuckets; 
        int myNumActiveParticles;
        // First and Last indices of GPU particle
        int FirstGPUParticleIndex;
        int LastGPUParticleIndex;
        int NumberOfGPUParticles;
        BucketActiveInfo *bucketActiveInfo;

        int getNumBuckets(){
        	return numBuckets;
        }

        void callFreeRemoteChunkMemory(int chunk);

        int getActiveRung(){ return activeRung; }
#ifdef HAPI_INSTRUMENT_WRS
        int getInstrumentId(){ return instrumentId; }
#endif
        // returns either all particles or only active particles,
        // depending on fraction of active particles to their
        // total count.
        int getDMNumParticles(){
          if(largePhase()){
            return myNumParticles;
          }
          else{
            return myNumActiveParticles; 
          }
        }

        int getNumActiveParticles(){
          return myNumActiveParticles;
        }

        void calculateNumActiveParticles(){ 
          myNumActiveParticles = 0;
          for(int i = 1; i <= myNumParticles; i++){
            if(myParticles[i].rung >= activeRung){
              myNumActiveParticles++;
            }
          }
        }

        bool largePhase(){
          return (1.0*myNumActiveParticles/myNumParticles) >= largePhaseThreshold;
        }

        void getDMParticles(CompactPartData *fillArray, int &fillIndex){
          NumberOfGPUParticles = 0;
          FirstGPUParticleIndex = fillIndex;//This is for the GPU Ewald
          if(largePhase()){
            for(int b = 0; b < numBuckets; b++){
              GenericTreeNode *bucket = bucketList[b];
              int buckstart = bucket->firstParticle;
              int buckend = bucket->lastParticle;
              GravityParticle *buckparts = bucket->particlePointer;
              bucket->bucketArrayIndex = fillIndex;
              for(int i = buckstart; i <= buckend; i++){
                fillArray[fillIndex] = buckparts[i-buckstart];
                fillIndex++;
              }
            }
          }
          else{
            for(int b = 0; b < numBuckets; b++){
              GenericTreeNode *bucket = bucketList[b];
              if(bucket->rungs < activeRung){
                continue;
              }
              BucketActiveInfo *binfo = &(bucketActiveInfo[b]);
              
              int buckstart = bucket->firstParticle;
              int buckend = bucket->lastParticle;
              GravityParticle *buckparts = bucket->particlePointer;

              binfo->start = fillIndex;
              for(int i = buckstart; i <= buckend; i++){
                if(buckparts[i-buckstart].rung >= activeRung){
                  fillArray[fillIndex] = buckparts[i-buckstart];
                  fillIndex++;
                }
              }
              binfo->size = fillIndex-binfo->start;
            }
          }
          //This is for the GPU Ewald
          if(FirstGPUParticleIndex == fillIndex){
            //This means no particle is on GPU
            FirstGPUParticleIndex = -1;
            LastGPUParticleIndex = -1;
            NumberOfGPUParticles = 0;
          }
          else{
            LastGPUParticleIndex = fillIndex - 1;
            NumberOfGPUParticles = LastGPUParticleIndex - FirstGPUParticleIndex + 1;
          }
        }

        bool isActive(int partNum){
          return myParticles[partNum].rung >= activeRung;
        }

        void clearMarkedBuckets(CkVec<GenericTreeNode *> &markedBuckets);
        void clearMarkedBucketsAll();

#ifdef HAPI_TRACE
        long long localNodeInteractions;
        long long localPartInteractions;
        long long remoteNodeInteractions;
        long long remotePartInteractions;
        long long remoteResumeNodeInteractions;
        long long remoteResumePartInteractions;
#endif

#ifdef HAPI_INSTRUMENT_WRS
        int instrumentId;

        double localNodeListConstructionTime;
        double remoteNodeListConstructionTime;
        double remoteResumeNodeListConstructionTime;
        double localPartListConstructionTime;
        double remotePartListConstructionTime;
        double remoteResumePartListConstructionTime;
        
        int nLocalNodeReqs;
        int nRemoteNodeReqs;
        int nRemoteResumeNodeReqs;
        int nLocalPartReqs;
        int nRemotePartReqs;
        int nRemoteResumePartReqs;

#endif

#endif

        void continueStartRemoteChunk(int chunk);
#ifdef CUDA
        void updateParticles(intptr_t data, int partIndex);
#endif
        void continueWrapUp();

#if INTERLIST_VER > 0
        void getBucketsBeneathBounds(GenericTreeNode *&node, int &start, int &end);
        void updateBucketState(int start, int end, int n, int chunk, State *state);
        void updateUnfinishedBucketState(int start, int end, int n, int chunk, State *state);
#endif

#if defined CHANGA_REFACTOR_WALKCHECK || defined CHANGA_REFACTOR_WALKCHECK_INTERLIST
        void addToBucketChecklist(int bucketIndex, NodeKey k){
          bucketcheckList[bucketIndex].insert(k);
          if(bucketIndex == TEST_BUCKET && thisIndex == TEST_TP)
            CkPrintf("[%d] add %ld\n", thisIndex, k);
        }
#endif

#if INTERLIST_VER > 0
        GenericTreeNode *getStartAncestor(int current, int previous, GenericTreeNode *dflt);
#endif
	/// \brief convert a key to a node using the nodeLookupTable
	inline GenericTreeNode *keyToNode(const Tree::NodeKey k){
          NodeLookupType::iterator iter = nodeLookupTable.find(k);
          if (iter != nodeLookupTable.end()) return iter->second;
          else return NULL;
        }

        GenericTreeNode *getBucket(int i){
          return bucketList[i];
        }

	/// Time read in from input file
	double dStartTime;

private:        
	// liveViz 
	liveVizRequestMsg * savedLiveVizMsg;

        LBStrategy foundLB;
        // jetley - saved first internal node
        Vector3D<float> savedCentroid;
        /// The phase for which we have just collected load balancing data.
        int iPrevRungLB;
        /// The phase for which we are about to do load balancing
        int iActiveRungLB;

	/// @brief Used to inform the mainchare that the requested operation has
	/// globally finished
	CkCallback callback;
	/// gravity globally finished
	CkCallback cbGravity;
	/// smooth globally finished
	CkCallback cbSmooth;
  CkCallback after_dd_callback;
	/// Total number of particles contained in this chare
	unsigned int myNumParticles;
	unsigned int numActiveParticles;
	/// Array with the particles in this chare
	GravityParticle* myParticles;
  int nbor_msgs_count_;
	/// Actual storage in the above array
	int nStore;
        /// Accelerations are initialized
        bool bBucketsInited;

  // Temporary location to hold the particles that have come from outside this
  // TreePiece. This is used in the case where we migrate the particles and
  // detect completion of migration using QD.
	std::vector<GravityParticle> myTmpShuffleParticle;
	std::vector<extraSPHData> myTmpShuffleSphParticle;
	std::vector<extraStarData> myTmpShuffleStarParticle;
  ParticleShuffleMsg* myShuffleMsg;
	/// Number of particles in my tree.  Can be different from
	/// myNumParticles when particles are created.
	int myTreeParticles;
 public:
	/// Total Particles in the simulation
	int64_t nTotalParticles;
	/// Total Gas Particles
	int64_t nTotalSPH;
	/// Total Dark Particles
	int64_t nTotalDark;
	/// Total Star Particles
	int64_t nTotalStar;
 private:
	/// Total gas in this chare
	unsigned int myNumSPH;
	/// Array with SPH particle data
	extraSPHData *mySPHParticles;
	/// Actual storage in the above array
	int nStoreSPH;
	/// Total stars in this chare
	unsigned int myNumStar;
	/// Array with star particle data
	extraStarData *myStarParticles;
	/// Actual storage in the above array
	int nStoreStar;
	/// MaxIOrder for output
	int64_t nMaxOrder;
        /// Start particle for reading
        int64_t nStartRead;
	/// particle count for output
	int myIOParticles;

	/// List of all the node-buckets in this TreePiece
	std::vector<GenericTreeNode *> bucketList;

	/// Array with sorted particles for domain decomposition (ORB)
	std::vector<GravityParticle> mySortedParticles;
	/// Array with sorted SPH data for domain decomposition (ORB)
	std::vector<extraSPHData> mySortedParticlesSPH;
	/// Array with sorted Star data for domain decomposition (ORB)
	std::vector<extraStarData> mySortedParticlesStar;
        /// Array with incoming particles messages for domain decomposition
	std::vector<ParticleShuffleMsg*> incomingParticlesMsg;
        /// How many particles have already arrived during domain decomposition
        int incomingParticlesArrived;
        /// Flag to acknowledge that the current TreePiece has already
        /// contributed to itself, prevents that acceptSortedParticles will
        /// change the TreePiece particles before the one belonging to someone
        /// else have been sent out
        bool incomingParticlesSelf;

	/// holds the total mass of the current TreePiece
	double piecemass;

	/// used to determine which coordinate we are outputting, while printing
	/// accelerations in ASCII format
	int cnt;
	/// used to determine if the x, y, z coordinates should be printed
	/// together while outputting accelerations in ASCII format
	int packed;

          // the index assigned by the CacheManager upon registration
          int localIndex;

	//unsigned int numSplitters;
        //SFC::Key* splitters;
	CProxy_TreePiece pieces;
	/// A proxy to the DataManager.
	//CProxy_DataManager dataManager; // unused...
	/// A local pointer to my DataManager.
	DataManager* dm;

	/// The counts of how many particles belonging to other
	/// TreePieces I currently hold
#ifndef REDUCTION_HELPER
	CkVec<int64_t> myBinCounts;
#endif
	std::vector<int> myBinCountsORB;
	/// My index in the responsibility array.
	int myPlace;
	/// The keys determined by the Sorter that separate me from my
	/// neighbors.
	SFC::Key leftSplitter, rightSplitter;

	std::string basefilename;
	// Bounding box of the entire simulation
	OrientedBox<float> boundingBox;
	unsigned iterationNo;
	/// The root of the global tree, always local to any chare
	GenericTreeNode* root;
	/// pool of memory to hold TreeNodes: makes allocation more efficient.
	NodePool *pTreeNodes;

	typedef std::map<NodeKey, CkVec<int>* >   MomentRequestType;
	/// Keep track of the requests for remote moments not yet satisfied.
	/// Used only during the tree construction.  This is a map
	/// from NodeKey to a vector of treepieces that have requested it.
	MomentRequestType momentRequests;

	/// Opening angle
	//double theta; -- moved as readonly
	/// Opening angle - monopole
	//double thetaMono; -- moved as readonly
        /// The current active mask for force computation in multistepping
        int activeRung;

	/// Periodic Boundary stuff
	int bPeriodic;
	int bComove;
	/// Background density of the Universe
	double dRhoFac;
	Vector3D<cosmoType> fPeriod;
	int nReplicas;
	int bEwald;		/* Perform Ewald */
	double fEwCut;
	double dEwhCut;
	EWT *ewt;
	int nMaxEwhLoop;
	int nEwhLoop;
#ifdef HEXADECAPOLE
	MOMC momcRoot;		/* complete moments of root */
#endif
        /// Have the Ewald h loop tables been calculated.
        bool bEwaldInited;

	int bGasCooling;
#ifndef COOLING_NONE
	clDerivsData *CoolData;
#endif
        /// indices of my newly formed particles in the starlog table
        std::vector<int> iSeTab;
	/// Setup for writing
	int nSetupWriteStage;
	int64_t nStartWrite;	// Particle number at which this piece starts
				// to write file.

	/// Map between Keys and TreeNodes, used to get a node from a key
	NodeLookupType nodeLookupTable;

	/// Number of nodes still missing before starting the real computation
	//u_int64_t prefetchWaiting;
	/// Array of keys that will be the root of the prefetching chunks
	Tree::NodeKey *prefetchRoots;
        /// Bounding box to use for prefetching
        OrientedBox<cosmoType> prefetchReq;

	/// number of chunks in which the tree will be chopped for prefetching
	int numChunks;

#if COSMO_STATS > 0
	u_int64_t myNumMACChecks;
	u_int64_t nodesOpenedLocal;
	u_int64_t nodesOpenedRemote;
	u_int64_t numOpenCriterionCalls;
#endif
	/// node interaction count for statistics
	u_int64_t nodeInterLocal;
	/// node interaction count for statistics
	u_int64_t *nodeInterRemote;
	/// particle interaction count for statistics
	u_int64_t particleInterLocal;
	/// particle interaction count for statistics
	u_int64_t *particleInterRemote;

	int nActive;		// number of particles that are active

	/// Size of bucketList, total number of buckets present
	unsigned int numBuckets;
#if INTERLIST_VER > 0
	/// Completed buckets, remote gravity walk
	int prevRemoteBucket;
	/// Completed buckets, local gravity walk
	int prevBucket;
#endif
	/// Used to start the Ewald computation for all buckets, one after the other
	unsigned int ewaldCurrentBucket;
	/// @brief Used as a placeholder while traversing the tree and computing
	/// forces. This should not exist in the cache version since all its
	/// information is duplicate from the correspondent TreeNode (in
	/// bucketList), and myParticles.
	/// @todo Eliminate the usage of this in the cache walk
	BucketGravityRequest *bucketReqs;

	/// Pointer to the instance of the local cache
        //CacheManager *localCache;

	// Entries used for the CacheManager
	EntryTypeGravityParticle gravityParticleEntry;
	EntryTypeSmoothParticle smoothParticleEntry;
	EntryTypeGravityNode gravityNodeEntry;

#if COSMO_DEBUG > 1 || defined CHANGA_REFACTOR_WALKCHECK || defined CHANGA_REFACTOR_WALKCHECK_INTERLIST
  ///A set for each bucket to check the correctness of the treewalk
  typedef std::vector< std::multiset<Tree::NodeKey> > DebugList;
  DebugList bucketcheckList;
#endif

  ///Variables added for ORB decomposition

  ///Used to denote the first call to evaluateParticleCounts within each phase
  bool firstTime;

  ///Temporary variable for the particle counts in each ORB box
  std::vector<int> tempBinCounts;

  ///Particles defining boundaries for ORB boxes
  std::list<GravityParticle *> orbBoundaries;

  ///Expected number of particles for each TreePiece after ORB decomposition
  int myExpectedCount;
  ///Expected number of SPH particles for each TreePiece after ORB
  ///decomposition
  int myExpectedCountSPH;
  ///Expected number of Star particles for each TreePiece after ORB
  ///decomposition
  int myExpectedCountStar;

  ///Level after which the local subtrees of all the TreePieces start
  unsigned int chunkRootLevel;

  ///Keeps track of the bounding boxes and dimensions along which they are splitted, starting from
  ///root to chunkRootLevel
  OrientedBox<float>* boxes;
  char* splitDims;

  ///Phase of ORB decomposition: the number of boxes double in each phase till they are equal to the number
  ///of TreePieces
  int phase;

  double myTotalMass;

 #if INTERLIST_VER > 0


 public:
#endif
  // called when a chunk has been used completely (chunkRemaining[chunk] == 0)
  void finishedChunk(int chunk);

  ///Array of comparison function pointers
  bool (*compFuncPtr[3])(GravityParticle,GravityParticle);
  double tmpTime;
  double totalTime;
 public:

#ifdef SPCUDA
  EwaldData *h_idata;
  CkCallback *cbEwaldGPU;
#endif
  void EwaldGPU(); 
  void EwaldGPUComplete();

#if COSMO_DEBUG > 1 || defined CHANGA_REFACTOR_WALKCHECK || defined CHANGA_REFACTOR_WALKCHECK_INTERLIST
  ///This function checks the correctness of the treewalk
  void checkWalkCorrectness();
  void combineKeys(Tree::NodeKey key,int bucket);
#endif

  /* DEBUGGING */
  void quiescence();
  /*END DEBUGGING */

  NodeLookupType &getNodeLookupTable() {
	  return nodeLookupTable;
  }

  int getNumNodes() {
	  return nodeLookupTable.size();
  }

  /// delete treenodes if allocated
  void deleteTree() {
    if(pTreeNodes != NULL) {
        delete pTreeNodes;
        pTreeNodes = NULL;
        root = NULL;
        nodeLookupTable.clear();
        }
    else {
        if (root != NULL) {
            root->fullyDelete();
            delete root;
            root = NULL;
            nodeLookupTable.clear();
            }
        }
    }

	/// Compute all the moments for the nodes that are NonLocal, so that
	/// during the tree traversal, they contain useful information to decide
	/// whether to open or not.
	void calculateRemoteMoments(GenericTreeNode* node);

	/// Checks that every particle is indeed included in its bucket node
	/// (may not be true due to truncation of the last two bits while
	/// generating the 63 bit keys.
	void checkTree(GenericTreeNode* node);

	/// Given a node, check who is the first owner and the last owner of it.
	/// It assumes that there are splitters, and that there is an ordering
	/// of them across chares. It works for SFC ordering.
	bool nodeOwnership(const Tree::NodeKey nkey, int &firstOwner, int &lastOwner);


	/// Initialize all the buckets for the tree walk
	/// @TODO: Eliminate this redundant copy!
	void initBuckets();
	template <class Tsmooth>
	void initBucketsSmooth(Tsmooth tSmooth);
	void smoothNextBucket();
	void reSmoothNextBucket();
	void markSmoothNextBucket();
	void smoothBucketComputation();
	/** @brief Start the treewalk for the next bucket among those belonging
	 * to me. The buckets are simply ordered in a vector.
	 */
	void startNextBucket();
	/** @brief Start a full step of bucket computation, it sends a message
	 * to trigger nextBucket() which will loop over all the buckets.
	 */
	void doAllBuckets();
	void reconstructNodeLookup(GenericTreeNode *node);
	//void rebuildSFCTree(GenericTreeNode *node,GenericTreeNode *parent,int *);

public:
 TreePiece() : pieces(thisArrayID), root(0),
            iPrevRungLB (-1), sTopDown(0), sGravity(0),
            sPrefetch(0), sLocal(0), sRemote(0), sPref(0), sSmooth(0), 
            nPrevActiveParts(0) {
	  dm = NULL;
	  foundLB = Null; 
	  iterationNo=0;
	  usesAtSync = true;
	  pTreeNodes = NULL;
	  bucketReqs=NULL;
	  nCacheAccesses = 0;
	  memWithCache = 0;
	  memPostCache = 0;
	  nNodeCacheEntries = 0;
	  nPartCacheEntries = 0;
	  completedActiveWalks = 0;
	  myPlace = -1;
	  nSetupWriteStage = -1;
    //openingDiffCount=0;
    chunkRootLevel=0;
    //splitters = NULL;
#if COSMO_STATS > 0
	  nodesOpenedLocal = 0;
	  nodesOpenedRemote = 0;
	  nodeInterLocal = 0;
	  particleInterLocal = 0;

	  numOpenCriterionCalls=0;

	  piecemass = 0.0;
#endif
	  packed=0;			/* output accelerations in x,y,z
					   packed format */
	  cnt=0;			/* Counter over x,y,z when not packed */
	  particleInterRemote = NULL;
	  nodeInterRemote = NULL;

#if INTERLIST_VER > 0
	  sInterListWalk = NULL;
#endif
#ifdef CUDA
          numActiveBuckets = -1;
#ifdef HAPI_TRACE
          localNodeInteractions = 0;
          localPartInteractions = 0;
          remoteNodeInteractions = 0;
          remotePartInteractions = 0;
          remoteResumeNodeInteractions = 0;
          remoteResumePartInteractions = 0;
#endif
#endif

	  tmpTime=0.0;
	  totalTime=0.0;
	  // temporarely set to -1, it will updated after the tree is built
	  numChunks=-1;
	  prefetchRoots = NULL;
	  ewt = NULL;
	  nMaxEwhLoop = 100;
          bEwaldInited = false;

          incomingParticlesMsg.clear();
          incomingParticlesArrived = 0;
          incomingParticlesSelf = false;

          myParticles = NULL;
          mySPHParticles = NULL;
          myStarParticles = NULL;
	  myNumParticles = myNumSPH = myNumStar = 0;
	  nStore = nStoreSPH = nStoreStar = 0;
          bBucketsInited = false;
	  myTreeParticles = -1;
	  orbBoundaries.clear();
	  boxes = NULL;
	  splitDims = NULL;
	  bGasCooling = 0;

#ifdef PUSH_GRAVITY
          createdSpanningTree = false;
#endif

          localTreeBuildComplete = false;
	}

    TreePiece(CkMigrateMessage* m): pieces(thisArrayID) {

	  usesAtSync = true;
	  //localCache = NULL;
	  dm = NULL;
	  bucketReqs = NULL;
	  numChunks=-1;
	  nCacheAccesses = 0;
	  memWithCache = 0;
	  memPostCache = 0;
	  nNodeCacheEntries = 0;
	  nPartCacheEntries = 0;
	  completedActiveWalks = 0;
	  prefetchRoots = NULL;
	  //remaining Chunk = NULL;
          ewt = NULL;
          bEwaldInited = false;
	  root = NULL;
	  pTreeNodes = NULL;

      sTopDown = 0;
	  sGravity = NULL;
	  sPrefetch = NULL;
	  sSmooth = NULL;
#if INTERLIST_VER > 0
	  sInterListWalk = NULL;
#endif

          incomingParticlesMsg.clear();
          incomingParticlesArrived = 0;
          incomingParticlesSelf = false;

	  cnt=0;
	  nodeInterRemote = NULL;
          particleInterRemote = NULL;

	  orbBoundaries.clear();
	  boxes = NULL;
	  splitDims = NULL;
          bBucketsInited = false;
	  myTreeParticles = -1;


          localTreeBuildComplete = false;
	}

        private:
        void freeWalkObjects();
	void allocateStars() {
	    nStoreStar = (int) (myNumStar*(1.0 + dExtraStore));
	    // Stars tend to form out of gas, so make sure there is
	    // enough room.
	    nStoreStar += 12 + (int) (myNumSPH*dExtraStore);
	    myStarParticles = new extraStarData[nStoreStar];
	    }

        public:
	~TreePiece() {
	  if (verbosity>1) ckout <<"Deallocating treepiece "<<thisIndex<<endl;
	  if(nStore > 0) delete[] myParticles;
	  if(nStoreSPH > 0) delete[] mySPHParticles;
	  if(nStoreStar > 0) delete[] myStarParticles;
	  delete[] nodeInterRemote;
	  delete[] particleInterRemote;
	  delete[] bucketReqs;
          delete[] ewt;

	  deleteTree();

	  if(boxes!= NULL ) delete[] boxes;
	  if(splitDims != NULL) delete[] splitDims;

#ifndef COOLING_NONE
	  if(bGasCooling)
	      CoolDerivsFinalize(CoolData);
#endif
          if (verbosity>1) ckout <<"Finished deallocation of treepiece "<<thisIndex<<endl;
	}

	void setPeriodic(int nReplicas, Vector3D<cosmoType> fPeriod, int bEwald,
			 double fEwCut, double fEwhCut, int bPeriod,
                         int bComove, double dRhoFac);
	void BucketEwald(GenericTreeNode *req, int nReps,double fEwCut);
	void EwaldInit();
	void calculateEwald(dummyMsg *m);
  void calculateEwaldUsingCkLoop(dummyMsg *msg, int yield_num);
  void callBucketEwald(int id);
  void doParallelNextBucketWork(int id, LoopParData* lpdata);
	void initCoolingData(const CkCallback& cb);
	// Scale velocities (needed to convert to canonical momenta for
	// comoving coordinates.)
	void velScale(double dScale, const CkCallback& cb);

	/// @brief Load I.C. from NChilada file
        /// @param dTuFac conversion factor from temperature to
        /// internal energy
        void loadNChilada(const std::string& filename, const double dTuFac,
                          const CkCallback& cb);
        void readFloatBinary(OutputParams& params, int bParaRead,
            const CkCallback& cb);
	/// @brief Load I.C. from Tipsy file
        /// @param filename tipsy file
        /// @param dTuFac conversion factor from temperature to
        /// internal energy
        /// @param bDoublePos Positions are in double precision
        /// @param bDoubleVel Velocities are in double precision
	void loadTipsy(const std::string& filename, const double dTuFac,
                       const bool bDoublePos,
                       const bool bDoubleVel,
		       const CkCallback& cb);
        /// @brief read a tipsy array file (binary or ascii)
        void readTipsyArray(OutputParams& params, const CkCallback& cb);
        void resetMetals(const CkCallback& cb);
        void getMaxIOrds(const CkCallback& cb);
        void RestartEnergy(double dTuFac, const CkCallback& cb);
        void findTotalMass(const CkCallback &cb);
        void recvTotalMass(CkReductionMsg *msg);

	// Write a Tipsy file
	void writeTipsy(Tipsy::TipsyWriter& w,
			const double dvFac, // scale velocities
			const double duTFac, // convert temperature
                        const bool bDoublePos,
                        const bool bDoubleVel,
			const int bCool);
	// Find position in the file to start writing
	void setupWrite(int iStage, u_int64_t iPrevOffset,
			const std::string& filename, const double dTime,
			const double dvFac, const double duTFac,
                        const bool bDoublePos,
                        const bool bDoubleVel,
			const int bCool, const CkCallback& cb);
	// control parallelism in the write
	void parallelWrite(int iPass, const CkCallback& cb,
			   const std::string& filename, const double dTime,
			   const double dvFac, // scale velocities
			   const double duTFac, // convert temperature
                           const bool bDoublePos,
                           const bool bDoubleVel,
			   const int bCool);
	// serial output
	void serialWrite(u_int64_t iPrevOffset, const std::string& filename,
			 const double dTime,
			 const double dvFac, const double duTFac,
                         const bool bDoublePos,
                         const bool bDoubleVel,
			 const int bCool, const CkCallback& cb);
	// setup for serial output
	void oneNodeWrite(int iIndex,
			       int iOutParticles,
			       int iOutSPH,
			       int iOutStar,
			       GravityParticle *particles, // particles to
						     // write
			       extraSPHData *pGas, // SPH data
			       extraStarData *pStar, // Star data
			       int *piSPH, // SPH data offsets
			       int *piStar, // Star data offsets
			       const u_int64_t iPrevOffset,
			       const std::string& filename,  // output file
			       const double dTime,      // time or expansion
			       const double dvFac,  // velocity conversion
			     const double duTFac, // temperature
						  // conversion
                           const bool bDoublePos,
                           const bool bDoubleVel,
			  const int bCool,
			  const CkCallback &cb);
	// Reorder for output
        void reOrder(int64_t nMaxOrder, const CkCallback& cb);
	// move particles around for output
	void ioShuffle(CkReductionMsg *msg);
	void ioAcceptSortedParticles(ParticleShuffleMsg *);
        /// @brief Set the load balancing data after a restart from
        /// checkpoint.
        void resetObjectLoad(const CkCallback& cb);
	// Assign keys after loading tipsy file and finding Bounding box
	void assignKeys(CkReductionMsg* m);
	void evaluateBoundaries(SFC::Key* keys, const int n, int isRefine, const CkCallback& cb);
	void unshuffleParticles(CkReductionMsg* m);
	void acceptSortedParticles(ParticleShuffleMsg *);
  void shuffleAfterQD();
  void unshuffleParticlesWoDD(const CkCallback& cb);
  void acceptSortedParticlesFromOther(ParticleShuffleMsg *);
  void setNumExpectedNeighborMsgs();

  /*****ORB Decomposition*******/
  void initORBPieces(const CkCallback& cb);
  void initBeforeORBSend(unsigned int myCount, unsigned int myCountGas,
			 unsigned int myCountStar,
			 const CkCallback& cb, const CkCallback& cback);
  void sendORBParticles();
  void acceptORBParticles(const GravityParticle* particles, const int n);
  void acceptORBParticles(const GravityParticle* particles, const int n,
			  const extraSPHData *pGas, const int nGasIn,
			  const extraStarData *pStar, const int nStarIn);
  void finalizeBoundaries(ORBSplittersMsg *splittersMsg);
  void evaluateParticleCounts(ORBSplittersMsg *splittersMsg);
  /*****************************/

  void kick(int iKickRung, double dDelta[MAXRUNG+1], int bClosing,
	    int bNeedVPred, int bGasIsothermal, double dMaxEnergy, double duDelta[MAXRUNG+1],
        double gammam1, double dThermalCondSatCoeff,
        double dMultiPhaseMaxTime, double dMultiPhaseMinTemp, double dEvapCoeff, const CkCallback& cb);
  void drift(double dDelta, int bNeedVPred, int bGasIsothermal, double dvDelta,
             double duDelta, int nGrowMass, bool buildTree, double dMaxEnergy,
	     const CkCallback& cb);
  void initAccel(int iKickRung, const CkCallback& cb);
#ifdef COOLING_MOLECULARH
  void distribLymanWerner(const CkCallback& cb);
#endif /*COOLING_MOLECULARH*/

  void applyFrameAcc(int iKickRung, Vector3D<double> frameAcc, const CkCallback& cb);
/**
 * @brief Apply an external gravitational force
 * @param activeRung The rung to apply the force.
 * @param exGravParams Parameters of the external force
 * @param cb Callback function
 */
  void externalGravity(int activeRung, const ExternalGravity exGrav,
                       const CkCallback& cb);
/**
 * Adjust timesteps of active particles.
 * @param iKickRung The rung we are on.
 * @param bEpsAccStep Use sqrt(eps/acc) timestepping
 * @param bGravStep Use sqrt(r^3/GM) timestepping
 * @param bSphStep Use Courant condition
 * @param bViscosityLimitdt Use viscosity in Courant condition
 * @param dEta Factor to use in determing timestep
 * @param dEtaCourant Courant factor to use in determing timestep
 * @param dEtauDot Factor to use in uDot based timestep
 * @param dDiffCoeff Diffusion coefficent
 * @param dEtaDiffusion Factor to use in diffusion based timestep
 * @param dDelta Base timestep
 * @param dAccFac Acceleration scaling for cosmology
 * @param dCosmoFac Cosmo scaling for Courant
 * @param dhMinOverSoft minimum smoothing parameter.
 * @param dResolveJeans multiple of Jeans length to be resolved.
 * @param bDoGas We are calculating gas forces.
 * @param cb Callback function reduces currrent maximum rung
 */
  void adjust(int iKickRung, int bEpsAccStep, int bGravStep,
	      int bSphStep, int bViscosityLimitdt,
	      double dEta, double dEtaCourant, double dEtauDot,
              double dDiffCoeff, double dEtaDiffusion,
	      double dDelta, double dAccFac,
	      double dCosmoFac, double dhMinOverSoft,
              double dResolveJeans,
	      int bDoGas,
	      const CkCallback& cb);
  /**
   * @brief Truncate the highest rung
   * @param iCurrMaxRung new maximum rung.
   * @param cb callback.
   */
  void truncateRung(int iCurrMaxRung, const CkCallback& cb);
  void rungStats(const CkCallback& cb);
  void countActive(int activeRung, const CkCallback& cb);
  /// @brief count total number of particles of given type
  void countType(int iType, const CkCallback& cb);
  void outputBlackHoles(const std::string& pszFileName, double dvFac,
                        long lFPos, const CkCallback &cb);
  /// @brief set sink type based on formation time.
  void SetSink(double dSinkMassMin, const CkCallback &cb);
  /// @brief set sink timesteps.
  void SinkStep(int iCurrSinkRung, int iKickRung, const CkCallback &cb);
  void formSinks(int bJeans, double dJConst2, int bDensity,
		 double dDensityCut, double dTime, int iKickRung, int bSimple,
		 const CkCallback &cb);
  void emergencyAdjust(int iRung, double dDelta, double dDeltaThresh,
		       const CkCallback &cb);
  void assignDomain(const CkCallback& cb);
  void starCenterOfMass(const CkCallback& cb);
  void calcEnergy(const CkCallback& cb);
  /// add new particle
  void newParticle(GravityParticle *p);
  void adjustTreePointers(GenericTreeNode *node, GravityParticle *newParts);
  /// Count add/deleted particles, and compact main particle storage.
  void colNParts(const CkCallback &cb);
  /// Assign iOrders to recently added particles.
  void newOrder(const NewMaxOrder *nStarts, const int n, const CkCallback &cb) ;
  
  /// Update total particle numbers
  void setNParts(int64_t _nTotalSPH, int64_t _nTotalDark,
		 int64_t _nTotalStar, const CkCallback &cb);
  /// @brief Set the gravitational softening on all particles
  void setSoft(const double dSoft, const CkCallback &cb);
	void physicalSoft(const double dSoftMax, const double dFac,
			  const int bSoftMaxMul, const CkCallback& cb);
	void growMass(int nGrowMass, double dDeltaM, const CkCallback& cb);
	void InitEnergy(double dTuFac, double z, double dTime, double gammam1,
			const CkCallback& cb);
	void updateuDot(int activeRung, double duDelta[MAXRUNG+1],
			double dStartTime[MAXRUNG+1], int bCool, int bAll,
			int bUpdateState, double gammam1, const CkCallback& cb);
	void ballMax(int activeRung, double dFac, const CkCallback& cb);
	void sphViscosityLimiter(int bOn, int activeRung, const CkCallback& cb);
    void getAdiabaticGasPressure(double gamma, double gammam1, double dTuFac, double dThermalCondCoeff,
        double dThermalCond2Coeff, double dThermalCondSatCoeff, double dThermalCond2SatCoeff,
        double dEvapMinTemp, double dDtCourantFac, const CkCallback &cb);
    void getCoolingGasPressure(double gamma, double gammam1, double dThermalCondCoeff,
        double dThermalCond2Coeff, double dThermalCondSatCoeff, double dThermalCond2SatCoeff,
        double dEvapMinTemp, double dDtCourantFac, double dResolveJeans, const CkCallback &cb);
#ifdef SPLITGAS
	void SplitGas(double dInitGasMass, const CkCallback& cb);
#endif
	inline COOL* Cool() {return dm->Cool;}
	/// @brief initialize random seed for star formation
	void initRand(int iRand, const CkCallback &cb);
	void FormStars(Stfm param, double dTime, double dDelta, double dCosmoFac,
		       const CkCallback& cb);
	void flushStarLog(const CkCallback& cb);
        void Feedback(const Fdbk &fb, double dTime, double dDelta,
                      const CkCallback& cb);
	void massMetalsEnergyCheck(int bPreDist, const CkCallback& cb);
	void SetTypeFromFileSweep(int iSetMask, char *file,
	   struct SortStruct *ss, int nss, int *pniOrder, int *pnSet);
	void setTypeFromFile(int iSetMask, char *file, const CkCallback& cb);
	void getCOM(const CkCallback& cb, int bLiveViz);
	void getCOMByType(int iType, const CkCallback& cb, int bLiveViz);
	void DumpFrame(InDumpFrame in, const CkCallback& cb, int liveVizDump) ;
	void liveVizDumpFrameInit(liveVizRequestMsg *msg);
	void setProjections(int bOn);

	/// \brief Charm entry point to build the tree (called by Main).
#ifdef PUSH_GRAVITY
	void buildTree(int bucketSize, const CkCallback& cb, bool merge);
#else
	void buildTree(int bucketSize, const CkCallback& cb);
#endif

	/// \brief Real tree build, independent of other TreePieces.
	void startOctTreeBuild(CkReductionMsg* m);
  void recvBoundary(SFC::Key key, NborDir dir);
	void recvdBoundaries(CkReductionMsg* m);

  /********ORB Tree**********/
  //void receiveBoundingBoxes(BoundingBoxes *msg);
  void startORBTreeBuild(CkReductionMsg* m);
  OrientedBox<float> constructBoundingBox(GenericTreeNode *node,int level, int numChild);
  void buildORBTree(GenericTreeNode * node, int level);
  /**************************/

  /// When the node is found to be NULL, forward the request
	bool sendFillReqNodeWhenNull(CkCacheRequestMsg<KeyType> *msg);
	/// Request the moments for this node.
	void requestRemoteMoments(const Tree::NodeKey key, int sender);
	void receiveRemoteMoments(const Tree::NodeKey key, Tree::NodeType type,
    int firstParticle, int numParticles, int remIdx,
    const MultipoleMoments& moments, const OrientedBox<double>& box,
            const OrientedBox<double>& boxBall,
            const unsigned int iParticleTypes, const int64_t nSPH);

	/// Entry point for the local computation: for each bucket compute the
	/// force that its particles see due to the other particles hosted in
	/// this TreePiece. The opening angle theta has already been passed
	/// through startGravity().  This function just calls doAllBuckets().
	void calculateGravityLocal();
	/// Do some minor preparation for the local walkk then
	/// calculateGravityLocal().
	void commenceCalculateGravityLocal();

	/// Entry point for the remote computation: for each bucket compute the
	/// force that its particles see due to the other particles NOT hosted
	/// by this TreePiece, and belonging to a subset of the global tree
	/// (specified by chunkNum).
	void calculateGravityRemote(ComputeChunkMsg *msg);
  void executeCkLoopParallelization(LoopParData *lpdata, int startbucket,
      int yield_num, int chunkNum, State* gravityState);
  int doBookKeepingForTargetActive(int curbucket, int end, int chunkNum, bool
    updatestate, State* gravityState);
  int doBookKeepingForTargetInactive(int chunkNum, bool updatestate,
    State* gravityState);

	/// As above but for the Smooth operation
	void calculateSmoothLocal();
	void calculateReSmoothLocal();
	void calculateMarkSmoothLocal();
	void nextBucketSmooth(dummyMsg *msg);
	void nextBucketReSmooth(dummyMsg *msg);
	void nextBucketMarkSmooth(dummyMsg *msg);
#if INTERLIST_VER > 0
#if 0
  void calculateForceRemoteBucket(int bucketIndex, int chunk);
  void calculateForceLocalBucket(int bucketIndex);
#endif
#endif

  /// @brief Start a tree based gravity computation.
  /// @param am the active rung for the computation
  /// @param theta the opening angle
  /// @param cb the callback to use after all the computation has finished
  void startGravity(int am, double myTheta, const CkCallback& cb);
  /// Setup utility function for all the smooths.  Initializes caches.
  void setupSmooth();
  /// Start a tree based smooth computation.
  /// @param p parameters, including the type, of the smooth.
  /// @param iLowhFix true if a minimum h/softening is used.
  /// @param nSmooth Number of particles to smooth over.
  /// @param dfBall2OverSoft2 square of minimum ratio of h to softening.
  /// @param cb the callback to use after all the computation has finished
  /// @reference SmoothParams
  void startSmooth(SmoothParams *p, int iLowhFix, int nSmooth,
		   double dfBall2OverSoft2, const CkCallback &cb);
  void startReSmooth(SmoothParams *p, const CkCallback& cb);
  void startMarkSmooth(SmoothParams *p, const CkCallback& cb);

  void finishNodeCache(const CkCallback& cb);

    /// @brief Retrieve the remote node, goes through the cache if present
    GenericTreeNode* requestNode(int remoteIndex, Tree::NodeKey lookupKey,
                                 int chunk, int reqID, int awi, void *source);
	/// @brief Receive a request for Nodes from a remote processor, copy the
	/// data into it, and send back a message.
	void fillRequestNode(CkCacheRequestMsg<KeyType> *msg);
	/** @brief Receive the node from the cache as following a previous
	 * request which returned NULL, and continue the treewalk of the bucket
	 * which requested it with this new node.
	*/
#if 0
	void receiveNode(GenericTreeNode &node, int chunk, unsigned int reqID);
#endif
	/// @brief Find the key in the KeyTable, and copy the node over the passed pointer
	const GenericTreeNode* lookupNode(Tree::NodeKey key);
	/// Find the particles starting at "begin", and return a pointer to it
	const GravityParticle* lookupParticles(int begin);

	/// @brief Check if we have done with the treewalk on a specific bucket,
	/// and if we have, check also if we are done with all buckets
	void finishBucket(int iBucket);

	/** @brief Routine which does the tree walk on non-local nodes. It is
	 * called back for every incoming node (which are those requested to the
	 * cache during previous treewalks), and continue the treewalk from
	 * where it had been interrupted. It will possibly made other remote
	 * requests.
	 */
	void cachedWalkBucketTree(GenericTreeNode* node, int chunk, int reqID);

#if INTERLIST_VER > 0
	void calculateForcesNode(OffsetNode node, GenericTreeNode *myNode,
				 int level,int chunk);
	void calculateForces(OffsetNode node, GenericTreeNode *myNode,
			     int level,int chunk);
#endif

    ExternalGravityParticle *requestParticles(Tree::NodeKey key,int chunk,
                                              int remoteIndex,int begin,
                                              int end,int reqID, int awi,
                                              void *source);
    GravityParticle *requestSmoothParticles(Tree::NodeKey key, int chunk,
                                            int remoteIndex, int begin,int end,
                                            int reqID, int awi, void *source);
	void fillRequestParticles(CkCacheRequestMsg<KeyType> *msg);
	void fillRequestSmoothParticles(CkCacheRequestMsg<KeyType> *msg);
	void flushSmoothParticles(CkCacheFillMsg<KeyType> *msg);
	void processReqSmoothParticles();

	void getParticleInfoForLB(int64_t active_part, int64_t total_part);
        void startlb(const CkCallback &cb, int activeRung);
  void setTreePieceLoad(int activeRung);
  void populateSavedPhaseData(int phase, double tpload, unsigned int activeparts);
  bool havePhaseData(int phase);
  void savePhaseData(std::vector<double> &loads, std::vector<unsigned int>
    &parts_per_phase, double* shuffleloads, unsigned int *shuffleparts,
    int shufflelen);
	void ResumeFromSync();

	void outputASCII(OutputParams& params, int bParaWrite,
			 const CkCallback& cb);
	void oneNodeOutVec(OutputParams& params, Vector3D<double>* avOut,
			   int nPart, int iIndex, int bDone,
			   const CkCallback& cb) ;
	void oneNodeOutArr(OutputParams& params, double* adOut,
			   int nPart, int iIndex, int bDone,
			   const CkCallback& cb) ;
	void oneNodeOutIntArr(OutputParams& params, int *aiOut,
                              int nPart, int iIndex, const CkCallback& cb);
        void outputBinaryStart(OutputParams& params, int64_t nStart,
                               const CkCallback& cb);
	void outputBinary(Ck::IO::Session, OutputParams& params);
        void minmaxNCOut(OutputParams& params, const CkCallback& cb);

	void outputStatistics(const CkCallback& cb);
	/// Collect the total statistics from the various chares
        void collectStatistics(const CkCallback &cb);

        /** @brief Entry method used to split the processing of all the buckets
         * in small pieces. It calls startNextBucket() _yieldPeriod number of
         * times, and then it returns to the scheduler after enqueuing a message
         * for itself.  After each startNextBucket() call, the state
         * is checked for completed walks (and forces calculated, and
         * finishBucket() is called to clean up.
	 */
        void nextBucket(dummyMsg *m);
  void nextBucketUsingCkLoop(dummyMsg *m);

	void report();
	void printTreeViz(GenericTreeNode* node, std::ostream& os);
	void printTree(GenericTreeNode* node, std::ostream& os);
	void pup(PUP::er& p);

        // jetley
        // need this in TreeWalk
        GenericTreeNode *getRoot() {return root;}
        // need this in Compute
	inline Vector3D<cosmoType> decodeOffset(int reqID) {
	    int offsetcode = reqID >> 22;
	    int x = (offsetcode & 0x7) - 3;
	    int y = ((offsetcode >> 3) & 0x7) - 3;
	    int z = ((offsetcode >> 6) & 0x7) - 3;

	    Vector3D<cosmoType> offset(x*fPeriod.x, y*fPeriod.y, z*fPeriod.z);

	    return offset;
	    }

        void receiveNodeCallback(GenericTreeNode *node, int chunk, int reqID, int awi, void *source);
        void receiveParticlesCallback(ExternalGravityParticle *egp, int num, int chunk, int reqID, Tree::NodeKey &remoteBucket, int awi, void *source);
        void receiveParticlesFullCallback(GravityParticle *egp, int num, int chunk, int reqID, Tree::NodeKey &remoteBucket, int awi, void *source);
        void doAtSync();

        void balanceBeforeInitialForces(const CkCallback &cb);

        // For merging of remote moment requests
        // before sending messages during tree building
        public:
        void sendRequestForNonLocalMoments(GenericTreeNode *pickedNode);
        void mergeNonLocalRequestsDone();
        //void addTreeBuildMomentsClient(GenericTreeNode *targetNode, TreePiece *client, GenericTreeNode *clientNode);
        std::map<NodeKey,NonLocalMomentsClientList>::iterator createTreeBuildMomentsEntry(GenericTreeNode *pickedNode);


        private:
        // XXX - hashtable instead of map
        std::map<NodeKey,NonLocalMomentsClientList> nonLocalMomentsClients;
        bool localTreeBuildComplete;
        int getResponsibleIndex(int first, int last);
        

        GenericTreeNode *boundaryParentReady(GenericTreeNode *parent);
        void accumulateMomentsFromChild(GenericTreeNode *parent, GenericTreeNode *child);

        void deliverMomentsToClients(GenericTreeNode *);
        void deliverMomentsToClients(const std::map<NodeKey,NonLocalMomentsClientList>::iterator &it);
        void treeBuildComplete();
        void processRemoteRequestsForMoments();
        void sendParticlesDuringDD(bool withqd);
        void mergeAllParticlesAndSaveCentroid();
        bool otherIdlePesAvail();

};

/// @brief Class for shadow arrays to avoid reduction conflicts in
/// overlapping liveViz, SPH and gravity.
class LvArray : public CBase_LvArray {
 public:
    LvArray() {}
    LvArray(CkMigrateMessage* m) {}
    } ;

int decodeReqID(int);
int encodeOffset(int reqID, int x, int y, int z);
bool bIsReplica(int reqID);
void printGenericTree(GenericTreeNode* node, std::ostream& os) ;
//bool compBucket(GenericTreeNode *ln,GenericTreeNode *rn);

#ifdef REDUCTION_HELPER

class TreePieceCounter : public CkLocIterator { 
  public:
    TreePieceCounter() { reset(); }
    void addLocation(CkLocation &loc);
    void reset();

  public:
    CkVec<TreePiece *> presentTreePieces;
};



class ReductionHelper : public CBase_ReductionHelper {
  public:
  ReductionHelper();
  ReductionHelper(CkMigrateMessage *);
  void pup(PUP::er &p);

  void countTreePieces(const CkCallback &cb);
  void reduceBinCounts(int nBins, int64_t *binCounts, const CkCallback &cb);
  void evaluateBoundaries(SFC::Key *keys, const int n, int isRefine, const CkCallback& cb);
  void evaluateBoundaries(const CkBitVector &binsToSplit, const CkCallback& cb);

  private:
  void senseLocalTreePieces();

  private:

  CkVec<int64_t> myBinCounts;
  int numTreePiecesCheckedIn;

  TreePieceCounter localTreePieces;
  std::vector<SFC::Key> splitters;
};
#endif

/// @brief Used to count non-empty treepieces on the local processor.
class NonEmptyTreePieceCounter : public CkLocIterator {            
  public:
    NonEmptyTreePieceCounter() { reset(); }
    void addLocation(CkLocation &loc);
    void reset();

  public:
    int count;
};

#endif //PARALLELGRAVITY_H
