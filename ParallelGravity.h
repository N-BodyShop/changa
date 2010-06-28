/** @file ParallelGravity.h
 */

#ifndef PARALLELGRAVITY_H
#define PARALLELGRAVITY_H

#ifndef CELL
#ifdef CELL_NODE
#undef CELL_NODE
#endif
#ifdef CELL_PART
#undef CELL_PART
#endif
#ifdef CELL_EWALD
#undef CELL_EWALD
#endif
#endif

#include "config.h"
#include <string>
#include <map>
#include <vector>
#include <algorithm>

#include "pup_stl.h"
#include "comlib.h"

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
#include "liveViz.h"

#include "codes.h"
#include "CacheInterface.h"

#ifdef SPCUDA
#include "EwaldCUDA.h"
#endif

#ifdef CUDA
#include "cuda_typedef.h"
#endif

PUPbytes(InDumpFrame);
PUPbytes(COOL);
PUPbytes(COOLPARAM);

#ifdef HPM_COUNTER
#include <libhpm.h>
#endif

using namespace Tree;

enum DomainsDec {
  SFC_dec=0,
  Oct_dec=1,
  ORB_dec=2,
  SFC_peano_dec=3,
  SFC_peano_dec_3D=4,
  SFC_peano_dec_2D=5
};

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

#include "MultistepLB.decl.h"          // jetley - needed for CkIndex_MultistepLB

class SmoothParams;

#include "InOutput.h"

#include "ParallelGravity.decl.h"

extern CProxy_Main mainChare;
extern int verbosity;
extern bool _cache;
extern int _nocache;
extern int _cacheLineDepth;
extern unsigned int _yieldPeriod;
extern DomainsDec domainDecomposition;
extern GenericTrees useTree;
extern CProxy_TreePiece treeProxy;
extern CProxy_LvArray lvProxy;	    // Proxy for the liveViz array
extern CProxy_LvArray smoothProxy;  // Proxy for smooth reduction
extern CProxy_TreePiece streamingProxy;
extern CProxy_DataManager dMProxy;
extern unsigned int numTreePieces;
extern unsigned int particlesPerChare;

extern CProxy_CkCacheManager cacheGravPart;
extern CProxy_CkCacheManager cacheSmoothPart;
extern CProxy_CkCacheManager cacheNode;

extern ComlibInstanceHandle cinst1, cinst2;

/// The group ID of your DataManager.  You must set this!
extern CkGroupID dataManagerID;

extern int boundaryEvaluationUE;
extern int weightBalanceUE;
extern int networkProgressUE;
extern int nodeForceUE;
extern int partForceUE;

extern int _prefetch;
extern int _randChunks;
extern int _numChunks;
extern unsigned int bucketSize;
extern int lbcomm_cutoff_msgs;

//jetley
extern int localNodesPerReq;
extern int remoteNodesPerReq;
extern int remoteResumeNodesPerReq;
extern int localPartsPerReq;
extern int remotePartsPerReq;
extern int remoteResumePartsPerReq;

extern double largePhaseThreshold;

extern double theta;
extern double thetaMono;

extern int nSmooth;

class dummyMsg : public CMessage_dummyMsg{
public:
int val;
/*#if INTERLIST_VER > 0
int level;
GenericTreeNode* startNode;
#endif*/
};

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

/********************************************
class piecedata : public CMessage_piecedata {
public:
	int CellInteractions;
	int ParticleInteractions;
	int MACChecks;
	double totalmass;
	CkCallback cb;

	piecedata():CellInteractions(0),ParticleInteractions(0),MACChecks(0),totalmass(0.0) { }
	void modifypiecedata(int cell,int particle,int mac,double mass){
		CellInteractions += cell;
		ParticleInteractions += particle;
		MACChecks += mac;
		totalmass += mass;
	}
	void reset(){ CellInteractions=0;ParticleInteractions=0;MACChecks=0;totalmass=0.0; }
	void setcallback(CkCallback& cback) { cb = cback; }
	CkCallback& getcallback() { return cb; }
};
********************************************/

class ComputeChunkMsg : public CMessage_ComputeChunkMsg {
  ComputeChunkMsg() {} // not available
 public:
  int chunkNum;

  ComputeChunkMsg(int i) : chunkNum(i) {
  }
};


class ORBSplittersMsg : public CMessage_ORBSplittersMsg{
public:
	int length;
  double *pos;
  char *dim;
  CkCallback cb;

  ORBSplittersMsg(int len, CkCallback callback): length (len), cb(callback) {}

};

class Compare{ //Defines the comparison operator on the map used in balancer
  int dim;
public:
  Compare() {}
  Compare(int i) : dim(i) {}

  void setDim(int i){ dim = i; }

  bool operator()(GravityParticle& p1, GravityParticle& p2) const {
    return p1.position[dim] < p2.position[dim];
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

inline int RungToSubsteps(int iRung) {
  return 1 << (MAXRUNG - iRung);
}

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

inline double RungToDt(double dDelta, int iRung) {
  return dDelta*RungToSubsteps(iRung)*MAXSUBSTEPS_INV;
}

class Main : public CBase_Main {
	CkArgMsg *args;
	std::string basefilename;
	CProxy_Sorter sorter;
	int nTotalParticles;
	int nTotalSPH;
	//double theta; -- moved to readonly
	double dTime;		/* Simulation time */
	double dEcosmo;		/* variables for integrating
				   Lazer-Irvine eq. */
	double dUOld;
	double dTimeOld;
	unsigned int printBinaryAcc;
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
	int bChkFirst;		/* alternate between 0 and 1 for checkpoint */
	double dSimStartTime;   // Start time for entire simulation
	int iStop;		/* indicate we're stopping the
				   simulation early */
	int iPhase;
	int nPhases;		/* number of walks that a node cache
				   will service */
	int nActiveGrav;
	int nActiveSPH;

#ifdef CUDA
          double localNodesPerReqDouble;
          double remoteNodesPerReqDouble;
          double remoteResumeNodesPerReqDouble;
          double localPartsPerReqDouble;
          double remotePartsPerReqDouble;
          double remoteResumePartsPerReqDouble;
#endif

public:

	Main(CkArgMsg* m);
	Main(CkMigrateMessage *m);

	void niceExit();
	void setupICs();
	void initialForces();
	void doSimulation();
	void restart();
        void advanceBigStep(int);
	int adjust(int iKickRung);
	void rungStats();
	void countActive(int activeRung);
	void calcEnergy(double, double, char *) ;
	void getStartTime();
	void getOutTimes();
	int bOutTime();
	void writeOutput(int iStep) ;
	void setTotalParticles(int n) {
	    nTotalParticles = n;
        }
	void updateSoft();
	void growMass(double dTime, double dDelta);
	void initSph();
	void initCooling();
	void doSph(int activeRung);
	int DumpFrameInit(double dTime, double dStep, int bRestart);
	void DumpFrame(double dTime, double dStep);
	int nextMaxRungIncDF(int nextMaxRung);
	void pup(PUP::er& p);
	void liveVizImagePrep(liveVizRequestMsg *msg);
};

/* IBM brain damage */
#undef hz
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
typedef CkQ<OffsetNode> CheckList;
typedef CkVec<OffsetNode> UndecidedList;
typedef CkVec<UndecidedList> UndecidedLists;
#endif

typedef struct particlesInfoR{
    ExternalGravityParticle* particles;
    int numParticles;
    Vector3D<double> offset;
#if defined CHANGA_REFACTOR_PRINT_INTERACTIONS || defined CHANGA_REFACTOR_WALKCHECK_INTERLIST || defined CUDA
    NodeKey key;
#endif
#if COSMO_DEBUG > 1
    GenericTreeNode *nd;
#endif
} RemotePartInfo;

 ///Local Particle Info structure
typedef struct particlesInfoL{
    GravityParticle* particles;
    int numParticles;
    Vector3D<double> offset;
#if defined CHANGA_REFACTOR_PRINT_INTERACTIONS || defined CHANGA_REFACTOR_WALKCHECK_INTERLIST || defined CUDA
    NodeKey key;
#endif
#if COSMO_DEBUG > 1 || defined CUDA
    GenericTreeNode *nd;
#endif
} LocalPartInfo;


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

   TreeWalk *sTopDown;
   TreeWalk *twSmooth;
#if INTERLIST_VER > 0
   TreeWalk *sInterListWalk;
   // clearable, used for resumed walks
   State *sInterListStateRemoteResume;
#endif
   Compute *sGravity, *sPrefetch;
   Compute *sSmooth;
   
   Opt *sLocal, *sRemote, *sPref;
   Opt *optSmooth;

   State *sPrefetchState;
   State *sLocalGravityState, *sRemoteGravityState, *sSmoothState;
   typedef std::map<CkCacheKey, CkVec<int>* > SmPartRequestType;
   // buffer of requests for smoothParticles.
   SmPartRequestType smPartRequests;

   CkVec<ActiveWalk> activeWalks;
   int completedActiveWalks; // XXX this should be part of the gravity
			     // walk state.
   int nCacheAccesses; // keep track of outstanding cache accesses to
		       // know when writebacks complete.  XXX this
		       // should be part of the smooth state

 public:
#if COSMO_PRINT_BK > 1
  State *getSRemoteGravityState(){ return sRemoteGravityState; }
  State *getSLocalGravityState(){ return sLocalGravityState; }
#endif
   void addActiveWalk(int iAwi, TreeWalk *tw, Compute *c, Opt *o, State *s);

        void markWalkDone();
	void finishWalk();
        void markSmoothWalkDone();
	void finishSmoothWalk();

        int getIndex() {
          return thisIndex;
        }

        /*
        int decPrefetchWaiting() {
          prefetchWaiting--;
          return prefetchWaiting;
        }

        int incPrefetchWaiting() {
          prefetchWaiting++;
          return prefetchWaiting;
        }
        */

        /*
        int addToremaining Chunk(int chunk, int howMuch){
          remaining Chunk[chunk] += howMuch;
          return remaining Chunk[chunk];
        }
        */

        void addToNodeInterRemote(int chunk, int howmany){
          nodeInterRemote[chunk] += howmany;
        }

        void addToParticleInterRemote(int chunk, int howmany){
          particleInterRemote[chunk] += howmany;
        }

        void addToNodeInterLocal(int howmany){
          nodeInterLocal += howmany;
        }

        void addToParticleInterLocal(int howmany){
          particleInterLocal += howmany;
        }

        /// Start a new remote computation upon prefetch finished
        void startRemoteChunk();

        /*
        int getCurrentRemote Bucket(){
        	return currentRemote Bucket;
        }
        */

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
        BucketActiveInfo *bucketActiveInfo;


        int getNumParticles(){
        	return myNumParticles;
        }

        int getNumBuckets(){
        	return numBuckets;
        }

        void callFreeRemoteChunkMemory(int chunk);

        int getActiveRung(){ return activeRung; }
#ifdef CUDA_INSTRUMENT_WRS
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
          if(largePhase()){
            for(int b = 0; b < numBuckets; b++){
              GenericTreeNode *bucket = bucketList[b];
              int buckstart = bucket->firstParticle;
              int buckend = bucket->lastParticle;
              GravityParticle *buckparts = bucket->particlePointer;
              bucket->bucketArrayIndex = fillIndex;
              for(int i = buckstart; i <= buckend; i++){
                fillArray[fillIndex] = buckparts[i-buckstart];
#if defined CUDA_EMU_KERNEL_NODE_PRINTS || defined CUDA_EMU_KERNEL_PART_PRINTS
                fillArray[fillIndex].tp = thisIndex;
                fillArray[fillIndex].id = i;
#endif
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
#if defined CUDA_EMU_KERNEL_NODE_PRINTS || defined CUDA_EMU_KERNEL_PART_PRINTS
                  fillArray[fillIndex].tp = thisIndex;
                  fillArray[fillIndex].id = i;
#endif
                  fillIndex++;
                }
              }
              binfo->size = fillIndex-binfo->start;
            }
          }
        }

        bool isActive(int partNum){
          return myParticles[partNum].rung >= activeRung;
        }

        void clearMarkedBuckets(CkVec<GenericTreeNode *> &markedBuckets);
        void clearMarkedBucketsAll();

#ifdef CUDA_STATS
        long long localNodeInteractions;
        long long localPartInteractions;
        long long remoteNodeInteractions;
        long long remotePartInteractions;
        long long remoteResumeNodeInteractions;
        long long remoteResumePartInteractions;
#endif

#ifdef CUDA_INSTRUMENT_WRS
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
	/// convert a key to a node using the nodeLookupTable
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

        // jetley - proxy for load balancer
        CkGroupID proxy;
        // jetley - whether proxy is valid or not
        CmiBool proxyValid;
        // jetley - saved first internal node
        Vector3D<float> savedCentroid;
        CmiBool proxySet;
        // jetley - multistep load balancing
        int prevLARung;

	/// @brief Used to inform the mainchare that the requested operation has
	/// globally finished
	CkCallback callback;
	/// smooth globally finished
	CkCallback cbSmooth;
	/// Total number of particles contained in this chare
	unsigned int myNumParticles;
	/// Array with the particles in this chare
	GravityParticle* myParticles;
 public:
	/// Total Particles in the simulation
	int64_t nTotalParticles;
	/// Total Gas Particles
	int64_t nTotalSPH;
 private:
	/// Total gas in this chare
	unsigned int myNumSPH;
	/// Array with SPH particle data
	extraSPHData *mySPHParticles;
	/// Total Star Particles
	int64_t nTotalStars;

	/// List of all the node-buckets in this TreePiece
	std::vector<GenericTreeNode *> bucketList;

	/// Array with sorted particles for domain decomposition (ORB)
	std::vector<GravityParticle> mySortedParticles;
        /// Array with incoming particles for domain decomposition (without stl)
        GravityParticle *incomingParticles;
        /// How many particles have already arrived during domain decomposition
        int incomingParticlesArrived;
        /// Flag to acknowledge that the current TreePiece has already
        /// contributed to itself, prevents that acceptSortedParticles will
        /// change the TreePiece particles before the one belonging to someone
        /// else have been sent out
        bool incomingParticlesSelf;
	std::vector<extraSPHData> *incomingGas;

	/// holds the total mass of the current TreePiece
	double piecemass;

	/// @if ALL

	/// used to determine which coordinate we are outputting, while printing
	/// accelerations in ASCII format
	int cnt;
	/// used to determine if the x, y, z coordinates should be printed
	/// together while outputting accelerations in ASCII format
	int packed;

	/// @endif

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
	CkVec<int> myBinCounts;
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

	typedef std::map<NodeKey, CkVec<int>* >   MomentRequestType;
	/// Keep track of the requests for remote moments not yet satisfied.
	/// Used only during the tree construction.
	MomentRequestType momentRequests;

	/// Opening angle
	//double theta; -- moved as readonly
	/// Opening angle - monopole
	//double thetaMono; -- moved as readonly
        /// The current active mask for force computation in multistepping
        int activeRung;

	/// Periodic Boundary stuff
	int bPeriodic;
	Vector3D<double> fPeriod;
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

	/// Setup for writing
	int nSetupWriteStage;
	int64_t nStartWrite;	// Particle number at which this piece starts
				// to write file.

	/// Map between Keys and TreeNodes, used to get a node from a key
	NodeLookupType nodeLookupTable;

	/// Number of particles which are still traversing the tree
	//u_int64_t myNumParticlesPending;

        /// Number of pending chenks
        // A chunk is pending wrt a TP until that TP has finished using
        // it completely (i.e. state->counterArrays[1][chunk] == 0)
        //int numPendingChunks;

	/// Number of nodes still missing before starting the real computation
	//u_int64_t prefetchWaiting;
	/// Current prefetching chunk in progress
        //int currentPrefetch;
	/// Array of keys that will be the root of the prefetching chunks
	Tree::NodeKey *prefetchRoots;
	/// Placeholder for particles used for prefetching
	OrientedBox<double> *prefetchReq;
	unsigned int numPrefetchReq;
	/// number of particles/buckets still remaining to compute for the chunk
	//int *remaining Chunk;

	/// number of chunks in which the tree will be chopped for prefetching
	int numChunks;

	//u_int64_t openingDiffCount;
    /// @if STATISTICS

#if COSMO_STATS > 0
	//u_int64_t myNumCellInteractions;
	//u_int64_t myNumParticleInteractions;
	u_int64_t myNumMACChecks;
	//u_int64_t myNumProxyCalls;
	//u_int64_t myNumProxyCallsBack;
	// Same as myNumCellInteractions, only restricted to cached nodes
	//int cachecellcount;
	u_int64_t nodesOpenedLocal;
	u_int64_t nodesOpenedRemote;
	u_int64_t numOpenCriterionCalls;
#endif
	u_int64_t nodeInterLocal;
	u_int64_t *nodeInterRemote;
	u_int64_t particleInterLocal;
	u_int64_t *particleInterRemote;
	int nActive;		// number of particles that are active

	/// @endif

	/// Size of bucketList, total number of buckets present
	unsigned int numBuckets;
	/// Used to start the remote computation for a particular chunk for all
	/// buckets, one after the other
	//unsigned int currentRemote Bucket;
#if INTERLIST_VER > 0
	int prevRemoteBucket;
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

  ///Level after which the local subtrees of all the TreePieces start
  unsigned int chunkRootLevel;

  ///Keeps track of the bounding boxes and dimensions along which they are splitted, starting from
  ///root to chunkRootLevel
  OrientedBox<float>* boxes;
  char* splitDims;

  ///Phase of ORB decomposition: the number of boxes double in each phase till they are equal to the number
  ///of TreePieces
  int phase;

 #if INTERLIST_VER > 0


 public:
#ifdef CELL
	friend void cellSPE_callback(void*);
	friend void cellSPE_ewald(void*);
#endif
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
#endif
  void EwaldGPU(); 
  void EwaldGPUComplete();
  int bLoaded;		/* Are particles loaded? */

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

  GenericTreeNode *get3DIndex();

	/// Recursive call to build the subtree with root "node", level
	/// specifies the level at which "node" resides inside the tree
	void buildOctTree(GenericTreeNode* node, int level);
#ifdef TREE_BREADTH_FIRST
	void growBottomUp(GenericTreeNode* node);
#endif

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
	/** @brief Initial walk through the tree. It will continue until local
	 * nodes are found (excluding those coming from the cache). When the
	 * treewalk is finished it stops and cachedWalkBucketTree will continue
	 * with the incoming nodes.
	 */
	void walkBucketTree(GenericTreeNode* node, int reqID);
	/** @brief Start the treewalk for the next bucket among those belonging
	 * to me. The buckets are simply ordered in a vector.
	 */
	void startNextBucket();
	/** @brief Start a full step of bucket computation, it sends a message
	 * to trigger nextBucket which will loop over all the buckets.
	 */
	void doAllBuckets();
	void reconstructNodeLookup(GenericTreeNode *node);
	//void rebuildSFCTree(GenericTreeNode *node,GenericTreeNode *parent,int *);

public:
 TreePiece() : pieces(thisArrayID), root(0), proxyValid(false),
	    proxySet(false), prevLARung (-1), sTopDown(0), sGravity(0),
	    sPrefetch(0), sLocal(0), sRemote(0), sPref(0), sSmooth(0) {
	  //CkPrintf("[%d] TreePiece created on proc %d\n",thisIndex, CkMyPe());
	  // ComlibDelegateProxy(&streamingProxy);
	  dm = NULL;
	  iterationNo=0;
	  usesAtSync=CmiTrue;
	  bucketReqs=NULL;
	  nCacheAccesses = 0;
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
#ifdef CUDA_STATS
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
	  numPrefetchReq = 0;
	  prefetchReq = NULL;
	  //remaining Chunk = NULL;
	  ewt = NULL;
	  nMaxEwhLoop = 100;

          incomingParticles = NULL;
          incomingParticlesArrived = 0;
          incomingParticlesSelf = false;
    orbBoundaries.clear();
	}

	TreePiece(CkMigrateMessage* m) {
          // jetley
          proxyValid = false;
          proxySet = false;

	  usesAtSync = CmiTrue;
	  //localCache = NULL;
	  dm = NULL;
	  bucketReqs = NULL;
	  numChunks=-1;
	  nCacheAccesses = 0;
	  completedActiveWalks = 0;
	  prefetchRoots = NULL;
	  //remaining Chunk = NULL;
          ewt = NULL;

      sTopDown = 0;
	  sGravity = NULL;
	  sPrefetch = NULL;
	  sSmooth = NULL;
#if INTERLIST_VER > 0
	  sInterListWalk = NULL;
#endif

          incomingParticles = NULL;
          incomingParticlesArrived = 0;
          incomingParticlesSelf = false;

	  nodeInterRemote = NULL;
          particleInterRemote = NULL;

	  orbBoundaries.clear();
	}

        private:
        void freeWalkObjects();

        public:
	~TreePiece() {
	  if (verbosity>1) ckout <<"Deallocating treepiece "<<thisIndex<<endl;
	  delete[] myParticles;
	  delete[] mySPHParticles;
	  //delete[] splitters;
	  //delete[] prefetchRoots;
	  //delete[] remaining Chunk;
	  delete[] nodeInterRemote;
	  delete[] particleInterRemote;
	  delete[] bucketReqs;
	  delete[] prefetchReq;
          delete[] ewt;

	  // recursively delete the entire tree
	  if (root != NULL) {
	    root->fullyDelete();
	    delete root;
	  }

          if (verbosity>1) ckout <<"Finished deallocation of treepiece "<<thisIndex<<endl;
	}

	void restart();

	void setPeriodic(int nReplicas, Vector3D<double> fPeriod, int bEwald,
			 double fEwCut, double fEwhCut, int bPeriod);
	void BucketEwald(GenericTreeNode *req, int nReps,double fEwCut);
	void EwaldInit();
	void calculateEwald(dummyMsg *m);
	// Scale velocities (needed to convert to canonical momenta for
	// comoving coordinates.)
	void velScale(double dScale);

	// Parse NChilada description file
	int parseNC(const std::string& fn);
	// Load from mass and position files
	void load(const std::string& fn, const CkCallback& cb);

	// Load from Tipsy file
	void loadTipsy(const std::string& filename, const double dTuFac,
		       const CkCallback& cb);
	// Write a Tipsy file
	void writeTipsy(const std::string& filename, const double dTime,
			const double dvFac, const double duTfac,
			const int bCool);
	// Find position in the file to start writing
	void setupWrite(int iStage, u_int64_t iPrevOffset,
			const std::string& filename, const double dTime,
			const double dvFac, const double duTFac,
			const int bCool, const CkCallback& cb);
	// serial output
	void serialWrite(u_int64_t iPrevOffset, const std::string& filename,
			 const double dTime,
			 const double dvFac, const double duTFac,
			 const int bCool, const CkCallback& cb);
	// setup for serial output
	void oneNodeWrite(int iIndex,
			       int iOutParticles,
			       int iOutSPH,
			       GravityParticle *particles, // particles to
						     // write
			       extraSPHData *pGas, // SPH data
			       int *piSPH, // SPH data offsets
			       const u_int64_t iPrevOffset,
			       const std::string& filename,  // output file
			       const double dTime,      // time or expansion
			       const double dvFac,  // velocity conversion
			     const double duTFac, // temperature
						  // conversion
			  const int bCool,
			  const CkCallback &cb);
	// Reorder for output
	void reOrder(CkCallback& cb);
	// move particles around for output
	void ioAcceptSortedParticles(const GravityParticle* particles,
				     const int n, const extraSPHData *pGas,
				     const int nGasIn);
	/** Inform the DataManager of my node that I'm here.
	 The callback will receive a CkReductionMsg containing no data.
	void registerWithDataManager(const CkGroupID& dataManagerID,
				     const CkCallback& cb);
	 */
	// Assign keys after loading tipsy file and finding Bounding box
	void assignKeys(CkReductionMsg* m);
	void evaluateBoundaries(SFC::Key* keys, const int n, int isRefine, const CkCallback& cb);
	void unshuffleParticles(CkReductionMsg* m);
	void acceptSortedParticles(const GravityParticle* particles,
				   const int n, const extraSPHData *pExtra,
				   const int nGas);
  /*****ORB Decomposition*******/
  void initORBPieces(const CkCallback& cb);
  void initBeforeORBSend(unsigned int myCount, const CkCallback& cb, const CkCallback& cback);
  void sendORBParticles();
  void acceptORBParticles(const GravityParticle* particles, const int n);
  void finalizeBoundaries(ORBSplittersMsg *splittersMsg);
  void evaluateParticleCounts(ORBSplittersMsg *splittersMsg);
  /*****************************/

  void kick(int iKickRung, double dDelta[MAXRUNG+1], int bClosing,
	    int bNeedVPred, int bGasIsothermal, double duDelta[MAXRUNG+1],
	    const CkCallback& cb);
  void drift(double dDelta, int bNeedVPred, int bGasIsothermal, double dvDelta,
	     double duDelta, int nGrowMass, const CkCallback& cb);
  void initAccel(int iKickRung, const CkCallback& cb);
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
 * @param dDelta Base timestep
 * @param dAccFac Acceleration scaling for cosmology
 * @param dCosmoFac Cosmo scaling for Courant
 * @param cb Callback function reduces currrent maximum rung
 */
  void adjust(int iKickRung, int bEpsAccStep, int bGravStep,
	      int bSphStep, int bViscosityLimitdt,
		       double dEta, double dEtaCourant, double dEtauDot,
		       double dDelta, double dAccFac,
		       double dCosmoFac, const CkCallback& cb);
  void rungStats(const CkCallback& cb);
  void countActive(int activeRung, const CkCallback& cb);
	void calcEnergy(const CkCallback& cb);
	void setSoft(const double dSoft);
	void physicalSoft(const double dSoftMax, const double dFac,
			  const int bSoftMaxMul, const CkCallback& cb);
	void growMass(int nGrowMass, double dDeltaM, const CkCallback& cb);
	// void initCooling(COOL inCool, COOLPARAM inParam, const CkCallback& cb);
	void InitEnergy(double dTuFac, double z, double dTime,
			const CkCallback& cb);
	void updateuDot(int activeRung, double duDelta[MAXRUNG+1],
			double dTime, double z, int bCool,
			int bUpdateState, const CkCallback& cb);
	void ballMax(int activeRung, double dFac, const CkCallback& cb);
	void sphViscosityLimiter(int bOn, int activeRung, const CkCallback& cb);
	void getAdiabaticGasPressure(double gamma, double gammam1,
				     const CkCallback &cb);
	void getCoolingGasPressure(double gamma, double gammam1,
				   const CkCallback &cb);
	void SetTypeFromFileSweep(int iSetMask, char *file,
	   struct SortStruct *ss, int nss, int *pniOrder, int *pnSet);
	void setTypeFromFile(int iSetMask, char *file, const CkCallback& cb);
	void getCOM(const CkCallback& cb, int bLiveViz);
	void getCOMByType(int iType, const CkCallback& cb, int bLiveViz);
	void DumpFrame(InDumpFrame in, const CkCallback& cb, int liveVizDump) ;
	void liveVizDumpFrameInit(liveVizRequestMsg *msg);

	/// Charm entry point to build the tree (called by Main), calls collectSplitters
	void buildTree(int bucketSize, const CkCallback& cb);

	/// Collect the boundaries of all TreePieces, and trigger the real treebuild
	//void collectSplitters(CkReductionMsg* m);
	/// Real tree build, independent of other TreePieces; calls the recursive buildTree
	void startOctTreeBuild(CkReductionMsg* m);

  /********ORB Tree**********/
  //void receiveBoundingBoxes(BoundingBoxes *msg);
  void startORBTreeBuild(CkReductionMsg* m);
  OrientedBox<float> constructBoundingBox(GenericTreeNode *node,int level, int numChild);
  void buildORBTree(GenericTreeNode * node, int level);
  /**************************/

	/// Request the TreePiece to send back later the moments for this node.
	void requestRemoteMoments(const Tree::NodeKey key, int sender);
	void receiveRemoteMoments(const Tree::NodeKey key, Tree::NodeType type, int firstParticle, int numParticles, const MultipoleMoments& moments, const OrientedBox<double>& box, const OrientedBox<double>& boxBall);

	/// Decide whether the node should be opened for the force computation
	/// of the given request. --- Moved outside TreePiece class
	/*bool openCriterionBucket(GenericTreeNode *node,
				 GenericTreeNode *bucketNode,
				 Vector3D<double> offset // Offset of node
				 );

	int openCriterionNode(GenericTreeNode *node,
			      GenericTreeNode *myNode,
			      Vector3D<double> offset // Offset of node
			      );
	 */

	/// Entry point for the local computation: for each bucket compute the
	/// force that its particles see due to the other particles hosted in
	/// this TreePiece. The opening angle theta has already been passed
	/// through "startIteration"
	void calculateGravityLocal();
	void commenceCalculateGravityLocal();

	/// Entry point for the remote computation: for each bucket compute the
	/// force that its particles see due to the other particles NOT hosted
	/// by this TreePiece, and belonging to a subset of the global tree
	/// (specified by chunkNum).
	void calculateGravityRemote(ComputeChunkMsg *msg);

	/// Temporary function to recurse over all the buckets like in
	/// walkBucketTree, only that NonLocal nodes are the only one for which
	/// forces are computed
	void walkBucketRemoteTree(GenericTreeNode *node, int chunk, int reqID, bool isRoot);

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

  /// Start a tree based gravity computation.
  /// @param am the active rung for the computation
  /// @param theta the opening angle
  /// @param cb the callback to use after all the computation has finished
  void startIteration(int am, double myTheta, const CkCallback& cb);
  /// As above, but for a smooth operation.
  void setupSmooth();
  void startIterationSmooth(SmoothParams *p, const CkCallback& cb);
  void startIterationReSmooth(SmoothParams *p, const CkCallback& cb);
  void startIterationMarkSmooth(SmoothParams *p, const CkCallback& cb);

  void finishNodeCache(int iPhases, const CkCallback& cb);
	/// Function called by the CacheManager to send out request for needed
	/// remote data, so that the later computation will hit.
	void prefetch(GenericTreeNode *node, int offsetID);
	void prefetch(ExternalGravityParticle *part);

	/// @brief Retrieve the remote node, goes through the cache if present
        //GenericTreeNode* requestNode(int remoteIndex, Tree::NodeKey lookupKey, int chunk, int reqID, bool isPrefetch=false);

        GenericTreeNode* requestNode(int remoteIndex, Tree::NodeKey lookupKey, int chunk, int reqID, int awi, void *source, bool isPrefetch);
	/// @brief Receive a request for Nodes from a remote processor, copy the
	/// data into it, and send back a message.
	void fillRequestNode(CkCacheRequestMsg *msg);
	/** @brief Receive the node from the cache as following a previous
	 * request which returned NULL, and continue the treewalk of the bucket
	 * which requested it with this new node.
	*/
#if 0
	void receiveNode(GenericTreeNode &node, int chunk, unsigned int reqID);
	/// Just and inline version of receiveNode
	void receiveNode_inline(GenericTreeNode &node, int chunk, unsigned int reqID);
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

        ExternalGravityParticle *requestParticles(Tree::NodeKey key,int chunk,int remoteIndex,int begin,int end,int reqID, int awi, void *source, bool isPrefetch=false);
	GravityParticle *requestSmoothParticles(Tree::NodeKey key, int chunk,
				    int remoteIndex, int begin,int end,
				    int reqID, int awi, void *source, bool isPrefetch);
	void fillRequestParticles(CkCacheRequestMsg *msg);
	void fillRequestSmoothParticles(CkCacheRequestMsg *msg);
	void flushSmoothParticles(CkCacheFillMsg *msg);
	void processReqSmoothParticles();

	//void startlb(CkCallback &cb);
	void startlb(CkCallback &cb, int activeRung);
	void ResumeFromSync();

	void outputAccelerations(OrientedBox<double> accelerationBox, const std::string& suffix, const CkCallback& cb);
	void outputASCII(OutputParams& params, int bParaWrite,
			 const CkCallback& cb);
	void oneNodeOutVec(OutputParams& params, Vector3D<double>* avOut,
			   int nPart, int iIndex, int bDone,
			   CkCallback& cb) ;
	void oneNodeOutArr(OutputParams& params, double* adOut,
			   int nPart, int iIndex, int bDone,
			   CkCallback& cb) ;
	
	void outputIOrderASCII(const std::string& suffix, const CkCallback& cb);
	void outputStatistics(const CkCallback& cb);

	/// Collect the total statistics from the various chares
	void collectStatistics(CkCallback &cb);
	//void getPieceValues(piecedata *totaldata);

        /** @brief Entry method used to split the processing of all the buckets
         * in small pieces. It call startNextBucket a _yieldPeriod number of
         * times, and then it returns to the scheduler after enqueuing a message
         * for itself.
	 */
        void nextBucket(dummyMsg *m);

	void report();
	void printTreeViz(GenericTreeNode* node, std::ostream& os);
	void printTree(GenericTreeNode* node, std::ostream& os);
	void pup(PUP::er& p);

        // jetley
        // need this in TreeWalk
        GenericTreeNode *getRoot() {return root;}
        // need this in Compute
	Vector3D<double> decodeOffset(int reqID);

        GenericTreeNode *nodeMissed(int reqID, int remoteIndex, Tree::NodeKey &key, int chunk, bool isPrefetch, int awi, void *source);

        ExternalGravityParticle *particlesMissed(Tree::NodeKey &key, int chunk, int remoteIndex, int firstParticle, int lastParticle, int reqID, bool isPrefetch, int awi, void *source);

        void receiveNodeCallback(GenericTreeNode *node, int chunk, int reqID, int awi, void *source);
        void receiveParticlesCallback(ExternalGravityParticle *egp, int num, int chunk, int reqID, Tree::NodeKey &remoteBucket, int awi, void *source);
        void receiveParticlesFullCallback(GravityParticle *egp, int num, int chunk, int reqID, Tree::NodeKey &remoteBucket, int awi, void *source);
        void receiveProxy(CkGroupID _proxy){ proxy = _proxy; proxySet = true; /*CkPrintf("[%d : %d] received proxy\n", CkMyPe(), thisIndex);*/}
        void doAtSync();
        GravityParticle *getParticles(){return myParticles;}

};

class LvArray : public CBase_LvArray {
 public:
    LvArray() {}
    LvArray(CkMigrateMessage* m) {}
    } ;

int decodeReqID(int);
int encodeOffset(int reqID, int x, int y, int z);
bool bIsReplica(int reqID);
void initNodeLock();
void printGenericTree(GenericTreeNode* node, std::ostream& os) ;
//bool compBucket(GenericTreeNode *ln,GenericTreeNode *rn);
#endif //PARALLELGRAVITY_H
