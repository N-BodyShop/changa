/** @file ParallelGravity.h
 */
 
#ifndef PARALLELGRAVITY_H
#define PARALLELGRAVITY_H

#include <string>
#include <map>
#include <vector>
#include <algorithm>

#include "pup_stl.h"
#include "comlib.h"

#include "Vector3D.h"
#include "tree_xdr.h"
#include "SFC.h"
#include "TreeNode.h"
#include "GenericTreeNode.h"
#include "Interval.h"

using namespace Tree;

enum DomainsDec {
  SFC_dec,
  Oct_dec,
  ORB_dec
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

#include "CacheManager.h"
#include "ParallelGravity.decl.h"

extern CProxy_Main mainChare;
extern int verbosity;
extern bool _cache;
extern int _cacheLineDepth;
extern unsigned int _yieldPeriod;
extern DomainsDec domainDecomposition;
extern GenericTrees useTree;
extern CProxy_TreePiece streamingProxy;

/// The group ID of your DataManager.  You must set this!
extern CkGroupID dataManagerID;


extern bool _prefetch;
extern int _numChunks;

class dummyMsg : public CMessage_dummyMsg{
public:
int val;
};

class TreePieceStatistics {
  u_int64_t nodesOpenedLocal;
  u_int64_t nodesOpenedRemote;
  u_int64_t nodeInterLocal;
  u_int64_t nodeInterRemote;
  u_int64_t particleInterLocal;
  u_int64_t particleInterRemote;

  TreePieceStatistics() : nodesOpenedLocal(0), nodesOpenedRemote(0), nodeInterLocal(0),
    nodeInterRemote(0), particleInterLocal(0), particleInterRemote(0) { }

 public:
  TreePieceStatistics(u_int64_t nol, u_int64_t nor, u_int64_t nil, u_int64_t nir,
		      u_int64_t pil, u_int64_t pir) :
    nodesOpenedLocal(nol), nodesOpenedRemote(nor), nodeInterLocal(nil),
    nodeInterRemote(nir), particleInterLocal(pil), particleInterRemote(pir) { }

  void printTo(CkOStream &os) {
    os << "  TreePiece: " << nodesOpenedLocal << " local nodes opened, ";
    os << nodesOpenedRemote << " remote" << endl;
    os << "  TreePiece: " << nodeInterLocal << " local particle-node interactions, ";
    os << nodeInterRemote << " remote" << endl;
    os << "  TreePiece: " << particleInterLocal << " local particle-particle interactions, ";
    os << particleInterRemote << " remote" << endl;
  }

  static CkReduction::reducerType sum;

  static CkReductionMsg *sumFn(int nMsg, CkReductionMsg **msgs) {
    TreePieceStatistics ret;
    for (int i=0; i<nMsg; ++i) {
      CkAssert(msgs[i]->getSize() == sizeof(TreePieceStatistics));
      TreePieceStatistics *data = (TreePieceStatistics *)msgs[i]->getData();
      ret.nodesOpenedLocal += data->nodesOpenedLocal;
      ret.nodesOpenedRemote += data->nodesOpenedRemote;
      ret.nodeInterLocal += data->nodeInterLocal;
      ret.nodeInterRemote += data->nodeInterRemote;
      ret.particleInterLocal += data->particleInterLocal;
      ret.particleInterRemote += data->particleInterRemote;
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

  ComputeChunkMsg(int i) : chunkNum(i) { }
};

class RequestNodeMsg : public CMessage_RequestNodeMsg {
 public:
  int retIndex;
  int depth;
  Tree::NodeKey key;
  unsigned int reqID;

  RequestNodeMsg(int r, int d, Tree::NodeKey k, unsigned int req) : retIndex(r), depth(d), key(k), reqID(req) {}
};

class RequestParticleMsg : public CMessage_RequestParticleMsg {
 public:
  int retIndex;
  int begin;
  int end;
  Tree::NodeKey key;
  unsigned int reqID;

  RequestParticleMsg(int r, int b, int e, Tree::NodeKey k, unsigned int req) : retIndex(r), begin(b), end(e), key(k), reqID(req) {}
};

class FillNodeMsg : public CMessage_FillNodeMsg {
 public:
  int owner;
  //int count;
  char *nodes;

  FillNodeMsg(int index) : owner(index) { }
};

class Main : public Chare {
	std::string basefilename;
	CProxy_TreePiece pieces;
	CProxy_DataManager dataManager;
	CProxy_Sorter sorter;
	unsigned int numTreePieces;
	double theta;	
	double dTimeStep;
	unsigned int bucketSize;
	unsigned int printBinaryAcc;
public:
		
	Main(CkArgMsg* m);

	void niceExit();
	/**
	 * Principal method which does all the coordination of the work.
	 * It is threaded, so it suspends while waiting for the TreePieces to
	 * perform the task. It calls TreePieces for: 1) load, 2) buildTree, 3)
	 * calculateGravity and startlb for the desired number of iterations.
	 */
	void nextStage();
	void calcEnergy() ;
};

class TreePiece : public CBase_TreePiece {
 private:
	unsigned int numTreePieces;
	/// @brief Used to inform the mainchare that the requested operation has
	/// globally finished
	CkCallback callback;
	/// Total number of particles contained in this chare
	unsigned int myNumParticles;
	/// Array with the particles in this chare
	GravityParticle* myParticles;

	/// Array with sorted particles for domain decomposition
	vector<GravityParticle> mySortedParticles;

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

	unsigned int numSplitters;
	SFC::Key* splitters;
	CProxy_TreePiece pieces;
	/// A proxy to the DataManager.
	CProxy_DataManager dataManager;
	/// A local pointer to my DataManager.
	DataManager* dm;

	/// The counts of how many particles belonging to other
	/// TreePieces I currently hold
	std::vector<int> myBinCounts;
	/// My index in the responsibility array.
	int myPlace;
	/// The keys determined by the Sorter that separate me from my
	/// neighbors.
	SFC::Key leftSplitter, rightSplitter;

	std::string basefilename;
	OrientedBox<float> boundingBox;
	FieldHeader fh;
	bool started;
	unsigned iterationNo;
	/// The root of the global tree, always local to any chare
	GenericTreeNode* root;

	typedef std::map<NodeKey, CkVec<int>* >   MomentRequestType;
	/// Keep track of the requests for remote moments not yet satisfied.
	/// Used only during the tree construction.
	MomentRequestType momentRequests;

	/// Opening angle
	double theta;

	/// Map between Keys and TreeNodes, used to get a node from a key
	NodeLookupType nodeLookupTable;

	/// Number of particles which are still traversing the tree
	u_int64_t myNumParticlesPending;

	/// Number of nodes still missing before starting the real computation
	u_int64_t prefetchWaiting;
	/// Current prefetching chunk in progress
	int currentPrefetch;
	/// Array of keys that will be the root of the prefetching chunks
	Tree::NodeKey *prefetchRoots;
	/// Placeholder for particles used for prefetching
	BucketGravityRequest prefetchReq[2];
	/// number of particles/buckets still remaining to compute for the chunk
	int *remainingChunk;

	/// Start a new remote computation upon prefetch finished
	void startRemoteChunk();

	/// number of chunks in which the tree will be chopped for prefetching
	int numChunks;

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
	u_int64_t nodeInterLocal;
	u_int64_t nodeInterRemote;
	u_int64_t particleInterLocal;
	u_int64_t particleInterRemote;
#endif

	/// @endif

	/// Size of bucketList, total number of buckets present
	unsigned int numBuckets;
	/// Used to start the computation for all buckets, one after the other
	unsigned int currentBucket;
	/// Used to start the remote computation for a particular chunk for all
	/// buckets, one after the other
	unsigned int currentRemoteBucket;
	/// List of all the node-buckets in this TreePiece
	std::vector<GenericTreeNode *> bucketList;
	/// @brief Used as a placeholder while traversing the tree and computing
	/// forces. This should not exist in the cache version since all its
	/// information is duplicate from the correspondent TreeNode (in
	/// bucketList), and myParticles.
	/// @todo Eliminate the usage of this in the cache walk
	BucketGravityRequest *bucketReqs;
	
	/// Pointer to the instance of the local cache
	CacheManager *localCache;

	/// convert a key to a node using the nodeLookupTable
	inline GenericTreeNode *keyToNode(const Tree::NodeKey);

#if COSMO_DEBUG > 1
  ///A set for each bucket to check the correctness of the treewalk
  typedef std::vector< std::set<Tree::NodeKey> > DebugList;
  DebugList bucketcheckList;
#endif
  
 public:

	int bLoaded;		/* Are particles loaded? */
	
#if COSMO_DEBUG > 1
  ///This function checks the correctness of the treewalk
  void checkWalkCorrectness();
  void combineKeys(Tree::NodeKey key,int bucket);
#endif

  /* DEBUGGING */
  void quiescence();
  /*END DEBUGGING */


	/// Recursive call to build the subtree with root "node", level
	/// specifies the level at which "node" resides inside the tree
	void buildOctTree(GenericTreeNode* node, int level);

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
	/** @brief Initial walk through the tree. It will continue until local
	 * nodes are found (excluding those coming from the cache). When the
	 * treewalk is finished it stops and cachedWalkBucketTree will continue
	 * with the incoming nodes.
	 */
	void walkBucketTree(GenericTreeNode* node, BucketGravityRequest& req);
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
	
	TreePiece(unsigned int numPieces) : numTreePieces(numPieces), pieces(thisArrayID), started(false), root(0) {
	  //CkPrintf("[%d] TreePiece created\n",thisIndex);
	  // ComlibDelegateProxy(&streamingProxy);
	  // the lookup of localCache is done in startOctTreeBuild, when we also markPresence
	  /*if(_cache){	
	    localCache = cacheManagerProxy.ckLocalBranch();
	    }*/
	  localCache = NULL;
	  iterationNo=0;
	  usesAtSync=CmiTrue;
	  bucketReqs=NULL;
	  myPlace = -1;
#if COSMO_STATS > 0
	  nodesOpenedLocal = 0;
	  nodesOpenedRemote = 0;
	  nodeInterLocal = 0;
	  nodeInterRemote = 0;
	  particleInterLocal = 0;
	  particleInterRemote = 0;

	  piecemass = 0.0;
	  packed=0;
	  cnt=0;
#endif

	  // temporarely set to -1, it will updated after the tree is built
	  numChunks=-1;
	  prefetchRoots = NULL;
	}
	
	TreePiece(CkMigrateMessage* m) {
	  usesAtSync = CmiTrue;
	  localCache = NULL;
	  bucketReqs = NULL;
	  prefetchRoots = NULL;
	}

	~TreePiece() {
	  delete[] myParticles;
	  delete[] splitters;
	  delete[] prefetchRoots;
	  delete[] remainingChunk;
	  delete[] bucketReqs;
	  // recursively delete the entire tree
	  root->fullyDelete();
	  delete root;
	}
	
	// Load from mass and position files
	void load(const std::string& fn, const CkCallback& cb);

	// Load from Tispy file
	void loadTipsy(const std::string& filename, const CkCallback& cb);
	/** Inform the DataManager of my node that I'm here.
	 The callback will receive a CkReductionMsg containing no data.
	 */
	void registerWithDataManager(const CkGroupID& dataManagerID,
				     const CkCallback& cb);
	// Assign keys after loading tipsy file and finding Bounding box
	void assignKeys(CkReductionMsg* m);
	void evaluateBoundaries(const CkCallback& cb);
	void unshuffleParticles(CkReductionMsg* m);
	void acceptSortedParticles(const GravityParticle* particles,
				   const int n);
	void kick(double dDelta, const CkCallback& cb);
	void drift(double dDelta, const CkCallback& cb);
	void calcEnergy(const CkCallback& cb);
	/// Charm entry point to build the tree (called by Main), calls collectSplitters
	void buildTree(int bucketSize, const CkCallback& cb);

	/// Collect the boundaries of all TreePieces, and trigger the real treebuild
	void collectSplitters(CkReductionMsg* m);
	/// Real tree build, independent of other TreePieces; calls the recursive buildTree
	void startOctTreeBuild(CkReductionMsg* m);

	/// Request the TreePiece to send back later the moments for this node.
	void requestRemoteMoments(const Tree::NodeKey key, int sender);
	void receiveRemoteMoments(const Tree::NodeKey key, Tree::NodeType type, int firstParticle, int numParticles, const MultipoleMoments& moments);

	/// Decide whether the node should be opened for the force computation
	/// of the given request.
	bool openCriterionBucket(GenericTreeNode *node,
				 BucketGravityRequest& req);

	/// Entry point for the local computation: for each bucket compute the
	/// force that its particles see due to the other particles hosted in
	/// this TreePiece. The opening angle theta has already been passed
	/// through "startIteration"
	void calculateGravityLocal();

	/// Entry point for the remote computation: for each bucket compute the
	/// force that its particles see due to the other particles NOT hosted
	/// by this TreePiece, and belonging to a subset of the global tree
	/// (specified by chunkNum).
	void calculateGravityRemote(ComputeChunkMsg *msg);

	/// Temporary function to recurse over all the buckets like in
	/// walkBucketTree, only that NonLocal nodes are the only one for which
	/// forces are computed
	void walkBucketRemoteTree(GenericTreeNode *node, int chunk, BucketGravityRequest &req, bool isRoot);

	/// Function called by the CacheManager to start a new iteration.
	/// @param t the opening angle
	/// @param cb the callback to use after all the computation has finished
	void startIteration(double t, const CkCallback& cb);

	/// Function called by the CacheManager to send out request for needed
	/// remote data, so that the later computation will hit.
	void prefetch(GenericTreeNode *node);
	void prefetch(GravityParticle *part);

	/// @brief Retrieve the remote node, goes through the cache if present
	GenericTreeNode* requestNode(int remoteIndex, Tree::NodeKey lookupKey, int chunk,
				 BucketGravityRequest& req, bool isPrefetch=false);
	/// @brief Receive a request for Nodes from a remote processor, copy the
	/// data into it, and send back a message.
	void fillRequestNode(RequestNodeMsg *msg);
	/** @brief Receive the node from the cache as following a previous
	 * request which returned NULL, and continue the treewalk of the bucket
	 * which requested it with this new node.
	*/
	void receiveNode(GenericTreeNode &node, int chunk, unsigned int reqID);
	/// Just and inline version of receiveNode
	void receiveNode_inline(GenericTreeNode &node, int chunk, unsigned int reqID);
	/// @brief Find the key in the KeyTable, and copy the node over the passed pointer
	/// @todo Could the node copy be avoided?
	const GenericTreeNode* lookupNode(Tree::NodeKey key);
	/// Find the particles starting at "begin", and return a pointer to it
	const GravityParticle* lookupParticles(int begin);

	/// @brief Check if we have done with the treewalk on a specific bucket,
	/// and if we have, check also if we are done with all buckets
	inline void finishBucket(int iBucket);

	/** @brief Routine which does the tree walk on non-local nodes. It is
	 * called back for every incoming node (which are those requested to the
	 * cache during previous treewalks), and continue the treewalk from
	 * where it had been interrupted. It will possibly made other remote
	 * requests.
	 */
	void cachedWalkBucketTree(GenericTreeNode* node, int chunk,
				  BucketGravityRequest& req);
	GravityParticle *requestParticles(const Tree::NodeKey &key,int chunk,int remoteIndex,int begin,int end,BucketGravityRequest &req, bool isPrefetch=false);
	void fillRequestParticles(RequestParticleMsg *msg);
	//void fillRequestParticles(Tree::NodeKey key,int retIndex, int begin,int end,
	//			  unsigned int reqID);
	void receiveParticles(GravityParticle *part,int num,int chunk,
			      unsigned int reqID);
	void receiveParticles_inline(GravityParticle *part,int num,int chunk,
				     unsigned int reqID);
			  
	void startlb(CkCallback &cb);
	void ResumeFromSync();

	void outputAccelerations(OrientedBox<double> accelerationBox, const std::string& suffix, const CkCallback& cb);
	void outputAccASCII(OrientedBox<double> accelerationBox, const std::string& suffix, const CkCallback& cb);
	void outputStatistics(Interval<unsigned int> macInterval, Interval<unsigned int> cellInterval, Interval<unsigned int> particleInterval, Interval<unsigned int> callsInterval, double totalmass, const CkCallback& cb);
	void outputRelativeErrors(Interval<double> errorInterval, const CkCallback& cb);

	/// Collect the total statistics from the various chares
	void collectStatistics(CkCallback &cb);
	//void getPieceValues(piecedata *totaldata);

        /** @brief Entry method used to split the processing of all the buckets
         * in small pieces. It call startNextBucket a _yieldPeriod number of
         * times, and then it returns to the scheduler after enqueuing a message
         * for itself.
	 */
        void nextBucket(dummyMsg *m);

	void report(const CkCallback& cb);
	void printTree(GenericTreeNode* node, ostream& os);	
	void pup(PUP::er& p);
};

void printTree(GenericTreeNode* node, ostream& os) ;
//bool compBucket(GenericTreeNode *ln,GenericTreeNode *rn);
#endif //PARALLELGRAVITY_H
