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

extern int verbosity;
extern bool _cache;
extern int _cacheLineDepth;
extern unsigned int _yieldPeriod;
extern DomainsDec domainDecomposition;
extern GenericTrees useTree;
extern CProxy_TreePiece streamingProxy;

class dummyMsg : public CMessage_dummyMsg{
public:
int val;
};

/********************************************/
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
/********************************************/

class RequestNodeMsg : public CMessage_RequestNodeMsg {
 public:
  int retIndex;
  int depth;
  Tree::NodeKey key;
  unsigned int reqID;

  RequestNodeMsg(int r, int d, Tree::NodeKey k, unsigned int req) : retIndex(r), depth(d), key(k), reqID(req) {}
};

class FillNodeMsg : public CMessage_FillNodeMsg {
 public:
  int owner;
  //int count;
  char *nodes;

  FillNodeMsg(int index) : owner(index) { }

  /*
  static FillNodeMsg::alloc(int msgnum, size_t sz, int *sizes, int pb) {
    int offsets[2];
    offsets[0] = ALIGN8(sz - 1); // -1 because we take out the "char nodes[1]"
    if(sizes==0) offsets[1] = offsets[0];
    else offsets[1] = offsets[0] + ALIGN8(sizes[0]);
    FillNodeMsg *newmsg = (FillNodeMsg *) CkAllocMsg(msgnum, offsets[1], pb);
    return (void *) newmsg;
  }
  */
};

class Main : public Chare {
	std::string basefilename;
	CProxy_TreePiece pieces;
	unsigned int numTreePieces;
	double theta;
	unsigned int bucketSize;
	unsigned int printBinaryAcc;
public:
		
	Main(CkArgMsg* m);

	/**
	 * Principal method which does all the coordination of the work.
	 * It is threaded, so it suspends while waiting for the TreePieces to
	 * perform the task. It calls TreePieces for: 1) load, 2) buildTree, 3)
	 * calculateGravity and startlb for the desired number of iterations.
	 */
	void nextStage();
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
	//GravityParticle* leftBoundary;
	//GravityParticle* rightBoundary;

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
	std::string basefilename;
	OrientedBox<float> boundingBox;
	FieldHeader fh;
	bool started;
	unsigned iterationNo;
	/// The root of the global tree, always local to any chare
	GenericTreeNode* root;

	/// Count for how many boundaries are still not received while building
	/// the tree
	unsigned int boundaryNodesPending;
	typedef std::map<NodeKey, CkVec<int>* >   MomentRequestType;
	/// Keep track of the requests for remote moments not yet satisfied.
	/// Used only during the tree construction.
	MomentRequestType momentRequests;

	/// Opening angle
	double theta;

	/// Map between Keys and TreeNodes, used to get a node from a key
	NodeLookupType nodeLookupTable;

	/// @if ALL

	/// unused - SEND VERSION
	unsigned int mySerialNumber;

	/// @endif

	/// Number of particles which are still traversing the tree
	u_int64_t myNumParticlesPending;

	/// @if STATISTICS

#if COSMO_STATS > 0
	u_int64_t myNumCellInteractions;
	u_int64_t myNumParticleInteractions;
	u_int64_t myNumMACChecks;
	u_int64_t myNumProxyCalls;
	u_int64_t myNumProxyCallsBack;
	/// Same as myNumCellInteractions, only restricted to cached nodes
	int cachecellcount;
	int countIntersects;
	int countHits;
#endif

	/// @endif

	/// @if ALL

	/// this was used when particles were traversing the tree (and not buckets)
	u_int64_t nextParticle;

	/// @endif

	/// Size of bucketList, total number of buckets present
	unsigned int numBuckets;
	/// Used to start the computation for all buckets, one after the other
	unsigned int currentBucket;
	/// List of all the node-buckets in this TreePiece
	std::vector<GenericTreeNode *> bucketList;
	/// @brief Used as a placeholder while traversing the tree and computing
	/// forces. This should not exist in the cache version since all its
	/// information is duplicate from the correspondent TreeNode (in
	/// bucketList), and myParticles.
	/// @todo Eliminate the usage of this in the cache walk
	BucketGravityRequest *bucketReqs;
	
	/// @if ALL

	// all these are used only in the SEND VERSION
	typedef std::map<unsigned int, GravityRequest> UnfilledRequestsType;
	UnfilledRequestsType unfilledRequests;
	typedef std::map<unsigned int, BucketGravityRequest> UnfilledBucketRequestsType;
	UnfilledBucketRequestsType unfilledBucketRequests;

	/// @endif

	/// Pointer to the instance of the local cache
	CacheManager *localCache;

	/*
	SFCTreeNode* lookupLeftChild(SFCTreeNode* node);
	SFCTreeNode* lookupRightChild(SFCTreeNode* node);
	*/

	/// convert a key to a node using the nodeLookupTable
	inline GenericTreeNode *keyToNode(const Tree::NodeKey);

 public:

	/* DEBUGGING
	void quiescence() { 
	  CkPrintf("[%d] quiescence detected, pending %d\n",thisIndex,myNumParticlesPending);
	  for (int i=0; i<numBuckets; ++i) {
	    if (bucketReqs[i].numAdditionalRequests != 0)
	      CkPrintf("[%d] requests for %d remaining %d\n",thisIndex,i,bucketReqs[i].numAdditionalRequests);
	  }
	  CkExit();
	}
	END DEBUGGING */

	// Recursive call to build the subtree with root "node", and SFC bounded by the two particles
	//void buildOctTree(GenericTreeNode* node, GravityParticle* leftParticle, GravityParticle* rightParticle);

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
	bool nodeOwnership(const GenericTreeNode *const node, int &firstOwner, int &lastOwner);
	//bool nodeOwnership(SFCTreeNode* node, unsigned int* designatedOwner = 0, unsigned int* numOwners = 0, unsigned int* firstOwner = 0, unsigned int* lastOwner = 0);

	/// @if ALL

	/// unused - SEND VERSION - Recursive tree walk over a tree
	void walkTree(GenericTreeNode* node, GravityRequest& req);
	/// unused - SEND VERSION - Sending the next particle for computation
	void startNextParticle();

	/// @endif

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
	/// Copy one node into the position passed be "tmp" (used to fill a request)
	//void copySFCTreeNode(SFCTreeNode &tmp,SFCTreeNode *node);
	/// Recursive function to copy nodes into the filling request (to be returned to the requester)
	//void prefixCopyNode(SFCTreeNode *node,Key lookupKey,Key *cacheKeys,SFCTreeNode *cacheNodes,int *count,int depth);
	void rebuildSFCTree(GenericTreeNode *node,GenericTreeNode *parent,int *);
	
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
#if COSMO_STATS > 0
	  countIntersects=0;
	  piecemass = 0.0;
	  packed=0;
	  cnt=0;
#endif
	}
	
	TreePiece(CkMigrateMessage* m) { 
		usesAtSync=CmiTrue;
	}
	~TreePiece() {
		delete[] myParticles;
		delete[] splitters;
	}
	
	void load(const std::string& fn, const CkCallback& cb);

	/// Charm entry point to build the tree (called by Main), calls collectSplitters
	void buildTree(int bucketSize, const CkCallback& cb);

	/// Collect the boundaries of all TreePieces, and trigger the real treebuild
	void collectSplitters(CkReductionMsg* m);
	/// Real tree build, independent of other TreePieces; calls the recursive buildTree
	void startOctTreeBuild(CkReductionMsg* m);

	/// @if ALL
	
	/// Receive a contribution for the multipole computation of a boundary
	/// node. This contribution is received only by one of the co-owners,
	/// which will later forward the result through acceptBoundaryNode.
	void acceptBoundaryNodeContribution(const Tree::NodeKey key, const int numParticles, const MultipoleMoments& moments);
	/// Receive the multipole for a particular boundary node, after all
	/// contributions by sharers have been added.
	void acceptBoundaryNode(const Tree::NodeKey key, const int numParticles, const MultipoleMoments& moments);

	/// @endif

	/// Request the TreePiece to send back later the moments for this node.
	/// Since it is [inline], it should do few work. The reply sent back
	/// cannot be [inline]
	void requestRemoteMoments(const Tree::NodeKey key, int sender);
	void receiveRemoteMoments(const Tree::NodeKey key, Tree::NodeType type, int firstParticle, int numParticles, const MultipoleMoments& moments);

	/// @if ALL

	/// unused - SEND VERSION
	bool openCriterion(GenericTreeNode *node, GravityRequest& req);

	/// @endif

	/// Decide whether the node should be opened for the force computation
	/// of the given request.
	bool openCriterionBucket(GenericTreeNode *node,
				 BucketGravityRequest& req);

	/// @if ALL

	/// unused (only for correctness check) - Used to perform n^2 (particle-particle) computation
	void calculateGravityDirect(const CkCallback& cb);
	/// unused (only for correctness check) - Filler for n^2 (particle-particle) computation
	void fillRequestDirect(GravityRequest req);
	/// unused (only for correctness check) - Receive incoming particle for computation n^2
	void receiveGravityDirect(const GravityRequest& req);

	/// unused - Treewalk one particle at a time
	void calculateGravityTree(double t, const CkCallback& cb);
	/// unused - Fill the request, almost identical to fillRequestBucketTree
	void fillRequestTree(GravityRequest req);	
	/// unused - Receive incoming data for particle tree walk
	void receiveGravityTree(const GravityRequest& req);

	/// unused - SEND VERSION - Return to the requester the data (typically the cache manager)
	void fillRequestBucketTree(BucketGravityRequest req);
	/// unused - SEND VERSION - callback where the data requested to another TreePiece will come back.
	void receiveGravityBucketTree(const BucketGravityRequest& req);

	/// @brief Main entry point to start gravity computation
	void calculateGravityBucketTree(double t, const CkCallback& cb);

	/// @endif

	/// Entry point for the local computation: for each bucket compute the
	/// force that its particles see due to the other particles hosted in
	/// this TreePiece. The opening angle theta has already been passed
	/// through "startIteration"
	void calculateGravityLocal();
	/// Entry point for the remote computation: for each bucket compute the
	/// force that its particles see due to the other particles NOT hosted
	/// by this TreePiece, and belonging to a subset of the global tree
	/// (specified by chunkNum).
	void calculateGravityRemote(int chunkNum);
	/// Temporary function to recurse over all the buckets like in
	/// walkBucketTree, only that NonLocal nodes are the only one for which
	/// forces are computed
	void walkBucketRemoteTree(GenericTreeNode *node, BucketGravityRequest &req);
	/// Function called by the CacheManager to start a new iteration.
	/// @param t the opening angle
	/// @param cb the callback to use after all the computation has finished
	void startIteration(double t, const CkCallback& cb);
	/// Function called by the CacheManager to send out request for needed
	/// remote data, so that the later computation will hit.
	void prefetch(GenericTreeNode *node);

	/// @brief Retrieve the remote node, goes through the cache if present
	GenericTreeNode* requestNode(int remoteIndex, Tree::NodeKey lookupKey,
				 BucketGravityRequest& req);
	/// @brief Receive a request for Nodes from a remote processor, copy the
	/// data into it, and send back a message.
	void fillRequestNode(RequestNodeMsg *msg);
	/** @brief Receive the node from the cache as following a previous
	 * request which returned NULL, and continue the treewalk of the bucket
	 * which requested it with this new node.
	*/
	void receiveNode(GenericTreeNode &node, unsigned int reqID);
	/// Just and inline version of receiveNode
	void receiveNode_inline(GenericTreeNode &node, unsigned int reqID);
	/// @brief Find the key in the KeyTable, and copy the node over the passed pointer
	/// @todo Could the node copy be avoided?
	const GenericTreeNode* lookupNode(Tree::NodeKey key);
	/// Find the particles starting at "begin", and return a pointer to it
	const GravityParticle* lookupParticles(int begin);

	/// @if ALL

	/// unused - highly inefficient version that requests one single particle
	GravityParticle* requestParticle(int remoteIndex, int iPart,
					 BucketGravityRequest& req);
	/// unused - highly inefficient version that returns one single particle
	void fillRequestParticle(int retIndex, int iPart,
				 BucketGravityRequest& req);
	/// unused - highly inefficient version that receives one single particle
	void receiveParticle(GravityParticle part, BucketGravityRequest& req);

	/// @endif

	/// @brief Check if we have done with the treewalk on a specific bucket,
	/// and if we have, check also if we are done with all buckets
	inline void finishBucket(int iBucket);
	/** @brief Routine which does the tree walk on non-local nodes. It is
	 * called back for every incoming node (which are those requested to the
	 * cache during previous treewalks), and continue the treewalk from
	 * where it had been interrupted. It will possibly made other remote
	 * requests.
	 */
	void cachedWalkBucketTree(GenericTreeNode* node,
				  BucketGravityRequest& req);
	GravityParticle *requestParticles(const Tree::NodeKey &key,int remoteIndex,int begin,int end,BucketGravityRequest &req);
	void fillRequestParticles(Tree::NodeKey key,int retIndex, int begin,int end,
				  unsigned int reqID);
	void receiveParticles(GravityParticle *part,int num,
			      unsigned int reqID);
	void receiveParticles_inline(GravityParticle *part,int num,
				     unsigned int reqID);
			  
	void startlb(CkCallback &cb);
	void ResumeFromSync();

	void outputAccelerations(OrientedBox<double> accelerationBox, const std::string& suffix, const CkCallback& cb);
	void outputAccASCII(OrientedBox<double> accelerationBox, const std::string& suffix, const CkCallback& cb);
	void outputStatistics(Interval<unsigned int> macInterval, Interval<unsigned int> cellInterval, Interval<unsigned int> particleInterval, Interval<unsigned int> callsInterval, double totalmass, const CkCallback& cb);
	void outputRelativeErrors(Interval<double> errorInterval, const CkCallback& cb);

	/// Collect the total statistics from the various chares
	void getPieceValues(piecedata *totaldata);

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
bool compBucket(GenericTreeNode *ln,GenericTreeNode *rn);
#endif //PARALLELGRAVITY_H
