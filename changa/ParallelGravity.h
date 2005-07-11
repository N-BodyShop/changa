/** @file ParallelGravity.h
 */
 
#ifndef PARALLELGRAVITY_H
#define PARALLELGRAVITY_H

#include <string>
#include <map>
#include <vector>
#include <algorithm>

#include "pup_stl.h"
#include "ComlibManager.h"

#include "Vector3D.h"
#include "tree_xdr.h"
#include "SFC.h"
#include "TreeNode.h"
#include "Interval.h"

using namespace Tree;

class BucketGravityRequest {
public:
		
	SFC::Key startingNode;
	unsigned int identifier;
	unsigned int requestingPieceIndex;
	OrientedBox<double> boundingBox;
	unsigned int numParticlesInBucket;
	double *softs;
	Vector3D<double>* positions;
	Vector3D<double>* accelerations;
	double *potentials;
	unsigned int numAdditionalRequests;
	int finished;
	
	BucketGravityRequest(unsigned int bucketSize = 0) : identifier(0),
	    numParticlesInBucket(bucketSize), numAdditionalRequests(0),
	    finished(0) {
		if(numParticlesInBucket) {
			softs = new double[numParticlesInBucket];
			positions = new Vector3D<double>[numParticlesInBucket];
			accelerations = new Vector3D<double>[numParticlesInBucket];
			potentials = new double[numParticlesInBucket];
		    } else {
			positions = accelerations = 0;
			softs = potentials = 0;
			}
	}
	
	BucketGravityRequest(const BucketGravityRequest& req) {
		startingNode = req.startingNode;
		identifier = req.identifier;
		requestingPieceIndex = req.requestingPieceIndex;
		boundingBox = req.boundingBox;
		numParticlesInBucket = req.numParticlesInBucket;
		finished = req.finished;
		if(numParticlesInBucket) {
			softs = new double[numParticlesInBucket];
			positions = new Vector3D<double>[numParticlesInBucket];
			accelerations = new Vector3D<double>[numParticlesInBucket];
			potentials = new double[numParticlesInBucket];
			for(unsigned int i = 0; i < numParticlesInBucket; ++i) {
				softs[i] = req.softs[i];
				positions[i] = req.positions[i];
				accelerations[i] = req.accelerations[i];
				potentials[i] = req.potentials[i];
			}
		    } else {
			positions = accelerations = 0;
			softs = potentials = 0;
			}
		numAdditionalRequests = req.numAdditionalRequests;
	}
	
	BucketGravityRequest& operator=(const BucketGravityRequest& req) {
		startingNode = req.startingNode;
		identifier = req.identifier;
		requestingPieceIndex = req.requestingPieceIndex;
		boundingBox = req.boundingBox;
		numParticlesInBucket = req.numParticlesInBucket;
		finished = req.finished;
		delete[] positions;
		delete[] accelerations;
		if(numParticlesInBucket) {
			softs = new double[numParticlesInBucket];
			positions = new Vector3D<double>[numParticlesInBucket];
			accelerations = new Vector3D<double>[numParticlesInBucket];
			potentials = new double[numParticlesInBucket];
			for(unsigned int i = 0; i < numParticlesInBucket; ++i) {
				softs[i] = req.softs[i];
				positions[i] = req.positions[i];
				accelerations[i] = req.accelerations[i];
				potentials[i] = req.potentials[i];
			}
		    } else {
			positions = accelerations = 0;
			softs = potentials = 0;
			}
		numAdditionalRequests = req.numAdditionalRequests;
		return *this;
	}
	
	~BucketGravityRequest() {
		delete[] softs;
		delete[] positions;
		delete[] accelerations;
		delete[] potentials;
	}
	
	void merge(const BucketGravityRequest& req) {
	    for(unsigned int i = 0; i < numParticlesInBucket; ++i) {
		    accelerations[i] += req.accelerations[i];
		    potentials[i] += req.potentials[i];
	    }
	}
	
	void pup(PUP::er &p) {
		p | startingNode;
		p | identifier;
		p | requestingPieceIndex;
		p | boundingBox;
		p | numParticlesInBucket;
		if(p.isUnpacking()) {
			if(numParticlesInBucket) {
				softs = new double[numParticlesInBucket];
				positions = new Vector3D<double>[numParticlesInBucket];
				accelerations = new Vector3D<double>[numParticlesInBucket];
				potentials = new double[numParticlesInBucket];
			}
		}
		for(unsigned int i = 0; i < numParticlesInBucket; ++i) {
			p | softs[i];
			p | positions[i];
			p | accelerations[i];
			p | potentials[i];
		}	
		p | numAdditionalRequests;
	}	
};

class GravityRequest {
public:
	SFC::Key startingNode;
	unsigned int identifier;
	unsigned int requestingPieceIndex;
	double soft;
	Vector3D<double> position;
	Vector3D<double> acceleration;
	double potential;
	unsigned int numAdditionalRequests;
	unsigned int numCellInteractions;
	unsigned int numParticleInteractions;
	unsigned int numMACChecks;
	unsigned int numEntryCalls;
	
	GravityRequest() : identifier(0), requestingPieceIndex(0), numAdditionalRequests(0), numCellInteractions(0), numParticleInteractions(0), numMACChecks(0), numEntryCalls(0) { }
	
	void merge(const GravityRequest& req) {
		acceleration += req.acceleration;
		potential += req.potential;
		numCellInteractions += req.numCellInteractions;
		numParticleInteractions += req.numParticleInteractions;
		numMACChecks += req.numMACChecks;
		numEntryCalls += req.numEntryCalls;
	}
	
	void pup(PUP::er &p) {
		p | startingNode;
		p | identifier;
		p | requestingPieceIndex;
		p | soft;
		p | position;
		p | acceleration;
		p | potential;
		p | numAdditionalRequests;
		p | numCellInteractions;
		p | numParticleInteractions;
		p | numMACChecks;
		p | numEntryCalls;
	}
};

class GravityParticle {
public:

	SFC::Key key;
	double mass;
	double soft;
	Vector3D<double> position;
	Vector3D<double> acceleration;
	Vector3D<double> treeAcceleration;
	double potential;
	unsigned int numCellInteractions;
	unsigned int numParticleInteractions;
	unsigned int numMACChecks;
	unsigned int numEntryCalls;

	/*********************************/
	double intcellmass;
	double intpartmass;
	double extcellmass;
	double extpartmass;
	/********************************/
	
	GravityParticle(SFC::Key k = 0) : key(k), mass(0), numCellInteractions(0), numParticleInteractions(0), numMACChecks(0), numEntryCalls(0) { }
	GravityParticle(float m, Vector3D<float> p, Vector3D<float> a) : mass(m), position(p), acceleration(a), numCellInteractions(0), numParticleInteractions(0), numMACChecks(0), numEntryCalls(0) { }
	
	inline bool operator<(const GravityParticle& p) const {
		return key < p.key;
	}
	
	void update(const GravityRequest& req) {
		treeAcceleration += req.acceleration;
		potential += req.potential;
		numCellInteractions += req.numCellInteractions;
		numParticleInteractions += req.numParticleInteractions;
		numMACChecks += req.numMACChecks;
		numEntryCalls += req.numEntryCalls;
	}

	void pup(PUP::er &p)
	    {
		p | position;
		p | mass;
		p | soft;
		}
	
};

#include "CacheManager.h"
#include "ParallelGravity.decl.h"

extern int verbosity;
extern bool _cache;
extern int _cacheLineDepth;
extern unsigned int _yieldPeriod;

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

class Main : public Chare {
	std::string basefilename;
	CProxy_TreePiece pieces;
	unsigned int numTreePieces;
	double theta;
	unsigned int bucketSize;
	unsigned int printBinaryAcc;
public:
		
	Main(CkArgMsg* m);

	void nextStage();
};

class TreePiece : public CBase_TreePiece {
	unsigned int numTreePieces;
	CkCallback callback;
	unsigned int myNumParticles;
	GravityParticle* myParticles;
	GravityParticle* leftBoundary;
	GravityParticle* rightBoundary;
	/***************************/
	double piecemass;
	int cnt;
	int packed;
	/*************************/
	unsigned int numSplitters;
	SFC::Key* splitters;
	CProxy_TreePiece pieces;
	CProxy_TreePiece streamingProxy;
	std::string basefilename;
	OrientedBox<float> boundingBox;
	FieldHeader fh;
	bool started;
	unsigned iterationNo;
	//TreeStuff::TreeNode* root;
	SFCTreeNode* root;
	unsigned int boundaryNodesPending;
	SFCTreeNode tempNode;
	/// Opening angle
	double theta;
	unsigned int mySerialNumber;
	u_int64_t myNumParticlesPending;
	u_int64_t myNumCellInteractions;
	u_int64_t myNumParticleInteractions;
	u_int64_t myNumMACChecks;
	u_int64_t myNumProxyCalls;
	u_int64_t myNumProxyCallsBack;
	u_int64_t nextParticle;

	unsigned int numBuckets;
	unsigned int currentBucket;
	int cachecellcount;
	std::vector<SFCTreeNode *> bucketList;
	BucketGravityRequest *bucketReqs;
	
	/** A table of the nodes in my tree, indexed by their keys.
	 @todo XXX: Make this lookup a hash table, so we get O(1) behavior instead of O(log N).
	 */
	typedef std::map<SFC::Key, SFCTreeNode *> NodeLookupType;
	NodeLookupType nodeLookup;

	typedef std::map<unsigned int, GravityRequest> UnfilledRequestsType;
	UnfilledRequestsType unfilledRequests;
	typedef std::map<unsigned int, BucketGravityRequest> UnfilledBucketRequestsType;
	UnfilledBucketRequestsType unfilledBucketRequests;
	/**
	pointer to the instance of the local cache
	*/
	CacheManager *localCache;
	int countIntersects;
	int countHits;
	
	SFCTreeNode* lookupLeftChild(SFCTreeNode* node);
	SFCTreeNode* lookupRightChild(SFCTreeNode* node);
	void buildTree(SFCTreeNode* node, GravityParticle* leftParticle, GravityParticle* rightParticle);
	void calculateRemoteMoments(SFCTreeNode* node);
	void checkTree(SFCTreeNode* node);
	bool nodeOwnership(SFCTreeNode* node, unsigned int* designatedOwner = 0, unsigned int* numOwners = 0, unsigned int* firstOwner = 0, unsigned int* lastOwner = 0);
	
	void walkTree(GravityTreeNode* node, GravityRequest& req);
	void startNextParticle();
	void walkBucketTree(GravityTreeNode* node, BucketGravityRequest& req);
	void startNextBucket();
	void doAllBuckets();
	void copySFCTreeNode(SFCTreeNode &tmp,SFCTreeNode *node);
	void prefixCopyNode(SFCTreeNode *node,Key lookupKey,Key *cacheKeys,SFCTreeNode *cacheNodes,int *count,int depth);
	void rebuildSFCTree(SFCTreeNode *node,SFCTreeNode *parent,int *);
	
public:
	
	TreePiece(unsigned int numPieces) : numTreePieces(numPieces), pieces(thisArrayID), streamingProxy(thisArrayID), started(false), root(0) {
	//CkPrintf("[%d] TreePiece created\n",thisIndex);
	    // ComlibDelegateProxy(&streamingProxy);
		if(_cache){	
			localCache = cacheManagerProxy.ckLocalBranch();
		}	
		iterationNo=0;
		usesAtSync=CmiTrue;
		countIntersects=0;
		/****************/
		piecemass = 0.0;
		packed=0;
		cnt=0;
		/****************/
	}
	
	TreePiece(CkMigrateMessage* m) { 
		usesAtSync=CmiTrue;
	}
	~TreePiece() {
		delete[] myParticles;
		delete[] splitters;
	}
	
	void load(const std::string& fn, const CkCallback& cb);
	
	void buildTree(int bucketSize, const CkCallback& cb);
	
	void collectSplitters(CkReductionMsg* m);
	void startTreeBuild(CkReductionMsg* m);
	
	void acceptBoundaryNodeContribution(const SFC::Key lookupKey, const u_int64_t numParticles, const MultipoleMoments& moments);
	void acceptBoundaryNode(const SFC::Key lookupKey, const u_int64_t numParticles, const MultipoleMoments& moments);

	bool openCriterion(GravityTreeNode *node, GravityRequest& req);
	bool openCriterionBucket(GravityTreeNode *node,
				 BucketGravityRequest& req);
	
	void calculateGravityDirect(const CkCallback& cb);
	void fillRequestDirect(GravityRequest req);
	void receiveGravityDirect(const GravityRequest& req);
	
	void calculateGravityTree(double t, const CkCallback& cb);
	void fillRequestTree(GravityRequest req);	
	void receiveGravityTree(const GravityRequest& req);
	
	void calculateGravityBucketTree(double t, const CkCallback& cb);
	void fillRequestBucketTree(BucketGravityRequest req);	
	void receiveGravityBucketTree(const BucketGravityRequest& req);
	
	SFCTreeNode* requestNode(int remoteIndex, Key lookupKey,
				 BucketGravityRequest& req);
	void fillRequestNode(int retIndex, Key lookupKey,
			     BucketGravityRequest& req);
	void receiveNode(SFCTreeNode node, BucketGravityRequest& req);
	void receiveNode_inline(SFCTreeNode node, BucketGravityRequest& req);
	void lookupNode(Key ,SFCTreeNode *);

	GravityParticle* requestParticle(int remoteIndex, int iPart,
					 BucketGravityRequest& req);
	void fillRequestParticle(int retIndex, int iPart,
				 BucketGravityRequest& req);
	void receiveParticle(GravityParticle part, BucketGravityRequest& req);
	void finishBucket(int iBucket);
	void cachedWalkBucketTree(GravityTreeNode* node,
				  BucketGravityRequest& req);
	GravityParticle *requestParticles(Key &key,int remoteIndex,int begin,int end,BucketGravityRequest& req);
	void fillRequestParticles(Key key,int retIndex, int begin,int end,
				    BucketGravityRequest& req);
	void receiveParticles(GravityParticle *part,int num,
				BucketGravityRequest& req);
	void receiveParticles_inline(GravityParticle *part,int num,
				BucketGravityRequest& req);
			  
	void startlb(CkCallback &cb);
	void ResumeFromSync();

	void outputAccelerations(OrientedBox<double> accelerationBox, const std::string& suffix, const CkCallback& cb);
	void outputAccASCII(OrientedBox<double> accelerationBox, const std::string& suffix, const CkCallback& cb);
	void outputStatistics(Interval<unsigned int> macInterval, Interval<unsigned int> cellInterval, Interval<unsigned int> particleInterval, Interval<unsigned int> callsInterval, double totalmass, const CkCallback& cb);
	void outputRelativeErrors(Interval<double> errorInterval, const CkCallback& cb);

/*******************ADDED*******************/
  void getPieceValues(piecedata *totaldata);
/*******************************************/	
        void nextBucket(dummyMsg *m);

	void report(const CkCallback& cb);
	
	void pup(PUP::er& p);
};

void printTree(SFCTreeNode* node, ostream& os) ;
bool compBucket(SFCTreeNode *ln,SFCTreeNode *rn);
#endif //PARALLELGRAVITY_H
