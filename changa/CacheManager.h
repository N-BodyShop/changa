#ifndef __CACHEMANAGER_H__
#define __CACHEMANAGER_H__
#include <sys/types.h>
#include <vector>
#include <map>
#include <set>
#include "SFC.h"
#include "GenericTreeNode.h"
#include "charm++.h"

/** NodeCacheEntry represents the entry for a remote 
node that is requested by the chares 
on a processor.
It stores the index of the remote chare from 
which node is to be requested and the local
chares that request it.***/

using namespace std;
//using namespace Tree;
//using namespace SFC;

typedef GenericTreeNode CacheNode;
typedef u_int64_t CacheKey;

class FillNodeMsg;

#include "CacheManager.decl.h"

class RequestorData {//: public CkPool<RequestorData, 128> {
 public:
  int arrayID;
  int reqID;

  RequestorData(int a, int r) {
    arrayID = a;
    reqID = r;
  }
};

class CacheEntry {
public:
	CacheKey requestID; // node or particle ID
	int home; // index of the array element that contains this node
	vector<RequestorData> requestorVec;
	//vector<int> requestorVec; // index of the array element that made the request
	//vector<BucketGravityRequest *> reqVec;  //request data structure that is different for each requestor

	bool requestSent;
	bool replyRecvd;
#if COSMO_STATS > 0
	int totalRequests;
	int hits;
#endif
	CacheEntry(){
		replyRecvd = false;
		requestSent=false;
		home = -1;
#if COSMO_STATS > 0
		totalRequests=0;
		hits=0;
#endif
	}

};

class NodeCacheEntry : public CacheEntry, public CkPool<NodeCacheEntry, 64> {
public:
	CacheNode *node;
	
	NodeCacheEntry():CacheEntry(){
		node = NULL;
	}
	~NodeCacheEntry(){
		delete node;
		CkAssert(requestorVec.empty());
		//reqVec.clear();
	}
	/// 
	void sendRequest(BucketGravityRequest *);
};	

class ParticleCacheEntry : public CacheEntry, public CkPool<ParticleCacheEntry, 64> {
public:
	GravityParticle *part;
	//int num;
	int begin;
	int end;
	ParticleCacheEntry():CacheEntry(){
		part = NULL;
		begin=0;
		end=0;
	}
	~ParticleCacheEntry(){
		delete []part;
		CkAssert(requestorVec.empty());
		//reqVec.clear();
	}
	void sendRequest(BucketGravityRequest *);
};

class MapKey {//: public CkPool<MapKey, 64> {
public:
	CacheKey k;
	int home;
	MapKey(){
		k=0;
		home=0;
	};
	MapKey(CacheKey _k,int _home){
		k = _k;
		home = _home;
	}
};
bool operator<(MapKey lhs,MapKey rhs);

class CacheManager : public CBase_CacheManager {
private:
	/* The Cache Table is fully associative 
	A hashtable can be used as well.*/

	map<CacheKey,NodeCacheEntry *> nodeCacheTable;
#if COSMO_STATS > 0
	int reqRecvd;
	int repRecvd;
	//bool allReqSent;
	/** counts the total number of nodes requested by all
	the chares on that processor***/
	int totalNodesRequested;
#endif
	/// used to generate new Nodes of the correct type (inheriting classes of CacheNode)
	CacheNode *prototype;
	/// list of TreePieces registered to this branch
	set<int> registeredChares;

	int storedNodes;
	unsigned int iterationNo;

	map<CacheKey,ParticleCacheEntry*> particleCacheTable;
	int storedParticles;
	//bool proxyInitialized; // checks if the streaming proxy has been delegated or not
	set<MapKey> outStandingRequests;
	
	/// Insert all nodes with root "node" coming from "from" into the nodeCacheTable.
	void addNodes(int from,CacheNode *node);
	//void addNode(CacheKey key,int from,CacheNode &node);
	/// @brief Check all the TreePiece buckets which requested a node, and call
	/// them back so that they can continue the treewalk.
	void processRequests(CacheNode *node,int from,int depth);
	/// @brief Fetches the Node from the correct TreePiece. If the TreePiece
	/// is in the same processor fetch it directly, otherwise send a message
	/// to the remote TreePiece::fillRequestNode
	CacheNode *sendNodeRequest(NodeCacheEntry *e,BucketGravityRequest *);
	GravityParticle *sendParticleRequest(ParticleCacheEntry *e,BucketGravityRequest *);

	public:
	CacheManager();
	~CacheManager(){};

	/** @brief Called by TreePiece to request for a particular node. It can
	 * return the node if it is already present in the cache, or call
	 * sendNodeRequest to get it. Returns null if the Node has to come from
	 * remote.
	*/
	CacheNode *requestNode(int ,int ,CacheKey ,BucketGravityRequest *);
	// Shortcut for the other recvNodes, this receives only one node
	//void recvNodes(CacheKey ,int ,CacheNode &);
	/** @brief Receive the nodes incoming from the remote
	 * TreePiece::fillRequestNode. It imports the nodes into the cache and
	 * process all pending requests.
	 */
	void recvNodes(FillNodeMsg *msg);
	
	GravityParticle *requestParticles(int ,const CacheKey ,int ,int ,int ,BucketGravityRequest *);
	void recvParticles(CacheKey ,GravityParticle *,int ,int);

	/** Invoked from the mainchare to start a new iteration. It calls the
	    prefetcher of all the chare elements residing on this processor, and
	    send a message to them to start the local computation. At this point
	    it is guaranteed that all the chares have already called
	    markPresence, since that is done in ResumeFromSync, and cacheSync is
	    called only after a reduction from ResumeFromSync.
	 */
	void cacheSync(double theta, const CkCallback& cb);
	/** Inform the CacheManager that the chare element will be resident on
	    this processor for the next iteration. It is done after treebuild in
	    the first iteration (when looking up localCache), and then it is
	    maintained until migration. Since this occur only with loadbalancer,
	    the unregistration and re-registration is done right before AtSync
	    and after ResumeFromSync.
	 */
	void markPresence(int index, GenericTreeNode*, int numChunks);
	void revokePresence(int index);

};

extern CProxy_CacheManager cacheManagerProxy;
#endif

