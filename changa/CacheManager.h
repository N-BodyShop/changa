#ifndef __CACHEMANAGER_H__
#define __CACHEMANAGER_H__
#include <sys/types.h>
#include <vector>
#include <map>
#include "SFC.h"
#include "GravityTreeNode.h"
#include "charm++.h"

/** CacheEntry represents the entry for a remote 
node that is requested by the chares 
on a processor.
It stores the index of the remote chare from 
which node is to be requested and the local
chares that request it.***/



using namespace std;
using namespace Tree;
using namespace SFC;

typedef SFCTreeNode CacheNode;

class ProtoCacheEntry{
public:
	Key requestNodeID; // nodeID
	int home; // index of the array element that contains this node
	vector<int> requestorVec; // index of the array element that made the request
	vector<BucketGravityRequest *> reqVec;  //request data structure that is different for each requestor

	bool requestSent;
	bool replyRecvd;
	int totalRequests;
	int hits;
	ProtoCacheEntry(){
		replyRecvd = false;
		home = -1;
		totalRequests=0;
		hits=0;
		requestSent=false;
	}

};

class CacheEntry : public ProtoCacheEntry{
public:
	CacheNode *node;
	
	CacheEntry():ProtoCacheEntry(){
		node = NULL;
	}
	~CacheEntry(){
		delete node;
		requestorVec.clear();
		reqVec.clear();
	}
	void sendRequest(BucketGravityRequest *);
};	

class ParticleCacheEntry : public ProtoCacheEntry{
public:
	GravityParticle *part;
	int num;
	int begin;
	int end;
	ParticleCacheEntry():ProtoCacheEntry(){
		part = NULL;
		num=0;
		begin=0;
		end=0;
	}
	~ParticleCacheEntry(){
		delete []part;
		requestorVec.clear();
		reqVec.clear();
	}
	void sendRequest(BucketGravityRequest *);
};

class MapKey {
public:
	Key k;
	int home;
	MapKey(){
		k=0;
		home=0;
	};
	MapKey(Key _k,int _home){
		k = _k;
		home = _home;
	}
};
bool operator<(MapKey lhs,MapKey rhs);

class CacheManager : public Group {
private:
	/* The Cache Table is fully associative 
	A hashtable can be used as well.*/

	map<MapKey,CacheEntry *> CacheTable;
	int reqRecvd;
	int repRecvd;
	bool allReqSent;
	/** counts the total number of nodes requested by all
	the chares on that processor***/
	int totalNodesRequested;
	
	int storedNodes;

	map<MapKey,ParticleCacheEntry*> particleCacheTable;
	int storedParticles;
	bool proxyInitialized; //checks if the streaming proxy has been delegated or not
	
	void addNode(Key key,int from,CacheNode node);
	void processRequests(Key key,int from);

	public:
	CacheManager();
	~CacheManager(){};

	CacheNode *requestNode(int ,int ,Key ,BucketGravityRequest *);
	void recvNodes(Key ,int ,CacheNode );
	void recvNodes(int ,Key *,CacheNode *,int );
	
	GravityParticle *requestParticles(int ,Key ,int ,int ,int ,BucketGravityRequest *);
	void recvParticles(Key ,GravityParticle *,int ,int);
	void cacheSync();

};

#include "CacheManager.decl.h"
extern CProxy_CacheManager cacheManagerProxy;
#endif

