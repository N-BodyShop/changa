#include "ParallelGravity.h"
#include "CacheManager.h"
#include <assert.h>
CProxy_CacheManager cacheManagerProxy;
extern CProxy_TreePiece treeProxy;
CpvDeclare(CProxy_TreePiece,streamingTreeProxy);

bool operator<(MapKey lhs,MapKey rhs){
	if(lhs.k<rhs.k){
		return true;
	}
	if(lhs.k>rhs.k){
		return false;
	}
	if(lhs.home < rhs.home){
		return true;
	}
	return false;
}


inline void CacheEntry::sendRequest(BucketGravityRequest *req){
	requestSent = true;
	treeProxy[home].fillRequestNode(CkMyPe(),requestNodeID,*req);
//	CpvAccess(streamingTreeProxy)[home].fillRequestNode(CkMyPe(),requestNodeID,*req);
};

inline void ParticleCacheEntry::sendRequest(BucketGravityRequest *req){
	requestSent = true;
	treeProxy[home].fillRequestParticles(requestNodeID,CkMyPe(),begin,end,*req);
//	CpvAccess(streamingTreeProxy)[home].fillRequestParticles(requestNodeID,CkMyPe(),begin,end,*req);
};


CacheManager::CacheManager(){
	reqRecvd=0;
	repRecvd=0;
	allReqSent = false;

	storedNodes = 0;
	totalNodesRequested = 0;
	proxyInitialized = false;
}

CacheNode *CacheManager::requestNode(int requestorIndex,int remoteIndex,Key key,BucketGravityRequest *req){
	if(!proxyInitialized){
		CpvInitialize(CProxy_TreePiece,streamingTreeProxy);
		CpvAccess(streamingTreeProxy)=treeProxy;
		ComlibDelegateProxy(&CpvAccess(streamingTreeProxy));
		proxyInitialized = true;
	}

	map<MapKey,CacheEntry *>::iterator p;
	p = CacheTable.find(MapKey(key,remoteIndex));
	CacheEntry *e;
	totalNodesRequested++;
	if(p != CacheTable.end()){
		e = p->second;
		assert(e->home == remoteIndex);
		e->totalRequests++;
		if(e->node != NULL){
			e->hits++;
			return e->node;
		}
		if(!e->requestSent){
		//sendrequest to the the tree
			e->sendRequest(req);
		}
	}else{
		e = new CacheEntry();
		e->requestNodeID = key;
		e->home = remoteIndex;
		e->totalRequests++;
		CacheTable.insert(pair<MapKey,CacheEntry *>(MapKey(key,remoteIndex),e));
		//sendrequest to the the tree
		e->sendRequest(req);
	}
	e->requestorVec.push_back(requestorIndex);
	e->reqVec.push_back(req);
	return NULL;
}

void CacheManager::addNode(Key key,int from,CacheNode node){
	map<MapKey,CacheEntry *>::iterator p;
	p = CacheTable.find(MapKey(key,from));
	if(p == CacheTable.end()){
		CacheEntry *e = new CacheEntry();
		e->requestNodeID = key;
		e->home = from;
		e->totalRequests++;
		e->requestSent = true;
		e->replyRecvd=true;
		e->node = new CacheNode;
		memcpy(e->node,&node,sizeof(CacheNode));
		CacheTable.insert(pair<MapKey,CacheEntry *>(MapKey(key,from),e));
		storedNodes++;
		return;
	}
	CacheEntry *e = p->second;
	e->node = new CacheNode;
	if(e->replyRecvd){
		CkPrintf("repeat request seen \n");
	}
	e->replyRecvd=true;
	memcpy(e->node,&node,sizeof(CacheNode));
	storedNodes++;
}

void CacheManager::processRequests(Key key,int from){
	map<MapKey,CacheEntry *>::iterator p;
	p = CacheTable.find(MapKey(key,from));
	assert(p != CacheTable.end());
	CacheEntry *e = p->second;

	vector<int>::iterator caller;
	vector<BucketGravityRequest *>::iterator callreq;
	for(caller = e->requestorVec.begin(),callreq = e->reqVec.begin();caller != e->requestorVec.end();caller++,callreq++){
		TreePiece *p = treeProxy[*caller].ckLocal();
		p->receiveNode(*(e->node),*(*callreq));
	}
	e->requestorVec.clear();
	e->reqVec.clear();

}

void CacheManager::recvNodes(Key key,int from,CacheNode node){
	addNode(key,from,node);
	processRequests(key,from);
}
void CacheManager::recvNodes(int num,Key *cacheKeys,CacheNode *cacheNodes,int index){
	for(int i=0;i<num;i++){
		addNode(cacheKeys[i],index,cacheNodes[i]);
	}
	for(int i=0;i<num;i++){
		processRequests(cacheKeys[i],index);
	}	
}

GravityParticle *CacheManager::requestParticles(int requestorIndex,Key key,int remoteIndex,int begin,int end,BucketGravityRequest *req){
	map<MapKey,ParticleCacheEntry *>::iterator p;
	p = particleCacheTable.find(MapKey(key,remoteIndex));
	ParticleCacheEntry *e;
	if(p != particleCacheTable.end()){
		e = p->second;
		assert(e->home == remoteIndex);
		e->totalRequests++;
		if(e->part != NULL){
			e->hits++;
			return e->part;
		}
		if(!e->requestSent){
			e->sendRequest(req);
		}
	}else{
		e = new ParticleCacheEntry();
		e->requestNodeID = key;
		e->home = remoteIndex;
		e->begin = begin;
		e->end=end;
		e->totalRequests++;
		particleCacheTable.insert(pair<MapKey,ParticleCacheEntry *>(MapKey(key,remoteIndex),e));
		e->sendRequest(req);
	}
	e->requestorVec.push_back(requestorIndex);
	e->reqVec.push_back(req);
	return NULL;
}

void CacheManager::recvParticles(Key key,GravityParticle *part,int num,int from){
	map<MapKey,ParticleCacheEntry *>::iterator p;
	p = particleCacheTable.find(MapKey(key,from));
	if(p == particleCacheTable.end()){
		printf("[%d] particle data received for a node that was never asked for \n",CkMyPe());
		return;
	}
	ParticleCacheEntry *e = p->second;
	e->part = new GravityParticle[num];
	memcpy(e->part,part,num*sizeof(GravityParticle));
	e->num = num;
	vector<int>::iterator caller;
	vector<BucketGravityRequest *>::iterator callreq;
	for(caller = e->requestorVec.begin(),callreq = e->reqVec.begin();caller != e->requestorVec.end();caller++,callreq++){
		TreePiece *p = treeProxy[*caller].ckLocal();
		p->receiveParticles(e->part,e->num,*(*callreq));
	}
	storedParticles+=num;
}


void CacheManager::cacheSync(){
	map<MapKey,CacheEntry *>::iterator p;
	for(p = CacheTable.begin();p != CacheTable.end();p++){
		printf("[%d] Key %llu  total number of requests %d hits %d\n",CkMyPe(),p->first.k,p->second->totalRequests,p->second->hits);
	}
	CacheTable.clear();
	particleCacheTable.clear();
	printf("[%d] Total number of requests %d storedNodes %d \n",CkMyPe(),totalNodesRequested,storedNodes);
}

