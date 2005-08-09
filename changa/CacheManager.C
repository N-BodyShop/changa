#include "ParallelGravity.h"
#include "CacheManager.h"
#include <assert.h>
CProxy_CacheManager cacheManagerProxy;
extern CProxy_TreePiece treeProxy;
CpvDeclare(CProxy_TreePiece,streamingTreeProxy);
extern int _cacheLineDepth;

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


inline void NodeCacheEntry::sendRequest(BucketGravityRequest *req){
	requestSent = true;
	treeProxy[home].fillRequestNode(new RequestNodeMsg(CkMyPe(),_cacheLineDepth,requestID,req->identifier));
//	CpvAccess(streamingTreeProxy)[home].fillRequestNode(CkMyPe(),requestNodeID,*req);
};

inline void ParticleCacheEntry::sendRequest(BucketGravityRequest *req){
	requestSent = true;
	treeProxy[home].fillRequestParticles(requestID,CkMyPe(),begin,end,req->identifier);
//	CpvAccess(streamingTreeProxy)[home].fillRequestParticles(requestNodeID,CkMyPe(),begin,end,*req);
};


CacheManager::CacheManager(){
	reqRecvd=0;
	repRecvd=0;
	allReqSent = false;

	storedNodes = 0;
	totalNodesRequested = 0;
	proxyInitialized = false;
	iterationNo=0;
}

CacheNode *CacheManager::requestNode(int requestorIndex,int remoteIndex,CacheKey key,BucketGravityRequest *req){
	if(!proxyInitialized){
		CpvInitialize(CProxy_TreePiece,streamingTreeProxy);
		CpvAccess(streamingTreeProxy)=treeProxy;
		ComlibDelegateProxy(&CpvAccess(streamingTreeProxy));
		proxyInitialized = true;
	}

	map<MapKey,NodeCacheEntry *>::iterator p;
	p = nodeCacheTable.find(MapKey(key,remoteIndex));
	NodeCacheEntry *e;
	totalNodesRequested++;
	if(p != nodeCacheTable.end()){
		e = p->second;
		assert(e->home == remoteIndex);
		e->totalRequests++;

		if(e->node != NULL){
			e->hits++;
			return e->node;
		}
		if(!e->requestSent && !e->replyRecvd){
		//sendrequest to the the tree, if it is local return the node 
			if(sendNodeRequest(e,req)){
				return e->node;
			}
		}
	}else{
		e = new NodeCacheEntry();
		e->requestID = key;
		e->home = remoteIndex;
		e->totalRequests++;
		nodeCacheTable.insert(pair<MapKey,NodeCacheEntry *>(MapKey(key,remoteIndex),e));
		//sendrequest to the the tree, if it is local return the node that is returned
		if(sendNodeRequest(e,req)){
			return e->node;
		}
	}
	e->requestorVec.push_back(requestorIndex);
	e->reqVec.push_back(req);
	return NULL;
}

#include "TreeNode.h"

CacheNode *CacheManager::sendNodeRequest(NodeCacheEntry *e,BucketGravityRequest *req){
	/*
		check if the array element to which the request is directed is local or not.
		If it is local then just fetch it store it and return it.
	*/
		TreePiece *p = treeProxy[e->home].ckLocal();
		if(p != NULL){
			e->requestSent = true;
			e->replyRecvd = true;
			e->node = prototype->createNew();
			const GenericTreeNode *nd = p->lookupNode(e->requestID);
#ifdef COSMO_PRINT
			if (nd==NULL) CkPrintf("Requested for node %s in %d!\n",(TreeStuff::keyBits(e->requestID,63)).c_str(),e->home);
#endif
			CkAssert(nd != NULL);
			if (e->node != NULL) e->node = nd->clone();
			if (e->node->getType() != NonLocal) e->node->remoteIndex = e->home;
			return e->node;
		}
	/*
		check if there is an outstanding request to the same array element that
		will also fetch this node because of the cacheLineDepth
	*/
		/// @TODO: despite the cacheLineDepth, the node can never be fetched if it is nonLocal to the TreePiece from which it should come!
	for(int i=0;i<_cacheLineDepth;i++){
		CacheKey k = e->requestID >> i;
		set<MapKey>::iterator pred = outStandingRequests.find(MapKey(k,e->home));
		if(pred != outStandingRequests.end()){
			return NULL;
		}
	}
	outStandingRequests.insert(MapKey(e->requestID,e->home));
	e->sendRequest(req);
	return NULL;
}


void CacheManager::addNodes(int from,CacheNode *node){
	map<MapKey,NodeCacheEntry *>::iterator p;
	p = nodeCacheTable.find(MapKey(node->getKey(),from));
	if (p == nodeCacheTable.end()) {
		NodeCacheEntry *e = new NodeCacheEntry();
		e->requestID = node->getKey();
		e->home = from;
		e->totalRequests++;
		e->requestSent = true;
		e->replyRecvd=true;
		if (node->getType() != NonLocal) node->remoteIndex = from;
		e->node = node;
		//memcpy(e->node,&node,sizeof(CacheNode));
		nodeCacheTable.insert(pair<MapKey,NodeCacheEntry *>(MapKey(node->getKey(),from),e));
		storedNodes++;
		//return;
	} else {
	  NodeCacheEntry *e = p->second;
	  if(e->replyRecvd){
	    //CkPrintf("[%d]repeat request seen for lookupKey %d\n",CkMyPe(),node->getKey());
		//return;
	  } else {
	    if (node->getType() != NonLocal) node->remoteIndex = from;
	    e->node = node; //new CacheNode;
	    e->replyRecvd=true;
	    //memcpy(e->node,&node,sizeof(CacheNode));
	    storedNodes++;
	    //debug code
	    if(node->getType() == Empty){
	      //e->node->key = key;
	    }
	  }
	}
	for (unsigned int i=0; i<node->numChildren(); ++i) {
	  if (node->getChildren(i) != NULL)
	    addNodes(from, node->getChildren(i));
	}
}

/*
void CacheManager::addNode(CacheKey key,int from,CacheNode *node){
	map<MapKey,NodeCacheEntry *>::iterator p;
	p = nodeCacheTable.find(MapKey(key,from));
	if(p == nodeCacheTable.end()){
		NodeCacheEntry *e = new NodeCacheEntry();
		e->requestID = key;
		e->home = from;
		e->totalRequests++;
		e->requestSent = true;
		e->replyRecvd=true;
		e->node = node;
		//memcpy(e->node,&node,sizeof(CacheNode));
		nodeCacheTable.insert(pair<MapKey,NodeCacheEntry *>(MapKey(key,from),e));
		storedNodes++;
		return;
	}
	NodeCacheEntry *e = p->second;
	if(e->replyRecvd){
		//CkPrintf("[%d]repeat request seen for lookupKey %d\n",CkMyPe(),key);
		return;
	}
	e->node = node; //new CacheNode;
	e->replyRecvd=true;
	//memcpy(e->node,&node,sizeof(CacheNode));
	storedNodes++;
	//debug code
	if(node.getType() == Empty){
		e->node->key = key;
	}
}
*/

void CacheManager::processRequests(CacheNode *node,int from){
	map<MapKey,NodeCacheEntry *>::iterator p;
	p = nodeCacheTable.find(MapKey(node->getKey(),from));
	CkAssert(p != nodeCacheTable.end());
	NodeCacheEntry *e = p->second;
	CkAssert(e->node != NULL);
	CkAssert(e->requestorVec.size() == e->reqVec.size());

	vector<int>::iterator caller;
	vector<BucketGravityRequest *>::iterator callreq;
	for(caller = e->requestorVec.begin(),callreq = e->reqVec.begin();caller != e->requestorVec.end();caller++,callreq++){
		TreePiece *p = treeProxy[*caller].ckLocal();
		
		p->receiveNode(*(e->node),(*callreq)->identifier);
		//treeProxy[*caller].receiveNode_inline(*(e->node),*(*callreq));
	}
	e->requestorVec.clear();
	e->reqVec.clear();
	// iterate over the children of the node
	for (unsigned int i=0; i<node->numChildren(); ++i) {
	  if (node->getChildren(i) != NULL)
	    processRequests(node->getChildren(i), from);
	}

}

/*
void CacheManager::recvNodes(Key key,int from,CacheNode &node){
  	addNode(key,from,node);
	outStandingRequests.erase(MapKey(key,from));
	processRequests(key,from);
}
*/
/* old version replaced by the one below!
void CacheManager::recvNodes(int num,Key *cacheKeys,CacheNode *cacheNodes,int index){
	for(int i=0;i<num;i++){
		addNode(cacheKeys[i],index,cacheNodes[i]);
	}
	outStandingRequests.erase(MapKey(cacheKeys[0],index));
	for(int i=0;i<num;i++){
		processRequests(cacheKeys[i],index);
	}
}
*/
void CacheManager::recvNodes(FillNodeMsg *msg){
  GenericTreeNode *node = prototype->createNew();
  PUP::fromMem p(msg->nodes);
  node->pup(p);
#ifdef COSMO_PRINT
  CkPrintf("Cache: Reveived nodes from %d\n",msg->owner);
#endif
  // recursively add all nodes in the tree rooted by "node"
  addNodes(msg->owner,node);
  outStandingRequests.erase(MapKey(node->getKey(),msg->owner));
  // recursively process all nodes just inserted in the cache
  processRequests(node,msg->owner);
}


/// @TODO If the TreePiece is local it should get them immediately without sending a full message
GravityParticle *CacheManager::requestParticles(int requestorIndex,CacheKey key,int remoteIndex,int begin,int end,BucketGravityRequest *req){
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
		e->requestID = key;
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

void CacheManager::recvParticles(CacheKey key,GravityParticle *part,int num,int from){
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

		p->receiveParticles(e->part,e->num,(*callreq)->identifier);
		//treeProxy[*caller].receiveParticles_inline(e->part,e->num,*(*callreq));
	}
	storedParticles+=num;
}


void CacheManager::cacheSync(unsigned int iter, GenericTreeNode *proto){
  prototype = proto;
	if(iter > iterationNo){
		iterationNo = iter;
		map<MapKey,NodeCacheEntry *>::iterator p;
		/*for(p = nodeCacheTable.begin();p != nodeCacheTable.end();p++){
			printf("[%d] Key %llu  total number of requests %d hits %d\n",CkMyPe(),p->first.k,p->second->totalRequests,p->second->hits);
		}*/
		CkPrintf("[%d] Total number of requests %d storedNodes %d outStandingRequests %d iterationNo %d \n",CkMyPe(),totalNodesRequested,storedNodes,outStandingRequests.size(),iterationNo);
		nodeCacheTable.clear();
		particleCacheTable.clear();
		outStandingRequests.clear();
		totalNodesRequested=0;
		storedNodes=0;
		reqRecvd=0;
		repRecvd=0;
		storedParticles=0;
		
	}	
}

