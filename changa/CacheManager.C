#include "ParallelGravity.h"
#include "CacheManager.h"

CProxy_CacheManager cacheManagerProxy;
extern CProxy_TreePiece treeProxy;
//CpvDeclare(CProxy_TreePiece,streamingTreeProxy);
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
  RequestNodeMsg *msg = new (32) RequestNodeMsg(CkMyPe(),_cacheLineDepth,requestID,req->identifier);
  *(int*)CkPriorityPtr(msg) = -110000000; 
  CkSetQueueing(msg, CK_QUEUEING_IFIFO);
  treeProxy[home].fillRequestNode(msg);
  //	CpvAccess(streamingTreeProxy)[home].fillRequestNode(CkMyPe(),requestNodeID,*req);
};

inline void ParticleCacheEntry::sendRequest(BucketGravityRequest *req){
  requestSent = true;
  RequestParticleMsg *msg = new (32) RequestParticleMsg(CkMyPe(),begin,end,requestID,req->identifier);
  *(int*)CkPriorityPtr(msg) = -100000000;
  CkSetQueueing(msg, CK_QUEUEING_IFIFO);
  treeProxy[home].fillRequestParticles(msg);
  //treeProxy[home].fillRequestParticles(requestID,CkMyPe(),begin,end,req->identifier);
  //	CpvAccess(streamingTreeProxy)[home].fillRequestParticles(requestNodeID,CkMyPe(),begin,end,*req);
};


CacheManager::CacheManager(){
  numChunks = 0;
  nodeCacheTable = NULL;
  particleCacheTable = NULL;
  storedNodes = 0;
  storedParticles = 0;
  //proxyInitialized = false;
  iterationNo=0;
#if COSMO_STATS > 0
  nodesArrived = 0;
  nodesMessages = 0;
  nodesDuplicated = 0;
  nodesMisses = 0;
  nodesLocal = 0;
  particlesArrived = 0;
  particlesTotalArrived = 0;
  particlesMisses = 0;
  particlesLocal = 0;
  totalNodesRequested = 0;
  totalParticlesRequested = 0;
#endif
}

CacheNode *CacheManager::requestNode(int requestorIndex,int remoteIndex,int chunk,CacheKey key,BucketGravityRequest *req,bool isPrefetch){
  /*
    if(!proxyInitialized){
    CpvInitialize(CProxy_TreePiece,streamingTreeProxy);
    CpvAccess(streamingTreeProxy)=treeProxy;
    ComlibDelegateProxy(&CpvAccess(streamingTreeProxy));
    proxyInitialized = true;
    }
  */

  map<CacheKey,NodeCacheEntry *>::iterator p;
  p = nodeCacheTable[chunk].find(key);
  NodeCacheEntry *e;
#if COSMO_STATS > 0
  totalNodesRequested++;
#endif
  if(p != nodeCacheTable[chunk].end()){
    e = p->second;
    //assert(e->home == remoteIndex); not anymore true
#if COSMO_STATS > 1
    e->totalRequests++;
#endif

    if(e->node != NULL){
      return e->node;
    }
    if(!e->requestSent && !e->replyRecvd){
      //sendrequest to the the tree, if it is local return the node 
      if(sendNodeRequest(chunk,e,req)){
	return e->node;
      }
    }
  }else{
    e = new NodeCacheEntry();
    e->requestID = key;
    e->home = remoteIndex;
#if COSMO_STATS > 1
    e->totalRequests++;
#endif
    nodeCacheTable[chunk].insert(pair<CacheKey,NodeCacheEntry *>(key,e));
    //sendrequest to the the tree, if it is local return the node that is returned
    if(sendNodeRequest(chunk,e,req)){
      return e->node;
    }
  }
  e->requestorVec.push_back(RequestorData(requestorIndex,req->identifier,isPrefetch));
  //e->reqVec.push_back(req);
#if COSMO_STATS > 1
  e->misses++;
#endif
#if COSMO_STATS > 0
  if (!isPrefetch) nodesMisses++;
#endif
  return NULL;
}

#include "TreeNode.h"

CacheNode *CacheManager::sendNodeRequest(int chunk, NodeCacheEntry *e,BucketGravityRequest *req){
  /*
    check if the array element to which the request is directed is local or not.
    If it is local then just fetch it store it and return it.
  */
  TreePiece *p = treeProxy[e->home].ckLocal();
  if(p != NULL){
    e->requestSent = true;
    e->replyRecvd = true;
    //e->node = prototype->createNew();
    const GenericTreeNode *nd = p->lookupNode(e->requestID);
#ifdef COSMO_PRINT
    if (nd==NULL) CkPrintf("Requested for node %s in %d!\n",(TreeStuff::keyBits(e->requestID,63)).c_str(),e->home);
#endif
    CkAssert(nd != NULL);
    CkAssert(e->node == NULL);
#if COSMO_STATS > 0
    nodesLocal++;
#endif
    e->node = nd->clone();
    switch (e->node->getType()) {
    case Bucket:
    case NonLocalBucket:
      e->node->setType(CachedBucket);
      break;
    case Empty:
      e->node->setType(CachedEmpty);
      break;
    case Invalid:
      break;
    default:
      e->node->setType(Cached);
    }
    return e->node;
  }
  /*
    check if there is an outstanding request to the same array element that
    will also fetch this node because of the cacheLineDepth
  */
  /// @TODO: despite the cacheLineDepth, the node can never be fetched if it is nonLocal to the TreePiece from which it should come!
  for(int i=0;i<_cacheLineDepth;i++){
    CacheKey k = e->requestID >> i;
    map<MapKey,int>::iterator pred = outStandingRequests.find(MapKey(k,e->home));
    if(pred != outStandingRequests.end()){
      return NULL;
    }
  }
  outStandingRequests[MapKey(e->requestID,e->home)] = chunk;
  e->sendRequest(req);
  return NULL;
}

void CacheManager::addNodes(int chunk,int from,CacheNode *node){
  map<CacheKey,NodeCacheEntry *>::iterator p;
  p = nodeCacheTable[chunk].find(node->getKey());
  NodeCacheEntry *e;
  CacheNode *oldnode = NULL;
#if COSMO_STATS > 0
  nodesArrived++;
#endif
  if (p == nodeCacheTable[chunk].end()) {
    // Completely new node, never seen in the cache
    e = new NodeCacheEntry();
    e->requestID = node->getKey();
    e->home = from;
    e->requestSent = true;
    e->replyRecvd = true;
    e->node = node;
    switch (node->getType()) {
    case Bucket:
    case NonLocalBucket:
      node->setType(CachedBucket);
      break;
    case Empty:
      node->setType(CachedEmpty);
      break;
    case Invalid:
      break;
    default:
      node->setType(Cached);
    }
    //memcpy(e->node,&node,sizeof(CacheNode));
    nodeCacheTable[chunk][node->getKey()] = e;

  } else {
    // The node placehoder is present, we substitute the node it contains with
    // the new one just arrived. The information this new node carries might be
    // newer or the same, but being the tree coherent across TreePieces, the
    // information is never wrong.
    // NOTE: inside the cache the type of a node is useless!!!
    e = p->second;
    oldnode = e->node;
    e->node = node;
    switch (node->getType()) {
    case Bucket:
    case NonLocalBucket:
      node->setType(CachedBucket);
      break;
    case Empty:
      node->setType(CachedEmpty);
      break;
    case Invalid:
      break;
    default:
      node->setType(Cached);
    }
    e->replyRecvd = true;
  }

  // recursively add the children. if a child is not present, check if it is
  // because the particular node has node (Empty, NonLocal, Bucket) or because
  // we stopped at cacheDepth. In this second case check if we have already a
  // child node stored in the cache, and if yes add it to the current node.
  map<CacheKey,NodeCacheEntry *>::iterator pchild;
  for (unsigned int i=0; i<node->numChildren(); ++i) {
    CacheNode *child = node->getChildren(i);
    if (child != NULL) addNodes(chunk,from, child);
    else {
      if (node->getType() != Empty && node->getType() != NonLocal && node->getType() != Bucket) {
	// find a child node into the cache, using the old node if present
	if (oldnode != NULL) child = oldnode->getChildren(i);
	else {
	  pchild = nodeCacheTable[chunk].find(node->getChildKey(i));
	  if (pchild != nodeCacheTable[chunk].end() && pchild->second->node != NULL) child = pchild->second->node;
	  else child = NULL;
	}
	if (child != NULL) {
	  // fix the parent/child pointers
	  if (child->getType() == Cached || child->getType() == CachedBucket || child->getType() == CachedEmpty) child->parent = e->node;
	  e->node->setChildren(i, child);
	}
      }
    }
  }
  // delete the old node (replaced) or count that we have inserted a new node
  if (oldnode != NULL) {
#if COSMO_STATS > 0
    nodesDuplicated++;
#endif
    delete oldnode;
  }
  else storedNodes++;
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

void CacheManager::processRequests(int chunk,CacheNode *node,int from,int depth){
  map<CacheKey,NodeCacheEntry *>::iterator p;
  p = nodeCacheTable[chunk].find(node->getKey());
  if (p == nodeCacheTable[chunk].end()) return; // this means the node is not stored in
					 // the cache, but is owned by some
					 // chare in this processor!
  NodeCacheEntry *e = p->second;
  if (e->node == NULL) return; // this is a strange condition where the node the
			       // parent pointed to is owned by some chare in
			       // this processor (like before), but the request
			       // for this node has been made by someone else
			       // (and the request has been sent out, not
			       // checking that a local chare had it)

  //CkAssert(e->requestorVec.size() == e->reqVec.size());

  vector<RequestorData>::iterator caller;
  //vector<BucketGravityRequest *>::iterator callreq;
  for(caller = e->requestorVec.begin(); caller != e->requestorVec.end(); caller++){
    TreePiece *p = treeProxy[caller->arrayID].ckLocal();

    if (caller->isPrefetch) p->prefetch(e->node);
    else {
      p->receiveNode(*(e->node),caller->reqID);
      //ckout <<"received node for computation"<<endl;
    }
    //treeProxy[*caller].receiveNode_inline(*(e->node),*(*callreq));
  }
  e->requestorVec.clear();
  //e->reqVec.clear();
  // iterate over the children of the node
  if (--depth <= 0) return;
  for (unsigned int i=0; i<node->numChildren(); ++i) {
    if (node->getChildren(i) != NULL)
      processRequests(chunk,node->getChildren(i), from, depth);
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
  GenericTreeNode *newnode = prototype->createNew();
  PUP::fromMem p(msg->nodes);
  newnode->pup(p);
#ifdef COSMO_PRINT
  CkPrintf("Cache: Reveived nodes from %d\n",msg->owner);
#endif
  map<MapKey,int>::iterator pchunk = outStandingRequests.find(MapKey(newnode->getKey(),msg->owner));
  CkAssert(pchunk != outStandingRequests.end());
  int chunk = pchunk->second;
  // recursively add all nodes in the tree rooted by "node"
#if COSMO_STATS > 0
  nodesMessages++;
#endif
  addNodes(chunk,msg->owner,newnode);
  map<CacheKey,NodeCacheEntry *>::iterator e = nodeCacheTable[chunk].find(newnode->getParentKey());
  if (e != nodeCacheTable[chunk].end() && e->second->node != NULL) {
    newnode->parent = e->second->node;
    if (e->second->node->getType() == Cached) {
      e->second->node->setChildren(e->second->node->whichChild(newnode->getKey()), newnode);
    }
  }
  outStandingRequests.erase(pchunk);
  // recursively process all nodes just inserted in the cache
  processRequests(chunk,newnode,msg->owner,_cacheLineDepth);
}


GravityParticle *CacheManager::requestParticles(int requestorIndex,int chunk,const CacheKey key,int remoteIndex,int begin,int end,BucketGravityRequest *req,bool isPrefetch){
  map<CacheKey,ParticleCacheEntry *>::iterator p;
  p = particleCacheTable[chunk].find(key);
  ParticleCacheEntry *e;
#if COSMO_STATS > 0
  totalParticlesRequested++;
#endif
  if(p != particleCacheTable[chunk].end()){
    e = p->second;
    CkAssert(e->home == remoteIndex);
    CkAssert(e->begin == begin);
    CkAssert(e->end == end);
#if COSMO_STATS > 1
    e->totalRequests++;
#endif
    if(e->part != NULL){
      return e->part;
    }
    if(!e->requestSent){
      if (sendParticleRequest(e,req)) {
	return e->part;
      }
    }
  }else{
    e = new ParticleCacheEntry();
    e->requestID = key;
    e->home = remoteIndex;
    e->begin = begin;
    e->end = end;
#if COSMO_STATS > 1
    e->totalRequests++;
#endif
    particleCacheTable[chunk][key] = e;
    if (sendParticleRequest(e, req)) {
      return e->part;
    }
  }

  e->requestorVec.push_back(RequestorData(requestorIndex,req->identifier,isPrefetch));
  //e->reqVec.push_back(req);
  outStandingParticleRequests[key] = chunk;
#if COSMO_STATS > 1
  e->misses++;
#endif
#if COSMO_STATS > 0
  if (!isPrefetch) particlesMisses++;
#endif
  return NULL;
}

GravityParticle *CacheManager::sendParticleRequest(ParticleCacheEntry *e, BucketGravityRequest *req) {
  TreePiece *p = treeProxy[e->home].ckLocal();
  if(p != NULL){
    e->requestSent = true;
    e->replyRecvd = true;
    const GravityParticle *gp = p->lookupParticles(e->begin);
    CkAssert(gp != NULL);
#if COSMO_STATS > 0
    particlesLocal++;
#endif
    e->part = new GravityParticle[e->end - e->begin + 1];
    memcpy(e->part, gp, (e->end - e->begin + 1)*sizeof(GravityParticle));
    return e->part;
  }
  e->sendRequest(req);
  return NULL;
}

void CacheManager::recvParticles(CacheKey key,GravityParticle *part,int num,int from){
  CkAssert(num>0);
  map<CacheKey,int>::iterator pchunk = outStandingParticleRequests.find(key);
  if(pchunk == outStandingParticleRequests.end()) {
    printf("[%d] particle data received for a node that was never asked for \n",CkMyPe());
#if COSMO_STATS > 0
    particlesError++;
#endif
    return;
  }
  int chunk = pchunk->second;
  CkAssert(chunk >=0 && chunk < numChunks);
  outStandingParticleRequests.erase(key);

  map<CacheKey,ParticleCacheEntry *>::iterator p;
  p = particleCacheTable[chunk].find(key);
  if(p == particleCacheTable[chunk].end()){
    printf("[%d] particle data received for a node that was never asked for in this chunk\n",CkMyPe());
#if COSMO_STATS > 0
    particlesError++;
#endif
    return;
  }
#if COSMO_STATS > 0
  particlesArrived++;
  particlesTotalArrived += num;
#endif
  ParticleCacheEntry *e = p->second;
  //CkAssert(e->from == from);
  e->part = new GravityParticle[num];
  memcpy(e->part,part,num*sizeof(GravityParticle));
  //e->num = num;
  vector<RequestorData>::iterator caller;
  //vector<BucketGravityRequest *>::iterator callreq;
  for(caller = e->requestorVec.begin(); caller != e->requestorVec.end(); caller++){
    TreePiece *p = treeProxy[caller->arrayID].ckLocal();

    //CkAssert(!caller->isPrefetch);
    if (caller->isPrefetch) p->prefetch(e->part);
    else p->receiveParticles(e->part,num,caller->reqID);
    //treeProxy[*caller].receiveParticles_inline(e->part,e->num,*(*callreq));
  }
  e->requestorVec.clear();
  storedParticles+=num;
}


void CacheManager::cacheSync(double theta, const CkCallback& cb) {
  //if(iter > iterationNo){
  iterationNo++;
  /*map<MapKey,NodeCacheEntry *>::iterator p;
    for(p = nodeCacheTable.begin();p != nodeCacheTable.end();p++){
    printf("[%d] Key %llu  total number of requests %d hits %d\n",CkMyPe(),p->first.k,p->second->totalRequests,p->second->hits);
    }*/
#if COSMO_STATS > 0
  //if (verbosity)
  //  CkPrintf("[%d] Total number of requests %d storedNodes %d outStandingRequests %d iterationNo %d \n",CkMyPe(),totalNodesRequested,storedNodes,outStandingRequests.size(),iterationNo);
  nodesArrived = 0;
  nodesDuplicated = 0;
  nodesMisses = 0;
  nodesLocal = 0;
  particlesArrived = 0;
  particlesTotalArrived = 0;
  particlesMisses = 0;
  particlesLocal = 0;
  totalNodesRequested = 0;
  totalParticlesRequested = 0;
#endif
  for (int chunk=0; chunk<numChunks; ++chunk) {
    map<CacheKey,NodeCacheEntry *>::iterator pn;
    for (pn = nodeCacheTable[chunk].begin(); pn != nodeCacheTable[chunk].end(); pn++) {
      NodeCacheEntry *e = pn->second;
      delete e;
    }
    nodeCacheTable[chunk].clear();
    map<CacheKey,ParticleCacheEntry *>::iterator pp;
    for (pp = particleCacheTable[chunk].begin(); pp != particleCacheTable[chunk].end(); pp++) {
      ParticleCacheEntry *e = pp->second;
      delete e;
    }
    particleCacheTable[chunk].clear();
  }

  CkAssert(outStandingRequests.empty());
  storedNodes=0;
  storedParticles=0;
  //}
  if (numChunks != newChunks) {
    delete[] nodeCacheTable;
    delete[] particleCacheTable;
    numChunks = newChunks;
    nodeCacheTable = new map<CacheKey,NodeCacheEntry *>[numChunks];
    particleCacheTable = new map<CacheKey,ParticleCacheEntry *>[numChunks];
  }

  // call the gravitational computation utility of each local chare element
  set<int>::iterator iter;
  for (iter = registeredChares.begin(); iter != registeredChares.end(); iter++) {
    TreePiece *p = treeProxy[*iter].ckLocal();
    CkAssert(p != NULL);
    p->startIteration(theta, cb);
  }
}

void CacheManager::markPresence(int index, GenericTreeNode *proto, int _numChunks) {
  prototype = proto;
  registeredChares.insert(index);
  newChunks = _numChunks;
}

void CacheManager::revokePresence(int index) {
  registeredChares.erase(index);
}

CkReduction::reducerType CacheStatistics::sum;

void CacheManager::collectStatistics(CkCallback& cb) {
#if COSMO_STATS > 0
  CacheStatistics cs(nodesArrived, nodesMessages, nodesDuplicated, nodesMisses,
		     nodesLocal, particlesArrived, particlesTotalArrived,
		     particlesMisses, particlesLocal, particlesError,
		     totalNodesRequested, totalParticlesRequested);
  contribute(sizeof(CacheStatistics), &cs, CacheStatistics::sum, cb);
#else
  CkAbort("Invalid call, only valid if COSMO_STATS is defined");
#endif
}
