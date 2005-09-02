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
  treeProxy[home].fillRequestNode(new RequestNodeMsg(CkMyPe(),_cacheLineDepth,requestID,req->identifier));
  //	CpvAccess(streamingTreeProxy)[home].fillRequestNode(CkMyPe(),requestNodeID,*req);
};

inline void ParticleCacheEntry::sendRequest(BucketGravityRequest *req){
  requestSent = true;
  treeProxy[home].fillRequestParticles(requestID,CkMyPe(),begin,end,req->identifier);
  //	CpvAccess(streamingTreeProxy)[home].fillRequestParticles(requestNodeID,CkMyPe(),begin,end,*req);
};


CacheManager::CacheManager(){
#if COSMO_STATS > 0
  reqRecvd=0;
  repRecvd=0;
  //allReqSent = false;
  totalNodesRequested = 0;
#endif

  storedNodes = 0;
  storedParticles = 0;
  //proxyInitialized = false;
  iterationNo=0;
}

CacheNode *CacheManager::requestNode(int requestorIndex,int remoteIndex,CacheKey key,BucketGravityRequest *req){
  /*
    if(!proxyInitialized){
    CpvInitialize(CProxy_TreePiece,streamingTreeProxy);
    CpvAccess(streamingTreeProxy)=treeProxy;
    ComlibDelegateProxy(&CpvAccess(streamingTreeProxy));
    proxyInitialized = true;
    }
  */

  map<CacheKey,NodeCacheEntry *>::iterator p;
  p = nodeCacheTable.find(key);
  NodeCacheEntry *e;
#if COSMO_STATS > 0
  totalNodesRequested++;
#endif
  if(p != nodeCacheTable.end()){
    e = p->second;
    //assert(e->home == remoteIndex); not anymore true
#if COSMO_STATS > 0
    e->totalRequests++;
#endif

    if(e->node != NULL){
#if COSMO_STATS > 0
      e->hits++;
#endif
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
#if COSMO_STATS > 0
    e->totalRequests++;
#endif
    nodeCacheTable.insert(pair<CacheKey,NodeCacheEntry *>(key,e));
    //sendrequest to the the tree, if it is local return the node that is returned
    if(sendNodeRequest(e,req)){
      return e->node;
    }
  }
  e->requestorVec.push_back(RequestorData(requestorIndex,req->identifier));
  //e->reqVec.push_back(req);
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
    //e->node = prototype->createNew();
    const GenericTreeNode *nd = p->lookupNode(e->requestID);
#ifdef COSMO_PRINT
    if (nd==NULL) CkPrintf("Requested for node %s in %d!\n",(TreeStuff::keyBits(e->requestID,63)).c_str(),e->home);
#endif
    CkAssert(nd != NULL);
    CkAssert(e->node == NULL);
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
  map<CacheKey,NodeCacheEntry *>::iterator p;
  p = nodeCacheTable.find(node->getKey());
  NodeCacheEntry *e;
  CacheNode *oldnode = NULL;
  if (p == nodeCacheTable.end()) {
    // Completely new node, never seen in the cache
    e = new NodeCacheEntry();
    e->requestID = node->getKey();
    e->home = from;
#if COSMO_STATS > 0
    e->totalRequests++;
#endif
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
    nodeCacheTable[node->getKey()] = e;

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
    if (child != NULL) addNodes(from, child);
    else {
      if (node->getType() != Empty && node->getType() != NonLocal && node->getType() != Bucket) {
	// find a child node into the cache, using the old node if present
	if (oldnode != NULL) child = oldnode->getChildren(i);
	else {
	  pchild = nodeCacheTable.find(node->getChildKey(i));
	  if (pchild != nodeCacheTable.end() && pchild->second->node != NULL) child = pchild->second->node;
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
  if (oldnode != NULL) delete oldnode;
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

void CacheManager::processRequests(CacheNode *node,int from,int depth){
  map<CacheKey,NodeCacheEntry *>::iterator p;
  p = nodeCacheTable.find(node->getKey());
  if (p == nodeCacheTable.end()) return; // this means the node is not stored in
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
		
    p->receiveNode(*(e->node),caller->reqID);
    //treeProxy[*caller].receiveNode_inline(*(e->node),*(*callreq));
  }
  e->requestorVec.clear();
  //e->reqVec.clear();
  // iterate over the children of the node
  if (--depth <= 0) return;
  for (unsigned int i=0; i<node->numChildren(); ++i) {
    if (node->getChildren(i) != NULL)
      processRequests(node->getChildren(i), from, depth);
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
  // recursively add all nodes in the tree rooted by "node"
  addNodes(msg->owner,newnode);
  map<CacheKey,NodeCacheEntry *>::iterator e = nodeCacheTable.find(newnode->getParentKey());
  if (e != nodeCacheTable.end() && e->second->node != NULL) {
    newnode->parent = e->second->node;
    if (e->second->node->getType() == Cached) {
      e->second->node->setChildren(e->second->node->whichChild(newnode->getKey()), newnode);
    }
  }
  outStandingRequests.erase(MapKey(newnode->getKey(),msg->owner));
  // recursively process all nodes just inserted in the cache
  processRequests(newnode,msg->owner,_cacheLineDepth);
}


GravityParticle *CacheManager::requestParticles(int requestorIndex,const CacheKey key,int remoteIndex,int begin,int end,BucketGravityRequest *req){
  map<CacheKey,ParticleCacheEntry *>::iterator p;
  p = particleCacheTable.find(key);
  ParticleCacheEntry *e;
  if(p != particleCacheTable.end()){
    e = p->second;
    CkAssert(e->home == remoteIndex);
    CkAssert(e->begin == begin);
    CkAssert(e->end == end);
#if COSMO_STATS > 0
    e->totalRequests++;
#endif
    if(e->part != NULL){
#if COSMO_STATS > 0
      e->hits++;
#endif
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
#if COSMO_STATS > 0
    e->totalRequests++;
#endif
    particleCacheTable[key] = e;
    if (sendParticleRequest(e, req)) {
      return e->part;
    }
  }
  e->requestorVec.push_back(RequestorData(requestorIndex,req->identifier));
  //e->reqVec.push_back(req);
  return NULL;
}

GravityParticle *CacheManager::sendParticleRequest(ParticleCacheEntry *e, BucketGravityRequest *req) {
  TreePiece *p = treeProxy[e->home].ckLocal();
  if(p != NULL){
    e->requestSent = true;
    e->replyRecvd = true;
    const GravityParticle *gp = p->lookupParticles(e->begin);
    CkAssert(gp != NULL);
    e->part = new GravityParticle[e->end - e->begin + 1];
    memcpy(e->part, gp, (e->end - e->begin + 1)*sizeof(GravityParticle));
    return e->part;
  }
  e->sendRequest(req);
  return NULL;
}

void CacheManager::recvParticles(CacheKey key,GravityParticle *part,int num,int from){
  CkAssert(num>0);
  map<CacheKey,ParticleCacheEntry *>::iterator p;
  p = particleCacheTable.find(key);
  if(p == particleCacheTable.end()){
    printf("[%d] particle data received for a node that was never asked for \n",CkMyPe());
    return;
  }
  ParticleCacheEntry *e = p->second;
  //CkAssert(e->from == from);
  e->part = new GravityParticle[num];
  memcpy(e->part,part,num*sizeof(GravityParticle));
  //e->num = num;
  vector<RequestorData>::iterator caller;
  //vector<BucketGravityRequest *>::iterator callreq;
  for(caller = e->requestorVec.begin(); caller != e->requestorVec.end(); caller++){
    TreePiece *p = treeProxy[caller->arrayID].ckLocal();

    p->receiveParticles(e->part,num,caller->reqID);
    //treeProxy[*caller].receiveParticles_inline(e->part,e->num,*(*callreq));
  }
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
  if (verbosity)
    CkPrintf("[%d] Total number of requests %d storedNodes %d outStandingRequests %d iterationNo %d \n",CkMyPe(),totalNodesRequested,storedNodes,outStandingRequests.size(),iterationNo);
  reqRecvd=0;
  repRecvd=0;
  totalNodesRequested=0;
#endif
  map<CacheKey,NodeCacheEntry *>::iterator pn;
  for (pn = nodeCacheTable.begin(); pn != nodeCacheTable.end(); pn++) {
    NodeCacheEntry *e = pn->second;
    delete e;
  }
  nodeCacheTable.clear();
  map<CacheKey,NodeCacheEntry *>::iterator pp;
  for (pp = nodeCacheTable.begin(); pp != nodeCacheTable.end(); pp++) {
    NodeCacheEntry *e = pp->second;
    delete e;
  }
  particleCacheTable.clear();

  CkAssert(outStandingRequests.empty());
  storedNodes=0;
  storedParticles=0;
  //}
  // call the gravitational computatino utility of each local chare element
  set<int>::iterator iter;
  for (iter = registeredChares.begin(); iter != registeredChares.end(); iter++) {
    TreePiece *p = treeProxy[*iter].ckLocal();
    CkAssert(p != NULL);
    p->startIteration(theta, cb);
  }
}

void CacheManager::markPresence(int index, GenericTreeNode *proto, int numChunks) {
  prototype = proto;
  registeredChares.insert(index);
  
}

void CacheManager::revokePresence(int index) {
  registeredChares.erase(index);
}
