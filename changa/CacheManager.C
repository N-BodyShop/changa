#include "ParallelGravity.h"
#include "CacheManager.h"

#include "lbdb.h"

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


CacheManager::CacheManager(int size){
  numChunks = 0;
  prefetchRoots = NULL;
  nodeCacheTable = NULL;
  particleCacheTable = NULL;
  chunkAck = NULL;
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

  maxSize = (u_int64_t)size * 1024 * 1024 / (sizeof(NodeCacheEntry) + sizeof(CacheNode));
  if (verbosity) CkPrintf("Cache: accepting at most %llu nodes\n",maxSize);
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
  CkAssert(chunkAck[chunk] > 0);
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
    storedNodes ++;
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
      LDObjHandle objHandle;
		  int objstopped = 0;
      objHandle = p->timingBeforeCall(&objstopped);
      p->receiveNode(*(e->node), chunk, caller->reqID);
      p->timingAfterCall(objHandle,&objstopped);
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
  CkAssert(chunkAck[chunk] > 0);
  // recursively add all nodes in the tree rooted by "node"
#if COSMO_STATS > 0
  nodesMessages++;
#endif
  addNodes(chunk,msg->owner,newnode);
  delete msg;
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
  CkAssert(chunkAck[chunk] > 0);
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
    storedParticles += (e->end - e->begin + 1);
    return e->part;
  }
  e->sendRequest(req);
  return NULL;
}

void CacheManager::recvParticles(CacheKey key,GravityParticle *part,int num,int from){
  CkAssert(num>0);
  map<CacheKey,int>::iterator pchunk = outStandingParticleRequests.find(key);
  if(pchunk == outStandingParticleRequests.end()) {
    printf("[%d] particle data received for a node that was never asked for (%016llx))\n",CkMyPe(),key);
#if COSMO_STATS > 0
    particlesError++;
#endif
    return;
  }
  int chunk = pchunk->second;
  CkAssert(chunk >=0 && chunk < numChunks);
  CkAssert(chunkAck[chunk] > 0);
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
    else {
      LDObjHandle objHandle;
		  int objstopped = 0;
      objHandle = p->timingBeforeCall(&objstopped);
      p->receiveParticles(e->part,num,chunk,caller->reqID);
      p->timingAfterCall(objHandle,&objstopped);
    }
    //treeProxy[*caller].receiveParticles_inline(e->part,e->num,*(*callreq));
  }
  e->requestorVec.clear();
  storedParticles+=num;
}

#ifdef CACHE_TREE
GenericTreeNode *CacheManager::buildProcessorTree(int n, GenericTreeNode **gtn) {
  //CkAssert(n > 1); // the recursion should have stopped before!
  int pick = -1;
  int count = 0;
  for (int i=0; i<n; ++i) {
    NodeType nt = gtn[i]->getType();
    if (nt == Internal || nt == Bucket) {
      // we can use this directly, noone else can have it other than NL
#if COSMO_DEBUG > 0
      (*ofs) << "cache "<<CkMyPe()<<": "<<keyBits(gtn[i]->getKey(),63)<<" using Internal node"<<endl;
#endif
      return gtn[i];
    } else if (nt == Boundary) {
      // let's count up how many boundaries we find
      pick = i;
      count++;
    } else {
      // here it can be NonLocal, NonLocalBucket or Empty. In all cases nothing to do.
    }
  }
  if (count == 0) {
    // only NonLocal (or Empty). any is good
#if COSMO_DEBUG > 0
    (*ofs) << "cache "<<CkMyPe()<<": "<<keyBits(gtn[0]->getKey(),63)<<" using NonLocal node"<<endl;
#endif
    return gtn[0];
  } else if (count == 1) {
    // only one single Boundary, all others are NonLocal, use this directly
#if COSMO_DEBUG > 0
    (*ofs) << "cache "<<CkMyPe()<<": "<<keyBits(gtn[pick]->getKey(),63)<<" using Boundary node"<<endl;
#endif
return gtn[pick];
  } else {
    // more than one boundary, need recursion
    GenericTreeNode *newNode = gtn[pick]->clone();
    // keep track if all the children are internal, in which case we have to
    // change this node type too from boundary to internal
    bool isInternal = true;
#if COSMO_DEBUG > 0
    (*ofs) << "cache "<<CkMyPe()<<": "<<keyBits(newNode->getKey(),63)<<" duplicating node"<<endl;
#endif
    nodeLookupTable[newNode->getKey()] = newNode;
    GenericTreeNode **newgtn = new GenericTreeNode*[count];
    for (int child=0; child<gtn[0]->numChildren(); ++child) {
      for (int i=0, j=0; i<n; ++i) {
	if (gtn[i]->getType() == Boundary) newgtn[j++]=gtn[i]->getChildren(child);
      }
      GenericTreeNode *ch = buildProcessorTree(count, newgtn);
      newNode->setChildren(child, ch);
      if (ch->getType() == Boundary || ch->getType() == NonLocal || ch->getType() == NonLocalBucket) isInternal = false;
    }
    delete[] newgtn;
    if (isInternal) {
      newNode->setType(Internal);
#if COSMO_DEBUG > 0
      (*ofs) << "cache "<<CkMyPe()<<": "<<keyBits(newNode->getKey(),63)<<" converting to Internal"<<endl;
#endif
    }
    return newNode;
  }
}

int CacheManager::createLookupRoots(GenericTreeNode *node, Tree::NodeKey *keys) {
  // assumes that the keys are ordered in tree depth first!
  if (node->getKey() == *keys) {
    // ok, found a chunk root, we can end the recursion
    //CkPrintf("mapping key %s\n",keyBits(*keys,63).c_str());
    chunkRootTable[*keys] = node;
    return 1;
  }
  // need to continue the recursion on the children
  int count = 0;
  for (int i=0; i<node->numChildren(); ++i) {
    GenericTreeNode *child = node->getChildren(i);
    int partial;
    if (child != NULL) {
      // child present, get the recursion going
      partial = createLookupRoots(node->getChildren(i), keys);
      keys += partial;
    } else {
      // the child does not exist, count the keys falling under it
      Tree::NodeKey childKey = node->getChildKey(i);
      for (partial=0; ; ++partial, ++keys) {
	int k;
	for (k=0; k<63; ++k) {
	  if (childKey == ((*keys)>>k)) break;
	}
	if (((*keys)|(~0 << k)) == ~0) break;
	//if (k==63) break;
      }
      // add the last key found to the count
      ++partial;
      ++keys;
      //CkPrintf("missed keys of %s, %d\n",keyBits(childKey,63).c_str(),partial);
    }
    count += partial;
  }
  return count;
}
#endif

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
  maxNodes = 0;
  maxParticles = 0;
#endif

  for (int chunk=0; chunk<numChunks; ++chunk) {
    CkAssert(nodeCacheTable[chunk].empty());
    CkAssert(particleCacheTable[chunk].empty());
    CkAssert(chunkAck[chunk]==0);
  }

  int i;
  map<int,GenericTreeNode*>::iterator iter;
#ifdef CACHE_TREE
  // build a local tree inside the cache. This will be an exact superset of all
  // the trees in this processor. Only the minimum number of nodes is duplicated
#if COSMO_DEBUG > 0
  char fout[100];
  sprintf(fout,"cache.%d.%d",CkMyPe(),iterationNo);
  ofs = new ofstream(fout);
#endif
  GenericTreeNode **gtn = new GenericTreeNode*[registeredChares.size()];
  for (i=0, iter = registeredChares.begin(); iter != registeredChares.end(); iter++, ++i) {
    gtn[i] = iter->second;
#if COSMO_DEBUG > 0
    *ofs << "registered "<<iter->first<<endl;
#endif
  }
  // delete old tree
  NodeLookupType::iterator nodeIter;
  for (nodeIter = nodeLookupTable.begin(); nodeIter != nodeLookupTable.end(); nodeIter++) {
    delete nodeIter->second;
  }
  nodeLookupTable.clear();
  root = buildProcessorTree(registeredChares.size(), gtn);
#if COSMO_DEBUG > 0
  printTree(root,*ofs);
  ofs->close();
  delete ofs;
#endif
#endif

  CkAssert(outStandingRequests.empty());
  storedNodes=0;
  storedParticles=0;
  //}

  // for now the number of chunks is constant throughout the computation, it is only updated dynamically
  //int newChunks = numChunks;
  if (numChunks == 0) {
    numChunks = _numChunks;
#if CACHE_TREE > 0
    root->getChunks(_numChunks, prefetchRoots);
#else
    prototype->getChunks(_numChunks, prefetchRoots);
#endif
    chunkWeight = new u_int64_t[numChunks];
    nodeCacheTable = new map<CacheKey,NodeCacheEntry *>[numChunks];
    particleCacheTable = new map<CacheKey,ParticleCacheEntry *>[numChunks];
    chunkAck = new int[numChunks];
  }

  //for (int i=0; i<numChunks; ++i) printf("%d chunk %d: %s\n",CkMyPe(), i, keyBits(prefetchRoots[i],63).c_str());

  /*
  if (numChunks != newChunks) {
    delete[] nodeCacheTable;
    delete[] particleCacheTable;
    delete[] chunkAck;
    numChunks = newChunks;
    nodeCacheTable = new map<CacheKey,NodeCacheEntry *>[numChunks];
    particleCacheTable = new map<CacheKey,ParticleCacheEntry *>[numChunks];
    chunkAck = new int[numChunks];
  }
  */

  // update the prefetchRoots with the collected weights
  // @TODO

#if CACHE_TREE > 0
  chunkRootTable.clear();
  int numMappedRoots = createLookupRoots(root, prefetchRoots);
  CkAssert(numMappedRoots == numChunks);
#endif

  for (i=0; i<numChunks; ++i) {
    chunkAck[i] = registeredChares.size();
    chunkWeight[i] = 0;
  }

  // call the gravitational computation utility of each local chare element
  for (iter = registeredChares.begin(); iter != registeredChares.end(); iter++) {
    TreePiece *p = treeProxy[iter->first].ckLocal();
    CkAssert(p != NULL);
    p->startIteration(theta, numChunks, prefetchRoots, cb);
  }
}

void CacheManager::markPresence(int index, GenericTreeNode *root) {
  prototype = root;
  registeredChares[index] = root;
  //newChunks = _numChunks;
}

void CacheManager::revokePresence(int index) {
  registeredChares.erase(index);
}

void CacheManager::finishedChunk(int num, u_int64_t weight) {
  CkAssert(chunkAck[num] > 0);
  if (--chunkAck[num] == 0) {
    // we can safely delete the chunk from the cache
    chunkWeight[num] += weight;
#ifdef COSMO_PRINT
    CkPrintf("Deleting chunk %d from processor %d\n",num,CkMyPe());
#endif
#if COSMO_STATS > 0
    if (maxNodes < storedNodes) maxNodes = storedNodes;
    if (maxParticles < storedParticles) maxParticles = storedParticles;
    int releasedNodes=0;
    int releasedParticles=0;
#endif
    //storedNodes -= nodeCacheTable[num].size();
    //storedParticles -= particleCacheTable[num].size();
    map<CacheKey,NodeCacheEntry *>::iterator pn;
    for (pn = nodeCacheTable[num].begin(); pn != nodeCacheTable[num].end(); pn++) {
      NodeCacheEntry *e = pn->second;
      storedNodes --;
#if COSMO_STATS > 0
      releasedNodes++;
#endif
      delete e;
    }
    nodeCacheTable[num].clear();
    map<CacheKey,ParticleCacheEntry *>::iterator pp;
    for (pp = particleCacheTable[num].begin(); pp != particleCacheTable[num].end(); pp++) {
      ParticleCacheEntry *e = pp->second;
      storedParticles -= (e->end - e->begin + 1);
#if COSMO_STATS > 0
      releasedParticles += (e->end - e->begin + 1);
#endif
      delete e;
    }
    particleCacheTable[num].clear();
    CkAssert(storedNodes >= 0);
    CkAssert(storedParticles >= 0);
#if COSMO_STATS > 0
    if (verbosity)
      CkPrintf(" Cache [%d]: in iteration %d chunk %d has %d nodes and %d particles, weight %llu\n",CkMyPe(),iterationNo,num,releasedNodes,releasedParticles,weight);
#endif
#ifdef COSMO_PRINT
    CkPrintf("%d After purging chunk %d left %d nodes and %d particles\n",CkMyPe(),num,storedNodes,storedParticles);
#endif
  }
}

CkReduction::reducerType CacheStatistics::sum;

void CacheManager::collectStatistics(CkCallback& cb) {
#if COSMO_STATS > 0
  CacheStatistics cs(nodesArrived, nodesMessages, nodesDuplicated, nodesMisses,
		     nodesLocal, particlesArrived, particlesTotalArrived,
		     particlesMisses, particlesLocal, particlesError,
		     totalNodesRequested, totalParticlesRequested,
		     maxNodes, maxParticles, CkMyPe());
  contribute(sizeof(CacheStatistics), &cs, CacheStatistics::sum, cb);
#else
  CkAbort("Invalid call, only valid if COSMO_STATS is defined");
#endif
}
