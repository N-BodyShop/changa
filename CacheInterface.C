/// @file CacheInterface.C
///
/// Implementation of the interfaces used by the CacheManager.
///
#include "CacheInterface.h"
#include "ParallelGravity.h"
#include "Opt.h"
#include "smooth.h"
#include "Compute.h"
#include "TreeWalk.h"

static const int PAD_reply = sizeof(NodeKey);  // Assume this is bigger
                                           // than a pointer

EntryTypeGravityParticle::EntryTypeGravityParticle() {
  CkCacheFillMsg<KeyType> msg(0);
}

/// @param idx Index of the TreePiece
/// @param key Key of the requested bucket.
///
/// Calls TreePiece::fillRequestParticles() to fullfill the request.
void * EntryTypeGravityParticle::request(CkArrayIndexMax& idx, KeyType key) {
  CkCacheRequest req(key, CkMyPe(), particleRequest);
  aggregator.ckLocalBranch()->insertData(req, *idx.data());
  return NULL;
}

/// @param msg Message containing requested data.
/// @param chunk chunk of cache
/// @param from Index of TreePiece which supplied the data
/// @return pointer to cached data
void * EntryTypeGravityParticle::unpack(CkCacheFillMsg<KeyType> *msg, int chunk, CkArrayIndexMax &from) {
  CacheParticle *data = (CacheParticle*) msg->data;
  data->msg = msg;
  return (void*) data;
}

void EntryTypeGravityParticle::writeback(CkArrayIndexMax& idx, KeyType k, void *data) { }

void EntryTypeGravityParticle::free(void *data) {
  CkFreeMsg(((CacheParticle*)data)->msg);
}

int EntryTypeGravityParticle::size(void * data) {
  CacheParticle *p = (CacheParticle *) data;
  return sizeof(CacheParticle) + (p->end - p->begin) * sizeof(ExternalGravityParticle);
}

/// @param requestorID ID of TreePiece CkArray
/// @param requestorIdx Index of requesting TreePiece
/// @param key Key of bucket requested
/// @param userData Encodes which bucket and which walk requested this
/// data.
/// @param data Pointer to the cached data
/// @param chunk Which chunk is this for.
void EntryTypeGravityParticle::callback(CkArrayID requestorID, CkArrayIndexMax &requestorIdx, KeyType key, CkCacheUserData &userData, void *data, int chunk) {
  CkArrayIndex1D idx(requestorIdx.data()[0]);
  CProxyElement_TreePiece elem(requestorID, idx);
  CacheParticle *cp = (CacheParticle *)data;
  int reqID = (int)(userData.d0 & 0xFFFFFFFF);
  int awi = userData.d0 >> 32;
  void *source = (void *)userData.d1;

  elem.receiveParticlesCallback(cp->part, cp->end - cp->begin + 1, chunk, reqID, key, awi, source);
}

void TreePiece::fillRequestParticles(const CkCacheRequest &req) {
  // the key used in the cache is shifted to the left of 1, this makes
  // a clear distinction between nodes and particles
  const GenericTreeNode *bucket = lookupNode(req.key >> 1);
  CkAssert(bucket != NULL);
  int total = sizeof(CacheParticle) + (bucket->lastParticle - bucket->firstParticle) * sizeof(ExternalGravityParticle);
  CkCacheFillMsg<KeyType> *reply = new (total) CkCacheFillMsg<KeyType>(req.key);
  CkAssert(reply != NULL);
  CacheParticle *data = (CacheParticle*)reply->data;
  CkAssert(data != NULL);
  data->begin = bucket->firstParticle;
  data->end = bucket->lastParticle;
  
  for (unsigned int i=0; i<bucket->particleCount; ++i) {
    data->part[i] = *((ExternalGravityParticle*)&myParticles[i+bucket->firstParticle]);
  }
  
  cacheGravPart[req.replyTo].recvData(reply);
}

// Methods for "combiner" cache

EntryTypeSmoothParticle::EntryTypeSmoothParticle() {
}

void * EntryTypeSmoothParticle::request(CkArrayIndexMax& idx, KeyType key) {
  CkCacheRequest req(key, CkMyPe(), smoothParticleRequest);
  aggregator.ckLocalBranch()->insertData(req, *idx.data());
  return NULL;
}

void * EntryTypeSmoothParticle::unpack(CkCacheFillMsg<KeyType> *msg, int chunk, CkArrayIndexMax &from) {
    // incoming data
    CacheSmoothParticle *cPartsIn = (CacheSmoothParticle*) msg->data;
    // Cached copy
    CacheSmoothParticle *cParts = new CacheSmoothParticle;
    
    cParts->begin = cPartsIn->begin;
    cParts->end = cPartsIn->end;
    cParts->key = msg->key;
    int nTotal = 1 + cParts->end - cParts->begin;
    cParts->partCached = new GravityParticle[nTotal];
    cParts->extraSPHCached = new extraSPHData[nTotal];
    // Expand External particles to full particles in cache
    for(int i = 0; i < nTotal; i++) {
	cParts->partCached[i].extraData = &cParts->extraSPHCached[i];
	cPartsIn->partExt[i].getParticle(&cParts->partCached[i]);
      	if(TYPETest(&(cParts->partCached[i]), globalSmoothParams->iType))
	   globalSmoothParams->initSmoothCache(&(cParts->partCached[i]));	// Clear cached copy
	}
    CkFreeMsg(msg);
    return (void*) cParts;
}

/// Send the message back to the original TreePiece.
void EntryTypeSmoothParticle::writeback(CkArrayIndexMax& idx, KeyType k, void *data) {
    CacheSmoothParticle *cPart = (CacheSmoothParticle *)data;
    int total = sizeof(CacheSmoothParticle)
	+ (cPart->end - cPart->begin)*sizeof(ExternalSmoothParticle);
    CkCacheFillMsg<KeyType> *reply = new (total) CkCacheFillMsg<KeyType>(cPart->key);
    CacheSmoothParticle *rdata = (CacheSmoothParticle*)reply->data;
    rdata->begin = cPart->begin;
    rdata->end = cPart->end;
  
    for (int i=0; i < 1 + cPart->end - cPart->begin; ++i) {
	rdata->partExt[i] = cPart->partCached[i].getExternalSmoothParticle();
	}
    treeProxy[*idx.data()].flushSmoothParticles(reply);
}

void EntryTypeSmoothParticle::free(void *data) {
    CacheSmoothParticle *cPart = (CacheSmoothParticle *)data;
    delete[] cPart->partCached;
    delete[] cPart->extraSPHCached;
    delete cPart;
}

int EntryTypeSmoothParticle::size(void * data) {
    CacheSmoothParticle *cPart = (CacheSmoothParticle *)data;
    return sizeof(CacheSmoothParticle)
	+ (cPart->end - cPart->begin)*sizeof(ExternalSmoothParticle);
}

void EntryTypeSmoothParticle::callback(CkArrayID requestorID, CkArrayIndexMax &requestorIdx, KeyType key, CkCacheUserData &userData, void *data, int chunk) {
  CkArrayIndex1D idx(requestorIdx.data()[0]);
  CProxyElement_TreePiece elem(requestorID, idx);
  int reqID = (int)(userData.d0 & 0xFFFFFFFF);
  int awi = userData.d0 >> 32;
  void *source = (void *)userData.d1;
  CacheSmoothParticle *cPart = (CacheSmoothParticle *)data;

  elem.receiveParticlesFullCallback(cPart->partCached,
				cPart->end - cPart->begin + 1, chunk, reqID,
				key, awi, source);
}

// satisfy buffered requests

void TreePiece::processReqSmoothParticles() {
    for(SmPartRequestType::iterator iter = smPartRequests.begin();
	iter != smPartRequests.end();) {
	KeyType bucketKey = iter->first;
	const GenericTreeNode *bucket = lookupNode(bucketKey >> 1);
	int total = sizeof(CacheSmoothParticle) + (bucket->lastParticle - bucket->firstParticle) * sizeof(ExternalSmoothParticle);

	CkVec<int> *vRec = iter->second;

	iter++;  // The current request gets deleted below, so
		 // increment first.
	CkCacheFillMsg<KeyType> *reply = new (total, 8*sizeof(int)) CkCacheFillMsg<KeyType>(bucketKey);
	CacheSmoothParticle *data = (CacheSmoothParticle*)reply->data;
	data->begin = bucket->firstParticle;
	data->end = bucket->lastParticle;
	for (unsigned int ip=0; ip<bucket->particleCount; ++ip) {
	    data->partExt[ip] = myParticles[ip+bucket->firstParticle].getExternalSmoothParticle();
	    }

	*(int*)CkPriorityPtr(reply) = -10000000;
	CkSetQueueing(reply, CK_QUEUEING_IFIFO);
	for(unsigned int i = 0; i < vRec->length(); ++i) {
	    nCacheAccesses++;
	    if(i < vRec->length() - 1) { // Copy message if there is
					 // more than one outstanding request.
		CkCacheFillMsg<KeyType> *replyCopy = (CkCacheFillMsg<KeyType> *) CkCopyMsg((void **) &reply);
		*(int*)CkPriorityPtr(replyCopy) = -10000000;
		CkSetQueueing(replyCopy, CK_QUEUEING_IFIFO);
		cacheSmoothPart[(*vRec)[i]].recvData(replyCopy);
		}
	    else
		cacheSmoothPart[(*vRec)[i]].recvData(reply);
	    }
	delete vRec;
	smPartRequests.erase(bucketKey);
	}
    }

void TreePiece::fillRequestSmoothParticles(const CkCacheRequest &req) {
    // buffer request if we are not ready for it.
    if(sSmooth == NULL) {
	CkVec<int> *vReq = smPartRequests[req.key];
	if(vReq == NULL) {
	    vReq = new CkVec<int>();
	    smPartRequests[req.key] = vReq;
	    }
	vReq->push_back(req.replyTo);
	return;
	}

    CkAssert(req.replyTo != CkMyPe());
  
  // the key used in the cache is shifted to the left of 1, this makes
  // a clear distinction between nodes and particles
  const GenericTreeNode *bucket = lookupNode(req.key >> 1);
  
  int total = sizeof(CacheSmoothParticle) + (bucket->lastParticle - bucket->firstParticle) * sizeof(ExternalSmoothParticle);
  CkCacheFillMsg<KeyType> *reply = new (total, 8*sizeof(int)) CkCacheFillMsg<KeyType>(req.key);
  CacheSmoothParticle *data = (CacheSmoothParticle*)reply->data;
  data->begin = bucket->firstParticle;
  data->end = bucket->lastParticle;
  
  for (unsigned int i=0; i<bucket->particleCount; ++i) {
      data->partExt[i] = myParticles[i+bucket->firstParticle].getExternalSmoothParticle();
  }
  
  nCacheAccesses++;
  
  *(int*)CkPriorityPtr(reply) = -10000000;
  CkSetQueueing(reply, CK_QUEUEING_IFIFO);
  cacheSmoothPart[req.replyTo].recvData(reply); 
}

/// Combine cached copies with the originals on the treepiece.
/// This function also decrements the count of outstanding cache
/// accesses and does a check to see if the smooth walk is finished.
void TreePiece::flushSmoothParticles(CkCacheFillMsg<KeyType> *msg) {
  
  CacheSmoothParticle *data = (CacheSmoothParticle*)msg->data;
  SmoothCompute *sc = sSmooth;
  
  CkAssert(sc != NULL);
  CkAssert(nCacheAccesses > 0);
  
  int j = 0;
  for(int i = data->begin; i <= data->end; i++) {
      if(TYPETest(&myParticles[i], sc->params->iType))
	  sc->params->combSmoothCache(&myParticles[i], &data->partExt[j]);
      j++;
      }
  
  nCacheAccesses--;
  delete msg;

  if(sSmoothState->bWalkDonePending)
      finishSmoothWalk();
}

// Node Cache methods

EntryTypeGravityNode::EntryTypeGravityNode() {
  BinaryTreeNode node;
  // save the virtual function table.
  // Note that this is compiler dependent; also note that it is unused
  // at the moment -- see unpackSingle() below.
  memcpy(&vptr, &node, sizeof(void*));
}

void * EntryTypeGravityNode::request(CkArrayIndexMax& idx, KeyType key) {
  CkCacheRequest req(key, CkMyPe(), nodeRequest);
  aggregator.ckLocalBranch()->insertData(req, *idx.data());
  return NULL;
}

void * EntryTypeGravityNode::unpack(CkCacheFillMsg<KeyType> *msg, int chunk, CkArrayIndexMax &from) {
  // recreate the entire tree inside this message
  Tree::BinaryTreeNode *node = (Tree::BinaryTreeNode *) (((char*)msg->data) + PAD_reply);
  node->unpackNodes();
  // recursively add all the nodes in this message to the cache
  // and link the leaves of this message to nodes in the cache (if present)
  unpackSingle(msg, node, chunk, from, true);
  
  // link node to its parent if present in the cache
  KeyType ckey(node->getParentKey());
  Tree::BinaryTreeNode *parent = (Tree::BinaryTreeNode *) cacheNode.ckLocalBranch()->requestDataNoFetch(ckey, chunk);
  if (parent != NULL) {
    node->parent = parent;
    parent->setChildren(parent->whichChild(node->getKey()), node);
  }

  return (void *) node;
}

void EntryTypeGravityNode::unpackSingle(CkCacheFillMsg<KeyType> *msg, Tree::BinaryTreeNode *node, int chunk, CkArrayIndexMax &from, bool isRoot) {

  // Store pointer to message in front of node storage so it can be
  // freed when we are done.  See free() method below.

  *(CkCacheFillMsg<KeyType> **) (((char*)node)-PAD_reply) = msg;

  // Overwrite virtual pointer table.  Something like this will be
  // needed for heterogeneous architectures.  Commented out for now
  // since it breaks on the PGI compiler.

  // memcpy(node, &vptr, sizeof(void*));

  if (!isRoot) CmiReference(UsrToEnv(msg));
  for (int i=0; i < 2; ++i) {
    if (node->children[i] != NULL) {
      unpackSingle(msg, node->children[i], chunk, from, false);
    } else {
      KeyType ckey = node->getChildKey(i);
      Tree::BinaryTreeNode *child = (Tree::BinaryTreeNode *) cacheNode.ckLocalBranch()->requestDataNoFetch(ckey, chunk);
      if (child != NULL) {
        child->parent = node;
        node->setChildren(node->whichChild(child->getKey()), child);
      }
    }
  }
  switch (node->getType()) {
  case Tree::Bucket:
  case Tree::NonLocalBucket:
    node->setType(Tree::CachedBucket);
    break;
  case Tree::Empty:
    node->setType(Tree::CachedEmpty);
    break;
  case Tree::Invalid:
    break;
  default:
    node->setType(Tree::Cached);
  }
  KeyType ckey(node->getKey());
  if (!isRoot) cacheNode.ckLocalBranch()->recvData(ckey, from, (EntryTypeGravityNode*)this, chunk, (void*)node);
}

void EntryTypeGravityNode::writeback(CkArrayIndexMax& idx, KeyType k, void *data) { }

void EntryTypeGravityNode::free(void *data) {
    // msg pointer is stored in front of the node data.
  CkFreeMsg(*(void **)(((char*)data)-PAD_reply));
}

int EntryTypeGravityNode::size(void * data) {
  return sizeof(Tree::BinaryTreeNode);
}

void EntryTypeGravityNode::callback(CkArrayID requestorID, CkArrayIndexMax &requestorIdx, KeyType key, CkCacheUserData &userData, void *data, int chunk) {
  CkArrayIndex1D idx(requestorIdx.data()[0]);
  CProxyElement_TreePiece elem(requestorID, idx);
  int reqID = (int)(userData.d0 & 0xFFFFFFFF);
  int awi = userData.d0 >> 32;
  void *source = (void *)userData.d1;
  elem.receiveNodeCallback((Tree::GenericTreeNode*)data, chunk, reqID, awi, source);
}

void TreePiece::process(const CkCacheRequest &req) {
  switch(req.type) {
  case nodeRequest:
    fillRequestNode(req);
    break;
  case particleRequest:
    fillRequestParticles(req);
    break;
  case smoothParticleRequest:
    fillRequestSmoothParticles(req);
    break;
  }
}

void TreePiece::fillRequestNode(CkCacheRequestMsg<KeyType> *msg) {
  const Tree::GenericTreeNode* node = lookupNode(msg->key);
  //GenericTreeNode tmp;
  if(node != NULL) {
    if(_cache) {
      int count = ((Tree::BinaryTreeNode*)node)->countDepth(_cacheLineDepth);
      // Extra bytes are allocated to store the msg pointer at the
      // beginning of the buffer.  See the free() and the
      // unpackSingle() method above.
      CkAssert(sizeof(msg) <= PAD_reply);  // be sure there is enough rooom
      //CkCacheFillMsg<KeyType> *reply = new (count * (sizeof(Tree::BinaryTreeNode)+PAD_reply), 8*sizeof(int)) CkCacheFillMsg<KeyType>(msg->key);
      CkCacheFillMsg<KeyType> *reply = new (count * ALIGN_DEFAULT(sizeof(Tree::BinaryTreeNode)+PAD_reply), 8*sizeof(int)) CkCacheFillMsg<KeyType>(msg->key);
      ((Tree::BinaryTreeNode*)node)->packNodes((Tree::BinaryTreeNode*)(reply->data+PAD_reply), _cacheLineDepth, PAD_reply);
      *(int*)CkPriorityPtr(reply) = -10000000;
      CkSetQueueing(reply, CK_QUEUEING_IFIFO);
      cacheNode[msg->replyTo].recvData(reply);
    } else {
      CkAbort("Non cached version not anymore supported, feel free to fix it!");
      //copySFCTreeNode(tmp,node);
      //streamingProxy[retIndex].receiveNode(tmp,msg->reqID);
    }
    delete msg;
  }
  else {// Handle NULL nodes
    if (!sendFillReqNodeWhenNull(msg)) {
      CkAbort("Ok, before it handled this, but why do we have a null pointer in the tree?!?");
    }
  }
}

void TreePiece::fillRequestNode(const CkCacheRequest &req) {
  const Tree::GenericTreeNode* node = lookupNode(req.key);
  //GenericTreeNode tmp;
  if(node != NULL) {
    if(_cache) {
      int count = ((Tree::BinaryTreeNode*)node)->countDepth(_cacheLineDepth);
      // Extra bytes are allocated to store the msg pointer at the
      // beginning of the buffer.  See the free() and the
      // unpackSingle() method above.
      CkCacheFillMsg<KeyType> *reply = new (count * ALIGN_DEFAULT(sizeof(Tree::BinaryTreeNode)+PAD_reply), 8*sizeof(int)) CkCacheFillMsg<KeyType>(req.key);
      ((Tree::BinaryTreeNode*)node)->packNodes((Tree::BinaryTreeNode*)(reply->data+PAD_reply), _cacheLineDepth, PAD_reply);
      *(int*)CkPriorityPtr(reply) = -10000000;
      CkSetQueueing(reply, CK_QUEUEING_IFIFO);
      cacheNode[req.replyTo].recvData(reply);
    } else {
      CkAbort("Non cached version not anymore supported, feel free to fix it!");
      //copySFCTreeNode(tmp,node);
      //streamingProxy[retIndex].receiveNode(tmp,req.reqID);
    }
  }
  else {	// Handle NULL nodes
    if (!sendFillReqNodeWhenNull(req)) {
      CkAbort("Ok, before it handled this, but why do we have a null pointer in the tree?!?");
    }
  }
}
