/// @file CacheInterface.C
///
/// Implementation of the interfaces used by the CacheManager.
///
#include "CacheInterface.h"
#include "ParallelGravity.h"
#include "Opt.h"
#include "smooth.h"

EntryTypeGravityParticle::EntryTypeGravityParticle() {
  CkCacheFillMsg msg(0);
}

/// @param idx Index of the TreePiece
/// @param key Key of the requested bucket.
///
/// Calls TreePiece::fillRequestParticles() to fullfill the request.
void * EntryTypeGravityParticle::request(CkArrayIndexMax& idx, CkCacheKey key) {
  CkCacheRequestMsg *msg = new (32) CkCacheRequestMsg(key, CkMyPe());
  // This is a high priority message
  *(int*)CkPriorityPtr(msg) = -100000000;
  CkSetQueueing(msg, CK_QUEUEING_IFIFO);
  treeProxy[*idx.data()].fillRequestParticles(msg);
  return NULL;
}

/// @param msg Message containing requested data.
/// @param chunk chunk of cache
/// @param from Index of TreePiece which supplied the data
/// @return pointer to cached data
void * EntryTypeGravityParticle::unpack(CkCacheFillMsg *msg, int chunk, CkArrayIndexMax &from) {
  CacheParticle *data = (CacheParticle*) msg->data;
  data->msg = msg;
  return (void*) data;
}

void EntryTypeGravityParticle::writeback(CkArrayIndexMax& idx, CkCacheKey k, void *data) { }

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
void EntryTypeGravityParticle::callback(CkArrayID requestorID, CkArrayIndexMax &requestorIdx, CkCacheKey key, CkCacheUserData &userData, void *data, int chunk) {
  CkArrayIndex1D idx(requestorIdx.data()[0]);
  CProxyElement_TreePiece elem(requestorID, idx);
  CacheParticle *cp = (CacheParticle *)data;
  int reqID = (int)(userData.d0 & 0xFFFFFFFF);
  int awi = userData.d0 >> 32;
  void *source = (void *)userData.d1;

  elem.ckLocal()->receiveParticlesCallback(cp->part, cp->end - cp->begin + 1, chunk, reqID, key, awi, source);
}


void TreePiece::fillRequestParticles(CkCacheRequestMsg *msg) {
  // the key used in the cache is shifted to the left of 1, this makes
  // a clear distinction between nodes and particles
  const GenericTreeNode *bucket = lookupNode(msg->key >> 1);
  CkAssert(bucket != NULL);
  int total = sizeof(CacheParticle) + (bucket->lastParticle - bucket->firstParticle) * sizeof(ExternalGravityParticle);
  CkCacheFillMsg *reply = new (total) CkCacheFillMsg(msg->key);
  CkAssert(reply != NULL);
  CacheParticle *data = (CacheParticle*)reply->data;
  CkAssert(data != NULL);
  data->begin = bucket->firstParticle;
  data->end = bucket->lastParticle;
  
  for (unsigned int i=0; i<bucket->particleCount; ++i) {
    data->part[i] = *((ExternalGravityParticle*)&myParticles[i+bucket->firstParticle]);
  }
  
  cacheGravPart[msg->replyTo].recvData(reply);
  
  delete msg;
}

// Methods for "combiner" cache

EntryTypeSmoothParticle::EntryTypeSmoothParticle() {
}

void * EntryTypeSmoothParticle::request(CkArrayIndexMax& idx, CkCacheKey key) {
  CkCacheRequestMsg *msg = new (32) CkCacheRequestMsg(key, CkMyPe());
  *(int*)CkPriorityPtr(msg) = -100000000;
  CkSetQueueing(msg, CK_QUEUEING_IFIFO);
  treeProxy[*idx.data()].fillRequestSmoothParticles(msg);
  return NULL;
}

void * EntryTypeSmoothParticle::unpack(CkCacheFillMsg *msg, int chunk, CkArrayIndexMax &from) {
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
void EntryTypeSmoothParticle::writeback(CkArrayIndexMax& idx, CkCacheKey k, void *data) {
    CacheSmoothParticle *cPart = (CacheSmoothParticle *)data;
    int total = sizeof(CacheSmoothParticle)
	+ (cPart->end - cPart->begin)*sizeof(ExternalSmoothParticle);
    CkCacheFillMsg *reply = new (total) CkCacheFillMsg(cPart->key);
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

void EntryTypeSmoothParticle::callback(CkArrayID requestorID, CkArrayIndexMax &requestorIdx, CkCacheKey key, CkCacheUserData &userData, void *data, int chunk) {
  CkArrayIndex1D idx(requestorIdx.data()[0]);
  CProxyElement_TreePiece elem(requestorID, idx);
  int reqID = (int)(userData.d0 & 0xFFFFFFFF);
  int awi = userData.d0 >> 32;
  void *source = (void *)userData.d1;
  CacheSmoothParticle *cPart = (CacheSmoothParticle *)data;

  elem.ckLocal()->receiveParticlesFullCallback(cPart->partCached,
				cPart->end - cPart->begin + 1, chunk, reqID,
				key, awi, source);
}

// satisfy buffered requests

void TreePiece::processReqSmoothParticles() {
    for(SmPartRequestType::iterator iter = smPartRequests.begin();
	iter != smPartRequests.end();) {
	CkCacheKey bucketKey = iter->first;
	const GenericTreeNode *bucket = lookupNode(bucketKey >> 1);
	int total = sizeof(CacheSmoothParticle) + (bucket->lastParticle - bucket->firstParticle) * sizeof(ExternalSmoothParticle);

	CkVec<int> *vRec = iter->second;

	iter++;  // The current request gets deleted below, so
		 // increment first.
	CkCacheFillMsg *reply = new (total) CkCacheFillMsg(bucketKey);
	CacheSmoothParticle *data = (CacheSmoothParticle*)reply->data;
	data->begin = bucket->firstParticle;
	data->end = bucket->lastParticle;
	for (unsigned int ip=0; ip<bucket->particleCount; ++ip) {
	    data->partExt[ip] = myParticles[ip+bucket->firstParticle].getExternalSmoothParticle();
	    }

	for(unsigned int i = 0; i < vRec->length(); ++i) {
	    nCacheAccesses++;
	    if(i < vRec->length() - 1) { // Copy message if there is
					 // more than one outstanding request.
		CkCacheFillMsg *replyCopy = (CkCacheFillMsg *) CkCopyMsg((void **) &reply);
		cacheSmoothPart[(*vRec)[i]].recvData(replyCopy);
		}
	    else
		cacheSmoothPart[(*vRec)[i]].recvData(reply);
	    }
	delete vRec;
	smPartRequests.erase(bucketKey);
	}
    }

void TreePiece::fillRequestSmoothParticles(CkCacheRequestMsg *msg) {
    // buffer request if we are not ready for it.
    if(sSmooth == NULL) {
	CkVec<int> *vReq = smPartRequests[msg->key];
	if(vReq == NULL) {
	    vReq = new CkVec<int>();
	    smPartRequests[msg->key] = vReq;
	    }
	vReq->push_back(msg->replyTo);
	delete msg;
	return;
	}

  CkAssert(msg->replyTo != CkMyPe());
  
  // the key used in the cache is shifted to the left of 1, this makes
  // a clear distinction between nodes and particles
  const GenericTreeNode *bucket = lookupNode(msg->key >> 1);
  
  int total = sizeof(CacheSmoothParticle) + (bucket->lastParticle - bucket->firstParticle) * sizeof(ExternalSmoothParticle);
  CkCacheFillMsg *reply = new (total) CkCacheFillMsg(msg->key);
  CacheSmoothParticle *data = (CacheSmoothParticle*)reply->data;
  data->begin = bucket->firstParticle;
  data->end = bucket->lastParticle;
  
  for (unsigned int i=0; i<bucket->particleCount; ++i) {
      data->partExt[i] = myParticles[i+bucket->firstParticle].getExternalSmoothParticle();
  }
  
  nCacheAccesses++;
  
  cacheSmoothPart[msg->replyTo].recvData(reply);
  
  delete msg;
}

/// Combine cached copies with the originals on the treepiece.
/// This function also decrements the count of outstanding cache
/// accesses and does a check to see if the smooth walk is finished.
void TreePiece::flushSmoothParticles(CkCacheFillMsg *msg) {
  
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

void * EntryTypeGravityNode::request(CkArrayIndexMax& idx, CkCacheKey key) {
  CkCacheRequestMsg *msg = new (32) CkCacheRequestMsg(key, CkMyPe());
  *(int*)CkPriorityPtr(msg) = -110000000;
  CkSetQueueing(msg, CK_QUEUEING_IFIFO);
  treeProxy[*idx.data()].fillRequestNode(msg);
  return NULL;
}

void * EntryTypeGravityNode::unpack(CkCacheFillMsg *msg, int chunk, CkArrayIndexMax &from) {
  // recreate the entire tree inside this message
  Tree::BinaryTreeNode *node = (Tree::BinaryTreeNode *) (((char*)msg->data) + 8);
  node->unpackNodes();
  // recursively add all the nodes in this message to the cache
  // and link the leaves of this message to nodes in the cache (if present)
  unpackSingle(msg, node, chunk, from, true);
  
  // link node to its parent if present in the cache
  CkCacheKey ckey(node->getParentKey());
  Tree::BinaryTreeNode *parent = (Tree::BinaryTreeNode *) ((CkCacheManager *)cacheNode.ckLocalBranch())->requestDataNoFetch(ckey, chunk);
  if (parent != NULL) {
    node->parent = parent;
    parent->setChildren(parent->whichChild(node->getKey()), node);
  }

  return (void *) node;
}

void EntryTypeGravityNode::unpackSingle(CkCacheFillMsg *msg, Tree::BinaryTreeNode *node, int chunk, CkArrayIndexMax &from, bool isRoot) {

  // Store pointer to message in front of node storage so it can be
  // freed when we are done.  See free() method below.

  *(CkCacheFillMsg **) (((char*)node)-8) = msg;

  // Overwrite virtual pointer table.  Something like this will be
  // needed for heterogeneous architectures.  Commented out for now
  // since it breaks on the PGI compiler.

  // memcpy(node, &vptr, sizeof(void*));

  if (!isRoot) CmiReference(UsrToEnv(msg));
  for (int i=0; i < 2; ++i) {
    if (node->children[i] != NULL) {
      unpackSingle(msg, node->children[i], chunk, from, false);
    } else {
      CkCacheKey ckey = node->getChildKey(i);
      Tree::BinaryTreeNode *child = (Tree::BinaryTreeNode *) ((CkCacheManager *)cacheNode.ckLocalBranch())->requestDataNoFetch(ckey, chunk);
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
  CkCacheKey ckey(node->getKey());
  if (!isRoot) ((CkCacheManager*)cacheNode.ckLocalBranch())->recvData(ckey, from, (EntryTypeGravityNode*)this, chunk, (void*)node);
}

void EntryTypeGravityNode::writeback(CkArrayIndexMax& idx, CkCacheKey k, void *data) { }

void EntryTypeGravityNode::free(void *data) {
    // msg pointer is stored in front of the node data.
  CkFreeMsg(*(void **)(((char*)data)-8));
}

int EntryTypeGravityNode::size(void * data) {
  return sizeof(Tree::BinaryTreeNode);
}

void EntryTypeGravityNode::callback(CkArrayID requestorID, CkArrayIndexMax &requestorIdx, CkCacheKey key, CkCacheUserData &userData, void *data, int chunk) {
  CkArrayIndex1D idx(requestorIdx.data()[0]);
  CProxyElement_TreePiece elem(requestorID, idx);
  int reqID = (int)(userData.d0 & 0xFFFFFFFF);
  int awi = userData.d0 >> 32;
  void *source = (void *)userData.d1;
  elem.ckLocal()->receiveNodeCallback((Tree::GenericTreeNode*)data, chunk, reqID, awi, source);
}


void TreePiece::fillRequestNode(CkCacheRequestMsg *msg) {
  const Tree::GenericTreeNode* node = lookupNode(msg->key);
  //GenericTreeNode tmp;
  if(node != NULL) {
    if(_cache) {
#if 1 || defined CACHE_BUFFER_MSGS
      int count = ((Tree::BinaryTreeNode*)node)->countDepth(_cacheLineDepth);
      //FillBinaryNodeMsg *reply = new (count, 0) FillBinaryNodeMsg(thisIndex);
      // 8 extra bytes are allocated to store the msg pointer at the
      // beginning of the buffer.  See the free() and the
      // unpackSingle() method above.
      CkCacheFillMsg *reply = new (count * (sizeof(Tree::BinaryTreeNode)+8)) CkCacheFillMsg(msg->key);
      //reply->magic[0] = 0xd98cb23a;
      //new (reply->data) BinaryTreeNode[count];
      //CkPrintf("calling packing function: starting %p, magic=%p\n",reply->nodes,reply->magic);
      ((Tree::BinaryTreeNode*)node)->packNodes((Tree::BinaryTreeNode*)(reply->data+8), _cacheLineDepth, 8);
      //CkAssert(reply->magic[0] == 0xd98cb23a);
#else
      PUP::sizer p1;
      node->pup(p1, msg->depth);
      FillNodeMsg *reply = new (p1.size(), 0) FillNodeMsg(thisIndex);

      /// @TODO: check that at destination of "remoteIndex" are correct
      PUP::toMem p2((void*)reply->nodes);
      node->pup(p2, msg->depth);
      //int count = node->copyTo(reply->nodes, msg->depth);
#endif
      cacheNode[msg->replyTo].recvData(reply);
    } else {
      CkAbort("Non cached version not anymore supported, feel free to fix it!");
      //copySFCTreeNode(tmp,node);
      //streamingProxy[retIndex].receiveNode(tmp,msg->reqID);
    }
  }
  else {	// Handle NULL nodes
    CkAbort("Ok, before it handled this, but why do we have a null pointer in the tree?!?");
  }
  delete msg;
}
