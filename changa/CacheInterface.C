#include "CacheInterface.h"
#include "ParallelGravity.h"

EntryTypeGravityParticle::EntryTypeGravityParticle() {
  CkCacheFillMsg msg(0);
  offset = ((char*)&msg.data) - ((char*)&msg);
}

void * EntryTypeGravityParticle::request(CkArrayIndexMax& idx, CkCacheKey key) {
  CkCacheRequestMsg *msg = new (32) CkCacheRequestMsg(key, CkMyPe());
  *(int*)CkPriorityPtr(msg) = -100000000;
  CkSetQueueing(msg, CK_QUEUEING_IFIFO);
  treeProxy[*idx.data()].fillRequestParticles(msg);
  return NULL;
}

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

void EntryTypeGravityParticle::callback(CkArrayID requestorID, CkArrayIndexMax &requestorIdx, CkCacheKey key, CkCacheUserData &userData, void *data, int chunk) {
  CkArrayIndex1D idx(requestorIdx.data()[0]);
  CProxyElement_TreePiece elem(requestorID, idx);
  CacheParticle *cp = (CacheParticle *)data;
  int reqID = (int)(userData.d0 & 0xFFFFFFFF);
  int awi = userData.d0 >> 32;
  void *source = (void *)userData.d1;

  elem.receiveParticlesCallback(cp->part, cp->end - cp->begin + 1, chunk, reqID, key, awi, source);
}


void TreePiece::fillRequestParticles(CkCacheRequestMsg *msg) {
  // the key used in the cache is shifted to the left of 1, this makes
  // a clear distinction between nodes and particles
  const GenericTreeNode *bucket = lookupNode(msg->key >> 1);
  
  int total = sizeof(CacheParticle) + (bucket->lastParticle - bucket->firstParticle) * sizeof(ExternalGravityParticle);
  CkCacheFillMsg *reply = new (total) CkCacheFillMsg(msg->key);
  CacheParticle *data = (CacheParticle*)reply->data;
  data->begin = bucket->firstParticle;
  data->end = bucket->lastParticle;
  
  for (int i=0; i<bucket->particleCount; ++i) {
    data->part[i] = *((ExternalGravityParticle*)&myParticles[i+bucket->firstParticle]);
  }
  
  cacheManagerProxy[msg->replyTo].recvData(reply);
  
  delete msg;
}

// Methods for "combiner" cache

EntryTypeSmoothParticle::EntryTypeSmoothParticle() {
  CkCacheFillMsg msg(0);
  offset = ((char*)&msg.data) - ((char*)&msg);
}

void * EntryTypeSmoothParticle::request(CkArrayIndexMax& idx, CkCacheKey key) {
  CkCacheRequestMsg *msg = new (32) CkCacheRequestMsg(key, CkMyPe());
  *(int*)CkPriorityPtr(msg) = -100000000;
  CkSetQueueing(msg, CK_QUEUEING_IFIFO);
  treeProxy[*idx.data()].fillRequestParticles(msg);
  return NULL;
}

void * EntryTypeSmoothParticle::unpack(CkCacheFillMsg *msg, int chunk, CkArrayIndexMax &from) {
  CacheParticle *data = (CacheParticle*) msg->data;
  data->msg = msg;
  return (void*) data;
}

void EntryTypeSmoothParticle::writeback(CkArrayIndexMax& idx, CkCacheKey k, void *data) {
    // Filippo: since it destroyed the data, should be moved to "free"
    // Send the message back to the original TreePiece.
    treeProxy[*idx.data()].flushSmoothParticles(((CacheParticle*)data)->msg);
}

void EntryTypeSmoothParticle::free(void *data) { }

int EntryTypeSmoothParticle::size(void * data) {
  CacheParticle *p = (CacheParticle *) data;
  return sizeof(CacheParticle) + (p->end - p->begin) * sizeof(ExternalGravityParticle);
}

void EntryTypeSmoothParticle::callback(CkArrayID requestorID, CkArrayIndexMax &requestorIdx, CkCacheKey key, CkCacheUserData &userData, void *data, int chunk) {
  CkArrayIndex1D idx(requestorIdx.data()[0]);
  CProxyElement_TreePiece elem(requestorID, idx);
  CacheParticle *cp = (CacheParticle *)data;
  int reqID = (int)(userData.d0 & 0xFFFFFFFF);
  int awi = userData.d0 >> 32;
  void *source = (void *)userData.d1;

  elem.receiveParticlesCallback(cp->part, cp->end - cp->begin + 1, chunk, reqID, key, awi, source);
}

void TreePiece::flushSmoothParticles(CkCacheFillMsg *msg) {
  // the key used in the cache is shifted to the left of 1, this makes
  // a clear distinction between nodes and particles
  const GenericTreeNode *bucket = lookupNode(msg->key >> 1);
  
  CacheParticle *data = (CacheParticle*)msg->data;
  
  for (int i=0; i<bucket->particleCount; ++i) {
      fcnCombine(&myParticles[i+bucket->firstParticle], &data->part[i]);
  }
  
  delete msg;
}

EntryTypeGravityNode::EntryTypeGravityNode() {
  BinaryTreeNode node;
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
  Tree::BinaryTreeNode *parent = (Tree::BinaryTreeNode *) cacheManagerProxy[CkMyPe()].requestDataNoFetch(ckey, chunk);
  if (parent != NULL) {
    node->parent = parent;
    parent->setChildren(parent->whichChild(node->getKey()), node);
  }

  return (void *) node;
}

void EntryTypeGravityNode::unpackSingle(CkCacheFillMsg *msg, Tree::BinaryTreeNode *node, int chunk, CkArrayIndexMax &from, bool isRoot) {
  *(CkCacheFillMsg **) (((char*)node)-8) = msg;
  memcpy(node, &vptr, sizeof(void*));
  if (!isRoot) CmiReference(UsrToEnv(msg));
  for (int i=0; i < 2; ++i) {
    if (node->children[i] != NULL) {
      unpackSingle(msg, node->children[i], chunk, from, false);
    } else {
      CkCacheKey ckey = node->getChildKey(i);
      Tree::BinaryTreeNode *child = (Tree::BinaryTreeNode *) cacheManagerProxy[CkMyPe()].requestDataNoFetch(ckey, chunk);
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
  if (!isRoot) cacheManagerProxy[CkMyPe()].recvData(ckey, from, (EntryTypeGravityNode*)this, chunk, (void*)node);
}

void EntryTypeGravityNode::writeback(CkArrayIndexMax& idx, CkCacheKey k, void *data) { }

void EntryTypeGravityNode::free(void *data) {
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
  elem.receiveNodeCallback((Tree::GenericTreeNode*)data, chunk, reqID, awi, source);
}


void TreePiece::fillRequestNode(CkCacheRequestMsg *msg) {
  const Tree::GenericTreeNode* node = lookupNode(msg->key);
  //GenericTreeNode tmp;
  if(node != NULL) {
    if(_cache) {
#if 1 || defined CACHE_BUFFER_MSGS
      int count = ((Tree::BinaryTreeNode*)node)->countDepth(_cacheLineDepth);
      //FillBinaryNodeMsg *reply = new (count, 0) FillBinaryNodeMsg(thisIndex);
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
      cacheManagerProxy[msg->replyTo].recvData(reply);
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
