#include "CacheInterface.h"

EntryTypeParticle::EntryTypeParticle() {
  FillDataMsg msg;
  offset = ((char*)&msg.data) - ((char*)&msg);
}

void *EntryTypeParticle::request(CkArrayIndexMax& idx, CacheKey key) {
  RequestDataMsg *msg = new (32) RequestDataMsg(key, CkMyPe());
  *(int*)CkPriorityPtr(msg) = -100000000;
  CkSetQueueing(msg, CK_QUEUEING_IFIFO);
  treeProxy[idx].fillRequestParticles(msg);
}

void *EntryTypeParticle::unpack(FillDataMsg *msg) {
  CacheParticle *data = (CacheParticle*) msg->data;
  data->msg = msg;
  return (void*) data;
}

void EntryTypeParticle::writeback(CkArrayIndexMax& idx, CacheKey k, void *data) {
  CkFreeMsg(((CacheParticle*)data)->msg);
}


void TreePiece::fillRequestParticles(RequestDataMsg *msg) {
  // the key used in the cache is shifted to the left of 1, this makes
  // a clear distinction between nodes and particles
  GenericTreeNode *bucket = lookupNode(msg->key >> 1);
  
  int total = sizeof(CacheParticle) + (bucket->lastParticle - bucket->firstParticle) * sizeof(ExternalGravityParticle);
  FillDataMsg *reply = new (total) FillDataMsg(msg->key);
  CacheParticle *data = (CacheParticle*)reply->data;
  data->begin = bucket->firstParticle;
  data->end = bucket->lastParticle;
  
  for (int i=0; i<bucket->particleCount; ++i) {
    data->part = *((ExternalGravityParticle*)&myParticles[i+bucket->firstParticle]);
  }
  
  cacheManagerProxy[msg->replyTo].recvData(reply);
  
  delete msg;
}
