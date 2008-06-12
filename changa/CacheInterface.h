#ifndef __CACHEINTERFACE_H__
#define __CACHEINTERFACE_H__

#include "CkCache.h"
#include "gravity.h"
#include "GenericTreeNode.h"

/*********************************************************
 * Gravity interface: Particles
 *********************************************************/

class CacheParticle {
public:
  CkCacheFillMsg *msg;
  int begin;
  int end;
  ExternalGravityParticle part[1];
};

class EntryTypeGravityParticle : public CkCacheEntryType {
  int offset;
public:
  EntryTypeGravityParticle();
  void * request(CkArrayIndexMax&, CkCacheKey);
  void * unpack(CkCacheFillMsg *, int, CkArrayIndexMax &);
  void writeback(CkArrayIndexMax&, CkCacheKey, void *);
  int size(void *);
  
  static void callback(CkArrayID, CkArrayIndexMax&, CkCacheKey, CmiUInt8, void*, int);
};

/*********************************************************
 * Gravity interface: Nodes
 *********************************************************/

class EntryTypeGravityNode : public CkCacheEntryType {
  void *vptr;
  void unpackSingle(CkCacheFillMsg *, Tree::BinaryTreeNode *, int, CkArrayIndexMax &, bool);
public:
  EntryTypeGravityNode();
  void * request(CkArrayIndexMax&, CkCacheKey);
  void * unpack(CkCacheFillMsg *, int, CkArrayIndexMax &);
  void writeback(CkArrayIndexMax&, CkCacheKey, void *);
  int size(void *);
  
  static void callback(CkArrayID, CkArrayIndexMax&, CkCacheKey, CmiUInt8, void*, int);
};

#endif


