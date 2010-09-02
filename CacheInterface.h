#ifndef __CACHEINTERFACE_H__
#define __CACHEINTERFACE_H__

#include "CkCache.h"
#include "config.h"
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
  void free(void *);
  int size(void *);
  
  static void callback(CkArrayID, CkArrayIndexMax&, CkCacheKey, CkCacheUserData &, void*, int);
};

/*********************************************************
 * Smooth interface: Particles
 *********************************************************/

class CacheSmoothParticle {
public:
    int begin; // Beginning particle number
    int end;	// ending Particle number
    CkCacheKey key;
    GravityParticle *partCached;
    extraSPHData *extraSPHCached;
  ExternalSmoothParticle partExt[1];
};

class EntryTypeSmoothParticle : public CkCacheEntryType {
    // N.B. can't have helpful attributes because of the static function.
public:
  EntryTypeSmoothParticle();
  void * request(CkArrayIndexMax&, CkCacheKey);
  void * unpack(CkCacheFillMsg *, int, CkArrayIndexMax &);
  void writeback(CkArrayIndexMax&, CkCacheKey, void *);
  void free(void *);
  int size(void *);
  
  static void callback(CkArrayID, CkArrayIndexMax&, CkCacheKey, CkCacheUserData &, void*, int);
};

/*********************************************************
 * Gravity interface: Nodes
 *********************************************************/

class EntryTypeGravityNode : public CkCacheEntryType {
  void *vptr; // For saving a copy of the virtual function table.
	      // It's use will be compiler dependent.
  void unpackSingle(CkCacheFillMsg *, Tree::BinaryTreeNode *, int, CkArrayIndexMax &, bool);
public:
  EntryTypeGravityNode();
  void * request(CkArrayIndexMax&, CkCacheKey);
  void * unpack(CkCacheFillMsg *, int, CkArrayIndexMax &);
  void writeback(CkArrayIndexMax&, CkCacheKey, void *);
  void free(void *);
  int size(void *);
  CmiUInt8 getUsedBy(void *data);
  
  static void callback(CkArrayID, CkArrayIndexMax&, CkCacheKey, CkCacheUserData &, void*, int);
};

#endif


