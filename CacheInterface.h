#ifndef __CACHEINTERFACE_H__
#define __CACHEINTERFACE_H__

/** @file CacheInterface.h
 *
 *  Declares the interfaces used by the CacheManager: the software
 *  cache for requesting off processor particle and node data.
 */

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

/// @brief Cache interface to particles for the gravity calculation.
/// This is a read-only cache of particles.
class EntryTypeGravityParticle : public CkCacheEntryType {
public:
  EntryTypeGravityParticle();
  /// @brief Request a bucket of particles from a TreePiece.
  void * request(CkArrayIndexMax&, CkCacheKey);
  /// @brief Return data from fufilled cache request.
  void * unpack(CkCacheFillMsg *, int, CkArrayIndexMax &);
  /// @brief Do nothing: this is a read-only cache.
  void writeback(CkArrayIndexMax&, CkCacheKey, void *);
  /// @brief free cached data.
  void free(void *);
  /// @brief return size of cached data.
  int size(void *);
  
  /// @brief callback to TreePiece after data is received.
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
  
  static void callback(CkArrayID, CkArrayIndexMax&, CkCacheKey, CkCacheUserData &, void*, int);
};

#endif


