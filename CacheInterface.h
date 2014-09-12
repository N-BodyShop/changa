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
#include "keytype.h"

/*********************************************************
 * Gravity interface: Particles
 *********************************************************/

class CacheParticle {
public:
  CkCacheFillMsg<KeyType> *msg;
  int begin;
  int end;
  ExternalGravityParticle part[1];
};

/// @brief Cache interface to particles for the gravity calculation.
/// This is a read-only cache of particles.
class EntryTypeGravityParticle : public CkCacheEntryType<KeyType> {
public:
  EntryTypeGravityParticle();
  /// @brief Request a bucket of particles from a TreePiece.
  void * request(CkArrayIndexMax&, KeyType);
  /// @brief Return data from fufilled cache request.
  void * unpack(CkCacheFillMsg<KeyType> *, int, CkArrayIndexMax &);
  /// @brief Do nothing: this is a read-only cache.
  void writeback(CkArrayIndexMax&, KeyType, void *);
  /// @brief free cached data.
  void free(void *);
  /// @brief return size of cached data.
  int size(void *);
  
  /// @brief callback to TreePiece after data is received.
  static void callback(CkArrayID, CkArrayIndexMax&, KeyType, CkCacheUserData &, void*, int);
};

/*********************************************************
 * Smooth interface: Particles
 *********************************************************/

/// @brief particle data in the smooth particle cache messages
class CacheSmoothParticle {
public:
    int begin; // Beginning particle number
    int end;	// ending Particle number
    KeyType key;
    GravityParticle *partCached;
    extraSPHData *extraSPHCached;
  ExternalSmoothParticle partExt[1];
};

/// @brief Cache interface to the particles for smooth calculations.
/// This cache is a writeback cache.
class EntryTypeSmoothParticle : public CkCacheEntryType<KeyType> {
    // N.B. can't have helpful attributes because of the static function.
public:
  EntryTypeSmoothParticle();
  /// @brief Request a bucket of particles from a TreePiece.
  void * request(CkArrayIndexMax&, KeyType);
  /// @brief Return data from fufilled cache request.
  void * unpack(CkCacheFillMsg<KeyType> *, int, CkArrayIndexMax &);
  void writeback(CkArrayIndexMax&, KeyType, void *);
  /// @brief free cached data.
  void free(void *);
  /// @brief return size of cached data.
  int size(void *);
  
  /// @brief callback to TreePiece after data is received.
  static void callback(CkArrayID, CkArrayIndexMax&, KeyType, CkCacheUserData &, void*, int);
};

/*********************************************************
 * Gravity interface: Nodes
 *********************************************************/

class EntryTypeGravityNode : public CkCacheEntryType<KeyType> {
  void *vptr; // For saving a copy of the virtual function table.
	      // It's use will be compiler dependent.
  void unpackSingle(CkCacheFillMsg<KeyType> *, Tree::BinaryTreeNode *, int, CkArrayIndexMax &, bool);
public:
  EntryTypeGravityNode();
  void * request(CkArrayIndexMax&, KeyType);
  void * unpack(CkCacheFillMsg<KeyType> *, int, CkArrayIndexMax &);
  void writeback(CkArrayIndexMax&, KeyType, void *);
  void free(void *);
  int size(void *);
  
  static void callback(CkArrayID, CkArrayIndexMax&, KeyType, CkCacheUserData &, void*, int);
};

#endif


