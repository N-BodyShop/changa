#ifndef __CACHEINTERFACE_H__
#define __CACHEINTERFACE_H__

#include "CacheManager.h"

/*********************************************************
 * Gravity interface: Particles
 *********************************************************/

class CacheParticle {
public:
  FillDataMsg *msg;
  int begin;
  int end;
  ExternalGravityParticle part[1];
};

class EntryTypeParticle : public Cache::EntryType {
  int offset;
public:
  EntryTypeParticle();
  void * request(CkArrayIndexMax&, CacheKey);
  void * unpack(FillDataMsg *);
  void writeback(CkArrayIndexMax&, CacheKey, void *);
};

/*********************************************************
 * Gravity interface: Nodes
 *********************************************************/

#endif


