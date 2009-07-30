#ifndef __SSEDEFS_H__
#define __SSEDEFS_H__

#ifdef COSMO_FLOAT
#define CONVERT_TO_COSMO_TYPE (float)  
typedef float cosmoType;
#define COSMO_CONST(val) val##f
#else
#define CONVERT_TO_COSMO_TYPE
typedef double cosmoType; 
#define COSMO_CONST(val) val 
#endif

#if defined(__SSE2__) && defined(COSMO_FLOAT)
#include "SSE-Float.h"
#define SSE_VECTOR_WIDTH 4
#define FORCE_INPUT_LIST_PAD 3
typedef SSEFloat SSEcosmoType; 
#define SSELoad(where, p1, p2) where(p1 p2, p1 +1 p2, p1 +2 p2, p1 +3 p2)
#define SSEStore(what, p1, p2) { \
  float p[4]; \
  storeu(p, what); \
  p1 p2 = p[0]; \
  p1 +1 p2 = p[1]; \
  p1 +2 p2 = p[2]; \
  p1 +3 p2 = p[3]; \
}
enum { cosmoMask=0x7 };
#elif defined(__SSE2__) && !defined(COSMO_FLOAT)
#include "SSE-Double.h"
#define SSE_VECTOR_WIDTH 2
#define FORCE_INPUT_LIST_PAD 1
typedef SSEDouble SSEcosmoType; 
#define SSELoad(where, p1, p2) where(p1 p2, p1 +1 p2)
#define SSEStore(what, p1, p2) storel(&p1 p2, what); storeh(&p1 +1 p2, what)
enum { cosmoMask=0x3 };
#endif

#endif
