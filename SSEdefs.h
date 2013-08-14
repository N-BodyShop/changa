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

#if CMK_USE_FMA4 && defined(__FMA4__)

#ifndef COSMO_FLOAT

#include "SSE-Double.h"
#define SSE_VECTOR_WIDTH 4
#define FORCE_INPUT_LIST_PAD 3
typedef SSEDouble SSEcosmoType; 
#define SSELoad(where, arr, idx, field) where(arr[idx]field, arr[idx+1]field, arr[idx+2]field, arr[idx+3]field)
#define SSEStore(what, arr, idx, field) { \
  double p[4];  \
  storeu(p, what); \
  arr[idx]field = p[0]; \
  arr[idx+1]field = p[1]; \
  arr[idx+2]field = p[2]; \
  arr[idx+3]field = p[3]; \
} 
enum { cosmoMask=0xf };

#endif /* COSMO_FLOAT */

#elif defined(__SSE2__)

#ifdef COSMO_FLOAT

#include "SSE-Float.h"
#define SSE_VECTOR_WIDTH 4
#define FORCE_INPUT_LIST_PAD 3
typedef SSEFloat SSEcosmoType; 
#define SSELoad(where, arr, idx, field) where(arr[idx]field, arr[idx+1]field, arr[idx+2]field, arr[idx+3]field)
#define SSEStore(what, arr, idx, field) { \
  float p[4]; \
  storeu(p, what); \
  arr[idx]field = p[0]; \
  arr[idx+1]field = p[1]; \
  arr[idx+2]field = p[2]; \
  arr[idx+3]field = p[3]; \
}
enum { cosmoMask=0xf };

#else

#include "SSE-Double.h"
#define SSE_VECTOR_WIDTH 2
#define FORCE_INPUT_LIST_PAD 1
typedef SSEDouble SSEcosmoType; 
#define SSELoad(where, arr, idx, field) where(arr[idx]field, arr[idx+1]field)
#define SSEStore(what, arr, idx, field) { \
  storel(&arr[idx]field, what); \
  storeh(&arr[idx+1]field, what); \
} 
enum { cosmoMask=0x3 };

#endif   /* end of COSMO_FLOAT under __SSE2__ */

#endif    /*  __SSE2__  */

#endif   /* end of __SSEDEFS_H__ */
