#ifndef __COSMOTYPE_H__
#define __COSMOTYPE_H__

#define COSMO_FLOAT
#ifdef COSMO_FLOAT
#define CONVERT_TO_COSMO_TYPE (float)
typedef float cosmoType;
#define COSMO_CONST(val) val##f
#else
#define CONVERT_TO_COSMO_TYPE
typedef double cosmoType;
#define COSMO_CONST(val) val
#endif

#endif
