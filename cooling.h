#ifndef COOLING_HINCLUDED
#define COOLING_HINCLUDED

/*
 * Cooling includes initially from GASOLINE
 */

#ifdef COOLING_NONE
#include "param.h"

typedef struct CoolingParametersStruct {
	double    dParam1;
	} COOLPARAM;

typedef struct CoolingPKDStruct COOL;

struct CoolingPKDStruct { 
   double     dTime;
};

inline COOL *CoolInit() {return NULL ;}
inline void CoolFinalize( COOL *cl ){};
inline void CoolAddParams( COOLPARAM *CoolParam, PRM ) {};
#else

#ifdef COOLING_DISK
#include "cooling_disk.h"
#else

#ifdef COOLING_PLANET
#include "cooling_planet.h"
#else

#ifdef COOLING_COSMO
#include "cooling_cosmo.h"
#else

#ifdef COOLING_METAL
#include "cooling_metal.h"
#else

#ifdef COOLING_BATE
#include "cooling_bate.h"
#else

#error "No valid cooling function specified"

#endif
#endif
#endif
#endif
#endif

#endif

#endif
