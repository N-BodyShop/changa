#ifdef VORONOI
#include "HydroParticle.h"

void HydroUtils::MapToCoolParticle( varVector &prim, COOLPARTICLE &cp, COOL *cl ) {
#ifdef COOLING_SROEOS
  cp.ye = fmax(0.1,prim[iye]);
#endif
#ifdef COOLING_MESA
//  cp.X = prim[iX];
  cp.X = cl->X; // cp.X = 0.0; // 0.7;
  cp.Z = cl->Z; // cp.Z = 1.0; // 0.02;
#endif
} 
#endif //VORONOI
