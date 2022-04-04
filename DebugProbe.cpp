#ifdef VORONOI
#include "DebugProbe.h"
#include "HydroParticle.h"

void DebugProbe::listInfo( GravityParticle *p) {
  double g2 = p->treeAcceleration.lengthSquared();
  double g = sqrt(g2);

#ifdef VORONOI
  dVector grav = p->hydroParticle().getParticleAcceleration();
  dVector tgrav( p->treeAcceleration.x, p->treeAcceleration.y, p->treeAcceleration.z);

  dVector delta = grav-tgrav;
  //printf( "%d %5.3e %5.3e %5.3e\n", p->iOrder, g, length(delta), p->mass);
#endif

  //printf( "%d %5.3e %5.3e\n", p->iOrder, g, p->mass);
} 

#endif //VORONOI
