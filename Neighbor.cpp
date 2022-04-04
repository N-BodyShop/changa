#ifdef VORONOI
#include "Neighbor.h"
#include <algorithm>

void VoronoiCompute::compute(HydroParticle& hp) {
  std::sort(neighborList.begin(), neighborList.end());

  double mrs = 0.; int i = 0; int lastP=0;
  for (std::vector<Neighbor>::iterator iter = neighborList.begin();
       iter != neighborList.end(); iter++, i++) {
    cell.nplane(iter->x, iter->y, iter->z, iter->id);

    // compute up to here
    mrs = iter->r2;
    if ((i % 5 == 0) &&  ((cell.max_radius_squared()) < mrs)) {
        break;
    }
  }
  
  bool validCell = (cell.max_radius_squared() < msr*msr);

  hp.setFacesAndNeighbors(cell, validCell);
}
#endif //VORONOI
