#ifndef NEIGHBOR_H
#define NEIGHBOR_H
#include <vector>
#include "voro++.hh"
#include "HydroParticle.h"

class Neighbor {
public:

  double x, y, z;
  int    id;
  double r2;

  Neighbor() {}

  Neighbor(int identifier, double xp, double yp, double zp) {
    id = identifier;
    x  = xp; y = yp; z = zp;
    r2 = x * x + y * y + z * z;
  }

  bool operator<(const Neighbor& rhs) const {
    return r2 < rhs.r2;
  }
};

class VoronoiCompute {
  std::vector<Neighbor> neighborList;
  voro::voronoicell_neighbor cell;
  double x, y, z;
  static const int NumNeighbors = 512;
  double msr;
public:

  VoronoiCompute(double xp, double yp, double zp, double r, double maxSearchRadius) {
    VoronoiCompute( xp, yp, zp, r, r, r, maxSearchRadius);
  }

  VoronoiCompute(double xp, double yp, double zp, double rx, double ry, double rz, double maxSearchRadius) {
    x = xp; y = yp; z = zp;
    msr = maxSearchRadius;
    double xmin = x - rx, xmax = x + rx, ymin = y - ry, ymax = y + ry, zmin = z - rz,
           zmax = z + rz;
    cell.init(xmin, xmax, ymin, ymax, zmin, zmax);
    neighborList.reserve(NumNeighbors);
  }

  void addNeighbor(int id, double xn, double yn, double zn) {
    Neighbor n(id, xn - x, yn - y, zn - z);

    neighborList.push_back(n);
  }

  void compute(HydroParticle& hp);
};
#endif // ifndef NEIGHBOR_H
