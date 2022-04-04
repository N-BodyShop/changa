#ifndef USERRADIATION_SOURCES_H
#define USERRADIATION_SOURCES_H
#include "UserSources.h"
#include "dVector.h"
//#include "param.h"
#define MAXRADSOURCES 2

class UserRadiationSources : public UserSources {
  double I0, F0;
  double radiusSource, mu0;
  double y0, deltaY;
  bool pointSource, slabSource, constantFlux;
  dVector normals[MAXRADSOURCES], position[MAXRADSOURCES];
  double tau;
  int numRadSources;

public:
  UserRadiationSources(){ numRadSources = MAXRADSOURCES; constantFlux = true;};

  void initialize(); 
  void initializeParameters( Parameters param);

  // apply the sources with a dt
  void applySource( double dt, dVector &pos, double volume, varVector &cons, varVector &deltas);

};
#endif
