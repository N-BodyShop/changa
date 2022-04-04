#ifndef FIXEDBC_H
#define FIXEDBC_H
#include "UserSources.h"
#include "dVector.h"
//#include "param.h"
#define MAXRADBOUNDSOURCES 2


class FixedBC : public UserSources {
  dVector boundaries[6];
  double norm2[6];
  double rho0[6];
  double tau;
  int numBoundaries;
  bool useRadiation[6];

  // spherical boundary
  dVector center;
  double rhoSph, radius, tempSph;

#ifdef RADIATION
  dVector normals[MAXRADBOUNDSOURCES];
  double mu0, I0;
  int numRadSources;
#endif

public:
  FixedBC(){ 
#ifdef RADIATION
    numRadSources = MAXRADBOUNDSOURCES;
#endif 
  }
  void initialize(); 
  void initializeParameters( Parameters param);

  // apply the sources with a dt
  void applySource( double dt, dVector &pos, double volume, varVector &cons, varVector &deltas);

};
#endif
