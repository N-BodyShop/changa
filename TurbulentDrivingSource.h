#ifndef TURB_SOURCES_H
#define TURB_SOURCES_H
#include "UserGravity.h"
#include "dVector.h"
//#include "param.h"
#define TURB_NX 100
#define TURB_NY 100
#define TURB_NZ 100
#define NSTIR 4
#define TWOPI 6.282

class TurbulentDrivingSource : public UserGravity {
  double klo, khi, kslope, amplitude;
  int seed;
  dVector turbAcc[TURB_NX][TURB_NY][TURB_NZ];

  void interpolate( dVector &pos, dVector &acc);
public:
  TurbulentDrivingSource(){klo = 0.; khi = 0.; kslope = 0.; amplitude = 0.;};
  void initializeParameters( Parameters param);
  void initialize();
  void getGravity( double t,double x, double y, double z, double& gx, double& gy, double& gz, double &dtg);
  int gravType() { return TURB_STIR;}
};
#endif
