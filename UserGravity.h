#ifndef USERGRAVITY_H
#define USERGRAVITY_H
#include <pup_stl.h>
#include <string>
#include "GravityParticle.h"
#include "parameters.h"
#define RAD_GRAVITY 1
#define TURB_STIR 2
#define Z_GRAVITY 3
#define PW_GRAVITY 4
#define MAX_INSTANCES 10

class UserGravity {
public:
  static double secUnit;
  static double cmUnit;
  static double xPeriod, yPeriod, zPeriod;

  virtual void getGravity( double t, double x, double y, double z, double& gx, double &gy, double &gz, double &dtg) = 0;
  virtual int gravType() = 0;

  static UserGravity *instance;
  static UserGravity *gravInstances[MAX_INSTANCES];
  static int numInstances;
  static int type;

  // initialized a source
  virtual void initialize( ) = 0;
  virtual void initializeParameters( Parameters param) = 0;

  static void initConstants( int type, double sec_unit, double cm_unit, double xp, double yp, double zp);

  static void initUserGravityParameters( Parameters param);
  static void initUserGravity();
};
#endif
