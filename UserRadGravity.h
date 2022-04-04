#ifndef USERRADGRAVITY_H
#define USERRADGRAVITY_H
#include <vector>
#include "UserGravity.h"

class UserRadGravity : public UserGravity {
  std::vector<double> r, m, g;
  char filename[256];
public:
  UserRadGravity() {}
  
  void initialize( );
  void initializeParameters( Parameters param);

  void getGravity( double t, double x, double y, double z, double& gx, double &gy, double &gz, double &dtg);
  int gravType() { return RAD_GRAVITY;}

};
#endif
