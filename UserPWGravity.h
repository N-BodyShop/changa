#ifndef USERPWGRAVITY_H
#define USERPWGRAVITY_H
#include <vector>
#include "UserGravity.h"

class UserPWGravity : public UserGravity {
public:
  double Mbh, rsoft, cspeed;
  double xbh, ybh, zbh;
  double vxbh, vybh, vzbh;
  UserPWGravity() {}
  void initialize( ){ }
  void initializeParameters( Parameters param);

  void getGravity( double t, double x, double y, double z, double& gx, double &gy, double &gz, double &dtg);
  int gravType() { return PW_GRAVITY;}
};
#endif
