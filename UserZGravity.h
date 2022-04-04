#ifndef USERZGRAVITY_H
#define USERZGRAVITY_H
#include <vector>
#include "UserGravity.h"

class UserZGravity : public UserGravity {
public:
  double zGravity, gravLengthScale;
  int direction;
  bool bGravityFixedAtOne;

  UserZGravity() {}
  void initialize( ){ direction=1; bGravityFixedAtOne = true;}
  void initializeParameters( Parameters param);

  void getGravity( double t, double x, double y, double z, double& gx, double &gy, double &gz, double &dtg);
  int gravType() { return Z_GRAVITY;}
};
#endif
