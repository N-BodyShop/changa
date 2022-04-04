#ifndef USERSOURCES_H
#define USERSOURCES_H
#define MAXUSERSOURCES 10
//#include "HydroParticle.h"
#include <pup_stl.h>
#include "GravityParticle.h"
#include "parameters.h"

class UserSources {
public:

  // initialized a source
  virtual void initialize( ) = 0; 
  virtual void initializeParameters( Parameters param) = 0;

  // apply the sources with a dt
  virtual void applySource( double dt, dVector &pos, double volume, varVector &cons, varVector &deltas) = 0;

  static int numSources;
  static UserSources *instances[MAXUSERSOURCES];// = new UserSources[UserSources::numSources];

  static void initSourceParameters( Parameters param); 
  static void registerSource( UserSources *userSources); 
};
#endif
