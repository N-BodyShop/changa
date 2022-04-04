#ifdef VORONOI
#include "UserSources.h"
#include <cstdlib>

int UserSources::numSources = 0;
UserSources *UserSources::instances[MAXUSERSOURCES] = {0};

void UserSources::initSourceParameters( Parameters param) {
  for( int i = 0; i < numSources; i++) {
    instances[i] -> initializeParameters( param);
    instances[i] -> initialize();
  }
}

void UserSources::registerSource( UserSources *userSource){
  UserSources::instances[UserSources::numSources++] = userSource; 
}
#endif //VORONOI
