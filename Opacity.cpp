#ifdef RADIATION
#include "Opacity.h"

int Opacity::numOpacities = 0;
Opacity *Opacity::instances[MAXOPACITIES] = {0};
bool Opacity::isTempDependent = false;

void Opacity::initOpacityParameters( Parameters param) {
  for( int i = 0; i < numOpacities; i++) {
    instances[i] -> initializeParameters( param);
    instances[i] -> initialize();
  }
}

void Opacity::registerOpacity( std::string nameWanted, Opacity *opacity){
  if( nameWanted.compare(opacity->name()) == 0) {
    Opacity::instances[Opacity::numOpacities++] = opacity; 
    isTempDependent = isTempDependent || opacity->isTemperatureDependent();
  }
}

bool Opacity::isTemperatureDependentOpacity() { return isTempDependent;}

void Opacity::computeTotalOpacity( double rho, double T, varVector &prim, double &kappaPlanck, 
                               double &kappaMean, double &kappaSca) {
  kappaPlanck = 0.; kappaMean = 0.; kappaSca = 0.;
  for(int i = 0; i < numOpacities; i++) {
    double tmpKappaPlanck = 0., tmpKappaMean = 0., tmpKappaSca = 0.;
    instances[i]->computeOpacity( rho, T, prim, tmpKappaPlanck, tmpKappaMean, tmpKappaSca);
    kappaPlanck += tmpKappaPlanck; kappaMean += tmpKappaMean; kappaSca += tmpKappaSca;
  }
}
#endif //VORONOI
