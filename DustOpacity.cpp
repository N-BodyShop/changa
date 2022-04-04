#ifdef RADIATION
#include "DustOpacity.h"

void DustOpacity::initialize() {}

bool DustOpacity::isTemperatureDependent() { return true;}

void DustOpacity::initializeParameters( Parameters param){}

std::string DustOpacity::name() { return "DustOpacity";}

void DustOpacity::computeOpacity( double rho, double T, varVector &prim, double &kappaPlanck, 
                       double &kappaMean, double &kappaSca) {
  //double T2 = (fmin(T,100.)/10.); T2 *= T2; 
  double T2 = 0.01*T*T; 
  kappaPlanck = 0.1*T2;
  kappaMean = 0.0316*T2;
  kappaSca = 0.;                        
}
#endif //VORONOI
