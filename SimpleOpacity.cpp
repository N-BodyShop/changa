#ifdef RADIATION
#include "SimpleOpacity.h"

void SimpleOpacity::initialize() {}

bool SimpleOpacity::isTemperatureDependent() { return false;}

std::string SimpleOpacity::name() { return "SimpleOpacity";}

void SimpleOpacity::initializeParameters( Parameters param){
  SimpleKappaPl = param.dSimpleKappaPl;
  SimpleKappaMean = param.dSimpleKappaMean;
  SimpleKappaSca = param.dSimpleKappaSca;
}

void SimpleOpacity::computeOpacity( double rho, double T, varVector &prim, double &kappaPlanck, 
                       double &kappaMean, double &kappaSca) {
  kappaPlanck = SimpleKappaPl;
  kappaMean = SimpleKappaMean;
  kappaSca =SimpleKappaSca;                        
}
#endif //RADIATION
