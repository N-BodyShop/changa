#ifndef DUSTOPACITY_H
#define DUSTOPACITY_H
#include "Opacity.h"

class DustOpacity : public Opacity {
public:
  DustOpacity(){}
  bool isTemperatureDependent();
  void initialize( ); 
  void initializeParameters( Parameters param);
  std::string name();

  void computeOpacity( double rho, double T, varVector &prim, double &kappaPlanck, 
                       double &kappaMean, double &kappaSca);
};

#endif
