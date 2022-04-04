#ifndef SIMPLEOPACITY_H
#define SIMPLEOPACITY_H
#include "Opacity.h"

class SimpleOpacity : public Opacity {
public:
  double SimpleKappaPl, SimpleKappaMean, SimpleKappaSca; 
  SimpleOpacity(){}
  bool isTemperatureDependent();
  void initialize( ); 
  void initializeParameters( Parameters param);
  std::string name();

  void computeOpacity( double rho, double T, varVector &prim, double &kappaPlanck, 
                       double &kappaMean, double &kappaSca);
};

#endif
