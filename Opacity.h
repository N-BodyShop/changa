#ifndef OPACITY_H
#define OPACITY_H
#define MAXOPACITIES 10
#include <pup_stl.h>
#include <string>
#include "GravityParticle.h"
#include "parameters.h"
#include "HydroParticle.h"

class Opacity {
public:
  // initialized opacities
  virtual bool isTemperatureDependent() = 0;
  virtual void initialize( ) = 0; 
  virtual void initializeParameters( Parameters param) = 0;
  virtual std::string name() = 0;

  virtual void computeOpacity( double rho, double T, varVector &prim, double &kappaPlanck, 
                               double &kappaMean, double &kappaSca) = 0;

  static int numOpacities;
  static Opacity *instances[MAXOPACITIES];// = new UserSources[UserSources::numSources];
  static bool isTempDependent;
  static void initOpacityParameters( Parameters param); 
  static void registerOpacity( std::string nameWanted, Opacity *opacity); 
  static bool isTemperatureDependentOpacity();
  static void computeTotalOpacity(double rho, double T, varVector &prim, double &kappaPlanck, 
                               double &kappaMean, double &kappaSca);
};



#endif
