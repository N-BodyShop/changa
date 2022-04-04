#ifndef RADIATION_SOURCES_H
#define RADIATION_SOURCES_H
#include "UserSources.h"
#include "dVector.h"
#include "HydroParticle.h"
#include "physconst.h"

//#include "param.h"

class RadiationSource : public UserSources {

private:
  int FourthPolyRoot(const double coef4, const double tconst, double &root);
  bool noLorentz, noHydro;
  double computeGamma( dVector v, dVector normal);
public:
  RadiationSource(){};

  void initialize(); 
  void initializeParameters( Parameters param);

  // apply the sources with a dt
  void applySource( double dt, 
                    dVector &pos, 
                    double volume, 
                    varVector &cons, 
                    varVector &deltas);

  double opacityConversion, ergsccConversion, densityConversion;
  inline double emissivity( double T) {return ARAD*T*T*T*T/(4.*M_PI)*ergsccConversion;} 
  //void computeJ( varVector &prim, double &J);
  void computeNewTgas( double dt,
                       double &J, 
                       double &Tgas, 
                       double &DeltaE,
                       varVector &prim,
                       double radWeightsCM[]);
  void labToCom( varVector &prim, varVector &primCM, double radWeights[]);
  void comToLab( varVector &primCM, varVector &prim);

  static RadiationSource* instance;
};
#endif
