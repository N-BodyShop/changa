#ifdef VORONOI
#include "UserProblemGenerator.h"
#include <cstdlib>
//#include <math>
#include <iostream>
#include "physconst.h"


UserProblemGenerator *UserProblemGenerator::instance = 0;

void UserProblemGenerator::initUserProblemGeneratorParameters( Parameters param) {
  instance -> initializeParameters( param);
  instance -> initialize();
}

void UserProblemGenerator::registerUserProblemGenerator( std::string name, UserProblemGenerator *userGenerator){
  if( userGenerator -> name().compare( name) == 0) {
    UserProblemGenerator::instance = userGenerator; 
  }
}

void WaveProblemGenerator::initializeParameters( Parameters param) {
  eigenVectorReal = 0.; eigenVectorImag = 0.; 
  dens              = param.dWaveRho/param.dGmPerCcUnit;
  // isothermal sound speed
  cs2               = param.dWaveCs*param.dWaveCs/param.dErgPerGmUnit;
  amplitude         = param.dWaveAmplitude;
  kx                = param.dWaveKxBox/param.vPeriod.x*2.*M_PI;
  eigenVectorReal[iRho] = param.dWaveRhoReal;
  eigenVectorImag[iRho] = param.dWaveRhoImag;  
  eigenVectorReal[ivx]  = param.dWaveVxReal;
  eigenVectorImag[ivx]  = param.dWaveVxImag;  
  eigenVectorReal[iE]   = param.dWaveIeReal;
  eigenVectorImag[iE]   = param.dWaveIeImag;
#ifdef RADIATION
  ErReal                = param.dWaveErReal;
  ErImag                = param.dWaveErImag;
  FrReal                = param.dWaveFrReal;
  FrImag                = param.dWaveFrImag; 
  // Compute T2 
  const double T = param.dWaveCs*param.dWaveCs/RGAS*HydroUtils::defaultMu;
  const double T4 = T*T*T*T;
  Iinit = ARAD*T4/(param.dErgPerGmUnit*param.dGmPerCcUnit)/(4.*M_PI);

#endif

}

void WaveProblemGenerator::applyGenerator( dVector &pos, varVector &prim) {

  const double fourPi = 4.*M_PI;
  // set background
  prim       = 0.;
  prim[iRho] = dens;
  prim[ivx]  = 0.;
  prim[iE]   = cs2/(HydroUtils::defaultGamma - 1.);

#ifdef RADIATION

  for( int ivar = IRADVARSTART; ivar < IRADVAREND; ivar++) {
	  prim[ivar] = Iinit; //HydroUtils::instance->initRadIntensity; 
  }

  // compute wmu * mu*mu
  double mu2norm = 0.;
  for( int iAng = 0; iAng < NANGVARS; iAng++) {
    int i =  iAng;
    double y2norm = HydroUtils::instance->radNormals[i][0]; 
    y2norm*=y2norm*HydroUtils::instance->radWeights[i];
    mu2norm += y2norm;
  }

#endif

  // Now set the perturbations
  prim[iRho] *= 1. + amplitude*(eigenVectorReal[iRho]*cos( kx*pos[0]) + eigenVectorImag[iRho]*sin(kx*pos[0]));
#ifdef RADIATION
  prim[iE] *= 1. + amplitude*(eigenVectorReal[iE]*cos( kx*pos[0]) + eigenVectorImag[iE]*sin(kx*pos[0]));
  prim[ivx] = sqrt(cs2)*amplitude*(eigenVectorReal[ivx]*cos( kx*pos[0]) + eigenVectorImag[ivx]*sin(kx*pos[0]));
#else
  prim[iE] *= 1. + amplitude*(HydroUtils::defaultGamma-1.)*(eigenVectorReal[iE]*cos( kx*pos[0]) + eigenVectorImag[iE]*sin(kx*pos[0]));
  prim[ivx] = sqrt(HydroUtils::defaultGamma*cs2)*amplitude*(eigenVectorReal[ivx]*cos( kx*pos[0]) + eigenVectorImag[ivx]*sin(kx*pos[0]));
#endif 

#ifdef RADIATION
  for( int iFreq = 0; iFreq < NFREQVARS; iFreq++) {
    double dEr = amplitude*(ErReal*cos(kx*pos[0]) + ErImag*sin(kx*pos[0]));
    double dFr = amplitude*(FrReal*cos(kx*pos[0]) + FrImag*sin(kx*pos[0]));
    for( int iAng = 0; iAng < NANGVARS; iAng++) {
      int i =  iAng + iFreq*NANGVARS;
      prim[i + IRADVARSTART] *= 1. + dEr + dFr*HydroUtils::instance->radNormals[i][0]*HydroUtils::instance->radWeights[i]/mu2norm;
    }
  }
#endif

}


void DynDiffusionProblemGenerator::initializeParameters( Parameters param) {
  dens              = param.dWaveRho/param.dGmPerCcUnit;
  // isothermal sound speed
  cs2               = param.dWaveCs*param.dWaveCs/param.dErgPerGmUnit;
  vx                = param.dvx*param.dSecUnit/(3.0857e21*param.dKpcUnit);
#ifdef RADIATION
  DeltaX            = param.dDeltaX;
  // Compute T2 
  const double T = param.dRadTeff;
  const double T4 = T*T*T*T;
  Iinit = ARAD*T4/(param.dErgPerGmUnit*param.dGmPerCcUnit)/(4.*M_PI);
#endif

}

void DynDiffusionProblemGenerator::applyGenerator( dVector &pos, varVector &prim) {

  const double fourPi = 4.*M_PI;
  // set background
  prim       = 0.;
  prim[iRho] = dens;
  prim[ivx]  = vx;
  prim[iE]   = cs2/(HydroUtils::defaultGamma - 1.);

#ifdef RADIATION

  // compute wmu * mu*mu
  double mu2norm = 0.;
  for( int iAng = 0; iAng < NANGVARS; iAng++) {
    int i =  iAng;
    double y2norm = HydroUtils::instance->radNormals[i][0]; 
    y2norm*=y2norm*HydroUtils::instance->radWeights[i];
    mu2norm += y2norm;
  }
  const double x = pos[0];
  const double Er = Iinit*exp(-fmin(x*x/(DeltaX*DeltaX), 10));
  const double Fr = Er*4.*vx/HydroUtils::instance->cInCodeUnits;
  for( int iFreq = 0; iFreq < NFREQVARS; iFreq++) {
    for( int iAng = 0; iAng < NANGVARS; iAng++) {
      int i =  iAng + iFreq*NANGVARS;
      prim[i + IRADVARSTART] = Er + Fr*HydroUtils::instance->radNormals[i][0]*HydroUtils::instance->radWeights[i]/mu2norm;
    }
  }
#endif
}
#endif //VORONOI