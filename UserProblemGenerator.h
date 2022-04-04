#ifndef USERPROBLEMGENERATOR_H
#define USERPROBLEMGENERATOR_H
#include "HydroParticle.h"
#include <pup_stl.h>
#include "GravityParticle.h"
#include "parameters.h"
#include <string>

class UserProblemGenerator {
public:

  // initialized a source
  virtual void initialize( ) = 0; 
  virtual void initializeParameters( Parameters param) = 0;
  virtual std::string name() = 0;

  // apply the sources with a dt
  virtual void applyGenerator( dVector &pos, varVector &prim) = 0;

  static UserProblemGenerator *instance;

  static void initUserProblemGeneratorParameters( Parameters param); 
  static void registerUserProblemGenerator( std::string name, UserProblemGenerator *userGenerator); 
};

class WaveProblemGenerator : public UserProblemGenerator {
private:
  varVector eigenVectorReal, eigenVectorImag;
#ifdef RADIATION  
  double ErReal, ErImag, FrReal, FrImag;
  double Iinit;
#endif
  double dens, cs2; 
  double amplitude, kx;
public:
  void initialize() {}
  void initializeParameters( Parameters param);
  std::string name() { std::string problemName("Wave"); return problemName;}
  void applyGenerator( dVector &pos, varVector &prim);
};

class DynDiffusionProblemGenerator : public UserProblemGenerator {
private:
#ifdef RADIATION  
  double Iinit, DeltaX;
#endif
  double dens, cs2, vx; 
public:
  void initialize() {}
  void initializeParameters( Parameters param);
  std::string name() { std::string problemName("DynamicDiffusion"); return problemName;}
  void applyGenerator( dVector &pos, varVector &prim);
};
#endif

