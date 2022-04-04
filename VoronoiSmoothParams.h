#ifndef VORONOISMOOTHPARAMS_H

#define VORONOISMOOTHPARAMS_H

#include "GravityParticle.h"
#include "smooth.h"
#define VSP_ACTIVEONLY 1
#define VSP_NEIGHBORONLY 2
#define VSP_ACTIVE_AND_NEIGHBOR 3
#define VSP_EVERYTHING 0


/// @brief Parameters and functions for the first SPH smooth: density
/// and velocity derivatives.
class VoronoiSmoothParams : public SmoothParams {
protected:

  double a, H; // Cosmological parameters
  int    bActiveOnly;
  bool   bUseDt;
//  double dt;
  int    hydroStep;
  bool   bReSearch;
  bool   bVoronoiICs;
  bool   bProblemGenerator;

  virtual void fcnSmooth(GravityParticle *p,
                         int              nSmooth,
                         pqSmoothNode    *nList);
  virtual int  isSmoothActive(GravityParticle *p);
  virtual void initTreeParticle(GravityParticle *p);
  virtual void postTreeParticle(GravityParticle *p) {}

  virtual void initSmoothParticle(GravityParticle *p);
  virtual void initSmoothCache(GravityParticle *p);
  virtual void combSmoothCache(GravityParticle        *p1,
                               ExternalSmoothParticle *p2);

public:

  VoronoiSmoothParams() {}

  /// @param _iType Type of particle to operate on
  /// @param am Active rung
  /// @param csm Cosmology information
  /// @param dTime Current time
  /// @param _bActiveOnly Only operate on active particles.
  /// @param _iType Type of particle to operate on
  /// @param am Active rung
  /// @param csm Cosmology information
  /// @param dTime Current time
  /// @param _bActiveOnly Only operate on active particles.

  VoronoiSmoothParams(int _iType, int am, CSM csm, double dTime,
                      int _bActiveOnly, double _dt, int _hydroStep, 
                      bool _bReSearch=false) {

    iType       = _iType;
    activeRung  = am;
  //  dt          = _dt;
    bActiveOnly = _bActiveOnly;
    bReSearch = _bReSearch;
    bVoronoiICs = false;
    bProblemGenerator = false;
    hydroStep   = _hydroStep;
  
    //if(bCheckSolution) { 
    //  initialize = false; bSolveRiemann = false;
    //}

    if (csm->bComove) {
      H = csmTime2Hub(csm, dTime);
      a = csmTime2Exp(csm, dTime);
    }
    else {
      H = 0.0;
      a = 1.0;
    }

    if ( hydroStep == INITIALIZE) ckout << " We are initializing\n ";
  }

  void useVoronoiICs() {
    bVoronoiICs = true; 
  }

  void useProblemGenerator() {
    bProblemGenerator = true;
  }
  
  static const int INITIALIZE = 0;
  static const int BUILD_MESH = 1;
  static const int HALF_STEP  = 2;
  static const int FULL_STEP  = 3;
  static const int CHECK_SOLUTION  = 4;
  static const int MAP_A_TO_B = 5;

  PUPable_decl(VoronoiSmoothParams);
  VoronoiSmoothParams(CkMigrateMessage *m) : SmoothParams(m) {}

  virtual void pup(PUP::er& p) {
    SmoothParams::pup(p); // Call base class
    p | a;
    p | H;
    p | bActiveOnly;
  //  p | dt;
    p | hydroStep;
    p | bReSearch;
    p | bVoronoiICs;
    p | bProblemGenerator;
  }
};

#endif // ifndef VORONOISMOOTHPARAMS_H
