#ifndef __PARTICLE_H__
#define __PARTICLE_H__

#include "Vector3D.h"
#include "SFC.h"

/*
 *  XXX this doesn't belong here.  It should be moved to a more
 *  appropriate place.
 */
class TreeWalk;
class Compute;
class Opt;
class State;

// Object to record a type of active walk. Contains pointers to TreeWalk/Compute/Opt (T/C/O) combinations
class ActiveWalk {
  public:
  TreeWalk *tw;
  Compute *c;
  Opt *o;
  State *s;
  
  ActiveWalk(TreeWalk *_tw, Compute *_c, Opt *_o, State *state) : 
      tw(_tw), c(_c), o(_o), s(state){}
  ActiveWalk(){}
};

// Object to bookkeep a Bucket Walk.
class BucketGravityRequest {
public:
		
	unsigned int numAdditionalRequests;
	int finished;

	BucketGravityRequest(unsigned int bucketSize = 0) :
	  numAdditionalRequests(0), finished(0) {
	}
	
};

/*
 * XXX The above doesn't belong here
 */

#ifdef STORED_PRECISION_SINGLE
typedef float storedType;
#else
typedef double storedType;
#endif

/* 
 * class needed to send gravity information to another processor
 */
class ExternalGravityParticle {
 private:
  storedType _soft;
  storedType _mass;
  Vector3D<storedType> _position;

 public:
  inline storedType& soft() {return _soft;}
  inline storedType& mass() {return _mass;}
  inline Vector3D<storedType>& position() {return _position;}
  void pup(PUP::er &p) {
    p | _position;
    p | _mass;
    p | _soft;
  }
};

/*
class ExternalDensityParticle {
  storedType mass;
  Vector3D<storedType> position;
};
*/
typedef ExternalGravityParticle ExternalDensityParticle;

/* 
 * class needed to send information for pressure calculation
 */
class ExternalPressureParticle {
 private:
  storedType _mass;
  Vector3D<storedType> _position;
  storedType _fBall;
  storedType _fDensity;
  unsigned int _iType;	// Bitmask to hold particle type information

 public:
  inline storedType& soft() {return _soft;}
  inline storedType& mass() {return _mass;}
  inline Vector3D<storedType>& position() {return _position;}
  inline storedType& fBall() {return _fBall;}
  inline storedType& fDensity() {return _fDensity;}
  inline unsigned int& iType() {return _iType;}
};

/*
 * This might be needed for a "send" version of gravity
 */
class GravityParticle {
 private:
  Vector3D<storedType> _treeAcceleration;
  storedType _potential;
  storedType _dtGrav;

 public:
  inline Vector3D<storedType>& treeAcceleration() {return _treeAcceleration;}
  inline storedType& potential() {return _potential;}
  inline storedType& dtGrav() {return _dtGrav;}
};

/*
 * Writeback information for symmetric Density
 */
class DensityParticle {
 private:
  storedType _fDensity;
  unsigned int _iType;	// Bitmask to hold particle type information

 public:
  inline storedType& fDensity() {return _fDensity;}
  inline unsigned int& iType() {return _iType;}
};

/*
 * Writeback information for symmetric Pressure
 */
class PressureParticle {
 private:
  Vector3D<storedType> _treeAcceleration;

 public:
  inline Vector3D<storedType>& treeAcceleration() {return _treeAcceleration;}
};

/*
 * Base particle class
 */
class Particle {
 private:
  storedType _soft;
  storedType _mass;
  Vector3D<storedType> _position;
  storedType _fBall;
  storedType _fDensity;
  unsigned int _iType;	// Bitmask to hold particle type information
  int _rung;  ///< the current rung (greater means faster)
  Vector3D<storedType> _treeAcceleration;
  storedType _potential;
  storedType _dtGrav;
  Vector3D<storedType> _velocity;
  SFC::Key _key;
  int _iOrder;                /* input order of particles */
  unsigned int _extraDataIdx;

 public:
  inline storedType& soft() {return _soft;}
  inline storedType& mass() {return _mass;}
  inline Vector3D<storedType>& position() {return _position;}
  inline storedType& fBall() {return _fBall;}
  inline storedType& fDensity() {return _fDensity;}
  inline unsigned int& iType() {return _iType;}
  inline int& rung() {return _rung;}
  inline Vector3D<storedType>& treeAcceleration() {return _treeAcceleration;}
  inline storedType& potential() {return _potential;}
  inline storedType& dtGrav() {return _dtGrav;}
  inline Vector3D<storedType>& velocity() {return _velocity;}
  inline SFC::Key& key() {return _key;}
  inline int& iOrder() {return _iOrder;}
  inline unsigned int& extraDataIdx() {return _extraDataIdx;}

  inline operator ExternalGravityParticle&() {return *(ExternalGravityParticle*)this;}
  //inline operator ExternalDensityParticle&() {return *(ExternalDensityParticle*)&this->_mass;}
  inline operator ExternalPressureParticle&() {return *(ExternalPressureParticle*)&this->_mass;}
  inline operator GravityParticle&() {return *(GravityParticle*)&this->_treeAcceleration;}
  inline operator DensityParticle&() {return *(DensityParticle*)&this->_fDensity;}
  inline operator PressureParticle&() {return *(PressureParticle*)&this->_treeAcceleration;}

  inline ExternalGravityParticle* getExternalGravityParticle() {return (ExternalGravityParticle*)this;}
  //inline ExternalDensityParticle* getExternalDensityParticle() {return (ExternalDensityParticle*)&this->_mass;}
  inline ExternalPressureParticle* getExternalPressureParticle() {return (ExternalPressureParticle*)&this->_mass;}
  inline GravityParticle* getGravityParticle() {return (GravityParticle*)&this->_treeAcceleration;}
  inline DensityParticle* getDensityParticle() {return (DensityParticle*)&this->_fDensity;}
  inline PressureParticle* getPressreParticle() {return (PressureParticle*)&this->_treeAcceleration;}
};


/*
 * Class to send information for a "smooth + pressure" calculation
 */
class ExternalSmoothPressureParticle {
 private:
  storedType _PoverRho2;          /* P/rho^2 */
  Vector3D<storedType> _vPred;    /* predicted velocity (time centered) */
  storedType _c;                  /* sound speed */
  storedType _mumax;              /* sound speed like viscosity term */

 public:
  inline storedType& PoverRho2() {return _PoverRho2;}
  inline Vector3D<storedType>& vPred() {return _vPred;}
  inline storedType& c() {return _c;}
  inline storedType& mumax() {return _mumax;}
};

/*
 * Writeback for "smooth + pressure" calculation
 */
class SmoothPressureParticle {
 private:
  storedType _mumax;              /* sound speed like viscosity term */
  storedType _PdV;                /* P dV heating (includes shocking) */
  storedType _divv;
  Vector3D<storedType> _curlv;

 public:
  inline storedType& mumax() {return _mumax;}
  inline storedType& PdV() {return _PdV;}
  inline storedType& divv() {return _divv;}
  inline Vector3D<storedType>& curlv() {return _curlv;}
};

/*
 * "extra data" class for SPH
 */
class SmoothParticle {
 private:
  storedType PoverRho2;          /* P/rho^2 */
  Vector3D<storedType> vPred;    /* predicted velocity (time centered) */
  storedType c;                  /* sound speed */
  storedType mumax;              /* sound speed like viscosity term */
  storedType PdV;                /* P dV heating (includes shocking) */
  storedType divv;
  Vector3D<storedType> curlv;
  storedType u;                  /* thermal energy */
  storedType uPred;              /* predicted thermal energy */
  storedType BalsaraSwitch;      /* Balsara viscosity reduction */
  storedType fMetals;

 public:
  inline storedType& PoverRho2() {return _PoverRho2;}
  inline Vector3D<storedType>& vPred() {return _vPred;}
  inline storedType& c() {return _c;}
  inline storedType& mumax() {return _mumax;}
  inline storedType& PdV() {return _PdV;}
  inline storedType& divv() {return _divv;}
  inline Vector3D<storedType>& curlv() {return _curlv;}
  inline storedType& u() {return _u;}
  inline storedType& uPred() {return _uPred;}
  inline storedType& BalsaraSwitch() {return _BalsaraSwitch;}
  inline storedType& fMetals() {return _fMetals;}

  inline operator ExternalSmoothPressureParticle&() {return *(ExternalSmoothPressureParticle*)this;}
  inline operator SmoothPressureParticle&() {return *(SmoothPressureParticle*)&this->mumax;}

  inline ExternalSmoothPressureParticle* getExternalSmoothPressureParticle() {return (ExternalSmoothPressureParticle*)this;}
  inline SmoothPressureParticle* getSmoothPressureParticle() {return (SmoothPressureParticle*)&this->mumax;}
};

#endif
