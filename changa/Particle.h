#ifndef __PARTICLE_H__
#define __PARTICLE_H__

#include "Vector3D.h"
#include "SFC.h"

#ifdef STORED_PRECISION_SINGLE
typedef float storedType;
#else
typedef double storedType;
#endif

class ExternalGravityParticle {
 private:
  storedType _soft;
  storedType _mass;
  Vector3D<storedType> _position;

 public:
  inline storedType& soft() {return _soft;}
  inline storedType& mass() {return _mass;}
  inline Vector3D<storedType>& position() {return _position;}
};

/*
class ExternalDensityParticle {
  storedType mass;
  Vector3D<storedType> position;
};
*/
typedef ExternalGravityParticle ExternalDensityParticle;

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

class DensityParticle {
 private:
  storedType _fDensity;
  unsigned int _iType;	// Bitmask to hold particle type information

 public:
  inline storedType& fDensity() {return _fDensity;}
  inline unsigned int& iType() {return _iType;}
};

class PressureParticle {
 private:
  Vector3D<storedType> _treeAcceleration;

 public:
  inline Vector3D<storedType>& treeAcceleration() {return _treeAcceleration;}
};

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
