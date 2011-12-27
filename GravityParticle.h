#ifndef GRAVITYPARTICLE_H
#define GRAVITYPARTICLE_H

#include "cooling.h"

#include "SFC.h"
#include <vector>

// Object to bookkeep a Bucket Walk.
class BucketGravityRequest {
public:
		
	int finished;

	BucketGravityRequest(unsigned int bucketSize = 0) : finished(0) {
	}
	
};

// Information needed to calculate gravity

class ExternalGravityParticle {
 public:

  double mass;
  double soft;
  Vector3D<double> position;

  void pup(PUP::er &p) {
    p | position;
    p | mass;
    p | soft;
  }
};

class GravityParticle;

class PushExternalGravityParticle : public ExternalGravityParticle {
  public:

  int rung;
  double dtGrav;
#ifdef CHANGESOFT
  double fSoft0;
#endif
#ifdef NEED_DT
  double dt;
#endif

  void pup(PUP::er &p){
    ExternalGravityParticle::pup(p);
    p|rung;
    p|dtGrav;
#ifdef CHANGESOFT
    p|fSoft0;
#endif
#ifdef NEED_DT
    p|dt;
#endif
  }

  PushExternalGravityParticle &operator=(const GravityParticle &p);
};

// Extra data needed for SPH
class extraSPHData 
{
 private:
    double _u;			/* Internal Energy */
    double _fMetals;		/* Metalicity */
    Vector3D<double> _vPred;	/* Predicted velocities for velocity
				   dependent forces */
    double _uPred;		/* Predicted internal energy */
    double _divv;		/* Diverence of the velocity */
    Vector3D<double> _curlv;	/* Curl of the velocity */
    double _mumax;		/* */
    double _PdV;
    double _c;			/* Speed of Sound */
    double _PoverRho2;		/* Pressure/rho^2 */
    double _BalsaraSwitch;	/* Pressure/rho^2 */
    double _fBallMax;		/* Radius for inverse neighbor finding */
#ifndef COOLING_NONE
    double _uDot;		/* Rate of change of u, for
				   predicting u */
    COOLPARTICLE _CoolParticle;	/* Abundances and any other cooling
				   internal variables */
#endif
    
 public:
    inline double& u() {return _u;}
    inline double& fMetals() {return _fMetals;}
    inline Vector3D<double>& vPred() {return _vPred;}
    inline double& uPred() {return _uPred;}
    inline double& divv() {return _divv;}
    inline Vector3D<double>& curlv() {return _curlv;}
    inline double& mumax() {return _mumax;}
    inline double& PdV() {return _PdV;}
    inline double& c() {return _c;}
    inline double& PoverRho2() {return _PoverRho2;}
    inline double& BalsaraSwitch() {return _BalsaraSwitch;}
    inline double& fBallMax() {return _fBallMax;}
#ifndef COOLING_NONE
    inline double& uDot() {return _uDot;}
    inline COOLPARTICLE& CoolParticle() {return _CoolParticle;}
#endif
    void pup(PUP::er &p) {
	p | _u;
	p | _fMetals;
	p | _vPred;
	p | _uPred;
	p | _divv;
	p | _curlv;
	p | _mumax;
	p | _PdV;
	p |  _c;
	p | _PoverRho2;
	p | _BalsaraSwitch;
	p | _fBallMax;
#ifndef COOLING_NONE
	p | _uDot;
	p((char *) &_CoolParticle, sizeof(_CoolParticle)); /* PUPs as bytes */
#endif
	}
    };

// Extra data needed for Stars
class extraStarData 
{
 private:
    double _fMetals;		/* Metalicity */
    double _fTimeForm;		/* Formation time */
    double _fMassForm;		/* Formation mass */
    int64_t _iGasOrder;		/* Gas from which this star formed */
 public:
    inline double& fMetals() {return _fMetals;}
    inline double& fTimeForm() {return _fTimeForm;}
    inline double& fMassForm() {return _fMassForm;}
    inline int64_t& iGasOrder() {return _iGasOrder;}
    void pup(PUP::er &p) {
	p | _fMetals;
	p | _fTimeForm;
	p | _fMassForm;
	p | _iGasOrder;
	}
    };

class GravityParticle;
int TYPETest(GravityParticle *a, unsigned int b);

class ExternalSmoothParticle;

// This class contains everything that a "dark matter" particle needs.
// Other classes of particles require this plus an "extra data" class.

class GravityParticle : public ExternalGravityParticle {
public:
	SFC::Key key;
	Vector3D<double> velocity;
	Vector3D<double> treeAcceleration;
	double potential;
	double dtGrav;
	double fBall;
	double fDensity;
	int iOrder;		/* input order of particles */
        int rung;  ///< the current rung (greater means faster)
	unsigned int iType;	// Bitmask to hold particle type information
#ifdef CHANGESOFT
	double fSoft0;
#endif
#ifdef NEED_DT
	double dt;
#endif
	void *extraData;	/* SPH or Star particle data */

#if COSMO_STATS > 1
	double intcellmass;
	double intpartmass;
	double extcellmass;
	double extpartmass;
#endif

        int interParticles;
	
	GravityParticle(SFC::Key k = 0) : ExternalGravityParticle() {
          key = k;
          rung = 0;
        }

	inline bool operator<(const GravityParticle& p) const {
		return key < p.key;
	}

	void pup(PUP::er &p) {
          ExternalGravityParticle::pup(p);
          p | key;
          p | velocity;
	  p | treeAcceleration;
	  p | fDensity;
	  p | fBall;
          p | iOrder;
          p | rung;
	  p | iType;
#ifdef CHANGESOFT
	  p | fSoft0;
#endif
#ifdef NEED_DT
	  p | dt;
#endif
        }
	// Access SPH quantities
	ExternalSmoothParticle getExternalSmoothParticle();
	inline double& u() { return (((extraSPHData*)extraData)->u());}
	inline double& fMetals() { return (((extraSPHData*)extraData)->fMetals());}
	inline Vector3D<double>& vPred() { return (((extraSPHData*)extraData)->vPred());}
	inline double& uPred() { return (((extraSPHData*)extraData)->uPred());}
	inline double& divv() { return (((extraSPHData*)extraData)->divv());}
	inline Vector3D<double>& curlv() { return (((extraSPHData*)extraData)->curlv());}
	inline double& mumax() { return (((extraSPHData*)extraData)->mumax());}
	inline double& PdV() { return (((extraSPHData*)extraData)->PdV());}
	inline double& c() { return (((extraSPHData*)extraData)->c());}
	inline double& PoverRho2() { return (((extraSPHData*)extraData)->PoverRho2());}
	inline double& BalsaraSwitch() { return (((extraSPHData*)extraData)->BalsaraSwitch());}
	inline double& fBallMax() { return (((extraSPHData*)extraData)->fBallMax());}
#ifndef COOLING_NONE
	inline double& uDot() { return (((extraSPHData*)extraData)->uDot());}
	inline COOLPARTICLE& CoolParticle() { return (((extraSPHData*)extraData)->CoolParticle());}
#endif
	// Access Star Quantities
	// XXX Beware overlaps with SPH; we could fix this by aligning
	// all common variables up at the start of the extraData structure.
	inline double& fStarMetals() { return (((extraStarData*)extraData)->fMetals());}
	inline double& fTimeForm() { return (((extraStarData*)extraData)->fTimeForm());}
	inline double& fMassForm() { return (((extraStarData*)extraData)->fMassForm());}
	inline int64_t& iGasOrder() { return (((extraStarData*)extraData)->iGasOrder());}

/* Particle Type Masks */

#define TYPE_GAS               (1<<0)
#define TYPE_DARK              (1<<1)
#define TYPE_STAR              (1<<2)

#define TYPE_DELETED           (1<<3)

#define TYPE_PHOTOGENIC        (1<<4)
#define TYPE_NbrOfACTIVE       (1<<5)

	inline bool isDark() { return TYPETest(this, TYPE_DARK);}
	inline bool isGas() { return TYPETest(this, TYPE_GAS);}
	inline bool isStar() { return TYPETest(this, TYPE_STAR);}

        GravityParticle &operator=(const ExternalGravityParticle &p){
          mass = p.mass;
          soft = p.soft;
          position = p.position;
        }

        GravityParticle &operator=(const PushExternalGravityParticle &p){
          /*
          mass = p.mass;
          soft = p.soft;
          position = p.position;
          */

          (*this) = *((const ExternalGravityParticle *)&p);

          rung = p.rung;
          dtGrav = p.dtGrav;
#ifdef CHANGESOFT
          fSoft0 = p.fSoft0;
#endif
#ifdef NEED_DT
          dt = p.dt;
#endif

        }
};

PushExternalGravityParticle &PushExternalGravityParticle::operator=(const GravityParticle &p){
  *((ExternalGravityParticle *)this) = p;

  rung = p.rung;
  dtGrav = p.dtGrav;
#ifdef CHANGESOFT
  fSoft0 = p.fSoft0;
#endif
#ifdef NEED_DT
  dt = p.dt;
#endif

}


inline int TYPETest(GravityParticle *a, unsigned int b) {
    return a->iType & b;
    }
inline int TYPESet(GravityParticle *a, unsigned int b) {
    return a->iType |= b;
    }
inline int TYPEReset(GravityParticle *a, unsigned int b) {
    return a->iType &= (~b);
    }

/// unmark particle as deleted
inline void unDeleteParticle(GravityParticle *p)
{
    CkAssert(TYPETest(p, TYPE_DELETED)); 

    TYPEReset(p, TYPE_DELETED); 
    }

/// mark particle as deleted
inline void deleteParticle(GravityParticle *p)
{
    TYPESet(p, TYPE_DELETED); 
    }

// Convert star particle to gas particle
// Note that new memory is allocated for the extradata.
inline GravityParticle StarFromGasParticle(GravityParticle *p) 
{
    GravityParticle starp = *p;
    starp.extraData = new extraStarData;
    starp.fStarMetals() = p->fMetals();
    return starp;
    }

// Class for cross processor data needed for smooth operations
class ExternalSmoothParticle {
 public:

  double mass;
  double fBall;
  double fDensity;
  Vector3D<double> position;
  Vector3D<double> velocity;
  unsigned int iType;	// Bitmask to hold particle type information
  int rung;
  Vector3D<double> vPred;
  Vector3D<double> treeAcceleration;
  double mumax;
  double PdV;
  double c;
  double PoverRho2;
  double BalsaraSwitch;
  double fBallMax;
  double u;
  double uPred;
  double uDot;
  double fMetals;

  ExternalSmoothParticle() {}

  ExternalSmoothParticle(GravityParticle *p) 
      {
	  mass = p->mass;
	  fBall = p->fBall;
	  fDensity = p->fDensity;
	  position = p->position;
	  velocity = p->velocity;
	  iType = p->iType;
	  rung = p->rung;
	  treeAcceleration = p->treeAcceleration;
	  if(TYPETest(p, TYPE_GAS)) {
	      vPred = p->vPred();
	      mumax = p->mumax();
	      PdV = p->PdV();
	      c = p->c();
	      PoverRho2 = p->PoverRho2();
	      BalsaraSwitch = p->BalsaraSwitch();
	      fBallMax = p->fBallMax();
	      u = p->u();
#ifndef COOLING_NONE
	      uDot = p->uDot();
	      uPred = p->uPred();
#endif
	      fMetals = p->fMetals();
	      }
	  }
  
  inline void getParticle(GravityParticle *tmp) { 
      tmp->mass = mass;
      tmp->fBall = fBall;
      tmp->fDensity = fDensity;
      tmp->position = position;
      tmp->velocity = velocity;
      tmp->iType = iType;
      tmp->rung = rung;
      tmp->treeAcceleration = treeAcceleration;
      if(TYPETest(tmp, TYPE_GAS)) {
	  tmp->vPred() = vPred;
	  tmp->mumax() = mumax;
	  tmp->PdV() = PdV;
	  tmp->c() = c;
	  tmp->PoverRho2() = PoverRho2;
	  tmp->BalsaraSwitch() = BalsaraSwitch;
	  tmp->fBallMax() = fBallMax;
	  tmp->u() = u;
#ifndef COOLING_NONE
	  tmp->uDot() = uDot;
	  tmp->uPred() = uPred;
#endif
	  tmp->fMetals() = fMetals;
	  }
      }
	  
  void pup(PUP::er &p) {
    p | position;
    p | velocity;
    p | mass;
    p | fBall;
    p | fDensity;
    p | iType;
    p | rung;
    p | treeAcceleration;
    p | vPred;
    p | mumax;
    p | PdV;
    p | c;
    p | PoverRho2;
    p | BalsaraSwitch;
    p | fBallMax;
    p | u;
    p | uPred;
    p | uDot;
    p | fMetals;
  }
};

inline ExternalSmoothParticle GravityParticle::getExternalSmoothParticle()
{ return ExternalSmoothParticle(this); }

#endif
