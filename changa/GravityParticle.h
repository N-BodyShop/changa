#ifndef GRAVITYPARTICLE_H
#define GRAVITYPARTICLE_H

#include "SFC.h"
#include <vector>

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
    void pup(PUP::er &p) {
	p | _u;
	p | _fMetals;
	}
    };

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
};

#define TYPE_GAS               (1<<0)
#define TYPE_DARK              (1<<1)
#define TYPE_STAR              (1<<2)
#define TYPE_PHOTOGENIC        (1<<3)
#define TYPE_NbrOfACTIVE       (1<<4)

inline int TYPETest(GravityParticle *a, unsigned int b) {
    return a->iType & b;
    }
inline int TYPESet(GravityParticle *a, unsigned int b) {
    return a->iType |= b;
    }


// Class for cross processor data needed for smooth operations
class ExternalSmoothParticle {
 public:

  double mass;
  double fBall;
  double fDensity;
  Vector3D<double> position;
  unsigned int iType;	// Bitmask to hold particle type information
  int rung;
  Vector3D<double> vPred;
  Vector3D<double> treeAcceleration;
  double mumax;
  double PdV;
  double c;
  double PoverRho2;
  double BalsaraSwitch;

  ExternalSmoothParticle() {}

  ExternalSmoothParticle(GravityParticle *p) 
      {
	  mass = p->mass;
	  fBall = p->fBall;
	  fDensity = p->fDensity;
	  position = p->position;
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
	      }
	  }
  
  inline void getParticle(GravityParticle *tmp) { 
      tmp->mass = mass;
      tmp->fBall = fBall;
      tmp->fDensity = fDensity;
      tmp->position = position;
      tmp->iType = iType;
      tmp->rung = rung;
      tmp->treeAcceleration = treeAcceleration;
      if(TYPETest(tmp, TYPE_GAS)) {
	  tmp->vPred() = vPred;
	  tmp->mumax() = mumax;
	  tmp->PdV() = PdV;
	  tmp->c() = c;
	  tmp->BalsaraSwitch() = BalsaraSwitch;
	  }
      }
	  
  void pup(PUP::er &p) {
    p | position;
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
  }
};

inline ExternalSmoothParticle GravityParticle::getExternalSmoothParticle()
{ return ExternalSmoothParticle(this); }

/* Particle Type Masks */

#endif
