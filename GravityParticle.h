/** @file GravityParticle.h
 * Defines the fundamental particle data structures.
 */
#ifndef GRAVITYPARTICLE_H
#define GRAVITYPARTICLE_H

#include "cooling.h"

#include "SFC.h"
#include <vector>

#ifdef DTADJUST
#define NEED_DT
#endif

/// @brief Object to bookkeep a Bucket Walk.
class BucketGravityRequest {
public:
		
	int finished;

	BucketGravityRequest(unsigned int bucketSize = 0) : finished(0) {
	}
	
};

/// @brief Information needed to calculate gravity
///
/// This is used in the CacheParticle class since it contains only the
/// information that an external procesor needs to calculate gravity.
///
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

/// @brief Extra data needed for SPH
class extraSPHData 
{
 private:
    double _u;			/* Internal Energy */
    double _fMetals;		/* Metalicity */
    double _fMFracOxygen;	/* Oxygen mass fraction  */
    double _fMFracIron;		/* Iron mass fraction  */
    double _fESNrate;		/* SN energy rate  */
    double _fTimeCoolIsOffUntil;/* time cooling is turned back on */
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
#ifdef CULLENALPHA
    double _CullenAlpha;        /* Alpha from Cullen & Dehnen 2010 */
    double _TimeDivV;           /* Time at which OldDivV was last updated */
    double _OldDivV;            /* Last active divv */
#endif
#ifdef DTADJUST
    double _dtNew;		/* New timestep from gas pressure */
#endif
#ifndef COOLING_NONE
    double _uDot;		/* Rate of change of u, for
				   predicting u */
    COOLPARTICLE _CoolParticle;	/* Abundances and any other cooling
				   internal variables */
#endif
#ifdef DIFFUSION
    double _diff;		/* Diffusion coefficient, based on Smagorinski  */
    double _fMetalsDot;
    double _fMetalsPred;
    double _fMFracOxygenDot;
    double _fMFracIronDot;
    double _fMFracOxygenPred;
    double _fMFracIronPred;
#endif
    
 public:
    inline double& u() {return _u;}
    inline double& fMetals() {return _fMetals;}
    inline double& fMFracOxygen() {return _fMFracOxygen;}
    inline double& fMFracIron() {return _fMFracIron;}
    inline double& fESNrate() {return _fESNrate;}
    inline double& fTimeCoolIsOffUntil() {return _fTimeCoolIsOffUntil;}
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
#ifdef CULLENALPHA
    inline double& CullenAlpha() {return _CullenAlpha;}
    inline double& TimeDivV() {return _TimeDivV;}
    inline double& OldDivV() {return _OldDivV;}
#endif
#ifdef DTADJUST
    inline double& dtNew() {return _dtNew;}
#endif
#ifndef COOLING_NONE
    inline double& uDot() {return _uDot;}
    inline COOLPARTICLE& CoolParticle() {return _CoolParticle;}
#endif
#ifdef DIFFUSION
    inline double& diff() {return _diff;}
    inline double& fMetalsDot() {return _fMetalsDot;}
    inline double& fMetalsPred() {return _fMetalsPred;}
    inline double& fMFracOxygenDot() {return _fMFracOxygenDot;}
    inline double& fMFracOxygenPred() {return _fMFracOxygenPred;}
    inline double& fMFracIronDot() {return _fMFracIronDot;}
    inline double& fMFracIronPred() {return _fMFracIronPred;}
#endif
    void pup(PUP::er &p) {
	p | _u;
	p | _fMetals;
	p | _fMFracIron;
	p | _fMFracOxygen;
	p | _fESNrate;
	p | _fTimeCoolIsOffUntil;
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
#ifdef CULLENALPHA
        p | _CullenAlpha;
        p | _TimeDivV;
        p | _OldDivV;
#endif 
#ifdef DTADJUST
        p | _dtNew;
#endif
#ifndef COOLING_NONE
	p | _uDot;
	p((char *) &_CoolParticle, sizeof(_CoolParticle)); /* PUPs as bytes */
#endif
#ifdef DIFFUSION
	p| _diff;
	p| _fMetalsDot;
	p| _fMetalsPred;
	p| _fMFracOxygenDot;
	p| _fMFracOxygenPred;
	p| _fMFracIronDot;
	p| _fMFracIronPred;
#endif
	}
    };

/// @brief Extra data needed for Stars
class extraStarData 
{
 private:
    double _fMetals;		/* Metalicity */
    double _fTimeForm;		/* Formation time */
    double _fMassForm;		/* Formation mass */
    double _fESNrate;		/* SN energy rate  */
    double _fNSN;               /* number of SN exploding */
    double _fMSN;               /* mass of feedback ejecta */
    double _fMFracOxygen;	/* Oxygen mass fraction  */
    double _fMFracIron;		/* Iron mass fraction  */
    double _fSNMetals;          /* Ejected metals from feedback */
    double _fMOxygenOut;        /* Ejected oxygen */
    double _fMIronOut;          /* Ejected iron */
    int64_t _iGasOrder;		/* Gas from which this star formed */
 public:
    inline double& fMetals() {return _fMetals;}
    inline double& fTimeForm() {return _fTimeForm;}
    inline double& fMassForm() {return _fMassForm;}
    inline double& fESNrate() {return _fESNrate;}
    inline double& fNSN() {return _fNSN;}
    inline double& fMSN() {return _fMSN;}
    inline double& fMFracOxygen() {return _fMFracOxygen;}
    inline double& fMFracIron() {return _fMFracIron;}
    inline double& fMIronOut() {return _fMIronOut;}
    inline double& fMOxygenOut() {return _fMOxygenOut;}
    inline double& fSNMetals() {return _fSNMetals;}
    inline int64_t& iGasOrder() {return _iGasOrder;}
    void pup(PUP::er &p) {
	p | _fMetals;
	p | _fTimeForm;
	p | _fMassForm;
	p | _fESNrate;
	p | _fNSN;    
	p | _fMSN;    
	p | _fMFracOxygen;
	p | _fMFracIron;
	p | _fSNMetals;
	p | _fMOxygenOut;
	p | _fMIronOut;
	p | _iGasOrder;
	}
    };

class GravityParticle;
int TYPETest(const GravityParticle *a, unsigned int b);

class ExternalSmoothParticle;

/// @brief Fundamental type for a particle
///
/// This class contains everything that a "dark matter" particle needs.
/// Other classes of particles require this plus an "extra data" class.

class GravityParticle : public ExternalGravityParticle {
public:
	SFC::Key key;
	Vector3D<double> velocity;
	Vector3D<double> treeAcceleration;
	double potential;
	double dtGrav;
	double fBall;
	double fDensity;
	int64_t iOrder;		/* input order of particles */
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

        double interMass;
	
        GravityParticle(SFC::Key k) : ExternalGravityParticle() {
            key = k;
            }
        GravityParticle() : ExternalGravityParticle() {
            }

	/// @brief Used to sort the particles into tree order.
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

// Debugging macros for the extra data fields.
// To enable, define GP_DEBUG_EXTRAS

#define GP_DEBUG_EXTRAS

#ifdef GP_DEBUG_EXTRAS
#define IMAGAS CkAssert(isGas())
#define IMASTAR CkAssert(isStar())
#else
#define IMAGAS
#define IMASTAR
#endif

	// Access SPH quantities
	/// @brief Get quantities needed for SPH smooths.
	ExternalSmoothParticle getExternalSmoothParticle();
	inline double& u() { IMAGAS; return (((extraSPHData*)extraData)->u());}
	inline double& fMetals() { IMAGAS; return (((extraSPHData*)extraData)->fMetals());}
	inline double& fMFracOxygen() {IMAGAS; return (((extraSPHData*)extraData)->fMFracOxygen());}
	inline double& fMFracIron() {IMAGAS; return (((extraSPHData*)extraData)->fMFracIron());}
	inline double& fESNrate() {IMAGAS; return (((extraSPHData*)extraData)->fESNrate());}
	inline double& fTimeCoolIsOffUntil() {IMAGAS; return (((extraSPHData*)extraData)->fTimeCoolIsOffUntil());}
	inline Vector3D<double>& vPred() { IMAGAS; return (((extraSPHData*)extraData)->vPred());}
	inline double& uPred() {IMAGAS;  return (((extraSPHData*)extraData)->uPred());}
	inline double& divv() { IMAGAS; return (((extraSPHData*)extraData)->divv());}
	inline Vector3D<double>& curlv() { IMAGAS; return (((extraSPHData*)extraData)->curlv());}
	inline double& mumax() { IMAGAS; return (((extraSPHData*)extraData)->mumax());}
	inline double& PdV() { IMAGAS; return (((extraSPHData*)extraData)->PdV());}
	inline double& c() { IMAGAS; return (((extraSPHData*)extraData)->c());}
	inline double& PoverRho2() { IMAGAS; return (((extraSPHData*)extraData)->PoverRho2());}
	inline double& BalsaraSwitch() { IMAGAS; return (((extraSPHData*)extraData)->BalsaraSwitch());}
	inline double& fBallMax() { IMAGAS; return (((extraSPHData*)extraData)->fBallMax());}
#ifdef CULLENALPHA
        inline double& CullenAlpha() {IMAGAS; return (((extraSPHData*)extraData)->CullenAlpha());}
        inline double& TimeDivV() {IMAGAS; return (((extraSPHData*)extraData)->TimeDivV());}
        inline double& OldDivV() {IMAGAS; return (((extraSPHData*)extraData)->OldDivV());}
#endif
#ifdef DTADJUST
        inline double& dtNew() { IMAGAS; return (((extraSPHData*)extraData)->dtNew());}
#endif
#ifndef COOLING_NONE
	inline double& uDot() { IMAGAS; return (((extraSPHData*)extraData)->uDot());}
	inline COOLPARTICLE& CoolParticle() { IMAGAS; return (((extraSPHData*)extraData)->CoolParticle());}
#endif
#ifdef DIFFUSION
	inline double& diff() { IMAGAS; return (((extraSPHData*)extraData)->diff());}
	inline double& fMetalsDot() { IMAGAS; return (((extraSPHData*)extraData)->fMetalsDot());}
	inline double& fMetalsPred() { IMAGAS; return (((extraSPHData*)extraData)->fMetalsPred());}
	inline double& fMFracOxygenDot() { IMAGAS; return (((extraSPHData*)extraData)->fMFracOxygenDot());}
	inline double& fMFracIronDot() { IMAGAS; return (((extraSPHData*)extraData)->fMFracIronDot());}
	inline double& fMFracOxygenPred() { IMAGAS; return (((extraSPHData*)extraData)->fMFracOxygenPred());}
	inline double& fMFracIronPred() { IMAGAS; return (((extraSPHData*)extraData)->fMFracIronPred());}
#endif
	// Access Star Quantities
	// XXX Beware overlaps with SPH; we could fix this by aligning
	// all common variables up at the start of the extraData structure.
	inline double& fStarMetals() { IMASTAR; return (((extraStarData*)extraData)->fMetals());}
	inline double& fStarMFracOxygen() {IMASTAR; return (((extraStarData*)extraData)->fMFracOxygen());}
	inline double& fStarMFracIron() {IMASTAR; return (((extraStarData*)extraData)->fMFracIron());}
	inline double& fTimeForm() { IMASTAR; return (((extraStarData*)extraData)->fTimeForm());}
	inline double& fMassForm() { IMASTAR; return (((extraStarData*)extraData)->fMassForm());}
	inline double& fStarESNrate() {IMASTAR; return (((extraStarData*)extraData)->fESNrate());}
	inline double& fNSN() {IMASTAR; return (((extraStarData*)extraData)->fNSN());}
	inline double& fMSN() {IMASTAR; return (((extraStarData*)extraData)->fMSN());}
	inline double& fMIronOut() {IMASTAR; return (((extraStarData*)extraData)->fMIronOut());}
	inline double& fMOxygenOut() {IMASTAR; return (((extraStarData*)extraData)->fMOxygenOut());}
	inline double& fSNMetals() {IMASTAR; return (((extraStarData*)extraData)->fSNMetals());}
	inline int64_t& iGasOrder() { IMASTAR; return (((extraStarData*)extraData)->iGasOrder());}

// See above debugging macros
#undef IMAGAS
#undef IMASTAR

/* Particle Type Masks */

#define TYPE_GAS               (1<<0)
#define TYPE_DARK              (1<<1)
#define TYPE_STAR              (1<<2)

#define TYPE_DELETED           (1<<3)

#define TYPE_PHOTOGENIC        (1<<4)
#define TYPE_NbrOfACTIVE       (1<<5)
#define TYPE_MAXTYPE           (1<<6)

	inline bool isDark() const { return TYPETest(this, TYPE_DARK);}
	inline bool isGas() const { return TYPETest(this, TYPE_GAS);}
	inline bool isStar() const { return TYPETest(this, TYPE_STAR);}

        GravityParticle &operator=(const ExternalGravityParticle &p){
          mass = p.mass;
          soft = p.soft;
          position = p.position;
	  return *this;
        }
};

inline int TYPETest(const GravityParticle *a, unsigned int b) {
    return a->iType & b;
    }
inline int TYPESet(GravityParticle *a, unsigned int b) {
    return a->iType |= b;
    }
inline int TYPEReset(GravityParticle *a, unsigned int b) {
    return a->iType &= (~b);
    }

/// @brief unmark particle as deleted
inline void unDeleteParticle(GravityParticle *p)
{
    CkAssert(TYPETest(p, TYPE_DELETED)); 

    TYPEReset(p, TYPE_DELETED); 
    }

/// @brief mark particle as deleted
inline void deleteParticle(GravityParticle *p)
{
    TYPESet(p, TYPE_DELETED); 
    }

/// @brief Create star particle from gas particle
/// Note that new memory is allocated for the extradata.
inline GravityParticle StarFromGasParticle(GravityParticle *p) 
{
    GravityParticle starp = *p;
    TYPESet(&starp, TYPE_STAR);
    starp.extraData = new extraStarData;
    starp.fStarMetals() = p->fMetals();
    starp.fStarMFracOxygen() = p->fMFracOxygen();
    starp.fStarMFracIron() = p->fMFracIron();
    return starp;
    }

/// @brief Class for cross processor data needed for smooth operations
class ExternalSmoothParticle {
 public:

  double mass;
  double fBall;
  double fDensity;
  Vector3D<double> position;
  Vector3D<double> velocity;
  unsigned int iType;	// Bitmask to hold particle type information
  int rung;
#ifdef DTADJUST
  double dt;
  double dtNew;
#endif
  Vector3D<double> vPred;
  Vector3D<double> treeAcceleration;
  double mumax;
  double PdV;
  double c;
  double PoverRho2;
  double BalsaraSwitch;
  double fBallMax;
#ifdef CULLENALPHA
  double CullenAlpha;
  double TimeDivV;
  double OldDivV;
#endif
  double u;
  double uPred;
  double uDot;
  double fESNrate;
  double fMetals;
  double fMFracOxygen;
  double fMFracIron;
  double fTimeCoolIsOffUntil;
  Vector3D<double> curlv;	/* Curl of the velocity */
#ifdef DIFFUSION
  double diff;
  double fMetalsDot;
  double fMFracOxygenDot;
  double fMFracIronDot;
#endif
  int iBucketOff;               /* Used by the Cache */

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
#ifdef CULLENALPHA
	      CullenAlpha = p->CullenAlpha();
              TimeDivV = p->TimeDivV();
              OldDivV = p->OldDivV();
#endif
	      curlv = p->curlv();
	      u = p->u();
#ifndef COOLING_NONE
	      uDot = p->uDot();
#endif
	      uPred = p->uPred();
	      fMetals = p->fMetals();
	      fESNrate = p->fESNrate();
	      fMFracOxygen = p->fMFracOxygen();
	      fMFracIron = p->fMFracIron();
	      fTimeCoolIsOffUntil = p->fTimeCoolIsOffUntil();
#ifdef DIFFUSION
	      diff = p->diff();
	      fMetalsDot = p->fMetalsDot();
	      fMFracOxygenDot = p->fMFracOxygenDot();
	      fMFracIronDot = p->fMFracIronDot();
#endif	      
#ifdef DTADJUST
              dt = p->dt;
              dtNew = p->dtNew();
#endif
	      }
	  }
  
  /// @brief Fill in a full gravity particle from this object.
  inline void getParticle(GravityParticle *tmp) const { 
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
#ifdef CULLENALPHA
	  tmp->CullenAlpha() = CullenAlpha;
          tmp->TimeDivV() = TimeDivV;
          tmp->OldDivV() = OldDivV;
#endif
	  tmp->curlv() = curlv;
	  tmp->u() = u;
#ifndef COOLING_NONE
	  tmp->uDot() = uDot;
#endif
	  tmp->uPred() = uPred;
	  tmp->fMetals() = fMetals;
	  tmp->fESNrate() = fESNrate;
	  tmp->fMFracOxygen() = fMFracOxygen;
	  tmp->fMFracIron() = fMFracIron;
	  tmp->fTimeCoolIsOffUntil() = fTimeCoolIsOffUntil;
#ifdef DIFFUSION
	  tmp->diff() = diff;
	  tmp->fMetalsDot() = fMetalsDot;
	  tmp->fMFracOxygenDot() = fMFracOxygenDot;
	  tmp->fMFracIronDot() = fMFracIronDot;
#endif
#ifdef DTADJUST
          tmp->dt = dt;
          tmp->dtNew() = dtNew;
#endif
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
#ifdef DTADJUST
    p | dt;
    p | dtNew;
#endif
    p | treeAcceleration;
    p | vPred;
    p | mumax;
    p | PdV;
    p | c;
    p | PoverRho2;
    p | BalsaraSwitch;
    p | fBallMax;
#ifdef CULLENALPHA
    p | CullenAlpha;
    p | TimeDivV;
    p | OldDivV;
#endif
    p | u;
    p | uPred;
    p | uDot;
    p | curlv;
    p | fMetals;
    p | fESNrate;
    p | fMFracOxygen;
    p | fMFracIron;
    p | fTimeCoolIsOffUntil;
#ifdef DIFFUSION
    p | diff;
    p | fMetalsDot;
    p | fMFracOxygenDot;
    p | fMFracIronDot;
#endif
    p | iBucketOff;
  }
};

inline ExternalSmoothParticle GravityParticle::getExternalSmoothParticle()
{ return ExternalSmoothParticle(this); }

#endif
