/** @file GravityParticle.h
 * Defines the fundamental particle data structures.
 */
#ifndef GRAVITYPARTICLE_H
#define GRAVITYPARTICLE_H

#include <charm.h>              /* for CkAssert */
#include "cooling.h"
#include "cosmoType.h"
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

  cosmoType mass;
  cosmoType soft;
  Vector3D<cosmoType> position;

#ifdef __CHARMC__
  void pup(PUP::er &p) {
    p | position;
    p | mass;
    p | soft;
  }
#endif
};

/// @brief Extra data needed for SPH
class extraSPHData 
{
 private:
    double _u;			/* Internal Energy */
    double _fMetals;		/* Metalicity */
#ifdef SPLITGAS
    int64_t _iSplitOrder;	/* Gas from which this particle split*/
#endif
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
    double _TimeDivV;           /* Time at which dvds was last updated */
    double _dvds;
    double _dvdsOnSFull;        ///< dvds for use in the Cullen R calculation
    double _dvds_old;           ///< Save dvds(OnSFull) for R calculation
#endif
#ifdef DTADJUST
    double _dtNew;		/* New timestep from gas pressure */
#endif
    double _dTimeFB;		/* Track feedback time */
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
#ifdef SUPERBUBBLE
    COOLPARTICLE _CoolParticleHot;	
    int _cpHotInit; /* Do we need to initialize the Hot Coolparticle? */
    double _uHot; /* Hot phase energy */
    double _uHotDot; /* Hot phase rate of energy change */
    double _uHotPred; /* Hot phase predicted energy */
    double _massHot; /* Hot phase mass*/
    double _fDensityU; /* Energy-scaled density */
    double _fThermalCond; /* Conduction rate */
    double _fThermalLength; /* Conduction length */
    double _fPromoteSum; /* Total evaporated mass */
    double _fPromoteSumuPred; /* Total evaporating energy */
    double _fPromoteuPredInit; /* Original energy pre-evaporation */
#endif
    
 public:
    inline double& u() {return _u;}
#ifdef SPLITGAS
    inline int64_t& iSplitOrder() {return _iSplitOrder;}
#endif
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
    inline const double CullenAlpha() const {return _CullenAlpha;}
    inline double& CullenAlpha() {return _CullenAlpha;}
    inline double& TimeDivV() {return _TimeDivV;}
    inline double& dvds() {return _dvds;}
    inline double& dvdsOnSFull() {return _dvdsOnSFull;}
    inline double& dvds_old() {return _dvds_old;}
#endif
#ifdef DTADJUST
    inline double& dtNew() {return _dtNew;}
#endif
    inline double& dTimeFB() {return _dTimeFB;}
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
#ifdef SUPERBUBBLE
    inline COOLPARTICLE& CoolParticleHot() {return _CoolParticleHot;}
    inline int& cpHotInit() {return _cpHotInit;}
    inline double& uHot() {return _uHot;}
    inline double& uHotPred() {return _uHotPred;}
    inline double& uHotDot() {return _uHotDot;}
    inline double& massHot() {return _massHot;}
    inline double& fDensityU() {return _fDensityU;}
    inline double& fThermalCond() {return _fThermalCond;}
    inline double& fThermalLength() {return _fThermalLength;}
    inline double& fPromoteSum() {return _fPromoteSum;}
    inline double& fPromoteSumuPred() {return _fPromoteSumuPred;}
    inline double& fPromoteuPredInit() {return _fPromoteuPredInit;}
#endif
#ifdef __CHARMC__
    void pup(PUP::er &p) {
	p | _u;
#ifdef SPLITGAS
	p | _iSplitOrder;
#endif
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
        p | _dvds;
        p | _dvdsOnSFull;
        p | _dvds_old;
#endif 
#ifdef DTADJUST
        p | _dtNew;
#endif
	p | _dTimeFB;
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
#ifdef SUPERBUBBLE
	p((char *) &_CoolParticleHot, sizeof(_CoolParticle)); /* PUPs as bytes */
    p| _cpHotInit;
    p| _uHot;
    p| _uHotDot;
    p| _uHotPred;
    p| _massHot;
    p| _fDensityU;
    p| _fThermalCond;
    p| _fThermalLength;
    p| _fPromoteSum;
    p| _fPromoteSumuPred;
    p| _fPromoteuPredInit;
#endif
	}
#endif
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
    int64_t _iEaterOrder;	/* iOrder for merging black holes */
    double _dMDot;		/* Accretion rate of black holes */
    double _dDeltaM;		/* Actual Mass Accreted on black holes */
#ifdef COOLING_MOLECULARH
    double _dStarLymanWerner;	/* Lyman Werner radiation emmited from star particles */
#endif /*COOLING_MOLECULARH*/
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
    inline int64_t& iEaterOrder() {return _iEaterOrder;}
    inline double& dMDot() {return _dMDot;}
    inline double& dDeltaM() {return _dDeltaM;}
#ifdef COOLING_MOLECULARH
    inline const double dStarLymanWerner() const {return _dStarLymanWerner;} 
    inline double& dStarLymanWerner() {return _dStarLymanWerner;} 
#endif /*COOLING_MOLECULARH*/
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
	p | _iEaterOrder;
	p | _dMDot;
	p | _dDeltaM;
#ifdef COOLING_MOLECULARH
	p | _dStarLymanWerner;
#endif /*COOLINg_MOLECULARH*/
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
	Vector3D<cosmoType> treeAcceleration;
	cosmoType potential;
        cosmoType dtGrav;       ///< timestep from gravity
        double fBall;           ///< Neighbor search radius for smoothing
	double fDensity;
        int64_t iOrder;	///< Input order of particles; unique particle ID
        int rung;  ///< the current rung (greater means faster)
        unsigned int iType;	///< Bitmask to hold particle type information
#ifdef SIDMINTERACT
        int iNSIDMInteractions; // SIDM number of interactions
#endif
#ifdef CHANGESOFT
	cosmoType fSoft0;
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

        cosmoType interMass;
	
        GravityParticle(SFC::Key k) : ExternalGravityParticle() {
            key = k;
            }
        GravityParticle() : ExternalGravityParticle() {
            }

	/// @brief Used to sort the particles into tree order.
	inline bool operator<(const GravityParticle& p) const {
		return key < p.key;
	}

#ifdef __CHARMC__
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
#ifdef SIDMINTERACT
          p | iNSIDMInteractions; // SIDM
#endif
#ifdef CHANGESOFT
	  p | fSoft0;
#endif
#ifdef NEED_DT
	  p | dt;
#endif
        }
#endif

// Debugging macros for the extra data fields.
// To enable, define GP_DEBUG_EXTRAS

#define GP_DEBUG_EXTRAS

#ifdef GP_DEBUG_EXTRAS
/// Debugging macro to be sure you are accessing gas properties from a
/// gas particle.
#define IMAGAS CkAssert(isGas())
/// Debugging macro to be sure you are accessing star properties from a
/// star particle.
#define IMASTAR CkAssert(isStar())
#else
#define IMAGAS
#define IMASTAR
#endif

	// Access SPH quantities
	/// @brief Get quantities needed for SPH smooths.
	ExternalSmoothParticle getExternalSmoothParticle();
	inline double& u() { IMAGAS; return (((extraSPHData*)extraData)->u());}
#ifdef SPLITGAS
	inline int64_t& iSplitOrder() { IMAGAS; return (((extraSPHData*)extraData)->iSplitOrder());}
#endif
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
        inline const double CullenAlpha() const {IMAGAS; return (((extraSPHData*)extraData)->CullenAlpha());}
        inline double& CullenAlpha() {IMAGAS; return (((extraSPHData*)extraData)->CullenAlpha());}
        inline double& TimeDivV() {IMAGAS; return (((extraSPHData*)extraData)->TimeDivV());}
        inline double& dvds() {IMAGAS; return (((extraSPHData*)extraData)->dvds());}
        inline double& dvdsOnSFull() {IMAGAS; return (((extraSPHData*)extraData)->dvdsOnSFull());}
        inline double& dvds_old() {IMAGAS; return (((extraSPHData*)extraData)->dvds_old());}
#endif
#ifdef DTADJUST
        inline double& dtNew() { IMAGAS; return (((extraSPHData*)extraData)->dtNew());}
#endif
	inline double& dTimeFB() { IMAGAS; return (((extraSPHData*)extraData)->dTimeFB());}
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
#ifdef SUPERBUBBLE
	inline COOLPARTICLE& CoolParticleHot() { IMAGAS; return (((extraSPHData*)extraData)->CoolParticleHot());}
	inline int& cpHotInit() { IMAGAS; return (((extraSPHData*)extraData)->cpHotInit());}
	inline double& uHot() { IMAGAS; return (((extraSPHData*)extraData)->uHot());}
	inline double& uHotPred() { IMAGAS; return (((extraSPHData*)extraData)->uHotPred());}
	inline double& uHotDot() { IMAGAS; return (((extraSPHData*)extraData)->uHotDot());}
	inline double& massHot() { IMAGAS; return (((extraSPHData*)extraData)->massHot());}
	inline double& fDensityU() { IMAGAS; return (((extraSPHData*)extraData)->fDensityU());}
	inline double& fThermalCond() { IMAGAS; return (((extraSPHData*)extraData)->fThermalCond());}
	inline double& fThermalLength() { IMAGAS; return (((extraSPHData*)extraData)->fThermalLength());}
	inline double& fPromoteSum() { IMAGAS; return (((extraSPHData*)extraData)->fPromoteSum());}
	inline double& fPromoteSumuPred() { IMAGAS; return (((extraSPHData*)extraData)->fPromoteSumuPred());}
	inline double& fPromoteuPredInit() { IMAGAS; return (((extraSPHData*)extraData)->fPromoteuPredInit());}
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
	inline int64_t& iEaterOrder() { IMASTAR; return (((extraStarData*)extraData)->iEaterOrder());}
	inline double& dDeltaM() { IMASTAR; return (((extraStarData*)extraData)->dDeltaM());}
	inline double& dMDot() { IMASTAR; return (((extraStarData*)extraData)->dMDot());}
#ifdef COOLING_MOLECULARH
	inline const double dStarLymanWerner() const { IMASTAR; return (((extraStarData*)extraData)->dStarLymanWerner());}	
	inline double& dStarLymanWerner() { IMASTAR; return (((extraStarData*)extraData)->dStarLymanWerner());} 
#endif /*COOLING_MOLECULARH*/

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
#define TYPE_SMOOTHACTIVE      (1<<6)

#define TYPE_SINK              (1<<7)
#define TYPE_SINKING           (1<<8)
#define TYPE_NEWSINKING        (1<<9)
#define TYPE_PROMOTED          (1<<10)
#define TYPE_FEEDBACK          (1<<11)
#define TYPE_MAXTYPE           (1<<12)

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

/// @brief Test for a type flag.
inline int TYPETest(const GravityParticle *a, unsigned int b) {
    return a->iType & b;
    }
/// @brief Set a type flag.
inline int TYPESet(GravityParticle *a, unsigned int b) {
    return a->iType |= b;
    }
/// @brief Unset a type flag.
inline int TYPEReset(GravityParticle *a, unsigned int b) {
    return a->iType &= (~b);
    }
inline int TYPEClear(GravityParticle *a) {
    return a->iType = 0;
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

  cosmoType mass;
  double fBall;
  double fDensity;
  Vector3D<cosmoType> position;
  Vector3D<double> velocity;
  int64_t iOrder;
  unsigned int iType;	// Bitmask to hold particle type information
  int rung;
#ifdef DTADJUST
  double dt;
  double dtNew;
#endif
  Vector3D<double> vPred;
  Vector3D<cosmoType> treeAcceleration;
  double mumax;
  double PdV;
  double c;
  double PoverRho2;
  double BalsaraSwitch;
  double fBallMax;
#ifdef CULLENALPHA
  double CullenAlpha;
  double TimeDivV;
  double dvds;
  double dvds_old;
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
#ifdef SUPERBUBBLE
  COOLPARTICLE CoolParticle;	
  double uHot;
  double uHotDot;
  double uHotPred;
  double massHot;
  double fDensityU;
  double fThermalCond;
  double fThermalLength;
  double fPromoteSum;
  double fPromoteSumuPred;
  double fPromoteuPredInit;
#endif
  double fNSN;
  int64_t iEaterOrder;
  double dTimeFB;
  int iBucketOff;               /* Used by the Cache */

  ExternalSmoothParticle() {}

  ExternalSmoothParticle(GravityParticle *p) 
      {
	  mass = p->mass;
	  fBall = p->fBall;
	  fDensity = p->fDensity;
	  position = p->position;
	  velocity = p->velocity;
	  iOrder = p->iOrder;
#ifdef SPLITGAS
	  iWriteOrder = p->iWriteOrder;
#endif
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
              dvds = p->dvds();
              dvds_old = p->dvds_old();
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
#ifdef SUPERBUBBLE
          CoolParticle = p->CoolParticle();
          uHot = p->uHot();
          uHotPred = p->uHotPred();
          uHotDot = p->uHotDot();
          massHot = p->massHot();
          fDensityU = p->fDensityU();
          fThermalCond = p->fThermalCond();
          fThermalLength = p->fThermalLength();
          fPromoteSum = p->fPromoteSum();
          fPromoteSumuPred = p->fPromoteSumuPred();
          fPromoteuPredInit = p->fPromoteuPredInit();
#endif
#ifdef DTADJUST
              dt = p->dt;
              dtNew = p->dtNew();
#endif
	      dTimeFB = p->dTimeFB();
	      }
	  if(TYPETest(p, TYPE_STAR)) {
	      fNSN = p->fNSN();
	      iEaterOrder = p->iEaterOrder();
	      }
	  }
  
  /// @brief Fill in a full gravity particle from this object.
  inline void getParticle(GravityParticle *tmp) const { 
      tmp->mass = mass;
      tmp->fBall = fBall;
      tmp->fDensity = fDensity;
      tmp->position = position;
      tmp->velocity = velocity;
      tmp->iOrder = iOrder;
#ifdef SPLITGAS
      tmp->iWriteOrder = iWriteOrder;
#endif
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
          tmp->dvds() = dvds;
          tmp->dvds_old() = dvds_old;
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
#ifdef SUPERBUBBLE
      tmp->CoolParticle() = CoolParticle;
      tmp->uHot() = uHot;
      tmp->uHotPred() = uHotPred;
      tmp->uHotDot() = uHotDot;
      tmp->massHot() = massHot;
      tmp->fDensityU() = fDensityU;
      tmp->fThermalCond() = fThermalCond;
      tmp->fThermalLength() = fThermalLength;
      tmp->fPromoteSum() = fPromoteSum;
      tmp->fPromoteSumuPred() = fPromoteSumuPred;
      tmp->fPromoteuPredInit() = fPromoteuPredInit;
#endif
#ifdef DTADJUST
          tmp->dt = dt;
          tmp->dtNew() = dtNew;
#endif
	  tmp->dTimeFB() = dTimeFB;
	  }
      if(TYPETest(tmp, TYPE_STAR)) {
	  tmp->fNSN() = fNSN;
	  tmp->iEaterOrder() = iEaterOrder;
	  }
      }
	  
#ifdef __CHARMC__
  void pup(PUP::er &p) {
    p | position;
    p | velocity;
    p | mass;
    p | fBall;
    p | fDensity;
    p | iOrder;
#ifdef SPLITGAS
    p | iWriteOrder;
#endif
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
    p | dvds;
    p | dvds_old;
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
#ifdef SUPERBUBBLE
	p((char *) &CoolParticle, sizeof(CoolParticle)); /* PUPs as bytes */
    p | uHot;
    p | uHotPred;
    p | uHotDot;
    p | massHot;
    p | fDensityU;
    p | fThermalCond;
    p | fThermalLength;
    p | fPromoteSum;
    p | fPromoteSumuPred;
    p | fPromoteuPredInit;
#endif
    p | fNSN;
    p | iEaterOrder;
    p | dTimeFB;
    p | iBucketOff;
  }
#endif
};

inline ExternalSmoothParticle GravityParticle::getExternalSmoothParticle()
{ return ExternalSmoothParticle(this); }

inline int TYPETest(ExternalSmoothParticle *a, unsigned int b) {
    return a->iType & b;
    }

#endif
