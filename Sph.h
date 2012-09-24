/* Classes to describe smooth functions needed for SPH */
#ifndef __SPH_H
#define __SPH_H

/// @brief Parameters and functions for the first SPH smooth: density
/// and velocity derivatives.
class DenDvDxSmoothParams : public SmoothParams
{
 protected:
    double a, H; // Cosmological parameters
    int bActiveOnly;
    
    virtual void fcnSmooth(GravityParticle *p, int nSmooth,
			   pqSmoothNode *nList);
    virtual int isSmoothActive(GravityParticle *p);
    virtual void initTreeParticle(GravityParticle *p);
    virtual void postTreeParticle(GravityParticle *p) {} 
    virtual void initSmoothParticle(GravityParticle *p);
    virtual void initSmoothCache(GravityParticle *p);
    virtual void combSmoothCache(GravityParticle *p1,
				 ExternalSmoothParticle *p2);
 public:
    DenDvDxSmoothParams() {}
    /// @param _iType Type of particle to operate on
    /// @param am Active rung
    /// @param csm Cosmology information
    /// @param dTime Current time
    /// @param _bActiveOnly Only operate on active particles.
    DenDvDxSmoothParams(int _iType, int am, CSM csm, double dTime,
			int _bActiveOnly) {
	iType = _iType;
	activeRung = am;
	bActiveOnly = _bActiveOnly;
	if(csm->bComove) {
	    H = csmTime2Hub(csm,dTime);
	    a = csmTime2Exp(csm,dTime);
	    }
	else {
	    H = 0.0;
	    a = 1.0;
	    }
    }
    PUPable_decl(DenDvDxSmoothParams);
    DenDvDxSmoothParams(CkMigrateMessage *m) : SmoothParams(m) {}
    virtual void pup(PUP::er &p) {
        SmoothParams::pup(p);//Call base class
	p|a;
	p|H;
	p|bActiveOnly;
	}
    };

/// @brief Get density and velocity derivatives of "Neighbor of
/// Active" particles
///
/// Like the above class, but only does non active particles marked
/// "Neighbor of Active" for the "fast gas" option.
/// Also, it doesn't mark any particles.  This is used in the
/// "FastGas" step.

class DenDvDxNeighborSmParams : public DenDvDxSmoothParams
{
    virtual void fcnSmooth(GravityParticle *p, int nSmooth,
			   pqSmoothNode *nList);
    virtual int isSmoothActive(GravityParticle *p);
    virtual void initTreeParticle(GravityParticle *p) {}
    virtual void postTreeParticle(GravityParticle *p) {} 
    virtual void initSmoothParticle(GravityParticle *p) {}
    virtual void initSmoothCache(GravityParticle *p) {}
    virtual void combSmoothCache(GravityParticle *p1,
				 ExternalSmoothParticle *p2) {}
 public:
    DenDvDxNeighborSmParams() {}
    /// @param _iType Type of particle to operate on
    /// @param am Active rung
    /// @param csm Cosmology information
    /// @param dTime Current time
    DenDvDxNeighborSmParams(int _iType, int am, CSM csm, double dTime)
	: DenDvDxSmoothParams(_iType, am, csm, dTime, 0) {}
    PUPable_decl(DenDvDxNeighborSmParams);
    DenDvDxNeighborSmParams(CkMigrateMessage *m) : DenDvDxSmoothParams(m) {}
    virtual void pup(PUP::er &p) {
        DenDvDxSmoothParams::pup(p);//Call base class
	}
    };

/// @brief Parameters for "Mark Smooth", used to find inverse nearest
/// neighbors.

class MarkSmoothParams : public SmoothParams
{
    virtual void fcnSmooth(GravityParticle *p, int nSmooth,
			   pqSmoothNode *nList) {}
    virtual int isSmoothActive(GravityParticle *p);
    virtual void initTreeParticle(GravityParticle *p) {}
    virtual void postTreeParticle(GravityParticle *p) {} 
    virtual void initSmoothParticle(GravityParticle *p) {}
    virtual void initSmoothCache(GravityParticle *p) {}
    virtual void combSmoothCache(GravityParticle *p1,
				 ExternalSmoothParticle *p2) {}
 public:
    MarkSmoothParams() {}
    /// @param _iType Type of particle to operate on
    /// @param am Active rung
    MarkSmoothParams(int _iType, int am) {
	iType = _iType;
	activeRung = am;
	}
    PUPable_decl(MarkSmoothParams);
    MarkSmoothParams(CkMigrateMessage *m) : SmoothParams(m) {}
    virtual void pup(PUP::er &p) {
        SmoothParams::pup(p);//Call base class
	}
    };

/// @brief Second pass in SPH: calculate pressure forces.
class PressureSmoothParams : public SmoothParams
{
    double a, H; // Cosmological parameters
    double alpha, beta; // SPH viscosity parameters
    
    virtual void fcnSmooth(GravityParticle *p, int nSmooth,
			   pqSmoothNode *nList);
    virtual int isSmoothActive(GravityParticle *p);
    virtual void initTreeParticle(GravityParticle *p) {}
    virtual void postTreeParticle(GravityParticle *p) {} 
    virtual void initSmoothParticle(GravityParticle *p);
    virtual void initSmoothCache(GravityParticle *p);
    virtual void combSmoothCache(GravityParticle *p1,
				 ExternalSmoothParticle *p2);
 public:
    PressureSmoothParams() {}
    /// @param _iType Type of particles to smooth
    /// @param am Active rung
    /// @param csm Cosmological parameters
    /// @param dTime Current time
    /// @param _alpha Artificial viscosity parameter
    /// @param _beta Artificial viscosity parameter
    PressureSmoothParams(int _iType, int am, CSM csm, double dTime,
			 double _alpha, double _beta) {
	iType = _iType;
	activeRung = am;
	if(csm->bComove) {
	    H = csmTime2Hub(csm,dTime);
	    a = csmTime2Exp(csm,dTime);
	    }
	else {
	    H = 0.0;
	    a = 1.0;
	    }
	alpha = _alpha;
	beta = _beta;
    }
    PUPable_decl(PressureSmoothParams);
    PressureSmoothParams(CkMigrateMessage *m) : SmoothParams(m) {}
    virtual void pup(PUP::er &p) {
        SmoothParams::pup(p);//Call base class
	p|a;
	p|H;
	p|alpha;
	p|beta;
	}
    };

///
/// @brief SmoothParams class for distributing deleted gas to neighboring
/// particles.
///
class DistDeletedGasSmoothParams : public SmoothParams
{
    virtual void fcnSmooth(GravityParticle *p, int nSmooth,
			   pqSmoothNode *nList);
    virtual int isSmoothActive(GravityParticle *p);
    virtual void initSmoothParticle(GravityParticle *p) {};
    virtual void initTreeParticle(GravityParticle *p) {}
    virtual void postTreeParticle(GravityParticle *p) {}
    virtual void initSmoothCache(GravityParticle *p);
    virtual void combSmoothCache(GravityParticle *p1,
				 ExternalSmoothParticle *p2);
 public:
    DistDeletedGasSmoothParams() {}
    /// @param _iType Type of particles to smooth
    /// @param am Active rung
    DistDeletedGasSmoothParams(int _iType, int am) {
	iType = _iType;
	activeRung = am;
	bUseBallMax = 0;
	}
    PUPable_decl(DistDeletedGasSmoothParams);
    DistDeletedGasSmoothParams(CkMigrateMessage *m) : SmoothParams(m) {}
    virtual void pup(PUP::er &p) {
        SmoothParams::pup(p);//Call base class
	}
    };

#endif
