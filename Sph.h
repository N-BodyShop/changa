/* Classes to describe smooth functions needed for SPH */
#ifndef __SPH_H
#define __SPH_H

#ifdef DIFFUSION

#if defined(FEEDBACKDIFFLIMIT) && !defined(DIFFUSIONHARMONIC)
#define DIFFUSIONHARMONIC
#endif

#endif

/// @brief Parameters and functions for the first SPH smooth: density
/// and velocity derivatives.
class DenDvDxSmoothParams : public SmoothParams
{
 protected:
    double a, H; // Cosmological parameters
    int bActiveOnly;
    int bConstantDiffusion;
    double dTime; 
    double dAlphaMax;           ///< Maximum SPH alpha
    int bStarting;              ///< We are starting (or restarting)
                                ///  the simulation
    int bHaveAlpha;             ///< Alpha has been read in.

    virtual void fcnSmooth(GravityParticle *p, int nSmooth,
			   pqSmoothNode *nList);
    virtual int isSmoothActive(GravityParticle *p);
    virtual void initTreeParticle(GravityParticle *p);
    virtual void postTreeParticle(GravityParticle *p); 
    virtual void initSmoothParticle(GravityParticle *p);
    virtual void initSmoothCache(GravityParticle *p);
    virtual void combSmoothCache(GravityParticle *p1,
				 ExternalSmoothParticle *p2);
 public:
    DenDvDxSmoothParams() {}
    /// @param _iType Type of particle to operate on
    /// @param am Active rung
    /// @param csm Cosmology information
    /// @param _dTime Current time
    /// @param _bActiveOnly Only operate on active particles.
    /// @param _bConstantDiffusion Fixed diffusion constant
    /// @param _bStarting Simulation is starting
    /// @param _bHaveAlpha No need to calculate alpha
    /// @param _dAlphaMax Maximum SPH alpha
    DenDvDxSmoothParams(int _iType, int am, CSM csm, double _dTime,
			int _bActiveOnly, int _bConstantDiffusion,
                        int _bStarting, int _bHaveAlpha, double _dAlphaMax) {
	iType = _iType;
	activeRung = am;
	bActiveOnly = _bActiveOnly;
	bConstantDiffusion = _bConstantDiffusion;
        bStarting = _bStarting;
        bHaveAlpha = _bHaveAlpha;
	dAlphaMax = _dAlphaMax;
        dTime = _dTime;
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
        p|dTime;
	p|a;
	p|H;
	p|bActiveOnly;
	p|bConstantDiffusion;
	p|bStarting;
	p|bHaveAlpha;
	p|dAlphaMax;
	}
    };

/// @brief Get density and velocity derivatives of "Neighbor of
/// Active" particles
///
/// Like the above class, but only does non active particles marked
/// "Neighbor of Active" for the "fast gas" option.
/// Also, marking of neighbors is not complete.  This is used in the
/// "FastGas" step.  It uses the same fcnSmooth() as DenDvDxSmoothParams.

class DenDvDxNeighborSmParams : public DenDvDxSmoothParams
{
    virtual int isSmoothActive(GravityParticle *p);
    virtual void initTreeParticle(GravityParticle *p) {}
    virtual void postTreeParticle(GravityParticle *p);
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
    /// @param bConstantDiffusion Fixed diffusion constant
    /// @param dAlphaMax Maximum SPH alpha
    /// This calls the DenDvDx constructor, and we assume bActiveOnly,
    /// bStarting, and bHaveAlpha are not set.
    DenDvDxNeighborSmParams(int _iType, int am, CSM csm, double dTime,
			    int bConstDiffusion, double dAlphaMax)
      : DenDvDxSmoothParams(_iType, am, csm, dTime, 0, bConstDiffusion,
                            0, 0, dAlphaMax) {}
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
				 ExternalSmoothParticle *p2);
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
    double dTime;
    double a, H; // Cosmological parameters
    double alpha, beta; // SPH viscosity parameters
    double dThermalDiffusionCoeff;
    double dMetalDiffusionCoeff;
    double dtFacCourant; // Courant timestep factor
    double dtFacDiffusion; // Diffusion timestep factor
    
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
    PressureSmoothParams(int _iType, int am, CSM csm, double _dTime,
			 double _alpha, double _beta,
                         double _dThermalDiff, double _dMetalDiff,
                         double dEtaCourant, double dEtaDiffusion) {
	iType = _iType;
	activeRung = am;
        dTime = _dTime;
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
	dThermalDiffusionCoeff = _dThermalDiff;
	dMetalDiffusionCoeff = _dMetalDiff;
	dtFacCourant = dEtaCourant*a*2.0/1.6;
	dtFacDiffusion = 2.0*dEtaDiffusion;
    }
    PUPable_decl(PressureSmoothParams);
    PressureSmoothParams(CkMigrateMessage *m) : SmoothParams(m) {}
    virtual void pup(PUP::er &p) {
        SmoothParams::pup(p);//Call base class
        p|dTime;
	p|a;
	p|H;
	p|alpha;
	p|beta;
	p|dThermalDiffusionCoeff;
	p|dMetalDiffusionCoeff;
	p|dtFacCourant;
	p|dtFacDiffusion;
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
