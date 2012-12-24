#ifndef SINKS_HINCLUDED
#define SINKS_HINCLUDED

class Sinks 
{
    //
    // Sink parameters
    //
 public:
    int bDoSinks;
    int bBHSink;
    int bSmallBHSmooth;
    int bBHTurnOffCooling;
    int bDoBHKick;
    double dDeltaStarForm;
    double dKmPerSecUnit;
    double dBHSinkEddEff;
    double dBHSinkFeedbackEff;
    double dBHSinkAlpha;
    double dBHSinkEddFactor;
    double dBHSinkFeedbackFactor;
    int bBHMindv;
    int bBHAccreteAll;
    int bDoSinksAtStart;
    int bSinkThermal;
    double dSinkRadius;
    double dSinkBoundOrbitRadius;
    double dSinkMustAccreteRadius;
    double dDeltaSink;
    double dSinkCurrentDelta;
    double dSinkMassMin;
    int iSinkRung;
    int iSinkCurrentRung;
    int bSinkForm;
    int nJeans;
    double dJeansConstant;
    int bSinkFormJeans;
    int bSinkFormDivV;
    double dSinkFormDivVCoeff;
    int bSinkFormDivAcc;
    double dSinkFormDivAccCoeff;
    int bSinkFormDV;
    int bSinkFormPotMin;
    double dSinkFormDensity;
    double dSinkTimeEligible;
    int bSinkFormSimple;
    int nSinkFormMin;
 public:
    void AddParams(PRM prm, struct parameters &param);
    void CheckParams(PRM prm, struct parameters &param);
    inline void pup(PUP::er &p);
    };

inline void Sinks::pup(PUP::er &p) {
    p|bDoSinks;
    p|bBHSink;
    p|bSmallBHSmooth;
    p|bBHTurnOffCooling;
    p|bDoBHKick;
    p|dDeltaStarForm;
    p|dKmPerSecUnit;
    p|dBHSinkEddEff;
    p|dBHSinkFeedbackEff;
    p|dBHSinkAlpha;
    p|dBHSinkEddFactor;
    p|dBHSinkFeedbackFactor;
    p|bBHMindv;
    p|bBHAccreteAll;
    p|bDoSinksAtStart;
    p|bSinkThermal;
    p|dSinkRadius;
    p|dSinkBoundOrbitRadius;
    p|dSinkMustAccreteRadius;
    p|dDeltaSink;
    p|dSinkCurrentDelta;
    p|dSinkMassMin;
    p|iSinkRung;
    p|iSinkCurrentRung;
    p|bSinkForm;
    p|nJeans;
    p|dJeansConstant;
    p|bSinkFormJeans;
    p|bSinkFormDivV;
    p|dSinkFormDivVCoeff;
    p|bSinkFormDivAcc;
    p|dSinkFormDivAccCoeff;
    p|bSinkFormDV;
    p|bSinkFormPotMin;
    p|dSinkFormDensity;
    p|dSinkTimeEligible;
    p|bSinkFormSimple;
    p|nSinkFormMin;
    }

class SinkFormTestSmoothParams : public SmoothParams
{
 protected:
    Sinks s;
    
    virtual void fcnSmooth(GravityParticle *p, int nSmooth,
			   pqSmoothNode *nList);
    virtual int isSmoothActive(GravityParticle *p);
    virtual void initTreeParticle(GravityParticle *p) { initSmoothCache(p); }
    virtual void postTreeParticle(GravityParticle *p) {} 
    virtual void initSmoothParticle(GravityParticle *p) { initSmoothCache(p); }
    virtual void initSmoothCache(GravityParticle *p);
    virtual void combSmoothCache(GravityParticle *p1,
				 ExternalSmoothParticle *p2);
 public:
    SinkFormTestSmoothParams() {}
    SinkFormTestSmoothParams(int _iType, int am, Sinks _s) {
	iType = _iType;
	activeRung = am;
	s = _s;
	bUseBallMax = 0;
    }
    PUPable_decl(SinkFormTestSmoothParams);
    SinkFormTestSmoothParams(CkMigrateMessage *m) : SmoothParams(m) {}
    virtual void pup(PUP::er &p) {
        SmoothParams::pup(p);//Call base class
	p|s;
	}
    };

class SinkFormSmoothParams : public SmoothParams
{
 protected:
    double dTime;
    double a, H; // Cosmological parameters
    Sinks s;
    
    virtual void fcnSmooth(GravityParticle *p, int nSmooth,
			   pqSmoothNode *nList);
    virtual int isSmoothActive(GravityParticle *p);
    virtual void initTreeParticle(GravityParticle *p) {}
    virtual void postTreeParticle(GravityParticle *p) {} 
    virtual void initSmoothParticle(GravityParticle *p) {}
    virtual void initSmoothCache(GravityParticle *p) {}
    virtual void combSmoothCache(GravityParticle *p1,
				 ExternalSmoothParticle *p2);
 public:
    SinkFormSmoothParams() {}
    SinkFormSmoothParams(int _iType, int am, CSM csm, double _dTime,
			 Sinks _s) {
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
	s = _s;
	bUseBallMax = 0;
    }
    PUPable_decl(SinkFormSmoothParams);
    SinkFormSmoothParams(CkMigrateMessage *m) : SmoothParams(m) {}
    virtual void pup(PUP::er &p) {
        SmoothParams::pup(p);//Call base class
	p|dTime;
	p|a;
	p|H;
	p|s;
	}
    };

class BHDensitySmoothParams : public SmoothParams
{
 protected:
    double dTime;
    double a, H; // Cosmological parameters
    Sinks s;
    
    virtual void fcnSmooth(GravityParticle *p, int nSmooth,
			   pqSmoothNode *nList);
    virtual int isSmoothActive(GravityParticle *p);
    virtual void initTreeParticle(GravityParticle *p);
    virtual void postTreeParticle(GravityParticle *p) {} 
    virtual void initSmoothParticle(GravityParticle *p) { initSmoothCache(p); } 
    virtual void initSmoothCache(GravityParticle *p);
    virtual void combSmoothCache(GravityParticle *p1,
				 ExternalSmoothParticle *p2);
 public:
    BHDensitySmoothParams() {}
    BHDensitySmoothParams(int _iType, int am, CSM csm, double _dTime,
			  Sinks _s) {
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
	s = _s;
    }
    PUPable_decl(BHDensitySmoothParams);
    BHDensitySmoothParams(CkMigrateMessage *m) : SmoothParams(m) {}
    virtual void pup(PUP::er &p) {
        SmoothParams::pup(p);//Call base class
	p|dTime;
	p|a;
	p|H;
	p|s;
	}
    };

class BHAccreteSmoothParams : public SmoothParams
{
 protected:
    double a, H; // Cosmological parameters
    double dTime;
    /// Big timestep to convert rungs into delta t.
    double dDelta;
    Sinks s;
    double gamma;	// Adiabatic index for pressure
    
    virtual void fcnSmooth(GravityParticle *p, int nSmooth,
			   pqSmoothNode *nList);
    virtual int isSmoothActive(GravityParticle *p);
    virtual void initTreeParticle(GravityParticle *p);
    virtual void postTreeParticle(GravityParticle *p); 
    virtual void initSmoothParticle(GravityParticle *p) {}
    virtual void initSmoothCache(GravityParticle *p);
    virtual void combSmoothCache(GravityParticle *p1,
				 ExternalSmoothParticle *p2);
 public:
    BHAccreteSmoothParams() {}
    BHAccreteSmoothParams(int _iType, int am, CSM csm, double _dTime,
			  double _dDelta, Sinks _s, double _gamma) {
	iType = _iType;
	activeRung = am;
	dTime = _dTime;
	dDelta = _dDelta;
	if(csm->bComove) {
	    H = csmTime2Hub(csm,dTime);
	    a = csmTime2Exp(csm,dTime);
	    }
	else {
	    H = 0.0;
	    a = 1.0;
	    }
	s = _s;
	gamma = _gamma;
    }
    PUPable_decl(BHAccreteSmoothParams);
    BHAccreteSmoothParams(CkMigrateMessage *m) : SmoothParams(m) {}
    virtual void pup(PUP::er &p) {
        SmoothParams::pup(p);//Call base class
	p|a;
	p|H;
	p|dTime;
	p|dDelta;
	p|s;
	p|gamma;
	}
    };

class BHIdentifySmoothParams : public SmoothParams
{
 protected:
    double a, H; // Cosmological parameters
    
    virtual void fcnSmooth(GravityParticle *p, int nSmooth,
			   pqSmoothNode *nList);
    virtual int isSmoothActive(GravityParticle *p);
    virtual void initTreeParticle(GravityParticle *p) {}
    virtual void postTreeParticle(GravityParticle *p) {} 
    virtual void initSmoothParticle(GravityParticle *p) {}
    virtual void initSmoothCache(GravityParticle *p);
    virtual void combSmoothCache(GravityParticle *p1,
				 ExternalSmoothParticle *p2);
 public:
    BHIdentifySmoothParams() {}
    BHIdentifySmoothParams(int _iType, int am, CSM csm, double dTime,
			  Sinks _s) {
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
    }
    PUPable_decl(BHIdentifySmoothParams);
    BHIdentifySmoothParams(CkMigrateMessage *m) : SmoothParams(m) {}
    virtual void pup(PUP::er &p) {
        SmoothParams::pup(p);//Call base class
	p|a;
	p|H;
	}
    };

class BHSinkMergeSmoothParams : public SmoothParams
{
 protected:
    double a, H; // Cosmological parameters
    double dTime;
    Sinks s;
    
    virtual void fcnSmooth(GravityParticle *p, int nSmooth,
			   pqSmoothNode *nList);
    virtual int isSmoothActive(GravityParticle *p);
    virtual void initTreeParticle(GravityParticle *p) {}
    virtual void postTreeParticle(GravityParticle *p) {} 
    virtual void initSmoothParticle(GravityParticle *p) {}
    virtual void initSmoothCache(GravityParticle *p) {}
    virtual void combSmoothCache(GravityParticle *p1,
				 ExternalSmoothParticle *p2);
    
 public:
    BHSinkMergeSmoothParams() {}
    BHSinkMergeSmoothParams(int _iType, int am, CSM csm, double _dTime,
			  Sinks _s, double _gamma) {
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
	s = _s;
    }
    PUPable_decl(BHSinkMergeSmoothParams);
    BHSinkMergeSmoothParams(CkMigrateMessage *m) : SmoothParams(m) {}
    virtual void pup(PUP::er &p) {
        SmoothParams::pup(p);//Call base class
	p|a;
	p|H;
	p|dTime;
	p|s;
	}
    };

class SinkAccreteTestSmoothParams : public SmoothParams
{
 protected:
    Sinks s;
    
    virtual void fcnSmooth(GravityParticle *p, int nSmooth,
			   pqSmoothNode *nList);
    virtual int isSmoothActive(GravityParticle *p);
    virtual void initTreeParticle(GravityParticle *p) {}
    virtual void postTreeParticle(GravityParticle *p) {} 
    virtual void initSmoothParticle(GravityParticle *p) {} 
    virtual void initSmoothCache(GravityParticle *p);
    virtual void combSmoothCache(GravityParticle *p1,
				 ExternalSmoothParticle *p2);
 public:
    SinkAccreteTestSmoothParams() {}
    SinkAccreteTestSmoothParams(int _iType, int am, double dTime,
			  Sinks _s) {
	iType = _iType;
	activeRung = am;
	s = _s;
	bUseBallMax = 0;
    }
    PUPable_decl(SinkAccreteTestSmoothParams);
    SinkAccreteTestSmoothParams(CkMigrateMessage *m) : SmoothParams(m) {}
    virtual void pup(PUP::er &p) {
        SmoothParams::pup(p);//Call base class
	p|s;
	}
    };

class SinkingAverageSmoothParams : public SmoothParams
{
 protected:
    Sinks s;
    
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
    SinkingAverageSmoothParams() {}
    SinkingAverageSmoothParams(int _iType, int am, double dTime,
			  Sinks _s) {
	iType = _iType;
	activeRung = am;
	s = _s;
    }
    PUPable_decl(SinkingAverageSmoothParams);
    SinkingAverageSmoothParams(CkMigrateMessage *m) : SmoothParams(m) {}
    virtual void pup(PUP::er &p) {
        SmoothParams::pup(p);//Call base class
	p|s;
	}
    };

class SinkAccreteSmoothParams : public SmoothParams
{
 protected:
    Sinks s;
    
    virtual void fcnSmooth(GravityParticle *p, int nSmooth,
			   pqSmoothNode *nList);
    virtual int isSmoothActive(GravityParticle *p);
    virtual void initTreeParticle(GravityParticle *p) {}
    virtual void postTreeParticle(GravityParticle *p) {} 
    virtual void initSmoothParticle(GravityParticle *p) {}
    virtual void initSmoothCache(GravityParticle *p);
    virtual void combSmoothCache(GravityParticle *p1,
				 ExternalSmoothParticle *p2);
 public:
    SinkAccreteSmoothParams() {}
    SinkAccreteSmoothParams(int _iType, int am, double dTime,
			  Sinks _s) {
	iType = _iType;
	activeRung = am;
	s = _s;
	bUseBallMax = 0;
    }
    PUPable_decl(SinkAccreteSmoothParams);
    SinkAccreteSmoothParams(CkMigrateMessage *m) : SmoothParams(m) {}
    virtual void pup(PUP::er &p) {
        SmoothParams::pup(p);//Call base class
	p|s;
	}
    };
#endif
