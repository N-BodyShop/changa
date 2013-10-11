
class SIDMSmoothParams : public SmoothParams
{
  protected:
  	double a, H; // Cosmological parameters
        double dTime;
        double dDelta;
        double dSIDMSigma;
  	virtual void fcnSmooth(GravityParticle *p, int nSmooth, pqSmoothNode *nList);
    	virtual int isSmoothActive(GravityParticle *p);
    	virtual void initTreeParticle(GravityParticle *p);
    	virtual void postTreeParticle(GravityParticle *p) {}
    	virtual void initSmoothParticle(GravityParticle *p);
    	virtual void initSmoothCache(GravityParticle *p);
    	virtual void combSmoothCache(GravityParticle *p1,ExternalSmoothParticle *p2);

  public:
    SIDMSmoothParams() {} //empty constructor
    SIDMSmoothParams(int _iType, int am, CSM csm, double _dTime, double _dSIDMSigma, double _dDelta) {
    iType = _iType;
    activeRung = am;
    dTime = _dTime;
    dDelta = _dDelta;
    dSIDMSigma= _dSIDMSigma;


        if(csm->bComove) {
            H = csmTime2Hub(csm,dTime);
            a = csmTime2Exp(csm,dTime);
            }
        else {
            H = 0.0;
            a = 1.0;
            }

    }
    PUPable_decl(SIDMSmoothParams);
    SIDMSmoothParams(CkMigrateMessage *m) : SmoothParams(m) {}
    virtual void pup(PUP::er &p) {
        SmoothParams::pup(p);//Call  base class
        p|dTime;
        p|dDelta;
        p|dSIDMSigma;
        p|a;
        p|H;
        }
    };

