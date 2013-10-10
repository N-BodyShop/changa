
class SIDMSmoothParams : public SmoothParams
{
  protected:
  	double a, H; // Cosmological parameters
  	virtual void fcnSmooth(GravityParticle *p, int nSmooth, pqSmoothNode *nList);
    	virtual int isSmoothActive(GravityParticle *p);
    	virtual void initTreeParticle(GravityParticle *p);
    	virtual void postTreeParticle(GravityParticle *p) {}
    	virtual void initSmoothParticle(GravityParticle *p);
    	virtual void initSmoothCache(GravityParticle *p);
    	virtual void combSmoothCache(GravityParticle *p1,ExternalSmoothParticle *p2);

  public:
    SIDMSmoothParams() {} //empty constructor
    SIDMSmoothParams(int _iType, int am, CSM csm, double dTime) {
    iType = _iType;

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
        SmoothParams::pup(p);//Call base class
        p|a;
        p|H;
        }
    };

