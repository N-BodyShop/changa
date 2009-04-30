#ifndef __INOUTPUT_H
#define __INOUTPUT_H

class OutputParams : public PUP::able 
{
 public:
    virtual double dValue(GravityParticle *p) = 0;
    virtual Vector3D<double> vValue(GravityParticle *p) = 0;
    int bVector;	// Is a vector, as opposed to a scalar
    std::string suffix;	// suffix of output file

    OutputParams() {}
    PUPable_abstract(OutputParams);
    OutputParams(CkMigrateMessage *m) : PUP::able(m) {}
    virtual void pup(PUP::er &p) {
        PUP::able::pup(p);//Call base class
        p|suffix;
        p|bVector;
	}
    };

class AccOutputParams : public OutputParams
{
 public:
    virtual double dValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual Vector3D<double> vValue(GravityParticle *p)
				{return p->treeAcceleration;}
    AccOutputParams() {}
    AccOutputParams(std::string _suffix) { bVector = 1; suffix = _suffix;}
    PUPable_decl(AccOutputParams);
    AccOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };
    
class DenOutputParams : public OutputParams
{
    virtual double dValue(GravityParticle *p) {return p->fDensity;}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
 public:
    DenOutputParams() {}
    DenOutputParams(std::string _suffix) { bVector = 0; suffix = _suffix;}
    PUPable_decl(DenOutputParams);
    DenOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };

class HsmOutputParams : public OutputParams
{
    virtual double dValue(GravityParticle *p) {return 0.5*p->fBall;}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
 public:
    HsmOutputParams() {}
    HsmOutputParams(std::string _suffix) { bVector = 0; suffix = _suffix;}
    PUPable_decl(HsmOutputParams);
    HsmOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };

class PresOutputParams : public OutputParams
{
    virtual double dValue(GravityParticle *p)
    {
	if (TYPETest(p, TYPE_GAS))
	    return p->fDensity*p->fDensity*p->PoverRho2();
	else
	    return 0.0;
	}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
 public:
    PresOutputParams() {}
    PresOutputParams(std::string _suffix) { bVector = 0; suffix = _suffix;}
    PUPable_decl(PresOutputParams);
    PresOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };

class DivVOutputParams : public OutputParams
{
    virtual double dValue(GravityParticle *p)
    {
	if (TYPETest(p, TYPE_GAS))
	    return p->divv();
	else
	    return 0.0;
	}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
 public:
    DivVOutputParams() {}
    DivVOutputParams(std::string _suffix) { bVector = 0; suffix = _suffix;}
    PUPable_decl(DivVOutputParams);
    DivVOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };

class DtOutputParams : public OutputParams
{
    virtual double dValue(GravityParticle *p)
    {
#ifdef NEED_DT
	return p->dt;
#else
	return 0.0;
#endif
	}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
 public:
    DtOutputParams() {}
    DtOutputParams(std::string _suffix) { bVector = 0; suffix = _suffix;}
    PUPable_decl(DtOutputParams);
    DtOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };
#endif
