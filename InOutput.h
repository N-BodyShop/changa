#ifndef __INOUTPUT_H
#define __INOUTPUT_H

class OutputParams;
#include "DataManager.h"

class OutputParams : public PUP::able 
{
 public:
    virtual double dValue(GravityParticle *p) = 0;
    virtual Vector3D<double> vValue(GravityParticle *p) = 0;
    int bVector;	// Is a vector, as opposed to a scalar
    std::string fileName;	// output file
    DataManager *dm;	// For extra state information (e.g. cooling)

    OutputParams() {dm = NULL;}
    PUPable_abstract(OutputParams);
    OutputParams(CkMigrateMessage *m) : PUP::able(m) {dm = NULL;}
    virtual void pup(PUP::er &p) {
        PUP::able::pup(p);//Call base class
        p|fileName;
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
    AccOutputParams(std::string _fileName) { bVector = 1; fileName = _fileName;}
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
    DenOutputParams(std::string _fileName) { bVector = 0; fileName = _fileName;}
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
    HsmOutputParams(std::string _fileName) { bVector = 0; fileName = _fileName;}
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
    PresOutputParams(std::string _fileName) { bVector = 0; fileName = _fileName;}
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
    DivVOutputParams(std::string _fileName) { bVector = 0; fileName = _fileName;}
    PUPable_decl(DivVOutputParams);
    DivVOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };

class PDVOutputParams : public OutputParams
{
    virtual double dValue(GravityParticle *p)
    {
	if (TYPETest(p, TYPE_GAS))
	    return p->PdV();
	else
	    return 0.0;
	}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
 public:
    PDVOutputParams() {}
    PDVOutputParams(std::string _fileName) { bVector = 0; fileName = _fileName;}
    PUPable_decl(PDVOutputParams);
    PDVOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };

class MuMaxOutputParams : public OutputParams
{
    virtual double dValue(GravityParticle *p)
    {
	if (TYPETest(p, TYPE_GAS))
	    return p->mumax();
	else
	    return 0.0;
	}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
 public:
    MuMaxOutputParams() {}
    MuMaxOutputParams(std::string _fileName) { bVector = 0; fileName = _fileName;}
    PUPable_decl(MuMaxOutputParams);
    MuMaxOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };

class BSwOutputParams : public OutputParams
{
    virtual double dValue(GravityParticle *p)
    {
	if (TYPETest(p, TYPE_GAS))
	    return p->BalsaraSwitch();
	else
	    return 0.0;
	}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
 public:
    BSwOutputParams() {}
    BSwOutputParams(std::string _fileName) { bVector = 0; fileName = _fileName;}
    PUPable_decl(BSwOutputParams);
    BSwOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };

class CsOutputParams : public OutputParams
{
    virtual double dValue(GravityParticle *p)
    {
	if (TYPETest(p, TYPE_GAS))
	    return p->c();
	else
	    return 0.0;
	}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
 public:
    CsOutputParams() {}
    CsOutputParams(std::string _fileName) { bVector = 0; fileName = _fileName;}
    PUPable_decl(CsOutputParams);
    CsOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };

class EDotOutputParams : public OutputParams
{
    virtual double dValue(GravityParticle *p)
    {
#ifndef COOLING_NONE
	CkAssert(dm != NULL);
	if (TYPETest(p, TYPE_GAS)) {
	    double r[3];  // For conversion to C
	    p->position.array_form(r);
	    return (COOL_EDOT(dm->Cool, &p->CoolParticle(), p->u(), p->fDensity, p->fMetals(), r));
	    }
	else
#endif
	    return 0.0;
	}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
 public:
    EDotOutputParams() {}
    EDotOutputParams(std::string _fileName) { bVector = 0; fileName = _fileName;}
    PUPable_decl(EDotOutputParams);
    EDotOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };

class Cool0OutputParams : public OutputParams
{
    virtual double dValue(GravityParticle *p)
    {
#ifndef COOLING_NONE
	if (TYPETest(p, TYPE_GAS))
	    return COOL_ARRAY0(unused1, &p->CoolParticle(), unused2);
	else
#endif
	    return 0.0;
	}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
 public:
    Cool0OutputParams() {}
    Cool0OutputParams(std::string _fileName) { bVector = 0; fileName = _fileName;}
    PUPable_decl(Cool0OutputParams);
    Cool0OutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };

class Cool1OutputParams : public OutputParams
{
    virtual double dValue(GravityParticle *p)
    {
#ifndef COOLING_NONE
	if (TYPETest(p, TYPE_GAS))
	    return COOL_ARRAY1(unused1, &p->CoolParticle(), unused2);
	else
#endif
	    return 0.0;
	}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
 public:
    Cool1OutputParams() {}
    Cool1OutputParams(std::string _fileName) { bVector = 0; fileName = _fileName;}
    PUPable_decl(Cool1OutputParams);
    Cool1OutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };

class Cool2OutputParams : public OutputParams
{
    virtual double dValue(GravityParticle *p)
    {
#ifndef COOLING_NONE
	if (TYPETest(p, TYPE_GAS))
	    return COOL_ARRAY2(unused1, &p->CoolParticle(), unused2);
	else
#endif
	    return 0.0;
	}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
 public:
    Cool2OutputParams() {}
    Cool2OutputParams(std::string _fileName) { bVector = 0; fileName = _fileName;}
    PUPable_decl(Cool2OutputParams);
    Cool2OutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };

class OxOutputParams : public OutputParams
{
    virtual double dValue(GravityParticle *p) {
	if (TYPETest(p, TYPE_STAR)) return p->fStarMFracOxygen();
	if (TYPETest(p, TYPE_GAS)) return p->fMFracOxygen();
	else return 0.0;
	}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
 public:
    OxOutputParams() {}
    OxOutputParams(std::string achFileName) { 
	bVector = 0; 
	fileName = achFileName+".OxMassFrac";
	}
    PUPable_decl(OxOutputParams);
    OxOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };

class FeOutputParams : public OutputParams
{
    virtual double dValue(GravityParticle *p) {
	if (TYPETest(p, TYPE_STAR)) return p->fStarMFracIron();
	if (TYPETest(p, TYPE_GAS)) return p->fMFracIron();
	else return 0.0;
	}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
 public:
    FeOutputParams() {}
    FeOutputParams(std::string achFileName) { 
	bVector = 0; 
	fileName = achFileName+".FeMassFrac";
	}
    PUPable_decl(FeOutputParams);
    FeOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };

class MetalsDotOutputParams : public OutputParams
{
    virtual double dValue(GravityParticle *p) {
#ifdef DIFFUSION
	if (p->isGas()) return p->fMetalsDot();
	else
	    return 0.0;
#endif
	}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
 public:
    MetalsDotOutputParams() {}
    MetalsDotOutputParams(std::string achFileName) { 
	bVector = 0; 
	fileName = achFileName+".Metalsdot";
	}
    PUPable_decl(MetalsDotOutputParams);
    MetalsDotOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };

class OxygenMassFracDotOutputParams : public OutputParams
{
    virtual double dValue(GravityParticle *p) {
#ifdef DIFFUSION
	if (p->isGas()) return p->fMFracOxygenDot();
	else
	    return 0.0;
#endif
	}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
 public:
    OxygenMassFracDotOutputParams() {}
    OxygenMassFracDotOutputParams(std::string achFileName) { 
	bVector = 0; 
	fileName = achFileName+".OxMassFracdot";
	}
    PUPable_decl(OxygenMassFracDotOutputParams);
    OxygenMassFracDotOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };

class IronMassFracDotOutputParams : public OutputParams
{
    virtual double dValue(GravityParticle *p) {
#ifdef DIFFUSION
	if (p->isGas()) return p->fMFracIronDot();
	else
	    return 0.0;
#endif
	}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
 public:
    IronMassFracDotOutputParams() {}
    IronMassFracDotOutputParams(std::string achFileName) { 
	bVector = 0; 
	fileName = achFileName+".FeMassFracdot";
	}
    PUPable_decl(IronMassFracDotOutputParams);
    IronMassFracDotOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };


class coolontimeOutputParams : public OutputParams
{
    virtual double dValue(GravityParticle *p) {
	if (TYPETest(p, TYPE_GAS)) return p->fTimeCoolIsOffUntil();
	else return 0.0;
	}
    virtual Vector3D<double> vValue(GravityParticle *p)
    {CkAssert(0); return 0.0;}
 public:
    coolontimeOutputParams() {}
    coolontimeOutputParams(std::string achFileName) { 
	bVector = 0; 
	fileName = achFileName+".coolontime";
	}
    PUPable_decl(coolontimeOutputParams);
    coolontimeOutputParams(CkMigrateMessage *m) {}
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
    DtOutputParams(std::string _fileName) { bVector = 0; fileName = _fileName;}
    PUPable_decl(DtOutputParams);
    DtOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };

class RungOutputParams : public OutputParams
{
    virtual double dValue(GravityParticle *p)
    {
	return (double) p->rung;
	}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
 public:
    RungOutputParams() {}
    RungOutputParams(std::string _fileName) { bVector = 0; fileName = _fileName;}
    PUPable_decl(RungOutputParams);
    RungOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };
#endif
