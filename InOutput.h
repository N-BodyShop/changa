/// @file InOutput.h
/// Declarations for I/O implemenatation.
#ifndef __INOUTPUT_H
#define __INOUTPUT_H

class OutputParams;
class OutputIntParams;
#include "DataManager.h"

/// @brief Base class for output parameters.
///
/// This is an abstract class from which an output parameter class can
/// be derived.  Derived classes need to implement dValue() or
/// vValue() which returns the value to be output for a given particle.
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

/// @brief Output accelerations.
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
    
/// @brief Output densities.
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

/// @brief Output smoothing lengths.
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

/// @brief Output gravitational softening.
class SoftOutputParams : public OutputParams
{
    virtual double dValue(GravityParticle *p) {return p->soft;}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
 public:
    SoftOutputParams() {}
    SoftOutputParams(std::string _fileName) { bVector = 0; fileName = _fileName;}
    PUPable_decl(SoftOutputParams);
    SoftOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };

/// @brief Output pressure.
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

/// @brief Output divergence of velocity.
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

/// @brief Output pressure times change in volume.
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

/// @brief Output artificial viscosity mumax.
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

/// @brief Output value of Balsara switch.
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

/// @brief Output sound speed.
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

/// @brief Output the cooling rate.
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

/// @brief Output the value in cool_array0.
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

/// @brief Output the value in cool_array1.
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

/// @brief Output the value in cool_array2.
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

class MFormOutputParams : public OutputParams
{
    virtual double dValue(GravityParticle *p) {
	if (TYPETest(p, TYPE_STAR)) return p->fMassForm();
	else return 0.0;
	}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
 public:
    MFormOutputParams() {}
    MFormOutputParams(std::string achFileName) { 
	bVector = 0; 
	fileName = achFileName+".massform";
	}
    PUPable_decl(MFormOutputParams);
    MFormOutputParams(CkMigrateMessage *m) {}
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

/// @brief Output Supernova heating rate
class ESNRateOutputParams : public OutputParams
{
    virtual double dValue(GravityParticle *p) {
	if (p->isGas()) return p->fESNrate();
	else if(p->isStar()) return p->fStarESNrate();
	else return 0.0;
	}
    virtual Vector3D<double> vValue(GravityParticle *p)
    {CkAssert(0); return 0.0;}
 public:
    ESNRateOutputParams() {}
    ESNRateOutputParams(std::string achFileName) { 
	bVector = 0; 
	fileName = achFileName+".ESNRate";
	}
    PUPable_decl(ESNRateOutputParams);
    ESNRateOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };

/// @brief Output timesteps.
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

/// @brief Output Keys.
class KeyOutputParams : public OutputParams
{
    virtual double dValue(GravityParticle *p)
    {
	return p->key;
	}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
 public:
    KeyOutputParams() {}
    KeyOutputParams(std::string _fileName) { bVector = 0; fileName = _fileName;}
    PUPable_decl(KeyOutputParams);
    KeyOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };

/// @brief Output Domains.
class DomainOutputParams : public OutputParams
{
    virtual double dValue(GravityParticle *p)
    {
	return p->interMass; // Hack: this gets assigned in assignDomain()
			     // just for this diagnostic.
	}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
 public:
    DomainOutputParams() {}
    DomainOutputParams(std::string _fileName) { bVector = 0; fileName = _fileName;}
    PUPable_decl(DomainOutputParams);
    DomainOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };

/// @brief Base class for Integer output parameters.
///
/// This is an abstract class from which an output parameter class can
/// be derived.  Derived classes need to implement iValue()
/// which returns the value to be output for a given particle.
class OutputIntParams : public PUP::able 
{
 public:
    virtual int iValue(GravityParticle *p) = 0;
    std::string fileName;	// output file

    OutputIntParams() {}
    PUPable_abstract(OutputIntParams);
    OutputIntParams(CkMigrateMessage *m) : PUP::able(m) {}
    virtual void pup(PUP::er &p) {
        PUP::able::pup(p);//Call base class
        p|fileName;
	}
    };

/// @brief Output iOrder.
class IOrderOutputParams : public OutputIntParams
{
    virtual int iValue(GravityParticle *p)
    {
	return p->iOrder;
	}
 public:
    IOrderOutputParams() {}
    IOrderOutputParams(std::string _fileName) { fileName = _fileName;}
    PUPable_decl(IOrderOutputParams);
    IOrderOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputIntParams::pup(p);//Call base class
	}
    };

/// @brief Output iGasOrder.
class IGasOrderOutputParams : public OutputIntParams
{
    virtual int iValue(GravityParticle *p)
    {
	if(p->isStar())
	    return p->iGasOrder();
	else
	    return 0;
	}
 public:
    IGasOrderOutputParams() {}
    IGasOrderOutputParams(std::string _fileName) { fileName = _fileName;}
    PUPable_decl(IGasOrderOutputParams);
    IGasOrderOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputIntParams::pup(p);//Call base class
	}
    };

/// @brief Output rungs.
class RungOutputParams : public OutputIntParams
{
    virtual int iValue(GravityParticle *p)
    {
	return p->rung;
	}
 public:
    RungOutputParams() {}
    RungOutputParams(std::string _fileName) { fileName = _fileName;}
    PUPable_decl(RungOutputParams);
    RungOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputIntParams::pup(p);//Call base class
	}
    };
#endif
