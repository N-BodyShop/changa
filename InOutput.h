/// @file InOutput.h
/// Declarations for I/O implemenatation.
#ifndef __INOUTPUT_H
#define __INOUTPUT_H

class OutputParams;
#include "DataManager.h"

int64_t ncGetCount(std::string typedir);

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
    virtual void setDValue(GravityParticle *p, double) = 0;
    virtual int64_t iValue(GravityParticle *p) = 0;
    virtual void setIValue(GravityParticle *p, int64_t iValue) = 0;
    int bFloat;         // Is a floating point number
    int bVector;	// Is a vector, as opposed to a scalar
    int iBinaryOut;     // Type of binary output
    double dTime;
    std::string fileName;	// output file
    std::string sTipsyExt;      // Extension for tipsy output
    std::string sNChilExt;      // file name for NChilada output
    unsigned int iType;         // mask of families containing this attribute
    unsigned int iTypeWriting;  // family being written in NC format
    DataManager *dm;	// For extra state information (e.g. cooling)

    OutputParams() {dm = NULL;}
    PUPable_abstract(OutputParams);
    OutputParams(CkMigrateMessage *m) : PUP::able(m) {dm = NULL;}
    virtual void pup(PUP::er &p) {
        PUP::able::pup(p);//Call base class
        p|fileName;
        p|sTipsyExt;
        p|sNChilExt;
        p|bFloat;
        p|bVector;
        p|iBinaryOut;
        p|dTime;
        p|iType;
        p|iTypeWriting;
	}
    };

/// @brief Output particle masses
class MassOutputParams : public OutputParams
{
 public:
    virtual double dValue(GravityParticle *p) {return p->mass;}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
    virtual void setDValue(GravityParticle *p, double val) {p->mass = val;}
    virtual int64_t iValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual void setIValue(GravityParticle *p, int64_t iValue) {CkAssert(0);}
    MassOutputParams() {}
    MassOutputParams(std::string _fileName, int _iBinaryOut, double _dTime) {
        bFloat = 1;
        bVector = 0; fileName = _fileName; iBinaryOut = _iBinaryOut;
        sTipsyExt = "mass"; sNChilExt = "mass";
        dTime = _dTime;
        iType = TYPE_GAS | TYPE_DARK | TYPE_STAR; }
    PUPable_decl(MassOutputParams);
    MassOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };

/// @brief Output particle positions
class PosOutputParams : public OutputParams
{
 public:
    virtual double dValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual Vector3D<double> vValue(GravityParticle *p)
				{return p->position;}
    virtual void setDValue(GravityParticle *p, double val) {CkAssert(0);}
    virtual int64_t iValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual void setIValue(GravityParticle *p, int64_t iValue) {CkAssert(0);}
    PosOutputParams() {}
    PosOutputParams(std::string _fileName, int _iBinaryOut, double _dTime) {
        bFloat = 1;
        bVector = 1; fileName = _fileName; iBinaryOut = _iBinaryOut;
        sTipsyExt = "pos"; sNChilExt = "pos";
        dTime = _dTime;
        iType = TYPE_GAS | TYPE_DARK | TYPE_STAR; }
    PUPable_decl(PosOutputParams);
    PosOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };

/// @brief Output particle velocities
class VelOutputParams : public OutputParams
{
 public:
    double dVFac;
    virtual double dValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual Vector3D<double> vValue(GravityParticle *p)
				{return dVFac*p->velocity;}
    virtual void setDValue(GravityParticle *p, double val) {CkAssert(0);}
    virtual int64_t iValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual void setIValue(GravityParticle *p, int64_t iValue) {CkAssert(0);}
    VelOutputParams() {}
    VelOutputParams(std::string _fileName, int _iBinaryOut, double _dTime,
                    double _dVFac) {
        bFloat = 1;
        bVector = 1; fileName = _fileName; iBinaryOut = _iBinaryOut;
        sTipsyExt = "vel"; sNChilExt = "vel";
        dTime = _dTime; dVFac = _dVFac;
        iType = TYPE_GAS | TYPE_DARK | TYPE_STAR; }
    PUPable_decl(VelOutputParams);
    VelOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
        p|dVFac;
	}
    };

/// @brief Output particle gravitational potential
class PotOutputParams : public OutputParams
{
 public:
    virtual double dValue(GravityParticle *p) {return p->potential;}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
    virtual void setDValue(GravityParticle *p, double val) {p->potential = val;}
    virtual int64_t iValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual void setIValue(GravityParticle *p, int64_t iValue) {CkAssert(0);}
    PotOutputParams() {}
    PotOutputParams(std::string _fileName, int _iBinaryOut, double _dTime) {
        bFloat = 1;
        bVector = 0; fileName = _fileName; iBinaryOut = _iBinaryOut;
        sTipsyExt = "pot"; sNChilExt = "pot";
        dTime = _dTime;
        iType = TYPE_GAS | TYPE_DARK | TYPE_STAR; }
    PUPable_decl(PotOutputParams);
    PotOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };

/// @brief Output particle gas density
class GasDenOutputParams : public OutputParams
{
 public:
    virtual double dValue(GravityParticle *p) {return p->fDensity;}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
    virtual void setDValue(GravityParticle *p, double val) {p->fDensity = val;}
    virtual int64_t iValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual void setIValue(GravityParticle *p, int64_t iValue) {CkAssert(0);}
    GasDenOutputParams() {}
    GasDenOutputParams(std::string _fileName, int _iBinaryOut, double _dTime) {
        bFloat = 1;
        bVector = 0; fileName = _fileName; iBinaryOut = _iBinaryOut;
        sTipsyExt = "gasden"; sNChilExt = "GasDensity";
        dTime = _dTime;
        iType = TYPE_GAS; }
    PUPable_decl(GasDenOutputParams);
    GasDenOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };

/// @brief Output particle gas temperature
class TempOutputParams : public OutputParams
{
 public:
    double duTFac;
    bool bGasCooling;
    virtual double dValue(GravityParticle *p) {
        if(bGasCooling) {
#ifndef COOLING_NONE
#ifdef COOLING_GRACKLE
            return CoolCodeEnergyToTemperature(dm->Cool, &p->CoolParticle(),
                p->u(), p->fDensity, p->fMetals());
#else
            return CoolCodeEnergyToTemperature(dm->Cool, &p->CoolParticle(),
                                               p->u(), p->fMetals());
#endif
#else
            CkAssert(0);
#endif
            }
        return duTFac*p->u();
        }
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
    virtual void setDValue(GravityParticle *p, double val) {CkAssert(0);}
    virtual int64_t iValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual void setIValue(GravityParticle *p, int64_t iValue) {CkAssert(0);}
    TempOutputParams() {}
    TempOutputParams(std::string _fileName, int _iBinaryOut, double _dTime,
                     bool _bGasCooling, double _duTFac) {
        bFloat = 1;
        bVector = 0; fileName = _fileName; iBinaryOut = _iBinaryOut;
        sTipsyExt = "temp"; sNChilExt = "temperature";
        dTime = _dTime; bGasCooling = _bGasCooling; duTFac = _duTFac;
        iType = TYPE_GAS; }
    PUPable_decl(TempOutputParams);
    TempOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
        p|bGasCooling;
        p|duTFac;
	}
    };


/// @brief Output accelerations.
class AccOutputParams : public OutputParams
{
 public:
    virtual double dValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual Vector3D<double> vValue(GravityParticle *p)
				{return p->treeAcceleration;}
    virtual void setDValue(GravityParticle *p, double val) {CkAssert(0);}
    virtual int64_t iValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual void setIValue(GravityParticle *p, int64_t iValue) {CkAssert(0);}
    AccOutputParams() {}
    AccOutputParams(std::string _fileName) { bFloat = 1; bVector = 1; fileName = _fileName;}
    AccOutputParams(std::string _fileName, int _iBinaryOut, double _dTime) {
        bFloat = 1;
        bVector = 1; fileName = _fileName; iBinaryOut = _iBinaryOut;
        sTipsyExt = "acc2"; sNChilExt = "acc";
        dTime = _dTime;
        iType = TYPE_GAS | TYPE_DARK | TYPE_STAR; }
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
    virtual void setDValue(GravityParticle *p, double val) {p->fDensity = val;}
    virtual int64_t iValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual void setIValue(GravityParticle *p, int64_t iValue) {CkAssert(0);}
 public:
    DenOutputParams() {}
    DenOutputParams(std::string _fileName) { bFloat = 1; bVector = 0; fileName = _fileName;}
    DenOutputParams(std::string _fileName, int _iBinaryOut, double _dTime) {
        bFloat = 1;
        bVector = 0; fileName = _fileName; iBinaryOut = _iBinaryOut;
        sTipsyExt = "den"; sNChilExt = "den";
        dTime = _dTime;
        iType = TYPE_GAS | TYPE_DARK | TYPE_STAR; }
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
    virtual void setDValue(GravityParticle *p, double val) {p->fBall = 2.0*val;}
    virtual int64_t iValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual void setIValue(GravityParticle *p, int64_t iValue) {CkAssert(0);}
 public:
    HsmOutputParams() {}
    HsmOutputParams(std::string _fileName, int _iBinaryOut, double _dTime) {
        bFloat = 1;
        bVector = 0; fileName = _fileName; iBinaryOut = _iBinaryOut;
        sTipsyExt = "smoothlength"; sNChilExt = "smoothlength";
        dTime = _dTime;
        iType = TYPE_GAS | TYPE_DARK | TYPE_STAR; }
    PUPable_decl(HsmOutputParams);
    HsmOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };

/// @brief Output particle gravitational softenings
class SoftOutputParams : public OutputParams
{
 public:
    virtual double dValue(GravityParticle *p) {
#ifdef CHANGESOFT
        return p->fSoft0;
#else
        return p->soft;
#endif
    }
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
    virtual void setDValue(GravityParticle *p, double val) {
#ifdef CHANGESOFT
	p->fSoft0 = val;
#else
	p->soft = val;
#endif
	}
    virtual int64_t iValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual void setIValue(GravityParticle *p, int64_t iValue) {CkAssert(0);}
    SoftOutputParams() {}
    SoftOutputParams(std::string _fileName, int _iBinaryOut, double _dTime) {
        bFloat = 1;
        bVector = 0; fileName = _fileName; iBinaryOut = _iBinaryOut;
        sTipsyExt = "soft"; sNChilExt = "soft";
        dTime = _dTime;
        iType = TYPE_GAS | TYPE_DARK | TYPE_STAR; }
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
    virtual void setDValue(GravityParticle *p, double val) {CkAssert(0);}
    virtual int64_t iValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual void setIValue(GravityParticle *p, int64_t iValue) {CkAssert(0);}
 public:
    PresOutputParams() {}
    PresOutputParams(std::string _fileName) { bFloat = 1; bVector = 0; fileName = _fileName;}
    PresOutputParams(std::string _fileName, int _iBinaryOut, double _dTime) {
        bFloat = 1;
        bVector = 0; fileName = _fileName; iBinaryOut = _iBinaryOut;
        sTipsyExt = "pres"; sNChilExt = "pres";
        dTime = _dTime;
        iType = TYPE_GAS; }
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
    virtual void setDValue(GravityParticle *p, double val) {p->divv() = val;}
    virtual int64_t iValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual void setIValue(GravityParticle *p, int64_t iValue) {CkAssert(0);}
 public:
    DivVOutputParams() {}
    DivVOutputParams(std::string _fileName, int _iBinaryOut, double _dTime) {
        bFloat = 1;
        bVector = 0; fileName = _fileName; iBinaryOut = _iBinaryOut;
        sTipsyExt = "divv"; sNChilExt = "divv";
        dTime = _dTime;
        iType = TYPE_GAS; }
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
    virtual void setDValue(GravityParticle *p, double val) {p->PdV() = val;}
    virtual int64_t iValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual void setIValue(GravityParticle *p, int64_t iValue) {CkAssert(0);}
 public:
    PDVOutputParams() {}
    PDVOutputParams(std::string _fileName) { bFloat = 1; bVector = 0; fileName = _fileName;}
    PDVOutputParams(std::string _fileName, int _iBinaryOut, double _dTime) {
        bFloat = 1;
        bVector = 0; fileName = _fileName; iBinaryOut = _iBinaryOut;
        sTipsyExt = "PdV"; sNChilExt = "PdV";
        dTime = _dTime;
        iType = TYPE_GAS; }
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
    virtual void setDValue(GravityParticle *p, double val) {p->mumax() = val;}
    virtual int64_t iValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual void setIValue(GravityParticle *p, int64_t iValue) {CkAssert(0);}
 public:
    MuMaxOutputParams() {}
    MuMaxOutputParams(std::string _fileName) { bFloat = 1; bVector = 0; fileName = _fileName;}
    MuMaxOutputParams(std::string _fileName, int _iBinaryOut, double _dTime) {
        bFloat = 1;
        bVector = 0; fileName = _fileName; iBinaryOut = _iBinaryOut;
        sTipsyExt = "mumax"; sNChilExt = "mumax";
        dTime = _dTime;
        iType = TYPE_GAS; }
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
    virtual void setDValue(GravityParticle *p, double val) {p->BalsaraSwitch() = val;}
    virtual int64_t iValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual void setIValue(GravityParticle *p, int64_t iValue) {CkAssert(0);}
 public:
    BSwOutputParams() {}
    BSwOutputParams(std::string _fileName, int _iBinaryOut, double _dTime) {
        bFloat = 1;
        bVector = 0; fileName = _fileName; iBinaryOut = _iBinaryOut;
        sTipsyExt = "BSw"; sNChilExt = "BSw";
        dTime = _dTime;
        iType = TYPE_GAS; }
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
    virtual void setDValue(GravityParticle *p, double val) {p->c() = val;}
    virtual int64_t iValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual void setIValue(GravityParticle *p, int64_t iValue) {CkAssert(0);}
 public:
    CsOutputParams() {}
    CsOutputParams(std::string _fileName, int _iBinaryOut, double _dTime) {
        bFloat = 1;
        bVector = 0; fileName = _fileName; iBinaryOut = _iBinaryOut;
        sTipsyExt = "c"; sNChilExt = "c";
        dTime = _dTime;
        iType = TYPE_GAS; }
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
    virtual void setDValue(GravityParticle *p, double val) {CkAssert(0);}
    virtual int64_t iValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual void setIValue(GravityParticle *p, int64_t iValue) {CkAssert(0);}
 public:
    EDotOutputParams() {}
    EDotOutputParams(std::string _fileName) { bFloat = 1; bVector = 0; fileName = _fileName;}
    EDotOutputParams(std::string _fileName, int _iBinaryOut, double _dTime) {
        bFloat = 1;
        bVector = 0; fileName = _fileName; iBinaryOut = _iBinaryOut;
        sTipsyExt = "eDot"; sNChilExt = "eDot";
        dTime = _dTime;
        iType = TYPE_GAS; }
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
	    return COOL_ARRAY0(dm->Cool, &p->CoolParticle(), p->fMetals());
	else
#endif
	    return 0.0;
	}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
    virtual void setDValue(GravityParticle *p, double val) {
#ifndef COOLING_NONE
	if (TYPETest(p, TYPE_GAS))
	    COOL_SET_ARRAY0(dm->Cool, &p->CoolParticle(), p->fMetals(), val);
#endif
	}
    virtual int64_t iValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual void setIValue(GravityParticle *p, int64_t iValue) {CkAssert(0);}
 public:
    Cool0OutputParams() {}
    Cool0OutputParams(std::string _fileName, int _iBinaryOut, double _dTime) {
        bFloat = 1;
        bVector = 0; fileName = _fileName; iBinaryOut = _iBinaryOut;
#ifndef COOLING_NONE
        sTipsyExt = COOL_ARRAY0_EXT; sNChilExt = COOL_ARRAY0_EXT;
#endif
        dTime = _dTime;
        iType = TYPE_GAS; }
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
	    return COOL_ARRAY1(dm->Cool, &p->CoolParticle(), p->fMetals());
	else
#endif
	    return 0.0;
	}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
    virtual void setDValue(GravityParticle *p, double val) {
#ifndef COOLING_NONE
        //XXX be sure metals has been set!
	if (TYPETest(p, TYPE_GAS))
	    COOL_SET_ARRAY1(dm->Cool, &p->CoolParticle(), p->fMetals(), val);
#endif
	}
    virtual int64_t iValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual void setIValue(GravityParticle *p, int64_t iValue) {CkAssert(0);}
 public:
    Cool1OutputParams() {}
    Cool1OutputParams(std::string _fileName, int _iBinaryOut, double _dTime) {
        bFloat = 1;
        bVector = 0; fileName = _fileName; iBinaryOut = _iBinaryOut;
#ifndef COOLING_NONE
        sTipsyExt = COOL_ARRAY1_EXT; sNChilExt = COOL_ARRAY1_EXT;
#endif
        dTime = _dTime;
        iType = TYPE_GAS; }
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
	    return COOL_ARRAY2(dm->Cool, &p->CoolParticle(), p->fMetals());
	else
#endif
	    return 0.0;
	}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
    virtual void setDValue(GravityParticle *p, double val) {
#ifndef COOLING_NONE
	if (TYPETest(p, TYPE_GAS))
	    COOL_SET_ARRAY2(dm->Cool, &p->CoolParticle(), p->fMetals(), val);
#endif
	}
    virtual int64_t iValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual void setIValue(GravityParticle *p, int64_t iValue) {CkAssert(0);}
 public:
    Cool2OutputParams() {}
    Cool2OutputParams(std::string _fileName, int _iBinaryOut, double _dTime) {
        bFloat = 1;
        bVector = 0; fileName = _fileName; iBinaryOut = _iBinaryOut;
#ifndef COOLING_NONE
        sTipsyExt = COOL_ARRAY2_EXT; sNChilExt = COOL_ARRAY2_EXT;
#endif
        dTime = _dTime;
        iType = TYPE_GAS; }
    PUPable_decl(Cool2OutputParams);
    Cool2OutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };

/// @brief Output the value in cool_array3.
class Cool3OutputParams : public OutputParams
{
    virtual double dValue(GravityParticle *p)
    {
#ifndef COOLING_NONE
	if (TYPETest(p, TYPE_GAS))
	    return COOL_ARRAY3(dm->Cool, &p->CoolParticle(), p->fMetals());
	else
#endif
	    return 0.0;
	}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
    virtual void setDValue(GravityParticle *p, double val) {
#ifndef COOLING_NONE
	if (TYPETest(p, TYPE_GAS))
	    CkAssert(0);
#endif
	}
    virtual int64_t iValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual void setIValue(GravityParticle *p, int64_t iValue) {CkAssert(0);}
 public:
    Cool3OutputParams() {}
    Cool3OutputParams(std::string _fileName, int _iBinaryOut, double _dTime) {
        bFloat = 1;
        bVector = 0; fileName = _fileName; iBinaryOut = _iBinaryOut;
#ifndef COOLING_NONE
        sTipsyExt = COOL_ARRAY3_EXT; sNChilExt = COOL_ARRAY3_EXT;
#endif
        dTime = _dTime;
        iType = TYPE_GAS; }
    PUPable_decl(Cool3OutputParams);
    Cool3OutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };

/// @brief Output Oxygen mass fraction.
class OxOutputParams : public OutputParams
{
    virtual double dValue(GravityParticle *p) {
	if (TYPETest(p, TYPE_STAR)) return p->fStarMFracOxygen();
	if (TYPETest(p, TYPE_GAS)) return p->fMFracOxygen();
	else return 0.0;
	}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
    virtual void setDValue(GravityParticle *p, double val) {
	if (TYPETest(p, TYPE_STAR)) p->fStarMFracOxygen() = val;
	if (TYPETest(p, TYPE_GAS)) p->fMFracOxygen() = val;
	}
    virtual int64_t iValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual void setIValue(GravityParticle *p, int64_t iValue) {CkAssert(0);}
 public:
    OxOutputParams() {}
    OxOutputParams(std::string _fileName, int _iBinaryOut, double _dTime) {
        bFloat = 1;
        bVector = 0; fileName = _fileName; iBinaryOut = _iBinaryOut;
        sTipsyExt = "OxMassFrac"; sNChilExt = "OxMassFrac";
        dTime = _dTime;
        iType = TYPE_GAS | TYPE_STAR; }
    PUPable_decl(OxOutputParams);
    OxOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };

/// @brief Output Iron mass fraction.
class FeOutputParams : public OutputParams
{
    virtual double dValue(GravityParticle *p) {
	if (TYPETest(p, TYPE_STAR)) return p->fStarMFracIron();
	if (TYPETest(p, TYPE_GAS)) return p->fMFracIron();
	else return 0.0;
	}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
    virtual void setDValue(GravityParticle *p, double val) {
	if (TYPETest(p, TYPE_STAR)) p->fStarMFracIron() = val;
	if (TYPETest(p, TYPE_GAS)) p->fMFracIron() = val;
	}
    virtual int64_t iValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual void setIValue(GravityParticle *p, int64_t iValue) {CkAssert(0);}
 public:
    FeOutputParams() {}
    FeOutputParams(std::string _fileName, int _iBinaryOut, double _dTime) {
        bFloat = 1;
        bVector = 0; fileName = _fileName; iBinaryOut = _iBinaryOut;
        sTipsyExt = "FeMassFrac"; sNChilExt = "FeMassFrac";
        dTime = _dTime;
        iType = TYPE_GAS | TYPE_STAR; }
    PUPable_decl(FeOutputParams);
    FeOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };

/// @brief Output metal mass fraction.
class MetalsOutputParams : public OutputParams
{
    virtual double dValue(GravityParticle *p) {
	if (TYPETest(p, TYPE_STAR)) return p->fStarMetals();
	if (TYPETest(p, TYPE_GAS)) return p->fMetals();
	else return 0.0;
	}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
    virtual void setDValue(GravityParticle *p, double val) {
	if (TYPETest(p, TYPE_STAR)) p->fStarMetals() = val;
	if (TYPETest(p, TYPE_GAS)) p->fMetals() = val;
        }
    virtual int64_t iValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual void setIValue(GravityParticle *p, int64_t iValue) {CkAssert(0);}
 public:
    MetalsOutputParams() {}
    MetalsOutputParams(std::string _fileName, int _iBinaryOut, double _dTime) {
        bFloat = 1;
        bVector = 0; fileName = _fileName; iBinaryOut = _iBinaryOut;
        sTipsyExt = "metals"; sNChilExt = "metals";
        dTime = _dTime;
        iType = TYPE_GAS | TYPE_STAR; }
    PUPable_decl(MetalsOutputParams);
    MetalsOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };

/// @brief Output mass at formation time.
class MFormOutputParams : public OutputParams
{
    virtual double dValue(GravityParticle *p) {
	if (TYPETest(p, TYPE_STAR)) return p->fMassForm();
	else return 0.0;
	}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
    virtual void setDValue(GravityParticle *p, double val) {
	if (TYPETest(p, TYPE_STAR)) p->fMassForm() = val;
	}
    virtual int64_t iValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual void setIValue(GravityParticle *p, int64_t iValue) {CkAssert(0);}
 public:
    MFormOutputParams() {}
    MFormOutputParams(std::string _fileName, int _iBinaryOut, double _dTime) {
        bFloat = 1;
        bVector = 0; fileName = _fileName; iBinaryOut = _iBinaryOut;
        sTipsyExt = "massform"; sNChilExt = "massform";
        dTime = _dTime;
        iType = TYPE_STAR; }
    PUPable_decl(MFormOutputParams);
    MFormOutputParams(CkMigrateMessage *m) {}
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
    virtual void setDValue(GravityParticle *p, double val) {
#ifdef DIFFUSION
	if (p->isGas()) p->fMetalsDot() = val;
#endif
	}
    virtual int64_t iValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual void setIValue(GravityParticle *p, int64_t iValue) {CkAssert(0);}
 public:
    MetalsDotOutputParams() {}
    MetalsDotOutputParams(std::string _fileName, int _iBinaryOut, double _dTime) {
        bFloat = 1;
        bVector = 0; fileName = _fileName; iBinaryOut = _iBinaryOut;
        sTipsyExt = "Metalsdot"; sNChilExt = "Metalsdot";
        dTime = _dTime;
        iType = TYPE_GAS; }
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
    virtual void setDValue(GravityParticle *p, double val) {
#ifdef DIFFUSION
	if (p->isGas()) p->fMFracOxygenDot() = val;
#endif
	}
    virtual int64_t iValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual void setIValue(GravityParticle *p, int64_t iValue) {CkAssert(0);}
 public:
    OxygenMassFracDotOutputParams() {}
    OxygenMassFracDotOutputParams(std::string _fileName, int _iBinaryOut, double _dTime) {
        bFloat = 1;
        bVector = 0; fileName = _fileName; iBinaryOut = _iBinaryOut;
        sTipsyExt = "OxMassFracdot"; sNChilExt = "OxMassFracdot";
        dTime = _dTime;
        iType = TYPE_GAS; }
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
    virtual void setDValue(GravityParticle *p, double val) {
#ifdef DIFFUSION
	if (p->isGas()) p->fMFracIronDot() = val;
#endif
	}
    virtual int64_t iValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual void setIValue(GravityParticle *p, int64_t iValue) {CkAssert(0);}
 public:
    IronMassFracDotOutputParams() {}
    IronMassFracDotOutputParams(std::string _fileName, int _iBinaryOut, double _dTime) {
        bFloat = 1;
        bVector = 0; fileName = _fileName; iBinaryOut = _iBinaryOut;
        sTipsyExt = "FeMassFracdot"; sNChilExt = "FeMassFracdot";
        dTime = _dTime;
        iType = TYPE_GAS; }
    PUPable_decl(IronMassFracDotOutputParams);
    IronMassFracDotOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };

/// @brief Output formation time.
class TimeFormOutputParams : public OutputParams
{
    virtual double dValue(GravityParticle *p) {
	if (TYPETest(p, TYPE_STAR)) return p->fTimeForm();
	else return 0.0;
	}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
    virtual void setDValue(GravityParticle *p, double val) {
	if (TYPETest(p, TYPE_STAR)) p->fTimeForm() = val;
        }
    virtual int64_t iValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual void setIValue(GravityParticle *p, int64_t iValue) {CkAssert(0);}
 public:
    TimeFormOutputParams() {}
    TimeFormOutputParams(std::string _fileName, int _iBinaryOut, double _dTime) {
        bFloat = 1;
        bVector = 0; fileName = _fileName; iBinaryOut = _iBinaryOut;
        sTipsyExt = "timeform"; sNChilExt = "timeform";
        dTime = _dTime;
        iType = TYPE_STAR; }
    PUPable_decl(TimeFormOutputParams);
    TimeFormOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };

class AgeOutputParams : public OutputParams
{
    virtual double dValue(GravityParticle *p) {
	if (TYPETest(p, TYPE_STAR)) return dTime - p->fTimeForm();
	else return 0.0;
	}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
    virtual void setDValue(GravityParticle *p, double val) {
	if (TYPETest(p, TYPE_STAR)) p->fTimeForm() = dTime - val;
        }
    virtual int64_t iValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual void setIValue(GravityParticle *p, int64_t iValue) {CkAssert(0);}
 public:
    AgeOutputParams() {}
    AgeOutputParams(std::string _fileName, int _iBinaryOut, double _dTime) {
        bVector = 0; fileName = _fileName; iBinaryOut = _iBinaryOut;
        sTipsyExt = "age"; sNChilExt = "age";
        dTime = _dTime;
        iType = TYPE_STAR; }
    PUPable_decl(AgeOutputParams);
    AgeOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };


/// @brief Output "cool on time" (time cooling is off until)
class coolontimeOutputParams : public OutputParams
{
    virtual double dValue(GravityParticle *p) {
	if (TYPETest(p, TYPE_GAS)) return p->fTimeCoolIsOffUntil();
	else return 0.0;
	}
    virtual Vector3D<double> vValue(GravityParticle *p)
    {CkAssert(0); return 0.0;}
    virtual void setDValue(GravityParticle *p, double val) {
	if (TYPETest(p, TYPE_GAS)) p->fTimeCoolIsOffUntil() = val;
	}
    virtual int64_t iValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual void setIValue(GravityParticle *p, int64_t iValue) {CkAssert(0);}
 public:
    coolontimeOutputParams() {}
    coolontimeOutputParams(std::string _fileName, int _iBinaryOut, double _dTime) {
        bFloat = 1;
        bVector = 0; fileName = _fileName; iBinaryOut = _iBinaryOut;
        sTipsyExt = "coolontime"; sNChilExt = "coolontime";
        dTime = _dTime;
        iType = TYPE_GAS; }
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
    virtual void setDValue(GravityParticle *p, double val) {
	if (p->isGas()) p->fESNrate() = val;
	else if(p->isStar()) p->fStarESNrate() = val;
	}
    virtual int64_t iValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual void setIValue(GravityParticle *p, int64_t iValue) {CkAssert(0);}
 public:
    ESNRateOutputParams() {}
    ESNRateOutputParams(std::string _fileName, int _iBinaryOut, double _dTime) {
        bFloat = 1;
        bVector = 0; fileName = _fileName; iBinaryOut = _iBinaryOut;
        sTipsyExt = "ESNRate"; sNChilExt = "ESNRate";
        dTime = _dTime;
        iType = TYPE_GAS | TYPE_STAR; }
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
    virtual void setDValue(GravityParticle *p, double val) {
#ifdef NEED_DT
	p->dt = val;
#endif
	}
    virtual int64_t iValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual void setIValue(GravityParticle *p, int64_t iValue) {CkAssert(0);}
 public:
    DtOutputParams() {}
    DtOutputParams(std::string _fileName, int _iBinaryOut, double _dTime) {
        bFloat = 1;
        bVector = 0; fileName = _fileName; iBinaryOut = _iBinaryOut;
        sTipsyExt = "dt"; sNChilExt = "dt";
        dTime = _dTime;
        iType = TYPE_GAS | TYPE_DARK | TYPE_STAR; }
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
    virtual void setDValue(GravityParticle *p, double val) {CkAssert(0);}
    virtual int64_t iValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual void setIValue(GravityParticle *p, int64_t iValue) {CkAssert(0);}
 public:
    KeyOutputParams() {}
    KeyOutputParams(std::string _fileName, int _iBinaryOut, double _dTime) {
        bFloat = 1;
        bVector = 0; fileName = _fileName; iBinaryOut = _iBinaryOut;
        sTipsyExt = "key"; sNChilExt = "key";
        dTime = _dTime;
        iType = TYPE_GAS | TYPE_DARK | TYPE_STAR; }
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
    virtual void setDValue(GravityParticle *p, double val) {CkAssert(0);}
    virtual int64_t iValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual void setIValue(GravityParticle *p, int64_t iValue) {CkAssert(0);}
 public:
    DomainOutputParams() {}
    DomainOutputParams(std::string _fileName) { bFloat = 1; bVector = 0; fileName = _fileName;}
    DomainOutputParams(std::string _fileName, int _iBinaryOut, double _dTime) {
        bFloat = 1;
        bVector = 0; fileName = _fileName; iBinaryOut = _iBinaryOut;
        sTipsyExt = "dom"; sNChilExt = "dom";
        dTime = _dTime;
        iType = TYPE_GAS | TYPE_DARK | TYPE_STAR; }
    PUPable_decl(DomainOutputParams);
    DomainOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };

/// @brief Output iOrder.
class IOrderOutputParams : public OutputParams
{
    virtual int64_t iValue(GravityParticle *p)
    {
	return p->iOrder;
	}
    virtual void setIValue(GravityParticle *p, int64_t iValue)
    {
	p->iOrder = iValue;
	}
    virtual double dValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
    virtual void setDValue(GravityParticle *p, double val) {CkAssert(0);}
 public:
    IOrderOutputParams() {}
    IOrderOutputParams(std::string _fileName, int _iBinaryOut, double _dTime) {
        bFloat = 0;
        bVector = 0;
        fileName = _fileName; iBinaryOut = _iBinaryOut;
        sTipsyExt = "iord"; sNChilExt = "iord";
        dTime = _dTime;
        iType = TYPE_GAS | TYPE_DARK | TYPE_STAR; }
    PUPable_decl(IOrderOutputParams);
    IOrderOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };

/// @brief Output iGasOrder.
class IGasOrderOutputParams : public OutputParams
{
    virtual int64_t iValue(GravityParticle *p)
    {
	if(p->isStar())
	    return p->iGasOrder();
	else
	    return 0;
	}
    virtual void setIValue(GravityParticle *p, int64_t iValue)
    {
        if(p->isStar()) p->iGasOrder() = iValue;
	}
    virtual double dValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
    virtual void setDValue(GravityParticle *p, double val) {CkAssert(0);}
 public:
    IGasOrderOutputParams() {}
    IGasOrderOutputParams(std::string _fileName, int _iBinaryOut, double _dTime) {
        bFloat = 0;
        bVector = 0;
        fileName = _fileName; iBinaryOut = _iBinaryOut;
        sTipsyExt = "igasorder"; sNChilExt = "igasorder";
        dTime = _dTime;
        iType = TYPE_STAR; }
    PUPable_decl(IGasOrderOutputParams);
    IGasOrderOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };

/// @brief Output rungs.
class RungOutputParams : public OutputParams
{
    virtual int64_t iValue(GravityParticle *p)
    {
	return p->rung;
	}
    virtual void setIValue(GravityParticle *p, int64_t iValue)
    {
        p->rung = iValue;
	}
    virtual double dValue(GravityParticle *p) {CkAssert(0); return 0.0;}
    virtual Vector3D<double> vValue(GravityParticle *p)
			    {CkAssert(0); return 0.0;}
    virtual void setDValue(GravityParticle *p, double val) {CkAssert(0);}
 public:
    RungOutputParams() {}
    RungOutputParams(std::string _fileName, int _iBinaryOut, double _dTime) {
        bFloat = 0;
        bVector = 0;
        fileName = _fileName; iBinaryOut = _iBinaryOut;
        sTipsyExt = "rung"; sNChilExt = "rung";
        dTime = _dTime;
        iType = TYPE_GAS| TYPE_DARK | TYPE_STAR; }
    PUPable_decl(RungOutputParams);
    RungOutputParams(CkMigrateMessage *m) {}
    virtual void pup(PUP::er &p) {
        OutputParams::pup(p);//Call base class
	}
    };
#endif
