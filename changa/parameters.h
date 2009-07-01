#ifndef PARAMETERS_HINCLUDED
#define PARAMETERS_HINCLUDED

#include "cosmo.h"

typedef struct parameters {
    /*
    ** Parameters for ParallelGravity.
    */
    int bDoDensity;
    int bDoGravity;
    double dSoft;
    int bPhysicalSoft;		/* Use physical softening in comoving coords */
    int bSoftMaxMul;		/* dSoftMax is a multiplier, not an
				   absolute limit */
    double dSoftMax;
    int nSteps;
    int iStartStep;
    int iWallRunTime;
    double dDelta;
    int bEpsAccStep;
    int bGravStep;
    double dEta;
    int iMaxRung;
    int bCannonical;
    int bKDK;
    int bPeriodic;
    int nReplicas;
    double fPeriod;
    Vector3D<double> vPeriod;
    int bEwald;
    double dEwCut;
    double dEwhCut;
    double dTheta2;
    int iOrder;
    CSM csm;			/* cosmo parameters */
    double dRedTo;
    /*
     * Gas parameters
     * Units: set by dMsolUnit and dKpcUnit
     */
    int bDoGas;
    int bGeometric;
    int bBulkViscosity;
    int bGasAdiabatic;
    int bGasIsothermal;
    double dhMinOverSoft;
    double dMsolUnit;
    double dKpcUnit;
    double dGasConst;
    double dConstAlpha;
    double dConstBeta;
    double dConstGamma;
    double dMeanMolWeight;
    double dErgPerGmUnit;
    double dGmPerCcUnit;
    double dSecUnit;
    double dComovingGmPerCcUnit;
    int bSphStep;
    int bFastGas;
    double dFracFastGas;
    int bViscosityLimiter;
    int iViscosityLimiter;
    int bViscosityLimitdt;
    double dEtaCourant;
    double dEtauDot;
    int bStandard;
    int bOverwrite;
    int bParaRead;
    int bParaWrite;
    char achInFile[256];
    char achOutName[256];
    int bStaticTest;
    int bBenchmark;
    int iOutInterval;
    int iCheckInterval;
    int iLogInterval;
    double dExtraStore;		/* Unused, here for PKDGRAV compatibility */
    double dDumpFrameStep;
    double dDumpFrameTime;
    int iDirector;
    int bLiveViz;
    } Parameters;

inline void operator|(PUP::er &p, Parameters &param) {
    p|param.bDoDensity;
    p|param.bDoGravity;
    p|param.dSoft;
    p|param.bPhysicalSoft;
    p|param.bSoftMaxMul;
    p|param.dSoftMax;
    p|param.nSteps;
    p|param.iStartStep;
    p|param.iWallRunTime;
    p|param.dDelta;
    p|param.bEpsAccStep;
    p|param.bGravStep;
    p|param.dEta;
    p|param.iMaxRung;
    p|param.bCannonical;
    p|param.bKDK;
    p|param.bPeriodic;
    p|param.nReplicas;
    p|param.fPeriod;
    p|param.vPeriod;
    p|param.bEwald;
    p|param.dEwCut;
    p|param.dEwhCut;
    p|param.dTheta2;
    p|param.iOrder;
    if(p.isUnpacking())
 	csmInitialize(&param.csm);
    p|*param.csm;
    p|param.dRedTo;
    p|param.bDoGas;
    p|param.bGeometric;
    p|param.bBulkViscosity;
    p|param.bGasAdiabatic;
    p|param.bGasIsothermal;
    p|param.bFastGas;
    p|param.dFracFastGas;
    p|param.bViscosityLimiter;
    p|param.iViscosityLimiter;
    p|param.dhMinOverSoft;
    p|param.dMsolUnit;
    p|param.dKpcUnit;
    p|param.dGasConst;
    p|param.dConstAlpha;
    p|param.dConstBeta;
    p|param.dConstGamma;
    p|param.dMeanMolWeight;
    p|param.dErgPerGmUnit;
    p|param.dGmPerCcUnit;
    p|param.dSecUnit;
    p|param.dComovingGmPerCcUnit;
    p|param.bSphStep;
    p|param.bViscosityLimitdt;
    p|param.dEtaCourant;
    p|param.dEtauDot;
    p|param.bStandard;
    p|param.bOverwrite;
    p|param.bParaRead;
    p|param.bParaWrite;
    p(param.achInFile, 256);
    p(param.achOutName, 256);
    p|param.bStaticTest;
    p|param.bBenchmark;
    p|param.iOutInterval;
    p|param.iCheckInterval;
    p|param.iLogInterval;
    p|param.dExtraStore;
    p|param.dDumpFrameStep;
    p|param.dDumpFrameTime;
    p|param.iDirector;
    p|param.bLiveViz;
    }

#endif
