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
    int nSteps;
    int iStartStep;
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
    int bEwald;
    double dEwCut;
    double dEwhCut;
    double dTheta2;
    int iOrder;
    CSM csm;			/* cosmo parameters */
    double dRedTo;
    /*
     * Units: set by dMsolUnit and dKpcUnit
     */
    double dMsolUnit;
    double dKpcUnit;
    double dGasConst;
    double dErgPerGmUnit;
    double dGmPerCcUnit;
    double dSecUnit;
    double dComovingGmPerCcUnit;
    int bStandard;
    int bOverwrite;
    int bParaRead;
    int bParaWrite;
    char achInFile[256];
    char achOutName[256];
    int bStaticTest;
    int iOutInterval;
    int iCheckInterval;
    int iLogInterval;
    double dExtraStore;		/* Unused, here for PKDGRAV compatibility */
    double dDumpFrameStep;
    double dDumpFrameTime;
    } Parameters;

inline void operator|(PUP::er &p, Parameters &param) {
    p|param.bDoDensity;
    p|param.bDoGravity;
    p|param.dSoft;
    p|param.nSteps;
    p|param.iStartStep;
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
    p|param.bEwald;
    p|param.dEwCut;
    p|param.dEwhCut;
    p|param.dTheta2;
    p|param.iOrder;
    if(p.isUnpacking())
 	csmInitialize(&param.csm);
    p|*param.csm;
    p|param.dRedTo;
    p|param.dMsolUnit;
    p|param.dKpcUnit;
    p|param.dGasConst;
    p|param.dErgPerGmUnit;
    p|param.dGmPerCcUnit;
    p|param.dSecUnit;
    p|param.dComovingGmPerCcUnit;
    p|param.bStandard;
    p|param.bOverwrite;
    p|param.bParaRead;
    p|param.bParaWrite;
    p(param.achInFile, 256);
    p(param.achOutName, 256);
    p|param.bStaticTest;
    p|param.iOutInterval;
    p|param.iCheckInterval;
    p|param.iLogInterval;
    p|param.dExtraStore;
    p|param.dDumpFrameStep;
    p|param.dDumpFrameTime;
    }

#if 0
PUPbytes(Parameters);
#endif
#endif
