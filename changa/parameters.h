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

#endif
