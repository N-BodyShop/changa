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
    int bStandard;
    int bOverwrite;
    char achInFile[256];
    char achOutName[256];
    int bStaticTest;
    int iOutInterval;
    int iCheckInterval;
    int iLogInterval;
    double dExtraStore;
    } Parameters;

#endif
