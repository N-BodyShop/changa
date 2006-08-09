#ifndef PARAMETERS_HINCLUDED
#define PARAMETERS_HINCLUDED

#include "cosmo.h"

typedef struct parameters {
    /*
    ** Parameters for ParallelGravity.
    */
    double dSoft;
    int nSteps;
    int iStartStep;
    double dTimeStep;
    double dEta;
    int bPeriodic;
    int nReplicas;
    double fPeriod;
    int bEwald;
    double dEwCut;
    double dEwhCut;
    CSM csm;			/* cosmo parameters */
    double dRedTo;
    char achInFile[256];
    char achOutName[256];
    int bStaticTest;
    int iOutInterval;
    int iLogInterval;
    } Parameters;

#endif
