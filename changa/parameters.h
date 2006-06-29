#ifndef PARAMETERS_HINCLUDED
#define PARAMETERS_HINCLUDED

typedef struct parameters {
    /*
    ** Parameters for ParallelGravity.
    */
    double dSoft;
    int nSteps;
    double dTimeStep;
    int bPeriodic;
    int nReplicas;
    double fPeriod;
    int bEwald;
    double dEwCut;
    double dEwhCut;
    char achInFile[256];
    char achOutName[256];
    int bStaticTest;
    int iOutInterval;
    int iLogInterval;
    } Parameters;

#endif
