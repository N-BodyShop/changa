#ifndef PARAMETERS_HINCLUDED
#define PARAMETERS_HINCLUDED

struct parameters {
    /*
    ** Parameters for ParallelGravity.
    */
    double dSoft;
    int nSteps;
    double dTimeStep;
    int bPeriodic;
    int nReplicas;
    char achInFile[256];
    char achOutName[256];
    int bStaticTest;
    int iOutInterval;
    int iLogInterval;
    };

#endif
