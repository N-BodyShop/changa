#ifndef PARAMETERS_HINCLUDED
#define PARAMETERS_HINCLUDED

struct parameters {
    /*
    ** Parameters for ParallelGravity.
    */
    double dSoft;
    int nSteps;
    double dTimeStep;
    char achInFile[256];
    int bStaticTest;
    };

#endif
