#ifndef SINKS_HINCLUDED
#define SINKS_HINCLUDED

class Sinks 
{
    //
    // Sink parameters
    //
    int bDoSinks;
    int bBHSink;
    int bSmallBHSmooth;
    int bBHTurnOffCooling;
    int bDoBHKick;
    double dBHSinkEddEff;
    double dBHSinkFeedbackEff;
    double dBHSinkAlpha;
    int bBHMindv;
    int bDoSinksAtStart;
    int bSinkThermal;
    double dSinkRadius;
    double dSinkBoundOrbitRadius;
    double dSinkMustAccreteRadius;
    double dDeltaSink;
    double dSinkMassMin;
    int iSinkRung;
    int bSinkForm;
    int nJeans;
    double dJeansConstant;
    int bSinkFormJeans;
    int bSinkFormDivV;
    double dSinkFormDivVCoeff;
    int bSinkFormDivAcc;
    double dSinkFormDivAccCoeff;
    int bSinkFormDV;
    int bSinkFormPotMin;
    double dSinkFormDenskty;
    double dSinkTimeEligible;
    int bSinkFormSimple;
    int nSinkFormMin;
    inline void pup(PUP::er &p);
    };

inline void Sinks::pup(PUP::er &p) {
    p|bDoSinks;
    p|bBHSink;
    p|bSmallBHSmooth;
    p|bBHTurnOffCooling;
    p|bDoBHKick;
    p|dBHSinkEddEff;
    p|dBHSinkFeedbackEff;
    p|dBHSinkAlpha;
    p|bBHMindv;
    p|bDoSinksAtStart;
    p|bSinkThermal;
    p|dSinkRadius;
    p|dSinkBoundOrbitRadius;
    p|dSinkMustAccreteRadius;
    p|dDeltaSink;
    p|dSinkMassMin;
    p|iSinkRung;
    p|bSinkForm;
    p|nJeans;
    p|dJeansConstant;
    p|bSinkFormJeans;
    p|bSinkFormDivV;
    p|dSinkFormDivVCoeff;
    p|bSinkFormDivAcc;
    p|dSinkFormDivAccCoeff;
    p|bSinkFormDV;
    p|bSinkFormPotMin;
    p|dSinkFormDensity;
    p|dSinkTimeEligible;
    p|bSinkFormSimple;
    p|nSinkFormMin;
    }


#endif
