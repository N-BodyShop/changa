#ifndef PARAMETERS_HINCLUDED
#define PARAMETERS_HINCLUDED

#include "cosmo.h"
#include "cooling.h"
#include "starform.h"
#include "feedback.h"
#include "sinks.h"

#include "externalGravity.h"

/** @brief Hold parameters of the run.
 */
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
    int nTruncateRung;
    int iMaxRung;
    int bCannonical;
    int bKDK;
    int bDtAdjust;
    int bPeriodic;
    int nReplicas;
    double fPeriod;
    Vector3D<double> vPeriod;
    int bEwald;
    double dEwCut;
    double dEwhCut;
    double dTheta;
    double dTheta2;
    double daSwitchTheta;
    int iOrder;
    int bConcurrentSph;
    double dFracNoDomainDecomp;
#ifdef PUSH_GRAVITY
    double dFracPushParticles;
#endif
    CSM csm;			/* cosmo parameters */
    double dRedTo;
    double dGlassDamper;
    /*
     * GrowMass parameters
     */
    int bDynGrowMass;
    int nGrowMass;
    double dGrowDeltaM;
    double dGrowStartT;
    double dGrowEndT;
    /*
     * Gas parameters
     * Units: set by dMsolUnit and dKpcUnit
     */
    int bDoGas;
    int bGeometric;
    int bBulkViscosity;
    int bGasAdiabatic;
    int bGasIsothermal;
    int bGasCooling;
    int nSmooth;
    COOLPARAM CoolParam;
    double dMaxEnergy;
    double dhMinOverSoft;
    double dResolveJeans;
    double dMsolUnit;
    double dKpcUnit;
    double ddHonHLimit;
    double dGasConst;
    double dConstAlpha;
    double dConstBeta;
    double dConstAlphaMax;
    double dConstGamma;
    double dMeanMolWeight;
    double dErgPerGmUnit;
    double dGmPerCcUnit;
    double dSecUnit;
    double dComovingGmPerCcUnit;
    double dThermalDiffusionCoeff;
    double dMetalDiffusionCoeff;
    int bConstantDiffusion;
    int bSphStep;
    int bFastGas;
    double dFracFastGas;
    int bViscosityLimiter;
    int iViscosityLimiter;
    int bViscosityLimitdt;
    double dEtaCourant;
    double dEtaDiffusion;
    double dEtauDot;
    int bStarForm;
    Stfm *stfm;
    int bFeedback;
    Fdbk *feedback;
#ifdef SUPERBUBBLE
    double dThermalCondCoeff;
    double dThermalCondSatCoeff;
    double dThermalCond2Coeff;
    double dThermalCond2SatCoeff;
    double dEvapMinTemp;
    double dEvapCoeff;
#endif
    int bDoExternalGravity;
    ExternalGravity externalGravity;
    int iRandomSeed;
    
    Sinks sinks;

    //SIDM
    int iSIDMSelect;
    double dSIDMSigma;
    double dSIDMVariable;
    
    //
    // Output parameters
    //
    int bStandard;
    int bDoublePos;
    int bDoubleVel;
    int bOverwrite;
    int bParaRead;
    int bParaWrite;
    int nIOProcessor;
    char achInFile[256];
    char achOutName[256];
    int bStaticTest;
    int bBenchmark;
    int iBinaryOut;
    int iOutInterval;
    int iCheckInterval;
    int iLogInterval;
    int iOrbitOutInterval;
    int bDoIOrderOutput;
    int bDoSoftOutput;
    int bDohOutput;
    int bDoCSound;
    int bDoStellarLW;
    int cacheLineDepth;
    double dExtraStore;
    double dMaxBalance;
    double dFracLoadBalance;
    double dDumpFrameStep;
    double dDumpFrameTime;
    int iDirector;
    int bLiveViz;
    int bUseCkLoopPar;
    int iVerbosity;
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
    p|param.nTruncateRung;
    p|param.iMaxRung;
    p|param.bCannonical;
    p|param.bKDK;
    p|param.bDtAdjust;
    p|param.bPeriodic;
    p|param.nReplicas;
    p|param.fPeriod;
    p|param.vPeriod;
    p|param.bEwald;
    p|param.dEwCut;
    p|param.dEwhCut;
    p|param.dTheta;
    p|param.dTheta2;
    p|param.daSwitchTheta;
    p|param.iOrder;
    p|param.bConcurrentSph;
    p|param.dFracNoDomainDecomp;
#ifdef PUSH_GRAVITY
    p|param.dFracPushParticles;
#endif
    if(p.isUnpacking())
 	csmInitialize(&param.csm);
    p|*param.csm;
    p|param.dGlassDamper;
    p|param.dRedTo;
    p|param.bDynGrowMass;
    p|param.nGrowMass;
    p|param.dGrowDeltaM;
    p|param.dGrowStartT;
    p|param.dGrowEndT;
    p|param.bDoGas;
    p|param.bGeometric;
    p|param.bBulkViscosity;
    p|param.bGasAdiabatic;
    p|param.bGasIsothermal;
    p|param.bGasCooling;
    p|param.nSmooth;
    p((char *)&param.CoolParam, sizeof(param.CoolParam));
    p|param.bFastGas;
    p|param.dFracFastGas;
    p|param.bViscosityLimiter;
    p|param.iViscosityLimiter;
    p|param.dMaxEnergy;
    p|param.dhMinOverSoft;
    p|param.dResolveJeans;
    p|param.dMsolUnit;
    p|param.dKpcUnit;
    p|param.ddHonHLimit;
    p|param.dGasConst;
    p|param.dConstAlpha;
    p|param.dConstBeta;
    p|param.dConstAlphaMax;
    p|param.dConstGamma;
    p|param.dMeanMolWeight;
    p|param.dErgPerGmUnit;
    p|param.dGmPerCcUnit;
    p|param.dSecUnit;
    p|param.dComovingGmPerCcUnit;
    p|param.dThermalDiffusionCoeff;
    p|param.dMetalDiffusionCoeff;
    p|param.bConstantDiffusion;
    p|param.bSphStep;
    p|param.bViscosityLimitdt;
    p|param.dEtaCourant;
    p|param.dEtaDiffusion;
    p|param.dEtauDot;
    p|param.bStarForm;
    if(p.isUnpacking())
 	param.stfm = new Stfm();
    p|*param.stfm;
    p|param.bFeedback;
    p|param.feedback;
#ifdef SUPERBUBBLE
    p|param.dThermalCondCoeff;
    p|param.dThermalCondSatCoeff;
    p|param.dThermalCond2Coeff;
    p|param.dThermalCond2SatCoeff;
    p|param.dEvapMinTemp;
    p|param.dEvapCoeff;
#endif
    p|param.bDoExternalGravity;
    p|param.externalGravity;
    p|param.iRandomSeed;
    p|param.sinks;
    p|param.dSIDMSigma;
    p|param.iSIDMSelect;
    p|param.dSIDMVariable;
    p|param.bStandard;
    p|param.bDoublePos;
    p|param.bDoubleVel;
    p|param.bOverwrite;
    p|param.bParaRead;
    p|param.bParaWrite;
    p|param.nIOProcessor;
    p(param.achInFile, 256);
    p(param.achOutName, 256);
    p|param.bStaticTest;
    p|param.bBenchmark;
    p|param.iBinaryOut;
    p|param.iOutInterval;
    p|param.iCheckInterval;
    p|param.iLogInterval;
    p|param.iOrbitOutInterval;
    p|param.bDoIOrderOutput;
    p|param.bDoSoftOutput;
    p|param.bDohOutput;
    p|param.bDoCSound;
    p|param.bDoStellarLW;
    p|param.cacheLineDepth;
    p|param.dExtraStore;
    p|param.dMaxBalance;
    p|param.dFracLoadBalance;
    p|param.dDumpFrameStep;
    p|param.dDumpFrameTime;
    p|param.iDirector;
    p|param.bLiveViz;
    p|param.bUseCkLoopPar;
    p|param.iVerbosity;
    }

#endif
