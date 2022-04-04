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
    double dDeltaMin;
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
    double dVoroKick;
    double dVoroMesh;
    double dHydroTheta;
    int  iRiemannSolver;
    double dSmallRho;
    double dBx, dBy, dBz;
    double dVPsi;
    int bUseEntropy;
    int bUseOldGravMethod;
    int bUseEOSSpeedup;
    int bUseVanLeer;
    double dMachThreshold;
    double dRedSpeedOfLight;
    int iParticleTrack;
    int iAngleTrack;
    int bFirstOrder;
    int bDoVoroStats;
    double dgrav;
    double dgravL;
    int iUserGravityType;
    double dKlo;
    double dKhi;
    double dKslope;
    double dTurbAmplitude;
    int iTurbSeed;        
    char rGravFile[256];
    double dMbh;
    double dxbh;
    double dybh;
    double dzbh;
    double dvxbh;
    double dvybh;
    double dvzbh;
    double drbhSoft;
    double dRadPosx1;
    double dRadPosy1;
    double dRadPosz1;
    double dRadPosx2;
    double dRadPosy2;
    double dRadPosz2;
    double dRadNormx1;
    double dRadNormy1;
    double dRadNormz1;
    double dRadNormx2;
    double dRadNormy2;
    double dRadNormz2;
    double dRadMu0;
    double dRadRadius;
    double dRadTeff;
    double dRadY0;
    double dRadDeltaY;
    double dRadTinit;
    double dRadMuinit;
    int    bRadPointSource;
    int    bRadSlabSource;
    int    bRadNoLorentz;
    int    bRadNoHydro;
    char   opacityName[256];
    double dSimpleKappaPl;
    double dSimpleKappaMean;
    double dSimpleKappaSca;
    Vector3D<double> dBoundaries[6];
    int    iNumBoundaries;
    int    bUseRadiation[6];
    double dRhoBoundary[6];
    Vector3D<double> dSphericalBC;
    double dSphericalBCRadius;
    double dRhoSpherical;
    double dTempSpherical;
    double dTauBoundary;
    char   problemGeneratorName[256];
    double dWaveRho;
    double dWaveCs;
    double dWaveAmplitude;
    double dWaveKxBox;
    double dWaveRhoReal;
    double dWaveRhoImag;
    double dWaveVxReal;
    double dWaveVxImag;
    double dWaveIeReal;
    double dWaveIeImag;
    double dWaveErReal;
    double dWaveErImag;
    double dWaveFrReal;
    double dWaveFrImag;
    int bDoCheckNaN;
    int bDoOutputPres;
    int bDoOutputGrav;
    double dDeltaX;
    double dvx; 
    int bStarForm;
    Stfm *stfm;
    int bFeedback;
    Fdbk *feedback;
    double dThermalCondCoeff;
    double dThermalCondSatCoeff;
    double dThermalCond2Coeff;
    double dThermalCond2SatCoeff;
    double dThermalCondCoeffCode;
    double dThermalCond2CoeffCode;
    double dEvapMinTemp;
    double dEvapCoeff;
    double dEvapCoeffCode;
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
    p|param.dDeltaMin;
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
    p|param.dVoroKick;
    p|param.dVoroMesh;
    p|param.dHydroTheta;
    p|param.iRiemannSolver;
    p|param.dSmallRho;
    p|param.dBx; p|param.dBy; p|param.dBz;
    p|param.dVPsi;
    p|param.bUseEntropy;
    p|param.bUseOldGravMethod;
    p|param.bUseEOSSpeedup;
    p|param.bUseVanLeer;
    p|param.dMachThreshold;
    p|param.dRedSpeedOfLight;
    p|param.iParticleTrack;
    p|param.iAngleTrack;
    p|param.bFirstOrder;
    p|param.bDoVoroStats;
    p|param.dgrav;
    p|param.dgravL;
    p|param.iUserGravityType;
    p|param.dKlo;
    p|param.dKhi;
    p|param.dKslope;
    p|param.dTurbAmplitude;
    p|param.iTurbSeed;   
    p(param.rGravFile, 256);    
    p|param.dMbh;
    p|param.dxbh;
    p|param.dybh;
    p|param.dzbh;
    p|param.dvxbh;
    p|param.dvybh;
    p|param.dvzbh;
    p|param.drbhSoft;
    p|param.dRadPosx1;
    p|param.dRadPosy1;
    p|param.dRadPosz1;
    p|param.dRadPosx2;
    p|param.dRadPosy2;
    p|param.dRadPosz2;
    p|param.dRadNormx1;
    p|param.dRadNormy1;
    p|param.dRadNormz1;
    p|param.dRadNormx2;
    p|param.dRadNormy2;
    p|param.dRadNormz2;
    p|param.dRadMu0;
    p|param.dRadRadius;
    p|param.dRadTeff;
    p|param.dRadY0;
    p|param.dRadDeltaY;
    p|param.dRadTinit;
    p|param.dRadMuinit;
    p|param.bRadPointSource;
    p|param.bRadSlabSource;
    p|param.bRadNoLorentz;
    p|param.bRadNoHydro;
    p(param.opacityName, 256);    
    p|param.dSimpleKappaPl;
    p|param.dSimpleKappaMean;
    p|param.dSimpleKappaSca;
    for( int i = 0; i < 6; i++) {
        p|param.dBoundaries[i];
        p|param.dRhoBoundary[i];
        p|param.bUseRadiation[i];
    }
    p|param.iNumBoundaries;
    p|param.dSphericalBC;
    p|param.dSphericalBCRadius;
    p|param.dRhoSpherical;
    p|param.dTempSpherical;
    p|param.dTauBoundary;
    p(param.problemGeneratorName, 256);
    p|param.dWaveRho;
    p|param.dWaveCs;
    p|param.dWaveAmplitude;
    p|param.dWaveKxBox;
    p|param.dWaveRhoReal;
    p|param.dWaveRhoImag;
    p|param.dWaveVxReal;
    p|param.dWaveVxImag;
    p|param.dWaveIeReal;
    p|param.dWaveIeImag;
    p|param.dWaveErReal;
    p|param.dWaveErImag;
    p|param.dWaveFrReal;
    p|param.dWaveFrImag; 
    p|param.bDoCheckNaN;
    p|param.bDoOutputPres;
    p|param.bDoOutputGrav;
    p|param.dDeltaX;
    p|param.dvx;
    p|param.bStarForm;
    if(p.isUnpacking())
 	param.stfm = new Stfm();
    p|*param.stfm;
    p|param.bFeedback;
    p|param.feedback;
    p|param.dThermalCondCoeff;
    p|param.dThermalCondSatCoeff;
    p|param.dThermalCond2Coeff;
    p|param.dThermalCond2SatCoeff;
    p|param.dThermalCondCoeffCode;
    p|param.dThermalCond2CoeffCode;
    p|param.dEvapMinTemp;
    p|param.dEvapCoeff;
    p|param.dEvapCoeffCode;
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
