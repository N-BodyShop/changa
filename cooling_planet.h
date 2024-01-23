#ifndef COOLING_PLANET_HINCLUDED
#define COOLING_PLANET_HINCLUDED
/*
 * Cooling code for planet formation simulations.
 * Originally written by James Wadsley, McMaster University for
 * GASOLINE.
 *
 * Updated for ChaNGa by Isaac Backus, University of Wasghinton
 */

#include "param.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Constants */
#define CL_Rgascode         8.2494e7
#define CL_Eerg_gm_degK     CL_Rgascode
#define CL_Eerg_gm_degK3_2  1.5*CL_Eerg_gm_degK
#define CL_NMAXBYTETABLE   56000

typedef struct CoolingParametersStruct {
    double    Y_Total;
    double dCoolingTmin;
    double dCoolingTmax;
    double dBetaCooling;
} COOLPARAM;

typedef struct CoolingParticleStruct {
    double Y_Total;
} COOLPARTICLE;

typedef struct {
    double Total;
} PERBARYON;

typedef struct CoolingPKDStruct COOL;

/* Per-thread data.  Not used for this cooling algorithm */
typedef struct clDerivsDataStruct {
} clDerivsData;

/* Heating Cooling Context */

struct CoolingPKDStruct {
    /* Cosmology hold-overs */
    double     z; /* Redshift */
    double     dTime;
    /* Units and conversion constants */
    double     dGmPerCcUnit;
    double     dComovingGmPerCcUnit;
    double     dErgPerGmUnit;
    double     dSecUnit;
    double     dErgPerGmPerSecUnit;
    double     diErgPerGmUnit;
    double     dKpcUnit;
    /* User parameters (see CoolAddParams) */
    double    Y_Total;
    double    Tmin;
    double    Tmax;
    double    beta;
    /* Star info */
    double dStarCenterOfMass[4];
    /**/
    clDerivsData       *DerivsData;

    int nTableRead; /* Internal Tables read from Files */

};

COOL *CoolInit( );

///Frees memory and deletes cl
void CoolFinalize( COOL *cl );

clDerivsData *CoolDerivsInit(COOL *cl);

void CoolDerivsFinalize(clDerivsData *cld ) ;

void clInitConstants( COOL *cl, double dGMPerCcunit,
                      double dComovingGmPerCcUnit, double dErgPerGmUnit,
                      double dSecUnit, double dKpcUnit, COOLPARAM CoolParam);

/* Doesn't do anything, needed by Sph.cpp */
void CoolInitRatesTable( COOL *cl, COOLPARAM CoolParam);

double clThermalEnergy( double Y_Total, double T );

double clTemperature( double Y_Total, double E );

void CoolAddParams( COOLPARAM *CoolParam, PRM );

// Saves CenterOfMass to cooling struct cl
void CoolSetStarCM(COOL *cl, double dCenterOfMass[4]);

/* Needed by InOutput.h */
#define COOL_ARRAY0_EXT  "Y_Total"
double COOL_ARRAY0(COOL *cl_, COOLPARTICLE *cp,double aa);
#define COOL_ARRAY0( cl_, cp, aa ) ((cp)->Y_Total)
double COOL_SET_ARRAY0(COOL *cl_, COOLPARTICLE *cp,double aa, double bb_val);
#define COOL_SET_ARRAY0( cl_, cp, aa, bb_val ) ((cp)->Y_Total = (bb_val))

#define COOL_ARRAY1_EXT  "NA"
double COOL_ARRAY1(COOL *cl_, COOLPARTICLE *cp,double aa);
#define COOL_ARRAY1( cl_, cp, aa ) (0)
double COOL_SET_ARRAY1(COOL *cl_, COOLPARTICLE *cp,double aa, double bb_val);
#define COOL_SET_ARRAY1( cl_, cp, aa, bb_val ) (0)

#define COOL_ARRAY2_EXT  "NA"
double COOL_ARRAY2(COOL *cl_, COOLPARTICLE *cp,double aa);
#define COOL_ARRAY2( cl_, cp, aa ) (0)
double COOL_SET_ARRAY2(COOL *cl_, COOLPARTICLE *cp,double aa, double bb_val);
#define COOL_SET_ARRAY2( cl_, cp, aa, bb_val ) (0)

#define COOL_ARRAY3_EXT  "H2"
double COOL_ARRAY3(COOL *cl, COOLPARTICLE *cp, double ZMetal);
#define COOL_ARRAY3(cl_, cp, aa ) (0)
double COOL_SET_ARRAY3(COOL *cl_, COOLPARTICLE *cp,double aa, double bb_val);
#define COOL_SET_ARRAY3( cl_, cp, aa, bb_val ) (0)

/// Not implemented, but required to keep compiling from crashing
double COOL_EDOT( COOL *cl_, COOLPARTICLE *cp_, double ECode_, double rhoCode_,
                  double ZMetal_, double *posCode_ );
#define COOL_EDOT( cl_, cp_, ECode_, rhoCode_, ZMetal_, posCode_) (0)

//// Not implemented, but required to keep compiling from crashing
double COOL_COOLING( COOL *cl_, COOLPARTICLE *cp_, double ECode_,
                     double rhoCode_, double ZMetal_, double *posCode_ );
#define COOL_COOLING( cl_, cp_, ECode_, rhoCode_, ZMetal_, posCode_) (0)

void CoolPARTICLEtoPERBARYON(COOL *cl_, PERBARYON *Y, COOLPARTICLE *cp,
                             double HTotal, double HeTotal);

#define CoolPARTICLEtoPERBARYON(cl_, Y, cp) { \
    (Y)->Total = (cp)->Y_Total; }

double CoolCodeEnergyToTemperature( COOL *Cool, COOLPARTICLE *cp, double E,
                                    double fMetal );

/* Note: nod to cosmology (z parameter) unavoidable unless we want to access
 * cosmo.[ch] from here */
void CoolSetTime( COOL *Cool, double dTime, double z );

/* Unit conversion routines */

double CoolCodeTimeToSeconds( COOL *Cool, double dCodeTime );

#define CoolCodeTimeToSeconds( Cool, dCodeTime ) ((Cool)->dSecUnit*(dCodeTime))

double CoolSecondsToCodeTime( COOL *Cool, double dTime );

#define CoolSecondsToCodeTime( Cool, dTime ) ((dTime)/(Cool)->dSecUnit)

double CoolCodeEnergyToErgPerGm( COOL *Cool, double dCodeEnergy );

#define CoolCodeEnergyToErgPerGm( Cool, dCodeEnergy ) \
    ((Cool)->dErgPerGmUnit*(dCodeEnergy))

double CoolErgPerGmToCodeEnergy( COOL *Cool, double dEnergy );

#define CoolErgPerGmToCodeEnergy( Cool, dEnergy ) \
    ((Cool)->diErgPerGmUnit*(dEnergy))

double CoolCodeWorkToErgPerGmPerSec( COOL *Cool, double dCodeWork );

#define CoolCodeWorkToErgPerGmPerSec( Cool, dCodeWork ) \
    ((Cool)->dErgPerGmPerSecUnit*(dCodeWork))

double CoolErgPerGmPerSecToCodeWork( COOL *Cool, double dWork );

#define CoolErgPerGmPerSecToCodeWork( Cool, dWork ) \
    ((dWork)/(Cool)->dErgPerGmPerSecUnit)

double CodeDensityToComovingGmPerCc( COOL *Cool, double dCodeDensity );

#define CodeDensityToComovingGmPerCc( Cool, dCodeDensity )  \
    ((Cool)->dComovingGmPerCcUnit*(dCodeDensity))

void CoolIntegrateEnergyCode(COOL *cl, clDerivsData *cData, COOLPARTICLE *cp,
                             double *E, double PdV, double rho, double ZMetal,
                             double *r, double tStep );

void CoolDefaultParticleData( COOLPARTICLE *cp );

void CoolInitEnergyAndParticleData( COOL *cl, COOLPARTICLE *cp, double *E,
    double dDensity, double dTemp, double fMetal);

/* Not implemented, but required to keep compiling from crashing */
double CoolEdotInstantCode(COOL *cl, COOLPARTICLE *cp, double ECode,
                           double rhoCode, double ZMetal, double *posCode );

void CoolCodePressureOnDensitySoundSpeed( COOL *cl, COOLPARTICLE *cp,
    double uPred, double fDensity, double gamma, double gammam1,
    double *PoverRho, double *c );
#define CoolCodePressureOnDensitySoundSpeed( cl__, cp__, uPred__, fDensity__, \
    gamma__, gammam1__, PoverRho__, c__ ) { \
  *(PoverRho__) = ((gammam1__)*(uPred__)); \
  *(c__) = sqrt((gamma__)*(*(PoverRho__))); }

/* Place holder function (not implemented but needed */
void CoolTableReadInfo( COOLPARAM *CoolParam, int cntTable, int *nTableColumns,
                        char *suffix );

/* Place holder function (not implemented but needed */
void CoolTableRead( COOL *Cool, int nData, void *vData);

#ifdef __cplusplus
}
#endif

#endif


