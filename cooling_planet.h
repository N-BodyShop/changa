#ifndef COOLING_PLANET_HINCLUDED
#define COOLING_PLANET_HINCLUDED

/*
 * Cooling code for planet formation simulations.
 * Originally written by James Wadsley, McMaster University for
 * GASOLINE.
 */

/* Global consts */

#include "param.h"

#ifdef __cplusplus
extern "C" {
#endif

#include "stiff.h"

/* Constants */
#define CL_B_gm         (6.022e23*(938.7830/931.494))
#define CL_k_Boltzmann  1.38066e-16
#define CL_eV_erg       1.60219e-12
#define CL_eV_per_K     (CL_k_Boltzmann/CL_eV_erg)
/* 
 * Work around for Dec ev6 flawed
 * treatment of sub-normal numbers 
 */
#define CL_MAX_NEG_EXP_ARG  -500.

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

typedef struct clDerivsDataStruct {
  COOL *cl;
  double rho,PdV,E,T,Y_Total,rFactor;
  double     dlnE;
  int        its;  /* Debug */
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

/* Doesn't do anything, needed by Sph.C */
void CoolInitRatesTable( COOL *cl, COOLPARAM CoolParam);

double clThermalEnergy( double Y_Total, double T );

double clTemperature( double Y_Total, double E );

void CoolAddParams( COOLPARAM *CoolParam, PRM );

/* The following are not implemented for cooling planet but are needed
 * by InOutput.h */
#define COOL_ARRAY0_EXT  "HI"
double COOL_ARRAY0(COOL *cl_, COOLPARTICLE *cp,double aa);
#define COOL_ARRAY0( cl_, cp, aa ) ((cp)->Y_Total)
double COOL_SET_ARRAY0(COOL *cl_, COOLPARTICLE *cp,double aa, double bb_val);
#define COOL_SET_ARRAY0( cl_, cp, aa, bb_val ) ((cp)->Y_Total = (bb_val))

#define COOL_ARRAY1_EXT  "HeI"
double COOL_ARRAY1(COOL *cl_, COOLPARTICLE *cp,double aa);
#define COOL_ARRAY1( cl_, cp, aa ) (0)
double COOL_SET_ARRAY1(COOL *cl_, COOLPARTICLE *cp,double aa, double bb_val);
#define COOL_SET_ARRAY1( cl_, cp, aa, bb_val ) (0)

#define COOL_ARRAY2_EXT  "HeII"
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
double COOL_EDOT( COOL *cl_, COOLPARTICLE *cp_, double ECode_, double rhoCode_, double ZMetal_, double *posCode_ );
#define COOL_EDOT( cl_, cp_, ECode_, rhoCode_, ZMetal_, posCode_) (CoolCodeWorkToErgPerGmPerSec( cl_, CoolEdotInstantCode( cl_, cp_, ECode_, rhoCode_, ZMetal_, posCode_ )))

//// Not implemented, but required to keep compiling from crashing
double COOL_COOLING( COOL *cl_, COOLPARTICLE *cp_, double ECode_, double rhoCode_, double ZMetal_, double *posCode_ );
#define COOL_COOLING( cl_, cp_, ECode_, rhoCode_, ZMetal_, posCode_) (CoolCodeWorkToErgPerGmPerSec( cl_, CoolEdotInstantCode( cl_, cp_, ECode_, rhoCode_, ZMetal_, posCode_ )))

void CoolPARTICLEtoPERBARYON(COOL *cl_, PERBARYON *Y, COOLPARTICLE *cp, double HTotal, double HeTotal);

#define CoolPARTICLEtoPERBARYON(cl_, Y, cp) { \
    (Y)->Total = (cp)->Y_Total; }

double CoolCodeEnergyToTemperature( COOL *Cool, COOLPARTICLE *cp, double E, double fMetal );

/* Note: nod to cosmology (z parameter) unavoidable unless we want to access cosmo.[ch] from here */
void CoolSetTime( COOL *Cool, double dTime, double z );

// Saves CenterOfMass to cooling struct cl
void CoolSetStarCM(COOL *cl, double dCenterOfMass[4]);

/* Unit conversion routines */

double CoolCodeTimeToSeconds( COOL *Cool, double dCodeTime );

#define CoolCodeTimeToSeconds( Cool, dCodeTime ) ((Cool)->dSecUnit*(dCodeTime))

double CoolSecondsToCodeTime( COOL *Cool, double dTime );

#define CoolSecondsToCodeTime( Cool, dTime ) ((dTime)/(Cool)->dSecUnit)

double CoolCodeEnergyToErgPerGm( COOL *Cool, double dCodeEnergy );

#define CoolCodeEnergyToErgPerGm( Cool, dCodeEnergy ) ((Cool)->dErgPerGmUnit*(dCodeEnergy))

double CoolErgPerGmToCodeEnergy( COOL *Cool, double dEnergy );

#define CoolErgPerGmToCodeEnergy( Cool, dEnergy ) ((Cool)->diErgPerGmUnit*(dEnergy))

double CoolCodeWorkToErgPerGmPerSec( COOL *Cool, double dCodeWork );

#define CoolCodeWorkToErgPerGmPerSec( Cool, dCodeWork ) ((Cool)->dErgPerGmPerSecUnit*(dCodeWork))

double CoolErgPerGmPerSecToCodeWork( COOL *Cool, double dWork );

#define CoolErgPerGmPerSecToCodeWork( Cool, dWork ) ((dWork)/(Cool)->dErgPerGmPerSecUnit)

double CodeDensityToComovingGmPerCc( COOL *Cool, double dCodeDensity );

#define CodeDensityToComovingGmPerCc( Cool, dCodeDensity )  ((Cool)->dComovingGmPerCcUnit*(dCodeDensity))

void CoolIntegrateEnergyCode(COOL *cl, clDerivsData *cData, COOLPARTICLE *cp, double *E, 
			     double PdV, double rho, double ZMetal, double *r, double tStep );

void CoolDefaultParticleData( COOLPARTICLE *cp );

void CoolInitEnergyAndParticleData( COOL *cl, COOLPARTICLE *cp, double *E, double dDensity, double dTemp, double fMetal );

/* Not implemented, but required to keep compiling from crashing */
double CoolEdotInstantCode(COOL *cl, COOLPARTICLE *cp, double ECode,
               double rhoCode, double ZMetal, double *posCode );

void CoolCodePressureOnDensitySoundSpeed( COOL *cl, COOLPARTICLE *cp, double uPred, double fDensity, double gamma, double gammam1, double *PoverRho, double *c );
#define CoolCodePressureOnDensitySoundSpeed( cl__, cp__, uPred__, fDensity__, gamma__, gammam1__, PoverRho__, c__ ) { \
  *(PoverRho__) = ((gammam1__)*(uPred__)); \
  *(c__) = sqrt((gamma__)*(*(PoverRho__))); }

/* Place holder function (not implemented but needed */
void CoolTableReadInfo( COOLPARAM *CoolParam, int cntTable, int *nTableColumns, char *suffix );

/* Place holder function (not implemented but needed */
void CoolTableRead( COOL *Cool, int nData, void *vData);

#ifdef __cplusplus
}
#endif

#endif


