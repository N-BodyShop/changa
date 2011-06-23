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
	double    dParam1;
	double    dParam2;
	double    dParam3;
	double    dParam4;
	double    Y_Total;
	double dCoolingTmin;
	double dCoolingTmax;
	} COOLPARAM;

/* How to do a dummy? --> can't I guess */
typedef struct CoolingParticleStruct {
	double Y_Total;
	} COOLPARTICLE;

typedef struct CoolingPKDStruct COOL;
    
typedef struct clDerivsDataStruct {
  STIFF *IntegratorContext;
  COOL *cl;
  double rho,PdV,E,T,Y_Total,rFactor;
  double     dlnE;
  int        its;  /* Debug */
} clDerivsData;

/* Heating Cooling Context */

struct CoolingPKDStruct { 
   double     z; /* Redshift */
   double     dTime;

   double     dGmPerCcUnit;
   double     dComovingGmPerCcUnit;
   double     dErgPerGmUnit;
   double     dSecUnit;
   double     dErgPerGmPerSecUnit;
   double     diErgPerGmUnit;
   double     dKpcUnit;

   double    dParam1;
   double    dParam2;
   double    dParam3;
   double    dParam4;
   double    Y_Total;
   double    Tmin;
   double    Tmax;
   clDerivsData       *DerivsData;

   int nTableRead; /* Internal Tables read from Files */

};


COOL *CoolInit( );
void CoolFinalize( COOL *cl );
clDerivsData *CoolDerivsInit(COOL *cl);
void CoolDerivsFinalize(clDerivsData *cld ) ;

void clInitConstants( COOL *cl, double dGMPerCcunit,
		      double dComovingGmPerCcUnit, double dErgPerGmUnit,
		      double dSecUnit, double dKpcUnit, COOLPARAM CoolParam);
void clInitUV(COOL *cl, int nTableColumns, int nTableRows, double *dTableData );
void clInitRatesTable( COOL *cl, double TMin, double TMax, int nTable );
void CoolInitRatesTable( COOL *cl, COOLPARAM CoolParam);

double clThermalEnergy( double Y_Total, double T );
double clTemperature( double Y_Total, double E );

double clEdotInstant( COOL *cl, double E, double T, double rho, double r );
void clIntegrateEnergy(COOL *cl, clDerivsData *clData,  double *E, double PdV,
		       double rho, double Y_Total, double radius, double tStep );

int clDerivs(double x, const double *y, double *dydx, void *Data) ;

int clJacobn(double x, const double y[], double dfdx[], double *dfdy, void *Data) ;
  
void CoolAddParams( COOLPARAM *CoolParam, PRM );
void CoolLogParams( COOLPARAM *CoolParam, FILE *fp );
void CoolOutputArray( COOLPARAM *CoolParam, int, int *, char * );

#define COOL_ARRAY0_EXT  "HI"
double COOL_ARRAY0(COOL *cl_, COOLPARTICLE *cp,double aa);
#define COOL_ARRAY0( cl_, cp, aa ) ((cp)->Y_Total)

#define COOL_ARRAY1_EXT  "HeI"
double COOL_ARRAY1(COOL *cl_, COOLPARTICLE *cp,double aa);
#define COOL_ARRAY1( cl_, cp, aa ) (0)

#define COOL_ARRAY2_EXT  "HeII"
double COOL_ARRAY2(COOL *cl_, COOLPARTICLE *cp,double aa);
#define COOL_ARRAY2( cl_, cp, aa ) (0)

#define COOL_ARRAY3_EXT  "H2"
double COOL_ARRAY3(COOL *cl, COOLPARTICLE *cp, double ZMetal);
#define COOL_ARRAY3(cl_, cp, aa ) (0)

double COOL_EDOT( COOL *cl_, COOLPARTICLE *cp_, double ECode_, double rhoCode_, double ZMetal_, double *posCode_ );
#define COOL_EDOT( cl_, cp_, ECode_, rhoCode_, ZMetal_, posCode_) (CoolCodeWorkToErgPerGmPerSec( cl_, CoolEdotInstantCode( cl_, cp_, ECode_, rhoCode_, ZMetal_, posCode_ )))

double COOL_COOLING( COOL *cl_, COOLPARTICLE *cp_, double ECode_, double rhoCode_, double ZMetal_, double *posCode_ );
#define COOL_COOLING( cl_, cp_, ECode_, rhoCode_, ZMetal_, posCode_) (CoolCodeWorkToErgPerGmPerSec( cl_, CoolEdotInstantCode( cl_, cp_, ECode_, rhoCode_, ZMetal_, posCode_ )))

double COOL_HEATING( COOL *cl_, COOLPARTICLE *cp_, double ECode_, double rhoCode_, double ZMetal_, double *posCode_ );
#define COOL_HEATING( cl_, cp_, ECode_, rhoCode_, ZMetal_, posCode_) (0)

double CoolCodeEnergyToTemperature( COOL *Cool, COOLPARTICLE *cp, double E, double fMetal );

/* Note: nod to cosmology (z parameter) unavoidable unless we want to access cosmo.[ch] from here */
void CoolSetTime( COOL *Cool, double dTime, double z );

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

void CoolIntegrateEnergy(COOL *cl, clDerivsData *cData, COOLPARTICLE *cp, double *E, 
		       double PdV, double rho, double ZMetal, double tStep );

void CoolIntegrateEnergyCode(COOL *cl, clDerivsData *cData, COOLPARTICLE *cp, double *E, 
			     double PdV, double rho, double ZMetal, double *r, double tStep );

void CoolDefaultParticleData( COOLPARTICLE *cp );

void CoolInitEnergyAndParticleData( COOL *cl, COOLPARTICLE *cp, double *E, double dDensity, double dTemp, double fMetal );

/* Deprecated */
double CoolHeatingRate( COOL *cl, COOLPARTICLE *cp, double E, double dDensity );

double CoolEdotInstantCode(COOL *cl, COOLPARTICLE *cp, double ECode, 
			   double rhoCode, double ZMetal, double *posCode );

void CoolCodePressureOnDensitySoundSpeed( COOL *cl, COOLPARTICLE *cp, double uPred, double fDensity, double gamma, double gammam1, double *PoverRho, double *c );
#define CoolCodePressureOnDensitySoundSpeed( cl__, cp__, uPred__, fDensity__, gamma__, gammam1__, PoverRho__, c__ ) { \
  *(PoverRho__) = ((gammam1__)*(uPred__)); \
  *(c__) = sqrt((gamma__)*(*(PoverRho__))); }

void CoolTableReadInfo( COOLPARAM *CoolParam, int cntTable, int *nTableColumns, char *suffix );

void CoolTableRead( COOL *Cool, int nData, void *vData);

#ifdef __cplusplus
}
#endif

#endif


