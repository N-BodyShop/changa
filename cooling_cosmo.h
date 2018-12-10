#ifndef COOLING_COSMO_HINCLUDED
#define COOLING_COSMO_HINCLUDED

/*
 * Cooling code for cosmology simulations.
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
#define CL_RT_FLOAT      float
#define CL_RT_MIN        1e-38
#define CL_RT_MIN        FLT_MIN
*/

#define CL_RT_FLOAT      double
#define CL_RT_MIN        1e-100

/*
#define CL_RT_MIN        DBL_MIN
*/
/* 
 * Work around for Dec ev6 flawed
 * treatment of sub-normal numbers 
 */
#define CL_MAX_NEG_EXP_ARG  -500.

#define CL_NMAXBYTETABLE   56000

/** @brief Input parameters for cooling */
typedef struct CoolingParametersStruct {
	int    bIonNonEqm;
	int    nCoolingTable;
	int    bUV;
	int    bUVTableUsesTime;
	int    bDoIonOutput;
        int    bLowTCool;
        int    bSelfShield;
	double dMassFracHelium;
	double dCoolingTmin;
	double dCoolingTmax;
	} COOLPARAM;

/** @brief per-particle cooling data */
typedef struct CoolingParticleStruct {
	double Y_HI,Y_HeI,Y_HeII;	/* Abundance of ions */
        double     dLymanWerner; /* Flux of Lyman Werner radiation at the gas particle.  Temporary CC */
	} COOLPARTICLE;

/** @brief abundance of various species in particles/baryon */
typedef struct { 
  double e,Total;
  double HI,HII,HeI,HeII,HeIII;
} PERBARYON;

/** @brief photoionization and heating rates from a uniform UV background */
typedef struct { 
  double   zTime;

  double   Rate_Phot_HI;
  double   Rate_Phot_HeI;
  double   Rate_Phot_HeII;

  double   Heat_Phot_HI;
  double   Heat_Phot_HeI;
  double   Heat_Phot_HeII;
} UVSPECTRUM;

/** @brief structure to hold Temperature independent cooling and heating
    rates */
typedef struct { 
  double   Rate_Phot_HI;
  double   Rate_Phot_HeI;
  double   Rate_Phot_HeII;

  double   Heat_Phot_HI;
  double   Heat_Phot_HeI;
  double   Heat_Phot_HeII;
 
  double   Cool_Coll_HI;
  double   Cool_Coll_HeI;
  double   Cool_Coll_HeII;
  double   Cool_Diel_HeII;
  double   Cool_Comp;
  double   Tcmb;
  double   Cool_LowTFactor;

} RATES_NO_T;

/** @brief structure to hold Temperature dependent cooling rates */
typedef struct { 
  CL_RT_FLOAT   Rate_Coll_HI;
  CL_RT_FLOAT   Rate_Coll_HeI;
  CL_RT_FLOAT   Rate_Coll_HeII;
  CL_RT_FLOAT   Rate_Radr_HII;
  CL_RT_FLOAT   Rate_Radr_HeII;
  CL_RT_FLOAT   Rate_Radr_HeIII;
  CL_RT_FLOAT   Rate_Diel_HeII;

  CL_RT_FLOAT   Cool_Brem_1;
  CL_RT_FLOAT   Cool_Brem_2;
  CL_RT_FLOAT   Cool_Radr_HII;
  CL_RT_FLOAT   Cool_Radr_HeII;
  CL_RT_FLOAT   Cool_Radr_HeIII;
  CL_RT_FLOAT   Cool_Line_HI;
  CL_RT_FLOAT   Cool_Line_HeI;
  CL_RT_FLOAT   Cool_Line_HeII;
  CL_RT_FLOAT   Cool_LowT;

} RATES_T;

typedef struct clDerivsDataStruct clDerivsData;
    
/**@brief  Heating/Cooling context: parameters and tables */

typedef struct CoolingPKDStruct { 
   double     z; /* Redshift */
   double     dTime;
 /* Rates independent of Temperature */ 
   RATES_NO_T  R;
 /* Table for Temperature dependent rates */ 
   int        nTable;
   double     TMin;
   double     TMax;
   double     TlnMin;
   double     TlnMax;
   double     rDeltaTln;
   RATES_T     *RT;

   int        nTableRead; /* number of Tables read from files */

   int        bUV;
   int        nUV;
   UVSPECTRUM *UV;
   int        bUVTableUsesTime;
   int        bUVTableLinear;
   int        bLowTCool;
   int        bSelfShield;

   double     dGmPerCcUnit;
   double     dComovingGmPerCcUnit;
   double     dErgPerGmUnit;
   double     dSecUnit;
   double     dErgPerGmPerSecUnit;
   double     diErgPerGmUnit;
   double     dKpcUnit;
   double     Y_H;
   double     Y_He;
   double     Y_eMAX;
/* Diagnostic */
   int        its;
} COOL;

/** @brief Rate information for a given particle */
typedef struct {
  double   T, Tln;
  double   Coll_HI;
  double   Coll_HeI;
  double   Coll_HeII;
  double   Radr_HII;
  double   Radr_HeII;
  double   Diel_HeII;
  double   Totr_HeII;
  double   Radr_HeIII;

  double   Phot_HI;
  double   Phot_HeI;
  double   Phot_HeII;
} RATE;

/** @brief return structure for clTestCool() */
typedef struct {
  double compton;
  double bremHII;
  double bremHeII;
  double bremHeIII;
  double radrecHII;
  double radrecHeII;
  double radrecHeIII;
  double collionHI; 
  double collionHeI;
  double collionHeII;
  double dielrecHeII;
  double lineHI;
  double lineHeI;
  double lineHeII;
  double lowT;
} COOL_ERGPERSPERGM;


/** @brief context for calculating cooling derivatives */
struct clDerivsDataStruct {
  STIFF *IntegratorContext;	/**< @brief Context for diff. eq. integrator */
  COOL *cl;			/**< @brief pointer to cooling context  */
  double rho,ExternalHeating,E,ZMetal;
  RATE Rate;
  PERBARYON Y;
  double     Y_Total0, Y_Total1;
  int        its;  /* Debug */
};


COOL *CoolInit( );
void CoolFinalize( COOL *cl );
clDerivsData *CoolDerivsInit(COOL *cl);
void CoolDerivsFinalize(clDerivsData *cld ) ;

void clInitConstants( COOL *cl, double dGMPerCcunit, double dComovingGmPerCcUnit,
					 double dErgPerGmUnit, double dSecUnit, double dKpcUnit, COOLPARAM CoolParam);
void clInitUV(COOL *cl, int nTableColumns, int nTableRows, double *dTableData );
void clInitRatesTable( COOL *cl, double TMin, double TMax, int nTable );
void CoolInitRatesTable( COOL *cl, COOLPARAM CoolParam);

void clRatesTableError( COOL *cl );
void clRatesRedshift( COOL *cl, double z, double dTime );
double clHeatTotal ( COOL *cl, PERBARYON *Y, RATE *Rate );
void clRates( COOL *cl, RATE *Rate, double T, double rho );
double clCoolTotal( COOL *cl, PERBARYON *Y, RATE *Rate, double rho, double ZMetal );
COOL_ERGPERSPERGM  clTestCool ( COOL *cl, PERBARYON *Y, RATE *Rate, double rho );
void clPrintCool( COOL *cl, PERBARYON *Y, RATE *Rate, double rho );
void clPrintCoolFile( COOL *cl, PERBARYON *Y, RATE *Rate, double rho, FILE *fp );

void clAbunds( COOL *cl, PERBARYON *Y, RATE *Rate, double rho);
double clThermalEnergy( double Y_Total, double T );
double clTemperature( double Y_Total, double E );
double clRateCollHI( double T );
double clRateCollHeI( double T );
double clRateCollHeII( double T );
double clRateRadrHII( double T );
double clRateRadrHeII( double T );
double clRateDielHeII( double T );
double clRateRadrHeIII( double T );
double clCoolBrem1( double T );
double clCoolBrem2( double T );
double clCoolRadrHII( double T );
double clCoolRadrHeII( double T );
double clCoolRadrHeIII( double T );
double clCoolLineHI( double T );
double clCoolLineHeI( double T );
double clCoolLineHeII( double T );
double clCoolLowT( double T );

double clEdotInstant ( COOL *cl, PERBARYON *Y, RATE *Rate, double rho,
		       double ZMetal, double *EdotHeat, double *EdotCool );
void clIntegrateEnergy(COOL *cl, clDerivsData *clData, PERBARYON *Y, double *E, 
		       double ExternalHeating, double rho, double ZMetal, double dt );
void clIntegrateEnergyDEBUG(COOL *cl, PERBARYON *Y, double *E, 
		       double ExternalHeating, double rho, double ZMetal, double dt );


void clDerivs(double x, const double *y, double *dHeat, double *dCool,
	     void *Data) ;

int clJacobn( double x, const double y[], double dfdx[], double *dfdy, void *Data) ;
  
void CoolAddParams( COOLPARAM *CoolParam, PRM );
void CoolLogParams( COOLPARAM *CoolParam, FILE *fp );
void CoolOutputArray( COOLPARAM *CoolParam, int, int *, char * );

#define COOL_ARRAY0_EXT  "HI"
double COOL_ARRAY0(COOL *cl_, COOLPARTICLE *cp,double aa);
#define COOL_ARRAY0( cl_, cp, aa ) ((cp)->Y_HI)
double COOL_SET_ARRAY0(COOL *cl_, COOLPARTICLE *cp,double aa, double bb_val);
#define COOL_SET_ARRAY0( cl_, cp, aa, bb_val ) ((cp)->Y_HI = (bb_val))

#define COOL_ARRAY1_EXT  "HeI"
double COOL_ARRAY1(COOL *cl_, COOLPARTICLE *cp,double aa);
#define COOL_ARRAY1( cl_, cp, aa ) ((cp)->Y_HeI)
double COOL_SET_ARRAY1(COOL *cl_, COOLPARTICLE *cp,double aa, double bb_val);
#define COOL_SET_ARRAY1( cl_, cp, aa, bb_val ) ((cp)->Y_HeI = (bb_val))

#define COOL_ARRAY2_EXT  "HeII"
double COOL_ARRAY2(COOL *cl_, COOLPARTICLE *cp,double aa);
#define COOL_ARRAY2( cl_, cp, aa ) ((cp)->Y_HeII)
double COOL_SET_ARRAY2(COOL *cl_, COOLPARTICLE *cp,double aa, double bb_val);
#define COOL_SET_ARRAY2( cl_, cp, aa, bb_val ) ((cp)->Y_HeII = (bb_val))

#define COOL_ARRAY3_EXT  "H2"
double COOL_ARRAY3(COOL *cl, COOLPARTICLE *cp, double ZMetal);
#define COOL_ARRAY3(cl_, cp, aa ) (0)
double COOL_SET_ARRAY3(COOL *cl_, COOLPARTICLE *cp,double aa, double bb_val);
#define COOL_SET_ARRAY3( cl_, cp, aa, bb_val ) (0)

double COOL_EDOT( COOL *cl_, COOLPARTICLE *cp_, double ECode_, double rhoCode_, double ZMetal_, double *posCode_ );
#define COOL_EDOT( cl_, cp_, ECode_, rhoCode_, ZMetal_, posCode_) (CoolCodeWorkToErgPerGmPerSec( cl_, CoolEdotInstantCode( cl_, cp_, ECode_, rhoCode_, ZMetal_, posCode_ )))

double COOL_COOLING( COOL *cl_, COOLPARTICLE *cp_, double ECode_, double rhoCode_, double ZMetal_, double *posCode_ );
#define COOL_COOLING( cl_, cp_, ECode_, rhoCode_, ZMetal_, posCode_) (CoolCodeWorkToErgPerGmPerSec( cl_, CoolCoolingCode( cl_, cp_, ECode_, rhoCode_, ZMetal_, posCode_ )))

double COOL_HEATING( COOL *cl_, COOLPARTICLE *cp_, double ECode_, double rhoCode_, double ZMetal_, double *posCode_ );
#define COOL_HEATING( cl_, cp_, ECode_, rhoCode_, ZMetal_, posCode_) (CoolCodeWorkToErgPerGmPerSec( cl_, CoolHeatingCode( cl_, cp_, ECode_, rhoCode_, ZMetal_, posCode_ )))

void CoolPARTICLEtoPERBARYON(COOL *cl_, PERBARYON *Y, COOLPARTICLE *cp, double HTotal, double HeTotal);

#define CoolPARTICLEtoPERBARYON(cl_, Y, cp) { \
    (Y)->HI = (cp)->Y_HI; \
    (Y)->HII = (cl_)->Y_H - (Y)->HI; \
    (Y)->HeI = (cp)->Y_HeI; \
    (Y)->HeII = (cp)->Y_HeII; \
    (Y)->HeIII = (cl_)->Y_He - (Y)->HeI - (Y)->HeII; \
    (Y)->e = (Y)->HII + (Y)->HeII + 2*(Y)->HeIII; \
    (Y)->Total = (Y)->e + (cl_)->Y_H + (cl_)->Y_He; }

void CoolPERBARYONtoPARTICLE(COOL *cl_, PERBARYON *Y, COOLPARTICLE *cp);

#define CoolPERBARYONtoPARTICLE(cl_, Y, cp) { \
    (cp)->Y_HI = (Y)->HI; \
    (cp)->Y_HeI = (Y)->HeI; \
    (cp)->Y_HeII = (Y)->HeII; }


double CoolEnergyToTemperature( COOL *Cool, COOLPARTICLE *cp, double E, double fMetal );
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
		       double ExternalHeating, double rho, double ZMetal, double tStep );

void CoolIntegrateEnergyCode(COOL *cl, clDerivsData *cData, COOLPARTICLE *cp, double *E, 
		       double ExternalHeating, double rho, double ZMetal, double *r, double tStep );

void CoolDefaultParticleData( COOLPARTICLE *cp );

void CoolInitEnergyAndParticleData( COOL *cl, COOLPARTICLE *cp, double *E, double dDensity, double dTemp, double fMetal );

/* Deprecated */
double CoolHeatingRate( COOL *cl, COOLPARTICLE *cp, double E, double dDensity, double ZMetal );

double CoolEdotInstantCode(COOL *cl, COOLPARTICLE *cp, double ECode, 
			   double rhoCode, double ZMetal, double *posCode );
double CoolCoolingCode(COOL *cl, COOLPARTICLE *cp, double ECode, 
			   double rhoCode, double ZMetal, double *posCode );
double CoolHeatingCode(COOL *cl, COOLPARTICLE *cp, double ECode, 
			   double rhoCode, double ZMetal, double *posCode );

void CoolCodePressureOnDensitySoundSpeed( COOL *cl, COOLPARTICLE *cp, double uPred, double fDensity, double gamma, double gammam1, double *PoverRho, double *c );

/* Note: gamma should be 5/3 for this to be consistent! */
#define CoolCodePressureOnDensitySoundSpeed( cl__, cp__, uPred__, fDensity__, gamma__, gammam1__, PoverRho__, c__ ) { \
  *(PoverRho__) = ((5./3.-1)*(uPred__)); \
  *(c__) = sqrt((5./3.)*(*(PoverRho__))); }

/*
double CoolCodePressureOnDensity( COOL *cl, COOLPARTICLE *cp, double uPred, double fDensity, double gammam1 );

#define CoolCodePressureOnDensity( cl, cp, uPred, fDensity, gammam1 ) ((gammam1)*(uPred))
*/

void CoolTableReadInfo( COOLPARAM *CoolParam, int cntTable, int *nTableColumns, char *suffix );

void CoolTableRead( COOL *Cool, int nData, void *vData);

#ifdef __cplusplus
}
#endif

#endif


