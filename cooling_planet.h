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

/** @brief structure containing the user set cooling parameters */
typedef struct CoolingParametersStruct {
    double    dCoolRateEff;
	double    Y_Total;
	double dCoolingTmin;
	double dCoolingTmax;
	} COOLPARAM;

/** @brief Defines the cooling particle structure.  The only
 * piece of information each particle needs is Y_total which
 * is equal the (# deg. freedom)/3
 */
typedef struct CoolingParticleStruct {
	double Y_Total;
	} COOLPARTICLE;

/** @brief abundance of various species in particles/baryon */
typedef struct { 
  double Total;
} PERBARYON;

typedef struct CoolingPKDStruct COOL;

/** @brief Cooling data for each particle */
typedef struct clDerivsDataStruct {
  COOL *cl;
  double rho,PdV,E,T,Y_Total,rFactor;
  double     dlnE;
  int        its;  /* Debug */
} clDerivsData;

/* Heating Cooling Context */

/** @brief Cooling data accessible by all threads */
struct CoolingPKDStruct { 
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
   double    dCoolRateEff;
   double    Y_Total;
   double    Tmin;
   double    Tmax;
   clDerivsData       *DerivsData;

   int nTableRead; /* Internal Tables read from Files */

};

/**
 * @brief CoolInit initializes a COOL data structure for storing
 *      cooling constants and parameters
 * @return a blank COOL data structure with nTableRead = 0
 */
COOL *CoolInit( );

///Frees memory and deletes cl
void CoolFinalize( COOL *cl );

/**
 * Per thread initialization of cooling data.
 * @param cl Initialized COOL structure.
 */
clDerivsData *CoolDerivsInit(COOL *cl);

/**
 * Deallocate memory for per-thread data.
 */
void CoolDerivsFinalize(clDerivsData *cld ) ;

/**
 * @brief clInitConstants sets physical constants and parameters
 *          for cooling
 * @param cl a COOL data structure (see CoolInit)
 * @param CoolParam an initialized COOLPARAM structure (see CoolAddParams)
 */
void clInitConstants( COOL *cl, double dGMPerCcunit,
		      double dComovingGmPerCcUnit, double dErgPerGmUnit,
		      double dSecUnit, double dKpcUnit, COOLPARAM CoolParam);

//void clInitUV(COOL *cl, int nTableColumns, int nTableRows, double *dTableData );
//void clInitRatesTable( COOL *cl, double TMin, double TMax, int nTable );

/* Doesn't do anything, needed by Sph.C */
void CoolInitRatesTable( COOL *cl, COOLPARAM CoolParam);

/**
 * @brief Calculates thermal energy from temperature and Y (deg. freedom/3)
 * @param Y_Total (deg. freedom)/3
 * @param T Temperature
 * @return Thermal (internal) energy per mass (CGS)
 */
double clThermalEnergy( double Y_Total, double T );

/**
 * @brief Calculates temperature from internal energy and Y (deg.freedom/3)
 * @param Y_Total (deg. freedom)/3
 * @param E Internal energy per mass (CGS units)
 * @return Temperature
 */
double clTemperature( double Y_Total, double E );

/**
 * @brief CoolAddParams Parses parameters and adds them to the COOLPARAM
 *          struct
 * @param CoolParam A COOLPARAM struct
 * @param prm The param struct
 *
 * The available parameters to be set are:
 *      dCoolRateEff : The "effective" cooling rate.  This is given by
 *          (r0^3/2)/t_cool where t_cool is the desired cooling rate
 *          at r0, in cgs units.
 *      dY_Total : (degrees of freedom)/3 for the gas
 *      dCoolingTmin : Minimum allowed temperature (in Kelvin)
 *      dCoolingTmax : Maximum allowed temperature (in Kelvin)
 */
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

//double COOL_HEATING( COOL *cl_, COOLPARTICLE *cp_, double ECode_, double rhoCode_, double ZMetal_, double *posCode_ );
//#define COOL_HEATING( cl_, cp_, ECode_, rhoCode_, ZMetal_, posCode_) (0)

void CoolPARTICLEtoPERBARYON(COOL *cl_, PERBARYON *Y, COOLPARTICLE *cp, double HTotal, double HeTotal);

#define CoolPARTICLEtoPERBARYON(cl_, Y, cp) { \
    (Y)->Total = (cp)->Y_Total; }

void CoolPARTICLEtoPERBARYON(COOL *cl_, PERBARYON *Y, COOLPARTICLE *cp, double HTotal, double HeTotal);

#define CoolPARTICLEtoPERBARYON(cl_, Y, cp) { \
    (Y)->Total = (cp)->Y_Total; }

double CoolCodeEnergyToTemperature( COOL *Cool, COOLPARTICLE *cp, double E, double fMetal );

/* Note: nod to cosmology (z parameter) unavoidable unless we want to access cosmo.[ch] from here */
void CoolSetTime( COOL *Cool, double dTime, double z );

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

//void CoolIntegrateEnergy(COOL *cl, clDerivsData *cData, COOLPARTICLE *cp, double *E,
//		       double PdV, double rho, double ZMetal, double tStep );

/**
 * Updates a particle's internal energy based on cooling and PdV work.
 * Cooling is done according to the equation:
 *      Edot = -rFactor*E + PdV
 * where PdV is the PdV work per time done on the particle and rFactor
 * is given by:
 *      rFactor = dCoolRateEff*r^3/2 = ( (r/r0)^-3/2 )/t_cool
 * where t_cool is the cooling time at r0
 * Solving the cooling equation over a time interval from 0 to t gives:
 *      E(t) = exp(-rFactor*t) (E0 - PdV/rFactor) + PdV/rFactor
 *
 * The Energy (Ecode) is updated, along with clData->rho,PdV,Y_total,rFactor,
 * E,T
 *
 * @brief CoolIntegrateEnergyCode Calculates particle's internal energy
 * due to cooling and PdV work
 * @param cl The COOL structure containing constants/params for the simulation
 * @param clData Per-thread cooling data
 * @param cp Cooling particle
 * @param ECode Particle internal energy in Code units
 * @param PdVCode PdV work per time in code units
 * @param rhoCode Particle density in code units
 * @param ZMetal Unused
 * @param posCode particle position in code units.  Assumes the disk center is at r=0
 * @param tStep Time step to integrate over, in seconds
 */
void CoolIntegrateEnergyCode(COOL *cl, clDerivsData *cData, COOLPARTICLE *cp, double *E, 
			     double PdV, double rho, double ZMetal, double *r, double tStep );
/**
 * @brief Initializes data and sets defaults for a cooling particle
 * @param cp a COOLPARTICLE structure (see cooling_planet.h)
 */
void CoolDefaultParticleData( COOLPARTICLE *cp );

/**
 * @brief Initializes data for a cooling particle, given a COOL struct
 * @param cl a COOL struct which has been initialized with constants defined
 * @param cp a cooling particle COOLPARTICLE struct
 * @param E Pointer to the particle's internal energy
 *
 * Initializes data for a cooling particle, given a COOL struct.
 * Sets Y_total and E in the cooling particle
 */
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


