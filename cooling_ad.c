#ifndef NOCOOLING

/*
 * Cooling code originally written by James Wadsley, McMaster
 * University for GASOLINE.
 *
 * Updated for ChaNGa by Isaac Backus, University of Washington
 */
/*
 * Cooling functions for a planetary disk where the cooling time is
 * proportional to the orbital time.
 *
 * ALGORITHM:
 *
 * Integrating over a time step (see CoolIntegrateEnergyCode) gives:
 *      E(t+dt) = exp(-dt/tcool) (E(t) - PdV*tcool) + PdV*tcool
 * where:
 *      tcool = beta/omega
 * beta is a dimensionless cooling time (see CoolAddParams) and omega is the
 * approximate Keplerian angular velocity, calculated as:
 *      omega = sqrt(G*M/R**3)
 * where M is the total mass of all star particles and R is the distance to
 * the center of mass of star particles.
 *
 * NOTES:
 *
 * The convention is that Cool...() functions are public and cl..()
 * functions are private (which could be enforced if we move this to
 * C++).   The COOL class will contain data which is constant across
 * all processors, and will be contained in a Charm++ NodeGroup.
 * The clDerivsData structure will contain data that changes from
 * particle to particle.
 *
 * By convention, all cooling is done in CGS.  The header file defines macros
 * which convert between code and CGS units [cooling_planet.h].  For example,
 * to convert energy into CGS, do :
 *      E = CoolCodeEnergyToErgPerGm(cl, E);
 * Conversion factors are stored in a COOL struct (see clInitConstants)
 *
  For a list of parameters to set, see CoolAddParams in cooling_planet.c
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "cooling.h"

/**
 * @brief CoolInit initializes a COOL data structure for storing
 *      cooling constants and parameters
 * @return a blank COOL data structure with nTableRead = 0
 */
COOL *CoolInit( )
{
    COOL *cl;
    cl = (COOL *) malloc(sizeof(COOL));
    assert(cl!=NULL);
    return cl;
}

/**
 * Per thread initialization of cooling data.
 * @param cl Initialized COOL structure.
 */
clDerivsData *CoolDerivsInit(COOL *cl)
{
    clDerivsData *Data;

    assert(cl != NULL);
    Data = malloc(sizeof(clDerivsData));
    assert(Data != NULL);

    return Data;
}

///Frees memory and deletes cl
void CoolFinalize(COOL *cl )
{
    free(cl);
}

/**
 * Deallocate memory for per-thread data.
 */
void CoolDerivsFinalize(clDerivsData *clData)
{
    free(clData);
}

/**
 * @brief clInitConstants sets physical constants and parameters
 *          for cooling
 * @param cl a COOL data structure (see CoolInit)
 * @param CoolParam an initialized COOLPARAM structure (see CoolAddParams)
 */
void clInitConstants(COOL *cl, double dGmPerCcUnit,
                     double dComovingGmPerCcUnit, double dErgPerGmUnit,
                     double dSecUnit, double dKpcUnit, COOLPARAM CoolParam)
{
    assert(cl!=NULL);
    /* Units */
    cl->dGmPerCcUnit = dGmPerCcUnit;
    cl->dComovingGmPerCcUnit = dComovingGmPerCcUnit;
    cl->dErgPerGmUnit = dErgPerGmUnit;
    cl->dSecUnit = dSecUnit;
    cl->dErgPerGmPerSecUnit = cl->dErgPerGmUnit / cl->dSecUnit;
    cl->diErgPerGmUnit = 1./dErgPerGmUnit;
    cl->dKpcUnit = dKpcUnit;

    cl->dSmallT = CoolParam.dSmallT;
    cl->dIsoT   = CoolParam.dIsoT;
    cl->bSetIsoT = CoolParam.bSetIsoT;
    cl->mu = CoolParam.dMu;
    cl->gamma = CoolParam.dGamma;
}

/**
 * @brief Calculates thermal energy from temperature and Y (deg. freedom/3)
 * @param Y_Total (deg. freedom)/3
 * @param T Temperature
 * @return Thermal (internal) energy per mass (CGS)
 */
double clThermalEnergy( COOL *cl, COOLPARTICLE *cp, double T, double dens) {
  return 1.0/(cp->gamma-1.0)*CL_Rgascode*T;
}

/* Module Interface routines */

/**
 * @brief CoolAddParams Parses parameters and adds them to the COOLPARAM
 *          struct
 * @param CoolParam A COOLPARAM struct
 * @param prm The param struct
 *
 * The available parameters to be set are:
 *      dBetaCooling : Effective cooling time = t_cool * omega.  Also
 *          equal to 2*pi*(t_cool/T_orbit)
 *      dY_Total : (degrees of freedom)/3 for the gas
 *      dCoolingTmin : Minimum allowed temperature (in Kelvin)
 *      dCoolingTmax : Maximum allowed temperature (in Kelvin)
 */
void CoolAddParams( COOLPARAM *CoolParam, PRM prm ) {
  CoolParam -> dSmallT = 1e0;// prm.dSmallT;
  prmAddParam(prm,"dSmallT",paramDouble,&CoolParam->dSmallT,
              sizeof(double),"smallT", "<Small T> = 1e0");
  CoolParam -> dIsoT = -1e0;// prm.dIsoT;
  prmAddParam(prm,"dIsoT",paramDouble,&CoolParam->dIsoT,
              sizeof(double),"IsoT", "<Iso T> = -1e0");
  CoolParam -> bSetIsoT = 0;// prm.dSmallT;
  prmAddParam(prm,"bSetIsoT",paramInt,&CoolParam->bSetIsoT,
              sizeof(int),"bSetIsoT", "<useIso> = 0");
}

/* Placeholder functions which just need to be defined to make
 * Sph.C and the compiler happy */

void CoolTableReadInfo( COOLPARAM *CoolParam, int cntTable, int *nTableColumns,
                        char *suffix ){ *nTableColumns = 0; }//no tables to be read

void CoolTableRead( COOL *Cool, int nData, void *vData) {}

void CoolInitRatesTable( COOL *cl, COOLPARAM CoolParam ){}

double CoolEdotInstantCode(COOL *cl, COOLPARTICLE *cp, double ECode,
                           double rhoCode, double ZMetal, double *posCode ) { return 0.;}

/* Initialization Routines */

/**
 * @brief Initializes data and sets defaults for a cooling particle
 * @param cp a COOLPARTICLE structure (see cooling_planet.h)
 */
void CoolDefaultParticleData( COOL *cl, COOLPARTICLE *cp ){
  cp->mu = cl->mu;
  cp->gamma = cl->gamma;
}

/**
 * @brief Initializes data for a cooling particle, given a COOL struct
 * @param cl a COOL struct which has been initialized with constants defined
 * @param cp a cooling particle COOLPARTICLE struct
 * @param E Pointer to the particle's internal energy
 *
 * Initializes data for a cooling particle, given a COOL struct.
 * Sets Y_total and E in the cooling particle
 */
void CoolInitEnergyAndParticleData( COOL *cl, COOLPARTICLE *cp, double *E,
                                    double dDensity, double dTemp, double fMetal)
{
    CoolDefaultParticleData( cl, cp);
    *E = 1.0/(cp->gamma - 1.0)*CL_Rgascode*fmax(dTemp, cl->dSmallT)/cp->mu*cl->diErgPerGmUnit;
}

///A hold over from cosmology.  Sets the time which is available to all threads
void CoolSetTime( COOL *cl, double dTime, double z ) {
    cl->z = z;
    cl->dTime = dTime;
}


///Calculates the temperature of a COOLPARTICLE given its energy (in Code units)
double CoolCodeEnergyToTemperature( COOL *cl, COOLPARTICLE *cp, double E,
                                    double fMetal) {
    CoolDefaultParticleData( cl, cp);
    if( cl->bSetIsoT && cl->dIsoT > 0.) 
        return cl->dIsoT; 
    double temp = E*cl->dErgPerGmUnit*(cp->gamma - 1.)*cp->mu/CL_Rgascode;
    return fmax(temp, cl->dSmallT);
}

/**
 * Updates a particle's internal energy based on cooling and PdV work.
 * Cooling is done according to the equation:
 *      Edot = -E/tcool + PdV
 * where PdV is the PdV work per time done on the particle and tcool
 * is given by:
 *      tcool = beta/omega
 * beta is a dimensionless cooling set by the runtime parameter dBetaCooling
 * and omega is the approximate keplerian orbital angular velocity
 * Solving the cooling equation over a time interval from 0 to t gives:
 *      E(t) = exp(-t/tcool) (E0 - PdV*tcool) + PdV*tcool
 *
 * The Energy (Ecode) is updated
 *
 * @brief CoolIntegrateEnergyCode Integrates a particle's internal energy
 * due to cooling and PdV work
 * @param cl The COOL structure containing constants/params for the simulation
 * @param clData [unused]
 * @param cp Cooling particle
 * @param ECode Particle internal energy in Code units
 * @param PdVCode PdV work per time in code units
 * @param rhoCode [unused]
 * @param ZMetal [Unused]
 * @param posCode particle position in code units
 * @param tStep Time step to integrate over, in seconds
 */


void CoolIntegrateEnergyCode(COOL *cl, clDerivsData *cData, COOLPARTICLE *cp, double *ECode,
			     double PdVCode, double rhoCode, double ZMetal, double *posCode, double tStep )
{
  double rho = rhoCode*cl->dGmPerCcUnit;

  double E = CoolCodeEnergyToErgPerGm( cl, *ECode );	/*Convert ergs to GM units, found in the header.*/
  double PdV = CoolCodeWorkToErgPerGmPerSec( cl, PdVCode );
  E += PdV*tStep;
  if( cl->bSetIsoT && cl->dIsoT > 0.) {
    CoolDefaultParticleData( cl, cp);
    double eIso = 1.0/(cp->gamma - 1.)*CL_Rgascode*cl->dIsoT/cp->mu;
    E = eIso;
  }
  *ECode = CoolErgPerGmToCodeEnergy(cl, E);
}

int CoolClean( COOL *cl, COOLPARTICLE *cp, double gamma, double *rhoCode, double *eCode, double *pCode, double *csCode) {

  CoolDefaultParticleData( cl, cp);
  double eFloor = 1.0/(gamma - 1.)*CL_Rgascode*cl->dSmallT/cp->mu*cl->diErgPerGmUnit;
  if( *eCode < eFloor) {
    *eCode = eFloor;
    *csCode = sqrt((gamma - 1.)* eFloor);
    *pCode = *rhoCode * (gamma - 1.)* eFloor/gamma;
    return 1;
  }
  return 0;
}

void CoolCodePressureOnDensitySoundSpeed( COOL *cl, COOLPARTICLE *cp,
    double uPred, double fDensity, double gamma, double gammam1,
    double *PoverRho, double *c ) {
  double E = CoolCodeEnergyToErgPerGm( cl, uPred);
  //double rho = fDensity*cl->dGmPerCcUnit;
  double cs = sqrt(gammam1*E);
  *PoverRho = gammam1*E*cl->diErgPerGmUnit; /*Calculate P/rho & convert*/
  *c = cs/(cl->dKpcUnit*3.0857e21/cl->dSecUnit);
}


#endif /* NOCOOLING */
