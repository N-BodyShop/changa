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
 * For a list of parameters to set, see CoolAddParams in cooling_planet.c
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
    /* Internal Tables read from Files */
    cl->nTableRead = 0;

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

    /* config parameters */
    cl->Y_Total = CoolParam.Y_Total;
    cl->Tmin = CoolParam.dCoolingTmin;
    cl->Tmax = CoolParam.dCoolingTmax;
    cl->beta = CoolParam.dBetaCooling;
}
/* Define physical constants */


/**
 * @brief Calculates thermal energy from temperature and Y (deg. freedom/3)
 * @param Y_Total (deg. freedom)/3
 * @param T Temperature
 * @return Thermal (internal) energy per mass (CGS)
 */
double clThermalEnergy( double Y_Total, double T ) {
    return Y_Total*CL_Eerg_gm_degK3_2*T;

}
/**
 * @brief Calculates temperature from internal energy and Y (deg.freedom/3)
 * @param Y_Total (deg. freedom)/3
 * @param E Internal energy per mass (CGS units)
 * @return Temperature
 */
double clTemperature( double Y_Total, double E ) {
    return E/(Y_Total*CL_Eerg_gm_degK3_2);
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
    CoolParam->dBetaCooling = 1.0;
    prmAddParam(prm, "dBetaCooling", paramDouble, &CoolParam->dBetaCooling,
                sizeof(double), "betacool",
                "<Effective cooling time (tCool*omega)> = 1");
    CoolParam->Y_Total = 2;
    prmAddParam(prm,"dY_Total",paramDouble,&CoolParam->Y_Total,
                sizeof(double),"Y_Total",
                "<Y_Total> = 2");
    CoolParam->dCoolingTmin = 0;
    prmAddParam(prm,"dCoolingTmin",paramDouble,&CoolParam->dCoolingTmin,
                sizeof(double),"ctmin",
                "<Minimum Temperature for Cooling> = 10K");
    CoolParam->dCoolingTmax = 1e9;
    prmAddParam(prm,"dCoolingTmax",paramDouble,&CoolParam->dCoolingTmax,
                sizeof(double),"ctmax",
                "<Maximum Temperature for Cooling> = 1e9K");
}

/* Placeholder functions which just need to be defined to make
 * Sph.C and the compiler happy */

void CoolTableReadInfo( COOLPARAM *CoolParam, int cntTable, int *nTableColumns, char *suffix )
{
    *nTableColumns = 0;
}
void CoolTableRead( COOL *Cool, int nData, void *vData) {}
void CoolInitRatesTable( COOL *cl, COOLPARAM CoolParam ) {}
double CoolEdotInstantCode(COOL *cl, COOLPARTICLE *cp, double ECode,
                           double rhoCode, double ZMetal, double *posCode ) {}

/* Initialization Routines */

/**
 * @brief Initializes data and sets defaults for a cooling particle
 * @param cp a COOLPARTICLE structure (see cooling_planet.h)
 */
void CoolDefaultParticleData( COOLPARTICLE *cp )
{
    cp->Y_Total = 2;
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
void CoolInitEnergyAndParticleData( COOL *cl, COOLPARTICLE *cp, double *E, double dDensity, double dTemp, double fMetal )
{
    cp->Y_Total = cl->Y_Total;
    *E = clThermalEnergy(cp->Y_Total,dTemp)*cl->diErgPerGmUnit;
}

///A hold over from cosmology.  Sets the time which is available to all threads
void CoolSetTime( COOL *cl, double dTime, double z ) {
    cl->z = z;
    cl->dTime = dTime;
}

/**
 * @brief CoolSetStarCM saves the center of mass of the star(s) to the COOL
 * struct cl
 * @param cl COOL struct, which stores cooling data/info
 * @param dCenterOfMass Array(length 4) which contains the star(s) center of
 * mass as the first 3 entries and the total star mass as the final entry
 */
void CoolSetStarCM(COOL *cl, double dCenterOfMass[4]) {
    int i;
    for(i=0; i<4; ++i) {
        cl->dStarCenterOfMass[i] = dCenterOfMass[i];
    }
}

///Calculates the temperature of a COOLPARTICLE given its energy (in Code units)
double CoolCodeEnergyToTemperature( COOL *Cool, COOLPARTICLE *cp, double E,
                                    double fMetal) {
    return clTemperature(cp->Y_Total, E*Cool->dErgPerGmUnit);
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
 * @param posCode particle position in code units.  Assumes the disk center is
 * at r=0
 * @param tStep Time step to integrate over, in seconds
 */
void CoolIntegrateEnergyCode(COOL *cl, clDerivsData *clData, COOLPARTICLE *cp,
                             double *ECode, double PdVCode, double rhoCode,
                             double ZMetal, double *posCode, double tStep )
{
    double radius = 0;
    double omega, tcool, PdV, T;

    /* This is used as a switch to disable cooling before a certain time*/
    if (tStep <= 0) return;
    /* Calculate the radius (in simulation units)
     * Defined as the distance to the center of mass of all the star
     * particles */
    int i=0;
    for(i=0; i<3; ++i) {
        radius += pow(posCode[i] - cl->dStarCenterOfMass[i], 2);
    }
    radius = sqrt(radius);
    /* Calculate approximate keplerian angular velocity */
    omega = sqrt(cl->dStarCenterOfMass[3] / pow(radius,3));
    /* Calculate cooling time (in seconds)*/
    tcool = cl->dSecUnit * cl->beta/omega;
    /* Convert the particle's internal energy to CGS */
    *ECode = CoolCodeEnergyToErgPerGm( cl, *ECode );
    /* Convert PdV work to CGS */
    PdV = CoolCodeWorkToErgPerGmPerSec(cl, PdVCode);
    /* Integrate energy */
    *ECode = exp(-tStep/tcool) * (*ECode - PdV*tcool) + PdV*tcool;
    /* apply minimum temperature cutoff */
    T = clTemperature(cp->Y_Total, *ECode);
    if (T < cl->Tmin) {
        T = cl->Tmin;
        *ECode = clThermalEnergy(cp->Y_Total, T);
    }
    /*Convert energy back to simulation units */
    *ECode = CoolErgPerGmToCodeEnergy(cl, *ECode);
}

#endif /* NOCOOLING */
