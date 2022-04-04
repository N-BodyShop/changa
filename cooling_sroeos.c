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
#include "nuc_eos_for_c.h"

#define K_TO_MEV 8.625e-11
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

    cl->heat0 = 1.544e20;
    cl->cool0 = 1.399e20;
    cl->Tnu = CoolParam.dTnu;
    cl->L0 = CoolParam.dLnu/1e52;
    cl->includeHeating = CoolParam.bIncludeHeating;
}

/**
 * @brief Calculates thermal energy from temperature and Y (deg. freedom/3)
 * @param Y_Total (deg. freedom)/3
 * @param T Temperature
 * @return Thermal (internal) energy per mass (CGS)
 */
double clThermalEnergy( COOL *cl, double T, double dens, double ye) {

  double cs, press,  gamma, ener;
  cl_sro_eos_rho_T(cl, ye, dens, T, &press, &ener, &cs, &gamma);

  return ener;
  //return (xeps+cl->dSroEosEnergyShift)*INVEPSGF;
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
void CoolAddParams( COOLPARAM *CoolParam, PRM prm )  { 
  CoolParam -> EosFileName[0] = '\0';
  prmAddParam(prm,"EosFileName",paramString,CoolParam->EosFileName,
              512, "eos_file", "eos file name");
  CoolParam -> dSmallT = 1e8;// prm.dSmallT;
  prmAddParam(prm,"dSmallT",paramDouble,&CoolParam->dSmallT,
              sizeof(double),"smallT", "<Small T> = 1e8");
  CoolParam -> dTnu = 4.;// in MeV;
  prmAddParam(prm,"dTnu",paramDouble,&CoolParam->dTnu,
              sizeof(double),"Tnu", "<Tnu> = 4 MeV");
  CoolParam -> dLnu = 1e52;// in ergs/s;
  prmAddParam(prm,"dLnu",paramDouble,&CoolParam->dLnu,
              sizeof(double),"Lnu", "<Lnu> = 1e52");
  CoolParam -> bIncludeHeating = 0;
  prmAddParam(prm,"bIncludeHeating",paramInt,&CoolParam->bIncludeHeating,
              sizeof(int),"IncludeHeating", "<IncludeHeating> = 0");
}

/* Placeholder functions which just need to be defined to make
 * Sph.C and the compiler happy */

void CoolTableReadInfo( COOLPARAM *CoolParam, int cntTable, int *nTableColumns,
                        char *suffix ){ *nTableColumns = 0; }//no tables to be read

void CoolTableRead( COOL *Cool, int nData, void *vData) {}

void CoolInitRatesTable( COOL *cl, COOLPARAM CoolParam ){
  nuc_eos_C_ReadTable_for_c( CoolParam.EosFileName);
  cl->dSroEosEnergyShift = nuc_eos_energy_shift_for_c();
  cl->dSroEosRhoMin = nuc_eos_rhomin_for_c()/RHOGF;
  cl->dSroEosTempMin = nuc_eos_tempmin_for_c()/K_TO_MEV;
  cl->dSroEosYeMax = nuc_eos_yemax_for_c();
}

double CoolEdotInstantCode(COOL *cl, COOLPARTICLE *cp, double ECode,
                           double rhoCode, double ZMetal, double *posCode ) { return 0.;}

/* Initialization Routines */

/**
 * @brief Initializes data and sets defaults for a cooling particle
 * @param cp a COOLPARTICLE structure (see cooling_planet.h)
 */
void CoolDefaultParticleData( COOLPARTICLE *cp ){
  cp->ye = 0.1;
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
                                    double dDensity, double dTemp, double ye)
{
    CoolDefaultParticleData( cp);
    cp->ye = fmin(ye, 0.5);
    *E = clThermalEnergy( cl, dTemp, dDensity*cl->dGmPerCcUnit, cp->ye) * cl->diErgPerGmUnit;
}

///A hold over from cosmology.  Sets the time which is available to all threads
void CoolSetTime( COOL *cl, double dTime, double z ) {
    cl->z = z;
    cl->dTime = dTime;
}


///Calculates the temperature of a COOLPARTICLE given its energy (in Code units)
double CoolCodeEnergyToTemperature( COOL *Cool, COOLPARTICLE *cp, double ECode, double rhoCode,
                                    double fMetal) {
  double E = CoolCodeEnergyToErgPerGm( Cool, ECode);
  double rho = rhoCode*Cool->dGmPerCcUnit;

  double cs, press, T,  gamma;

  cl_sro_eos_rho_ener(Cool, cp->ye, rho, E, &press, &T, &cs, &gamma);
  
  return T; 
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
  double heating = 0.;
  double cooling = 0.;
  if( cl->includeHeating) {
    double tau_nu = rho/1e11; 
    double expTau = exp(-fmin(tau_nu,20.));
    double T2 = (cl->Tnu/4.)*(cl->Tnu/4.);
    double cm = cl->dKpcUnit*3.0857e21;
    double r2 = posCode[0]*posCode[0] + posCode[1]*posCode[1] + posCode[2]*posCode[2];
    r2 *= cm*cm; 
    double TMeV = CoolCodeEnergyToTemperature( cl, cp, *ECode, rhoCode, -1.)*K_TO_MEV;
    cooling = cl->cool0*TMeV*TMeV/4. *expTau;
    heating = cl->heat0*cl->L0*T2*(1e7*1e7)/r2*expTau;
  }
  // compute the temperature
  E += PdV*tStep + (heating - cooling)*tStep;
  *ECode = CoolErgPerGmToCodeEnergy(cl, E);
}

int CoolClean( COOL *cl, COOLPARTICLE *cp, double gamma, double *rhoCode, double *eCode, double *pCode, double *csCode) {
  double rho = *rhoCode * cl->dGmPerCcUnit;
  double cs, press, T,  gam1, ener;
  T = cl->dSmallT;
  cl_sro_eos_rho_T(cl, cp->ye, rho, T, &press, &ener, &cs, &gam1);

  double eFloorCode = ener*cl->diErgPerGmUnit;
  if( *eCode < eFloorCode) { 
    *eCode = eFloorCode; /*Calculate P/rho & convert*/
    *csCode = cs/(cl->dKpcUnit*3.0857e21/cl->dSecUnit);
    *pCode = press/(cl->dGmPerCcUnit*cl->dErgPerGmUnit);
    return 1;
  }
  return 0;
}

void CoolCodePressureDensity( COOL *cl, COOLPARTICLE *cp,
    double press, double fDensity, double *ie, double *c, double *gamma) {

  double rho = fDensity*cl->dGmPerCcUnit;
  double pressure = press*cl->dGmPerCcUnit*cl->dErgPerGmUnit;
  double ienergy = pressure/rho/(*gamma - 1.);
  double soundspeed = sqrt(*gamma*pressure/rho);
  *ie = ienergy*cl->dGmPerCcUnit*cl->dErgPerGmUnit;
  *c = soundspeed/(cl->dKpcUnit*3.0857e21/cl->dSecUnit);
}

void CoolCodePressureOnDensitySoundSpeed( COOL *cl, COOLPARTICLE *cp,
    double uPred, double fDensity, double gam, double gamm1,
    double *PoverRho, double *c ) {
  double rho = fDensity*cl->dGmPerCcUnit;
  double cs, press, T,  gamma;
  const double Efloor = clThermalEnergy( cl, cl->dSmallT, rho, cp->ye);
  const double Ein = CoolCodeEnergyToErgPerGm( cl, uPred);
  double E = fmax(Efloor, Ein);

  cl_sro_eos_rho_ener(cl, cp->ye, rho, E, &press, &T, &cs, &gamma);
  *PoverRho = press/rho*cl->diErgPerGmUnit; /*Calculate P/rho & convert*/
  *c = cs/(cl->dKpcUnit*3.0857e21/cl->dSecUnit);
  //if( rho > 1e13) {
  //   printf("press: %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e\n", rho, gamma, press, Ein, Efloor, *PoverRho);
  //}
  if( isnan(gamma)) {
     printf("cooling_sroeos.c: gamma is nan %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %.3e %.3e %.3e\n", gamma, rho, press, cp->ye, T, E, Efloor, Ein, cs);
     gamma = 2.;
     *c = sqrt(gamma* (*PoverRho));
  }
}

void cl_sro_eos_rho_T(COOL *cl, const double ye, const double rho, double T, double *press, double *E, double *cs, double *gamma) {

  if( rho < 1.2*cl->dSroEosRhoMin || ye >= cl->dSroEosYeMax) {
    cl_sro_eos_low( 1, cl, rho, press, &T, E, cs);
    return;
  }
  int n = 1, keyerr = 0, anyerr = 0;
  double xrho = rho * RHOGF;
  double xeps, xprs,xent,xcs2,xdedt,xdpderho,xdpdrhoe,xmunu;

  double xtemp = T*K_TO_MEV;
  double xprec = 1e-4;
  const double xye = fmin(ye, 0.5);

  nuc_eos_m_kt1_short_for_c(&n,&xrho,&xtemp,&xye,
		      &xeps,&xprs,&xent,&xcs2,&xdedt,
		      &xdpderho,&xdpdrhoe,&xmunu,
		      &keyerr,&anyerr);

  const double cspeed = sqrt(xcs2)*CL_SpeedOfLight;
  const double p = xprs * INVPRESSGF;
  const double gam1 = cspeed*cspeed/(p/rho);

  *E = (xeps+cl->dSroEosEnergyShift)*INVEPSGF;
  *cs = cspeed;
  *press = p;
  *gamma = gam1;

}


void cl_sro_eos_rho_ener(COOL *cl, const double ye, const double rho, const double E, double *press, double *T, double *cs, double *gamma) {

  //printf("Min Rho = %.3e\n", cl->dSroEosRhoMin);
  if( rho < 1.2*cl->dSroEosRhoMin || ye >= cl->dSroEosYeMax) {
    double ener = E;
    cl_sro_eos_low( 0, cl, rho, press, T, &ener, cs);
    return;
  }
  int n = 1, keyerr = 0, anyerr = 0;
  double xrho = rho * RHOGF;
  double xprs,xent,xcs2,xdedt,xdpderho,xdpdrhoe,xmunu;

  double T_guess = 2.0*cl->dSmallT;
  double xtemp = T_guess*K_TO_MEV;
  double xprec = 1e-4;
  const double xye = fmin(ye, 0.5);

  double xeps = E*EPSGF-cl->dSroEosEnergyShift;
  
  nuc_eos_m_kt0_short_for_c(&n,&xrho,&xtemp,&xye,
		      &xeps,&xprs,&xent,&xcs2,&xdedt,
		      &xdpderho,&xdpdrhoe,&xmunu,&xprec,
		      &keyerr,&anyerr);

  const double cspeed = sqrt(xcs2)*CL_SpeedOfLight;
  const double p = xprs * INVPRESSGF;
  const double gam1 = cspeed*cspeed/(p/rho);
  const double temp = xtemp/K_TO_MEV;
  if( gam1 < 1e-3) {
     printf("cooling_sroeos.c: gamma = %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %.3e\n", gam1, ye, E, temp, cspeed, p, rho);
     printf("%d %d \n", keyerr, anyerr);
  }

  *cs = cspeed;
  *press = p;
  *gamma = gam1;
  *T = temp;
}

void cl_sro_eos_low( int mode, COOL *cl, const double rho, double *press, double *T, double *ener, double *cs) {
  const double ONE_HALF_FOUR_THIRDS = 0.3968511799;
  const double idealK1 =  1.2435e15 * ONE_HALF_FOUR_THIRDS;
  const double idealGamma = 1.41;
  double p, E;
  if(mode == 1) {
    //    energy wanted
    p=idealK1*pow(rho, idealGamma);
    E=p/rho/(idealGamma-1.);
  }
  else {
    E = *ener;
    p = (idealGamma - 1.0) * rho * E;
    *T = cl->dSmallT;
  }

  *cs=sqrt(idealGamma* (p)/rho);
  *ener = E;
  *press = p;
}

#endif /* NOCOOLING */
