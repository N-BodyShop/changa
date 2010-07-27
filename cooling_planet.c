#ifndef NOCOOLING

/*
 * Cooling code originally written by James Wadsley, McMaster
 * University for GASOLINE.
 */
/*
 * Cooling functions for a planetary disk where the cooling time is
 * proportional to the orbital time.
 *
 * The convention is that Cool...() functions are public and cl..()
 * functions are private (which could be enforced if we move this to
 * C++).   The COOL class will contain data which is constant across
 * all processors, and will be contained in a Charm++ NodeGroup.
 * The clDerivsData structure will contain data that changes from
 * particle to particle.  This includes both particle data and the
 * Integrator context.  The integrator context has to be per Treepiece
 * because each will be separately integrating a particle producing
 * separate intermediate data.
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

/* General accuracy target */
#define EPS 1e-5
#define MAXABUNDITERATIONS 20
/* Accuracy for Temperature estimation for E,rho assuming eqm abundances */
#define EPSTEMP 1e-5
#define ETAMINTIMESTEP 1e-4

#include "gsl/gsl_errno.h"
#include "gsl/gsl_odeiv.h"
#include "cooling.h"

/* Integrator constants */

/* When to use to just do a first order step and quit */
/* Good value to use I think */
#define ECHANGEFRACSMALL 1e-4  

/* Accuracy target for intergrators */
#define EPSINTEG  1e-5
#define MAXINTEGITS 20000

#define EMUL (1.01)

COOL *CoolInit( )
{
  COOL *cl;
  cl = (COOL *) malloc(sizeof(COOL));
  assert(cl!=NULL);

  cl->nTableRead = 0; /* Internal Tables read from Files */

  cl->DerivsData = malloc(sizeof(clDerivsData));
  assert(cl->DerivsData != NULL);
  cl->DerivsData->IntegratorContext = 
    StiffInit( EPSINTEG, 1, cl->DerivsData, clDerivs, clJacobn );
  
  return cl;
}

void CoolFinalize(COOL *cl ) 
{
  StiffFinalize(cl->DerivsData->IntegratorContext );
  free(cl->DerivsData);
  free(cl);
}

void clInitConstants(COOL *cl, double dGmPerCcUnit,
		     double dComovingGmPerCcUnit, double dErgPerGmUnit,
		     double dSecUnit, double dKpcUnit, COOLPARAM CoolParam) 
{
  assert(cl!=NULL);
  cl->dGmPerCcUnit = dGmPerCcUnit;
  cl->dComovingGmPerCcUnit = dComovingGmPerCcUnit;
  cl->dErgPerGmUnit = dErgPerGmUnit;
  cl->dSecUnit = dSecUnit;
  cl->dErgPerGmPerSecUnit = cl->dErgPerGmUnit / cl->dSecUnit;
  cl->diErgPerGmUnit = 1./dErgPerGmUnit;
  cl->dKpcUnit = dKpcUnit;
  
  cl->dParam1 = CoolParam.dParam1;
  cl->dParam2 = CoolParam.dParam2;
  cl->dParam3 = CoolParam.dParam3;
  cl->dParam4 = CoolParam.dParam4;
  cl->Y_Total = CoolParam.Y_Total;
  cl->Tmin = CoolParam.dCoolingTmin;
  cl->Tmax = CoolParam.dCoolingTmax;

  /* Derivs Data Struct */
  {
    clDerivsData *Data = cl->DerivsData;

    assert(Data != NULL);

    Data->cl = cl;
    Data->dlnE = (log(EMUL)-log(1/EMUL));
  }
}

#define CL_Rgascode         8.2494e7
#define CL_Eerg_gm_degK     CL_Rgascode
#define CL_ev_degK          1.0/1.1604e4
#define CL_Eerg_gm_ev       CL_Eerg_gm_degK/CL_ev_degK
#define CL_Eerg_gm_degK3_2  1.5*CL_Eerg_gm_degK

/* 
 * Though 13.6eV is lost to the Gas as radiation during H recombination, calculating the
 * Energy using u = E per unit mass = 3/2 n/rho k T requires we don't subtract it there.
 * This formulation is useful because then pressure = 2/3 u rho.
 * Instead we subtract the 13.6eV for H collisional ionization which actually
 * removes no energy from the Gas ( similarly for Helium ) 
 * It also means photoionization doesn't add the 13.6eV, only the excess.
 */
double clThermalEnergy( double Y_Total, double T ) {
  return Y_Total*CL_Eerg_gm_degK3_2*T;

}

double clTemperature( double Y_Total, double E ) {
  return E/(Y_Total*CL_Eerg_gm_degK3_2);
}

double clEdotInstant( COOL *cl, double E, double T, double rho, double rFactor )
{
	double Edot;

	Edot = (T < cl->Tmin ? 0 : -E*rFactor);

	return Edot;
	}

/*
 *  We solve the implicit equation:  
 *  Eout = Ein + dt * Cooling ( Yin, Yout, Ein, Eout )
 *
 *  E erg/gm, PdV erg/gm/s, rho gm/cc, dt sec
 * 
 */

int clDerivs(double x, double *y, double *dydx, void *Data) {
  clDerivsData *d = Data;

  d->E = y[0];
  d->T = clTemperature( d->Y_Total, d->E );
  dydx[0] = clEdotInstant( d->cl, d->E, d->T, d->rho, d->rFactor ) + d->PdV;
  return GSL_SUCCESS;
}

int clJacobn(double x, double y[], double dfdx[], double *dfdy, void *Data) {
  clDerivsData *d = Data;
  double E = y[0],dE;
  const int ndim = 1;

  dfdx[0] = 0;

  /* Approximate dEdt/dE */
  d->E = E*(EMUL);
  d->T = clTemperature( d->Y_Total, d->E );
  dE = clEdotInstant( d->cl, d->E, d->T, d->rho, d->rFactor );

  d->E = E*(1/EMUL);
  d->T = clTemperature( d->Y_Total, d->E );
  dE -= clEdotInstant( d->cl, d->E, d->T, d->rho, d->rFactor );

  dfdy[0*ndim + 0] = dE/(E*d->dlnE);
  return GSL_SUCCESS;
}

void clIntegrateEnergy(COOL *cl, double *E, double PdV, double rho,
		       double Y_Total, double radius, double tStep ) {

  double dEdt,dE,Ein = *E,EMin;
  double t=0,dtused,dtnext,tstop = tStep*(1-1e-8),dtEst;
  clDerivsData *d = cl->DerivsData;
  STIFF *sbs = d->IntegratorContext;
  
  if (tStep <= 0) return;

  d->rho = rho;
  d->PdV = PdV;
  d->Y_Total = Y_Total;
  d->rFactor = cl->dParam1*pow(radius,-3./2.);
  
  EMin = clThermalEnergy( d->Y_Total, cl->Tmin );

  dtnext = tStep;
  gsl_odeiv_evolve_reset(sbs->evolver);
  {
    int its = 0;
    while (t<tstop) {
      double Eold;
      if (its++ > MAXINTEGITS) assert(0);
      d->E = *E;
	  d->T = clTemperature( d->Y_Total, d->E );
	  clDerivs(t, E, &dEdt, d);
      if (fabs(dEdt) > 0) {
		  dtEst = fabs(*E/dEdt);

      /* 
	 Since there is no time dependence and the function is smooth
	 if the changes become very small it must have reached a saddle or
	 an equilibrium.  I'll put my money on Equilibrium and abort 
      */
      /*
      if (tStep-t < ECHANGEFRACSMALL*dtEst) {
	fprintf(stderr,"Aborting -- changes too small\n");
	*E += (tStep-t)*dEdt;
	break;
      }
      if (dtEst < tStep-t) sbs->hMin = dtEst*ETAMINTIMESTEP;
      else sbs->hMin = (tStep-t)*ETAMINTIMESTEP;
      */

	if (dtnext > 0.5*dtEst) dtnext = 0.5*dtEst;
      }
      if (dtnext >= tStep-t) dtnext = tStep-t;
      StiffStep( sbs, E, &dEdt,  &t, dtnext, &Ein, &dtused, &dtnext );
      Eold = *E;
#ifdef ASSERTENEG      
      assert(*E > 0.0);
#else
      if (*E < EMin) {
	*E = EMin;
	break;
      }
#endif    
    }
  }
  /* 
     Note Stored abundances are not necessary with
     this eqm integrator therefore the following operations
     could be moved to output routines 
  */

  d->E = *E;
  d->T = clTemperature( d->Y_Total, d->E );
  if (d->T < cl->Tmin ) {
	  d->T = cl->Tmin;
	  *E = clThermalEnergy( d->Y_Total, d->T );
	  }
}

/* Module Interface routines */

void CoolAddParams( COOLPARAM *CoolParam, PRM prm ) {
	CoolParam->dParam1 = 0;
	prmAddParam(prm,"CooldParam1",paramDouble,&CoolParam->dParam1,
				sizeof(double),"dparam1",
				"<Param1> = 0.0");
	CoolParam->dParam2 = 0;
	prmAddParam(prm,"CooldParam2",paramDouble,&CoolParam->dParam2,
				sizeof(double),"dParam2",
				"<Param2> = 0.0");
	CoolParam->dParam3 = 0;
	prmAddParam(prm,"CooldParam3",paramDouble,&CoolParam->dParam3,
				sizeof(double),"dParam3",
				"<Param3> = 0.0");
	CoolParam->dParam4 = 0;
	prmAddParam(prm,"CooldParam4",paramDouble,&CoolParam->dParam4,
				sizeof(double),"dParam4",
				"<Param4> = 0.0");
	CoolParam->Y_Total = 2;
	prmAddParam(prm,"dY_Total",paramDouble,&CoolParam->Y_Total,
				sizeof(double),"Y_Total",
				"<Y_Total> = 2");
	CoolParam->dCoolingTmin = 10;
	prmAddParam(prm,"dCoolingTmin",paramDouble,&CoolParam->dCoolingTmin,
				sizeof(double),"ctmin",
				"<Minimum Temperature for Cooling> = 10K");
	CoolParam->dCoolingTmax = 1e9;
	prmAddParam(prm,"dCoolingTmax",paramDouble,&CoolParam->dCoolingTmax,
				sizeof(double),"ctmax",
				"<Maximum Temperature for Cooling> = 1e9K");
	}
	
void CoolOutputArray( COOLPARAM *CoolParam, int cnt, int *type, char *suffix ) {
    /* *type = OUT_NULL; */

	}

/* Initialization Routines */

void CoolTableReadInfo( COOLPARAM *CoolParam, int cntTable, int *nTableColumns, char *suffix )
{
   int localcntTable = 0;

   *nTableColumns = 0;
   }

void CoolTableRead( COOL *Cool, int nData, void *vData)
{
   fprintf(stderr," Attempt to initialize non-exitent table in cooling\n");
   assert(0);

   }

void CoolDefaultParticleData( COOLPARTICLE *cp )
{
	cp->Y_Total = 2;
	}

void CoolInitEnergyAndParticleData( COOL *cl, COOLPARTICLE *cp, double *E, double dDensity, double dTemp, double fMetal )
{
	cp->Y_Total = cl->Y_Total;
	*E = clThermalEnergy(cp->Y_Total,dTemp)*cl->diErgPerGmUnit;
	}

void CoolInitRatesTable( COOL *cl, COOLPARAM CoolParam ) {
	}

void CoolSetTime( COOL *cl, double dTime, double z ) {
	cl->z = z;
	cl->dTime = dTime;
	}

/* Output Conversion Routines */

double CoolEnergyToTemperature( COOL *Cool, COOLPARTICLE *cp, double E ) {
	return clTemperature( cp->Y_Total, E );
	}

double CoolCodeEnergyToTemperature( COOL *Cool, COOLPARTICLE *cp, double E,
				    double fMetal) {
	return CoolEnergyToTemperature( Cool, cp, E*Cool->dErgPerGmUnit );
	}

/* Integration Routines */

#define CONVERT_CMPERKPC (3.0857e21)

void CoolIntegrateEnergyCode(COOL *cl, COOLPARTICLE *cp, double *ECode, 
		       double PdVCode, double rhoCode, double ZMetal, double *posCode, double tStep ) {
	double radius;
	
	radius= sqrt(posCode[0]*posCode[0]+posCode[1]*posCode[1]+posCode[2]*posCode[2])
		*cl->dKpcUnit*CONVERT_CMPERKPC;

	*ECode = CoolCodeEnergyToErgPerGm( cl, *ECode );
	clIntegrateEnergy(cl,  ECode, CoolCodeWorkToErgPerGmPerSec( cl, PdVCode ), 
					  CodeDensityToComovingGmPerCc(cl, rhoCode), cp->Y_Total, radius, tStep);
	*ECode = CoolErgPerGmToCodeEnergy(cl, *ECode);

	}

/* Star form function -- do not use */
double CoolHeatingRate( COOL *cl, COOLPARTICLE *cp, double T, double dDensity ) {
	assert(0);
	return 0;
	}

/* Not implemented */
double CoolEdotInstantCode(COOL *cl, COOLPARTICLE *cp, double ECode, 
			  double rhoCode, double ZMetal, double *posCode ) {
    double T,E,rho,Edot;

    assert(0);
    return CoolErgPerGmPerSecToCodeWork( cl, Edot );
    }

#endif /* NOCOOLING */
