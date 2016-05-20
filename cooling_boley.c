#ifndef NOCOOLING

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

/*
 * Planetary disk cooling code as described in Boley 2009, ApJ 695,L53
 * and Boley et al 2010, Icarus 207, 509.  The cooling is calculated
 * from \grad \cdot F \sim (T^4 - T^4_{irr})/(\Delta \tau + 1/\Delta
 * \tau).  This approximates an irradiated disk.
 *
 * This implementation is taken from the GASOLINE implementation by
 * Hayfield.
 */

/* General accuracy target */
#define EPS 1e-5
#define MAXABUNDITERATIONS 20
/* Accuracy for Temperature estimation for E,rho assuming eqm abundances */
#define EPSTEMP 1e-5
#define ETAMINTIMESTEP 1e-4

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

  cl->rossTab = NULL;
  cl->plckTab = NULL;

  cl->nTableRead = 0; /* Internal Tables read from Files */

  return cl;
}

/**
 * Per thread initialization of cooling
 * @param cl Initialized COOL structure.
 */
clDerivsData *CoolDerivsInit(COOL *cl)
{
    clDerivsData *Data;
    double dEMin;

    assert(cl != NULL);
    Data = malloc(sizeof(clDerivsData));
    assert(Data != NULL);
    Data->IntegratorContext = StiffInit( EPSINTEG, 1, Data, clDerivs);
    Data->cl = cl;
    dEMin =  clThermalEnergy(cl->Y_Total, cl->Tmin);
    StiffSetYMin(Data->IntegratorContext, &dEMin);
    return Data;
  }
  

void CoolFinalize(COOL *cl ) 
{
  free(cl->DerivsData);

  if(cl->rossTab != NULL){
    free(cl->rossTab[0]);
    free(cl->rossTab);
  }

  if(cl->plckTab != NULL){
    free(cl->plckTab[0]);
    free(cl->plckTab);
  }

  free(cl);
}

/**
 * Deallocate memory for per-thread data.
 */
void CoolDerivsFinalize(clDerivsData *clData)
{
    StiffFinalize(clData->IntegratorContext );
    free(clData);
    }

void clInitConstants( COOL *cl, double dGmPerCcUnit, double dComovingGmPerCcUnit, 
		double dErgPerGmUnit, double dSecUnit, double dKpcUnit, COOLPARAM CoolParam) 
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
  cl->dRhoCutoff = CoolParam.dRhoCutoff;
  cl->dParam3 = CoolParam.dParam3;
  cl->dParam4 = CoolParam.dParam4;
  cl->Y_Total = CoolParam.Y_Total;
  cl->Tmin = CoolParam.dCoolingTmin;
  cl->Tmax = CoolParam.dCoolingTmax;

  cl->rossTab = NULL;
  cl->plckTab = NULL;

  if((CoolParam.achRossName[0] != '\0') && (CoolParam.achPlckName[0] != '\0')){
    int i,j, fstat;
    double * tt = NULL;
    FILE * fp;

    cl->rossTab = (double **)malloc(sizeof(double *)*CL_TABRC);
    assert(cl->rossTab != NULL);

    tt = (double *)malloc(sizeof(double)*CL_TABRC2);
    assert(tt != NULL);

    for(i = 0; i < CL_TABRC; i++){
      cl->rossTab[i] = tt+i*CL_TABRC;
    }

    cl->plckTab = (double **)malloc(sizeof(double *)*CL_TABRC);
    assert(cl->plckTab != NULL);

    tt = (double *)malloc(sizeof(double)*CL_TABRC2);
    assert(tt != NULL);

    for(i = 0; i < CL_TABRC; i++){
      cl->plckTab[i] = tt+i*CL_TABRC;
    }

    // read tables - rosseland
    assert((fp = fopen(CoolParam.achRossName, "r")) != NULL);
    
    for(j = 0; j < CL_TABRC; j++)
      for(i = 0; i < CL_TABRC; i++){
        fstat = fscanf(fp, "%lf", cl->rossTab[j]+i);
        assert(fstat != EOF);
      }
    
    
    fclose(fp);

    // read tables - planck
    assert((fp = fopen(CoolParam.achPlckName, "r")) != NULL);
    
    for(j = 0; j < CL_TABRC; j++)
      for(i = 0; i < CL_TABRC; i++){
        fstat = fscanf(fp, "%lf", cl->plckTab[j]+i);
        assert(fstat != EOF);
      }

    fclose(fp);
  }
  else {
      fprintf(stderr,"Need to specify achRossName and achPlkName\n");
      assert(0);
      }
  
  cl->t4i = cl->Tmin*cl->Tmin;
  cl->t4i = cl->t4i*cl->t4i;
}

#define CL_Rgascode         8.2494e7
#define CL_Eerg_gm_degK     CL_Rgascode
#define CL_ev_degK          1.0/1.1604e4
#define CL_Eerg_gm_ev       CL_Eerg_gm_degK/CL_ev_degK
#define CL_Eerg_gm_degK3_2  1.5*CL_Eerg_gm_degK
#define CL_SIGMA 5.6704e-05

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

double clEdotInstant( COOL *cl, double E, double T, double rho, double rFactor,
                      double *dDeltaTau, 
                      double *dEdotHeat, double *dEdotCool)
{ 
  // I'm not sure this is 'good' but tries to handle negative E, T without crashing
  // Also apparently need to handle nan's?

  double Edot;

  // NB that we may need to use the cooling timestep 
  // along with the scheme
  // assume here that the gas ideal
  double kr, kp, t1, t2;

  const double ilt = 1.0/log(10.0);
  const double gfac = pow(36.0*M_PI, 1.0/3.0);
  const double newT = (T <= cl->Tmin) ? cl->Tmin : T;
  const double p = cl->Y_Total*CL_Rgascode*rho*newT;

  *dEdotHeat = 0.0;

  double tpl = ilt*log(p);
  double ttl = ilt*log(newT);
  
  tpl = (tpl - CL_TABPMIN)/CL_TABDP;
  ttl = (ttl - CL_TABTMIN)/CL_TABDT;

  if((tpl < 0.0) || (tpl > CL_TABRCM2) || (tpl != tpl) ||
     (ttl < 0.0) || (ttl > CL_TABRCM2) || (ttl != ttl) ||
     (rho > cl->dRhoCutoff)){
    *dEdotCool = 0.0;
    return 0.0;
  }

  //  tpl = (tpl < 0) ? 0.0 : 
  //    ((tpl > CL_TABRCM2) ? CL_TABRCM2 : tpl);

  //  ttl = (ttl < 0) ? 0.0 : 
  //    ((ttl > CL_TABRCM2) ? CL_TABRCM2 : ttl);

  const int tipl = (int)tpl;
  const int titl = (int)ttl;
  const double dp = tpl-tipl;
  const double dt = ttl-titl;

  t1 = cl->rossTab[tipl][titl] + 
    (cl->rossTab[tipl][titl+1] - cl->rossTab[tipl][titl])*dt;

  t2 = cl->rossTab[tipl+1][titl] + 
    (cl->rossTab[tipl+1][titl+1] - cl->rossTab[tipl+1][titl])*dt;

  kr = pow(10.0, t1 + (t2 - t1)*dp);

  t1 = cl->plckTab[tipl][titl] + 
    (cl->plckTab[tipl][titl+1] - cl->plckTab[tipl][titl])*dt;

  t2 = cl->plckTab[tipl+1][titl] + 
    (cl->plckTab[tipl+1][titl+1] - cl->plckTab[tipl+1][titl])*dt;

  kp = pow(10.0, t1 + (t2 - t1)*dp);

  // now t1 is delta tau planck, t2 delta tau ross
  // rFactor is (m/rho)^(1/3)

  t1 = rho*kp*rFactor;
  /*
   * See Boley 2009, ApJ 695, L53, following eq. 2.  where 
   * the local optical depth across a cell is calculated as
   * \Delta \tau = \rho(\kappa_R(1 - exp(-2\Delta \tau_P))
   *                    + \kappa_P exp(-2 \Delta \tau_P)) s
   */
  t2 = rho*(kp*exp(-2.0*t1) - kr*expm1(-2.0*t1))*rFactor;
  *dDeltaTau = t2;
  // finally t2 = f(tau)
#ifdef FTAU_EXP
  t2 *= 1.36*exp(-t2);
#else //FTAU_EXP
  t2 += 1.0/t2;
  t2 = 1.0/t2;
#endif //FTAU_EXP

  // t1 is now T^4
  t1 = newT*newT;
  t1 = t1*t1;

  Edot = -gfac*CL_SIGMA*(t1-cl->t4i)/(rho*rFactor)*t2;
  *dEdotCool = -Edot;

  return Edot;
}

/*
 *  We solve the implicit equation:  
 *  Eout = Ein + dt * Cooling ( Yin, Yout, Ein, Eout )
 *
 *  E erg/gm, PdV erg/gm/s, rho gm/cc, dt sec
 * 
 */

void clDerivs(double x, const double *y, double *dHeat, double *dCool,
              void *Data) {
  clDerivsData *d = Data;

  d->E = y[0];
  d->T = clTemperature( d->Y_Total, d->E );
  clEdotInstant( d->cl, d->E, d->T, d->rho, d->rFactor, &d->dDeltaTau,
                 dHeat, dCool );
  if(d->PdV > 0.0)
      *dHeat += d->PdV;
  else
      *dCool -= d->PdV;
}

#if 0
int clJacobn( double x, const double y[], double dfdx[], double *dfdy,
	       void *Data) {
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
#endif

void clIntegrateEnergy(COOL *cl, clDerivsData *clData, double *E, 
		       double PdV, double rho, double Y_Total,
                       double radius, /* radius of particle */
                       double tStep ) {

  double EMin;
  double t=0;
  clDerivsData *d = clData;
  STIFF *sbs = d->IntegratorContext;
  
  if (tStep <= 0) return;

  d->rho = rho;
  d->PdV = PdV;
  d->Y_Total = Y_Total;

  d->rFactor = radius;

  EMin = clThermalEnergy( d->Y_Total, cl->Tmin );

  {
      StiffStep( sbs, E, t, tStep);
#ifdef ASSERTENEG      
      assert(*E > 0.0);
#else
      if (*E < EMin) {
	*E = EMin;
      }
#endif
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
  // cooling time (orbital frequency at 1 cm) or fixed cooling frequency
	CoolParam->dParam1 = 0;
	prmAddParam(prm,"CooldParam1",paramDouble,&CoolParam->dParam1,
				sizeof(double),"dparam1",
				"<Param1> = 0.0");
  // cooling cutoff density
        CoolParam->dRhoCutoff = 1e38;
        prmAddParam(prm,"dCoolRhoCutoff",paramDouble,&CoolParam->dRhoCutoff,
                    sizeof(double),"dCoolRhoCutoff",
                    "<no cooling above this density (gm/cc)> = HUGE");
  // cooling gamma
	CoolParam->dParam3 = 5.0/3.0;
#if 0
        // Unimplemented, and fixed at 5/3 for now
	prmAddParam(prm,"CooldParam3",paramDouble,&CoolParam->dParam3,
				sizeof(double),"dParam3",
				"<Param3> = 5/3");
#endif
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

	CoolParam->achRossName[0] = '\0';
	prmAddParam(prm,"achRossName",paramString,CoolParam->achRossName,
				256,"ross",
				"<Input filename for Rosseland opacities>");

	CoolParam->achPlckName[0] = '\0';
	prmAddParam(prm,"achPlckName",paramString,CoolParam->achPlckName,
				256,"plck",
				"<Input filename for Planck opacities>");

	}
	
void CoolOutputArray( COOLPARAM *CoolParam, int cnt, int *type, char *suffix ) {
#if 0
	*type = OUT_NULL;
#endif
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

void CoolInitEnergyAndParticleData( COOL *cl, COOLPARTICLE *cp, double *E, double dDensity, double dTemp, double dMetals )
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
				    double fMetals ) {
	return CoolEnergyToTemperature( Cool, cp, E*Cool->dErgPerGmUnit );
	}

/* Integration Routines */

#define CONVERT_CMPERKPC (3.0857e21)

void CoolIntegrateEnergyCode(COOL *cl, clDerivsData *clData, COOLPARTICLE *cp,
			     double *ECode, double PdVCode, double rhoCode,
			     double ZMetal, double *posCode, double tStep ) {
	double radius;          /* radius of a particle to calculate tau */
	
	radius = cp->mrho*cl->dKpcUnit*CONVERT_CMPERKPC;

	*ECode = CoolCodeEnergyToErgPerGm( cl, *ECode );
	clIntegrateEnergy(cl, clData, ECode, CoolCodeWorkToErgPerGmPerSec( cl, PdVCode ), 
					  CodeDensityToComovingGmPerCc(cl, rhoCode), cp->Y_Total, radius, tStep);
	*ECode = CoolErgPerGmToCodeEnergy(cl, *ECode);
        cp->dDeltaTau = clData->dDeltaTau;
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
