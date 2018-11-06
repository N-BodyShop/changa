#ifndef NOCOOLING

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

/*
 * Cooling code originally written by James Wadsley, McMaster
 * University for GASOLINE.
 */
/* Usage/Interfaces:
   Functions starting with

   Cool are public: intended to be used by outside callers

   cl are private: only intended for using by cooling routine itself
   * The COOL class will contain data which is constant across
 * all processors, and will be contained in a Charm++ NodeGroup.
 * The clDerivsData structure will contain data that changes from
 * particle to particle.  This includes both particle data and the
 * Integrator context.  The integrator context has to be per Treepiece
 * because each will be separately integrating a particle producing
 * separate intermediate data.
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

#define NDOT_MIN  1e-60
  
#define CL_eHI     (13.60*CL_eV_erg)
#define CL_eHeI    (24.59*CL_eV_erg)
#define CL_eHeII   (54.42*CL_eV_erg)
#define CL_E2HeII  (3.0*13.6*CL_eV_erg)

#define EMUL (1.0001)

COOL *CoolInit( )
{
  COOL *cl;
  cl = (COOL *) malloc(sizeof(COOL));
  assert(cl!=NULL);

  cl->nUV = 0;
  cl->UV = NULL; 

  cl->nTable = 0;
  cl->RT = NULL;

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
    Data->Y_Total0 = (cl->Y_H+cl->Y_He)*.9999; /* neutral */
    Data->Y_Total1 = (cl->Y_eMAX+cl->Y_H+cl->Y_He)*1.0001; /* Full Ionization */
    dEMin =  clThermalEnergy(Data->Y_Total0, cl->TMin);
    StiffSetYMin(Data->IntegratorContext, &dEMin);
    return Data;
    }

void CoolFinalize(COOL *cl ) 
{
  if (cl->UV != NULL) free(cl->UV);
  if (cl->RT != NULL) free(cl->RT);
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
  
  cl->Y_H = 1.0-CoolParam.dMassFracHelium;
  cl->Y_He = CoolParam.dMassFracHelium/4;
  cl->Y_eMAX = cl->Y_H + cl->Y_He*2;

  cl->bUV = CoolParam.bUV;
  cl->bUVTableUsesTime = CoolParam.bUVTableUsesTime;
  cl->bUVTableLinear = CoolParam.bUVTableUsesTime; /* Linear if using time */
  cl->bLowTCool = CoolParam.bLowTCool;
  cl->bSelfShield = CoolParam.bSelfShield;

}

void clInitUV(COOL *cl, int nTableColumns, int nTableRows, double *dTableData )
{
    int i;

    assert(cl!=NULL);
    assert(cl->UV == NULL);

	assert(nTableColumns == 7);

    cl->nUV = nTableRows;
    cl->UV = (UVSPECTRUM *) malloc(nTableRows*sizeof(UVSPECTRUM));    
	assert(cl->UV!=NULL);
    
    for (i=0;i<nTableRows;i++) {
		(cl->UV)[i].zTime = dTableData[i*nTableColumns];
		
		(cl->UV)[i].Rate_Phot_HI = dTableData[i*nTableColumns+1];
		(cl->UV)[i].Rate_Phot_HeI = dTableData[i*nTableColumns+2];
		(cl->UV)[i].Rate_Phot_HeII = dTableData[i*nTableColumns+3];
		
		(cl->UV)[i].Heat_Phot_HI = dTableData[i*nTableColumns+4];
		(cl->UV)[i].Heat_Phot_HeI = dTableData[i*nTableColumns+5];
		(cl->UV)[i].Heat_Phot_HeII = dTableData[i*nTableColumns+6];

		/* Make sure the heating is in units of ergs per ionization */
		assert( (cl->UV)[i].Heat_Phot_HI>1e-15 && (cl->UV)[i].Heat_Phot_HI<1e-10);

		if (i) assert (((cl->UV)[i-1].zTime > (cl->UV)[i].zTime 
						&& !cl->bUVTableUsesTime) ||
					   ((cl->UV)[i-1].zTime < (cl->UV)[i].zTime 
						&& cl->bUVTableUsesTime));
		}
    }
   
void clInitRatesTable( COOL *cl, double TMin, double TMax, int nTable ) {
/*-----------------------------------------------------------------
 *      Function of ln(T) tables:
 * cl assumed to be predefined (allocated)
 * A table spacing of 1.23e-3 in log(T) eg. nTable=15001 10->1e9 K
 * Maximum 1% errors except at discontinuities in the functions.
 * Storage for such a table is (15001*15*4|8b) => 0.9|1.7 Mb
 *-----------------------------------------------------------------*/
  int i;
  double DeltaTln, Tln, T;

  assert(cl!=NULL);
  assert(cl->RT == NULL);

  cl->R.Cool_Coll_HI = CL_eHI*CL_B_gm;
  cl->R.Cool_Coll_HeI = CL_eHeI*CL_B_gm;
  cl->R.Cool_Coll_HeII = CL_eHeII*CL_B_gm;
  cl->R.Cool_Diel_HeII = (CL_E2HeII+CL_eHeI)*CL_B_gm;

  cl->nTable = nTable;
  cl->TMin = TMin;
  cl->TMax = TMax;
  cl->TlnMin = log( TMin );
  cl->TlnMax = log( TMax );
  DeltaTln = ( cl->TlnMax-cl->TlnMin )/( nTable - 1 );
  cl->rDeltaTln = 1./DeltaTln;
  cl->RT = (RATES_T *) malloc( nTable * sizeof(RATES_T) );
  assert(cl->RT != NULL);

  for ( i=0; i<nTable; i++ ) {
    Tln = cl->TlnMin + DeltaTln*(i-1);
    T=exp(Tln);

    (cl->RT+i)->Rate_Coll_HI = clRateCollHI( T );
    if ( (cl->RT+i)->Rate_Coll_HI < CL_RT_MIN ) (cl->RT+i)->Rate_Coll_HI = CL_RT_MIN;
    (cl->RT+i)->Rate_Coll_HeI = clRateCollHeI( T );
    if ( (cl->RT+i)->Rate_Coll_HeI < CL_RT_MIN ) (cl->RT+i)->Rate_Coll_HeI = CL_RT_MIN;
    (cl->RT+i)->Rate_Coll_HeII = clRateCollHeII( T );
    if ( (cl->RT+i)->Rate_Coll_HeII < CL_RT_MIN ) (cl->RT+i)->Rate_Coll_HeII = CL_RT_MIN;

    (cl->RT+i)->Rate_Radr_HII = clRateRadrHII( T );
    (cl->RT+i)->Rate_Radr_HeII = clRateRadrHeII( T );
    (cl->RT+i)->Rate_Radr_HeIII = clRateRadrHeIII( T );
    (cl->RT+i)->Rate_Diel_HeII = clRateDielHeII( T );

    (cl->RT+i)->Cool_Brem_1 = clCoolBrem1( T );
    (cl->RT+i)->Cool_Brem_2 = clCoolBrem2( T );
    (cl->RT+i)->Cool_Radr_HII = clCoolRadrHII( T );
    (cl->RT+i)->Cool_Radr_HeII = clCoolRadrHeII( T );
    (cl->RT+i)->Cool_Radr_HeIII = clCoolRadrHeIII( T );
    (cl->RT+i)->Cool_Line_HI = clCoolLineHI( T );
    (cl->RT+i)->Cool_Line_HeI = clCoolLineHeI( T );
    (cl->RT+i)->Cool_Line_HeII = clCoolLineHeII( T );
    (cl->RT+i)->Cool_LowT = clCoolLowT( T );
  }    
}

void clRatesTableError( COOL *cl ) {
  /* This estimates the error for a table of half the size of
   * the current one.
   * The collisional, dielectric and line cooling rates can all go to
   * zero (minimum) and will even for double precision (exp(-1e5/T)) 
   * the critical values that must never get down to zero
   * are the radiative recombination rates: Check this before using
   * CL_RT_FLOAT float for the table.
   */
  int i,j,ierr[15];
  double min[15],max[15];
  double err,maxerr[15];
  CL_RT_FLOAT *p,*p1,*p2;

  for (j=0;j<15;j++) {
    maxerr[j]=0.0;
    min[j]=1e300;
    max[j]=-1e300;
  }
  for ( i=1; i<cl->nTable-1; i+=2 ) {
    p = (CL_RT_FLOAT *) &((cl->RT+i)->Rate_Coll_HI);
    p1 = (CL_RT_FLOAT *) &((cl->RT+i-1)->Rate_Coll_HI);
    p2 = (CL_RT_FLOAT *) &((cl->RT+i+1)->Rate_Coll_HI);
    if (i==10001) {
      printf(" Comp %i %i %e %e\n",i,(int) sizeof(CL_RT_FLOAT),p[1],(cl->RT+i)->Rate_Coll_HeI);
    }
    for (j=0;j<15;j++) {
      if (p[j] < min[j]) min[j]=p[j];
      if (p[j] > max[j]) max[j]=p[j];
      err= fabs((0.5*(p1[j]+p2[j])-p[j])/(p[j]+1e-30));
      if (err > maxerr[j] ) {
	maxerr[j]=err;
	ierr[j]=i;
      }
    }
  }  
  for (j=0;j<15;j++) {
    printf("Col %i  Max Error: %e at T=%e dlnT=%e (min %e max %e)\n",j,
	 maxerr[j],exp(cl->TlnMin+ierr[j]/cl->rDeltaTln),2/cl->rDeltaTln,min[j],max[j]);
  }
}

#define CL_Ccomp0 0.565e-9
#define CL_Tcmb0  2.735
#define CL_Ccomp  (CL_Ccomp0*CL_Tcmb0)

void clRatesRedshift( COOL *cl, double zIn, double dTimeIn ) {
  int i;
  double xx;
  double zTime;
  UVSPECTRUM *UV,*UV0;

  cl->z = zIn;
  cl->dTime = dTimeIn;
  cl->dComovingGmPerCcUnit = cl->dGmPerCcUnit*pow(1.+zIn,3.);
     
  cl->R.Cool_Comp = pow((1+zIn)*CL_Ccomp,4.0)*CL_B_gm;
  cl->R.Tcmb = CL_Tcmb0*(1+zIn);
  cl->R.Cool_LowTFactor = (cl->bLowTCool ? CL_B_gm*cl->Y_H*cl->Y_H/0.001 : 0 );

  /* Photo-Ionization rates */

  UV = cl->UV;

  if (cl->bUV) {
	  assert( UV != NULL );
	  if (cl->bUVTableUsesTime) {
		  /*
		   ** Table in order of increasing time
		   */
		  zTime = dTimeIn;
		  for ( i=0; i < cl->nUV && zTime >= UV->zTime ; i++,UV++ );
		  }
	  else {
		  /*
		   ** Table in order of high to low redshift 
		   */
		  zTime = zIn;
		  for ( i=0; i < cl->nUV && zTime <= UV->zTime ; i++,UV++ );
		  }
	  }

  if (!cl->bUV || i==0) {
	  cl->R.Rate_Phot_HI = CL_RT_MIN;
	  cl->R.Rate_Phot_HeI = CL_RT_MIN;
	  cl->R.Rate_Phot_HeII = CL_RT_MIN;
	  
	  cl->R.Heat_Phot_HI = 0.0;
	  cl->R.Heat_Phot_HeI = 0.0;
	  cl->R.Heat_Phot_HeII = 0.0;
	  return;
	  }
  
  UV0=UV-1;
  if (i == cl->nUV ) {
	  cl->R.Rate_Phot_HI = UV0->Rate_Phot_HI;
	  cl->R.Rate_Phot_HeI = UV0->Rate_Phot_HeI;
	  cl->R.Rate_Phot_HeII = UV0->Rate_Phot_HeII;
	  
	  cl->R.Heat_Phot_HI = UV0->Heat_Phot_HI*CL_B_gm;
	  cl->R.Heat_Phot_HeI = UV0->Heat_Phot_HeI*CL_B_gm;
	  cl->R.Heat_Phot_HeII = UV0->Heat_Phot_HeII*CL_B_gm;
	  }
  else {
	  if (cl->bUVTableLinear) { /* use Linear interpolation */	
		  xx = (zTime - UV0->zTime)/(UV->zTime - UV0->zTime);
		  cl->R.Rate_Phot_HI = UV0->Rate_Phot_HI*(1-xx)+UV->Rate_Phot_HI*xx;
		  cl->R.Rate_Phot_HeI = UV0->Rate_Phot_HeI*(1-xx)+UV->Rate_Phot_HeI*xx;
		  cl->R.Rate_Phot_HeII = UV0->Rate_Phot_HeII*(1-xx)+UV->Rate_Phot_HeII*xx;
		  
		  cl->R.Heat_Phot_HI = (UV0->Heat_Phot_HI*(1-xx)+UV->Heat_Phot_HI*xx)*CL_B_gm;
		  cl->R.Heat_Phot_HeI = (UV0->Heat_Phot_HeI*(1-xx)+UV->Heat_Phot_HeI*xx)*CL_B_gm;
		  cl->R.Heat_Phot_HeII = (UV0->Heat_Phot_HeII*(1-xx)+UV->Heat_Phot_HeII*xx)*CL_B_gm;
		  }
	  else { /* use Log interpolation with 1+zTime */
		  xx = log((1+zTime)/(1+UV0->zTime))/log((1+UV->zTime)/(1+UV0->zTime));
		  cl->R.Rate_Phot_HI = pow(UV0->Rate_Phot_HI,1-xx)*pow(UV->Rate_Phot_HI,xx);
		  cl->R.Rate_Phot_HeI = pow(UV0->Rate_Phot_HeI,1-xx)*pow(UV->Rate_Phot_HeI,xx);
		  cl->R.Rate_Phot_HeII = pow(UV0->Rate_Phot_HeII,1-xx)*pow(UV->Rate_Phot_HeII,xx);
		  
		  cl->R.Heat_Phot_HI = pow(UV0->Heat_Phot_HI,1-xx)*pow(UV->Heat_Phot_HI,xx)*CL_B_gm;
		  cl->R.Heat_Phot_HeI = pow(UV0->Heat_Phot_HeI,1-xx)*pow(UV->Heat_Phot_HeI,xx)*CL_B_gm;
		  cl->R.Heat_Phot_HeII = pow(UV0->Heat_Phot_HeII,1-xx)*pow(UV->Heat_Phot_HeII,xx)*CL_B_gm;
		  }
	  }
  if (cl->R.Rate_Phot_HI < CL_RT_MIN) cl->R.Rate_Phot_HI = CL_RT_MIN;
  if (cl->R.Rate_Phot_HeI < CL_RT_MIN) cl->R.Rate_Phot_HeI = CL_RT_MIN;
  if (cl->R.Rate_Phot_HeII < CL_RT_MIN) cl->R.Rate_Phot_HeII = CL_RT_MIN;

  return;
  }

/* The following are self shielding correction tables based on the
   results of Pontzen, Governato et al 2008.  Based solely on local
   density, the HI, HeI and HeII photoionization rates from the UV
   background reduced by the factors in these tables.  The array of
   densities is unused, but documents the density dependence of the
   factor arrays.
 */

double AP_log_den_mp_percm3[] = { -10.25, -9.75, -9.25, -8.75, -8.25, -7.75, -7.25,
			    -6.75, -6.25, -5.75, -5.25, -4.75, -4.25, -3.75, 
			    -3.25, -2.75, -2.25,
			    -1.75, -1.25, -0.75, -0.25, 0.25, 0.75, 1.25, 1.75, 2.25 };
 
 
double AP_Gamma_HI_factor[] = { 0.99805271764596307, 0.99877911567687988, 0.99589340865612034,
			  0.99562060764857702, 0.99165170359332663, 0.9900889877822455,
			  0.98483276828954668, 0.97387675312245325, 0.97885673164000397,
			  0.98356305803821331, 0.96655786672182487, 0.9634906824933207,
			  0.95031917373653985, 0.87967606627349137, 0.79917533618355074,
			  0.61276011113763151, 0.16315185162187529, 0.02493663181368239,
			  0.0044013580765645335, 0.00024172553511936628, 1.9576102058649783e-10,
			  0.0, 0.0, 0.0, 0.0, 0.0 };

double AP_Gamma_HeI_factor[] = { 0.99284882980782224, 0.9946618686265097, 0.98641914356740497,
			   0.98867015777574848, 0.96519214493135597, 0.97188336387980656,
			   0.97529866247535113 , 0.97412477991428936, 0.97904139838765991,
			   0.98368372768570034, 0.96677432842215549, 0.96392622083382651,
			   0.95145730833093178, 0.88213871255482879, 0.80512823597731886,
			   0.62474472739578646, 0.17222786134467002, 0.025861959933038869,
			   0.0045265030237581529, 0.00024724339438128221, 1.3040144591221284e-08,
			   0.0, 0.0, 0.0, 0.0, 0.0};
 
 
double AP_Gamma_HeII_factor[] = { 0.97990208047216765, 0.98606251822654412,
			0.97657215632444849, 0.97274858503068629, 0.97416108746560681,
			0.97716929017896703, 0.97743607605974214, 0.97555305775319012,
			0.97874250764784809, 0.97849791914637996, 0.95135572977973504,
			0.92948461312852582, 0.89242272355549912, 0.79325512242742746 ,
			0.6683745597121028, 0.51605924897038324, 0.1840253816147828,
			0.035905775349044489, 0.0045537756654992923, 0.00035933897136804514,
			1.2294426136470751e-06, 0.0, 0.0, 0.0, 0.0, 0.0 };

void clRates( COOL *cl, RATE *Rate, double T, double rho ) {
  double Tln;
  double xTln,wTln0,wTln1;
  RATES_T *RT0,*RT1;
  int iTln;

  if (T >= cl->TMax) T=cl->TMax*(1.0 - EPS);   
  if (T < cl->TMin) T=cl->TMin;
  Tln = log(T);

  Rate->T = T;
  Rate->Tln = Tln; 

  xTln = (Tln-cl->TlnMin)*cl->rDeltaTln;
  iTln = xTln;
  assert(iTln >= 0);
  assert(iTln < cl->nTable - 1);
  RT0 = (cl->RT+iTln);
  RT1 = RT0+1; 
  wTln1 = xTln-iTln;
  wTln0 = 1-wTln1;
  
  Rate->Coll_HI = (wTln0*RT0->Rate_Coll_HI+wTln1*RT1->Rate_Coll_HI);
  Rate->Coll_HeI = (wTln0*RT0->Rate_Coll_HeI+wTln1*RT1->Rate_Coll_HeI);
  Rate->Coll_HeII = (wTln0*RT0->Rate_Coll_HeII+wTln1*RT1->Rate_Coll_HeII);

  Rate->Radr_HII = (wTln0*RT0->Rate_Radr_HII+wTln1*RT1->Rate_Radr_HII);
  Rate->Radr_HeII = (wTln0*RT0->Rate_Radr_HeII+wTln1*RT1->Rate_Radr_HeII);
  Rate->Diel_HeII = (wTln0*RT0->Rate_Diel_HeII+wTln1*RT1->Rate_Diel_HeII);
  Rate->Totr_HeII = Rate->Radr_HeII + Rate->Diel_HeII;
  Rate->Radr_HeIII = (wTln0*RT0->Rate_Radr_HeIII+wTln1*RT1->Rate_Radr_HeIII);

  Rate->Phot_HI = cl->R.Rate_Phot_HI;
  Rate->Phot_HeI = cl->R.Rate_Phot_HeI;
  Rate->Phot_HeII = cl->R.Rate_Phot_HeII;
  if (cl->bSelfShield) {
      /* Apply self shielding correction.  Note that density to array
         index calculation needs to match the AP_log_den_mp_percm3 array
         above. */
      double logen_B;
      logen_B = log10(rho*CL_B_gm);
      if (logen_B > 2.2499) {
	  Rate->Phot_HI = 0;
	  Rate->Phot_HeI = 0;
	  Rate->Phot_HeII = 0;
	  }
      else if (logen_B > -10.25) {
	  double x = (logen_B+10.25)*2.0;
	  int ix;
	  ix = floor(x);
	  x -= ix;
	  Rate->Phot_HI *= (AP_Gamma_HI_factor[ix]*(1-x)+AP_Gamma_HI_factor[ix+1]*x);
	  Rate->Phot_HeI *= (AP_Gamma_HeI_factor[ix]*(1-x)+AP_Gamma_HeI_factor[ix+1]*x);
	  Rate->Phot_HeII *= (AP_Gamma_HeII_factor[ix]*(1-x)+AP_Gamma_HeII_factor[ix+1]*x);
	  }
      }
}

/* Deprecated except for testing: use EdotInstant */
/* Need density in here to make this work with Self-Shielding */
double clHeatTotal ( COOL *cl, PERBARYON *Y, RATE *Rate ) {
  /* erg /gram /sec
     Note: QQ_* premultiplied by (CL_B_gm*erg_ev) */
    return Y->HI   * cl->R.Heat_Phot_HI * Rate->Phot_HI +
	Y->HeI  * cl->R.Heat_Phot_HeI * Rate->Phot_HeI +
	Y->HeII * cl->R.Heat_Phot_HeII * Rate->Phot_HeII;
}

/* Deprecated except for testing: use EdotInstant */
double clCoolTotal ( COOL *cl, PERBARYON *Y, RATE *Rate, double rho, double ZMetal ) {
  /* Assumes clRates called previously */
  /* erg /gram /sec */

  double en_B=rho*CL_B_gm;
  double xTln,wTln0,wTln1,LowTCool;
  RATES_T *RT0,*RT1;
  int iTln;

  xTln = (Rate->Tln-cl->TlnMin)*cl->rDeltaTln;
  iTln = xTln;
  assert(iTln >= 0);
  assert(iTln < cl->nTable - 1);
  RT0 = (cl->RT+iTln);
  RT1 = RT0+1; 
  wTln1 = xTln-iTln;
  wTln0 = 1-wTln1;

  if (Rate->T > cl->R.Tcmb)
      LowTCool = (wTln0*RT0->Cool_LowT+wTln1*RT1->Cool_LowT)*cl->R.Cool_LowTFactor*en_B*ZMetal;
  else
      LowTCool = 0;

  /* PUT INTO erg/gm/sec */
  return Y->e * ( 
#ifndef NOCOMPTON
    cl->R.Cool_Comp * ( Rate->T - cl->R.Tcmb ) + 
#endif
    en_B * (
    (wTln0*RT0->Cool_Brem_1+wTln1*RT1->Cool_Brem_1) * ( Y->HII + Y->HeII ) +
    (wTln0*RT0->Cool_Brem_2+wTln1*RT1->Cool_Brem_2) * Y->HeIII +

    (wTln0*RT0->Cool_Radr_HII+wTln1*RT1->Cool_Radr_HII) * Y->HII * Rate->Radr_HII +
    (wTln0*RT0->Cool_Radr_HeII+wTln1*RT1->Cool_Radr_HeII) * Y->HeII * Rate->Radr_HeII +
    (wTln0*RT0->Cool_Radr_HeIII+wTln1*RT1->Cool_Radr_HeIII) * Y->HeIII * Rate->Radr_HeIII +
 
    cl->R.Cool_Coll_HI * Y->HI * Rate->Coll_HI +
    cl->R.Cool_Coll_HeI * Y->HeI * Rate->Coll_HeI +
    cl->R.Cool_Coll_HeII * Y->HeII * Rate->Coll_HeII +

    cl->R.Cool_Diel_HeII * Y->HeII * Rate->Diel_HeII +

    (wTln0*RT0->Cool_Line_HI+wTln1*RT1->Cool_Line_HI) * Y->HI +
    (wTln0*RT0->Cool_Line_HeI+wTln1*RT1->Cool_Line_HeI) * Y->HeI +
    (wTln0*RT0->Cool_Line_HeII+wTln1*RT1->Cool_Line_HeII) * Y->HeII ) ) + LowTCool;
}

COOL_ERGPERSPERGM  clTestCool ( COOL *cl, PERBARYON *Y, RATE *Rate, double rho ) {
  /* Assumes clRates called previously */
  /* erg /gram /sec */

  double en_B=rho*CL_B_gm;
  double xTln,wTln0,wTln1;
  RATES_T *RT0,*RT1;
  int iTln;
  COOL_ERGPERSPERGM ret;

  xTln = (Rate->Tln-cl->TlnMin)*cl->rDeltaTln;
  iTln = xTln;
  assert(iTln >= 0);
  assert(iTln < cl->nTable - 1);
  RT0 = (cl->RT+iTln);
  RT1 = RT0+1; 
  wTln1 = xTln-iTln;
  wTln0 = 1-wTln1;

  /* PUT INTO erg/gm/sec */
  ret.compton = Y->e * ( 
    cl->R.Cool_Comp * ( Rate->T - cl->R.Tcmb ));
  ret.bremHII = Y->e * en_B * (
    (wTln0*RT0->Cool_Brem_1+wTln1*RT1->Cool_Brem_1) * ( Y->HII ));
  ret.bremHeII = Y->e * en_B * (
    (wTln0*RT0->Cool_Brem_1+wTln1*RT1->Cool_Brem_1) * ( Y->HeII ));
  ret.bremHeIII = Y->e * en_B * (
    (wTln0*RT0->Cool_Brem_2+wTln1*RT1->Cool_Brem_2) * Y->HeIII );
  ret.radrecHII = Y->e * en_B * 
    (wTln0*RT0->Cool_Radr_HII+wTln1*RT1->Cool_Radr_HII) * Y->HII * Rate->Radr_HII;
  ret.radrecHeII = Y->e * en_B * 
    (wTln0*RT0->Cool_Radr_HeII+wTln1*RT1->Cool_Radr_HeII) * Y->HeII * Rate->Radr_HeII;
  ret.radrecHeIII = Y->e * en_B * 
    (wTln0*RT0->Cool_Radr_HeIII+wTln1*RT1->Cool_Radr_HeIII) * Y->HeIII * Rate->Radr_HeIII;
  ret.collionHI = Y->e * en_B * 
    cl->R.Cool_Coll_HI * Y->HI * Rate->Coll_HI;
  ret.collionHeI = Y->e * en_B * 
    cl->R.Cool_Coll_HeI * Y->HeI * Rate->Coll_HeI;
  ret.collionHeII = Y->e * en_B * 
    cl->R.Cool_Coll_HeII * Y->HeII * Rate->Coll_HeII;
  ret.dielrecHeII = Y->e * en_B * 
    cl->R.Cool_Diel_HeII * Y->HeII * Rate->Diel_HeII;
  ret.lineHI = Y->e * en_B * 
    (wTln0*RT0->Cool_Line_HI+wTln1*RT1->Cool_Line_HI) * Y->HI;
  ret.lineHeI = Y->e * en_B * 
    (wTln0*RT0->Cool_Line_HeI+wTln1*RT1->Cool_Line_HeI) * Y->HeI;
  ret.lineHeII = Y->e * en_B * 
    (wTln0*RT0->Cool_Line_HeII+wTln1*RT1->Cool_Line_HeII) * Y->HeII;
  ret.lowT = en_B * 
    (wTln0*RT0->Cool_LowT+wTln1*RT1->Cool_LowT)*cl->R.Cool_LowTFactor*0.001; /* assume metallicity 0.001 */

  return ret;
}

void clPrintCool ( COOL *cl, PERBARYON *Y, RATE *Rate, double rho ) {
  /* Assumes clRates called previously */
  /* erg /gram /sec */

  double en_B=rho*CL_B_gm;
  double xTln,wTln0,wTln1;
  RATES_T *RT0,*RT1;
  int iTln;

  xTln = (Rate->Tln-cl->TlnMin)*cl->rDeltaTln;
  iTln = xTln;
  assert(iTln >= 0);
  assert(iTln < cl->nTable - 1);
  RT0 = (cl->RT+iTln);
  RT1 = RT0+1; 
  wTln1 = xTln-iTln;
  wTln0 = 1-wTln1;

  /* PUT INTO erg/gm/sec */
  printf("Compton:  %e\n",
    Y->e * ( 
    cl->R.Cool_Comp * ( Rate->T - cl->R.Tcmb )));
  printf("Cool Brem HII    %e\n",
    Y->e * en_B * (
    (wTln0*RT0->Cool_Brem_1+wTln1*RT1->Cool_Brem_1) * ( Y->HII )) );
  printf("Cool Brem HeII   %e\n",
    Y->e * en_B * (
    (wTln0*RT0->Cool_Brem_1+wTln1*RT1->Cool_Brem_1) * ( Y->HeII )) );
  printf("Cool Brem HeIII  %e\n",
    Y->e * en_B * (
    (wTln0*RT0->Cool_Brem_2+wTln1*RT1->Cool_Brem_2) * Y->HeIII ) );
  printf("Radiative Recombination  HII    %e\n",
    Y->e * en_B * 
    (wTln0*RT0->Cool_Radr_HII+wTln1*RT1->Cool_Radr_HII) * Y->HII * Rate->Radr_HII );
  printf("Radiative Recombination  HeII   %e\n",
    Y->e * en_B * 
    (wTln0*RT0->Cool_Radr_HeII+wTln1*RT1->Cool_Radr_HeII) * Y->HeII * Rate->Radr_HeII);
  printf("Radiative Recombination  HeIII  %e\n",
    Y->e * en_B * 
    (wTln0*RT0->Cool_Radr_HeIII+wTln1*RT1->Cool_Radr_HeIII) * Y->HeIII * Rate->Radr_HeIII);
  printf("Collisional Ionization  HI    %e\n",
    Y->e * en_B * 
    cl->R.Cool_Coll_HI * Y->HI * Rate->Coll_HI);
  printf("Collisional Ionization  HeI   %e\n",
    Y->e * en_B * 
    cl->R.Cool_Coll_HeI * Y->HeI * Rate->Coll_HeI);
  printf("Collisional Ionization  HeII  %e\n",
    Y->e * en_B * 
    cl->R.Cool_Coll_HeII * Y->HeII * Rate->Coll_HeII);
  printf("Dielectric Recombination HeII %e\n",
    Y->e * en_B * cl->R.Cool_Diel_HeII * Y->HeII * Rate->Diel_HeII);

  printf("Line cooling HI   %e\n",
    Y->e * en_B * 
    (wTln0*RT0->Cool_Line_HI+wTln1*RT1->Cool_Line_HI) * Y->HI);
  printf("Line cooling HeI  %e\n",
    Y->e * en_B * 
    (wTln0*RT0->Cool_Line_HeI+wTln1*RT1->Cool_Line_HeI) * Y->HeI);
  printf("Line cooling HeII  %e\n",
    Y->e * en_B * 
    (wTln0*RT0->Cool_Line_HeII+wTln1*RT1->Cool_Line_HeII) * Y->HeII );
  printf("Low T cooling (Z=0.001)  %e\n",
    en_B * 
    (wTln0*RT0->Cool_LowT+wTln1*RT1->Cool_LowT)*cl->R.Cool_LowTFactor*0.001);
}

void clPrintCoolFile( COOL *cl, PERBARYON *Y, RATE *Rate, double rho, FILE *fp ) {
  /* Assumes clRates called previously */
  /* erg /gram /sec */

  double en_B=rho*CL_B_gm;
  double xTln,wTln0,wTln1;
  RATES_T *RT0,*RT1;
  int iTln;

  xTln = (Rate->Tln-cl->TlnMin)*cl->rDeltaTln;
  iTln = xTln;
  assert(iTln >= 0);
  assert(iTln < cl->nTable - 1);
  RT0 = (cl->RT+iTln);
  RT1 = RT0+1; 
  wTln1 = xTln-iTln;
  wTln0 = 1-wTln1;

  /* PUT INTO erg/gm/sec */
  fprintf(fp,"%e ",
    Y->e * ( 
    cl->R.Cool_Comp * ( Rate->T - cl->R.Tcmb )));
  fprintf(fp,"%e ",
    Y->e * en_B * (
    (wTln0*RT0->Cool_Brem_1+wTln1*RT1->Cool_Brem_1) * ( Y->HII )) );
  fprintf(fp,"%e ",
    Y->e * en_B * (
    (wTln0*RT0->Cool_Brem_1+wTln1*RT1->Cool_Brem_1) * ( Y->HeII )) );
  fprintf(fp,"%e ",
    Y->e * en_B * (
    (wTln0*RT0->Cool_Brem_2+wTln1*RT1->Cool_Brem_2) * Y->HeIII ) );
  fprintf(fp,"%e ",
    Y->e * en_B * 
    (wTln0*RT0->Cool_Radr_HII+wTln1*RT1->Cool_Radr_HII) * Y->HII * Rate->Radr_HII );
  fprintf(fp,"%e ",
    Y->e * en_B * 
    (wTln0*RT0->Cool_Radr_HeII+wTln1*RT1->Cool_Radr_HeII) * Y->HeII * Rate->Radr_HeII);
  fprintf(fp,"%e ",
    Y->e * en_B * 
    (wTln0*RT0->Cool_Radr_HeIII+wTln1*RT1->Cool_Radr_HeIII) * Y->HeIII * Rate->Radr_HeIII);
  fprintf(fp,"%e ",
    Y->e * en_B * 
    cl->R.Cool_Coll_HI * Y->HI * Rate->Coll_HI);
  fprintf(fp,"%e ",
    Y->e * en_B * 
    cl->R.Cool_Coll_HeI * Y->HeI * Rate->Coll_HeI);
  fprintf(fp,"%e ",
    Y->e * en_B * 
    cl->R.Cool_Coll_HeII * Y->HeII * Rate->Coll_HeII);
  fprintf(fp,"%e ",
    Y->e * en_B * cl->R.Cool_Diel_HeII * Y->HeII * Rate->Diel_HeII);

  fprintf(fp,"%e ",
    Y->e * en_B * 
    (wTln0*RT0->Cool_Line_HI+wTln1*RT1->Cool_Line_HI) * Y->HI);
  fprintf(fp,"%e ",
    Y->e * en_B * 
    (wTln0*RT0->Cool_Line_HeI+wTln1*RT1->Cool_Line_HeI) * Y->HeI);
  fprintf(fp,"%e ",
    Y->e * en_B * 
    (wTln0*RT0->Cool_Line_HeII+wTln1*RT1->Cool_Line_HeII) * Y->HeII );
  fprintf(fp,"%e\n",
    en_B * 
    (wTln0*RT0->Cool_LowT+wTln1*RT1->Cool_LowT)*cl->R.Cool_LowTFactor*0.001);
}

void clAbunds( COOL *cl, PERBARYON *Y, RATE *Rate, double rho ) {
  double en_B=rho*CL_B_gm;
  double rcirrHI   = (Rate->Coll_HI)/(Rate->Radr_HII);
  double rcirrHeI  = (Rate->Coll_HeI)/(Rate->Totr_HeII);
  double rcirrHeII = (Rate->Coll_HeII)/(Rate->Radr_HeIII);
  double rpirrHI   = (Rate->Phot_HI)/(Rate->Radr_HII * en_B);
  double rpirrHeI  = (Rate->Phot_HeI)/(Rate->Totr_HeII * en_B);
  double rpirrHeII = (Rate->Phot_HeII)/(Rate->Radr_HeIII * en_B);
  double yHI = 0.;
  double yHeI = 0.;
  double yHeII = 0.;
  double yH = cl->Y_H;
  double yHe = cl->Y_He;
  double yeMAX = cl->Y_eMAX;
  double rye,ye;
  double fHI,fHeI,fHeII,rfHe,yHI_old,yHeII_old;
  int i;  

  for ( i=0 ; i<MAXABUNDITERATIONS ; i++ ) {
    yHI_old = yHI;
    yHeII_old = yHeII;

    ye = (yeMAX-(yHI + 2 * yHeI + yHeII));
    if (ye <= 0) {
      Y->e = 0;
      Y->HI = yH;
      Y->HII = 0;
      Y->HeI = yHe;
      Y->HeII = 0;
      Y->HeIII = 0;
      Y->Total = yH + yHe;
      return;
    }
    rye = 1/ye;
    fHI = rcirrHI + rpirrHI * rye;
    fHeI = rcirrHeI + rpirrHeI * rye;
    rfHe = 1 / ( 1 + fHeI * (1+rcirrHeII+rpirrHeII*rye) );

    yHI = yH / (1.0+fHI);
    yHeI = yHe * rfHe;
    yHeII = yHe * fHeI * rfHe;

    if ( fabs(yHeII_old-yHeII) < EPS * yHeII && fabs(yHI_old-yHI) < EPS * yHI ) break;
  }

  Y->e = ye;
  Y->HI = yHI;
  Y->HII = yH / (1.0/fHI+1.0);
  Y->HeI = yHeI;
  Y->HeII = yHeII;
  fHeII = rcirrHeII + rpirrHeII*rye;
  Y->HeIII = yHe / ((1.0/fHeI+1.0)/fHeII+1.0);
  Y->Total = Y->e + yH + yHe;
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

double clTemperaturePrimordial( COOL *cl, double Y_HI, double Y_HeI, double Y_HeII, double E ) {
	return clTemperature( 2*cl->Y_H - Y_HI + 3*cl->Y_He - 2*Y_HeI - Y_HeII,  E );
	}

/*-----------------------------------------------------------------
 *     Collisional Ionization rates
 *-----------------------------------------------------------------*/
/*     H + e- -> H+ + 2e-  Janev et al. 1987 (Abel 1996) */
double clRateCollHI( double T ) {
  double TL,arg;
    
  TL = log(T*CL_eV_per_K);
  arg = -32.713967867 + TL*(13.536556     + TL*(-5.73932875 +
      TL*(1.56315498 +
      TL*(-0.2877056     + TL*(3.48255977e-2 + TL*(-2.63197617e-3 +
      TL*(1.11954395e-4  + TL*(-2.03914985e-6))))))));
  if (arg < CL_MAX_NEG_EXP_ARG) return 0;
  return exp ( arg );
}
  
/*     He + e- -> He+ + 2e-  Janev et al. 1987 (Abel 1996) */
double clRateCollHeI( double T ) {
  double TL,arg;
    
  TL = log(T*CL_eV_per_K);
  arg = -44.09864886  + TL*(23.91596563   + TL*(-10.7532302 + 
      TL*(3.05803875 + 
      TL*(-0.56851189    + TL*(6.79539123e-2 + TL*(-5.00905610e-3 +
      TL*(2.06723616e-4  + TL*(-3.64916141e-6))))))));
  if (arg < CL_MAX_NEG_EXP_ARG) return 0;
  return exp( arg );
}

/*     He+ + e- -> He++ + 2e- Aladdin Database 1989 (Abel 1996) */
double clRateCollHeII( double T ) {
  double TL,arg;

  TL = log(T*CL_eV_per_K);
  arg = -68.71040990  + TL*(43.93347633   + TL*(-18.4806699 + 
      TL*(4.70162649 +
      TL*(-0.76924663    + TL*(8.113042e-2   + TL*(-5.32402063e-3 + 
      TL*(1.97570531e-4  + TL*(-3.16558106e-6))))))));
  if (arg < CL_MAX_NEG_EXP_ARG) return 0;
  return exp( arg );
}

/*-----------------------------------------------------------------
 *     Radiative Recombination rates
 *-----------------------------------------------------------------*/
/*     H+ + e- -> H + gam  Verner & Ferland 1996 */
double clRateRadrHII( double T ) {
  double Tsq = sqrt(T);

  return 7.982e-11/( Tsq*0.563615 *
		pow(1+Tsq*0.563615,0.252) * pow(1+Tsq*1.192167e-3,1.748));
}

/*     He+ + e- -> He + gam  radiative  Verner & Ferland 1996 */
double clRateRadrHeII( double T ) {
  /* 
   * Note that these functions do not meet perfectly at 1e6 -- 2% difference 
   * The derivatives are different there also: So the apparent error is large 
   */
  double Tsq = sqrt(T);
  if (T < 1e6)
    return  3.294e-11/( Tsq*0.253673 *
             pow(1+Tsq*0.253673,0.309) * pow(1+Tsq*1.649348e-4,1.691));
  else
    return  9.356e-10/( Tsq*4.841607 *
             pow(1+Tsq*4.841607,0.2108) * pow(1+Tsq*4.628935e-4,1.7892));
}

/*     He+ + e- -> He + gam  dielectronic  Aldovandi&Pequignot 1973 (Black 1981) */
double clRateDielHeII( double T ) {
  double T_inv = 1.0/T,arg;

  arg = -4.7e5*T_inv;
  if (arg < CL_MAX_NEG_EXP_ARG) return 0;
  return 1.9e-3*pow(T,-1.5)*exp(arg)*(1+0.3*exp(-9.4e4*T_inv));
}

/*     He++ + e- -> He+ + gam  Verner & Ferland 1996 */
double clRateRadrHeIII( double T ) {
  double Tsq = sqrt(T);

  return 1.891e-10/( Tsq*0.326686 *
          pow(1+Tsq*0.326686,0.2476) * pow(1+Tsq*6.004084e-4,1.7524));
}

/*-----------------------------------------------------------------
 *     Bremsstrahlung   
 *-----------------------------------------------------------------*/
#define CL_Cbremss1 1.426e-27
#define CL_al       0.79464
#define CL_bl       0.1243
#define CL_ar       2.13164
#define CL_br       (-0.1240)

double clCoolBrem1( double T ) {
  double Tlog10, Tsq;

  Tlog10 = log10(T);
  Tsq = sqrt(T);
  if (T < 3.2e5) 
     return Tsq*CL_Cbremss1*(CL_al+CL_bl*Tlog10)*CL_B_gm;
  else   
     return Tsq*CL_Cbremss1*(CL_ar+CL_br*Tlog10)*CL_B_gm;
}

#define CL_alog4   0.602059991
#define CL_alII    (4.0*(CL_al-CL_bl*CL_alog4))
#define CL_blII    (4.0*CL_bl)
#define CL_arII    (4.0*(CL_ar-CL_br*CL_alog4))
#define CL_brII    (4.0*CL_br)

double clCoolBrem2( double T ) {
  double Tlog10, Tsq;

  Tlog10 = log10(T);
  Tsq = sqrt(T);

  if (T<12.8e5) 
    return Tsq*CL_Cbremss1*(CL_alII+CL_blII*Tlog10)*CL_B_gm;
  else
    return Tsq*CL_Cbremss1*(CL_arII+CL_brII*Tlog10)*CL_B_gm;
}

/*-----------------------------------------------------------------
 *     Cooling multiplier for radiative recombination
 *-----------------------------------------------------------------*/
#define CL_aHII  0.0215964
#define CL_b     0.270251

double clCoolRadrHII( double T ) {

  double Tpow;
    
  Tpow=pow(T,CL_b);
  /* return CL_B_gm*(CL_eHI+exp(-CL_aHII*Tpow)*CL_k_Boltzmann*T); */
  /* Though 13.6eV is lost to the Gas as radiation, calculating the
   * Energy using u = 3/2 k T requires we don't subtract it here.
   */
  return CL_B_gm*(exp(-CL_aHII*Tpow)*CL_k_Boltzmann*T);
}
 
double clCoolRadrHeII( double T ) {

  double Tpow;
    
  Tpow=pow(T,CL_b);
  return CL_B_gm*(exp(-(CL_aHII*pow(13.6/24.59,CL_b))*Tpow)*CL_k_Boltzmann*T);
}

double clCoolRadrHeIII( double T ) {
  double Tpow;
    
  Tpow=pow(T,CL_b);
  return CL_B_gm*(exp(-(CL_aHII*pow(13.6/54.42,CL_b))*Tpow)*CL_k_Boltzmann*T);
}

/*-----------------------------------------------------------------
 *     Line Cooling
 *-----------------------------------------------------------------*/
/*      CEN (1992, Ap.J.Suppl 78,341) ADVOCATES MULTIPLYING EACH OF 
 *      THESE RATES BY Cen_correctn - HE CLAIMS THIS GIVES THE RIGHT
 *      HIGH T LIMIT FOR PROCESSES INVOLVING A FREE EL INTERACTING 
 *      WITH AN ORBITAL ELECTRON ?? */
#define CL_aHI  7.5e-19
#define CL_bHI  1.18348e05

double clCoolLineHI( double T ) {
  double T_inv, arg;
  double Cen_correctn = 1.0/(1.0+sqrt(T/1.0e5));

  T_inv=1.0/T;
  arg = -CL_bHI*T_inv;
  if (arg < CL_MAX_NEG_EXP_ARG) return 0;
  return CL_B_gm*CL_aHI*exp( arg )*Cen_correctn;
}

#define CL_aHeI   9.10e-27
#define CL_bHeI   1.3179e04
#define CL_p_HeI  0.1687

double clCoolLineHeI( double T ) {
  double T_inv,arg;
  double Cen_correctn = 1.0/(1.0+sqrt(T/1.0e5));

  T_inv=1.0/T;
  arg = -CL_bHeI*T_inv;
  if (arg < CL_MAX_NEG_EXP_ARG) return 0;
  return CL_B_gm*CL_aHeI*exp(-CL_bHeI*T_inv)*pow(T_inv,CL_p_HeI)*Cen_correctn;
}

#define CL_aHeII   5.54e-17
#define CL_bHeII   4.73638e05
#define CL_p_HeII  0.397

double clCoolLineHeII( double T ) {

  double T_inv,arg;
  double Cen_correctn = 1.0/(1.0+sqrt(T/1.0e5));

  T_inv=1.0/T;
  arg = -CL_bHeII*T_inv;
  if (arg < CL_MAX_NEG_EXP_ARG) return 0;
  return CL_B_gm*CL_aHeII*exp(-CL_bHeII*T_inv)*pow(T_inv,CL_p_HeII)*Cen_correctn;
}

double clCoolLowT( double T ) {
    double x;
    /* Cooling Rate for low T, fit from Bromm et al. MNRAS, 328, 969 (Figure 1). by Maschenko */
    /* Fit for metallicity Z = 0.001 -- scales linearly with Z */
    /* Returns cooling in erg cm^-3 s^-1 (after multiplied by n_H ^2) */
    /* Code uses erg g^-1 s^-1 so need to multiply the return value by Y_H^2 n_B * B_gm */
    
    if (T > 1e4 || T <= 10.001) return 0;
    x = log10(log10(log10(T)));
    return pow(10.0,-27.81 + 2.928*x - 0.6982*x*x);
    }

/* Returns Heating - Cooling excluding External Heating, units of ergs s^-1 g^-1 
   Public interface CoolEdotInstantCode */
double clEdotInstant( COOL *cl, PERBARYON *Y, RATE *Rate, double rho,
		    double ZMetal,  double *dEdotHeat, double *dEdotCool)
{
  double en_B = rho*CL_B_gm;
  double xTln,wTln0,wTln1;
  RATES_T *RT0,*RT1;
  int iTln;

  double ne,LowTCool;

  ne = Y->e*en_B;

  xTln = (Rate->Tln-cl->TlnMin)*cl->rDeltaTln;
  iTln = xTln;
  assert(iTln >= 0);
  assert(iTln < cl->nTable - 1);
  RT0 = (cl->RT+iTln);
  RT1 = RT0+1; 
  wTln1 = (xTln-iTln);
  wTln0 = (1-wTln1);

#define DTFRACLOWTCOOL 0.25
  if (Rate->T > cl->R.Tcmb*(1+DTFRACLOWTCOOL))
      LowTCool = (wTln0*RT0->Cool_LowT+wTln1*RT1->Cool_LowT)*cl->R.Cool_LowTFactor*en_B*ZMetal;
  else if (Rate->T < cl->R.Tcmb*(1-DTFRACLOWTCOOL))
      LowTCool = -(wTln0*RT0->Cool_LowT+wTln1*RT1->Cool_LowT)*cl->R.Cool_LowTFactor*en_B*ZMetal;
  else {
      double x = (Rate->T/cl->R.Tcmb-1)*(1./DTFRACLOWTCOOL);
      LowTCool = -(wTln0*RT0->Cool_LowT+wTln1*RT1->Cool_LowT)*cl->R.Cool_LowTFactor*en_B*ZMetal
	  *x*(3-3*fabs(x)+x*x);
      }

  *dEdotCool = 
#ifndef NOCOMPTON
      Y->e * cl->R.Cool_Comp * ( Rate->T - cl->R.Tcmb ) 
#endif
    + 
    ne * ((wTln0*RT0->Cool_Brem_1+wTln1*RT1->Cool_Brem_1) * ( Y->HII + Y->HeII ) +
	  (wTln0*RT0->Cool_Brem_2+wTln1*RT1->Cool_Brem_2) * Y->HeIII +
	  
	  cl->R.Cool_Diel_HeII * Y->HeII * Rate->Diel_HeII +

	  (wTln0*RT0->Cool_Line_HI+wTln1*RT1->Cool_Line_HI) * Y->HI +
	  (wTln0*RT0->Cool_Line_HeI+wTln1*RT1->Cool_Line_HeI) * Y->HeI +
	  (wTln0*RT0->Cool_Line_HeII+wTln1*RT1->Cool_Line_HeII) * Y->HeII +

	  (wTln0*RT0->Cool_Radr_HII+wTln1*RT1->Cool_Radr_HII) * Y->HII * Rate->Radr_HII  +
	  (wTln0*RT0->Cool_Radr_HeII+wTln1*RT1->Cool_Radr_HeII) * Y->HeII * Rate->Radr_HeII +
	  (wTln0*RT0->Cool_Radr_HeIII+wTln1*RT1->Cool_Radr_HeIII) * Y->HeIII * Rate->Radr_HeIII +
	    
	  cl->R.Cool_Coll_HI * Y->HI * Rate->Coll_HI +
	  cl->R.Cool_Coll_HeI * Y->HeI * Rate->Coll_HeI + 
	  cl->R.Cool_Coll_HeII * Y->HeII * Rate->Coll_HeII )
    + 
      LowTCool;
  *dEdotHeat = 
    Y->HI   * cl->R.Heat_Phot_HI * Rate->Phot_HI +
    Y->HeI  * cl->R.Heat_Phot_HeI * Rate->Phot_HeI +
    Y->HeII * cl->R.Heat_Phot_HeII * Rate->Phot_HeII;
  return *dEdotHeat - *dEdotCool;
}

/*
 *  We solve the implicit equation:  
 *  Eout = Ein + dt * Cooling ( Yin, Yout, Ein, Eout )
 *
 *  E erg/gm, ExternalHeating erg/gm/s, rho gm/cc, dt sec
 * 
 */

double clfTemp(void *Data, double T) 
{
  clDerivsData *d = Data; 
  
  d->its++;
  clRates( d->cl, &d->Rate, T, d->rho );
  clAbunds( d->cl, &d->Y, &d->Rate, d->rho );

  return d->E-clThermalEnergy( d->Y.Total, T );
}

void clTempIteration( clDerivsData *d )
{
 double T,TA,TB;
 d->its = 0;

 if (d->E <= 0) T=d->cl->TMin;
 else {

   TB = clTemperature( d->Y_Total0, d->E );
   TA = clTemperature( d->Y_Total1, d->E );
   if (TA > TB) { T=TA; TA=TB; TB=T; }

   T = RootFind(clfTemp, d, TA, TB, EPSTEMP*TA ); 
 } 
 d->its++;
 clRates( d->cl, &d->Rate, T, d->rho );
 clAbunds( d->cl, &d->Y, &d->Rate, d->rho );
}

void clDerivs(double x, const double *y, double *dHeat, double *dCool,
	     void *Data) {
  clDerivsData *d = Data;

  d->E = y[0];
  clTempIteration( d );
  clEdotInstant( d->cl, &d->Y, &d->Rate, d->rho, d->ZMetal, dHeat, dCool );
  if(d->ExternalHeating > 0.0)
      *dHeat += d->ExternalHeating;
  else
      *dCool -= d->ExternalHeating;
}

#if 0
/* Unused */
int clJacobn(double x, const double y[], double dfdx[], double *dfdy, void *Data) {
  clDerivsData *d = Data;
  double E = y[0],dE;
  const int ndim = 1;

  dfdx[0] = 0;

  /* Approximate dEdt/dE */
  d->E = E*(EMUL);
  clTempIteration( d );
  dE = clEdotInstant( d->cl, &d->Y, &d->Rate, d->rho, d->ZMetal );

  d->E = E*(1/EMUL);
  clTempIteration( d );
  dE -= clEdotInstant( d->cl, &d->Y, &d->Rate, d->rho, d->ZMetal );

  dfdy[0*ndim + 0] = dE/(E*d->dlnE);
  return GSL_SUCCESS;
}
#endif

void clIntegrateEnergy(COOL *cl, clDerivsData *clData, PERBARYON *Y, double *E,
		       double ExternalHeating, double rho, double ZMetal, double tStep ) {

  double EMin;
  double t=0;
  clDerivsData *d = clData;
  STIFF *sbs = d->IntegratorContext;
  
  if (tStep == 0) return;

  d->rho = rho;
  d->ExternalHeating = ExternalHeating;
  d->ZMetal = ZMetal;
  
  clRates( d->cl, &d->Rate, cl->TMin, d->rho );
  clAbunds( d->cl, &d->Y, &d->Rate, d->rho );

  EMin = clThermalEnergy( d->Y.Total, cl->TMin );
   if (tStep < 0)  {   /* Don't integrate energy equation */
      tStep = fabs(tStep);
      *E += ExternalHeating*tStep;
      if (*E < EMin) *E=EMin;
      d->E = *E;
      clTempIteration( d );
      *Y = d->Y;
      return;
      };  

  {
      StiffStep( sbs, E, t, tStep);
#ifdef ASSERTENEG      
      assert(*E > 0.0);
#else
      if (*E < EMin) {
	*E = EMin;
      }
#endif    
      cl->its = 1;
    }
  /* 
     Note Stored abundances are not necessary with
     this eqm integrator therefore the following operations
     could be moved to output routines 
  */

  d->E = *E;
  clTempIteration( d );
  *Y = d->Y;
}

/* Module Interface routines */

void CoolAddParams( COOLPARAM *CoolParam, PRM prm ) {
	CoolParam->bIonNonEqm = 0;
	prmAddParam(prm,"bIonNonEqm",paramBool,&CoolParam->bIonNonEqm,
				sizeof(int),"IonNonEqm",
				"<Gas is Cooling Non-Equilibrium Abundances> = +IonNonEqm");
	CoolParam->bUV = 1;
	prmAddParam(prm,"bUV",paramBool,&CoolParam->bUV,sizeof(int),"UV",
				"read in an Ultra Violet file = +UV");
	CoolParam->bUVTableUsesTime = 0;
	prmAddParam(prm,"bUVTableUsesTime",paramBool,&CoolParam->bUVTableUsesTime,sizeof(int),"UVTableUsesTime",
				"Ultra Violet Table uses time = +UVTableUsesTime");
	CoolParam->dMassFracHelium = 0.25;
	prmAddParam(prm,"dMassFracHelium",paramDouble,&CoolParam->dMassFracHelium,
				sizeof(double),"hmf",
				"<Primordial Helium Fraction (by mass)> = 0.25");
	CoolParam->dCoolingTmin = 10;
	prmAddParam(prm,"dCoolingTmin",paramDouble,&CoolParam->dCoolingTmin,
				sizeof(double),"ctmin",
				"<Minimum Temperature for Cooling> = 10K");
	CoolParam->dCoolingTmax = 1e9;
	prmAddParam(prm,"dCoolingTmax",paramDouble,&CoolParam->dCoolingTmax,
				sizeof(double),"ctmax",
				"<Maximum Temperature for Cooling> = 1e9K");
	CoolParam->nCoolingTable = 15001;
	prmAddParam(prm,"nCoolingTable",paramInt,&CoolParam->nCoolingTable,
				sizeof(int),"nctable","<# Cooling table elements> = 15001");
	CoolParam->bDoIonOutput = 1;
	prmAddParam(prm,"bDoIonOutput",paramBool,&CoolParam->bDoIonOutput,sizeof(int),
				"Iout","enable/disable Ion outputs (cooling only) = +Iout");
	CoolParam->bLowTCool = 0;
	prmAddParam(prm,"bLowTCool",paramBool,&CoolParam->bLowTCool,sizeof(int),
				"ltc","enable/disable low T cooling = +ltc");
	CoolParam->bSelfShield = 0;
	prmAddParam(prm,"bSelfShield",paramBool,&CoolParam->bSelfShield,sizeof(int),
				"ssc","enable/disable Self Shielded Cooling = +ssc");
	}
	
void CoolOutputArray( COOLPARAM *CoolParam, int cnt, int *type, char *suffix ) {
	}

/* Output Conversion Routines */

double CoolEnergyToTemperature( COOL *Cool, COOLPARTICLE *cp, double E, double fMetal ) {
	return clTemperature(2*(Cool)->Y_H - cp->Y_HI + 
						 3*(Cool)->Y_He - 2*cp->Y_HeI - cp->Y_HeII, E );
	}

double CoolCodeEnergyToTemperature( COOL *Cool, COOLPARTICLE *cp, double E, double fMetal ) {
	return CoolEnergyToTemperature( Cool, cp, E*Cool->dErgPerGmUnit, fMetal );
	}

/* Initialization Routines */

void CoolTableReadInfo( COOLPARAM *CoolParam, int cntTable, int *nTableColumns, char *suffix )
{
   int localcntTable = 0;

   *nTableColumns = 0;

   if (CoolParam->bUV) {
	   if (localcntTable == cntTable) {
		   *nTableColumns = 7;
		   sprintf(suffix,"UV");
		   return;
		   }
	   localcntTable++;
	   }
   }

void CoolTableRead( COOL *Cool, int nData, void *vData)
{
   double *dTableData = vData;
   int nTableColumns; 
   int nTableRows;
   int localcntTable = 0;

   if (Cool->bUV) {
	   if (localcntTable == Cool->nTableRead) {
		   nTableColumns = 7;
		   nTableRows = nData/(sizeof(double)*nTableColumns);
		   assert( nData == sizeof(double)*nTableColumns*nTableRows );
		   clInitUV( Cool, nTableColumns, nTableRows, dTableData );
		   Cool->nTableRead++;
		   return;
		   }
	   localcntTable++;
	   }
   
   fprintf(stderr," Attempt to initialize non-exitent table in cooling\n");
   assert(0);

   }

void CoolDefaultParticleData( COOLPARTICLE *cp )
{
	cp->Y_HI = 0.75;
	cp->Y_HeI = 0.0625;
	cp->Y_HeII = 0.0;
	}

void CoolInitEnergyAndParticleData( COOL *cl, COOLPARTICLE *cp, double *E, double dDensity, double dTemp, double fMetal )
{
	PERBARYON Y;
    RATE r;

	cp->Y_HI = cl->Y_H;
	cp->Y_HeI = cl->Y_He;
	cp->Y_HeII = 0.0;
	CoolPARTICLEtoPERBARYON(cl, &Y, cp);
	clRates(cl,&r,dTemp,CodeDensityToComovingGmPerCc(cl,dDensity));
	clAbunds(cl,&Y,&r,CodeDensityToComovingGmPerCc(cl,dDensity));
	CoolPERBARYONtoPARTICLE(cl, &Y,cp);
	*E = clThermalEnergy(Y.Total,dTemp)*cl->diErgPerGmUnit;
	}

void CoolInitRatesTable( COOL *cl, COOLPARAM CoolParam ) {
	clInitRatesTable( cl, CoolParam.dCoolingTmin, CoolParam.dCoolingTmax, CoolParam.nCoolingTable );
	}

void CoolSetTime( COOL *cl, double dTime, double z ) {
	clRatesRedshift( cl, z, dTime );
	}

/* Integration Routines */

/**
 * Physical units
 */
void CoolIntegrateEnergy(COOL *cl, clDerivsData *clData, COOLPARTICLE *cp,
			 double *E,
                       double PdV, double rho, double ZMetal, double tStep ) {
        PERBARYON Y;

        CoolPARTICLEtoPERBARYON(cl, &Y,cp);
        clIntegrateEnergy(cl, clData, &Y, E, PdV, rho, ZMetal, tStep);
        CoolPERBARYONtoPARTICLE(cl, &Y, cp);
        }

/**
 * Code Units except for dt
 */
void CoolIntegrateEnergyCode(COOL *cl, clDerivsData *clData, COOLPARTICLE *cp,
			     double *ECode, double ExternalHeatingCode,
			     double rhoCode, double ZMetal, double *posCode,
			     double tStep ) {
	PERBARYON Y;
	double E;
	
	E = CoolCodeEnergyToErgPerGm( cl, *ECode );
	CoolPARTICLEtoPERBARYON(cl, &Y, cp);
	clIntegrateEnergy(cl, clData, &Y, &E,
			  CoolCodeWorkToErgPerGmPerSec( cl, ExternalHeatingCode ), 
			  CodeDensityToComovingGmPerCc(cl, rhoCode), ZMetal, tStep);
	CoolPERBARYONtoPARTICLE(cl, &Y, cp);
	*ECode = CoolErgPerGmToCodeEnergy(cl, E);

	}

/* Deprecated: */
double CoolHeatingRate( COOL *cl, COOLPARTICLE *cp, double T, double dDensity, double ZMetal ) {
    PERBARYON Y;
    RATE Rate;

    CoolPARTICLEtoPERBARYON(cl, &Y, cp);
    clRates(cl, &Rate, T, dDensity);

    return (-clCoolTotal(cl, &Y, &Rate, dDensity, ZMetal ) + clHeatTotal(cl, &Y, &Rate));
    }

/* Code heating - cooling rate excluding external heating (PdV, etc..) */
double CoolEdotInstantCode(COOL *cl, COOLPARTICLE *cp, double ECode, 
			  double rhoCode, double ZMetal, double *posCode ) {
    PERBARYON Y;
    RATE Rate;
    double T,E,rho,Edot;
    double dHeat, dCool;

    E = CoolCodeEnergyToErgPerGm( cl, ECode );
    T = CoolEnergyToTemperature( cl, cp, E, ZMetal );
    rho = CodeDensityToComovingGmPerCc(cl,rhoCode );
    CoolPARTICLEtoPERBARYON(cl, &Y, cp);
    clRates(cl, &Rate, T, rho);
    
    Edot = clEdotInstant( cl, &Y, &Rate, rho, ZMetal, &dHeat, &dCool );

    return CoolErgPerGmPerSecToCodeWork( cl, Edot );
    }

/* Code heating - cooling rate excluding external heating (PdV, etc..) */
double CoolCoolingCode(COOL *cl, COOLPARTICLE *cp, double ECode, 
			  double rhoCode, double ZMetal, double *posCode ) {
    PERBARYON Y;
    RATE Rate;
    double T,E,rho,Edot;

    E = CoolCodeEnergyToErgPerGm( cl, ECode );
    T = CoolEnergyToTemperature( cl, cp, E, ZMetal );
    rho = CodeDensityToComovingGmPerCc(cl,rhoCode );
    CoolPARTICLEtoPERBARYON(cl, &Y, cp);
    clRates(cl, &Rate, T, rho);
    
    Edot = clCoolTotal(cl, &Y, &Rate, rho, ZMetal );

    return CoolErgPerGmPerSecToCodeWork( cl, Edot );
    }

/* Code heating due to atomic/radiative processes only */
double CoolHeatingCode(COOL *cl, COOLPARTICLE *cp, double ECode, 
			  double rhoCode, double ZMetal, double *posCode ) {
    PERBARYON Y;
    RATE Rate;
    double T,E,rho,Edot;

    E = CoolCodeEnergyToErgPerGm( cl, ECode );
    T = CoolEnergyToTemperature( cl, cp, E, ZMetal );
    rho = CodeDensityToComovingGmPerCc(cl,rhoCode );
    CoolPARTICLEtoPERBARYON(cl, &Y, cp);
    clRates(cl, &Rate, T, rho);
    
    Edot = clHeatTotal ( cl, &Y, &Rate );

    return CoolErgPerGmPerSecToCodeWork( cl, Edot );
    }

#ifdef TEST_COOLING

int main(int argc, char **argv)
{
    COOL *cl;
    double rate;
    PERBARYON Y;
    double rho = 1.0; /* 1 gm/cc */
    double ZMetal = 0.0;
    double lgT, dlgT;
    clDerivsData *d;
    double dCoolingTmin = 10.0;
    double dCoolingTmax = 1e9;
    double nCoolingTable = 10000;
    
    cl = CoolInit();
    cl->dGmPerCcUnit = 1.0;
    cl->dComovingGmPerCcUnit = 1.0;
    cl->dErgPerGmUnit = 1.0;
    cl->dSecUnit = 1.0;
    cl->dErgPerGmPerSecUnit = cl->dErgPerGmUnit / cl->dSecUnit;
    cl->diErgPerGmUnit = 1./cl->dErgPerGmUnit;
    cl->dKpcUnit = 0.0;
  
    cl->Y_H = 1.0-0.24;
    cl->Y_He = 0.24/4;
    cl->Y_eMAX = cl->Y_H + cl->Y_He*2;

    cl->bUV = 0;
    cl->bLowTCool = 0;
    cl->bSelfShield = 0;
    clInitRatesTable( cl, dCoolingTmin, dCoolingTmax, nCoolingTable );
    d = CoolDerivsInit(cl);
    
    for(lgT = 3.4; lgT < 8.0; lgT += .05) {
	double T = exp(lgT*log(10.0));
	d->rho = 1.67e-24;
	d->E = clfTemp(d, T);
	clRates( d->cl, &d->Rate, T, d->rho );
	clAbunds( d->cl, &d->Y, &d->Rate, d->rho );
 	d->E=clThermalEnergy( d->Y.Total, T );
	rate = clEdotInstant( cl, &d->Y, &d->Rate, rho, ZMetal );
	printf("%g %g %g %g\n", T, d->rho*d->rho*rate, d->rho*d->E, d->E/rate);
	}
    return 0;
    }

#endif
#endif /* NOCOOLING */

