/*#define ASSERTENEG*/

#ifndef NOCOOLING

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
/* #include <rpc/xdr.h> */

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

/* This min is to prevent Y_e = 0 which shuts off all collisional processes (cooling, recomb or ionization)
   and prevents unphysical negative Y_e (resulting from roundoff)
   This is important for cases with no cosmic rays or UV 
   */
#define Y_EMIN 1e-7

#define TESTRATE 1e-3
/* #define CUBICTABLEINTERP   */
/* #define USETABLE */

#ifdef USETABLE
#define CLRATES( _cl, _Rate, _T, _rho, _ZMetal, _columnL, _Rate_Phot_H2_stellar) clRates_Table( _cl, _Rate, _T, _rho, _ZMetal, _columnL, _Rate_Phot_H2_stellar)
#define CLEDOTINSTANT( _cl, _Y, _Rate, _rho, _ZMetal, _Heat, _Cool ) clEdotInstant_Table( _cl, _Y, _Rate, _rho, _ZMetal, _Heat, _Cool)

#else
#define CLRATES( _cl, _Rate, _T, _rho, _ZMetal, _columnL, _Rate_Phot_H2_stellar)   clRates( _cl, _Rate, _T, _rho, _ZMetal, _columnL, _Rate_Phot_H2_stellar)
#define CLEDOTINSTANT( _cl, _Y, _Rate, _rho, _ZMetal, _Heat, _Cool ) clEdotInstant( _cl, _Y, _Rate, _rho, _ZMetal, _Heat, _Cool)

#endif

#ifdef CUBICTABLEINTERP
#define TABLEFACTOR 2
#else 
#define TABLEFACTOR 1
#endif

#ifdef CUBICTABLEINTERP
#define TABLEINTERP( _rname ) (wTln0*RT0->_rname+wTln1*RT1->_rname+wTln0d*RT0d->_rname+wTln1d*RT1d->_rname)
#else
#define TABLEINTERP( _rname ) (wTln0*RT0->_rname+wTln1*RT1->_rname)
#endif

#include "cooling.h"

/* Integrator constants */

/* Accuracy target for integrators */
#define EPSINTEG  1e-3
#define WARNINTEGITS 100
#define MAXINTEGITS 10000

/* Debugging Information */
#define PARTICLEIORD    100

#define CL_eH2     (4.476*CL_eV_erg) /*Energy lost during H2 dissociation, Shapira & Kang (1987) through Abel 1997, page 194 */
#define CL_eHI     (13.60*CL_eV_erg)
#define CL_eHeI    (24.59*CL_eV_erg)
#define CL_eHeII   (54.42*CL_eV_erg)
#define CL_E2HeII  (3.0*13.6*CL_eV_erg)

#define EMUL (1.000001)

#define M_H      1.672e-24
#define SUM_Metal 0.0184707   
COOL *CoolInit( )
{
  COOL *cl;
  cl = (COOL *) malloc(sizeof(COOL));
  assert(cl!=NULL);
  
  cl->nUV = 0;
  cl->UV = NULL; 

  cl->nTable = 0;
  cl->RT = NULL;
  cl->MetalCoolln = NULL; 
  cl->MetalHeatln = NULL; 

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

    assert(cl != NULL);
    Data = malloc(sizeof(clDerivsData));
    assert(Data != NULL);
    Data->IntegratorContext = StiffInit(EPSINTEG, 5, Data, clDerivs); /*Change array length to 5 to include H2*/
    Data->cl = cl;

    return Data;
    }

void CoolFinalize(COOL *cl ) 
{
  if (cl->UV != NULL) free(cl->UV);
  if (cl->RT != NULL) free(cl->RT);
  if (cl->MetalCoolln != NULL) free(cl->MetalCoolln); 
  if (cl->MetalHeatln != NULL) free (cl->MetalHeatln); 
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
  cl->dMassFracHelium = CoolParam.dMassFracHelium;

  cl->bUV = CoolParam.bUV;
  cl->bUVTableUsesTime = CoolParam.bUVTableUsesTime;
  cl->bUVTableLinear = CoolParam.bUVTableUsesTime; /* Linear if using time */
  cl->bLowTCool = CoolParam.bLowTCool;
  cl->bSelfShield = CoolParam.bSelfShield;
  cl->bMetal = CoolParam.bMetal; 
  cl->bShieldHI = CoolParam.bShieldHI;
  cl->dClump = CoolParam.dClump;
  cl->dLymanWernerFrac = CoolParam.dLymanWernerFrac;

}

/* caculate the number fraction YH, YHe as a function of metalicity. Cosmic 
   production rate of helium relative to metals (in mass)  
   delta Y/delta Z = 2.1 and primordial He Yp = 0.236 (Jimenez et al. 2003,Science 
   299, 5612. 
   piecewise linear
   Y = Yp + dY/dZ*ZMetal up to ZMetal = 0.1, then linear decrease to 0 at Z=1)  

  SUM_Metal = sum(ni/nH *mi),it is a fixed number for cloudy abundance. 
  Massfraction fmetal = Z*SUM_metal/(1 + 4*nHe/nH + Z*SUM_metal) (1)
  4*nHe/nH = mHe/mH = fHe/fH 
  also fH + fHe + fMetal = 1  (2)
  if fHe specified, combining the 2 eq above will solve for 
  fH and fMetal 
  
  
  a = CoolParam.dMassFracHelium; 
  c = SUM_Metal*ZMetal;
  fH = (-(2.0*a + a*c -1) + sqrt(pow((2.*a+a*c-1),2.)-4.0*a*(1+c)*(a-1)))/(2.*(1+c)); 
  fMetal = 1.0 - fH - CoolParam.dMassFracHelium; 

 */ 

void clSetAbundanceTotals(COOL *cl, double ZMetal, double *pY_H, double *pY_He, double *pY_eMax) {
    double Y_H, Y_He;
    
    if (ZMetal <= 0.1) {
        Y_He = (0.236 + 2.1*ZMetal)/4.0;
	}
    else {
        Y_He = (-0.446*(ZMetal - 0.1)/0.9 + 0.446)/4.0;
	}
    Y_H = 1.0 - Y_He*4.0 - ZMetal; 

    *pY_H = Y_H;
    *pY_He = Y_He;
    *pY_eMax = Y_H+ Y_He*2; /* Ignoring any electrons from metals */

}

void CoolPARTICLEtoPERBARYON(COOL *cl, PERBARYON *Y, COOLPARTICLE *cp, double ZMetal) {
    double Y_H, Y_He, Y_eMax;
    clSetAbundanceTotals(cl,ZMetal,&Y_H,&Y_He,&Y_eMax);
    Y->HI = cp->f_HI*Y_H;
    Y->H2 = cp->f_H2/2.0*Y_H; /*Number of H2 molecules total from fraction of H atoms in H2*/
    Y->HII = Y_H - Y->HI - 2.0*Y->H2;
    Y->HeI = cp->f_HeI*Y_He;
    Y->HeII = cp->f_HeII*Y_He;
    Y->HeIII = Y_He - Y->HeI - Y->HeII;
    Y->e = Y->HII + Y->HeII + 2*Y->HeIII;
    Y->Total = Y->e + Y_H + Y_He + ZMetal/MU_METAL; 
}

void CoolPERBARYONtoPARTICLE(COOL *cl, PERBARYON *Y, COOLPARTICLE *cp, double ZMetal) {
    double Y_H, Y_He, Y_eMax;
    clSetAbundanceTotals(cl,ZMetal,&Y_H,&Y_He,&Y_eMax);
    cp->f_HI = Y->HI/Y_H;
    cp->f_H2 = 2.0*(Y->H2)/Y_H; /*Fraction of H atoms that are in H2 from number of H2 molecules CC*/
    cp->f_HeI = Y->HeI/Y_He;
    cp->f_HeII = Y->HeII/Y_He;
}


/*Returns baryonic fraction for a given species*/
double COOL_ARRAY0(COOL *cl, COOLPARTICLE *cp, double ZMetal) {
    double Y_H, Y_He, Y_eMax;
    clSetAbundanceTotals(cl,ZMetal,&Y_H,&Y_He,&Y_eMax);
    return (cp->f_HI*Y_H);
}

double COOL_ARRAY1(COOL *cl, COOLPARTICLE *cp, double ZMetal) {
    double Y_H, Y_He, Y_eMax;
    clSetAbundanceTotals(cl,ZMetal,&Y_H,&Y_He,&Y_eMax);
    return (cp->f_HeI*Y_He);
}

double COOL_ARRAY2(COOL *cl, COOLPARTICLE *cp, double ZMetal) {
    double Y_H, Y_He, Y_eMax;
    clSetAbundanceTotals(cl,ZMetal,&Y_H,&Y_He,&Y_eMax);
    return (cp->f_HeII*Y_He);
}

double COOL_ARRAY3(COOL *cl, COOLPARTICLE *cp, double ZMetal) {
    double Y_H, Y_He, Y_eMax;
    clSetAbundanceTotals(cl,ZMetal,&Y_H,&Y_He,&Y_eMax);
    return (cp->f_H2/2.0*Y_H);
}

void clReadMetalTable(COOL *cl, COOLPARAM clParam)
{
  int i,j,k,l, nt, nnH, nz; 
  double dt, dnH, dz, tminlog, tmaxlog, nHminlog, nHmaxlog,zmin,zmax; 
  FILE *fp;  
  XDR xdrs; 
  
  /* fp = fopen(clParam.CoolInFile, "r"); */
  fp = fopen("cooltable_xdr", "r"); 
  assert(fp != NULL); 
  fscanf(fp, "%d %lf %lf %lf \n",&nz, &zmin, &zmax, &dz);
  fscanf(fp, "%d %lf %lf %lf \n",&nnH, &nHminlog ,&nHmaxlog,&dnH);
  fscanf(fp, "%d %lf %lf %lf \n",&nt, &tminlog, &tmaxlog, &dt);

  cl->bMetal = clParam.bMetal; 
  cl->nzMetalTable = nz;
  cl->nnHMetalTable = nnH;
  cl->nTMetalTable = nt; 
  
  cl->MetalTMin = pow(10.0, tminlog); 
  cl->MetalTMax = pow(10.0, tmaxlog);
  cl->MetalnHMin = pow(10.0, nHminlog);   /* Note: here nH means mass density */ 
  cl->MetalnHMax = pow(10.0, nHmaxlog); 
  cl->MetalzMin = zmin;
  cl->MetalzMax = zmax;
  
  cl->MetalTlogMin = tminlog;
  cl->MetalTlogMax = tmaxlog;
  cl->MetalnHlogMin = nHminlog; 
  cl->MetalnHlogMax = nHmaxlog;

  cl->rDeltaTlog = 1./dt; 
  cl->rDeltanHlog = 1./dnH; 
  cl->rDeltaz = 1./dz;  
  
  cl->nUV = nz -1;
  cl->UV = (UVSPECTRUM *)malloc(sizeof(UVSPECTRUM)*(cl->nUV));
  
  /* 3D table, f = f(T, rho, z), in redshift the order is from high z-->low z, 
     density and temperature are in the order low --> high */
 
  cl->MetalCoolln = (float ***)malloc(nz*sizeof(float ***)); 
  for (i=0; i<=nz-1;i++){
    cl->MetalCoolln[i] = (float **) malloc(nnH*sizeof(float **));
  }
  
  for (i=0; i<=nz-1;i++){
    for(j=0; j<=nnH-1; j++){
     cl->MetalCoolln[i][j] = (float *) malloc(nt*sizeof(float *));
    }
  }

  cl->MetalHeatln = (float ***)malloc(nz*sizeof(float ***)); 
  for (i=0; i<=nz-1;i++){
    cl->MetalHeatln[i] = (float **) malloc(nnH*sizeof(float **));
  }
  
  for (i=0; i<=nz-1;i++){
    for(j=0; j<=nnH-1; j++){
     cl->MetalHeatln[i][j] = (float *) malloc(nt*sizeof(float *));
    }
  }

  /*  for(;;){
    fscanf(fp, "%s", comments);
    if (comments[0] eq "$") break;  
    } */
  xdrstdio_create(&xdrs, fp, XDR_DECODE);
  for(j=0; j<=nnH-1; j++){
    for(k=0; k<=nt-1; k++){
      xdr_float(&xdrs, &cl->MetalCoolln[0][j][k]); /*redshift: i, density: j, temperature: k */
      xdr_float(&xdrs, &cl->MetalHeatln[0][j][k]);
    }
   }     
  
  for(i=nz-2, l=0; i>=0; i--,l++){
    xdr_double(&xdrs, &cl->UV[l].zTime); 
    xdr_double(&xdrs, &cl->UV[l].Rate_Phot_HI);
    xdr_double(&xdrs, &cl->UV[l].Rate_Phot_HeI);
    xdr_double(&xdrs, &cl->UV[l].Rate_Phot_HeII);

    xdr_double(&xdrs, &cl->UV[l].Heat_Phot_HI);
    xdr_double(&xdrs, &cl->UV[l].Heat_Phot_HeI);
    xdr_double(&xdrs, &cl->UV[l].Heat_Phot_HeII);

    for(j=0; j<=nnH-1; j++){
      for(k=0; k<=nt-1; k++){
	xdr_float(&xdrs, &cl->MetalCoolln[l+1][j][k]);
	xdr_float(&xdrs, &cl->MetalHeatln[l+1][j][k]);
      }
    } 
  }
 
  fclose(fp);
  return;  
}


void clInitUV(COOL *cl, int nTableColumns, int nTableRows, double *dTableData )
{
    int i;
    printf("call clInit UV...\n"); 
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
		(cl->UV)[i].Rate_Phot_H2_cosmo = dTableData[i*nTableColumns+4];

		(cl->UV)[i].Heat_Phot_HI = dTableData[i*nTableColumns+5];
		(cl->UV)[i].Heat_Phot_HeI = dTableData[i*nTableColumns+6];
		(cl->UV)[i].Heat_Phot_HeII = dTableData[i*nTableColumns+7];
		(cl->UV)[i].Heat_Phot_H2 = dTableData[i*nTableColumns+8];

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
  int i,j;
  double DeltaTln, Tln, T;

  assert(cl!=NULL);
  assert(cl->RT == NULL);

  cl->R.Cool_Coll_HI = CL_eHI*CL_B_gm; 
  cl->R.Cool_Coll_HeI = CL_eHeI*CL_B_gm;
  cl->R.Cool_Coll_HeII = CL_eHeII*CL_B_gm;
  cl->R.Cool_Diel_HeII = (CL_E2HeII+CL_eHeI)*CL_B_gm;
#ifdef NOMOLECULARHCOOLING
  cl->R.Cool_Coll_H2 = 0;
#else
  cl->R.Cool_Coll_H2 = CL_eH2*CL_B_gm;
#endif

  cl->nTable = nTable;
  cl->TMin = TMin;
  cl->TMax = TMax;
  cl->TlnMin = log( TMin );
  cl->TlnMax = log( TMax );
  DeltaTln = ( cl->TlnMax-cl->TlnMin )/( nTable - 1 );
  cl->rDeltaTln = 1./DeltaTln;
  cl->RT = (RATES_T *) malloc( nTable * sizeof(RATES_T) * TABLEFACTOR );
  assert(cl->RT != NULL);

  for ( j=0; j<nTable; j++ ) {
    Tln = cl->TlnMin + DeltaTln*j;
    T=exp(Tln);

    i = j*TABLEFACTOR;
    (cl->RT+i)->Rate_Coll_HI = clRateCollHI( T );
    if ( (cl->RT+i)->Rate_Coll_HI < CL_RT_MIN ) (cl->RT+i)->Rate_Coll_HI = CL_RT_MIN;
    (cl->RT+i)->Rate_Coll_HeI = clRateCollHeI( T );
    if ( (cl->RT+i)->Rate_Coll_HeI < CL_RT_MIN ) (cl->RT+i)->Rate_Coll_HeI = CL_RT_MIN;
    (cl->RT+i)->Rate_Coll_HeII = clRateCollHeII( T );
    if ( (cl->RT+i)->Rate_Coll_HeII < CL_RT_MIN ) (cl->RT+i)->Rate_Coll_HeII = CL_RT_MIN;
    (cl->RT+i)->Rate_Coll_e_H2 = clRateColl_e_H2( T );
    if ( (cl->RT+i)->Rate_Coll_e_H2 < CL_RT_MIN ) (cl->RT+i)->Rate_Coll_e_H2 = CL_RT_MIN;
    (cl->RT+i)->Rate_Coll_HI_H2 = clRateColl_HI_H2( T );
    if ( (cl->RT+i)->Rate_Coll_HI_H2 < CL_RT_MIN ) (cl->RT+i)->Rate_Coll_HI_H2 = CL_RT_MIN;
    (cl->RT+i)->Rate_Coll_H2_H2 = clRateColl_H2_H2( T );
    if ( (cl->RT+i)->Rate_Coll_H2_H2 < CL_RT_MIN ) (cl->RT+i)->Rate_Coll_H2_H2 = CL_RT_MIN;
    (cl->RT+i)->Rate_Coll_Hm_e = clRateColl_Hm_e( T );/*gas phase form of H2 */
    if ( (cl->RT+i)->Rate_Coll_Hm_e < CL_RT_MIN ) (cl->RT+i)->Rate_Coll_Hm_e = CL_RT_MIN;
    (cl->RT+i)->Rate_Coll_Hm_HII = clRateColl_Hm_HII( T );/*gas phase form of H2 */
    if ( (cl->RT+i)->Rate_Coll_Hm_HII < CL_RT_MIN ) (cl->RT+i)->Rate_Coll_Hm_HII = CL_RT_MIN;
    (cl->RT+i)->Rate_Coll_HII_H2 = clRateColl_HII_H2( T );/*gas phase form of H2 */
    if ( (cl->RT+i)->Rate_Coll_HII_H2 < CL_RT_MIN ) (cl->RT+i)->Rate_Coll_HII_H2 = CL_RT_MIN;
    (cl->RT+i)->Rate_HI_e = clRateHI_e( T );/*gas phase form of H2 */
    if ( (cl->RT+i)->Rate_HI_e < CL_RT_MIN ) (cl->RT+i)->Rate_HI_e = CL_RT_MIN;
    (cl->RT+i)->Rate_HI_Hm = clRateHI_Hm( T );/*gas phase form of H2 */
    if ( (cl->RT+i)->Rate_HI_Hm < CL_RT_MIN ) (cl->RT+i)->Rate_HI_Hm = CL_RT_MIN;

    (cl->RT+i)->Rate_Radr_HII = clRateRadrHII( T );
    (cl->RT+i)->Rate_Radr_HeII = clRateRadrHeII( T );
    (cl->RT+i)->Rate_Radr_HeIII = clRateRadrHeIII( T );
    (cl->RT+i)->Rate_Diel_HeII = clRateDielHeII( T );
    (cl->RT+i)->Rate_Chtr_HeII = clRateChtrHeII( T );
   
    (cl->RT+i)->Cool_Brem_1 = clCoolBrem1( T );
    (cl->RT+i)->Cool_Brem_2 = clCoolBrem2( T );
    (cl->RT+i)->Cool_Radr_HII = clCoolRadrHII( T );
    (cl->RT+i)->Cool_Radr_HeII = clCoolRadrHeII( T );
    (cl->RT+i)->Cool_Radr_HeIII = clCoolRadrHeIII( T );
    (cl->RT+i)->Cool_Line_H2_H = clCoolLineH2_HI( T ); /*H2 Rot-Vib transitions out of collisionally-induced excited states. H2-H collisions*/
    (cl->RT+i)->Cool_Line_H2_H2 = clCoolLineH2_H2( T ); /*H2 Rot-Vib transitions out of collisionally-induced excited states. H2-H2 collisions*/
    (cl->RT+i)->Cool_Line_H2_He = clCoolLineH2_He( T ); /*H2 Rot-Vib transitions out of collisionally-induced excited states. H2-He collisions*/
    (cl->RT+i)->Cool_Line_H2_e = clCoolLineH2_e( T ); /*H2 Rot-Vib transitions out of collisionally-induced excited states. H2-e collisions*/
    (cl->RT+i)->Cool_Line_H2_HII = clCoolLineH2_HII( T ); /*H2 Rot-Vib transitions out of collisionally-induced excited states. H2-p collisions*/
    (cl->RT+i)->Cool_Line_HI = clCoolLineHI( T );
    (cl->RT+i)->Cool_Line_HeI = clCoolLineHeI( T );
    (cl->RT+i)->Cool_Line_HeII = clCoolLineHeII( T );
    (cl->RT+i)->Cool_LowT = clCoolLowT( T );

#ifdef CUBICTABLEINTERP
    {
    double Tup = exp(Tln+0.01*DeltaTln);
    double Tdn = exp(Tln-0.01*DeltaTln);
    double Tfac = 1./(0.02*DeltaTln)*DeltaTln; /* Table contains dy/dlnT * DeltaTln */
    
    i++;
    (cl->RT+i)->Rate_Coll_HI = (clRateCollHI( Tup )-clRateCollHI( Tdn ))*Tfac;
    if ( (cl->RT+i-1)->Rate_Coll_HI < CL_RT_MIN ) (cl->RT+i)->Rate_Coll_HI = 0;
    (cl->RT+i)->Rate_Coll_HeI = (clRateCollHeI( Tup )-clRateCollHeI( Tdn ))*Tfac;
    if ( (cl->RT+i-1)->Rate_Coll_HeI < CL_RT_MIN ) (cl->RT+i)->Rate_Coll_HeI = 0;
    (cl->RT+i)->Rate_Coll_HeII = (clRateCollHeII( Tup )-clRateCollHeII( Tdn ))*Tfac;
    if ( (cl->RT+i-1)->Rate_Coll_HeII < CL_RT_MIN ) (cl->RT+i)->Rate_Coll_HeII = 0;
    (cl->RT+i)->Rate_Coll_e_H2 = (clRateColl_e_H2( Tup )-clRateColl_e_H2( Tdn ))*Tfac;
    if ( (cl->RT+i-1)->Rate_Coll_e_H2 < CL_RT_MIN ) (cl->RT+i)->Rate_Coll_e_H2 = 0; 
    (cl->RT+i)->Rate_Coll_HI_H2 = (clRateColl_HI_H2( Tup )-clRateColl_HI_H2( Tdn ))*Tfac;
    if ( (cl->RT+i-1)->Rate_Coll_HI_H2 < CL_RT_MIN ) (cl->RT+i)->Rate_Coll_HI_H2 = 0; 
    (cl->RT+i)->Rate_Coll_H2_H2 = (clRateColl_H2_H2( Tup )-clRateColl_H2_H2( Tdn ))*Tfac;
    if ( (cl->RT+i-1)->Rate_Coll_H2_H2 < CL_RT_MIN ) (cl->RT+i)->Rate_Coll_H2_H2 = 0; 
    (cl->RT+i)->Rate_Coll_Hm_e = (clRateColl_Hm_e( Tup )-clRateColl_Hm_e( Tdn ))*Tfac;/*gas phase form of H2 */
    if ( (cl->RT+i)->Rate_Coll_Hm_e < CL_RT_MIN ) (cl->RT+i)->Rate_Coll_Hm_e = 0;
    (cl->RT+i)->Rate_Coll_Hm_HII = (clRateColl_Hm_HII( Tup )-clRateColl_Hm_HII( Tdn ))*Tfac;/*gas phase form of H2 */
    if ( (cl->RT+i)->Rate_Coll_Hm_HII < CL_RT_MIN ) (cl->RT+i)->Rate_Coll_Hm_HII = 0;
    (cl->RT+i)->Rate_Coll_HII_H2 = (clRateColl_HII_H2( Tup )-clRateColl_HII_H2( Tdn ))*Tfac;/*gas phase form of H2 */
    if ( (cl->RT+i)->Rate_Coll_HII_H2 < CL_RT_MIN ) (cl->RT+i)->Rate_Coll_HII_H2 = 0;
    (cl->RT+i)->Rate_HI_e = (clRate_HI_e( Tup )-clRate_HI_e( Tdn ))*Tfac;/*gas phase form of H2 */
    if ( (cl->RT+i)->Rate_HI_e < CL_RT_MIN ) (cl->RT+i)->Rate_HI_e = 0;
    (cl->RT+i)->Rate_HI_Hm = (clRate_HI_Hm( Tup )-clRate_HI_Hm( Tdn ))*Tfac;/*gas phase form of H2 */
    if ( (cl->RT+i)->Rate_HI_Hm < CL_RT_MIN ) (cl->RT+i)->Rate_HI_Hm = 0;

    (cl->RT+i)->Rate_Radr_HII = ( clRateRadrHII( Tup )-clRateRadrHII( Tdn ))*Tfac;
    (cl->RT+i)->Rate_Radr_HeII = ( clRateRadrHeII( Tup )-clRateRadrHeII( Tdn ))*Tfac;
    (cl->RT+i)->Rate_Radr_HeIII = ( clRateRadrHeIII( Tup )-clRateRadrHeIII( Tdn ))*Tfac;
    (cl->RT+i)->Rate_Diel_HeII = ( clRateDielHeII( Tup )-clRateDielHeII( Tdn ))*Tfac;
    (cl->RT+i)->Rate_Chtr_HeII = ( clRateChtrHeII( Tup )-clRateChtrHeII( Tdn ))*Tfac;
    
    (cl->RT+i)->Cool_Brem_1 = ( clCoolBrem1( Tup )-clCoolBrem1( Tdn ))*Tfac;
    (cl->RT+i)->Cool_Brem_2 = ( clCoolBrem2( Tup )-clCoolBrem2( Tdn ))*Tfac;
    (cl->RT+i)->Cool_Radr_HII = ( clCoolRadrHII( Tup )-clCoolRadrHII( Tdn ))*Tfac;
    (cl->RT+i)->Cool_Radr_HeII = ( clCoolRadrHeII( Tup )-clCoolRadrHeII( Tdn ))*Tfac;
    (cl->RT+i)->Cool_Radr_HeIII = ( clCoolRadrHeIII( Tup )-clCoolRadrHeIII( Tdn ))*Tfac;
    (cl->RT+i)->Cool_Line_HI = ( clCoolLineHI( Tup )-clCoolLineHI( Tdn ))*Tfac;
    (cl->RT+i)->Cool_Line_HeI = ( clCoolLineHeI( Tup )-clCoolLineHeI( Tdn ))*Tfac; 
    (cl->RT+i)->Cool_Line_HeII = ( clCoolLineHeII( Tup )-clCoolLineHeII( Tdn ))*Tfac;
    /* (cl->RT+i)->Cool_Line_H2 =  ( clCoolLineH2_table( Tup )-clCoolLineH2_table( Tdn ))*Tfac;*/ /*H2 Rot-Vib transitions CC */
    (cl->RT+i)->Cool_Line_H2_H = ( clCoolLineH2_HI( Tup ) - clCoolLineH2_HI( Tdn ))*Tfac; /*H2 Rot-Vib transitions out of collisionally-induced excited states. H2-H collisions*/
    (cl->RT+i)->Cool_Line_H2_H2 = ( clCoolLineH2_H2( Tup ) - clCoolLineH2_H2( Tdn ))*Tfac; /*H2 Rot-Vib transitions out of collisionally-induced excited states. H2-H2 collisions*/
    (cl->RT+i)->Cool_Line_H2_He = ( clCoolLineH2_He( Tup ) - clCoolLineH2_He( Tdn ))*Tfac; /*H2 Rot-Vib transitions out of collisionally-induced excited states. H2-He collisions*/
    (cl->RT+i)->Cool_Line_H2_e = ( clCoolLineH2_e( Tup ) - clCoolLineH2_e( Tdn ))*Tfac; /*H2 Rot-Vib transitions out of collisionally-induced excited states. H2-e collisions*/
    (cl->RT+i)->Cool_Line_H2_HII = ( clCoolLineH2_HII( Tup ) - clCoolLineH2_HII( Tdn ))*Tfac; /*H2 Rot-Vib transitions out of collisionally-induced excited states. H2-p collisions*/

    (cl->RT+i)->Cool_LowT = ( clCoolLowT( Tup )-clCoolLowT( Tdn ))*Tfac;
    }
#endif
  }    

#if (0)
#ifdef CUBICTABLEINTERP
  for ( j=0; j<nTable; j++ ) {
    double xTln,wTln0,wTln1,wTln0d,wTln1d;
    RATES_T *RT0,*RT1,*RT0d,*RT1d;
    int iTln;
    
    Tln = cl->TlnMin + DeltaTln*j;
    T=exp(Tln);

    i = j*TABLEFACTOR+1;

    xTln = (Tln+0.5/cl->rDeltaTln-cl->TlnMin)*cl->rDeltaTln; /* midpoint */
    iTln = xTln;
    RT0 = (cl->RT+iTln*TABLEFACTOR);
    RT1 = RT0+TABLEFACTOR; 
    xTln = xTln-iTln;
    RT0d = RT0+1;
    RT1d = RT1+1;
	{
	double x2 = xTln*xTln;
	wTln1 = x2*(3-2*xTln);
	wTln0 = 1-wTln1;
	wTln0d = xTln*(1+xTln*(xTln-2));
	wTln1d = x2*(xTln-1);
/*	wTln1 = xTln;
	wTln0 = 1-xTln;
	wTln0d = 0;
	wTln1d = 0;*/
	}

	/*    if ((j%50)==20) {
	printf("%f  %e %e  %e _%e_ %e\n",T,((cl->RT+i-1)->Rate_Radr_HII-(cl->RT+i-3)->Rate_Radr_HII)/(T),(cl->RT+i-2)->Rate_Radr_HII, (cl->RT+i-1)->Rate_Radr_HII, TABLEINTERP( Rate_Radr_HII ), (cl->RT+i+1)->Rate_Radr_HII );
	}*/
      }
#endif
#endif

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
  double Y_H, Y_He, Y_eMax;
  double ten = 10.0,output, expon;

  /* printf("Redshift: %f \n", zIn); */ 
  /* cl->z = 0.0; */ 
  cl->z = zIn; 
  cl->dTime = dTimeIn;
  cl->dComovingGmPerCcUnit = cl->dGmPerCcUnit*pow(1.+zIn,3.);
  cl->dExpand = 1.0/(1.0+zIn);

  cl->R.Cool_Comp = pow((1+zIn)*CL_Ccomp,4.0)*CL_B_gm; 
  cl->R.Tcmb = CL_Tcmb0*(1+zIn);
  clSetAbundanceTotals(cl,0.0,&Y_H,&Y_He,&Y_eMax); /* Hack to estimate Y_H */
  cl->R.Cool_LowTFactor = (cl->bLowTCool ? CL_B_gm*Y_H*Y_H/0.001 : 0 );

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
	  cl->R.Rate_Phot_H2_cosmo = CL_RT_MIN; 
	  
	  cl->R.Heat_Phot_HI = 0.0;
	  cl->R.Heat_Phot_HeI = 0.0;
	  cl->R.Heat_Phot_HeII = 0.0;
	  cl->R.Heat_Phot_H2 = 0.0; 
	  return;
	  }
  
  UV0=UV-1;  
  if (i == cl->nUV ) {
	  cl->R.Rate_Phot_HI = UV0->Rate_Phot_HI;
	  cl->R.Rate_Phot_HeI = UV0->Rate_Phot_HeI;
	  cl->R.Rate_Phot_HeII = UV0->Rate_Phot_HeII;
	  expon = 0.90632725*zTime - 0.16790918*zTime*zTime + 0.010241484*zTime*zTime*zTime - 12.518825;
	  output = pow(ten,expon); /*haardt_madau gal+quasar, z = 0*/ 
	  cl->R.Rate_Phot_H2_cosmo = output;
	  cl->R.Heat_Phot_HI = UV0->Heat_Phot_HI*CL_B_gm;
	  cl->R.Heat_Phot_HeI = UV0->Heat_Phot_HeI*CL_B_gm;
	  cl->R.Heat_Phot_HeII = UV0->Heat_Phot_HeII*CL_B_gm;
	  cl->R.Heat_Phot_H2 = 6.4e-13*CL_B_gm;
	  }
  else {
	  if (cl->bUVTableLinear) { /* use Linear interpolation */	
		  xx = (zTime - UV0->zTime)/(UV->zTime - UV0->zTime);
		  cl->R.Rate_Phot_HI = UV0->Rate_Phot_HI*(1-xx)+UV->Rate_Phot_HI*xx;
		  cl->R.Rate_Phot_HeI = UV0->Rate_Phot_HeI*(1-xx)+UV->Rate_Phot_HeI*xx;
		  cl->R.Rate_Phot_HeII = UV0->Rate_Phot_HeII*(1-xx)+UV->Rate_Phot_HeII*xx;
		  expon = 0.90632725*zTime - 0.16790918*zTime*zTime + 0.010241484*zTime*zTime*zTime - 12.518825;
		  cl->R.Rate_Phot_H2_cosmo = pow(ten,expon); /*haardt_madau gal+quasar, z = 0*/ 
		  cl->R.Heat_Phot_HI = (UV0->Heat_Phot_HI*(1-xx)+UV->Heat_Phot_HI*xx)*CL_B_gm;
		  cl->R.Heat_Phot_HeI = (UV0->Heat_Phot_HeI*(1-xx)+UV->Heat_Phot_HeI*xx)*CL_B_gm;
		  cl->R.Heat_Phot_HeII = (UV0->Heat_Phot_HeII*(1-xx)+UV->Heat_Phot_HeII*xx)*CL_B_gm;
		  cl->R.Heat_Phot_H2 = 6.4e-13*CL_B_gm;
		  }
	  else { /* use Log interpolation with 1+zTime */
		  xx = log((1+zTime)/(1+UV0->zTime))/log((1+UV->zTime)/(1+UV0->zTime));
		  cl->R.Rate_Phot_HI = pow(UV0->Rate_Phot_HI,1-xx)*pow(UV->Rate_Phot_HI,xx);
		  cl->R.Rate_Phot_HeI = pow(UV0->Rate_Phot_HeI,1-xx)*pow(UV->Rate_Phot_HeI,xx);
		  cl->R.Rate_Phot_HeII = pow(UV0->Rate_Phot_HeII,1-xx)*pow(UV->Rate_Phot_HeII,xx);
		  expon = 0.90632725*zTime - 0.16790918*zTime*zTime + 0.010241484*zTime*zTime*zTime - 12.518825;
		  cl->R.Rate_Phot_H2_cosmo = pow(ten,expon);  /*haardt_madau gal+quasar, z = 0*/ 
		  
		  cl->R.Heat_Phot_HI = pow(UV0->Heat_Phot_HI,1-xx)*pow(UV->Heat_Phot_HI,xx)*CL_B_gm;
		  cl->R.Heat_Phot_HeI = pow(UV0->Heat_Phot_HeI,1-xx)*pow(UV->Heat_Phot_HeI,xx)*CL_B_gm;
		  cl->R.Heat_Phot_HeII = pow(UV0->Heat_Phot_HeII,1-xx)*pow(UV->Heat_Phot_HeII,xx)*CL_B_gm;
		  cl->R.Heat_Phot_H2 = 6.4e-13*CL_B_gm;
		  }
	  }
  if (cl->R.Rate_Phot_HI < CL_RT_MIN) cl->R.Rate_Phot_HI = CL_RT_MIN;
  if (cl->R.Rate_Phot_HeI < CL_RT_MIN) cl->R.Rate_Phot_HeI = CL_RT_MIN;
  if (cl->R.Rate_Phot_HeII < CL_RT_MIN) cl->R.Rate_Phot_HeII = CL_RT_MIN;
  if (cl->R.Rate_Phot_H2_cosmo < CL_RT_MIN) cl->R.Rate_Phot_H2_cosmo = CL_RT_MIN; 

  return;
  }

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
			1.2294426136470751e-6, 0.0, 0.0, 0.0, 0.0, 0.0 };

double AP_Gamma_H2_factor[] = {1.0, 1.0, 1.0,
			  1.0, 1.0, 1.0,
			  1.0, 1.0, 1.0,
			  1.0, 1.0, 1.0,
			  1.0, 1.0, 1.0,
			  1.0, 1.0, 1.0,
			  1.0, 1.0, 1.0,
			  1.0, 1.0, 1.0, 1.0, 1.0}; 
/*This relates the density of the surroundings to the likelihood of photodissociation
we have another way to do shielding for H2 so this is filled in only to prevent a crash

It's also untested with H2
*/

void clRates( COOL *cl, RATE *Rate, double T, double rho, double ZMetal, double columnL, double Rate_Phot_H2_stellar) {
  double Tln;

  if (T >= cl->TMax) T=cl->TMax*(1.0 - EPS);   
  if (T < cl->TMin) T=cl->TMin;
  Tln = log(T); /* Deprecated but log's, sqrt's etc... used in raw rate functions */

  Rate->T = T;
  Rate->Tln = Tln; 
  Rate->Coll_HI = clRateCollHI( T );
  Rate->Coll_HeI = clRateCollHeI( T );
  Rate->Coll_HeII = clRateCollHeII( T );
  Rate->Coll_e_H2 = clRateColl_e_H2( T );
  Rate->Coll_HI_H2 = clRateColl_HI_H2( T );
  Rate->Coll_H2_H2 = clRateColl_H2_H2( T );
  Rate->Coll_Hm_e = clRateColl_Hm_e(T);          /*gas phase formation of H2 */
  Rate->Coll_Hm_HII = clRateColl_Hm_HII(T);           /*gas phase formation of H2 */
  Rate->Coll_HII_H2 = clRateColl_HII_H2(T);          /*gas phase formation of H2 */
  Rate->HI_e = clRateHI_e(T);          /*gas phase formation of H2 */
  Rate->HI_Hm = clRateHI_Hm(T);          /*gas phase formation of H2 */

  Rate->Radr_HII = clRateRadrHII( T );
  Rate->Radr_HeII = clRateRadrHeII( T );
  Rate->Diel_HeII = clRateDielHeII( T );
  Rate->Chtr_HeII = clRateChtrHeII( T );  
  Rate->Totr_HeII = Rate->Radr_HeII + Rate->Diel_HeII + Rate->Chtr_HeII;
  Rate->Radr_HeIII = clRateRadrHeIII( T );
  Rate->DustForm_H2 = clRateDustFormH2( ZMetal, cl->dClump  ); /*H2 Formation on dust*/
  Rate->Phot_HI = cl->R.Rate_Phot_HI;
  Rate->Phot_HeI = cl->R.Rate_Phot_HeI;
  Rate->Phot_HeII = cl->R.Rate_Phot_HeII;

  Rate->Phot_H2 = cl->R.Rate_Phot_H2_cosmo + pow(10,Rate_Phot_H2_stellar - 58.9079)*cl->dLymanWernerFrac/cl->dKpcUnit/cl->dKpcUnit/cl->dExpand/cl->dExpand;
  Rate->CorreLength = columnL * cl->dKpcUnit * 3.08568025e21 * cl->dExpand; /* From system units to cm*/
  Rate->LymanWernerCode = Rate_Phot_H2_stellar;
  if (cl->bSelfShield) {
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

#define TABLEINTERPLIN( _rname ) (wTln0*RT0->_rname+wTln1*RT1->_rname)
void clRates_Table_Lin( COOL *cl, RATE *Rate, double T, double rho, double ZMetal, double columnL, double Rate_Phot_H2_stellar) {
  double Tln;
  double xTln,wTln0,wTln1;/*,wTln0d,wTln1d;*/
  RATES_T *RT0,*RT1; /**RT0d,*RT1d;*/
  int iTln;

  if (T >= cl->TMax) T=cl->TMax*(1.0 - EPS);   
  if (T < cl->TMin) T=cl->TMin;
  Tln = log(T);

  Rate->T = T;
  Rate->Tln = Tln; 

  xTln = (Tln-cl->TlnMin)*cl->rDeltaTln;
  iTln = xTln;
  RT0 = (cl->RT+iTln*TABLEFACTOR);
  RT1 = RT0+TABLEFACTOR; 
  xTln = xTln-iTln;
  wTln1 = xTln;
  wTln0 = 1-xTln;
  Rate->Coll_HI = TABLEINTERPLIN( Rate_Coll_HI );
  Rate->Coll_e_H2 = TABLEINTERPLIN( Rate_Coll_e_H2 ); 
  Rate->Coll_HI_H2 = TABLEINTERPLIN( Rate_Coll_HI_H2 ); 
  Rate->Coll_H2_H2 = TABLEINTERPLIN( Rate_Coll_H2_H2 ); 
  Rate->Coll_HeI = TABLEINTERPLIN( Rate_Coll_HeI );
  Rate->Coll_HeII = TABLEINTERPLIN( Rate_Coll_HeII );
  Rate->Coll_Hm_e = TABLEINTERPLIN(Rate_Coll_Hm_e);          /*gas phase form of H2 */
  Rate->Coll_Hm_HII = TABLEINTERPLIN(Rate_Coll_Hm_HII);           /*gas phase form of H2 */
  Rate->Coll_HII_H2 = TABLEINTERPLIN(Rate_Coll_HII_H2);         
  Rate->HI_e = TABLEINTERPLIN(Rate_HI_e);          /*gas phase form of H2 */
  Rate->HI_Hm = TABLEINTERPLIN(Rate_HI_Hm);          /*gas phase form of H2 */

  Rate->Radr_HII = TABLEINTERPLIN( Rate_Radr_HII );
  Rate->Radr_HeII = TABLEINTERPLIN( Rate_Radr_HeII );
  Rate->Diel_HeII = TABLEINTERPLIN( Rate_Diel_HeII );
  Rate->Chtr_HeII = TABLEINTERPLIN( Rate_Chtr_HeII ); 
  Rate->Totr_HeII = Rate->Radr_HeII + Rate->Diel_HeII + Rate->Chtr_HeII;
  Rate->Radr_HeIII = TABLEINTERPLIN( Rate_Radr_HeIII );
  Rate->DustForm_H2 = clRateDustFormH2( ZMetal, cl->dClump );

  Rate->Phot_HI = cl->R.Rate_Phot_HI;
  Rate->Phot_HeI = cl->R.Rate_Phot_HeI;
  Rate->Phot_HeII = cl->R.Rate_Phot_HeII;
  Rate->Phot_H2 = cl->R.Rate_Phot_H2_cosmo + pow(10,Rate_Phot_H2_stellar - 58.9079)*cl->dLymanWernerFrac/cl->dKpcUnit/cl->dKpcUnit/cl->dExpand/cl->dExpand; 
  Rate->CorreLength = columnL * cl->dKpcUnit * 3.08568025e21 * cl->dExpand; /* From system units to cm*/
  Rate->LymanWernerCode = Rate_Phot_H2_stellar;
  if (cl->bSelfShield) {
      double logen_B;
      logen_B = log10(rho*CL_B_gm);
      if (logen_B > 2.2499) {
	  Rate->Phot_HI = 0;
	  Rate->Phot_HeI = 0;
	  Rate->Phot_HeII = 0;
	  Rate->Phot_H2 = 0; 
	  }
      else if (logen_B > -10.25) {
	  double x = (logen_B+10.25)*2.0;
	  int ix;
	  ix = floor(x);
	  x -= ix;
	  Rate->Phot_HI *= (AP_Gamma_HI_factor[ix]*(1-x)+AP_Gamma_HI_factor[ix+1]*x);
	  Rate->Phot_HeI *= (AP_Gamma_HeI_factor[ix]*(1-x)+AP_Gamma_HeI_factor[ix+1]*x);
	  Rate->Phot_HeII *= (AP_Gamma_HeII_factor[ix]*(1-x)+AP_Gamma_HeII_factor[ix+1]*x);
	  Rate->Phot_H2 += (AP_Gamma_H2_factor[ix]*(1-x)+AP_Gamma_H2_factor[ix+1]*x);
	  }
      }
}


void clRates_Table( COOL *cl, RATE *Rate, double T, double rho, double ZMetal, double columnL, double Rate_Phot_H2_stellar) {
  double Tln;
  double xTln,wTln0,wTln1;/*,wTln0d,wTln1d;*/
  RATES_T *RT0,*RT1;/*,*RT0d,*RT1d;*/
  int iTln;

#ifdef TESTRATE
  RATE test;
  clRates( cl, &test, T, rho, ZMetal, columnL, Rate_Phot_H2_stellar);
#endif

  if (T >= cl->TMax) T=cl->TMax*(1.0 - EPS);   
  if (T < cl->TMin) T=cl->TMin;
  Tln = log(T);

  Rate->T = T;
  Rate->Tln = Tln; 

  xTln = (Tln-cl->TlnMin)*cl->rDeltaTln;
  iTln = xTln;
  RT0 = (cl->RT+iTln*TABLEFACTOR);
  RT1 = RT0+TABLEFACTOR; 
  xTln = xTln-iTln;
#ifdef CUBICTABLEINTERP
  RT0d = RT0+1;
  RT1d = RT1+1;
  {
  double x2 = xTln*xTln;
  wTln1 = x2*(3-2*xTln);
  wTln0 = 1-wTln1;
  wTln0d = xTln*(1+xTln*(xTln-2));
  wTln1d = x2*(xTln-1);
/*  wTln1 = xTln;
  wTln0 = 1-xTln;
  wTln0d = 0;
  wTln1d = 0;*/
  }
#else  
  wTln1 = xTln;
  wTln0 = 1-xTln;
#endif
  Rate->Coll_HI = TABLEINTERP( Rate_Coll_HI );
  Rate->Coll_HeI = TABLEINTERP( Rate_Coll_HeI );
  Rate->Coll_HeII = TABLEINTERP( Rate_Coll_HeII );
  Rate->Coll_e_H2 = TABLEINTERP( Rate_Coll_e_H2 );
  Rate->Coll_HI_H2 = TABLEINTERP( Rate_Coll_HI_H2 );
  Rate->Coll_H2_H2 = TABLEINTERP( Rate_Coll_H2_H2 );
  Rate->Coll_HII_H2 = TABLEINTERPLIN(Rate_Coll_HII_H2);        
  Rate->Coll_Hm_e = TABLEINTERPLIN(Rate_Coll_Hm_e);          /*gas phase form of H2 */
  Rate->Coll_Hm_HII = TABLEINTERPLIN(Rate_Coll_Hm_HII);           /*gas phase form of H2 */
  Rate->HI_e = TABLEINTERPLIN(Rate_HI_e);          /*gas phase form of H2 */
  Rate->HI_Hm = TABLEINTERPLIN(Rate_HI_Hm);          /*gas phase form of H2 */
  Rate->Radr_HII = TABLEINTERP( Rate_Radr_HII );
  Rate->Radr_HeII = TABLEINTERP( Rate_Radr_HeII );
  Rate->Diel_HeII = TABLEINTERP( Rate_Diel_HeII );
  Rate->Chtr_HeII = TABLEINTERP( Rate_Chtr_HeII );
  Rate->Totr_HeII = Rate->Radr_HeII + Rate->Diel_HeII + Rate->Chtr_HeII;
  Rate->Radr_HeIII = TABLEINTERP( Rate_Radr_HeIII );
  Rate->DustForm_H2 =  clRateDustFormH2( ZMetal, cl->dClump );

  Rate->Phot_HI = cl->R.Rate_Phot_HI;
  Rate->Phot_HeI = cl->R.Rate_Phot_HeI;
  Rate->Phot_HeII = cl->R.Rate_Phot_HeII;
  Rate->Phot_H2 = cl->R.Rate_Phot_H2_cosmo + pow(10,Rate_Phot_H2_stellar - 58.9079)*cl->dLymanWernerFrac/cl->dKpcUnit/cl->dKpcUnit/cl->dExpand/cl->dExpand;
  Rate->CorreLength = columnL * cl->dKpcUnit*3.08568025e21 * cl->dExpand; /* From system units to cm*/
  Rate->LymanWernerCode = Rate_Phot_H2_stellar;
  if (cl->bSelfShield) {
      double logen_B;
      logen_B = log10(rho*CL_B_gm);
      if (logen_B > 2.2499) {
	  Rate->Phot_HI = 0;
	  Rate->Phot_HeI = 0;
	  Rate->Phot_HeII = 0;
	  Rate->Phot_H2 = 0;
	  }
      else if (logen_B > -10.25) {
	  double x = (logen_B+10.25)*2.0;
	  int ix;
	  ix = floor(x);
	  x -= ix;
	  Rate->Phot_HI *= (AP_Gamma_HI_factor[ix]*(1-x)+AP_Gamma_HI_factor[ix+1]*x);
	  Rate->Phot_HeI *= (AP_Gamma_HeI_factor[ix]*(1-x)+AP_Gamma_HeI_factor[ix+1]*x);
	  Rate->Phot_HeII *= (AP_Gamma_HeII_factor[ix]*(1-x)+AP_Gamma_HeII_factor[ix+1]*x);
	  Rate->Phot_H2 *= (AP_Gamma_H2_factor[ix]*(1-x)+AP_Gamma_H2_factor[ix+1]*x); /* Untested */
	  }
      }

#ifdef TESTRATE
#define RATEVAR( _name ) (fabs((test._name - Rate->_name)/(test._name)) > TESTRATE)
  if ( (T <.5e6 && rho > 1e3*1e-31) && ((RATEVAR(Coll_HI)) || RATEVAR(Coll_HeI) || RATEVAR(Coll_HeII) || 
       RATEVAR(Radr_HII) || RATEVAR(Radr_HeII) || RATEVAR(Diel_HeII) || RATEVAR(Totr_HeII) || RATEVAR(Radr_HeIII) || 
					RATEVAR(Phot_HI) || RATEVAR(Phot_HeI) || RATEVAR(Phot_HeII) || RATEVAR(Phot_H2) )) { 
    printf("Bad interpolated rates at T=%e rho=%e\n",T,rho);
#ifdef CUBICTABLEINTERP
    printf("%e %e %e %e\n %e %e %e %e\n", RT0->Rate_Coll_HI, RT1->Rate_Coll_HI, RT0d->Rate_Coll_HI, RT1d->Rate_Coll_HI,wTln0, wTln1, wTln0d, wTln1d );
#endif
    printf("%12f %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e %e12 %e12 %e12\n", 
               T,  test.Coll_HI,  test.Coll_e_H2,  test.Coll_HI_H2,   test.Coll_H2_H2,  test.Coll_HeI,   test.Coll_HeII,   test.Radr_HII,   test.Radr_HeII,   test.Diel_HeII,  test.Totr_HeII,  test.Radr_HeIII,  test.Phot_HI,      test.Phot_H2,    test.Phot_HeI,  test.Phot_HeII);
    printf("%12f %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e \n", 
               T, Rate->Coll_HI,  Rate->Coll_e_H2, Rate->Coll_HI_H2,  Rate->Coll_H2_H2, Rate->Coll_HeI,  Rate->Coll_HeII,  Rate->Radr_HII,  Rate->Radr_HeII,  Rate->Diel_HeII, Rate->Chtr_HeII, Rate->Totr_HeII,  Rate->Radr_HeIII,  Rate->Phot_HI,   Rate->Phot_HI, Rate->Phot_HeI,  Rate->Phot_HeII );
    clRates_Table_Lin( cl, &test, T, rho, ZMetal, columnL, Rate_Phot_H2_stellar);
    printf("%12f %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e %e12 %12e %e12\n", T, test.Coll_HI,  test.Coll_e_H2, test.Coll_HI_H2, test.Coll_H2_H2, test.Coll_HeI,  test.Coll_HeII,  test.Radr_HII,  test.Radr_HeII,  test.Diel_HeII,  test.Totr_HeII,  test.Radr_HeIII,  test.Phot_HI,  test.Phot_H2, test.Phot_HeI,  test.Phot_HeII );
    clRates( cl, &test, T*1.00123, rho, ZMetal, columnL, Rate_Phot_H2_stellar);
    printf("%12f %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e %e12 %12e %e12\n", T*1.00123, test.Coll_HI, test.Coll_e_H2, test.Coll_HI_H2, test.Coll_H2_H2, test.Coll_HeI,  test.Coll_HeII,  test.Radr_HII,  test.Radr_HeII,  test.Diel_HeII,  test.Totr_HeII,  test.Radr_HeIII,  test.Phot_HI, test.Phot_H2, test.Phot_HeI,  test.Phot_HeII );
    clRates( cl, &test, T/1.00123, rho, ZMetal, columnL, Rate_Phot_H2_stellar);
    printf("%12f %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e %12e %e12 %12e %12e\n", T/1.00123, test.Coll_HI, test.Coll_e_H2, test.Coll_HI_H2, test.Coll_H2_H2, test.Coll_HeI,  test.Coll_HeII,  test.Radr_HII,  test.Radr_HeII,  test.Diel_HeII,  test.Totr_HeII,  test.Radr_HeIII,  test.Phot_HI, test.Phot_H2, test.Phot_HeI,  test.Phot_HeII );
    assert(0);
  }
#endif
}


void clRateMetalTable(COOL *cl, RATE *Rate, double T, double rho, double Y_H, double ZMetal)
{
  double tempT, tempnH,tempz, nH;
  double Tlog, nHlog; 
  double xTlog, wTlog0, wTlog1, xz, wz0, wz1, xnHlog, wnHlog0, wnHlog1; 
  int    iTlog, iz, inHlog; 
  double  Cool000, Cool010, Cool100, Cool110, Cool001, Cool011, Cool101, Cool111; 
  double  Cool00, Cool01, Cool10, Cool11, Cool0, Cool1, Cool;
  
  double  Heat000, Heat010, Heat100, Heat110, Heat001, Heat011, Heat101, Heat111; 
  double  Heat00, Heat01, Heat10, Heat11, Heat0, Heat1, Heat;

  
  if(!cl->bMetal) {
    Rate->Cool_Metal = 0.0; 
    Rate->Heat_Metal = 0.0; 
    return; 
  }
 
  nH = rho*Y_H/M_H;
  
  tempT = T; 
  tempnH = nH;
  tempz = cl->z; 

  if (T >= cl->MetalTMax) tempT = cl->MetalTMax*(1.0-EPS);
  if (T < cl->MetalTMin) tempT=cl->MetalTMin;

  if (nH >= cl->MetalnHMax) tempnH = cl->MetalnHMax*(1.0-EPS);
  if (nH < cl->MetalnHMin) tempnH = cl->MetalnHMin; 
 
  if (cl->z <= cl->MetalzMin) tempz = cl->MetalzMin+EPS;
  /* if redshift is too high or no UV, use the no UV metal cooling table*/
  if (cl->z > cl->MetalzMax || !cl->bUV) tempz = cl->MetalzMax;   

  Tlog = log10(tempT); 
  nHlog = log10(tempnH); 
  
  xz = cl->nzMetalTable -1 - (tempz - cl->MetalzMin)*cl->rDeltaz; 
  iz = xz;   

  xTlog = (Tlog - cl->MetalTlogMin)*cl->rDeltaTlog; 
  assert(xTlog >= 0.0);
  iTlog = xTlog; 

  xnHlog = (nHlog - cl->MetalnHlogMin)*cl->rDeltanHlog; 
  inHlog = xnHlog;
  if (inHlog == cl->nnHMetalTable - 1) inHlog == cl->nnHMetalTable - 2; /*CC; To prevent running over the table.  Should not be used*/
  
  Cool000 = cl->MetalCoolln[iz][inHlog][iTlog];
  Cool001 = cl->MetalCoolln[iz][inHlog][iTlog+1];
  Cool010 = cl->MetalCoolln[iz][inHlog+1][iTlog];
  Cool011 = cl->MetalCoolln[iz][inHlog+1][iTlog+1];
  Cool100 = cl->MetalCoolln[iz+1][inHlog][iTlog];
  Cool101 = cl->MetalCoolln[iz+1][inHlog][iTlog+1];
  Cool110 = cl->MetalCoolln[iz+1][inHlog+1][iTlog];
  Cool111 = cl->MetalCoolln[iz+1][inHlog+1][iTlog+1];

  Heat000 = cl->MetalHeatln[iz][inHlog][iTlog];
  Heat001 = cl->MetalHeatln[iz][inHlog][iTlog+1];
  Heat010 = cl->MetalHeatln[iz][inHlog+1][iTlog];
  Heat011 = cl->MetalHeatln[iz][inHlog+1][iTlog+1];
  Heat100 = cl->MetalHeatln[iz+1][inHlog][iTlog];
  Heat101 = cl->MetalHeatln[iz+1][inHlog][iTlog+1];
  Heat110 = cl->MetalHeatln[iz+1][inHlog+1][iTlog];
  Heat111 = cl->MetalHeatln[iz+1][inHlog+1][iTlog+1];

  xz = xz - iz; 
  wz1 = xz; 
  wz0 = 1-xz; 
  
  Cool00 = wz0*Cool000 + wz1*Cool100;
  Cool01 = wz0*Cool001 + wz1*Cool101; 
  Cool10 = wz0*Cool010 + wz1*Cool110; 
  Cool11 = wz0*Cool011 + wz1*Cool111;

  Heat00 = wz0*Heat000 + wz1*Heat100;
  Heat01 = wz0*Heat001 + wz1*Heat101; 
  Heat10 = wz0*Heat010 + wz1*Heat110; 
  Heat11 = wz0*Heat011 + wz1*Heat111;


  xnHlog = xnHlog - inHlog; 
  wnHlog1 = xnHlog; 
  wnHlog0 = 1-xnHlog; 

  Cool0 = wnHlog0*Cool00 + wnHlog1*Cool10; 
  Cool1 = wnHlog0*Cool01 + wnHlog1*Cool11; 

  
  Heat0 = wnHlog0*Heat00 + wnHlog1*Heat10; 
  Heat1 = wnHlog0*Heat01 + wnHlog1*Heat11; 
  
  xTlog = xTlog - iTlog; 
  wTlog1 = xTlog; 
  wTlog0 = 1 - xTlog; 
  Cool = wTlog0*Cool0 + wTlog1*Cool1; 
  Heat = wTlog0*Heat0 + wTlog1*Heat1; 
    /* convert unit to erg/g/sec, time a factor of nH^2/nH, also scale with metalicity */ 
  Rate->Cool_Metal = exp(Cool)*nH*Y_H/M_H * ZMetal/ZSOLAR; 
  Rate->Heat_Metal = exp(Heat)*nH*Y_H/M_H * ZMetal/ZSOLAR;   

}

/* Deprecated except for testing: use EdotInstant */
/* Need density in here to make this work with Self-Shielding */
double clHeatTotal ( COOL *cl, PERBARYON *Y, RATE *Rate, double rho, double ZMetal ) {
  /* erg /gram /sec
     Note: QQ_* premultiplied by (CL_B_gm*erg_ev) */
  double heating;
  double en_B = rho*CL_B_gm;
  double s_dust, s_self, Rate_Phot_HI;
  s_dust = clDustShield(Y->HI*en_B, Y->H2*en_B, ZMetal, Rate->CorreLength);
  s_self = clSelfShield(Y->H2*en_B, Rate->CorreLength);
  if (cl->bShieldHI) Rate_Phot_HI = Rate->Phot_HI*s_dust;
  else Rate_Phot_HI = Rate->Phot_HI;

  heating = 
    Y->HI   * cl->R.Heat_Phot_HI * Rate_Phot_HI +
    Y->HeI  * cl->R.Heat_Phot_HeI * Rate->Phot_HeI +
    Y->HeII * cl->R.Heat_Phot_HeII * Rate->Phot_HeII
#ifndef NOMOLECULARHCOOLING
    + Y->H2   * cl->R.Heat_Phot_H2 * Rate->Phot_H2*s_dust*s_self
#endif
#ifndef NOMETALCOOLING
    + Rate->Heat_Metal
#endif
    ;

  return heating; 
}

/* Deprecated except for testing: use EdotInstant */
double clCoolTotal ( COOL *cl, PERBARYON *Y, RATE *Rate, double rho, double ZMetal ) {
  /* Assumes clRates called previously */
  /* erg /gram /sec */

  double en_B=rho*CL_B_gm;
  double LowTCool;
  double s_dust, s_self;
  s_dust = clDustShield(Y->HI*en_B, Y->H2*en_B, ZMetal, Rate->CorreLength);
  s_self = clSelfShield(Y->H2*en_B, Rate->CorreLength);

  if (Rate->T > cl->R.Tcmb)
      LowTCool = clCoolLowT(Rate->T)*cl->R.Cool_LowTFactor*en_B*ZMetal;
  else
      LowTCool = 0;

  /* PUT INTO erg/gm/sec */
  return Y->e * ( 
#ifndef NOCOMPTON
    cl->R.Cool_Comp * ( Rate->T - cl->R.Tcmb ) + 
#endif
    en_B * (
    clCoolBrem1(Rate->T) * ( Y->HII + Y->HeII ) +
    clCoolBrem2(Rate->T) * Y->HeIII +

    clCoolRadrHII(Rate->T) * Y->HII * Rate->Radr_HII +
    clCoolRadrHeII(Rate->T) * Y->HeII * Rate->Radr_HeII +
    clCoolRadrHeIII(Rate->T) * Y->HeIII * Rate->Radr_HeIII +
 
    cl->R.Cool_Coll_HI * Y->HI * Rate->Coll_HI +
    cl->R.Cool_Coll_HeI * Y->HeI * Rate->Coll_HeI +
    cl->R.Cool_Coll_HeII * Y->HeII * Rate->Coll_HeII +
    cl->R.Cool_Diel_HeII * Y->HeII * Rate->Diel_HeII +
#ifndef NOMOLECULARHCOOLING
    clCoolLineH2_e(Rate->T) * Y->H2 * CL_B_gm /*Re-added CL_B_gm CC 9/16/10 */  +
    cl->R.Cool_Coll_H2 * Y->H2 * Rate->Coll_e_H2  + /*remove shielding of collisions CC 11/18/10*/
#endif
    clCoolLineHI(Rate->T) * Y->HI +
    clCoolLineHeI(Rate->T) * Y->HeI +
    clCoolLineHeII(Rate->T) * Y->HeII ) ) +

#ifndef NOMOLECULARHCOOLING
    clCoolLineH2_HI(Rate->T) * en_B * Y->H2 * Y->HI * CL_B_gm /*Re-added CL_B_gm CC 9/16/10:
 clCoolLineH2_HI(Rate->T) [erg cm^3 s^-1 H2atom^-1 HIatom^-1] * en_B [atom/cm^3] * Y->H2 [H2atom atom^-1] * Y->HI [HIatom atom^-1] CL_B_gm [atom g^-1] => [ergs s^-1 g^-1]*/ +
    clCoolLineH2_H2(Rate->T) * en_B * Y->H2 * Y->H2 * CL_B_gm  +
    clCoolLineH2_He(Rate->T) * en_B * Y->H2 * Y->HeI * CL_B_gm  + 
    clCoolLineH2_HII(Rate->T) * en_B * Y->H2 * Y->HII * CL_B_gm  +
    en_B * Y->HI * cl->R.Cool_Coll_H2 * Y->H2 * Rate->Coll_HI_H2  + /*remove shielding of collisions CC 11/18/10*/
    en_B * Y->H2 * cl->R.Cool_Coll_H2 * Y->H2 * Rate->Coll_H2_H2  + /*remove shielding of collisions CC 11/18/10*/
    en_B * Y->HII * cl->R.Cool_Coll_H2 * Y->H2 * Rate->Coll_HII_H2 + /*remove shielding of collisions CC 11/18/10*/
#endif
#ifndef NOMETALCOOLING
    Rate->Cool_Metal +
#endif
    LowTCool;  
 
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
  ret.collion_e_H2 = Y->e * en_B * 
    cl->R.Cool_Coll_H2 * Y->H2 * Rate->Coll_e_H2; 
  ret.collion_H_H2 = Y->HI * en_B * 
    cl->R.Cool_Coll_H2 * Y->H2 * Rate->Coll_HI_H2; 
  ret.collion_H2_H2 = Y->H2 * en_B * 
    cl->R.Cool_Coll_H2 * Y->H2 * Rate->Coll_H2_H2; 
  ret.collion_HII_H2 = Y->H2 * en_B *
    cl->R.Cool_Coll_H2 * Y->H2 * Rate->Coll_HII_H2;
  ret.collionHI = Y->e * en_B * 
    cl->R.Cool_Coll_HI * Y->HI * Rate->Coll_HI;
  ret.collionHeI = Y->e * en_B * 
    cl->R.Cool_Coll_HeI * Y->HeI * Rate->Coll_HeI;
  ret.collionHeII = Y->e * en_B * 
    cl->R.Cool_Coll_HeII * Y->HeII * Rate->Coll_HeII;
  ret.dielrecHeII = Y->e * en_B * 
    cl->R.Cool_Diel_HeII * Y->HeII * Rate->Diel_HeII;
  ret.lineH2 = Y->e * en_B * 
    (wTln0*RT0->Cool_Line_H2_H+wTln1*RT1->Cool_Line_H2_H) * Y->H2;    /*Assumes that the only collision that inspires line cooling is H-H2. */
  ret.lineHI = Y->e * en_B * 
    (wTln0*RT0->Cool_Line_HI+wTln1*RT1->Cool_Line_HI) * Y->HI;
  ret.lineHeI = Y->e * en_B * 
    (wTln0*RT0->Cool_Line_HeI+wTln1*RT1->Cool_Line_HeI) * Y->HeI;
  ret.lineHeII = Y->e * en_B * 
    (wTln0*RT0->Cool_Line_HeII+wTln1*RT1->Cool_Line_HeII) * Y->HeII;
  ret.lowT = en_B * 
    (wTln0*RT0->Cool_LowT+wTln1*RT1->Cool_LowT)*cl->R.Cool_LowTFactor*0.001; /* assume metallicity 0.001 */
  ret.NetMetalCool = Rate->Cool_Metal - Rate->Heat_Metal; /* assume clRateMetalTable called first */
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
  RT0 = (cl->RT+iTln);
  RT1 = RT0+1; 
  wTln1 = xTln-iTln;
  wTln0 = 1-wTln1;

  /* PUT INTO erg/gm/sec */
  printf("Density:  %e \n", rho);
  printf("Temperature: %e \n", Rate->T); 
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
  printf("Collisional Electron Dissociation  H2    %e\n",
    Y->e * en_B * 
    cl->R.Cool_Coll_H2 * Y->H2 *  Rate->Coll_e_H2);
  printf("Collisional H Dissociation  H2    %e\n",  
    Y->HI * en_B *
    cl->R.Cool_Coll_H2 * Y->H2 * Rate->Coll_HI_H2);
  printf("Collisional H2 Dissociation  H2    %e\n",  
    Y->H2 * en_B *
    cl->R.Cool_Coll_H2 * Y->H2 * Rate->Coll_H2_H2);
  printf("Collisional HII Dissociation  H2    %e\n",
    Y->H2 * en_B *
	 cl->R.Cool_Coll_H2 * Y->H2 * Rate->Coll_HII_H2);
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
    Y->e * en_B * 
    cl->R.Cool_Diel_HeII * Y->HeII * Rate->Diel_HeII);
  printf("Line cooling H2   %e\n",
    Y->e * en_B * 
	 (wTln0*RT0->Cool_Line_H2_H+wTln1*RT1->Cool_Line_H2_H) * Y->H2);  /*Assumes that the only collision that inspires line cooling is H2-H*/ 
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
  printf("Net Metal Heating %e \n", Rate->Heat_Metal); 
  printf("Net Metal Cooling %e \n", Rate->Cool_Metal);
    }  


/* printf("%e\n",
    Y->e * ( 
    cl->R.Cool_Comp * ( Rate->T - cl->R.Tcmb ))); 
  printf("%e\n",
    Y->e * en_B * (
    (wTln0*RT0->Cool_Brem_1+wTln1*RT1->Cool_Brem_1) * ( Y->HII )) );
  printf("%e\n",
    Y->e * en_B * (
    (wTln0*RT0->Cool_Brem_1+wTln1*RT1->Cool_Brem_1) * ( Y->HeII )) );
  printf("%e\n",
    Y->e * en_B * (
    (wTln0*RT0->Cool_Brem_2+wTln1*RT1->Cool_Brem_2) * Y->HeIII ) );
  printf("%e\n",
    Y->e * en_B * 
    (wTln0*RT0->Cool_Radr_HII+wTln1*RT1->Cool_Radr_HII) * Y->HII * Rate->Radr_HII );
  printf("%e\n",
    Y->e * en_B * 
    (wTln0*RT0->Cool_Radr_HeII+wTln1*RT1->Cool_Radr_HeII) * Y->HeII * Rate->Radr_HeII);
  printf("%e\n",
    Y->e * en_B * 
    (wTln0*RT0->Cool_Radr_HeIII+wTln1*RT1->Cool_Radr_HeIII) * Y->HeIII * Rate->Radr_HeIII);
  printf("%e\n",
    Y->e * en_B * 
    cl->R.Cool_Coll_HI * Y->HI * Rate->Coll_HI);
  printf("%e\n",
    Y->e * en_B * 
    cl->R.Cool_Coll_HeI * Y->HeI * Rate->Coll_HeI);
  printf("%e\n",
    Y->e * en_B * 
    cl->R.Cool_Coll_HeII * Y->HeII * Rate->Coll_HeII);
  printf("%e\n",
    Y->e * en_B * cl->R.Cool_Diel_HeII * Y->HeII * Rate->Diel_HeII);

  printf("%e\n",
    Y->e * en_B * 
    (wTln0*RT0->Cool_Line_HI+wTln1*RT1->Cool_Line_HI) * Y->HI);
  printf("%e\n",
    Y->e * en_B * 
    (wTln0*RT0->Cool_Line_HeI+wTln1*RT1->Cool_Line_HeI) * Y->HeI);
  printf("%e\n",
    Y->e * en_B * 
    (wTln0*RT0->Cool_Line_HeII+wTln1*RT1->Cool_Line_HeII) * Y->HeII );
  printf("%e\n",
    en_B * 
    (wTln0*RT0->Cool_LowT+wTln1*RT1->Cool_LowT)*cl->R.Cool_LowTFactor*0.001); 
    } */
void clPrintCoolFile( COOL *cl, PERBARYON *Y, RATE *Rate, double rho, double ZMetal, FILE *fp ) {
  /* Assumes clRates called previously */
  /* erg /gram /sec */

  double en_B=rho*CL_B_gm;
  double ne = Y->e*en_B;
  double xTln,wTln0,wTln1;
  RATES_T *RT0,*RT1;
  int iTln;
  double s_dust, s_self;
  s_dust = clDustShield(Y->HI*en_B, Y->H2*en_B, ZMetal, Rate->CorreLength);
  s_self = clSelfShield(Y->H2*en_B, Rate->CorreLength);

  xTln = (Rate->Tln-cl->TlnMin)*cl->rDeltaTln;
  iTln = xTln;
  RT0 = (cl->RT+iTln);
  RT1 = RT0+1; 
  wTln1 = xTln-iTln;
  wTln0 = 1-wTln1;

  /* PUT INTO erg/gm/sec */
  fprintf(fp,"Density:  %e \n", rho*CL_B_gm);
  fprintf(fp, "Temperature: %g \n", Rate->T); 
  fprintf(fp,"Compton:  %e\n",
    Y->e * ( 
    cl->R.Cool_Comp * ( Rate->T - cl->R.Tcmb )));
  /*  fprintf(fp,"Cool Brem HII    %e\n",
    Y->e * en_B * (
    (wTln0*RT0->Cool_Brem_1+wTln1*RT1->Cool_Brem_1) * ( Y->HII )) );
  fprintf(fp,"Cool Brem HeII   %e\n",
    Y->e * en_B * (
    (wTln0*RT0->Cool_Brem_1+wTln1*RT1->Cool_Brem_1) * ( Y->HeII )) );
  fprintf(fp,"Cool Brem HeIII  %e\n",
    Y->e * en_B * (
    (wTln0*RT0->Cool_Brem_2+wTln1*RT1->Cool_Brem_2) * Y->HeIII ) );
  fprintf(fp,"Radiative Recombination  HII    %e\n",
    Y->e * en_B * 
    (wTln0*RT0->Cool_Radr_HII+wTln1*RT1->Cool_Radr_HII) * Y->HII * Rate->Radr_HII );
  fprintf(fp,"Radiative Recombination  HeII   %e\n",
    Y->e * en_B * 
    (wTln0*RT0->Cool_Radr_HeII+wTln1*RT1->Cool_Radr_HeII) * Y->HeII * Rate->Radr_HeII);
  fprintf(fp,"Radiative Recombination  HeIII  %e\n",
    Y->e * en_B * 
    (wTln0*RT0->Cool_Radr_HeIII+wTln1*RT1->Cool_Radr_HeIII) * Y->HeIII * Rate->Radr_HeIII);
  fprintf(fp,"Collisional Electron Dissociation  H2    %e\n",
    Y->e * en_B * 
	  cl->R.Cool_Coll_H2 * Y->H2 * Rate->Coll_e_H2); 
  fprintf(fp,"Collisional H Dissociation  H2    %e\n",  
    Y->HI * en_B *
	  cl->R.Cool_Coll_H2 * Y->H2 * Rate->Coll_HI_H2);
  fprintf(fp,"Collisional H2 Dissociation  H2    %e\n",  
    Y->H2 * en_B *
	  cl->R.Cool_Coll_H2 * Y->H2 * Rate->Coll_H2_H2);
   fprintf(fp,"Collisional HII Dissociation  H2    %e\n",  
    Y->HII * en_B *
          cl->R.Cool_Coll_H2 * Y->H2 * Rate->Coll_HII_H2);
  fprintf(fp,"Collisional Ionization  HI    %e\n",
    Y->e * en_B * 
	  cl->R.Cool_Coll_HI * Y->HI * Rate->Coll_HI);
  fprintf(fp,"Collisional Ionization  HeI   %e\n",
    Y->e * en_B * 
    cl->R.Cool_Coll_HeI * Y->HeI * Rate->Coll_HeI);
  fprintf(fp,"Collisional Ionization  HeII  %e\n",
    Y->e * en_B * 
    cl->R.Cool_Coll_HeII * Y->HeII * Rate->Coll_HeII);
  fprintf(fp,"Dielectric Recombination HeII %e\n",
    Y->e * en_B * 
    cl->R.Cool_Diel_HeII * Y->HeII * Rate->Diel_HeII);
  fprintf(fp,"Line cooling H2   %e\n",
    Y->e * en_B * 
    (wTln0*RT0->Cool_Line_H2+wTln1*RT1->Cool_Line_H2) * Y->H2);
  fprintf(fp,"Line cooling HI   %e\n",
    Y->e * en_B * 
	  (wTln0*RT0->Cool_Line_HI+wTln1*RT1->Cool_Line_HI) * Y->HI);
  fprintf(fp,"Line cooling HeI  %e\n",
    Y->e * en_B * 
    (wTln0*RT0->Cool_Line_HeI+wTln1*RT1->Cool_Line_HeI) * Y->HeI);
  fprintf(fp,"Line cooling HeII  %e\n",
    Y->e * en_B * 
    (wTln0*RT0->Cool_Line_HeII+wTln1*RT1->Cool_Line_HeII) * Y->HeII );
  fprintf(fp,"Low T cooling (Z=0.001)  %e\n",
    en_B * 
    (wTln0*RT0->Cool_LowT+wTln1*RT1->Cool_LowT)*cl->R.Cool_LowTFactor*0.001);
  fprintf(fp,"Net Metal Heating %e \n", Rate->Heat_Metal); 
  fprintf(fp,"Net Metal Cooling %e \n", Rate->Cool_Metal);*/

  fprintf(fp,"Cool Brem HII    %e\n", Y->e * en_B * (clCoolBrem1(Rate->T) * Y->HII));
  fprintf(fp,"Cool Brem HeII   %e\n", Y->e * en_B * (clCoolBrem1(Rate->T) * Y->HeII));
  fprintf(fp,"Cool Brem HeIII  %e\n",Y->e * en_B * (clCoolBrem2(Rate->T) * Y->HeIII));
  fprintf(fp,"Radiative Recombination  HII    %e\n",Y->e * en_B * clCoolRadrHII(Rate->T) * Y->HII * Rate->Radr_HII );
  fprintf(fp,"Radiative Recombination  HeII   %e\n",Y->e * en_B * clCoolRadrHeII(Rate->T) * Y->HeII * Rate->Radr_HeII);
  fprintf(fp,"Radiative Recombination  HeIII  %e\n",Y->e * en_B * clCoolRadrHeIII(Rate->T) * Y->HeIII * Rate->Radr_HeIII);
  fprintf(fp,"Collisional Electron Dissociation  H2    %e %e\n", Y->e * en_B * cl->R.Cool_Coll_H2 * Y->H2 * Rate->Coll_e_H2, Y->e * en_B * cl->R.Cool_Coll_H2 * Y->H2 * Rate->Coll_e_H2);/** s_self * s_dust);  removing shielding of collisions 9/22/10*/
  fprintf(fp,"Collisional H Dissociation  H2    %e %e\n",  Y->HI * en_B *cl->R.Cool_Coll_H2 * Y->H2 * Rate->Coll_HI_H2, Y->e * en_B * cl->R.Cool_Coll_H2 * Y->H2 * Rate->Coll_HI_H2); /** s_self * s_dust); removing shielding of collisions 9/22/10*/
  fprintf(fp,"Collisional H2 Dissociation  H2    %e %e\n",  Y->H2 * en_B * cl->R.Cool_Coll_H2 * Y->H2 * Rate->Coll_H2_H2, Y->e * en_B * cl->R.Cool_Coll_H2 * Y->H2);/* * Rate->Coll_H2_H2* s_self * s_dust); removing shielding of collisions 9/22/10*/
  fprintf(fp,"Collisional HII Dissociation  H2    %e %e\n",  Y->HII * en_B * cl->R.Cool_Coll_H2 * Y->H2 * Rate->Coll_HII_H2, Y->e * en_B * cl->R.Cool_Coll_H2 * Y->H2 * Rate->Coll_HII_H2);/** s_self * s_dust); removing shielding of collisions 9/22/10*/
  fprintf(fp,"Collisional Ionization  HI    %e\n", Y->e * en_B *  cl->R.Cool_Coll_HI * Y->HI * Rate->Coll_HI);
  fprintf(fp,"Collisional Ionization  HeI   %e\n", Y->e * en_B *  cl->R.Cool_Coll_HeI * Y->HeI * Rate->Coll_HeI);
  fprintf(fp,"Collisional Ionization  HeII  %e\n", Y->e * en_B *  cl->R.Cool_Coll_HeII * Y->HeII * Rate->Coll_HeII);
  fprintf(fp,"Dielectric Recombination HeII %e\n", Y->e * en_B * cl->R.Cool_Diel_HeII * Y->HeII * Rate->Diel_HeII);
  /*fprintf(fp,"Line cooling H2   %e\n",Y->e * en_B * clCoolLineH2(Rate->T, Y->HI + Y->HII + 2.0*Y->H2) * Y->H2);*/
  fprintf(fp,"Line cooling H2 (H2)   %e\n",clCoolLineH2_H2(Rate->T) * Y->H2 * Y->H2); /*Glover and Abel 08, reaction 10*/
  fprintf(fp,"Line cooling H2 (HI)   %e\n",clCoolLineH2_HI(Rate->T) * Y->H2 * Y->HI); /*Glover and Abel 08, reaction 9*/
  fprintf(fp,"Line cooling H2 (HeI)   %e\n",clCoolLineH2_He(Rate->T) * Y->H2 * Y->HeI); /*Glover and Abel 08, reaction 11*/
  fprintf(fp,"Line cooling H2 (HII)   %e\n",clCoolLineH2_HII(Rate->T) * Y->H2 * Y->HII);
  fprintf(fp,"Line cooling H2 (e)   %e\n",clCoolLineH2_e(Rate->T) * Y->H2 * Y->e); /*Glover and Abel 08, reaction 8*/
  fprintf(fp,"Line cooling HI   %e\n", Y->e * en_B * clCoolLineHI(Rate->T) * Y->HI);
  fprintf(fp,"Line cooling HeI  %e\n",Y->e * en_B * clCoolLineHeI(Rate->T) * Y->HeI);
  fprintf(fp,"Line cooling HeII  %e\n",Y->e * en_B * clCoolLineHeII(Rate->T) * Y->HeII);
  fprintf(fp,"Low T cooling  %e\n",en_B * 
    (wTln0*RT0->Cool_LowT+wTln1*RT1->Cool_LowT)*cl->R.Cool_LowTFactor*0.001);
  fprintf(fp,"Photon Heat HI %e %e\n",Y->HI*cl->R.Heat_Phot_HI*Rate->Phot_HI, Y->HI*cl->R.Heat_Phot_HI*Rate->Phot_HI*s_dust);
  fprintf(fp,"Photon Heat HeI %e\n",Y->HeI   * cl->R.Heat_Phot_HeI   * Rate->Phot_HeI);
  fprintf(fp,"Photon Heat HeII %e \n",Y->HeII   * cl->R.Heat_Phot_HeII   * Rate->Phot_HeII);
  fprintf(fp,"Photon Heat H2 %e %e\n",Y->H2*cl->R.Heat_Phot_H2*Rate->Phot_H2, Y->H2*cl->R.Heat_Phot_H2*Rate->Phot_H2*s_self *s_dust);
  fprintf(fp,"Net Metal Heating %e \n", Rate->Heat_Metal); 
  fprintf(fp,"Net Metal Cooling %e \n", Rate->Cool_Metal);
  printf("Radr: %#2e\n",-1.0*ne*(clCoolRadrHII(Rate->T) * Y->HII * Rate->Radr_HII + clCoolRadrHeII(Rate->T) * Y->HeII * Rate->Radr_HeII + clCoolRadrHeIII(Rate->T) * Y->HeIII * Rate->Radr_HeIII));
  printf("Brems: %#2e\n",-1.0*ne*(clCoolBrem1(Rate->T) * ( Y->HII + Y->HeII ) + clCoolBrem2(Rate->T) * Y->HeIII));
  printf("Photon Heat HI %e %e \n",Y->HI*cl->R.Heat_Phot_HI*Rate->Phot_HI, Y->HI*cl->R.Heat_Phot_HI*Rate->Phot_HI*s_dust);
  printf("Photon Heat HeI %e   \n",Y->HeI*cl->R.Heat_Phot_HeI*Rate->Phot_HeI);
  printf("Photon Heat HeII %e  \n",Y->HeII*cl->R.Heat_Phot_HeII*Rate->Phot_HeII);
  printf("Photon Heat H2 %e %e \n",Y->H2*cl->R.Heat_Phot_H2*Rate->Phot_H2, Y->H2*cl->R.Heat_Phot_H2*Rate->Phot_H2*s_self*s_dust);
  printf("Net Metal Heating %e \n", Rate->Heat_Metal); 
  printf("Net Metal Cooling %e \n", -1.0*Rate->Cool_Metal);

}

void clAbunds( COOL *cl, PERBARYON *Y, RATE *Rate, double rho, double ZMetal) {
  double en_B =rho*CL_B_gm;

  double s_dust, s_self, Rate_Phot_HI; 

  /*Coll. dissos./Rad. Recomb*/
  double rcirrHeI  = (Rate->Coll_HeI)/(Rate->Totr_HeII);
  double rcirrHeII = (Rate->Coll_HeII)/(Rate->Radr_HeIII);

  /*Photon dissos./Rad. Recomb*/
  double rpirrHeI  = (Rate->Phot_HeI)/(Rate->Totr_HeII * en_B);
  double rpirrHeII = (Rate->Phot_HeII)/(Rate->Radr_HeIII * en_B);

  double yH;
  double yH2 = 0; 
  double yHI = 0; 
  double yHII = 0; 

  double yHe; 
  double yHeI = 0.;
  double yHeII = 0.;
  double yHeIII = 0;

  double yeMax;
  double rye,ye;
  double fHI,fHII,fHeI,fHeII,rfH,rfHe,yHI_old,yHeII_old,yHII_old; 
  double yH2_old;
  int i;  

  clSetAbundanceTotals(cl,ZMetal,&yH,&yHe,&yeMax);

  s_dust = clDustShield(yH*en_B, 0, ZMetal, Rate->CorreLength);
  s_self = 1.0; /* set to zero so that when setting fHI, I will not be dividing by zero */

  Rate_Phot_HI = Rate->Phot_HI;

  for ( i=0 ; i<MAXABUNDITERATIONS ; i++ ) {
    yHI_old   = yHI;
    yHeII_old = yHeII;
    yH2_old   = yH2;
    yHII_old = yH - yHI - 2.0*yH2;

    ye = (yeMax-(yHI + 2 * yH2 + 2 * yHeI + yHeII)); /*Free electrons*/
    if (ye <= 0) {
      ye = 0;
      yHII = 0;
      yHeI = yHe;
      yHeII = 0;
      yHeIII = 0;
      if (s_dust == 0 || s_self == 0) {
	yHI = 0;
	yHII = 0;
	yH2 = yH/2.0;
      }
      else {
	fHI = 2.0*(Rate->DustForm_H2*en_B*(yHI_old + 2.0*yH2_old))/
	           (Rate->Coll_HI_H2*en_B*yHI_old +
	           Rate->Coll_H2_H2*en_B*yH2_old +
		     Rate->Coll_HII_H2*en_B*yHII_old) +
		   (Rate->Phot_H2)*s_dust*s_self;
	yHI = yH/(1 + fHI);
	yH2 = (yH - yHI)/2.0;
      }
      if ( fabs(yHeII_old-yHeII) < EPS * yHeII && fabs(yHI_old-yHI) < EPS * yHI ) break;
      s_dust = clDustShield(yHI*en_B, yH2*en_B, ZMetal, Rate->CorreLength);
      s_self = clSelfShield(yH2*en_B, Rate->CorreLength);
      if (cl->bShieldHI) Rate_Phot_HI = Rate->Phot_HI*s_dust;
    }
    else {
      rye = 1/ye;

      if (s_dust == 0 || s_self == 0) {
	yHI = 0;
	yHII = 0;
	yH2 = yH/2.0;
      }
      else {
	fHII = (Rate->Radr_HII*en_B*ye)/
	       (Rate_Phot_HI*s_dust + 
	        Rate->Coll_HI*en_B*ye);  
	fHI = 2.0*(Rate->DustForm_H2*en_B*(yHI_old + 2.0*yH2_old))/
	         (Rate->Coll_e_H2*en_B*ye + 
	           Rate->Coll_HI_H2*en_B*yHI_old + 
	           Rate->Coll_H2_H2*en_B*yH2_old +
		  Rate->Coll_HII_H2*en_B*yHII_old +  
		  (Rate->Phot_H2)*s_dust*s_self);
	rfH  =  1 / ( 1 + fHII * (1 + fHI) ); 
	yHII =  yH * rfH;
	yHI  =  yH * rfH * fHII;
	yH2  = (yH - yHI - yHII)/2.0;
      }
      s_dust = clDustShield(yHI*en_B, yH2*en_B, ZMetal, Rate->CorreLength);
      s_self = clSelfShield(yH2*en_B, Rate->CorreLength);
      if (cl->bShieldHI) Rate_Phot_HI = Rate->Phot_HI*s_dust;
      fHeI  = rcirrHeI + rpirrHeI * rye;
      fHeII = rcirrHeII + rpirrHeII * rye;
      rfHe  = 1 / ( 1 + fHeI * (1 + fHeII) );
      yHeI  = yHe * rfHe;
      yHeII = yHe * fHeI * rfHe;
      yHeIII = yHe / ((1.0/fHeI+1.0)/fHeII+1.0);

      if ( fabs(yHeII_old-yHeII) < EPS * yHeII && fabs(yHI_old-yHI) < EPS * yHI ) break;
    }
  }

  Y->e = ye;
  Y->HI = yHI;
  Y->HII = yHII;
  Y->H2 = 2.0*yH2;
  Y->HeI = yHeI;
  Y->HeII = yHeII;
  Y->HeIII = yHeIII;
  Y->Total = Y->e + yH + yHe + ZMetal/MU_METAL - Y->H2;/*Don't want to double count hydrogen atoms in molecules when finding the total number of particles */
#ifdef COOLDEBUG
  if (0) {
    FILE *fpdebug;
    fpdebug = fopen("lte.txt","a");
    fprintf(fpdebug,"%#2f, %#2e, %g, %g, %g, %g, %g, %g, %g, %#2e, %#2e, %#2e, %#2e, %#2e, %#2e, %#2e, %#2e, %#2e, %#2e, %#2e, %#2e, %#2e, %#2e, %#2e, %#2e, %#2e, %#2e, %#2e\n",
	    en_B, ye, Y->H2, Y->HI, Y->HII, Y->HeI, Y->HeII, Y->HeIII, ZMetal, Rate->CorreLength, s_dust, s_self,
	    Rate->Coll_HI,
	    Rate->Coll_HeI,
	    Rate->Coll_HeII,
	    Rate->Coll_e_H2,
	    Rate->Coll_HI_H2,
	    Rate->Coll_H2_H2,
	    Rate->Coll_HII_H2,
	    Rate->Radr_HII,
	    Rate->Totr_HeII,
	    Rate->Radr_HeII,
	    Rate->Radr_HeIII,
	    Rate->DustForm_H2,
	    Rate->Phot_HI,
	    Rate->Phot_HeI,
	    Rate->Phot_HeII,
	    Rate->Phot_H2);
    fclose(fpdebug);
  }
#endif
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
    double Y_H, Y_He, Y_eMax;
    clSetAbundanceTotals(cl,0.0,&Y_H,&Y_He,&Y_eMax);/* no metals */
    return clTemperature( 2*Y_H - Y_HI + 3*Y_He - 2*Y_HeI - Y_HeII,  E );
    }

double clSelfShield (double yH2, double h) {
  double x, column_denH2, omega_H2 = 0.2; 
      if (yH2 <= 0) return 1.0;
      column_denH2 = h*yH2;
      x = column_denH2/5e14;
      return (1 - omega_H2)/(1 + x)/(1 + x) + omega_H2/sqrt(1 + x)*exp(-0.00085*sqrt(1 + x));
}

double clDustShield (double yHI, double yH2, double z, double h) {
  double column_denHI, column_denH2, sigmad = 2e-21; /*4e-21;*/
      if (yHI < 0) column_denHI = 0;
      else column_denHI = h*yHI;
      if (yH2 < 0) column_denH2 = 0;
      else column_denH2 = h*yH2;
      return exp(-1.0*sigmad*z/ZSOLAR*(column_denHI + 2.0*column_denH2));
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

/*H2 Collision*/
/*Lepp & Shull, 1983*/
/* Should I update this to reaction 10 from Glober and Abel 08 (Martin, Keogh & Mandy (1998); Shapiro & Kang (1987))?  CC */
double clRateColl_H2_H2( double T){
  double ratecoll_H2_H2;
  if (T < 7291) ratecoll_H2_H2 = 5.22e-14*exp(-3.22e4/T);
  else ratecoll_H2_H2 = 3.17e-15*exp(-4060./T - (7500./T)*(7500./T));
  return ratecoll_H2_H2;
}

/*Ionized H*/
/*Abel 1997, k11*/
/* Should I update this to reaction 8 from Glober and Abel 08 (Trevisan & Tennyson (2002a))?  CC*/
double clRateColl_HII_H2(double T){
  double ratecoll_HII_H2, TL = log(T*CL_eV_per_K);
  ratecoll_HII_H2 = exp(-24.24914687731536
		   + 3.400824447095291*TL
		   - 3.898003964650152*pow(TL,2)
		   + 2.045587822403071*pow(TL,3)
		   - 0.5416182856220388*pow(TL,4)
		   + 0.0841077503763412*pow(TL,5)
		   - 0.007879026154483455*pow(TL,6)
		   + 0.0004138398421504563*pow(TL,7)
		   - 9.36345888928611e-6*pow(TL,8));
  return ratecoll_HII_H2;
}

/*Electron Collision*/
/*Donahue & Shull 1991, Abel 1997, k12*/
/* Should I update this to reaction 8 from Glober and Abel 08 (Trevisan & Tennyson (2002a))? CC*/
double clRateColl_e_H2( double T){
  double ratecoll_e_H2;
  ratecoll_e_H2 = 5.6e-11*pow(T,0.5)*exp(-102124.0/T); /* Temperature listed in kelvin */ 
  return ratecoll_e_H2;
}

/*Neutral H*/
/*Dove and Mandy 1986, Abel 1997, k13*/
/* Should I update this to reaction 9 from Glober and Abel 08 (Mac Low & Shull (1986))?*/
double clRateColl_HI_H2(double T){
  double ratecoll_HI_H2;
  /*Donahue & Shull  ratecollH2 = 6.11e-14*exp(-4.48*CL_eV_per_K/T);*/
  ratecoll_HI_H2 = 1.067e-10*pow(CL_eV_per_K*T,2.012)*exp(-1.0*(4.463/T/CL_eV_per_K)*pow(1+0.2472*CL_eV_per_K*T,3.512));
  return ratecoll_HI_H2;
}

/*Helium collisions*/
/*Reaction 11 from Glober and Abel 08 (Dove et al. (1987)) -- Update code to include?*/
/*double clRateColl_He_H2(double T){
  double ratecollH2;
  ratecollH2 = log10(-27.029 + 3.801*log(T) - 29487.0/T);
  return ratecollH2;
  }*/

/*Abel 1997, k14*/
double clRateColl_Hm_e(double T){
  double ratecoll_Hm_e, LT = log(T*CL_eV_per_K);
  ratecoll_Hm_e = exp(-18.01849334273
		      + 2.360852208681*LT
		      - 0.2827443061704*pow(LT,2)
		      + 0.01623316639567*pow(LT,3)
		      - 0.03365012031362999*pow(LT,4)
		      + 0.01178329782711*pow(LT,5)
		      - 0.001656194699504*pow(LT,6)
		      + 0.0001068275202678*pow(LT,7)
		      - 2.631285809207e-6*pow(LT,8));
  return ratecoll_Hm_e;
}

/*Abel 1997, k16*/
double clRateColl_Hm_HII(double T){
  double ratecoll_Hm_HII = 7.0e-8/sqrt(T/100.0);
  return ratecoll_Hm_HII;
}

/*-----------------------------------------------------------------
 *     More Formation Paths for H2 in gas phase
 *-----------------------------------------------------------------*/
/*     H + e- = H- + gamma, Abel 1997, k7*/
double clRateHI_e(double T){
  double rateHI_e, Tlog10 = log10(T), temp = 4.0415e-5*pow(Tlog10,6) - 5.447e-3*pow(Tlog10,4);
  if (T < 6e7) rateHI_e = 1.429e-18*pow(T,0.762)*pow(T,0.1523*Tlog10)*pow(T,-3.274e-2*Tlog10*Tlog10);
  else rateHI_e = 3.802e-17*pow(T,0.1998*Tlog10)*pow(10.0,temp);
  return rateHI_e;
}


/*     H + H- = H2 + e-, Abel 1997, k8*/
double clRateHI_Hm(double T){
  double rateHI_Hm, Tev = T*CL_eV_per_K, TL = log(T*CL_eV_per_K);
  if (Tev > 0.1)
    rateHI_Hm = exp(-20.06913897587003
		   + 0.2289800603272916*TL
		   + 0.03599837721023835*pow(TL,2)
		   - 0.004555120027032095*pow(TL,3)
		   - 0.0003105115447124016*pow(TL,4)
		   + 0.0001073294010367247*pow(TL,5)
		   - 8.36671960467864e-6*pow(TL,6)
		   + 2.238306228891639e-7*pow(TL,7));
  else rateHI_Hm = 1.428e-9;
  return rateHI_Hm;
}

/*-----------------------------------------------------------------
 *     Radiative Recombination rates
 *-----------------------------------------------------------------*/
/*     H+ + e- -> H + gam  Verner & Ferland 1996 */
double clRateRadrHII( double T ) {
  double Tsq = sqrt(T); 

   return 7.982e-11/( Tsq*0.563615 *     
   pow(1+Tsq*0.563615,0.252) * pow(1+Tsq*1.192167e-3,1.748));  
  /*  double lgT = log10(T); 
  double lg_rec; 
  lg_rec = -9.77199 -1.20330  * lgT + 0.891480 *pow(lgT, 2.0) -0.677734 * pow(lgT, 3.0) + 0.290476 * pow(lgT, 4.0) - 0.0746242 * pow(lgT, 5.0) + 0.0115600 * pow(lgT, 6.0) -0.00105994 * pow(lgT, 7.0) + 5.30408e-05 * pow(lgT, 8.0)  -1.11692e-06 * pow(lgT, 9.0);
  return exp10(lg_rec);  */
}

/*     He+ + e- -> He + gam  radiative  Verner & Ferland 1996 */
double clRateRadrHeII( double T ) {
  /* 
   * Note that these functions do not meet perfectly at 1e6 -- 2% difference 
   * The derivatives are different there also: So the apparent error is large 
   */
  /*  double Tsq = sqrt(T);
  if (T < 1e6)
    return  3.294e-11/( Tsq*0.253673 *
             pow(1+Tsq*0.253673,0.309) * pow(1+Tsq*1.649348e-4,1.691));
  else
    return  9.356e-10/( Tsq*4.841607 *
    pow(1+Tsq*4.841607,0.2108) * pow(1+Tsq*4.628935e-4,1.7892)); */

  /* fix the the 2% difference in at 1e6 K solving the bump problem in the cooling curve.*/
  /* method: take the data from verner & ferland (1996) and fit it with 9th order polynomial */  
  double lgT = log10(T); 
  double lg_rec;
  lg_rec = -10.3217 + 0.749031* lgT - 1.66439*pow(lgT, 2.0) + 1.13114 * pow(lgT, 3.0) -0.459420 * pow(lgT, 4.0) + 0.114337 * pow(lgT, 5.0) -0.0174612 * pow(lgT, 6.0) + 0.00158680 * pow(lgT, 7.0) -7.86098e-05 * pow(lgT, 8.0) +  1.63398e-06 * pow(lgT, 9.0); 
  return pow(10.0, lg_rec); 
  /* fit the curve from data generated from cloudy */  
  
  /* double lgT = log10(T); 
  double lg_rec; 

  lg_rec = - 8.63240 -2.55916 * lgT + 2.36054 *pow(lgT, 2.0) -1.48299 * pow(lgT, 3.0) + 0.541563 * pow(lgT, 4.0) -0.120970 * pow(lgT, 5.0) + 0.0166637 * pow(lgT, 6.0) -0.00138250 * pow(lgT, 7.0) + 6.33987e-05 * pow(lgT, 8.0) - 1.23508e-06  * pow(lgT, 9.0);
  return exp10(lg_rec);  */
}

/*     He+ + e- -> He + gam  dielectronic  Aldovandi&Pequignot 1973 (Black 1981) */
double clRateDielHeII( double T ) {
  double T_inv = 1.0/T,arg;

  arg = -4.7e5*T_inv;
  if (arg < CL_MAX_NEG_EXP_ARG) return 0;
  return 1.9e-3*pow(T,-1.5)*exp(arg)*(1+0.3*exp(-9.4e4*T_inv));
}

#define CHTR_a 7.47e-15 

double clRateChtrHeII(double T) {
  /*double T_4 = T/1e4; 
    if (T_4 < 0.6) 
    return CHTR_a * pow(0.6, 2.06)*(1. + 9.93* exp(-3.89 *0.6));
  else if (T_4 >= 10.) 
    return CHTR_a * pow(10., 2.06)*(1. + 9.93* exp(-3.89 *10.));
  else
  return CHTR_a * pow(T_4, 2.06)*(1. + 9.93* exp(-3.89 *T_4));  */
  return 0. ; 
}

/*     He++ + e- -> He+ + gam  Verner & Ferland 1996 */
double clRateRadrHeIII( double T ) {
  double Tsq = sqrt(T);

  return 1.891e-10/( Tsq*0.326686 *
          pow(1+Tsq*0.326686,0.2476) * pow(1+Tsq*6.004084e-4,1.7524));
}

double clRateDustFormH2( double z, double clump ) {
  double Rate_dust = 0;
  /* clump = 10.0; Ranges from 2-10 to 30-100, Gendin et al 2008 CC*/ 
  Rate_dust = 3.5e-17*z/ZSOLAR*clump; /*Formation rate coefficient of molecular hydrogen on dust, Gnedin et al 2008, Wolfire 2008, unit of cc per s, divide metallicity by solar metallicity (was 0.0177 in first iterations of code) CC*/  
  return Rate_dust;
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

#ifdef CLOUDY 
#define CL_aHI   -1275.29
#define CL_bHI   374.535
#define CL_cHI   -42.8662
#define CL_dHI   2.19186
#define CL_fHI   -0.0422634
#define CL_bHI_old  1.18348e05

double clCoolLineHI( double T ) {
  double lnT, ln_cool, arg;
  arg = -CL_bHI_old/T;
  lnT = log(T); 
  if (arg < CL_MAX_NEG_EXP_ARG) return 0;
  ln_cool = CL_aHI +  CL_bHI*lnT +  CL_cHI * pow(lnT, 2.0) +  CL_dHI * pow(lnT, 3.0) +  CL_fHI*pow(lnT, 4.0); 
  return CL_B_gm*exp(ln_cool);
}

#else 
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
#endif 

#ifdef CLOUDY
#define CL_aHeI   -3133.67
#define CL_bHeI   980.557
#define CL_cHeI   -117.344
#define CL_dHeI   6.26898
#define CL_fHeI   -0.125976
#define CL_bHeI_old  1.3179e04
double clCoolLineHeI( double T ) {
  double lnT, ln_cool, arg;
  arg = -CL_bHeI_old/T;
  lnT = log(T); 
  if (arg < CL_MAX_NEG_EXP_ARG) return 0;
  ln_cool = CL_aHeI +  CL_bHeI*lnT +  CL_cHeI * pow(lnT, 2.0) +  CL_dHeI * pow(lnT, 3.0) +  CL_fHeI*pow(lnT, 4.0); 
  return CL_B_gm*exp(ln_cool);
}
#else
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
#endif

#ifdef CLOUDY 

#define CL_aHeII  -4002.63
#define CL_bHeII  1197.95 
#define CL_cHeII  -136.703 
#define CL_dHeII  6.96923 
#define CL_fHeII  -0.133835
#define CL_bHeII_old  4.73638e05

 double clCoolLineHeII( double T ) {
  double lnT, ln_cool, arg;
  arg = -CL_bHeII_old/T;
  lnT = log(T); 
  if (arg < CL_MAX_NEG_EXP_ARG) return 0;
  ln_cool = CL_aHeII +  CL_bHeII*lnT +  CL_cHeII * pow(lnT, 2.0) +  CL_dHeII * pow(lnT, 3.0) +  CL_fHeII*pow(lnT, 4.0); 
  return CL_B_gm*exp(ln_cool);
}    

#else 
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
#endif

double clCoolLineH2_HI( double T){ /* H2-H collisionally-induced cooling, Glover & Abel 08, Table 8 */
  double a00 = -16.818342,
    a10 = 37.383713,
    a20 = 58.145166,
    a30 = 48.656103,
    a40 = 20.159831,
    a50 = 3.8479610;
  double a01 = -24.311209,
    a11 = 3.5692468,
    a21 = -11.332860,
    a31 = -27.850082,
    a41 = -21.328264,
    a51 = -4.2519023;
  double a02 = -24.311209,
    a12 = 4.6450521,
    a22 = -3.7209846,
    a32 = 5.9369081,
    a42 = -5.5108047,
    a52 = 1.5538288;
  double xint = 6000, slope = 2.10095, yint = 1.86368e-22;

  if (T <= 100) return pow(10.0, a00 + 
                                 a10*log10(T/1000.0) + 
                                 a20*log10(T/1000.0)*log10(T/1000.0) + 
                                 a30*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0) + 
                                 a40*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0) + 
                                 a50*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0));
  else if (T <= 1000) return pow(10.0, a01 + 
                                 a11*log10(T/1000.0) + 
                                 a21*log10(T/1000.0)*log10(T/1000.0) + 
                                 a31*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0) + 
                                 a41*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0) + 
                                 a51*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0));
  else if (T <= xint) return pow(10.0, a02 + 
                                 a12*log10(T/1000.0) + 
                                 a22*log10(T/1000.0)*log10(T/1000.0) + 
                                 a32*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0) + 
                                 a42*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0) + 
                                 a52*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)) ;
  else return pow(10.0,slope*(log10(T/1000.0) - log10(xint/1000.0)) + log10(yint));
}

double clCoolLineH2_H2( double T){  /* H2-H2 collisionally-induced cooling, Glover & Abel 08, Table 8*/
  double a0 = -23.962112,
    a1 = 2.09433740,
    a2 = -0.77151436,
    a3 = 0.43693353,
    a4 = -0.14913216,
    a5 = -0.033638326;
 double xint = 6000, slope = 1.34460, yint = 2.19802e-23;
  
 if (T <= xint) return pow(10.0, a0 + 
	     a1*log10(T/1000.0) + 
	     a2*log10(T/1000.0)*log10(T/1000.0) + 
	     a3*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0) + 
	     a4*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0) + 
	     a5*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0));
 else return pow(10.0,slope*(log10(T/1000.0) - log10(xint/1000.0)) + log10(yint));
}

double clCoolLineH2_He( double T){ /* H2-He collisionally-induced cooling, Glover & Abel 08, Table 8 */
  double a0 = -23.689237,
    a1 = 2.1892372,
    a2 = -0.81520438,
    a3 = 0.29036281,
    a4 = -0.16596184,
    a5 = 0.19191375;
  double xint = 6000, slope = 1.48703, yint = 4.48145e-23;
  
  if (T <= xint) return pow(10.0, a0 + 
	     a1*log10(T/1000.0) + 
	     a2*log10(T/1000.0)*log10(T/1000.0) + 
	     a3*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0) + 
	     a4*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0) + 
	     a5*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0));
  else return pow(10.0,slope*(log10(T/1000.0) - log10(xint/1000.0)) + log10(yint));
}

double clCoolLineH2_HII( double T){ /*H2-HII collisionally-induced cooling , Glover & Abel 08, Table 8 */
  double a0 = -21.716699,
    a1 = 1.3865783,
    a2 = -0.37915285,
    a3 = 0.11453688,
    a4 = -0.23214154,
    a5 = 0.058538864;
  double xint = 10000, slope = 0.336011, yint = 1.70474e-21;
  
  if (T <= xint) return pow(10.0, a0 + 
	     a1*log10(T/1000.0) + 
	     a2*log10(T/1000.0)*log10(T/1000.0) + 
	     a3*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0) + 
	     a4*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0) + 
	     a5*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0));
  else return pow(10.0,slope*(log10(T/1000.0) - log10(xint/1000.0)) + log10(yint));
}

double clCoolLineH2_e( double T){  /*Cooling based on H2-e collisionally-induced dissociation, Glover & Abel 08, Table 8*/
  double a00 = -34.286155,
    a10 = -48.537163,
    a20 = -77.121176,
    a30 = -51.352459,
    a40 = -15.169160,
    a50 = -0.98120322;
  double a01 = -22.190316,
    a11 = 1.5728955,
    a21 = -0.21335100,
    a31 = 0.96149759,
    a41 = -0.91023495,
    a51 = 0.13749749;
 double xint = 10000, slope = 1.07723, yint = 2.28029e-21;

  if (T <= 200) return pow(10.0, a00 + 
                                 a10*log10(T/1000.0) + 
                                 a20*log10(T/1000.0)*log10(T/1000.0) + 
                                 a30*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0) + 
                                 a40*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0) + 
                                 a50*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0));
  else if (T <= xint) return pow(10.0, a01 + 
                                 a11*log10(T/1000.0) + 
                                 a21*log10(T/1000.0)*log10(T/1000.0) + 
                                 a31*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0) + 
                                 a41*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0) + 
                                 a51*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0)*log10(T/1000.0));
  else return pow(10.0,slope*(log10(T/1000.0) - log10(xint/1000.0)) + log10(yint));
}

/*Replaced by data from Glover & Abel 08
Cooling from rot-vib transitions of H2 CC -- Martin, Schwarz & Mandy, 1996
Note that this only includes excitations from H2-H collisions and may be a problem when most of the gas has been turned to H2*/
/*
double clCoolLineH2( double T, double YH){
  if (T > 45000 || YH < 1e-3) return 0;

  double alpha1 =-1.058656e3,
    alpha2 = 6.412920e2,
    alpha3 =-1.330667e2,
    alpha4 = 9.285717,
    alpha5 = 9.160106e1,
    alpha6 = 2.680075e3,
    alpha7 = 9.500433e3,
    alpha8 =-1.253746e3,
    alpha9 = 7.792906e2,
    alpha10=-1.628687e2,
    alpha11= 1.145088e1,
    alpha12= 1.057438e2,
    alpha13= 2.889762e3,
    alpha14= 1.382858e4,
    alpha15=-1.175497e2,
    alpha16= 7.886144e1,
    alpha17=-1.777082e1,
    alpha18= 1.338843,
    alpha19= 3.406424e3;
  double logt, logyh, lograd, a, b, f, w;

  logt=log10(T);
  logyh=log10(YH);

  a = -1.0*logyh + alpha1 + alpha2*logt + alpha3*logt*logt + alpha4*logt*logt*logt + alpha5*log10(1.0 + alpha6/T) + alpha7/T;
  b = alpha8 + alpha9*logt + alpha10*logt*logt + alpha11*logt*logt*logt + alpha12*log10(1.0 + alpha13/T) + alpha14/T;
  f = logyh - b; 
  w = alpha15 + alpha16*logt + alpha17*logt*logt + alpha18*logt*logt*logt + alpha19/T;
  lograd = a + 0.5*(f - sqrt(f*f+2.*w*w));
  return pow(10.0,lograd);
  }*/

/*Replaced by data from Glover & Abel 08*/
/*Cooling from rot-vib transitions or H2 CC -- Martin, Schwarz & Mandy, 1996*/
/*This function should rely on the density so it will return incorrect results, use clCoolLineH2 if possible*/
/*
double clCoolLineH2_table( double T){
  if (T > 45000) return 0;

  double alpha1 =-1.058656e3,
    alpha2 = 6.412920e2,
    alpha3 =-1.330667e2,
    alpha4 = 9.285717,
    alpha5 = 9.160106e1,
    alpha6 = 2.680075e3,
    alpha7 = 9.500433e3,
    alpha8 =-1.253746e3,
    alpha9 = 7.792906e2,
    alpha10=-1.628687e2,
    alpha11= 1.145088e1,
    alpha12= 1.057438e2,
    alpha13= 2.889762e3,
    alpha14= 1.382858e4,
    alpha15=-1.175497e2,
    alpha16= 7.886144e1,
    alpha17=-1.777082e1,
    alpha18= 1.338843,
    alpha19= 3.406424e3;
  double logt, logyh, lograd, a, b, f, w;

  logt=log10(T);
  logyh=0; Assumed Density

  a = -1.0*logyh + alpha1 + alpha2*logt + alpha3*logt*logt + alpha4*logt*logt*logt + alpha5*log10(1.0 + alpha6/T) + alpha7/T;
  b = alpha8 + alpha9*logt + alpha10*logt*logt + alpha11*logt*logt*logt + alpha12*log10(1.0 + alpha13/T) + alpha14/T;
  f = logyh - b; 
  w = alpha15 + alpha16*logt + alpha17*logt*logt + alpha18*logt*logt*logt + alpha19/T;
  lograd = a + 0.5*(f - sqrt(f*f+2.*w*w));
  return pow(10.0,lograd);
}*/

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
double clEdotInstant_Table( COOL *cl, PERBARYON *Y, RATE *Rate, double rho, 
			    double ZMetal, double *dEdotHeat, double *dEdotCool )
{
  double en_B = rho*CL_B_gm;
  double xTln,wTln0,wTln1;/*,wTln0d,wTln1d;*/
  RATES_T *RT0,*RT1;/*,*RT0d,*RT1d;*/
  int iTln;

  double ne,LowTCool;
  double s_dust, s_self, Rate_Phot_HI;
  s_dust = clDustShield(Y->HI*en_B, Y->H2*en_B, ZMetal, Rate->CorreLength);
  s_self = clSelfShield(Y->H2*en_B, Rate->CorreLength);
  if (cl->bShieldHI) Rate_Phot_HI = Rate->Phot_HI*s_dust;
  else Rate_Phot_HI = Rate->Phot_HI;

  ne = Y->e*en_B;

  xTln = (Rate->Tln-cl->TlnMin)*cl->rDeltaTln;
  iTln = xTln;
  RT0 = (cl->RT+iTln*TABLEFACTOR);
  RT1 = RT0+TABLEFACTOR; 
  xTln = xTln-iTln;
#ifdef CUBICTABLEINTERP
  RT0d = RT0+1;
  RT1d = RT1+1;
  {
  double x2 = xTln*xTln;
  wTln1 = x2*(3-2*xTln);
  wTln0 = 1-wTln1;
  wTln0d = xTln*(1+xTln*(xTln-2));
  wTln1d = x2*(xTln-1);
/*  wTln1 = xTln;
  wTln0 = 1-xTln;
  wTln0d = 0;
  wTln1d = 0;*/
  }
#else  
  wTln1 = xTln;
  wTln0 = 1-xTln;
#endif

#define DTFRACLOWTCOOL 0.25
  if (Rate->T > cl->R.Tcmb*(1+DTFRACLOWTCOOL))
      LowTCool = TABLEINTERP( Cool_LowT )*cl->R.Cool_LowTFactor*en_B*ZMetal;
  else if (Rate->T < cl->R.Tcmb*(1-DTFRACLOWTCOOL))
      LowTCool = -TABLEINTERP( Cool_LowT )*cl->R.Cool_LowTFactor*en_B*ZMetal;
  else {
      double x = (Rate->T/cl->R.Tcmb-1)*(1./DTFRACLOWTCOOL);
      LowTCool = -TABLEINTERP( Cool_LowT )*cl->R.Cool_LowTFactor*en_B*ZMetal
	  *x*(3-3*fabs(x)+x*x);
      }

  *dEdotCool = 
#ifndef NOCOMPTON
      Y->e * cl->R.Cool_Comp * ( Rate->T - cl->R.Tcmb ) 
#endif
    + 
    ne * (TABLEINTERP( Cool_Brem_1 ) * ( Y->HII + Y->HeII ) +
	  TABLEINTERP( Cool_Brem_2 ) * Y->HeIII +
	  
	  cl->R.Cool_Diel_HeII * Y->HeII * Rate->Diel_HeII +

	  TABLEINTERP( Cool_Line_HI ) * Y->HI +
	  TABLEINTERP( Cool_Line_HeI ) * Y->HeI +
	  TABLEINTERP( Cool_Line_HeII ) * Y->HeII +

	  TABLEINTERP( Cool_Radr_HII ) * Y->HII * Rate->Radr_HII  +
	  TABLEINTERP( Cool_Radr_HeII ) * Y->HeII * Rate->Radr_HeII +
	  TABLEINTERP( Cool_Radr_HeIII ) * Y->HeIII * Rate->Radr_HeIII +

	  clCoolLineH2_e(Rate->T) * Y->H2  * CL_B_gm /*Re-added CL_B_gm 9/16/10 */+
	  cl->R.Cool_Coll_H2 * Y->H2 * Rate->Coll_e_H2 /* s_dust * s_self*/  /*shielding of collisions 11/18/10*/ +
	  cl->R.Cool_Coll_HI * Y->HI * Rate->Coll_HI +
	  cl->R.Cool_Coll_HeI * Y->HeI * Rate->Coll_HeI + 
	  cl->R.Cool_Coll_HeII * Y->HeII * Rate->Coll_HeII )
    +
    clCoolLineH2_HI(Rate->T) /* en_B * Y->H2 * Y->HI * CL_B_gm*/ /*Re-added CL_B_gm 9/16/10:
 clCoolLineH2_HI(Rate->T) [erg cm^3 s^-1 H2atom^-1 HIatom^-1] * en_B [atom/cm^3] * Y->H2 [H2atom atom^-1] * Y->HI [HIatom atom^-1] CL_B_gm [atom g^-1] => [ergs s^-1 g^-1]*/
    +
    clCoolLineH2_H2(Rate->T) /* en_B * Y->H2 * Y->H2 * CL_B_gm */
    +
    clCoolLineH2_He(Rate->T) /* en_B * Y->H2 * Y->HeI * CL_B_gm */ 
    +
    clCoolLineH2_HII(Rate->T) /* en_B * Y->H2 * Y->HII * CL_B_gm */
    +
    cl->R.Cool_Coll_H2 * Y->H2 * Rate->Coll_HI_H2 * Y->HI *en_B 
    +
    cl->R.Cool_Coll_H2 * Y->H2 * Rate->Coll_H2_H2 * Y->H2 *en_B 
    +
    cl->R.Cool_Coll_H2 * Y->H2 * Rate->Coll_HII_H2 * Y->HII *en_B 
    + 
      LowTCool
#ifndef NOMETALCOOLING
    +
    Rate->Cool_Metal
#endif
    ;

  *dEdotHeat =
#ifndef NOMETALCOOLING
      Rate->Heat_Metal
    +
#endif
    Y->H2 * cl->R.Heat_Phot_H2 * Rate->Phot_H2*s_dust*s_self + /*photon heating and dissociation CC*/
    Y->HI   * cl->R.Heat_Phot_HI * Rate_Phot_HI +
    Y->HeI  * cl->R.Heat_Phot_HeI * Rate->Phot_HeI +
    Y->HeII * cl->R.Heat_Phot_HeII * Rate->Phot_HeII;

  return *dEdotHeat - *dEdotCool;
}

/* Returns Heating - Cooling excluding External Heating, units of ergs s^-1 g^-1 
   Public interface CoolEdotInstantCode */
double clEdotInstant( COOL *cl, PERBARYON *Y, RATE *Rate, double rho, 
		      double ZMetal, double *dEdotHeat, double *dEdotCool )
{
  double en_B = rho*CL_B_gm;
  double ne,LowTCool;
  double s_dust, s_self, Rate_Phot_HI;

  s_dust = clDustShield(Y->HI*en_B, Y->H2*en_B, ZMetal, Rate->CorreLength);
  s_self = clSelfShield(Y->H2*en_B, Rate->CorreLength);
  if (cl->bShieldHI) Rate_Phot_HI = Rate->Phot_HI*s_dust;
  else Rate_Phot_HI = Rate->Phot_HI;

  ne = Y->e*en_B;
#define DTFRACLOWTCOOL 0.25
  if (Rate->T > cl->R.Tcmb*(1+DTFRACLOWTCOOL))
      LowTCool = clCoolLowT(Rate->T)*cl->R.Cool_LowTFactor*en_B*ZMetal;
  else if (Rate->T < cl->R.Tcmb*(1-DTFRACLOWTCOOL))
      LowTCool = -clCoolLowT(Rate->T)*cl->R.Cool_LowTFactor*en_B*ZMetal;
  else {
      double x = (Rate->T/cl->R.Tcmb-1)*(1./DTFRACLOWTCOOL);
      LowTCool = -clCoolLowT(Rate->T)*cl->R.Cool_LowTFactor*en_B*ZMetal
	  *x*(3-3*fabs(x)+x*x);
      }
  *dEdotCool = 
#ifndef NOCOMPTON
     Y->e * cl->R.Cool_Comp * ( Rate->T - cl->R.Tcmb ) 
#endif
    +
    ne * (clCoolBrem1(Rate->T) * ( Y->HII + Y->HeII ) +
	  clCoolBrem2(Rate->T) * Y->HeIII +
	  
	  cl->R.Cool_Diel_HeII * Y->HeII * Rate->Diel_HeII +

	  clCoolRadrHII(Rate->T) * Y->HII * Rate->Radr_HII  +
	  clCoolRadrHeII(Rate->T) * Y->HeII * Rate->Radr_HeII +
	  clCoolRadrHeIII(Rate->T) * Y->HeIII * Rate->Radr_HeIII +

	  clCoolLineHI(Rate->T) * Y->HI +
	  clCoolLineHeI(Rate->T) * Y->HeI +
	  clCoolLineHeII(Rate->T) * Y->HeII +
	  
#ifndef NOMOLECULARHCOOLING /* Compile flag to take out molecular cooling for debugging purposes*/
	  clCoolLineH2_e(Rate->T) * Y->H2 * CL_B_gm + /*Re-added CL_B_gm 9/16/10 */
	  cl->R.Cool_Coll_H2 * Y->H2 * Rate->Coll_e_H2 +
#endif
	  cl->R.Cool_Coll_HI * Y->HI * Rate->Coll_HI +
	  cl->R.Cool_Coll_HeI * Y->HeI * Rate->Coll_HeI + 
	  cl->R.Cool_Coll_HeII * Y->HeII * Rate->Coll_HeII )
#ifndef NOMOLECULARHCOOLING /* Compile flag to take out molecular cooling for debugging purposes*/
    +
    clCoolLineH2_HI(Rate->T) * en_B * Y->H2 * Y->HI * CL_B_gm /*Re-added CL_B_gm 9/16/10:
    clCoolLineH2_HI(Rate->T) [erg cm^3 s^-1 H2atom^-1 HIatom^-1] * en_B [atom/cm^3] * Y->H2 [H2atom atom^-1] * Y->HI [HIatom atom^-1] CL_B_gm [atom g^-1] => [ergs s^-1 g^-1]*/
    +
    clCoolLineH2_H2(Rate->T) * en_B * Y->H2 * Y->H2 * CL_B_gm
    +
    clCoolLineH2_He(Rate->T) * en_B * Y->H2 * Y->HeI * CL_B_gm 
    +
    clCoolLineH2_HII(Rate->T) * en_B * Y->H2 * Y->HII  * CL_B_gm
    +
    cl->R.Cool_Coll_H2 * Y->H2 * Rate->Coll_HI_H2  * Y->HI * en_B 
    +
    cl->R.Cool_Coll_H2 * Y->H2 * Rate->Coll_H2_H2 * Y->H2 * en_B 
    +
    cl->R.Cool_Coll_H2 * Y->H2 * Rate->Coll_HII_H2 * Y->HII * en_B   
#endif
    +
      LowTCool
#ifndef NOMETALCOOLING
    + 
    Rate->Cool_Metal
#endif
    ;

  *dEdotHeat =
#ifndef NOMETALCOOLING
      Rate->Heat_Metal
    +
#endif
#ifndef NOMOLECULARHCOOLING
    Y->H2   * cl->R.Heat_Phot_H2 * Rate->Phot_H2*s_dust*s_self +
#endif
    Y->HI   * cl->R.Heat_Phot_HI * Rate_Phot_HI +/* This produces very large amounts of HI shielding but is consistent with everything else CC*/
    Y->HeI  * cl->R.Heat_Phot_HeI  * Rate->Phot_HeI +
    Y->HeII * cl->R.Heat_Phot_HeII * Rate->Phot_HeII;

#ifdef COOLDEBUG
  /*Debugging information on cooling and heating */
  if (cl->p->iOrder == PARTICLEIORD) { 
    if (0) { 
      printf("\nEdot-> Total: %#2e, T: %#2f; rho: %#2f; Shield: %#2e\n",*dEdotHeat - *dEdotCool,Rate->T,en_B,s_dust*s_self);
      printf("Edot-> HI: %#2e; H2: %#2e; HeI: %#2e, HeII: %#2e\n",Y->HI,Y->H2,Y->HeI,Y->HeII);
      printf("Edot-> Line: %#2e, Coll %#2e, %#2e, Phot: %#2e, %#2e\n",
	     ne*clCoolLineH2_e(Rate->T) * Y->H2 * CL_B_gm /*Re-added CL_B_gm 9/16/10, CC */ + 
	     clCoolLineH2_HI(Rate->T) * en_B * Y->H2 * Y->HI * CL_B_gm  /*Re-added CL_B_gm 9/16/10, CC:
 clCoolLineH2_H(Rate->T) [erg cm^3 s^-1 H2atom^-1 HIatom^-1] * en_B [atom/cm^3] * Y->H2 [H2atom atom^-1] * Y->HI [HIatom atom^-1] CL_B_gm [atom g^-1] => [ergs s^-1 g^-1]*/ + 
	     clCoolLineH2_H2(Rate->T) * en_B * Y->H2 * Y->H2 * CL_B_gm + 
	     clCoolLineH2_He(Rate->T) * en_B * Y->H2 * Y->HeI * CL_B_gm + 
	     clCoolLineH2_HII(Rate->T) * en_B * Y->H2 * Y->HII * CL_B_gm ,
	     cl->R.Cool_Coll_H2*Y->H2*Rate->Coll_e_H2*ne + cl->R.Cool_Coll_H2*Y->H2*Rate->Coll_HI_H2*Y->HI*en_B + cl->R.Cool_Coll_H2*Y->H2*Rate->Coll_H2_H2*Y->H2*en_B + cl->R.Cool_Coll_H2*Y->H2*Rate->Coll_HII_H2*Y->HII*en_B,
	     (cl->R.Cool_Coll_H2*Y->H2*Rate->Coll_e_H2*ne + cl->R.Cool_Coll_H2*Y->H2*Rate->Coll_HI_H2*Y->HI*en_B + cl->R.Cool_Coll_H2*Y->H2*Rate->Coll_H2_H2*Y->H2*en_B + cl->R.Cool_Coll_H2*Y->H2*Rate->Coll_HII_H2*Y->HII*en_B), /* *s_dust*s_self,removing shielding of collisions 9/22/10, CC*/
	     Y->H2*cl->R.Heat_Phot_H2*Rate->Phot_H2,
	     Y->H2*cl->R.Heat_Phot_H2*Rate->Phot_H2*s_self*s_dust);
      printf("Edot-> Comp: %#2e, Brems: %#2e, Diel: %#2e, Radr: %#2e, Line: %#2e, Coll: %#2e\n",
	     Y->e * cl->R.Cool_Comp * ( Rate->T - cl->R.Tcmb),
	     ne*(clCoolBrem1(Rate->T) * ( Y->HII + Y->HeII ) + clCoolBrem2(Rate->T) * Y->HeIII),
	     ne*(cl->R.Cool_Diel_HeII * Y->HeII * Rate->Diel_HeII),
	     ne*(clCoolRadrHII(Rate->T) * Y->HII * Rate->Radr_HII + clCoolRadrHeII(Rate->T) * Y->HeII * Rate->Radr_HeII + clCoolRadrHeIII(Rate->T) * Y->HeIII * Rate->Radr_HeIII),
	     ne*(clCoolLineHI(Rate->T) * Y->HI + clCoolLineHeI(Rate->T) * Y->HeI + clCoolLineHeII(Rate->T) * Y->HeII),
	     ne*(cl->R.Cool_Coll_HI * Y->HI * Rate->Coll_HI +cl->R.Cool_Coll_HeI * Y->HeI * Rate->Coll_HeI + cl->R.Cool_Coll_HeII * Y->HeII * Rate->Coll_HeII)
	     );
      printf("Edot-> LowTCool: %#2e, MetalCool: %#2e, Sum: %#2e\n",
	     LowTCool,
	     Rate->Cool_Metal,
	     Y->e * cl->R.Cool_Comp*(Rate->T-cl->R.Tcmb) + 
	     ne*(clCoolBrem1(Rate->T)*(Y->HII + Y->HeII)+clCoolBrem2(Rate->T)*Y->HeIII) +  
	     ne*(cl->R.Cool_Diel_HeII * Y->HeII * Rate->Diel_HeII) + 
	     ne*(clCoolRadrHII(Rate->T) * Y->HII * Rate->Radr_HII + clCoolRadrHeII(Rate->T) * Y->HeII * Rate->Radr_HeII + clCoolRadrHeIII(Rate->T) * Y->HeIII * Rate->Radr_HeIII) + 
	     ne*(clCoolLineHI(Rate->T)*Y->HI+clCoolLineHeI(Rate->T)*Y->HeI+clCoolLineHeII(Rate->T)*Y->HeII) +
	     ne*(cl->R.Cool_Coll_HI*Y->HI*Rate->Coll_HI+cl->R.Cool_Coll_HeI*Y->HeI*Rate->Coll_HeI+cl->R.Cool_Coll_HeII*Y->HeII*Rate->Coll_HeII) + 
	     LowTCool + Rate->Cool_Metal
	     );
      printf("Edot-> Phot: %#2e, Metal_Heat: %#2e",
	     Y->HI*cl->R.Heat_Phot_HI*Rate->Phot_HI  + Y->HeI*cl->R.Heat_Phot_HeI*Rate->Phot_HeI + Y->HeII*cl->R.Heat_Phot_HeII*Rate->Phot_HeII + Y->H2*cl->R.Heat_Phot_H2*Rate->Phot_H2*s_dust*s_self,Rate->Heat_Metal);
      printf("\n");
    }
  }
#endif


  return *dEdotHeat - *dEdotCool;
}

/*
 *  We solve the equation:  
 *  Eout = Ein + dt * Cooling ( Yin, Yout, Ein, Eout )
 *
 *  E erg/gm, ExternalHeating erg/gm/s, rho gm/cc, dt sec
 * 
 */

double clfTemp(void *Data, double T) 
{
  clDerivsData *d = Data; 
  d->its++;
  CLRATES( d->cl, &d->Rate, T, d->rho, d->ZMetal, d->Rate.CorreLength/(d->cl->dKpcUnit * 3.08568025e21 * d->cl->dExpand), d->Rate.LymanWernerCode );
  clRateMetalTable(d->cl, &d->Rate, T, d->rho, d->Y_H, d->ZMetal); 
  clAbunds( d->cl, &d->Y, &d->Rate, d->rho, d->ZMetal);

  return d->E-clThermalEnergy( d->Y.Total, T );
}

void clTempIteration( clDerivsData *d )
{
 double T,TA,TB;
 double Y_Total0 = (d->Y_H+d->Y_He+d->Y.H2)*.9999; 
 double Y_Total1 = (d->Y_eMax+d->Y_H+d->Y_He)*1.0001;  
 d->its = 0;
 if (d->E <= 0) T=d->cl->TMin;
 else {

   TB = clTemperature( Y_Total0, d->E );
   TA = clTemperature( Y_Total1, d->E );
   if (TA > TB) { T=TA; TA=TB; TB=T; }

   T = RootFind(clfTemp, d, TA, TB, EPSTEMP*TA );
 } 
 d->its++;
 CLRATES( d->cl, &d->Rate, T, d->rho, d->ZMetal, d->Rate.CorreLength/(d->cl->dKpcUnit * 3.08568025e21 * d->cl->dExpand), d->Rate.LymanWernerCode );
 clRateMetalTable(d->cl, &d->Rate, T, d->rho, d->Y_H, d->ZMetal); 
 clAbunds( d->cl, &d->Y, &d->Rate, d->rho, d->ZMetal);
}

void clDerivs(double x, const double *y, double *dGain, double *dLoss,
	     void *Data) {
  clDerivsData *d = Data;
  double T,ne,nHI, nHII, internalheat = 0, externalheat = 0;
  double internalcool = 0;
  double en_B = d->rho*CL_B_gm;
  double s_dust, s_self, nH2, nHminus, Rate_Phot_HI;
  
  d->E = y[0];
  d->Y.HI = y[1];
  d->Y.H2 = y[4];
  d->Y.HII = d->Y_H - d->Y.HI - d->Y.H2*2.0; /*Don't forget about molec H CC*/
  if(d->Y.HII < 0) d->Y.HII = 0;  
 
  d->Y.HeI = y[2];
  d->Y.HeII = y[3];
  d->Y.HeIII = d->Y_He - d->Y.HeI - d->Y.HeII; 
  if(d->Y.HeIII < 0) d->Y.HeIII = 0; 
 
  d->Y.e = d->Y.HII + d->Y.HeII + 2*d->Y.HeIII;
#ifdef Y_EMIN
  if (d->Y.e < Y_EMIN) d->Y.e = Y_EMIN;
#endif

  d->Y.Total = d->Y.e + d->Y_H + d->Y_He + d->ZMetal/MU_METAL - d->Y.H2;  /* H total from cl now -- in future from particle */ 
  T = clTemperature( d->Y.Total, d->E );
  CLRATES( d->cl, &d->Rate, T, d->rho, d->ZMetal, d->columnL, d->dLymanWerner);
  s_dust = clDustShield(d->Y.HI*en_B, d->Y.H2*en_B, d->ZMetal, d->Rate.CorreLength); /* dust shielding*/
  s_self = clSelfShield(d->Y.H2*en_B, d->Rate.CorreLength); /* H2 self shielding*/
  if (d->cl->bShieldHI) Rate_Phot_HI = d->Rate.Phot_HI*s_dust;/* If bShieldHI was set in parameter file, modify the photoionizing flux for HI by the dust shielding*/
  else Rate_Phot_HI = d->Rate.Phot_HI;

  externalheat = d->ExternalHeating;
  if (d->bCool) {
    clRateMetalTable(d->cl, &d->Rate, T, d->rho, d->Y_H, d->ZMetal);
    CLEDOTINSTANT( d->cl, &d->Y, &d->Rate, d->rho, d->ZMetal, &internalheat, &internalcool );
    
    dGain[0] = internalheat;
    dLoss[0] = internalcool;
  }
  if(externalheat > 0.0)
      dGain[0] += externalheat;
  else
      dLoss[0] -= externalheat;

  ne =  en_B*d->Y.e;
  nHI = en_B*d->Y.HI;
  nH2 = en_B*d->Y.H2;
  nHII = en_B*d->Y.HII;
  nHminus = d->Rate.HI_e* d->Y.HI*d->Y.e/(d->Rate.HI_Hm*d->Y.HI + d->Rate.Coll_Hm_HII*d->Y.HI + d->Rate.Coll_Hm_e*d->Y.e);

#ifdef COOLDEBUG
  if (0) {
  if (d->cl->p->iOrder == PARTICLEIORD) printf("%d, x: %e, en_B: %e, T: %e, HI: %e, H2: ,%e\n",d->cl->p->iOrder,x,en_B,T,d->Y.HI,d->Y.H2);
    if (d->cl->p->iOrder == PARTICLEIORD) printf("GasForm: %e,  PhotH2: %e = %e * %e,  Coll_e: %e, Coll_H2: %e, Coll_HI: %e, Coll_HII: %e, Shield: %e, Shielded: %e, Dust: %e \n",	      
	   d->Y.HI*nHminus*en_B*d->Rate.HI_Hm,
          -1.0*d->Y.H2*d->Rate.Phot_H2,d->Y.H2,d->Rate.Phot_H2,
          -1.0*d->Y.H2*d->Rate.Coll_e_H2*ne,
          -1.0*d->Y.H2*d->Rate.Coll_H2_H2*nH2,
          -1.0*d->Y.H2*d->Rate.Coll_HI_H2*nHI,
          -1.0*d->Y.H2*d->Rate.Coll_HII_H2*nHI,
               s_dust*s_self,
       (d->Y.HI*nHminus*en_B*d->Rate.HI_Hm - /*gas phase formation, adel 97 */
	     d->Y.H2*d->Rate.Phot_H2*s_dust*s_self -
	     d->Y.H2*d->Rate.Coll_e_H2*ne - 
	     d->Y.H2*d->Rate.Coll_H2_H2*nH2 -
	     d->Y.H2*d->Rate.Coll_HI_H2*nHI -
             d->Y.H2*d->Rate.Coll_HII_H2*nHII),
       (d->Y.HI + 2.0*d->Y.H2)*d->Rate.DustForm_H2*nHI);
  }
#endif 
  dGain[4] = d->Y.HI*nHminus*en_B*d->Rate.HI_Hm + /*gas phase formation of H2, adel 97 */
            (d->Y.HI + 2.0*d->Y.H2)*d->Rate.DustForm_H2*nHI;
  dLoss[4] = d->Y.H2*d->Rate.Phot_H2*s_dust*s_self +
	     d->Y.H2*d->Rate.Coll_e_H2*ne +
	     d->Y.H2*d->Rate.Coll_H2_H2*nH2 +
	     d->Y.H2*d->Rate.Coll_HI_H2*nHI +
             d->Y.H2*d->Rate.Coll_HII_H2*nHII;
  dGain[1] = ne*d->Y.HII*d->Rate.Radr_HII +
             2.0*dLoss[4]; /*Photo or collisionally dissociated H2*/
  dLoss[1] = ne*d->Y.HI*d->Rate.Coll_HI + 
             d->Y.HI*Rate_Phot_HI + /*Adding in molec H and shielding, should possibly add shielding for others*/
             2.0*dGain[4]; /*Formation of H2*/
  dGain[2] = ne*d->Y.HeII*d->Rate.Totr_HeII;
  dLoss[2] = ne*d->Y.HeI*d->Rate.Coll_HeI + 
             d->Y.HeI*d->Rate.Phot_HeI;
  dGain[3] = ne*d->Y.HeIII*d->Rate.Radr_HeIII +
             dLoss[2];
  dLoss[3] = ne*d->Y.HeII*d->Rate.Coll_HeII + 
             d->Y.HeII*d->Rate.Phot_HeII +
             dGain[2];
}

void clSetyscale( COOL *cl, double Y_H, double Y_He, double *y, double *yscale) {
    double YHII, YHeIII;

    yscale[0] = y[0]; /*e */
    yscale[1] = y[1]; /*Y->HI*/;
    yscale[2] = y[2]; /* Y->HeI*/;
    yscale[3] = y[3]; /* Y->HeII*/;
    yscale[4] = y[4]; /*Y->H2*/;

    /* MW 1e10 Msun in 10kpc radius => nH = 1e22 cm^-2
       If nHI/nH = 1e-12 then nHI = 1e10 cm^-2 -- undetectable */
#define YSCALEMIN 1e-10
#define YeSCALEMIN 1e-4

    /* Make sure error in YHII and YHeIII is controlled too
       |delta YHI| = |delta YHII|
       |delta YHI| = |delta YH2| = |delta YHII| 
       |delta YHeIII| < |delta YHeI| + |delta YHeII| so control those 
       to control delta YHeIII
       NB: Also gives error control on Y_e which was absent before! */
    
#define CONTROLYHII

    /* HI */
#ifdef CONTROLYHII
    if ((YHII = Y_H - y[1] - 2.0*y[4]) < YeSCALEMIN) YHII = YeSCALEMIN; /*Adding in molec H*/
    if (yscale[1] > YHII) yscale[1] = YHII; 
    else 
#endif
      if (yscale[1] < YSCALEMIN) yscale[1] = YSCALEMIN; 

    /* H2 */
#ifdef CONTROLYHII
    if (yscale[4] > YHII) yscale[4] = YHII;
    else
#endif
      if (yscale[4] < YSCALEMIN) yscale[4] = YSCALEMIN;

#ifdef CONTROLYHII
    if ((YHeIII = 0.5*(Y_He - y[2] - y[3])) <  (0.5*YeSCALEMIN)) YHeIII = (0.5*YeSCALEMIN);
    if (yscale[2] > YHeIII) yscale[2] =  YHeIII; 
    else 
#endif
      if (yscale[2] < YSCALEMIN) yscale[2] = YSCALEMIN;


#ifdef CONTROLYHII
    if (yscale[3] > YHeIII) yscale[3] =  YHeIII; 
    else 
#endif
      if (yscale[3] < YSCALEMIN) yscale[3] = YSCALEMIN; 

} 

void clIntegrateEnergy(COOL *cl, clDerivsData *clData, PERBARYON *Y, double *E, 
		       double ExternalHeating, double rho, double ZMetal, double tStep, double columnL, double dLymanWerner  ) {

/* Make sure no values are near roundoff ~ (1e-14)*Y_i */
#define YHMIN 1e-12
#define YHeMIN 1e-13
#define YH2MIN 1e-12/2.0 
  
  const int array_length=5; /*Arrays expanded for H2*/

  double y[array_length-1],yin[array_length-1],EMin,YTotal ;
  double dHeat, dCool;
  double t=0;
  double  s_dust, s_self;
  clDerivsData *d = clData;
  STIFF *sbs = d->IntegratorContext;
  int its = 0;
  FILE *fp; 
  double T;
 
  if (tStep == 0) return;
  d->bCool = 1;
  if (tStep < 0) {
      tStep = fabs(tStep);
      d->bCool = 0;
      }

#ifdef Y_EMIN
  if (Y->e < Y_EMIN) Y->e = Y_EMIN;
#endif

  y[0] = yin[0] = *E;
  y[1] = yin[1] = Y->HI;
  y[2] = yin[2] = Y->HeI;
  y[3] = yin[3] = Y->HeII;
  y[4] = yin[4] = Y->H2; /*Inititalyzing the number of H2 molecules*/ 

  d->rho = rho;
  d->ExternalHeating = ExternalHeating;
  d->ZMetal = ZMetal; 
  d->dLymanWerner = dLymanWerner;
  d->columnL = columnL;
  clSetAbundanceTotals( cl, ZMetal, &d->Y_H, &d->Y_He, &d->Y_eMax );
  
/* H, He total from cl now -- in future from particle */
  YTotal = Y->HII + Y->HeII + 2*Y->HeIII + d->Y_H + d->Y_He + d->ZMetal/MU_METAL -  Y->H2; 

  EMin = clThermalEnergy( YTotal, cl->TMin );

#ifdef COOLDEBUG
  if (cl->p->iOrder == PARTICLEIORD) { 
    printf("\ny[0]: %e, y[1]: %e, y[2]: %e, y[3]: %e, y[4] %e\n",y[0],y[1],y[2],y[3],y[4]);
    /*    clDerivs( d, t, yin, &dHeat, &dCool ); *//*, sqrt(cl->p->fBall2/2.0) );*/
    printf("tStep %g \n", tStep);
    printf("%d: rho %g \n", cl->p->iOrder,rho*CL_B_gm); 
    T = clTemperature( Y->Total, *E );
    printf("Temperature %g \n", T);
    printf("Y e:%g Total:%g H2:%g HI:%g HII:%g HeI:%g HeII:%g HeIII:%g\n",Y->e, YTotal, Y->H2, Y->HI, Y->HII, Y->HeI, Y->HeII, Y->HeIII);
    /*   printf("dHeat[0]: %e, dHeat[1]: %e, dHeat[2]: %e, dHeat[3]: %e, dHeat[4] %e\n",dHeat[0],dHeat[1],dHeat[2],dHeat[3],dHeat[4]);
	 printf("dCool[0]: %e, dCool[1]: %e, dCool[2]: %e, dCool[3]: %e, dCool[4] %e\n",dCool[0],dCool[1],dCool[2],dCool[3],dCool[4]);*/
/*  printf("Cooling p %i: %f %f %g %g %g %g %f %f %g %g %g %g \n",cl->p->iOrder,cl->z,d->Rate.T,rho,cl->p->fMass,cl->p->fBall2,*E,ExternalHeating, ZMetal, d->Rate.Cool_Metal, d->Rate.Heat_Metal, clCoolTotal(cl, Y, &d->Rate, rho, ZMetal), clHeatTotal(cl, Y, &d->Rate, rho) ); 
    printf("Cooling p %i \n", cl->p->iOrder); 
    printf("temperature %g \n", d->Rate.T);
    printf("redshift %g rho %g \n",cl->z, rho*CL_B_gm);
    printf("PdV %e \n",ExternalHeating); 
    printf("Y_H %g, Y_He %g \n", d->Y_H, d->Y_He); 
    clPrintCool(cl, Y, &d->Rate, rho);  
    printf("delt %e \n", tStep);  
    printf("Ein %e \n", *E);
    printf("dydt %e \n", (dydt-1)[1]); 
    printf("totalcool %g \n",clCoolTotal(cl, Y, &d->Rate, rho, ZMetal)); 
    printf("totalheat %g \n",clHeatTotal(cl, Y, &d->Rate, rho));
    printf("metalcool %g \n", d->Rate.Cool_Metal); 
    printf("metalheat %g \n", d->Rate.Heat_Metal); */
    printf("Internal %e \n",CLEDOTINSTANT( d->cl, &d->Y, &d->Rate, d->rho, d->ZMetal, &dHeat, &dCool ));
    printf("PdV %e \n",ExternalHeating); 
    /*    printf("Edot %g, %g \n", dHeat[0], dCool[0]);*/
    printf("dbCool %d \n",d->bCool);
    printf("E %g \n", *E); 

    
    FILE *fpdebug;
    fpdebug = fopen("cooldebug.txt","a");
    fprintf(fpdebug, "\ntStep %g \n", tStep); 
    fprintf(fpdebug, "Y e:%g Total:%g H2:%g HI:%g HII:%g HeI:%g HeII:%g HeIII:%g\n",Y->e, Y->Total, Y->H2, Y->HI, Y->HII, Y->HeI, Y->HeII, Y->HeIII); 
    clPrintCoolFile(cl, Y, &d->Rate, rho,d->ZMetal, fpdebug);       
    fprintf(fpdebug, "Metalicity %g \n", d->ZMetal);
    fprintf(fpdebug, "redshift %g \n", cl->z);
    fprintf(fpdebug, "rho %g \n", rho*CL_B_gm); 
    fprintf(fpdebug, "temperature %g \n", T);
    fprintf(fpdebug, "Y %g %g %g %g %g %g %g %g\n",Y->e, Y->Total, Y->H2, Y->HI, Y->HII, Y->HeI, Y->HeII, Y->HeIII); 
    fprintf(fpdebug, "Internal %e \n",CLEDOTINSTANT( d->cl, &d->Y, &d->Rate, d->rho, d->ZMetal, &dHeat, &dCool ));
    fprintf(fpdebug, "PdV %e \n",ExternalHeating);
    /*    fprintf(fpdebug, "Edot %e, %e \n", dHeat[0], dCool[0] );  */
    fprintf(fpdebug, "E %g \n", *E); 
    /*fprintf(fpdebug, "units %e %e %e %e %e \n", cl->dGmPerCcUnit, cl->dComovingGmPerCcUnit, cl->dErgPerGmUnit, cl->dSecUnit, cl->dKpcUnit );*/
    fclose(fpdebug);
  }
 #endif

 {
     double ymin[array_length];
     ymin[1] = YHMIN;
     ymin[2] = YHeMIN;
     ymin[3] = YHeMIN;
     ymin[4] = YH2MIN;

      YTotal = (d->Y_H - ymin[4]) + (d->Y_H - ymin[1] - 2.0*ymin[4]) + d->Y_He + ymin[3] + 
	2.0*(d->Y_He - ymin[2] - ymin[3]) + d->ZMetal/MU_METAL;
      EMin = clThermalEnergy( YTotal, cl->TMin );
      ymin[0] = EMin;
      StiffSetYMin(sbs, ymin);
      
      StiffStep( sbs, y, t, tStep);
      if(y[1] > d->Y_H) y[1] = d->Y_H;
      if(y[4] > d->Y_H/2.0) y[4] = d->Y_H/2.0;
      if(y[2] > d->Y_He) y[2] = d->Y_He;
      if(y[3] > d->Y_He) y[3] = d->Y_He;
      if (fabs(y[2])+fabs(y[3]) > d->Y_He) {
	  if(y[2] > y[3]) y[3] = d->Y_He - y[2]; 
	  else y[2] = d->Y_He - y[3]; 
	  } 

      if (fabs(y[1])+fabs(y[4]*2.0) > d->Y_H) { 
        if(y[1] > y[4]*2) y[4] = (d->Y_H - y[1])/2.0; 
          else y[1] = d->Y_H - y[4]*2.0; 
      }

      if(y[4] < YH2MIN) y[4] = YH2MIN;
      if(y[1] < YHMIN) y[1] = YHMIN;
      if (d->Y_H - y[1] - 2.0*y[4] < YHMIN) {
        if (y[1] >= 2.0*y[4]) y[1] = d->Y_H - 2.0*y[4] - YHMIN;
        else                  y[4] =(d->Y_H - y[1] - YHMIN)/2.0;

      }
/*    if(y[4] < YH2MIN) {
	y[4] = YH2MIN;
	if (d->Y_H - y[1] - 2.0*y[4] < YHMIN) y[4] = (d->Y_H - y[1] - YHMIN)/2.0;
      }
      if(y[1] < YHMIN) {
	y[1] = YHMIN;
	if (d->Y_H - y[1] - 2.0*y[4] < YHMIN) y[1] = d->Y_H - 2.0*y[4] - YHMIN;
	}*/

      if(y[2] < YHeMIN) y[2] = YHeMIN;
      if(y[3] < YHeMIN) y[3] = YHeMIN;
      if (d->Y_He - y[2] - y[3] < YHeMIN) {
        if (y[2] >= y[3]) y[2] = d->Y_He - y[3] - YHeMIN;
        else              y[3] = d->Y_He - y[2] - YHeMIN; 
      }

/*    if(y[2] < YHeMIN) {
	y[2] = YHeMIN;
	if (d->Y_He - y[2] - y[3] < YHeMIN) y[3] = d->Y_He - y[2] - YHeMIN;
      }
      if(y[3] < YHeMIN) {
	y[3] = YHeMIN;
	if (d->Y_He - y[2] - y[3] < YHeMIN) y[2] = d->Y_He - y[3] - YHeMIN;
	}*/

      YTotal = (d->Y_H) - y[4] + (d->Y_H - y[1] - 2*y[4]) + d->Y_He + y[3] +
	2.0*(d->Y_He - y[2] - y[3]) + d->ZMetal/MU_METAL;
      EMin = clThermalEnergy( YTotal, cl->TMin );

#ifdef ASSERTENEG      
      assert(*y > 0.0);
#else
      if (y[0] < EMin) {
	y[0] = EMin;
/*	y[1] = YHMIN;
	y[2] = d->Y_He;
	y[3] = YHeMIN;
	y[4] = d->Y_H/2.0;
	printf("Emin: %e\n", EMin);
	break;*/
      }
#endif   
   cl->its = 1;
   }

#ifdef COOLDEBUG
      if(cl->p->iOrder == PARTICLEIORD){
	printf("%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
	       t,
	       d->Rate.T, 
	       d->rho*CL_B_gm, 
	       d->ZMetal, 
	       d->columnL, 
	       d->dLymanWerner,
	       Y->H2,
	       y[4],
	       0, /*dydt[4],*/
          -1.0*d->Y.H2*d->Rate.Phot_H2,
          -1.0*d->Y.H2*d->Rate.Coll_e_H2*rho*CL_B_gm*Y->e,
          -1.0*d->Y.H2*d->Rate.Coll_H2_H2*rho*CL_B_gm*Y->H2,
          -1.0*d->Y.H2*d->Rate.Coll_HI_H2*rho*CL_B_gm*Y->HI,
          -1.0*d->Y.H2*d->Rate.Coll_HII_H2*rho*CL_B_gm*Y->HII,
               s_dust*s_self,
               (d->Y.HI + 2.0*d->Y.H2)*d->Rate.DustForm_H2*rho*CL_B_gm*Y->HI);
      }
#endif

   *E = y[0];
   d->E = y[0];
   Y->HI = y[1];
   Y->HeI = y[2];
   Y->HeII = y[3];
   Y->HeIII = d->Y_He - Y->HeI - Y->HeII;
   Y->H2 = y[4];
   Y->HII = d->Y_H - Y->HI - 2.0*Y->H2;
   Y->e = Y->HII + Y->HeII + 2*Y->HeIII;
   Y->Total =  Y->e + d->Y_H + d->Y_He + d->ZMetal/MU_METAL - Y->H2; /*as two hydrogen atoms make one molecular hydrogen particle, subtract off number of molecules to avoid double counting particles CC */

#ifdef COOLDEBUG
   /* Table of cooling and heating terms that can be read into an idl function to create a nice cooling curve*/
   if (cl->p->iOrder == PARTICLEIORD) {
     FILE *cooldebug;
     cooldebug = fopen("cooldebug_table.txt","a");
     double en_B = d->rho*CL_B_gm, Edot, LowTCool;
     double ne = en_B*Y->e;
     RATE *Rate = &d->Rate;

     T = clTemperature( Y->Total, *E );
     CLRATES( d->cl, &d->Rate, T, d->rho, d->ZMetal, d->columnL, d->dLymanWerner);
     s_dust = clDustShield(Y->HI*en_B, Y->H2*en_B, d->ZMetal, Rate->CorreLength);
     s_self = clSelfShield(Y->H2*en_B, Rate->CorreLength);
#define DTFRACLOWTCOOL 0.25
     if (Rate->T > cl->R.Tcmb*(1+DTFRACLOWTCOOL))
       LowTCool = clCoolLowT(Rate->T)*cl->R.Cool_LowTFactor*en_B*ZMetal;
     else if (Rate->T < cl->R.Tcmb*(1-DTFRACLOWTCOOL))
       LowTCool = -clCoolLowT(Rate->T)*cl->R.Cool_LowTFactor*en_B*ZMetal;
     else {
       double x = (Rate->T/cl->R.Tcmb-1)*(1./DTFRACLOWTCOOL);
       LowTCool = -clCoolLowT(Rate->T)*cl->R.Cool_LowTFactor*en_B*ZMetal
	 *x*(3-3*fabs(x)+x*x);
      }

     Edot = clEdotInstant( cl, &d->Y, &d->Rate, d->rho, d->ZMetal, &dHeat, &dCool );
     /*T [K]  Density  e  H2  HI  HII  HeI  HeII  HeIII  shield  Edot  Cool_Comp  Brem  Diel  Radr  LineHI  LineHeI  LineHeII  LineH2_H  LineH2_H2  LineH2_He  LineH2_e  LineH2_HII  Coll  CollH2  LowT  MCool  Cool  Phot  PhotH2  MHeat     */
     fprintf(cooldebug,"%15f, %15e, %15f, %15e, %15g, %15g, %15g, %15g, %15g, %15g, %15e, %15e, %15e, %15e, %15e, %15e, %15e, %15e, %15e, %15e, %15e, %15e, %15e, %15e, %15e, %15e, %15e, %15e, %15e, %15e, %15e, %15e, %15e, %15e\n",
             T, y[0], en_B, ne, Y->H2, Y->HI, Y->HII, Y->HeI, Y->HeII, Y->HeIII, s_dust*s_self, Edot,
	     -1.0*Y->e*cl->R.Cool_Comp*(T - cl->R.Tcmb),
	     -1.0*ne*(clCoolBrem1(T) * ( Y->HII + Y->HeII ) + clCoolBrem2(T) * Y->HeIII),
	     -1.0*ne*(cl->R.Cool_Diel_HeII * Y->HeII * Rate->Diel_HeII),
	     -1.0*ne*(clCoolRadrHII(T)*Y->HII*Rate->Radr_HII + 
                      clCoolRadrHeII(T)*Y->HeII*Rate->Radr_HeII + 
                      clCoolRadrHeIII(T)*Y->HeIII*Rate->Radr_HeIII),
	     -1.0*ne*(clCoolLineHI(T)*Y->HI),
	     -1.0*ne*(clCoolLineHeI(T)*Y->HeI),
	     -1.0*ne*(clCoolLineHeII(T)*Y->HeII),
	     -1.0*clCoolLineH2_HI(Rate->T)*en_B*  Y->HI*Y->H2* CL_B_gm, 
	     -1.0*clCoolLineH2_H2(Rate->T)*en_B* Y->H2*Y->H2* CL_B_gm, 
	     -1.0*clCoolLineH2_He(Rate->T)*en_B* Y->HeI*Y->H2*CL_B_gm,
	     -1.0*clCoolLineH2_e(Rate->T)*ne*          Y->H2* CL_B_gm,
	     -1.0*clCoolLineH2_HII(Rate->T)*en_B*Y->HII*Y->H2*CL_B_gm,
	     -1.0*ne*(cl->R.Cool_Coll_HI*Y->HI*Rate->Coll_HI + cl->R.Cool_Coll_HeI*Y->HeI*Rate->Coll_HeI + cl->R.Cool_Coll_HeII*Y->HeII*Rate->Coll_HeII),
	     -1.0*cl->R.Cool_Coll_H2*Y->H2*Rate->Coll_e_H2 *        ne, /**s_dust*s_self, removing shielding of collisions 9/22/10, CC*/
	     -1.0*cl->R.Cool_Coll_H2*Y->H2*Rate->Coll_HI_H2 *Y->HI*en_B, /**s_dust*s_self, removing shielding of collisions 9/22/10, CC*/
	     -1.0*cl->R.Cool_Coll_H2*Y->H2*Rate->Coll_H2_H2*Y->H2*en_B, /**s_dust*s_self, removing shielding of collisions 9/22/10 CC*/
	     LowTCool,
	     -1.0*Rate->Cool_Metal,
	     -1.0*(   Y->e*cl->R.Cool_Comp*(T-cl->R.Tcmb) + 
	              ne*(clCoolBrem1(T)*(Y->HII + Y->HeII) + clCoolBrem2(T)*Y->HeIII) +  
	              ne*(cl->R.Cool_Diel_HeII*Y->HeII*Rate->Diel_HeII) + 
	              ne*(clCoolRadrHII(T)*Y->HII*Rate->Radr_HII + clCoolRadrHeII(T)*Y->HeII*Rate->Radr_HeII + clCoolRadrHeIII(T)*Y->HeIII*Rate->Radr_HeIII) + 
	              ne*(clCoolLineHI(T)*Y->HI + clCoolLineHeI(T)*Y->HeI + clCoolLineHeII(T)*Y->HeII) +
	              ne*(cl->R.Cool_Coll_HI*Y->HI*Rate->Coll_HI + cl->R.Cool_Coll_HeI*Y->HeI*Rate->Coll_HeI + cl->R.Cool_Coll_HeII*Y->HeII*Rate->Coll_HeII) + 
		      CL_B_gm*(ne*clCoolLineH2_e(Rate->T)*Y->H2 + clCoolLineH2_HI(Rate->T)*en_B*Y->H2*Y->HI + clCoolLineH2_H2(Rate->T)*en_B*Y->H2*Y->H2 + clCoolLineH2_He(Rate->T)*en_B*Y->H2*Y->HeI + clCoolLineH2_HII(Rate->T)*en_B*Y->H2*Y->HII) + 
		      cl->R.Cool_Coll_H2*Y->H2*Rate->Coll_e_H2 + cl->R.Cool_Coll_H2*Y->H2*Rate->Coll_HI_H2*Y->HI*en_B + cl->R.Cool_Coll_H2*Y->H2*Rate->Coll_H2_H2*Y->H2*en_B +/*removing shielding of collisions 9/22/10 CC*/
		      LowTCool + Rate->Cool_Metal),
	     Y->HI*cl->R.Heat_Phot_HI*Rate->Phot_HI + 
                       Y->HeI*cl->R.Heat_Phot_HeI*Rate->Phot_HeI + 
                       Y->HeII*cl->R.Heat_Phot_HeII*Rate->Phot_HeII,
	     Y->H2*cl->R.Heat_Phot_H2*Rate->Phot_H2*s_dust*s_self,
	     Rate->Heat_Metal);
    fclose(cooldebug);
  }
#endif
 }

/* Module Interface routines */

void CoolAddParams( COOLPARAM *CoolParam, PRM prm ) {
	CoolParam->bIonNonEqm = 1;
	prmAddParam(prm,"bIonNonEqm",paramBool,&CoolParam->bIonNonEqm,
				sizeof(int),"IonNonEqm",
				"<Gas is Cooling Non-Equilibrium Abundances> = +IonNonEqm");
	CoolParam->bUV = 1;
	prmAddParam(prm,"bUV",paramBool,&CoolParam->bUV,sizeof(int),"UV",
				"read in an Ultra Violet file = +UV");
	CoolParam->bUVTableUsesTime = 0;
	prmAddParam(prm,"bUVTableUsesTime",paramBool,&CoolParam->bUVTableUsesTime,sizeof(int),"UVTableUsesTime",
				"Ultra Violet Table uses time = +UVTableUsesTime");
	CoolParam->dMassFracHelium = 0.236;
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
	CoolParam->bMetal = 1; 
	prmAddParam(prm,"bMetal",paramBool,&CoolParam->bMetal,sizeof(int),
				"mtc","enable/disable Metal heating/cooling = +mtc");
	CoolParam->bShieldHI = 1; /*Whether or not to shield HI by dust*/
	prmAddParam(prm,"bShieldHI",paramBool,&CoolParam->bShieldHI,sizeof(int),
				"shI","enable/disable Dust Shielding of HI = +shI");
	CoolParam->dClump = 10;
	prmAddParam(prm,"dClump",paramDouble,&CoolParam->dClump, /* Sub-grid clumping factor used to calculate H2 formation on dust */
				sizeof(double),"clump",
				"<Clumping Factor for H2 Formation>  = 10");
	CoolParam->dLymanWernerFrac = 0.5; /* Multipler for the H2 photo-dissociation by stellar LW photons (i.e., what fraction of the LW radiation escapes the birth cloud) */
	prmAddParam(prm,"dLymanWernerFrac",paramDouble,&CoolParam->dLymanWernerFrac,sizeof(double),
				"lwf","coeff times stellar lyman-werner photons");
	/*	strcpy(CoolParam->CoolInFile,"cooltable_xdr\0");
	prmAddParam(prm,"CoolInFile",3,&CoolParam->CoolInFile,256,"coolin",
	"<cooling table file> (file in xdr binary format)"); */

	}
	
void CoolOutputArray( COOLPARAM *CoolParam, int cnt, int *type, char *suffix ) {
 #if 0
	*type = OUT_NULL;

	switch (cnt) {
	case 0:
		if (!CoolParam->bDoIonOutput) return;
		*type = OUT_COOL_ARRAY0;
		sprintf(suffix,".HI");
		return;
	case 1:
		if (!CoolParam->bDoIonOutput) return;
		*type = OUT_COOL_ARRAY1;
		sprintf(suffix,".HeI");
		return;
	case 2:
		if (!CoolParam->bDoIonOutput) return;
		*type = OUT_COOL_ARRAY2;
		sprintf(suffix,".HeII");
		return;
        case 3:
		if (!CoolParam->bDoIonOutput) return;
		*type = OUT_COOL_ARRAY3;
		sprintf(suffix,".H2");
		return;
	}
#endif
}

/* Output Conversion Routines */
double CoolEnergyToTemperature( COOL *cl, COOLPARTICLE *cp, double E, double ZMetal ) {
    double Y_H, Y_He, Y_eMax;
    clSetAbundanceTotals(cl,ZMetal,&Y_H,&Y_He,&Y_eMax);
    return clTemperature(2*Y_H - cp->f_HI*Y_H 
			 -  cp->f_H2*Y_H - (cp->f_H2*Y_H)/2.0  /*Total - electrons in H2 - half the number of atoms in H2 */
			 + 3*Y_He - 2*cp->f_HeI*Y_He - cp->f_HeII*Y_He + ZMetal/MU_METAL, E );
    }

double CoolCodeEnergyToTemperature( COOL *cl, COOLPARTICLE *cp, double E, double ZMetal ) {
    return CoolEnergyToTemperature( cl, cp, E*cl->dErgPerGmUnit, ZMetal );
    }

/* Initialization Routines */

void CoolTableReadInfo( COOLPARAM *CoolParam, int cntTable, int *nTableColumns, char *suffix )
{
   int localcntTable = 0;


   *nTableColumns = 0;
#if (0)
   if (CoolParam->bMetal) { 
       /* Use self-consistent UV from metal data file for metal cooling */
       return; 
       }
   /* Overwrite UV in the metal data file */
   if (!CoolParam->bMetal && CoolParam->bUV) {
	   if (localcntTable == cntTable) {
		   *nTableColumns = 7;
		   sprintf(suffix,"UV");
		   return;
		   }
	   localcntTable++;
	   }
#endif
   }

void CoolTableRead( COOL *Cool, int nData, void *vData)
{
   double *dTableData = vData;
   int nTableColumns; 
   int nTableRows;
   int localcntTable = 0;
   printf("calling CoolTableRead...\n");
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
   
   fprintf(stderr," Attempt to initialize non-existent table in cooling\n");
   assert(0);

}

void CoolDefaultParticleData( COOLPARTICLE *cp )
{
 
	cp->f_HI = 0.75;
	cp->f_HeI = 0.06;
	cp->f_HeII = 0.0;
	cp->f_H2 = 0.0;
}

void CoolInitEnergyAndParticleData( COOL *cl, COOLPARTICLE *cp, double *E, double dDensity, double dTemp, double ZMetal)
{
	PERBARYON Y;
	RATE r;
	double  columnL = 1e20/cl->dKpcUnit/3.08568025e21/cl->dExpand; /*assumption to get things started before we know the energy to calculate the sound speed)*/
	
	cp->f_HI = 1.0;
	cp->f_HeI = 1.0;
	cp->f_HeII = 0.0;
	cp->f_H2 = 0.0;
	double dLymanWerner = cp->dLymanWerner;
	CoolPARTICLEtoPERBARYON(cl, &Y, cp, ZMetal);
	CLRATES(cl,&r,dTemp,CodeDensityToComovingGmPerCc(cl,dDensity), ZMetal, columnL, dLymanWerner);
	clAbunds(cl,&Y,&r,CodeDensityToComovingGmPerCc(cl,dDensity), ZMetal);
	CoolPERBARYONtoPARTICLE(cl, &Y,cp, ZMetal);
	*E = clThermalEnergy(Y.Total,dTemp)*cl->diErgPerGmUnit;

#ifdef COOLDEBUG
       if (cl->p->iOrder == 1) {
	FILE *cooldebug;
	/* Output heading for table of coooling of the gas particle*/
	cooldebug = fopen("cooldebug_table.txt","w");
	/*	fprintf(cooldebug,"%15s",T [K]     energy        Density   e            H2  HI           HII     HeI    HeII         HeIII    shield        Edot            Cool_Comp      Brem           Diel          Radr           LineHI   LineHeI   LineHeII  LineH2_H       LineH2_H2       LineH2_He     LineH2_e       LineH2_HII    Coll          CollH2         LowT            MCool         Cool           Phot           PhotH2       MHeat);*/
	fprintf(cooldebug,"%15s, %15s, %15s, %15s, %15s, %15s, %15s, %15s, %15s, %15s, %15s, %15s, %15s, %15s, %15s, %15s, %15s, %15s, %15s, %15s, %15s, %15s, %15s, %15s, %15s, %15s, %15s, %15s, %15s, %15s, %15s, %15s, %15s, %15s\n",
		"T [K]","Energy","Density","e","H2","HI","HII","HeI","HeII","HeIII","shield","Edot","Cool_Comp","Brem","Diel","Radr","LineHI","LineHeI","LineHeII","LineH2_H","LineH2_H2","LineH2_He","LineH2_e","LineH2_HII","Coll","CollH2_e","CollH2_H","CollH2_H2","LowT","MCool","Cool","Phot","PhotH2","MHeat");
	fclose(cooldebug);
       }

       if (cl->p->iOrder == 1) {
	FILE *fpdebug;
	fpdebug = fopen("cooldebug.txt","w");
	fprintf(fpdebug,"\n");
	fclose(fpdebug);
       }
#endif
}

void CoolInitRatesTable( COOL *cl, COOLPARAM CoolParam ) {
	clInitRatesTable( cl, CoolParam.dCoolingTmin, CoolParam.dCoolingTmax, CoolParam.nCoolingTable );
	clReadMetalTable(cl, CoolParam);
	
}

void CoolSetTime( COOL *cl, double dTime, double z ) {
	clRatesRedshift( cl, z, dTime );
}

/* ***************************************** Integration Routines ******************************* */

/**
 * Physical units
 */
void CoolIntegrateEnergy(COOL *cl, clDerivsData *clData, COOLPARTICLE *cp, 
			 double *E,
			 double ExternalHeating, double rho, double ZMetal, double tStep, double columnL  ) {
        PERBARYON Y;
	double dLymanWerner = cp->dLymanWerner;
        CoolPARTICLEtoPERBARYON(cl, &Y,cp, ZMetal);
        clIntegrateEnergy(cl, clData, &Y, E, ExternalHeating, rho, ZMetal, tStep, columnL, dLymanWerner);
        CoolPERBARYONtoPARTICLE(cl, &Y, cp, ZMetal);
        }

/**
 * Code Units except for dt
 */
void CoolIntegrateEnergyCode(COOL *cl, clDerivsData *clData, COOLPARTICLE *cp, 
			     double *ECode, double ExternalHeatingCode, 
			     double rhoCode, double ZMetal, double *posCode, 
			     double tStep, double columnL ) {
	PERBARYON Y;
	double T,E,rho;

	double dLymanWerner = cp->dLymanWerner;

	E = CoolCodeEnergyToErgPerGm( cl, *ECode );
	T = CoolEnergyToTemperature( cl, cp, E, ZMetal);
	rho = CodeDensityToComovingGmPerCc(cl,rhoCode );
	CoolPARTICLEtoPERBARYON(cl, &Y, cp, ZMetal);
#ifdef NOEXTHEAT
	ExternalHeatingCode = 0;
#endif
	clIntegrateEnergy(cl, clData, &Y, &E,
			  CoolCodeWorkToErgPerGmPerSec( cl, ExternalHeatingCode ), 
			  CodeDensityToComovingGmPerCc(cl, rhoCode), ZMetal, tStep, columnL,dLymanWerner );
	CoolPERBARYONtoPARTICLE(cl, &Y, cp, ZMetal);
	*ECode = CoolErgPerGmToCodeEnergy(cl, E);
	}

/* Deprecated: */
double CoolHeatingRate( COOL *cl, COOLPARTICLE *cp, double T, double dDensity, double ZMetal, double columnL ) {
    PERBARYON Y;
    RATE Rate;
    double Y_H, Y_He, Y_eMax;
    double dLymanWerner = cp->dLymanWerner;
    clSetAbundanceTotals(cl,ZMetal,&Y_H,&Y_He,&Y_eMax);
    CoolPARTICLEtoPERBARYON(cl, &Y, cp, ZMetal);
    CLRATES(cl, &Rate, T, dDensity, ZMetal, columnL, dLymanWerner);
    clRateMetalTable(cl, &Rate, T, dDensity, Y_H, ZMetal); 
    return (-clCoolTotal(cl, &Y, &Rate, dDensity, ZMetal ) + clHeatTotal(cl, &Y, &Rate, dDensity, ZMetal));
    }

/* Code heating - cooling rate excluding external heating (PdV, etc..) */
double CoolEdotInstantCode(COOL *cl, COOLPARTICLE *cp, double ECode, 
			     double rhoCode, double ZMetal, double *posCode, double columnL ) {
    PERBARYON Y;
    RATE Rate;
    double T,E,rho,Edot;
    double Y_H, Y_He, Y_eMax;
    double dHeat, dCool;
    double dLymanWerner = cp->dLymanWerner;
    clSetAbundanceTotals(cl,ZMetal,&Y_H,&Y_He,&Y_eMax);

    E = CoolCodeEnergyToErgPerGm( cl, ECode );
    T = CoolEnergyToTemperature( cl, cp, E, ZMetal);
    rho = CodeDensityToComovingGmPerCc(cl,rhoCode );
    CoolPARTICLEtoPERBARYON(cl, &Y, cp, ZMetal);
    CLRATES(cl, &Rate, T, rho, ZMetal, columnL, dLymanWerner);
    clRateMetalTable(cl, &Rate, T, rho, Y_H, ZMetal); 
    Edot = CLEDOTINSTANT( cl, &Y, &Rate, rho, ZMetal, &dHeat, &dCool );

    return CoolErgPerGmPerSecToCodeWork( cl, Edot );
    }

/* Code heating - cooling rate excluding external heating (PdV, etc..) */
double CoolCoolingCode(COOL *cl, COOLPARTICLE *cp, double ECode, 
		       double rhoCode, double ZMetal, double *posCode, double columnL ) {
    PERBARYON Y;
    RATE Rate;
    double T,E,rho,Edot;
    double Y_H, Y_He, Y_eMax;
    double dLymanWerner = cp->dLymanWerner;
    clSetAbundanceTotals(cl,ZMetal,&Y_H,&Y_He,&Y_eMax);

    E = CoolCodeEnergyToErgPerGm( cl, ECode );
    T = CoolEnergyToTemperature( cl, cp, E, ZMetal );
    rho = CodeDensityToComovingGmPerCc(cl,rhoCode );
    CoolPARTICLEtoPERBARYON(cl, &Y, cp, ZMetal);
    CLRATES(cl, &Rate, T, rho, ZMetal, columnL, dLymanWerner);
    clRateMetalTable(cl, &Rate, T, rho, Y_H, ZMetal);

    Edot = clCoolTotal(cl, &Y, &Rate, rho, ZMetal );

    return CoolErgPerGmPerSecToCodeWork( cl, Edot );
    }

/* Code heating due to atomic/radiative processes only */
double CoolHeatingCode(COOL *cl, COOLPARTICLE *cp, double ECode, 
		       double rhoCode, double ZMetal, double *posCode, double columnL ) {
    PERBARYON Y;
    RATE Rate;
    double T,E,rho,Edot;
    double Y_H, Y_He, Y_eMax;
    double dLymanWerner = cp->dLymanWerner;
    clSetAbundanceTotals(cl,ZMetal,&Y_H,&Y_He,&Y_eMax);

    E = CoolCodeEnergyToErgPerGm( cl, ECode );
    T = CoolEnergyToTemperature( cl, cp, E, ZMetal );
    rho = CodeDensityToComovingGmPerCc(cl,rhoCode );
    CoolPARTICLEtoPERBARYON(cl, &Y, cp, ZMetal);
    CLRATES(cl, &Rate, T, rho, ZMetal, columnL, dLymanWerner);
    clRateMetalTable(cl, &Rate, T, rho, Y_H, ZMetal);

    Edot = clHeatTotal ( cl, &Y, &Rate, rho, ZMetal );

    return CoolErgPerGmPerSecToCodeWork( cl, Edot );
    }

#endif /* NOCOOLING */
