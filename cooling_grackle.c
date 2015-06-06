#ifndef NOCOOLING

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
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

#include "cooling.h"

#define mh     1.67262171e-24   
#define kboltz 1.3806504e-16

COOL *CoolInit( )
    {
    COOL *cl;
    cl = (COOL *) malloc(sizeof(COOL));
    assert(cl!=NULL);
#ifdef CONFIG_BFLOAT_8
    assert(sizeof(gr_float)==8);
#else
#ifdef CONFIG_BFLOAT_4
    assert(sizeof(gr_float)==4);
#else
    fprintf(stderr,"Cooling Grackle: gr_float type not defined\n");
    assert(0);
#endif
#endif
    
    cl->pgrackle_data = &grackle_data;
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
    return Data;
    }

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

void clInitConstants( COOL *cl, double dGmPerCcUnit, double dComovingGmPerCcUnit, 
		      double dErgPerGmUnit, double dSecUnit, double dKpcUnit, COOLPARAM CoolParam) 
    {
    assert(cl!=NULL);
	grackle_verbose = CoolParam.grackle_verbose;

    cl->my_units.comoving_coordinates = CoolParam.bComoving; // 1 if cosmological sim, 0 if not
    cl->my_units.density_units = dGmPerCcUnit;
    cl->my_units.length_units = dKpcUnit*3.0857e21; // cm
    cl->my_units.time_units = dSecUnit;
    cl->my_units.velocity_units = cl->my_units.length_units / cl->my_units.time_units;
    cl->my_units.a_units = 1.0; // units for the expansion factor

    /* Erg per Gm unit is calculated as velocity units^2 */
    cl->dErgPerGmUnit = dErgPerGmUnit; // too useful not to keep
    cl->diErgPerGmUnit = 1/dErgPerGmUnit;
    cl->dSecUnit = dSecUnit;
    cl->dComovingGmPerCcUnit = dComovingGmPerCcUnit;
    cl->dErgPerGmPerSecUnit = dErgPerGmUnit/dSecUnit;

    // Second, create a chemistry object for parameters and rate data.
    if (set_default_chemistry_parameters() == 0) {
        fprintf(stderr, "Grackle Error in set_default_chemistry_parameters.\n");
        assert(0);
        }
    
    cl->pgrackle_data->use_grackle = CoolParam.use_grackle;            // chemistry on
    cl->pgrackle_data->with_radiative_cooling = CoolParam.with_radiative_cooling; // cooling on
    cl->pgrackle_data->primordial_chemistry = CoolParam.primordial_chemistry;   // 0-3 molecular network with H, He, D
    if (CoolParam.primordial_chemistry > GRACKLE_PRIMORDIAL_CHEMISTRY_MAX) {
        fprintf(stderr,"Must compile so that GRACKLE_PRIMORDIAL_CHEMISTRY_MAX >= primordial_chemistry parameter used\n");
        assert(GRACKLE_PRIMORDIAL_CHEMISTRY_MAX >= CoolParam.primordial_chemistry);
        };
    cl->pgrackle_data->metal_cooling = CoolParam.metal_cooling;          // metal cooling on
    cl->pgrackle_data->UVbackground = CoolParam.UVbackground;           // UV background on
    strncpy( cl->grackle_data_file, CoolParam.grackle_data_file, MAXPATHLEN ); // Permanent local copy
    cl->pgrackle_data->grackle_data_file = cl->grackle_data_file; // hdf5 cloudy data file (pointer to permanent copy)
    
        {
        double initial_redshift = 0.;
        double a_value = 1. / (1. + initial_redshift);
        
        // Finally, initialize the chemistry object.
        if (initialize_chemistry_data(&cl->my_units, a_value) == 0) {
            fprintf(stderr, "Grackle Error in initialize_chemistry_data.\n");
            assert(0);
            }
        }


    }

/*Returns baryonic fraction for a given species*/
double COOL_ARRAY0(COOL *cl, COOLPARTICLE *cp, double ZMetal) {
#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX>=1)
    return (cp->HI);
#else
    return 0;
#endif
}

double COOL_ARRAY1(COOL *cl, COOLPARTICLE *cp, double ZMetal) {
#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX>=1)
    return (cp->HII);
#else
    return 0;
#endif
}

double COOL_ARRAY2(COOL *cl, COOLPARTICLE *cp, double ZMetal) {
#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX>=1)
    return (cp->HeI);
#else
    return 0;
#endif
}

double COOL_ARRAY3(COOL *cl, COOLPARTICLE *cp, double ZMetal) {
#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX>=1)
    return (cp->HeII);
#else
    return 0;
#endif
}

double COOL_ARRAY4(COOL *cl, COOLPARTICLE *cp, double ZMetal) {
#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX>=1)
    return (cp->HeIII);
#else
    return 0;
#endif
}

double COOL_ARRAY5(COOL *cl, COOLPARTICLE *cp, double ZMetal) {
#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX>=1)
    return (cp->e);
#else
    return 0;
#endif
}

double COOL_ARRAY6(COOL *cl, COOLPARTICLE *cp, double ZMetal) {
#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX>=2)
    return (cp->HM);
#else
    return 0;
#endif
}

double COOL_ARRAY7(COOL *cl, COOLPARTICLE *cp, double ZMetal) {
#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX>=2)
    return (cp->H2I);
#else
    return 0;
#endif
}

double COOL_ARRAY8(COOL *cl, COOLPARTICLE *cp, double ZMetal) {
#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX>=2)
    return (cp->H2II);
#else
    return 0;
#endif
}

double COOL_ARRAY9(COOL *cl, COOLPARTICLE *cp, double ZMetal) {
#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX>=3)
    return (cp->DI);
#else
    return 0;
#endif
}

double COOL_ARRAY10(COOL *cl, COOLPARTICLE *cp, double ZMetal) {
#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX>=3)
    return (cp->DII);
#else
    return 0;
#endif
}

double COOL_ARRAY11(COOL *cl, COOLPARTICLE *cp, double ZMetal) {
#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX>=3)
    return (cp->HDI);
#else
    return 0;
#endif
}


void CoolTableReadInfo( COOLPARAM *CoolParam, int cntTable, int *nTableColumns, char *suffix ) {
    *nTableColumns = 0;
    }


void CoolTableRead( COOL *Cool, int nData, void *vData)  {
    assert(0);
}

void CoolInitRatesTable( COOL *cl, COOLPARAM CoolParam ) {
    } 

void CoolSetTime( COOL *cl, double dTime, double z ) {
    double a_value = 1. / (1. + z);

    cl->dTime = dTime;
    cl->z = z;
    cl->a = a_value;
/*
    if (initialize_chemistry_data(&cl->my_units, a_value) == 0) {
        fprintf(stderr, "Grackle Error in initialize_chemistry_data.\n");
        assert(0);
        }
*/
    }


void CoolDefaultParticleData( COOLPARTICLE *cp ) 
{
    /* Never used I think - just to set values */
    double tiny_number = 1.e-20;
#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX>=1)
//    gr_float HI, HII, HeI, HeII, HeII, e;
    cp->HI = 0.76;
    cp->HII = tiny_number;
    cp->HeI = 0.24;
    cp->HeII = tiny_number;
    cp->HeIII = tiny_number;
    cp->e = tiny_number;
#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX>=2)
//    gr_float HM, H2I, H2II;
    cp->HM = tiny_number;
    cp->H2I = tiny_number;
    cp->H2II = tiny_number;
#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX>=3)
//    gr_float DI, DII, HDI
    cp->DI = 2.0 * 3.4e-5;
    cp->DII = tiny_number;
    cp->HDI = tiny_number;
#endif
#endif
#endif
}

void CoolInitEnergyAndParticleData( COOL *cl, COOLPARTICLE *cp, double *E, double dDensity, double dTemp, double ZMetal) {
    double tiny_number = 1.e-20;

/* Ionization fractions arbitrary -- should be set to eqm or some early universe model (e ~ 1e-5) */
#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX>=1)
//    gr_float HI, HII, HeI, HeII, HeII, e;
    double fNonMetal = 1-ZMetal; // see make_consistent_g in solve_rate_cool.F
    cp->HI = fNonMetal * cl->pgrackle_data->HydrogenFractionByMass;
    cp->HII = tiny_number;
    cp->HeI = fNonMetal * (1.0 - cl->pgrackle_data->HydrogenFractionByMass); //(from c_example.c -- fails to do He increase w/ Z)
    cp->HeII = tiny_number;
    cp->HeIII = tiny_number;
    cp->e = tiny_number;
#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX>=2)
//    gr_float HM, H2I, H2II;
    cp->HM = tiny_number * dDensity;
    cp->H2I = tiny_number * dDensity;
    cp->H2II = tiny_number * dDensity;
#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX>=3)
//    gr_float DI, DII, HDI
    cp->DI = 2.0 * 3.4e-5 * dDensity;
    cp->DII = tiny_number * dDensity;
    cp->HDI = tiny_number * dDensity;
#endif
#endif
#endif
    // solar metallicity
//    metal_density[i] = grackle_data.SolarMetalFractionByMass * density[i];

    // No routine in Grackle to use to set up energy sensibly
    // Energy: erg per gm --> code units   assumes neutral gas (as above)
    // This is not called after checkpoints so should be ok
    *E = dTemp*(kboltz*1.5/(1.2*mh))*cl->diErgPerGmUnit;

//    printf("Grackle %d \n",GRACKLE_PRIMORDIAL_CHEMISTRY_MAX);
//    printf("%g %g   %g %g %g %g %g %g %g\n",dDensity,*E,cp->HI,cp->HII,cp->HeI,cp->HeII,cp->HeIII,cp->e);

    double Tnew=0;
    int its=0;
    while (its++ < 20 && fabs(Tnew/dTemp-1) > 1e-2) {
        Tnew = CoolCodeEnergyToTemperature( cl, cp, *E, dDensity, ZMetal );
        *E = *E*dTemp/Tnew; 
        }
    assert(its<20);
//    printf("dTemp %g (K) Tnew %g (K) Energy %g (erg/g)  mu %g \n",dTemp,Tnew,*E*cl->dErgPerGmUnit,Tnew/dTemp*1.2);
    }


void CoolIntegrateEnergy(COOL *cl, clDerivsData *clData, COOLPARTICLE *cp, double *E,
			 double ExternalHeat, double rho, double ZMetal, double *rp, double tStep ) {
    assert(0);
    }

void CoolIntegrateEnergyCode(COOL *cl, clDerivsData *clData, COOLPARTICLE *cp, double *E,
        double ExternalHeat, double rho, double ZMetal, double *rp, double tStep ) {

    double dt = tStep/cl->dSecUnit; // back into code units
    int zero[]={0,0,0},one[]={1,1,1};
    gr_float density = rho, energy = *E,
        x_velocity=0, y_velocity=0, z_velocity=0,
        HI_density, HII_density, HM_density,
        HeI_density, HeII_density, HeIII_density,
        H2I_density, H2II_density,
        DI_density, DII_density, HDI_density,
        e_density, metal_density;

    metal_density = ZMetal*density;

    *E += 0.5*ExternalHeat*dt;  /* Gnedin suggestion */

    if (cl->pgrackle_data->primordial_chemistry==0) {
/*
    int solve_chemistry_table(code_units *my_units, double a_value, double dt_value, 
        int grid_rank, int *grid_dimension, int *grid_start, int *grid_end, 
        gr_float *density, gr_float *internal_energy, 
        gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity, gr_float *metal_density); 
*/
  
        if (solve_chemistry_table(&(cl->my_units), cl->a, dt,
            1, one, zero, zero,
            &density, &energy,
            &x_velocity, &y_velocity, &z_velocity,
            &metal_density)== 0) {
            fprintf(stderr, "Grackle Error in solve_chemistry.\n");
            assert(0);
            }
        }
    else {
#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX>=1)
        HI_density = cp->HI*density;
        HII_density = cp->HII*density;
        HeI_density = cp->HeI*density;
        HeII_density = cp->HeII*density;
        HeIII_density = cp->HeIII*density;
        e_density = cp->e*density;
#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX>=2)
        HM_density = cp->HM*density;
        H2I_density = cp->H2I*density;
        H2II_density = cp->H2II*density;
#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX>=3)
        DI_density = cp->DI*density;
        DII_density = cp->DII*density;
        HDI_density = cp->HDI*density;
#endif
#endif
#endif
/*
  solve_chemistry(&(cl->my_units),
        a_value, dt,
        grid_rank, grid_dimension,
        grid_start, grid_end,
        density, energy,
        x_velocity, y_velocity, z_velocity,
        HI_density, HII_density, HM_density,
        HeI_density, HeII_density, HeIII_density,
        H2I_density, H2II_density,
        DI_density, DII_density, HDI_density,
        e_density, metal_density) */
        if (solve_chemistry(&(cl->my_units),
            cl->a, dt,
            1, one, zero, zero,
            &density, &energy,
            &x_velocity, &y_velocity, &z_velocity,
            &HI_density, &HII_density, &HM_density,
            &HeI_density, &HeII_density, &HeIII_density,
            &H2I_density, &H2II_density,
            &DI_density, &DII_density, &HDI_density,
                &e_density, &metal_density)== 0) {
    
            fprintf(stderr, "Grackle Error in solve_chemistry.\n");
            assert(0);
            }
        }
    
    *E = energy;
    *E += 0.5*ExternalHeat*dt;  /* Gnedin suggestion -- far from self-consistent */

#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX>=1)
    float dinv = 1./density;
    cp->HI = HI_density*dinv;
    cp->HII = HII_density*dinv;
    cp->HeI = HeI_density*dinv;
    cp->HeII = HeII_density*dinv;
    cp->e = e_density*dinv;
#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX>=2)
    cp->HM = HM_density*dinv;
    cp->H2I = H2I_density*dinv;
    cp->H2II = H2II_density*dinv;
#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX>=3)
    cp->DI = DI_density*dinv;
    cp->DII = DII_density*dinv;
    cp->HDI = HDI_density*dinv;
#endif
#endif
#endif

    }
    
double CoolHeatingRate( COOL *cl, COOLPARTICLE *cp, double T, double dDensity, double ZMetal, double rkpc) {
    assert(0);
    }

double CoolCoolingCode(COOL *cl, COOLPARTICLE *cp, double ECode, 
		       double rhoCode, double ZMetal, double *posCode ) {
    assert(0);
    }

double CoolHeatingCode(COOL *cl, COOLPARTICLE *cp, double ECode, 
		       double rhoCode, double ZMetal, double *posCode ) {
    assert(0);
    }

void CoolAddParams( COOLPARAM *CoolParam, PRM prm ) {
	CoolParam->bDoIonOutput = 1;
	prmAddParam(prm,"bDoIonOutput",paramBool,&CoolParam->bDoIonOutput,sizeof(int),
				"Iout","enable/disable Ion outputs (cooling only) = +Iout");


	CoolParam->grackle_verbose = 0;
	prmAddParam(prm,"grackle_verbose",paramBool,&CoolParam->grackle_verbose,sizeof(int),"grackle_verbose",
				"on =  +grackle_verbose [off]");
	CoolParam->use_grackle = 1;
	prmAddParam(prm,"use_grackle",paramBool,&CoolParam->use_grackle,sizeof(int),"use_grackle",
				"on = +use_grackle");
	CoolParam->with_radiative_cooling = 1;
	prmAddParam(prm,"with_radiative_cooling",paramBool,&CoolParam->with_radiative_cooling,sizeof(int),"with_radiative_cooling",
				"on = +with_radiative_cooling");
	CoolParam->primordial_chemistry = 0;
	prmAddParam(prm,"primordial_chemistry",paramInt,&CoolParam->primordial_chemistry,sizeof(int),"primordial_chemistry",
				"-primordial_chemistry=0 [values 0,1,2,3]");
	CoolParam->metal_cooling = 1;
	prmAddParam(prm,"metal_cooling",paramBool,&CoolParam->metal_cooling,sizeof(int),"metal_cooling",
				"on = +metal_cooling");
	CoolParam->UVbackground = 1;
	prmAddParam(prm,"UVbackground",paramBool,&CoolParam->UVbackground,sizeof(int),"UVbackground",
				"on = +UVbackground");
	CoolParam->bComoving = 1;
	prmAddParam(prm,"bComoving",paramBool,&CoolParam->bComoving,sizeof(int),"bComoving",
				"on = +bComoving");
	strcpy(CoolParam->grackle_data_file,"CloudyData_UVB=HM2012.h5\0");
	prmAddParam(prm,"grackle_data_file",paramString,&CoolParam->grackle_data_file,256,"grackle_data_file",
	"<cooling table file> (file in hdf5 format, e.g. CloudyData_UVB=HM2012.h5)"); 

	}
	

#if 0
/* Not needed for ChaNGa */
void CoolLogParams( COOLPARAM *CoolParam, LOGGER *lgr) {
    LogParams(lgr, "COOLING", "bDoIonOutput: %d",CoolParam->bDoIonOutput); 

    LogParams(lgr, "COOLING", "grackle_verbose: %d",CoolParam->grackle_verbose); 
    LogParams(lgr, "COOLING", "use_grackle: %d",CoolParam->use_grackle); 
    LogParams(lgr, "COOLING", "with_radiative_cooling: %d",CoolParam->with_radiative_cooling);
    LogParams(lgr, "COOLING", "primorial_chemistry: %d (MAX %d) ",CoolParam->primordial_chemistry,GRACKLE_PRIMORDIAL_CHEMISTRY_MAX); 
    LogParams(lgr, "COOLING", "metal_cooling: %d",CoolParam->metal_cooling); 
    LogParams(lgr, "COOLING", "UVbackground: %d",CoolParam->UVbackground); 
    LogParams(lgr, "COOLING", "bComoving: %d",CoolParam->bComoving); 
    LogParams(lgr, "COOLING", "grackle_data_file: %s",CoolParam->grackle_data_file); 

}
#endif

void CoolOutputArray( COOLPARAM *CoolParam, int cnt, int *type, char *suffix ) {
#if 0
    char *extensions[]= { "HI", "HII", "HeI", "HeII", "HeIII", "e",
                          "HM", "H2I", "H2II",
                          "DI", "DII", "HDI" };
	*type = OUT_NULL;
    if (!CoolParam->bDoIonOutput) return;

/*
#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX<1)
    float dummy;
#endif
#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX>=1)
    gr_float HI, HII, HeI, HeII, HeII, e;
#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX>=2)
    gr_float HM, H2I, H2II;
#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX>=3)
    gr_float DI, DII, HDI
#endif
#endif
#endif
*/
	switch (cnt) {
	case 11:
    case 10:
    case 9:
		if (CoolParam->primordial_chemistry<3) return;
    case 8:
    case 7:
    case 6:
		if (CoolParam->primordial_chemistry<2) return;
    case 5:
    case 4:
    case 3:
    case 2:
    case 1:
    case 0:
		if (CoolParam->primordial_chemistry<1) return;
		*type = OUT_COOL_ARRAY0+cnt;
		sprintf(suffix,extensions[cnt]);
		return;
	}
#endif
}

double CoolCodeEnergyToTemperature( COOL *cl, COOLPARTICLE *cp, double E, double rho, double ZMetal ) {
    int one[]={1,1,1};
    gr_float density = rho, energy = E,
        x_velocity=0, y_velocity=0, z_velocity=0,
        HI_density, HII_density, HM_density,
        HeI_density, HeII_density, HeIII_density,
        H2I_density, H2II_density,
        DI_density, DII_density, HDI_density,
        e_density, metal_density, temperature;

    metal_density = ZMetal*density;

    if (cl->pgrackle_data->primordial_chemistry==0) {
/*
    calculate_temperature_table(code_units *my_units, int grid_rank, int *grid_dimension, 
    gr_float *density, gr_float *internal_energy, gr_float *metal_density, gr_float *temperature);*/
        if (calculate_temperature_table(&cl->my_units,
                1, one,
                &density, &energy, &metal_density, &temperature) == 0) {
            fprintf(stderr, "Grackle Error in calculate_temperature.\n");
            assert(0);
            }
        }
    else {
#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX>=1)
        HI_density = cp->HI*density;
        HII_density = cp->HII*density;
        HeI_density = cp->HeI*density;
        HeII_density = cp->HeII*density;
        HeIII_density = cp->HeIII*density;
        e_density = cp->e*density;
#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX>=2)
        HM_density = cp->HM*density;
        H2I_density = cp->H2I*density;
        H2II_density = cp->H2II*density;
#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX>=3)
        DI_density = cp->DI*density;
        DII_density = cp->DII*density;
        HDI_density = cp->HDI*density;
#endif
#endif
#endif
        if (calculate_temperature(&cl->my_units,
                1, one,
                &density, &energy,
                &HI_density, &HII_density, &HM_density,
                &HeI_density, &HeII_density, &HeIII_density,
                &H2I_density, &H2II_density,
                &DI_density, &DII_density, &HDI_density,
                &e_density, &metal_density, 
                &temperature) == 0) {
            fprintf(stderr, "Grackle Error in calculate_temperature.\n");
            assert(0);
            }
        }
    return ((double) temperature);
    }

// Currently a MACRO -- assumes gamma=5/3
// Should use: Grackle calculate_pressure if H2 > 0
/* 
void CoolCodePressureOnDensitySoundSpeed( COOL *cl, COOLPARTICLE *cp, double uPred, double fDensity, double gamma, double gammam1, double *PoverRho, double *c ) {


}
*/

/* Code heating - cooling rate excluding external heating (PdV, etc..) */
double CoolEdotInstantCode(COOL *cl, COOLPARTICLE *cp, double ECode, 
			     double rhoCode, double ZMetal, double *posCode ) {
    int zero[]={0,0,0},one[]={1,1,1};
    gr_float density = rhoCode, energy = ECode,
        x_velocity=0, y_velocity=0, z_velocity=0,
        HI_density, HII_density, HM_density,
        HeI_density, HeII_density, HeIII_density,
        H2I_density, H2II_density,
        DI_density, DII_density, HDI_density,
        e_density, metal_density, cooling_time;

    metal_density = ZMetal*density;

    if (cl->pgrackle_data->primordial_chemistry==0) {
/*
        int calculate_cooling_time_table(code_units *my_units, double a_value, 
            int grid_rank, int *grid_dimension, int *grid_start, int *grid_end, 
            gr_float *density, gr_float *internal_energy, 
            gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity, gr_float *metal_density, gr_float *cooling_time); */
        if (calculate_cooling_time_table(&cl->my_units, cl->a,
                1, one, zero, zero,
                &density, &energy, 
                &x_velocity, &y_velocity, &z_velocity, &metal_density, &cooling_time) == 0) {
            fprintf(stderr, "Grackle Error in calculate_cooling_time_table.\n");
            assert(0);
            }
        }
    else {
#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX>=1)
        HI_density = cp->HI*density;
        HII_density = cp->HII*density;
        HeI_density = cp->HeI*density;
        HeII_density = cp->HeII*density;
        HeIII_density = cp->HeIII*density;
        e_density = cp->e*density;
#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX>=2)
        HM_density = cp->HM*density;
        H2I_density = cp->H2I*density;
        H2II_density = cp->H2II*density;
#if (GRACKLE_PRIMORDIAL_CHEMISTRY_MAX>=3)
        DI_density = cp->DI*density;
        DII_density = cp->DII*density;
        HDI_density = cp->HDI*density;
#endif
#endif
#endif
/*
        calculate_cooling_time(code_units *my_units, double a_value, 
            int grid_rank, int *grid_dimension, int *grid_start, int *grid_end, 
            gr_float *density, gr_float *internal_energy, 
            gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity, 
            gr_float *HI_density, gr_float *HII_density, gr_float *HM_density, 
            gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density, 
            gr_float *H2I_density, gr_float *H2II_density, gr_float *DI_density, gr_float *DII_density, 
            gr_float *HDI_density, gr_float *e_density, gr_float *metal_density, gr_float *cooling_time); 
*/
        if (calculate_cooling_time(&cl->my_units, cl->a,
            1, one, zero, zero,
            &density, &energy,
            &x_velocity, &y_velocity, &z_velocity,
            &HI_density, &HII_density, &HM_density,
            &HeI_density, &HeII_density, &HeIII_density,
            &H2I_density, &H2II_density,
            &DI_density, &DII_density, &HDI_density,
            &e_density, &metal_density, 
            &cooling_time) == 0) {
            fprintf(stderr, "Grackle Error in calculate_cooling_time.\n");
            assert(0);
            }
        }

        return (ECode/cooling_time); /* Edot (erg/g/s)   undoes code in cool_multi_time_g.F */
    }



#endif /* NOCOOLING */
