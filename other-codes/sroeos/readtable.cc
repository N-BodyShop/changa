#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define H5_USE_16_API 1
#include <hdf5.h>
#include "nuc_eos.hh"

// Catch HDF5 errors
#define HDF5_ERROR(fn_call)                                              \
  do {                                                                   \
    int _error_code = fn_call;                                           \
    if (_error_code < 0) {	       				\
      fprintf(stderr, \
                  "HDF5 call '%s' returned error code %d",               \
                  #fn_call, _error_code);                                \
      abort(); }							\
  } while (0)

static int file_is_readable(const char* filename);
static int file_is_readable(const char* filename)
{
    FILE* fp = NULL;
    fp = fopen(filename, "r");
    if(fp != NULL)
    {
        fclose(fp);
        return 1;
    }
    return 0;
}


// define the variables
namespace nuc_eos {
  int nrho;
  int ntemp;
  int nye;

  double *alltables;
  double *epstable;
  double *logrho; 
  double *logtemp;
  double temp0, temp1;
  double dlintemp, dlintempi;
  double drholintempi;
  double dlintempyei;
  double drholintempyei;
  double *yes;
  double energy_shift;
  double dtemp, dtempi;
  double drho, drhoi;
  double dye, dyei;
  double drhotempi;
  double drhoyei;
  double dtempyei;
  double drhotempyei;
  
  double eos_rhomax, eos_rhomin;
  double eos_tempmin, eos_tempmax;
  double eos_yemin, eos_yemax;
  
  double c2p_tempmin;
  double c2p_tempmax;

  int ivs_short[8];
  int ivs_presseps[2];
}



namespace nuc_eos {

void nuc_eos_readtable_(char* nuceos_table_name)
{
  using namespace nuc_eos;
  nuc_eos_C_ReadTable( nuceos_table_name);

}
// Cactus calls this function. It reads in the table and calls a fortran
// function to setup values for the fortran eos module
void nuc_eos_C_ReadTable(char* nuceos_table_name)
{
  using namespace nuc_eos;

  fprintf(stdout,"*******************************\n");
  fprintf(stdout,"Reading nuc_eos table file:\n");
  fprintf(stdout,"%s\n",nuceos_table_name);
  fprintf(stdout,"*******************************\n");

  hid_t file;
  if (!file_is_readable(nuceos_table_name)) {
    fprintf(stderr,
               "Could not read nuceos_table_name %s \n",
               nuceos_table_name);
    abort();
  }
  HDF5_ERROR(file = H5Fopen(nuceos_table_name, H5F_ACC_RDONLY, H5P_DEFAULT));

// Use these two defines to easily read in a lot of variables in the same way
// The first reads in one variable of a given type completely
#define READ_EOS_HDF5(NAME,VAR,TYPE,MEM)                                      \
  do {                                                                        \
    hid_t dataset;                                                            \
    HDF5_ERROR(dataset = H5Dopen(file, NAME));                                \
    HDF5_ERROR(H5Dread(dataset, TYPE, MEM, H5S_ALL, H5P_DEFAULT, VAR));       \
    HDF5_ERROR(H5Dclose(dataset));                                            \
  } while (0)
// The second reads a given variable into a hyperslab of the alltables_temp array
#define READ_EOSTABLE_HDF5(NAME,OFF)                                     \
  do {                                                                   \
    hsize_t offset[2]     = {OFF,0};                                     \
    H5Sselect_hyperslab(mem3, H5S_SELECT_SET, offset, NULL, var3, NULL); \
    READ_EOS_HDF5(NAME,alltables_temp,H5T_NATIVE_DOUBLE,mem3);                \
  } while (0)

  // Read size of tables
  READ_EOS_HDF5("pointsrho",  &nrho,  H5T_NATIVE_INT, H5S_ALL);
  READ_EOS_HDF5("pointstemp", &ntemp, H5T_NATIVE_INT, H5S_ALL);
  READ_EOS_HDF5("pointsye",   &nye,   H5T_NATIVE_INT, H5S_ALL);


  // Allocate memory for tables
  double* alltables_temp;
  if (!(alltables_temp = (double*)malloc(nrho * ntemp * nye * NTABLES * sizeof(double)))) {
    fprintf(stderr, "Cannot allocate memory for EOS table\n");
    abort();
  }
  if (!(logrho = (double*)malloc(nrho * sizeof(double)))) {
    fprintf(stderr, "Cannot allocate memory for EOS table\n");
    abort();
  }
  if (!(logtemp = (double*)malloc(ntemp * sizeof(double)))) {
    fprintf(stderr, "Cannot allocate memory for EOS table\n");
    abort();
  }
  if (!(yes = (double*)malloc(nye * sizeof(double)))) {
    fprintf(stderr, "Cannot allocate memory for EOS table\n");
    abort();
  }

  // Prepare HDF5 to read hyperslabs into alltables_temp
  hsize_t table_dims[2] = {NTABLES, nrho * ntemp * nye};
  hsize_t var3[2]       = { 1, nrho * ntemp * nye};
  hid_t mem3 =  H5Screate_simple(2, table_dims, NULL);

  // Read alltables_temp
  READ_EOSTABLE_HDF5("logpress",  0);
  READ_EOSTABLE_HDF5("logenergy", 1);
  READ_EOSTABLE_HDF5("entropy",   2);
  READ_EOSTABLE_HDF5("munu",      3);
  READ_EOSTABLE_HDF5("cs2",       4);
  READ_EOSTABLE_HDF5("dedt",      5);
  READ_EOSTABLE_HDF5("dpdrhoe",   6);
  READ_EOSTABLE_HDF5("dpderho",   7);
  // chemical potentials
  READ_EOSTABLE_HDF5("muhat",     8);
  READ_EOSTABLE_HDF5("mu_e",      9);
  READ_EOSTABLE_HDF5("mu_p",     10);
  READ_EOSTABLE_HDF5("mu_n",     11);
  // compositions
  READ_EOSTABLE_HDF5("Xa",       12);
  READ_EOSTABLE_HDF5("Xh",       13);
  READ_EOSTABLE_HDF5("Xn",       14);
  READ_EOSTABLE_HDF5("Xp",       15);
  // average nucleus
  READ_EOSTABLE_HDF5("Abar",     16);
  READ_EOSTABLE_HDF5("Zbar",     17);
  // Gamma
  READ_EOSTABLE_HDF5("gamma",    18);

  // Read additional tables and variables
  READ_EOS_HDF5("logrho",       logrho,        H5T_NATIVE_DOUBLE, H5S_ALL);
  READ_EOS_HDF5("logtemp",      logtemp,       H5T_NATIVE_DOUBLE, H5S_ALL);
  READ_EOS_HDF5("ye",           yes,            H5T_NATIVE_DOUBLE, H5S_ALL);
  READ_EOS_HDF5("energy_shift", &energy_shift, H5T_NATIVE_DOUBLE, H5S_ALL);

  HDF5_ERROR(H5Sclose(mem3));

  // check if we already have the relativistic cs2
  int have_rel_cs2;
  {
    hid_t lapl_id = 0;
    htri_t check_have_rel_cs2 = H5Lexists(file,"have_rel_cs2",lapl_id);
    if (check_have_rel_cs2) {
      READ_EOS_HDF5("have_rel_cs2",&have_rel_cs2,H5T_NATIVE_INT,H5S_ALL);
    } else {
      have_rel_cs2 = 0;
    }
  }

  HDF5_ERROR(H5Fclose(file));

  // change ordering of alltables array so that
  // the table kind is the fastest changing index
  if (!(alltables = (double*)malloc(nrho * ntemp * nye * NTABLES 
				    * sizeof(double)))) {
    fprintf(stderr, "Cannot allocate memory for EOS table\n");
    abort();
  }
  for(int iv = 0;iv<NTABLES;iv++) 
    for(int k = 0; k<nye;k++) 
      for(int j = 0; j<ntemp; j++) 
	for(int i = 0; i<nrho; i++) {
	  int indold = i + nrho*(j + ntemp*(k + nye*iv));
	  int indnew = iv + NTABLES*(i + nrho*(j + ntemp*k));
	  alltables[indnew] = alltables_temp[indold];
	}

  // free memory of temporary array
  free(alltables_temp);

  // convert units, convert logs to natural log
  // The latter is great, because exp() is way faster than pow()
  // pressure
  energy_shift = energy_shift * EPSGF;
  for(int i=0;i<nrho;i++) {
    logrho[i] = log(pow(10.0,logrho[i]) * RHOGF);
  }

  for(int i=0;i<ntemp;i++) {
    logtemp[i] = log(pow(10.0,logtemp[i]));
  }

  // allocate epstable; a linear-scale eps table
  // that allows us to extrapolate to negative eps
  if (!(epstable = (double*)malloc(nrho * ntemp * nye  
				    * sizeof(double)))) {
    fprintf(stderr, "Cannot allocate memory for eps table\n");
    abort();
  }

  // convert units
  for(int i=0;i<nrho*ntemp*nye;i++) {

    { // pressure
      int idx = 0 + NTABLES*i;
      alltables[idx] = log(pow(10.0,alltables[idx])*PRESSGF);
    }

    { // eps
      int idx = 1 + NTABLES*i;
      alltables[idx] = log(pow(10.0,alltables[idx])*EPSGF);
      epstable[i] = exp(alltables[idx]);
    }

    { // cs2
      int idx = 4 + NTABLES*i;
      // unit conversion
      alltables[idx] *= LENGTHGF*LENGTHGF/TIMEGF/TIMEGF;
    }

    { // dedT
      int idx = 5 + NTABLES*i;
      alltables[idx] *= EPSGF;
    }

    { // dpdrhoe
      int idx = 6 + NTABLES*i;
      alltables[idx] *= PRESSGF/RHOGF;
    }

    { // dpderho
      int idx = 7 + NTABLES*i;
      alltables[idx] *= PRESSGF/EPSGF;
    }

  }

  temp0 = exp(logtemp[0]);
  temp1 = exp(logtemp[1]);

  // set up some vars
  dtemp  = (logtemp[ntemp-1] - logtemp[0]) / (1.0*(ntemp-1));
  dtempi = 1.0/dtemp;

  dlintemp = temp1-temp0;
  dlintempi = 1.0/dlintemp;

  drho  = (logrho[nrho-1] - logrho[0]) / (1.0*(nrho-1));
  drhoi = 1.0/drho;

  dye  = (yes[nye-1] - yes[0]) / (1.0*(nye-1));
  dyei = 1.0/dye;

  drhotempi   = drhoi * dtempi;
  drholintempi = drhoi * dlintempi;
  drhoyei     = drhoi * dyei;
  dtempyei    = dtempi * dyei;
  dlintempyei = dlintempi * dyei;
  drhotempyei = drhoi * dtempi * dyei;
  drholintempyei = drhoi * dlintempi * dyei;

  eos_rhomax = exp(logrho[nrho-1]);
  eos_rhomin = exp(logrho[0]);
  
  eos_tempmax = exp(logtemp[ntemp-1]);
  eos_tempmin = exp(logtemp[0]);

  eos_yemax = yes[nye-1];
  eos_yemin = yes[0];

  //set up index vectors for interpolations
  for(int i=0;i<8;i++) {
    ivs_short[i] = i;
  }

  // set up index vector for getting pressure and spec. int. energy
  ivs_presseps[0] = 0;
  ivs_presseps[1] = 1;

  
  if(!have_rel_cs2) {
    // in the case of legacy O'Connor & Ott 2010
    // stellarcollapse.org EOS tables, convert
    // to relativistic speed of sound by dividing
    // by relativistic enthalpy.

    fprintf(stdout,"Old style table -- converting cs2 to relativistic cs2!\n");

    int iv_cs2 = 4;
    int iv_p = 0;
    
    for(int k = 0; k<nye;k++)
      for(int j = 0; j<ntemp; j++)
        for(int i = 0; i<nrho; i++) {

	  int idx_cs2 = iv_cs2 + NTABLES*(i + nrho*(j + ntemp*k));
	  int idx_p = iv_p + NTABLES*(i + nrho*(j + ntemp*k)); 
	  int idx_eps = (i + nrho*(j + ntemp*k));

	  double xcs2 = alltables[idx_cs2];
	  double xp = exp(alltables[idx_p]);
	  double xeps = epstable[idx_eps] - energy_shift; 
	  double xrho = exp(logrho[i]);

	  double new_cs2 = CLIGHT*CLIGHT * xcs2 / 
	    (CLIGHT*CLIGHT + xeps + xp/xrho);

	  alltables[idx_cs2] = new_cs2;

	}
  }

}

} // namespace nuc_eos
