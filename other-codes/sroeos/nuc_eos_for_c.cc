#include "nuc_eos_for_c.h"
#include "nuc_eos.hh"

extern "C" 
double nuc_eos_energy_shift_for_c() {
  return nuc_eos::energy_shift;
}

extern "C" 
double nuc_eos_rhomin_for_c() {
  return nuc_eos::eos_rhomin;
}

extern "C" 
double nuc_eos_yemax_for_c() {
  return nuc_eos::eos_yemax;
}

extern "C"
double nuc_eos_tempmin_for_c() { 
  return nuc_eos::eos_tempmin;
}

extern "C"
void nuc_eos_C_ReadTable_for_c(char* nuceos_table_name) {
  nuc_eos::nuc_eos_C_ReadTable( nuceos_table_name);
}

extern "C"
void nuc_eos_m_kt1_short_for_c(const int *n_in,
			   const double *rho, 
			   const double *temp,
			   const double *ye,
			   double *eps,
			   double *prs,
			   double *ent,
			   double *cs2,
			   double *dedt,
			   double *dpderho,
			   double *dpdrhoe,
			   double *munu,
			   int *keyerr,
			   int *anyerr) {
  nuc_eos::nuc_eos_m_kt1_short( n_in, rho, temp, ye, eps, prs, ent, cs2, 
				dedt, dpderho, dpdrhoe, munu, keyerr, anyerr);
}

extern "C"
void nuc_eos_m_kt0_short_for_c(const int *n_in,
			 const double *rho, 
			 double *temp,
			 const double *ye,
			 const double *eps,
			 double *prs,
			 double *ent,
			 double *cs2,
			 double *dedt,
			 double *dpderho,
			 double *dpdrhoe,
			 double *munu,
			 const double *prec,
			 int *keyerr,
			 int *anyerr) {

  nuc_eos::nuc_eos_m_kt0_short( n_in, rho, temp, ye, eps, prs, ent, cs2, 
				dedt, dpderho, dpdrhoe, munu,prec,keyerr, anyerr);
}

