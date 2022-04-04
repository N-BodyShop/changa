#ifndef NUC_EOS_FOR_C_H 
#define NUC_EOS_FOR_C_H

#define LENGTHGF 6.77269222552442e-06
#define TIMEGF 2.03040204956746e05
#define RHOGF 1.61887093132742e-18
#define PRESSGF 1.80123683248503e-39
#define EPSGF 1.11265005605362e-21
#define INVRHOGF 6.17714470405638e17
#define INVEPSGF 8.98755178736818e20
#define INVPRESSGF 5.55174079257738e38
#ifdef __cplusplus
extern "C" {
#endif
  
  double nuc_eos_energy_shift_for_c();
  double nuc_eos_rhomin_for_c();
  double nuc_eos_yemax_for_c();
  double nuc_eos_tempmin_for_c();

  void nuc_eos_C_ReadTable_for_c(char* nuceos_table_name);

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
			   int *anyerr);

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
			 int *anyerr);
#ifdef __cplusplus
}
#endif

#endif

