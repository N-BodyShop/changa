#ifndef NUC_EOS_HH
#define NUC_EOS_HH

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define NTABLES 19
#define LENGTHGF 6.77269222552442e-06
#define TIMEGF 2.03040204956746e05
#define RHOGF 1.61887093132742e-18
#define PRESSGF 1.80123683248503e-39
#define EPSGF 1.11265005605362e-21
#define INVRHOGF 6.17714470405638e17
#define INVEPSGF 8.98755178736818e20
#define INVPRESSGF 5.55174079257738e38
#define CLIGHT 1.0
#define restrict __restrict__

namespace nuc_eos {

  extern int nrho;
  extern int ntemp;
  extern int nye;

  extern double *alltables;
  extern double *epstable;
  extern double *logrho; 
  extern double *logtemp;
  extern double temp0, temp1;
  extern double dlintemp,dlintempi;
  extern double drholintempi;
  extern double dlintempyei;
  extern double drholintempyei;
  extern double *yes;
  extern double energy_shift;
  extern double dtemp, dtempi;
  extern double drho, drhoi;
  extern double dye, dyei;
  extern double drhotempi;
  extern double drhoyei;
  extern double dtempyei;
  extern double drhotempyei;

// min and max values

  extern double eos_rhomax, eos_rhomin;
  extern double eos_tempmin, eos_tempmax;
  extern double eos_yemin, eos_yemax;
  
  extern double c2p_tempmin;
  extern double c2p_tempmax;

// table key
// 0 logpress 
// 1 logenergy
// 2 entropy
// 3 munu
// 4 cs2
// 5 dedt
// 6 dpdrhoe
// 7 dpderho
// 8 muhat
// 9 mu_e
// 10 mu_p
// 11 mu_n
// 12 Xa
// 13 Xh
// 14 Xn
// 15 Xp
// 16 Abar
// 17 Zbar
// 18 Gamma

// some vectors for selecting variables for more
// efficient interpolation
  extern int ivs_short[8];
  extern int ivs_presseps[2];

extern "C"
  void nuc_eos_C_ReadTable(char* nuceos_table_name);

extern "C"
  void nuc_eos_readtable_(char* nuceos_table_name);

extern "C"
  void nuc_eos_m_kt1_press_eps(const int *restrict n_in,
			       const double *restrict rho, 
			       const double *restrict temp,
			       const double *restrict ye,
			       double *restrict eps,
			       double *restrict prs,
			       int *restrict keyerr,
			       int *restrict anyerr);
extern "C"
  void nuc_eos_m_kt0_press(const int *restrict n_in,
			   const double *restrict rho, 
			   double *restrict temp,
			   const double *restrict ye,
			   const double *restrict eps,
			   double *restrict prs,
			   const double *restrict prec,
			   int *restrict keyerr,
			   int *restrict anyerr);


  // call for derivatives needed in C2P; onlye keytemp=0 implemented
extern "C"
  void nuc_eos_m_kt0_dpdrhoe_dpderho(const int *restrict n_in,
				     const double *restrict rho, 
				     double *restrict temp,
				     const double *restrict ye,
				     const double *restrict eps,
				     double *restrict dpdrhoe,
				     double *restrict dpderho,
				     const double *restrict prec,
				     int *restrict keyerr,
				     int *restrict anyerr);
extern "C"
  void nuc_eos_m_kt1_press_eps_cs2(const int *restrict n_in,
				   const double *restrict rho, 
				   const double *restrict temp,
				   const double *restrict ye,
				   double *restrict eps,
				   double *restrict prs,
				   double *restrict cs2,
				   int *restrict keyerr,
				   int *restrict anyerr);
extern "C"
  void nuc_eos_m_kt0_press_cs2(const int *restrict n_in,
			       const double *restrict rho, 
			       double *restrict temp,
			       const double *restrict ye,
			       const double *restrict eps,
			       double *restrict prs,
			       double *restrict cs2,
			       const double *restrict prec,
			       int *restrict keyerr,
			       int *restrict anyerr);

extern "C"
  void nuc_eos_m_kt1_short(const int *restrict n_in,
			   const double *restrict rho, 
			   const double *restrict temp,
			   const double *restrict ye,
			   double *restrict eps,
			   double *restrict prs,
			   double *restrict ent,
			   double *restrict cs2,
			   double *restrict dedt,
			   double *restrict dpderho,
			   double *restrict dpdrhoe,
			   double *restrict munu,
			   int *restrict keyerr,
			   int *restrict anyerr);

extern "C"
void nuc_eos_m_kt0_short(const int *restrict n_in,
			 const double *restrict rho, 
			 double *restrict temp,
			 const double *restrict ye,
			 const double *restrict eps,
			 double *restrict prs,
			 double *restrict ent,
			 double *restrict cs2,
			 double *restrict dedt,
			 double *restrict dpderho,
			 double *restrict dpdrhoe,
			 double *restrict munu,
			 const double *restrict prec,
			 int *restrict keyerr,
			 int *restrict anyerr);

extern "C"
void nuc_eos_m_kt2_short(const int *restrict n_in,
			 const double *restrict rho, 
			 double *restrict temp,
			 const double *restrict ye,
			 double *restrict eps,
			 double *restrict prs,
			 const double *restrict ent,
			 double *restrict cs2,
			 double *restrict dedt,
			 double *restrict dpderho,
			 double *restrict dpdrhoe,
			 double *restrict munu,
			 const double *restrict prec,
			 int *restrict keyerr,
			 int *restrict anyerr);

extern "C"
  void nuc_eos_m_kt1_full(const int *restrict n_in,
			  const double *restrict rho, 
			  const double *restrict temp,
			  const double *restrict ye,
			  double *restrict eps,
			  double *restrict prs,
			  double *restrict ent,
			  double *restrict cs2,
			  double *restrict dedt,
			  double *restrict dpderho,
			  double *restrict dpdrhoe,
			  double *restrict xa,
			  double *restrict xh,
			  double *restrict xn,
			  double *restrict xp,
			  double *restrict abar,
			  double *restrict zbar,
			  double *restrict mue,
			  double *restrict mun,
			  double *restrict mup,
			  double *restrict muhat,
			  int *restrict keyerr,
			  int *restrict anyerr);

extern "C"
  void nuc_eos_m_kt0_full(const int *restrict n_in,
			  const double *restrict rho, 
			  double *restrict temp,
			  const double *restrict ye,
			  const double *restrict eps,
			  double *restrict prs,
			  double *restrict ent,
			  double *restrict cs2,
			  double *restrict dedt,
			  double *restrict dpderho,
			  double *restrict dpdrhoe,
			  double *restrict xa,
			  double *restrict xh,
			  double *restrict xn,
			  double *restrict xp,
			  double *restrict abar,
			  double *restrict zbar,
			  double *restrict mue,
			  double *restrict mun,
			  double *restrict mup,
			  double *restrict muhat,
			  const double *restrict prec,
			  int *restrict keyerr,
			  int *restrict anyerr);

  
  void testeos();

/////void nuc_eos_C_linterp_many(double x, double y, double z,
/////			    double* f, double* ft, 
/////			    int nx, int ny, int nz, int nvars,
/////			    double* xt,double*yt, double* zt);
/////
/////void nuc_eos_C_linterp_some(double x, double y, double z,
/////			    double* f, double* ft, 
/////			    int* ivs,
/////			    int nx, int ny, int nz, int nvars,
/////			    double* xt,double*yt, double* zt);
/////
/////void nuc_eos_C_linterp_for_temp(double x, double y, double z,
/////				double* f, double* ft, 
/////				int nx, int ny, int nz, 
/////				double* xt, double*yt, double* zt,
/////				double* linterp_for_temp);
/////
/////void nuc_eos_C_findtemp(double lr, double lt0, double ye, 
/////			double leps, double prec, double *lt,
/////			int *keyerr);
/////
/////void nuc_eos_C_testing();
/////
/////void dXdT(double* table, double* temp,
/////	  int nrho, int ntemp, int nye, int ntables,
/////	  int itable, double* diff);
/////
/////void dXdrho(double* table, double* rho,
/////	    int nrho, int ntemp, int nye, int ntables,
/////	    int itable, double* diff);
/////
/////
/////void dXdTforward(double* table, double* temp,
/////	  int nrho, int ntemp, int nye, int ntables,
/////	  int itable, double* diff);
/////
/////void dXdrhoforward(double* table, double* rho,
/////	    int nrho, int ntemp, int nye, int ntables,
/////	    int itable, double* diff);
/////

}

#endif // NUC_EOS_HH
