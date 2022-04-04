#include "../../HydroParticleCool.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef MHD
#define NWAVE 6
#else // HYDRO
#define NWAVE 5
#endif
#ifndef MIN
#define MIN(a,b) ( ((a) < (b)) ? (a) : (b) )
#endif
#ifndef MAX
#define MAX(a,b) ( ((a) > (b)) ? (a) : (b) )
#endif
#define SIGN(a) ( ((a) < 0.) ? -1. : 1. )
#define SQR(x) ( (x)*(x) )
#define TINY_NUMBER 1.0e-30

typedef struct Cons1D_s{
  double d;			/*!< density */
  double Mx;			/*!< momentum density in X,Y,Z; where X is    */
  double My;                      /*!< direction longitudinal to 1D slice; which*/
  double Mz;                      /*!< can be in any dimension: 1,2,or 3        */
  double E;			/*!< total energy density */
  double Gamma;
#ifdef MHD
  double By;			/*!< cell centered magnetic fields in Y */
  double Bz;			/*!< cell centered magnetic fields in Z */
#endif /* MHD */
#if (NSCALARS > 0)
  double s[NSCALARS];             /*!< passively advected scalars */
#endif
}Cons1DS;

/*----------------------------------------------------------------------------*/
/*! \struct Prim1DS
 *  \brief Primitive variables in 1D (does not contain Bx).
 *  IMPORTANT!! The order of the elements in Prim1DS CANNOT be changed.
 */
typedef struct Prim1D_s{
  double d;			/*!< density */
  double Vx;			/*!< velocity in X-direction */
  double Vy;			/*!< velocity in Y-direction */
  double Vz;			/*!< velocity in Z-direction */
  double P;			/*!< pressure */
  double Gamma;
#ifdef MHD
  double By;			/*!< cell centered magnetic fields in Y-dir */
  double Bz;			/*!< cell centered magnetic fields in Z-dir */
#endif /* MHD */
#if (NSCALARS > 0)
  double r[NSCALARS];             /*!< density-normalized advected scalars */
#endif
}Prim1DS;

int fluxes(const Cons1DS Ul, const Cons1DS Ur,
            const Prim1DS Wl, const Prim1DS Wr,
            const double Bxi, Cons1DS *pF, double *MaxMach);

int hllFluxes(const Cons1DS Ul, const Cons1DS Ur,
            const Prim1DS Wl, const Prim1DS Wr,
            const double Bxi, Cons1DS *pF, double *MaxMach);

Prim1DS Cons1D_to_Prim1D(const Cons1DS *pU, const double *pBx);

void esys_roe_adb_hyd(const double v1, const double v2, const double v3,
  const double h, const double Gamma, double eigenvalues[],
  double right_eigenmatrix[][5], double left_eigenmatrix[][5]);

void esys_roe_adb_mhd(const double d, const double v1, const double v2, const double v3,
  const double h, const double b1, const double b2, const double b3, 
  const double x, const double y, const double Gamma,
  double eigenvalues[],
  double right_eigenmatrix[][7], double left_eigenmatrix[][7]);

#ifdef __cplusplus
}
#endif
