#include "copyright.h"
/*============================================================================*/
/*! \file hlld.c
 *  \brief Computes 1D fluxes using the HLLD Riemann solver.
 *
 * PURPOSE: Computes 1D fluxes using the HLLD Riemann solver, an extension of
 *   the HLLE solver to MHD.  Only works for MHD problems.  SEPARATE code
 *   blocks for adiabatic and isothermal equations of state.
 *
 * REFERENCES:
 * - T. Miyoshi & K. Kusano, "A multi-state HLL approximate Riemann solver
 *   for ideal MHD", JCP, 208, 315 (2005)
 * - A. Mignone, "A simple and accurate Riemann solver for isothermal MHD",
 *   JPC, 225, 1427 (2007)
 *
 * HISTORY: Adiabatic version written by Brian Biskeborn, May 8, 2006,
 *             COS sophmore project.
 *          Isothermal version written by Nicole Lemaster, May 1, 2008.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 * - fluxes() - all Riemann solvers in Athena must have this function name and
 *              use the same argument list as defined in rsolvers/prototypes.h*/
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "athena_defs.h"

#define SMALL_NUMBER 1e-30

#ifndef MHD
#error : The HLLD flux only works for mhd.
#endif /* MHD */

/*----------------------------------------------------------------------------*/
/*! \fn void fluxes(const Cons1DS Ul, const Cons1DS Ur,
 *           const Prim1DS Wl, const Prim1DS Wr, const double Bxi, Cons1DS *pFlux)
 *  \brief Compute 1D fluxes
 * Input Arguments:
 * - Bxi = B in direction of slice at cell interface
 * - Ul,Ur = L/R-states of CONSERVED variables at cell interface
 *
 * Output Arguments:
 * - Flux = fluxes of CONSERVED variables at cell interface
 */

int fluxes(const Cons1DS Ul, const Cons1DS Ur,
            const Prim1DS Wl, const Prim1DS Wr, const double Bxi, Cons1DS *pFlux,
          double *MaxMach)
{
  Cons1DS Ulst,Uldst,Urdst,Urst;       /* Conserved variable for all states */
  Prim1DS Wlst,Wrst;                   /* Primitive variables for all states */
  Cons1DS Fl,Fr;                       /* Fluxes for left & right states */
  double spd[5];                        /* signal speeds, left to right */
/*  double maxspd; */
  double sdl,sdr,sdml,sdmr;             /* S_i-u_i, S_i-S_M (i=L or R) */
  double pbl,pbr;                       /* Magnetic pressures */
  double cfl,cfr,cfmax;                 /* Cf (left & right), max(cfl,cfr) */
  double gpl,gpr,gpbl,gpbr;             /* gamma*P, gamma*P + B */
  double sqrtdl,sqrtdr;                 /* sqrt of the L* & R* densities */
  double invsumd;                       /* 1/(sqrtdl + sqrtdr) */
  double ptl,ptr,ptst;                  /* total pressures */
  double vbstl,vbstr;                   /* v_i* dot B_i* for i=L or R */
  double Bxsig;                         /* sign(Bx) = 1 for Bx>0, -1 for Bx<0 */
  double Bxsq;                          /* Bx^2 */
  double tmp;                      /* Temp variable for repeated calculations */
#if (NSCALARS > 0)
  int n;
#endif
  int isOk = 1;


/*--- Step 1. ------------------------------------------------------------------
 * Convert left- and right- states in conserved to primitive variables.
 */

/*
  pbl = Cons1D_to_Prim1D(&Ul,&Wl,&Bxi);
  pbr = Cons1D_to_Prim1D(&Ur,&Wr,&Bxi);
*/

/*--- Step 2. ------------------------------------------------------------------
 * Compute left & right wave speeds according to Miyoshi & Kusano, eqn. (67)
 */

  Bxsq = Bxi*Bxi;
  pbl = 0.5*(Bxsq + SQR(Wl.By) + SQR(Wl.Bz));
  pbr = 0.5*(Bxsq + SQR(Wr.By) + SQR(Wr.Bz));
  gpl  = Wl.Gamma * Wl.P;
  gpr  = Wr.Gamma * Wr.P;
  gpbl = gpl + 2.0*pbl;
  gpbr = gpr + 2.0*pbr;

  cfl = sqrt((gpbl + sqrt(SQR(gpbl)-4.0*gpl*Bxsq))/(2.0*Wl.d));
  cfr = sqrt((gpbr + sqrt(SQR(gpbr)-4.0*gpr*Bxsq))/(2.0*Wr.d));
  cfmax = MAX(cfl,cfr);

  if(Wl.Vx <= Wr.Vx) {
    spd[0] = Wl.Vx - cfmax;
    spd[4] = Wr.Vx + cfmax;
  }
  else {
    spd[0] = Wr.Vx - cfmax;
    spd[4] = Wl.Vx + cfmax;
  }

/*  maxspd = MAX(fabs(spd[0]),fabs(spd[4])); */

 *MaxMach = MAX(fabs(spd[0]/cfmax),fabs(spd[4]/cfmax));
/*--- Step 3. ------------------------------------------------------------------
 * Compute L/R fluxes
 */

  /* total pressure */
  ptl = Wl.P + pbl;
  ptr = Wr.P + pbr;

  Fl.d  = Ul.Mx;
  Fl.Mx = Ul.Mx*Wl.Vx + ptl - Bxsq;
  Fl.My = Ul.d*Wl.Vx*Wl.Vy - Bxi*Ul.By;
  Fl.Mz = Ul.d*Wl.Vx*Wl.Vz - Bxi*Ul.Bz;
  Fl.E  = Wl.Vx*(Ul.E + ptl - Bxsq) - Bxi*(Wl.Vy*Ul.By + Wl.Vz*Ul.Bz);
  Fl.By = Ul.By*Wl.Vx - Bxi*Wl.Vy;
  Fl.Bz = Ul.Bz*Wl.Vx - Bxi*Wl.Vz;

  Fr.d  = Ur.Mx;
  Fr.Mx = Ur.Mx*Wr.Vx + ptr - Bxsq;
  Fr.My = Ur.d*Wr.Vx*Wr.Vy - Bxi*Ur.By;
  Fr.Mz = Ur.d*Wr.Vx*Wr.Vz - Bxi*Ur.Bz;
  Fr.E  = Wr.Vx*(Ur.E + ptr - Bxsq) - Bxi*(Wr.Vy*Ur.By + Wr.Vz*Ur.Bz);
  Fr.By = Ur.By*Wr.Vx - Bxi*Wr.Vy;
  Fr.Bz = Ur.Bz*Wr.Vx - Bxi*Wr.Vz;

#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++) {
    Fl.s[n] = Fl.d*Wl.r[n];
    Fr.s[n] = Fr.d*Wr.r[n];
  }
#endif

/*--- Step 4. ------------------------------------------------------------------
 * Return upwind flux if flow is supersonic
 */

  if(spd[0] >= 0.0){
    *pFlux = Fl;
    return isOk;
  }

  if(spd[4] <= 0.0){
    *pFlux = Fr;
    return isOk;
  }

/*--- Step 5. ------------------------------------------------------------------
 * Compute middle and Alfven wave speeds
 */

  sdl = spd[0] - Wl.Vx;
  sdr = spd[4] - Wr.Vx;

  /* S_M: eqn (38) of Miyoshi & Kusano */
  spd[2] = (sdr*Wr.d*Wr.Vx - sdl*Wl.d*Wl.Vx - ptr + ptl) /
           (sdr*Wr.d-sdl*Wl.d);

  sdml   = spd[0] - spd[2];
  sdmr   = spd[4] - spd[2];
  /* eqn (43) of Miyoshi & Kusano */
  Ulst.d = Ul.d * sdl/sdml;
  Urst.d = Ur.d * sdr/sdmr;
  sqrtdl = sqrt(Ulst.d);
  sqrtdr = sqrt(Urst.d);

  /* eqn (51) of Miyoshi & Kusano */
  spd[1] = spd[2] - fabs(Bxi)/sqrtdl;
  spd[3] = spd[2] + fabs(Bxi)/sqrtdr;

/*--- Step 6. ------------------------------------------------------------------
 * Compute intermediate states
 */

  ptst = ptl + Ul.d*sdl*(sdl-sdml);
  if( ptst < 0. ) isOk = 0;

/* Ul* */
  /* eqn (39) of M&K */
  Ulst.Mx = Ulst.d * spd[2];
//   if((fabs(spd[2]/Wl.Vx-1.0)<SMALL_NUMBER) ||
//      (fabs(spd[2])/fabs(spd[0]) <= SMALL_NUMBER &&
//       fabs(Wl.Vx)/fabs(spd[0]) <= SMALL_NUMBER)) {
//     Ulst.My = Ulst.d * Wl.Vy;
//     Ulst.Mz = Ulst.d * Wl.Vz;
//
//     Ulst.By = Ul.By;
//     Ulst.Bz = Ul.Bz;
//   }
  if (fabs(Ul.d*sdl*sdml-Bxsq) < SMALL_NUMBER*ptst) {
    /* Degenerate case */
    Ulst.My = Ulst.d * Wl.Vy;
    Ulst.Mz = Ulst.d * Wl.Vz;

    Ulst.By = Ul.By;
    Ulst.Bz = Ul.Bz;
  }
  else {
    /* eqns (44) and (46) of M&K */
    tmp = Bxi*(sdl-sdml)/(Ul.d*sdl*sdml-Bxsq);
    Ulst.My = Ulst.d * (Wl.Vy - Ul.By*tmp);
    Ulst.Mz = Ulst.d * (Wl.Vz - Ul.Bz*tmp);
//     if(Ul.By == 0.0 && Ul.Bz == 0.0) {
//       Ulst.By = 0.0;
//       Ulst.Bz = 0.0;
//     }
//     else {
//       /* eqns (45) and (47) of M&K */
//       tmp = (Ul.d*SQR(sdl)-Bxsq)/(Ul.d*sdl*sdml - Bxsq);
//       Ulst.By = Ul.By * tmp;
//       Ulst.Bz = Ul.Bz * tmp;
//     }

    /* eqns (45) and (47) of M&K */
    tmp = (Ul.d*SQR(sdl)-Bxsq)/(Ul.d*sdl*sdml - Bxsq);
    Ulst.By = Ul.By * tmp;
    Ulst.Bz = Ul.Bz * tmp;
  }
  vbstl = (Ulst.Mx*Bxi+Ulst.My*Ulst.By+Ulst.Mz*Ulst.Bz)/Ulst.d;
  /* eqn (48) of M&K */
  Ulst.E = (sdl*Ul.E - ptl*Wl.Vx + ptst*spd[2] +
            Bxi*(Wl.Vx*Bxi+Wl.Vy*Ul.By+Wl.Vz*Ul.Bz - vbstl))/sdml;
  Wlst = Cons1D_to_Prim1D(&Ulst,&Bxi);


/* Ur* */
  /* eqn (39) of M&K */
  Urst.Mx = Urst.d * spd[2];
//   if((fabs(spd[2]/Wr.Vx-1.0)<SMALL_NUMBER) ||
//      (fabs(spd[2])/fabs(spd[4]) <= SMALL_NUMBER &&
//       fabs(Wr.Vx)/fabs(spd[4]) <= SMALL_NUMBER)) {
//     Urst.My = Urst.d * Wr.Vy;
//     Urst.Mz = Urst.d * Wr.Vz;
//
//     Urst.By = Ur.By;
//     Urst.Bz = Ur.Bz;
//   }
  if (fabs(Ur.d*sdr*sdmr-Bxsq) < SMALL_NUMBER*ptst) {
    /* Degenerate case */
    Urst.My = Urst.d * Wr.Vy;
    Urst.Mz = Urst.d * Wr.Vz;

    Urst.By = Ur.By;
    Urst.Bz = Ur.Bz;
  }
  else {
    /* eqns (44) and (46) of M&K */
    tmp = Bxi*(sdr-sdmr)/(Ur.d*sdr*sdmr-Bxsq);
    Urst.My = Urst.d * (Wr.Vy - Ur.By*tmp);
    Urst.Mz = Urst.d * (Wr.Vz - Ur.Bz*tmp);

//     if(Ur.By == 0.0 && Ur.Bz == 0.0) {
//       Urst.By = 0.0;
//       Urst.Bz = 0.0;
//     }
//     else {
//       /* eqns (45) and (47) of M&K */
//       tmp = (Ur.d*SQR(sdr)-Bxsq)/(Ur.d*sdr*sdmr - Bxsq);
//       Urst.By = Ur.By * tmp;
//       Urst.Bz = Ur.Bz * tmp;
//     }

    /* eqns (45) and (47) of M&K */
    tmp = (Ur.d*SQR(sdr)-Bxsq)/(Ur.d*sdr*sdmr - Bxsq);
    Urst.By = Ur.By * tmp;
    Urst.Bz = Ur.Bz * tmp;
  }
  vbstr = (Urst.Mx*Bxi+Urst.My*Urst.By+Urst.Mz*Urst.Bz)/Urst.d;
  /* eqn (48) of M&K */
  Urst.E = (sdr*Ur.E - ptr*Wr.Vx + ptst*spd[2] +
            Bxi*(Wr.Vx*Bxi+Wr.Vy*Ur.By+Wr.Vz*Ur.Bz - vbstr))/sdmr;
  Wrst = Cons1D_to_Prim1D(&Urst,&Bxi);


/* Ul** and Ur** - if Bx is zero, same as *-states */
//   if(Bxi == 0.0) {
  if(0.5*Bxsq < SMALL_NUMBER*ptst) {
    Uldst = Ulst;
    Urdst = Urst;
  }
  else {
    invsumd = 1.0/(sqrtdl + sqrtdr);
    if(Bxi > 0.0) Bxsig =  1.0;
    else          Bxsig = -1.0;

    Uldst.d = Ulst.d;
    Urdst.d = Urst.d;

    Uldst.Mx = Ulst.Mx;
    Urdst.Mx = Urst.Mx;

    /* eqn (59) of M&K */
    tmp = invsumd*(sqrtdl*Wlst.Vy + sqrtdr*Wrst.Vy + Bxsig*(Urst.By-Ulst.By));
    Uldst.My = Uldst.d * tmp;
    Urdst.My = Urdst.d * tmp;

    /* eqn (60) of M&K */
    tmp = invsumd*(sqrtdl*Wlst.Vz + sqrtdr*Wrst.Vz + Bxsig*(Urst.Bz-Ulst.Bz));
    Uldst.Mz = Uldst.d * tmp;
    Urdst.Mz = Urdst.d * tmp;

    /* eqn (61) of M&K */
    tmp = invsumd*(sqrtdl*Urst.By + sqrtdr*Ulst.By +
                   Bxsig*sqrtdl*sqrtdr*(Wrst.Vy-Wlst.Vy));
    Uldst.By = Urdst.By = tmp;

    /* eqn (62) of M&K */
    tmp = invsumd*(sqrtdl*Urst.Bz + sqrtdr*Ulst.Bz +
                   Bxsig*sqrtdl*sqrtdr*(Wrst.Vz-Wlst.Vz));
    Uldst.Bz = Urdst.Bz = tmp;

    /* eqn (63) of M&K */
    tmp = spd[2]*Bxi + (Uldst.My*Uldst.By + Uldst.Mz*Uldst.Bz)/Uldst.d;
    Uldst.E = Ulst.E - sqrtdl*Bxsig*(vbstl - tmp);
    Urdst.E = Urst.E + sqrtdr*Bxsig*(vbstr - tmp);
  }

/*--- Step 7. ------------------------------------------------------------------
 * Compute flux
 */

  if(spd[1] >= 0.0) {
/* return Fl* */
    pFlux->d  = Fl.d  + spd[0]*(Ulst.d  - Ul.d);
    pFlux->Mx = Fl.Mx + spd[0]*(Ulst.Mx - Ul.Mx);
    pFlux->My = Fl.My + spd[0]*(Ulst.My - Ul.My);
    pFlux->Mz = Fl.Mz + spd[0]*(Ulst.Mz - Ul.Mz);
    pFlux->E  = Fl.E  + spd[0]*(Ulst.E  - Ul.E);
    pFlux->By = Fl.By + spd[0]*(Ulst.By - Ul.By);
    pFlux->Bz = Fl.Bz + spd[0]*(Ulst.Bz - Ul.Bz);
  }
  else if(spd[2] >= 0.0) {
/* return Fl** */
    tmp = spd[1] - spd[0];
    pFlux->d  = Fl.d  - spd[0]*Ul.d  - tmp*Ulst.d  + spd[1]*Uldst.d;
    pFlux->Mx = Fl.Mx - spd[0]*Ul.Mx - tmp*Ulst.Mx + spd[1]*Uldst.Mx;
    pFlux->My = Fl.My - spd[0]*Ul.My - tmp*Ulst.My + spd[1]*Uldst.My;
    pFlux->Mz = Fl.Mz - spd[0]*Ul.Mz - tmp*Ulst.Mz + spd[1]*Uldst.Mz;
    pFlux->E  = Fl.E  - spd[0]*Ul.E  - tmp*Ulst.E  + spd[1]*Uldst.E;
    pFlux->By = Fl.By - spd[0]*Ul.By - tmp*Ulst.By + spd[1]*Uldst.By;
    pFlux->Bz = Fl.Bz - spd[0]*Ul.Bz - tmp*Ulst.Bz + spd[1]*Uldst.Bz;
  }
  else if(spd[3] > 0.0) {
/* return Fr** */
    tmp = spd[3] - spd[4];
    pFlux->d  = Fr.d  - spd[4]*Ur.d  - tmp*Urst.d  + spd[3]*Urdst.d;
    pFlux->Mx = Fr.Mx - spd[4]*Ur.Mx - tmp*Urst.Mx + spd[3]*Urdst.Mx;
    pFlux->My = Fr.My - spd[4]*Ur.My - tmp*Urst.My + spd[3]*Urdst.My;
    pFlux->Mz = Fr.Mz - spd[4]*Ur.Mz - tmp*Urst.Mz + spd[3]*Urdst.Mz;
    pFlux->E  = Fr.E  - spd[4]*Ur.E  - tmp*Urst.E  + spd[3]*Urdst.E;
    pFlux->By = Fr.By - spd[4]*Ur.By - tmp*Urst.By + spd[3]*Urdst.By;
    pFlux->Bz = Fr.Bz - spd[4]*Ur.Bz - tmp*Urst.Bz + spd[3]*Urdst.Bz;
  }
  else {
/* return Fr* */
    pFlux->d  = Fr.d  + spd[4]*(Urst.d  - Ur.d);
    pFlux->Mx = Fr.Mx + spd[4]*(Urst.Mx - Ur.Mx);
    pFlux->My = Fr.My + spd[4]*(Urst.My - Ur.My);
    pFlux->Mz = Fr.Mz + spd[4]*(Urst.Mz - Ur.Mz);
    pFlux->E  = Fr.E  + spd[4]*(Urst.E  - Ur.E);
    pFlux->By = Fr.By + spd[4]*(Urst.By - Ur.By);
    pFlux->Bz = Fr.Bz + spd[4]*(Urst.Bz - Ur.Bz);
  }

/* Fluxes of passively advected scalars, computed from density flux */
#if (NSCALARS > 0)
  if (pFlux->d >= 0.0) {
    for (n=0; n<NSCALARS; n++) pFlux->s[n] = pFlux->d*Wl.r[n];
  } else {
    for (n=0; n<NSCALARS; n++) pFlux->s[n] = pFlux->d*Wr.r[n];
  }
#endif

  return isOk;
}
