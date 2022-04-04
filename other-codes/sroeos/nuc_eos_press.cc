#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nuc_eos.hh"
#include "helpers.hh"

namespace nuc_eos {

extern "C"
void nuc_eos_m_kt1_press_eps(const int *restrict n_in,
			     const double *restrict rho, 
			     const double *restrict temp,
			     const double *restrict ye,
			     double *restrict eps,
			     double *restrict prs,
			     int *restrict keyerr,
			     int *restrict anyerr)
{

  using namespace nuc_eos;
  
  const int n = *n_in;

  *anyerr = 0;
  for(int i=0;i<n;i++) {
    // check if we are fine
    keyerr[i] = checkbounds(rho[i], temp[i], ye[i]);
    if(keyerr[i] != 0) {
      *anyerr = 1;
    }
  }

  // Abort if there is any error in checkbounds.
  // This should never happen and the program should abort with
  // a fatal error anyway. No point in doing any further EOS calculations.
  if(*anyerr) return;

  for(int i=0;i<n;i++) {
    int idx[8];
    double delx,dely,delz;
    const double xrho = log(rho[i]);
    const double xtemp = log(temp[i]);

    get_interp_spots(xrho,xtemp,ye[i],&delx,&dely,&delz,idx);

    {
      const int iv = 0;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(prs[i]),iv);
    }

    {
      const int iv = 1;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(eps[i]),iv);
    }
  }

  // now get rid of ln:
  for(int i=0;i<n;i++) {
    prs[i] = exp(prs[i]);
    eps[i] = exp(eps[i]) - energy_shift;
  }


  return;
}


void nuc_eos_m_kt0_press(const int *restrict n_in,
			 const double *restrict rho, 
			 double *restrict temp,
			 const double *restrict ye,
			 const double *restrict eps,
			 double *restrict prs,
			 const double *restrict prec,
			 int *restrict keyerr,
			 int *restrict anyerr)
{

  using namespace nuc_eos;

  const int n = *n_in;

  *anyerr = 0;

  for(int i=0;i<n;i++) {
    // check if we are fine
    // Note that this code now requires that the
    // temperature guess be within the table bounds
    keyerr[i] = checkbounds_kt0_noTcheck(rho[i], ye[i]);
    if(keyerr[i] != 0) {
      *anyerr = 1;
    }
  }

  // Abort if there is any error in checkbounds.
  // This should never happen and the program should abort with
  // a fatal error anyway. No point in doing any further EOS calculations.
  if(*anyerr) return;

  // first must find the temperature
  for(int i=0;i<n;i++) {
    const double lr = log(rho[i]);
    const double lt = log(MIN(MAX(temp[i],eos_tempmin),eos_tempmax));
    double ltout;
    const double epstot = eps[i]+energy_shift;
    if(epstot>0.0e0) {
      // this is the standard scenario; eps is larger than zero
      // and we can operate with logarithmic tables
      const double lxeps = log(epstot);
#if DEBUG
      fprintf(stderr,"%d %15.6E %15.6E %15.6E %15.6E\n",i,lr,lt,ye[i],lxeps);
      fprintf(stderr,"%d %15.6E %15.6E %15.6E %15.6E\n",i,
              exp(lr),exp(lt),ye[i],exp(lxeps));
#endif
      nuc_eos_findtemp(lr,lt,ye[i],lxeps,*prec,
                       (double *restrict)(&ltout),&keyerr[i]);
    } else {
      // will be overwritten further down, only marks error
      keyerr[i] = 667;
    } // epstot > 0.0

    if(keyerr[i] != 0) {
      // now try negative temperature treatment
      double eps0, eps1;
      int idx[8];
      double delx,dely,delz;

      get_interp_spots_linT_low_eps(lr,temp1,ye[i],&delx,&dely,&delz,idx);
      nuc_eos_C_linterp_one_linT_low_eps(idx,delx,dely,delz,&(eps1));

      get_interp_spots_linT_low_eps(lr,temp0,ye[i],&delx,&dely,&delz,idx);
      nuc_eos_C_linterp_one_linT_low_eps(idx,delx,dely,delz,&(eps0));

      temp[i] = (epstot-eps0) * (temp1-temp0)/(eps1-eps0) + temp0;

      // set error codes
      *anyerr = 1;
      keyerr[i] = 668;
      // get pressure
      {
        const int iv = 0;
        get_interp_spots_linT_low(lr,temp[i],ye[i],&delx,&dely,&delz,idx);
        nuc_eos_C_linterp_one_linT_low(idx,delx,dely,delz,&(prs[i]),iv);
      }
      
    } else {
      temp[i] = exp(ltout);
      int idx[8];
      double delx,dely,delz;
      get_interp_spots(lr,ltout,ye[i],&delx,&dely,&delz,idx);
      {
        const int iv = 0;
        nuc_eos_C_linterp_one(idx,delx,dely,delz,&(prs[i]),iv);
      }
    }
  } // loop i<n

  // now get rid of ln:
  for(int i=0;i<n;i++) {
    prs[i] = exp(prs[i]);
  }

  return;
}

} // namespace nuc_eos
