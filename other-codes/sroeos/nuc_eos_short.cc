#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nuc_eos.hh"
#include "helpers.hh"

namespace nuc_eos {
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
    // get prs
    {
      const int iv = 0;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(prs[i]),iv);
    }
    {
      const int iv = 1;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(eps[i]),iv);
    }
    // get entropy
    {
      const int iv = 2;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(ent[i]),iv);
    }
    // get munu
    {
      const int iv = 3;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(munu[i]),iv);
    }
    // get cs2
    {
      const int iv = 4;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(cs2[i]),iv);
    }
    // get dedT
    {
      const int iv = 5;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(dedt[i]),iv);
    }
    // get dpdrhoe
    {
      const int iv = 6;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(dpdrhoe[i]),iv);
    }
    // get dpderho
    {
      const int iv = 7;
      nuc_eos_C_linterp_one(idx,delx,dely,delz,&(dpderho[i]),iv);
    }
  }
  
  for(int i=0;i<n;i++) {
    prs[i] = exp(prs[i]);
    eps[i] = exp(eps[i]) - energy_shift;
  }

  return;
}


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
      *anyerr = 1;
    } else {
      temp[i] = exp(ltout);
      int idx[8];
      double delx,dely,delz;
      get_interp_spots(lr,ltout,ye[i],&delx,&dely,&delz,idx);
      // get prs
      {
	const int iv = 0;
	nuc_eos_C_linterp_one(idx,delx,dely,delz,&(prs[i]),iv);
      }
      // get entropy
      {
	const int iv = 2;
	nuc_eos_C_linterp_one(idx,delx,dely,delz,&(ent[i]),iv);
      }
      // get munu
      {
	const int iv = 3;
	nuc_eos_C_linterp_one(idx,delx,dely,delz,&(munu[i]),iv);
      }
      // get cs2
      {
	const int iv = 4;
	nuc_eos_C_linterp_one(idx,delx,dely,delz,&(cs2[i]),iv);
      }
      // get dedT
      {
	const int iv = 5;
	nuc_eos_C_linterp_one(idx,delx,dely,delz,&(dedt[i]),iv);
      }
      // get dpdrhoe
      {
	const int iv = 6;
	nuc_eos_C_linterp_one(idx,delx,dely,delz,&(dpdrhoe[i]),iv);
      }
      // get dpderho
      {
	const int iv = 7;
	nuc_eos_C_linterp_one(idx,delx,dely,delz,&(dpderho[i]),iv);
      }
    }
  }

  for(int i=0;i<n;i++) {
    prs[i] = exp(prs[i]);
  }

  return;
}


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

    for(int i=0;i<n;i++) {
      const double lr = log(rho[i]);
      const double lt = log(MIN(MAX(temp[i],eos_tempmin),eos_tempmax));
      double ltout;
      nuc_eos_findtemp_entropy(lr,lt,ye[i],ent[i],*prec,
			       (double *restrict)(&ltout),&keyerr[i]);

      if(keyerr[i] != 0) {
        *anyerr = 1;
      } else {
	temp[i] = exp(ltout);
	int idx[8];
	double delx,dely,delz;
	get_interp_spots(lr,ltout,ye[i],&delx,&dely,&delz,idx);
	// get prs
	{
	  const int iv = 0;
	  nuc_eos_C_linterp_one(idx,delx,dely,delz,&(prs[i]),iv);
	}
	// get eps
	{
	  const int iv = 1;
	  nuc_eos_C_linterp_one(idx,delx,dely,delz,&(eps[i]),iv);
	}
	// get munu
	{
	  const int iv = 3;
	  nuc_eos_C_linterp_one(idx,delx,dely,delz,&(munu[i]),iv);
	}
	// get cs2
	{
	  const int iv = 4;
	  nuc_eos_C_linterp_one(idx,delx,dely,delz,&(cs2[i]),iv);
	}
	// get dedT
	{
	  const int iv = 5;
	  nuc_eos_C_linterp_one(idx,delx,dely,delz,&(dedt[i]),iv);
	}
	// get dpdrhoe
	{
	  const int iv = 6;
	  nuc_eos_C_linterp_one(idx,delx,dely,delz,&(dpdrhoe[i]),iv);
	}
	// get dpderho
	{
	  const int iv = 7;
	  nuc_eos_C_linterp_one(idx,delx,dely,delz,&(dpderho[i]),iv);
	}
      }
    }

    for(int i=0;i<n;i++) {
      prs[i] = exp(prs[i]);
      eps[i] = exp(eps[i]) - energy_shift;
    }


  return;
}



} // namespace nuc_eos
