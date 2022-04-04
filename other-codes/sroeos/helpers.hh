#include "nuc_eos.hh"
#include <cstdlib>
#define DEBUG 0

namespace nuc_eos {

static inline __attribute__((always_inline))  
int checkbounds(const double xrho, 
		const double xtemp, 
		const double xye) {
  
  using namespace nuc_eos;

  // keyerr codes:
  // 101 -- Y_e too high
  // 102 -- Y_e too low
  // 103 -- temp too high (if keytemp = 1)
  // 104 -- temp too low (if keytemp = 1)
  // 105 -- rho too high
  // 106 -- rho too low

  if(xrho > eos_rhomax) {
    return 105;
  }
  if(xrho < eos_rhomin) {
    return 106;
  }
  if(xye > eos_yemax) {
    return 101;
  }
  if(xye < eos_yemin) {
    return 102;
  }
  if(xtemp > eos_tempmax) {
    return 103;
  }
  if(xtemp < eos_tempmin) {
    return 104;
  }
  return 0;
}

static inline __attribute__((always_inline))  
int checkbounds_kt0_noTcheck(const double xrho, 
		const double xye) {
  
  using namespace nuc_eos;

  // keyerr codes:
  // 101 -- Y_e too high
  // 102 -- Y_e too low
  // 105 -- rho too high
  // 106 -- rho too low

  if(xrho > eos_rhomax) {
    return 105;
  }
  if(xrho < eos_rhomin) {
    return 106;
  }
  if(xye > eos_yemax) {
    return 101;
  }
  if(xye < eos_yemin) {
    return 102;
  }
  return 0;
}


static inline __attribute__((always_inline))  
void get_interp_spots(const double x,
		      const double y,
		      const double z,
		      double* __restrict__ delx,
		      double* __restrict__ dely,
		      double* __restrict__ delz,
		      int* __restrict__ idx) 
{
  using namespace nuc_eos;

  int ix = 1 + (int)( (x - logrho[0] - 1.0e-10) * drhoi );
  int iy = 1 + (int)( (y - logtemp[0] - 1.0e-10) * dtempi );
  int iz = 1 + (int)( (z - yes[0] - 1.0e-10) * dyei );

  ix = MAX( 1, MIN( ix, nrho-1 ) );
  iy = MAX( 1, MIN( iy, ntemp-1 ) );
  iz = MAX( 1, MIN( iz, nye-1 ) );

  idx[0] = NTABLES*(ix + nrho*(iy + ntemp*iz));
  idx[1] = NTABLES*((ix-1) + nrho*(iy + ntemp*iz));
  idx[2] = NTABLES*(ix + nrho*((iy-1) + ntemp*iz));
  idx[3] = NTABLES*(ix + nrho*(iy + ntemp*(iz-1)));
  idx[4] = NTABLES*((ix-1) + nrho*((iy-1) + ntemp*iz));
  idx[5] = NTABLES*((ix-1) + nrho*(iy + ntemp*(iz-1)));
  idx[6] = NTABLES*(ix + nrho*((iy-1) + ntemp*(iz-1)));
  idx[7] = NTABLES*((ix-1) + nrho*((iy-1) + ntemp*(iz-1)));

   // set up aux vars for interpolation
  *delx = logrho[ix] - x;
  *dely = logtemp[iy] - y;
  *delz = yes[iz] - z;

  return;
}

static inline __attribute__((always_inline))  
void get_interp_spots_linT_low(const double x,
			       const double y,
			       const double z,
			       double* __restrict__ delx,
			       double* __restrict__ dely,
			       double* __restrict__ delz,
			       int* __restrict__ idx) 
{
  using namespace nuc_eos;

  int ix = 1 + (int)( (x - logrho[0] - 1.0e-10) * drhoi );
  int iy = 1;
  int iz = 1 + (int)( (z - yes[0] - 1.0e-10) * dyei );

  ix = MAX( 1, MIN( ix, nrho-1 ) );
  iz = MAX( 1, MIN( iz, nye-1 ) );

  idx[0] = NTABLES*(ix + nrho*(iy + ntemp*iz));
  idx[1] = NTABLES*((ix-1) + nrho*(iy + ntemp*iz));
  idx[2] = NTABLES*(ix + nrho*((iy-1) + ntemp*iz));
  idx[3] = NTABLES*(ix + nrho*(iy + ntemp*(iz-1)));
  idx[4] = NTABLES*((ix-1) + nrho*((iy-1) + ntemp*iz));
  idx[5] = NTABLES*((ix-1) + nrho*(iy + ntemp*(iz-1)));
  idx[6] = NTABLES*(ix + nrho*((iy-1) + ntemp*(iz-1)));
  idx[7] = NTABLES*((ix-1) + nrho*((iy-1) + ntemp*(iz-1)));

   // set up aux vars for interpolation
  *delx = logrho[ix] - x;
  *dely = temp1 - y;
  *delz = yes[iz] - z;

  return;
}

static inline __attribute__((always_inline))  
void get_interp_spots_linT_low_eps(const double x,
				   const double y,
				   const double z,
				   double* __restrict__ delx,
				   double* __restrict__ dely,
				   double* __restrict__ delz,
				   int* __restrict__ idx) 
{
  using namespace nuc_eos;

  int ix = 1 + (int)( (x - logrho[0] - 1.0e-10) * drhoi );
  int iy = 1;
  int iz = 1 + (int)( (z - yes[0] - 1.0e-10) * dyei );

  ix = MAX( 1, MIN( ix, nrho-1 ) );
  iz = MAX( 1, MIN( iz, nye-1 ) );

  idx[0] = (ix + nrho*(iy + ntemp*iz));
  idx[1] = ((ix-1) + nrho*(iy + ntemp*iz));
  idx[2] = (ix + nrho*((iy-1) + ntemp*iz));
  idx[3] = (ix + nrho*(iy + ntemp*(iz-1)));
  idx[4] = ((ix-1) + nrho*((iy-1) + ntemp*iz));
  idx[5] = ((ix-1) + nrho*(iy + ntemp*(iz-1)));
  idx[6] = (ix + nrho*((iy-1) + ntemp*(iz-1)));
  idx[7] = ((ix-1) + nrho*((iy-1) + ntemp*(iz-1)));

   // set up aux vars for interpolation
  *delx = logrho[ix] - x;
  *dely = temp1 - y;
  *delz = yes[iz] - z;

  return;
}


static inline __attribute__((always_inline))  
void nuc_eos_C_linterp_one(const int* __restrict__ idx, 
			   const double delx, 
			   const double dely,
			   const double delz,
			   double* __restrict__ f,
			   const int iv)
{
  using namespace nuc_eos;

  // helper variables
  double fh[8], a[8];

  fh[0] = alltables[iv+idx[0]];
  fh[1] = alltables[iv+idx[1]];
  fh[2] = alltables[iv+idx[2]];
  fh[3] = alltables[iv+idx[3]];
  fh[4] = alltables[iv+idx[4]];
  fh[5] = alltables[iv+idx[5]];
  fh[6] = alltables[iv+idx[6]];
  fh[7] = alltables[iv+idx[7]];

  // set up coeffs of interpolation polynomical and
  // evaluate function values
  a[0] = fh[0];
  a[1] = drhoi *   ( fh[1] - fh[0] );
  a[2] = dtempi *   ( fh[2] - fh[0] );
  a[3] = dyei *   ( fh[3] - fh[0] );
  a[4] = drhotempi *  ( fh[4] - fh[1] - fh[2] + fh[0] );
  a[5] = drhoyei *  ( fh[5] - fh[1] - fh[3] + fh[0] );
  a[6] = dtempyei *  ( fh[6] - fh[2] - fh[3] + fh[0] );
  a[7] = drhotempyei * ( fh[7] - fh[0] + fh[1] + fh[2] + 
			 fh[3] - fh[4] - fh[5] - fh[6] );

  *f = a[0] + a[1] * delx
    + a[2] * dely
    + a[3] * delz
    + a[4] * delx * dely
    + a[5] * delx * delz
    + a[6] * dely * delz
    + a[7] * delx * dely * delz;

  return;
}

static inline __attribute__((always_inline))  
void nuc_eos_C_linterp_one_linT_low(const int* __restrict__ idx, 
				    const double delx, 
				    const double dely,
				    const double delz,
				    double* __restrict__ f,
				    const int iv)
{
  using namespace nuc_eos;

  // helper variables
  double fh[8], a[8];

  fh[0] = alltables[iv+idx[0]];
  fh[1] = alltables[iv+idx[1]];
  fh[2] = alltables[iv+idx[2]];
  fh[3] = alltables[iv+idx[3]];
  fh[4] = alltables[iv+idx[4]];
  fh[5] = alltables[iv+idx[5]];
  fh[6] = alltables[iv+idx[6]];
  fh[7] = alltables[iv+idx[7]];

  // set up coeffs of interpolation polynomical and
  // evaluate function values
  a[0] = fh[0];
  a[1] = drhoi *   ( fh[1] - fh[0] );
  a[2] = dlintempi *   ( fh[2] - fh[0] );
  a[3] = dyei *   ( fh[3] - fh[0] );
  a[4] = drholintempi *  ( fh[4] - fh[1] - fh[2] + fh[0] );
  a[5] = drhoyei *  ( fh[5] - fh[1] - fh[3] + fh[0] );
  a[6] = dlintempyei *  ( fh[6] - fh[2] - fh[3] + fh[0] );
  a[7] = drholintempyei * ( fh[7] - fh[0] + fh[1] + fh[2] + 
			 fh[3] - fh[4] - fh[5] - fh[6] );

  *f = a[0] + a[1] * delx
    + a[2] * dely
    + a[3] * delz
    + a[4] * delx * dely
    + a[5] * delx * delz
    + a[6] * dely * delz
    + a[7] * delx * dely * delz;

  return;
}


static inline __attribute__((always_inline))  
void nuc_eos_C_linterp_one_linT_low_eps(const int* __restrict__ idx, 
					const double delx, 
					const double dely,
					const double delz,
					double* __restrict__ f)

{
  using namespace nuc_eos;

  // helper variables
  double fh[8], a[8];

  fh[0] = epstable[idx[0]];
  fh[1] = epstable[idx[1]];
  fh[2] = epstable[idx[2]];
  fh[3] = epstable[idx[3]];
  fh[4] = epstable[idx[4]];
  fh[5] = epstable[idx[5]];
  fh[6] = epstable[idx[6]];
  fh[7] = epstable[idx[7]];

  // set up coeffs of interpolation polynomical and
  // evaluate function values
  a[0] = fh[0];
  a[1] = drhoi *   ( fh[1] - fh[0] );
  a[2] = dlintempi *   ( fh[2] - fh[0] );
  a[3] = dyei *   ( fh[3] - fh[0] );
  a[4] = drholintempi *  ( fh[4] - fh[1] - fh[2] + fh[0] );
  a[5] = drhoyei *  ( fh[5] - fh[1] - fh[3] + fh[0] );
  a[6] = dlintempyei *  ( fh[6] - fh[2] - fh[3] + fh[0] );
  a[7] = drholintempyei * ( fh[7] - fh[0] + fh[1] + fh[2] + 
			 fh[3] - fh[4] - fh[5] - fh[6] );

  *f = a[0] + a[1] * delx
    + a[2] * dely
    + a[3] * delz
    + a[4] * delx * dely
    + a[5] * delx * delz
    + a[6] * dely * delz
    + a[7] * delx * dely * delz;

  return;
}




static inline __attribute__((always_inline))
double linterp2D(const double *restrict xs, 
		 const double *restrict ys, 
		 const double *restrict fs, 
		 const double x, 
		 const double y)
{

  //  2     3 
  //
  //  0     1
  //
  // first interpolate in x between 0 and 1, 2 and 3
  // then interpolate in y
  // assume rectangular grid
  
  double t1 = (fs[1]-fs[0])/(xs[1]-xs[0]) * (x - xs[0]) + fs[0];
  double t2 = (fs[3]-fs[2])/(xs[1]-xs[0]) * (x - xs[0]) + fs[2];

  return (t2 - t1)/(ys[1]-ys[0]) * (y-ys[0]) + t1;
}

static inline __attribute__((always_inline))
void bisection(const double lr, 
	       const double lt0,
	       const double ye,
	       const double leps0,
	       const double prec,
	       double *restrict ltout,
	       const int iv,
	       int *restrict keyerrt) {
  // iv is the index of the variable we do the bisection on

  using namespace nuc_eos;

  int bcount = 0; 
  int maxbcount = 80;
  int itmax = 50;

  const double dlt0p = log(1.1);
  const double dlt0m = log(0.9);
  const double dltp = log(1.2);
  const double dltm = log(0.8);

  // temporary local vars
  double lt, lt1, lt2;
  double ltmin = logtemp[0];
  double ltmax = logtemp[ntemp-1];
  double f1,f2,fmid,dlt,ltmid;
  double f1a = 0.0;
  double f2a = 0.0;
  double delx,dely,delz;
  int idx[8];
  

  // prepare
  lt = lt0;
  lt1 = MIN(lt0 + dlt0p,ltmax);
  lt2 = MAX(lt0 + dlt0m,ltmin);

  get_interp_spots(lr,lt1,ye,&delx,&dely,&delz,idx);
  nuc_eos_C_linterp_one(idx,delx,dely,delz,&f1a,iv);  

  get_interp_spots(lr,lt2,ye,&delx,&dely,&delz,idx);
  nuc_eos_C_linterp_one(idx,delx,dely,delz,&f2a,iv);  

  f1=f1a-leps0;
  f2=f2a-leps0;

  // iterate until we bracket the right eps, but enforce
  // dE/dt > 0, so eps(lt1) > eps(lt2)
  while(f1*f2 >= 0.0) {
    lt1 = MIN(lt1 + dltp,ltmax);
    lt2 = MAX(lt2 + dltm,ltmin);
    get_interp_spots(lr,lt1,ye,&delx,&dely,&delz,idx);
    nuc_eos_C_linterp_one(idx,delx,dely,delz,&f1a,iv);  
    
    get_interp_spots(lr,lt2,ye,&delx,&dely,&delz,idx);
    nuc_eos_C_linterp_one(idx,delx,dely,delz,&f2a,iv);  

    f1=f1a-leps0;
    f2=f2a-leps0;

#if DEBUG
    fprintf(stderr,"bisection bracketing it %d, f1: %15.6E, f2: %15.6E, lt1: %15.6E, lt2: %15.6E, f1a: %18.11E, f2a: %18.11E leps0: %18.11E\n",
	    bcount,f1,f2,lt1,lt2,f1a,f2a,leps0);
#endif

    bcount++;
    if(bcount >= maxbcount) {
      *keyerrt = 667;
      return;
    }
  } // while

  if(f1 < 0.0) {
    lt = lt1;
    dlt = lt2 - lt1;
  } else {
    lt = lt2;
    dlt = lt1 - lt2;
  }

#if DEBUG
    fprintf(stderr,"bisection step 2 it -1, fmid: %15.6E ltmid: %15.6E dlt: %15.6E\n",
	    f2,lt,dlt);
    fprintf(stderr,"ltmax: %15.6E\n",ltmax);
#endif

  int it;
  for(it=0;it<itmax;it++) {
    dlt = dlt * 0.5;
    ltmid = lt + dlt;
    get_interp_spots(lr,ltmid,ye,&delx,&dely,&delz,idx);
    nuc_eos_C_linterp_one(idx,delx,dely,delz,&f2a,iv); 
    
    fmid=f2a-leps0;
    if(fmid <= 0.0) lt=ltmid;
#if DEBUG
    fprintf(stderr,"bisection step 2 it %d, fmid: %15.6E f2a: %15.6E lt: %15.6E ltmid: %15.6E dlt: %15.6E\n",
	    it,fmid,f2a,lt,ltmid,dlt);
#endif

    if(fabs(1.0-f2a/leps0) <= prec) {
      *ltout = ltmid;
      return;
    }
  } // for it = 0

  if(it >= itmax-1) {
    *keyerrt = 667;
    return;
  }
} // bisection



static inline __attribute__((always_inline))
void nuc_eos_findtemp(const double lr, 
		      const double lt0,
		      const double ye,
		      const double lepsin,
		      const double prec,
		      double *restrict ltout,
		      int *keyerrt) {

  using namespace nuc_eos;

  // local variables
  const int itmax = 200; // use at most 10 iterations, then go to bisection
  double dlepsdlt; // derivative dlogeps/dlogT
  double ldt;
  double leps,leps0; // temp vars for eps
  double ltn, lt; // temp vars for temperature
  double oerr;
  const double ltmax = logtemp[ntemp-1]; // max temp
  const double ltmin = logtemp[0]; // min temp
  const int iv = 1;
  int it = 0;

  // setting up some vars
  *keyerrt = 0;
  leps0 = lepsin;
  lt = lt0;

  // step 1: do we already have the right temperature
  int idx[8];
  double delx,dely,delz;
  get_interp_spots(lr,lt,ye,&delx,&dely,&delz,idx);
  nuc_eos_C_linterp_one(idx,delx,dely,delz,&leps,iv);
#if DEBUG
  fprintf(stderr,"it: %d t: %15.6E leps: %15.6E eps0: %15.6E del: %15.6E\n",
	  it,lt, leps,leps0,fabs(leps-leps0)/(fabs(leps0)));
#endif
  if(fabs(leps-leps0) < prec*fabs(leps0)) {
    *ltout = lt0;
    return;
  }

  oerr = 1.0e90;
  double fac = 1.0;
  const int irho = MIN(MAX(1 + (int)(( lr - logrho[0] - 1.0e-12) * drhoi),1),nrho-1); 
  const int iye = MIN(MAX(1 + (int)(( ye - yes[0] - 1.0e-12) * dyei),1),nye-1); 

  while(it < itmax) {
    it++;

    // step 2: check if the two bounding values of the temperature
    //         give eps values that enclose the new eps.
    const int itemp = MIN(MAX(1 + (int)(( lt - logtemp[0] - 1.0e-12) * dtempi),1),ntemp-1); 

    double epst1, epst2;
    // lower temperature
    {
      // get data at 4 points
      double fs[4];
      // point 0
      int ifs = 1 + NTABLES*(irho-1 + nrho*((itemp-1) + ntemp*(iye-1)));
      fs[0] = alltables[ifs];
      // point 1 
      ifs = 1 + NTABLES*(irho + nrho*((itemp-1) + ntemp*(iye-1)));
      fs[1] = alltables[ifs];
      // point 2 
      ifs = 1 + NTABLES*(irho-1 + nrho*((itemp-1) + ntemp*(iye)));
      fs[2] = alltables[ifs];
      // point 3
      ifs = 1 + NTABLES*(irho + nrho*((itemp-1) + ntemp*(iye)));
      fs[3] = alltables[ifs];
      
      epst1 = linterp2D(&logrho[irho-1],&yes[iye-1], fs, lr, ye);
    }
    // upper temperature
    {
      // get data at 4 points
      double fs[4];
      // point 0
      int ifs = 1 + NTABLES*(irho-1 + nrho*((itemp) + ntemp*(iye-1)));
      fs[0] = alltables[ifs];
      // point 1 
      ifs = 1 + NTABLES*(irho + nrho*((itemp) + ntemp*(iye-1)));
      fs[1] = alltables[ifs];
      // point 2 
      ifs = 1 + NTABLES*(irho-1 + nrho*((itemp) + ntemp*(iye)));
      fs[2] = alltables[ifs];
      // point 3
      ifs = 1 + NTABLES*(irho + nrho*((itemp) + ntemp*(iye)));
      fs[3] = alltables[ifs];
      
      epst2 = linterp2D(&logrho[irho-1],&yes[iye-1], fs, lr, ye);
    }

    // Check if we are already bracketing the input internal
    // energy. If so, interpolate for new T.
    if( (leps0 >= epst1 && leps0 <= epst2) 
	|| (leps0 <= epst1 && leps0 >=epst2)) {
      
      *ltout = (logtemp[itemp]-logtemp[itemp-1]) / (epst2 - epst1) * 
	(leps0 - epst1) + logtemp[itemp-1];
     
      return;
    }

    // well, then do a Newton-Raphson step
    // first, guess the derivative
    dlepsdlt = (epst2-epst1)/(logtemp[itemp]-logtemp[itemp-1]);
    ldt = -(leps - leps0) / dlepsdlt * fac;

    ltn = MIN(MAX(lt + ldt,ltmin),ltmax);
    lt = ltn;

    get_interp_spots(lr,lt,ye,&delx,&dely,&delz,idx);
    nuc_eos_C_linterp_one(idx,delx,dely,delz,&leps,iv);

#if DEBUG
    fprintf(stderr,"%d %d %d\n",irho,itemp,iye);
    fprintf(stderr,"it: %d t: %15.6E leps: %15.6E eps0: %15.6E del: %15.6E\n",
	    it,lt, leps,leps0,fabs(leps-leps0)/(fabs(leps0)));
#endif
    // drive the thing into the right direction
    double err = fabs(leps-leps0);
    if(oerr < err) fac *= 0.9;
    oerr = err;

    if(err < prec*fabs(leps0)) {
      *ltout = lt;
      return;
    }



  } // while(it < itmax)

  if(it >= itmax - 1) {
    // try bisection
    // bisection(lr, lt0, ye, leps0, ltout, 1, prec, keyerrt);
#if DEBUG
    fprintf(stderr, "Failed to converge. This is bad. Trying bisection!\n");
#endif
    bisection(lr,lt0,ye,leps0,prec,ltout,1,keyerrt);
#if DEBUG
    if(*keyerrt==667) {
      fprintf(stderr,"This is worse. Bisection failed!\n");
      abort();
    }      
#endif
  }


  return;
}

static inline __attribute__((always_inline))
void nuc_eos_findtemp_entropy(const double lr, 
		      const double lt0,
		      const double ye,
		      const double entin,
		      const double prec,
		      double *restrict ltout,
		      int *keyerrt) {

  using namespace nuc_eos;

  // local variables
  const int itmax = 200; // use at most 10 iterations, then go to bisection
  double dentdlt; // derivative dentropy/dlogT
  double ldt;
  double ent,ent0; // temp vars for eps
  double ltn, lt; // temp vars for temperature
  double oerr;
  const double ltmax = logtemp[ntemp-1]; // max temp
  const double ltmin = logtemp[0]; // min temp
  const int iv = 2;
  int it = 0;

  // setting up some vars
  *keyerrt = 0;
  ent0 = entin;
  lt = lt0;

  // step 1: do we already have the right temperature
  int idx[8];
  double delx,dely,delz;
  get_interp_spots(lr,lt,ye,&delx,&dely,&delz,idx);
  nuc_eos_C_linterp_one(idx,delx,dely,delz,&ent,iv);
#if DEBUG
  fprintf(stderr,"it: %d t: %15.6E leps: %15.6E eps0: %15.6E del: %15.6E\n",
	  it,lt,ent,ent0,fabs(ent-ent0)/(fabs(ent0)));
#endif
  if(fabs(ent-ent0) < prec*fabs(ent0)) {
    *ltout = lt0;
    return;
  }

  oerr = 1.0e90;
  double fac = 1.0;
  const int irho = MIN(MAX(1 + (int)(( lr - logrho[0] - 1.0e-12) * drhoi),1),nrho-1); 
  const int iye = MIN(MAX(1 + (int)(( ye - yes[0] - 1.0e-12) * dyei),1),nye-1); 

  while(it < itmax) {
    it++;

    // step 2: check if the two bounding values of the temperature
    //         give eps values that enclose the new eps.
    const int itemp = MIN(MAX(1 + (int)(( lt - logtemp[0] - 1.0e-12) * dtempi),1),ntemp-1); 

    double ent1, ent2;
    // lower temperature
    {
      // get data at 4 points
      double fs[4];
      // point 0
      int ifs = 2 + NTABLES*(irho-1 + nrho*((itemp-1) + ntemp*(iye-1)));
      fs[0] = alltables[ifs];
      // point 1 
      ifs = 2 + NTABLES*(irho + nrho*((itemp-1) + ntemp*(iye-1)));
      fs[1] = alltables[ifs];
      // point 2 
      ifs = 2 + NTABLES*(irho-1 + nrho*((itemp-1) + ntemp*(iye)));
      fs[2] = alltables[ifs];
      // point 3
      ifs = 2 + NTABLES*(irho + nrho*((itemp-1) + ntemp*(iye)));
      fs[3] = alltables[ifs];
      
      ent1 = linterp2D(&logrho[irho-1],&yes[iye-1], fs, lr, ye);
    }
    // upper temperature
    {
      // get data at 4 points
      double fs[4];
      // point 0
      int ifs = 2 + NTABLES*(irho-1 + nrho*((itemp) + ntemp*(iye-1)));
      fs[0] = alltables[ifs];
      // point 1 
      ifs = 2 + NTABLES*(irho + nrho*((itemp) + ntemp*(iye-1)));
      fs[1] = alltables[ifs];
      // point 2 
      ifs = 2 + NTABLES*(irho-1 + nrho*((itemp) + ntemp*(iye)));
      fs[2] = alltables[ifs];
      // point 3
      ifs = 2 + NTABLES*(irho + nrho*((itemp) + ntemp*(iye)));
      fs[3] = alltables[ifs];
      
      ent2 = linterp2D(&logrho[irho-1],&yes[iye-1], fs, lr, ye);
    }

    // Check if we are already bracketing the input internal
    // energy. If so, interpolate for new T.
    if( (ent0 >= ent1 && ent0 <= ent2) 
	|| (ent0 <= ent1 && ent0 >= ent2)) {
      
      *ltout = (logtemp[itemp]-logtemp[itemp-1]) / (ent2 - ent1) * 
	(ent0 - ent1) + logtemp[itemp-1];
     
      return;
    }

    // well, then do a Newton-Raphson step
    // first, guess the derivative
    dentdlt = (ent2-ent1)/(logtemp[itemp]-logtemp[itemp-1]);
    ldt = -(ent - ent0) / dentdlt * fac;

    ltn = MIN(MAX(lt + ldt,ltmin),ltmax);
    lt = ltn;

    get_interp_spots(lr,lt,ye,&delx,&dely,&delz,idx);
    nuc_eos_C_linterp_one(idx,delx,dely,delz,&ent,iv);

#if DEBUG
    fprintf(stderr,"%d %d %d\n",irho,itemp,iye);
    fprintf(stderr,"it: %d t: %15.6E ent: %15.6E ent0: %15.6E del: %15.6E\n",
	    it,lt, ent,ent0,fabs(ent-ent0)/(fabs(ent0)));
#endif
    // drive the thing into the right direction
    double err = fabs(ent-ent0);
    if(oerr < err) fac *= 0.9;
    oerr = err;

    if(err < prec*fabs(ent0)) {
      *ltout = lt;
      return;
    }



  } // while(it < itmax)

  if(it >= itmax - 1) {
    // try bisection
#if DEBUG
    fprintf(stderr, "Failed to converge. This is bad. Trying bisection!\n");
#endif
    bisection(lr,lt0,ye,ent0,prec,ltout,2,keyerrt);
#if DEBUG
    if(*keyerrt==667) {
      fprintf(stderr,"This is worse. Bisection failed!\n");
      abort();
    }      
#endif
  }


  return;
}

} // namespace nuc_eos
