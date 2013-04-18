#include <math.h>
#include "ParallelGravity.h"
#include "feedback.h"
#include "starlifetime.h"

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))

/* MAXIMUM MASS: 109 M_sol, lifetime= 3.26 Myr */
void Padova::CoefInit (double dMetals)
{
  double logZ;
  double Z = dMetals;
  
  Z = min(Z, zmax);
  Z = max(Z, zmin);
  logZ = log10(Z);
  a0 = a00 + a01*logZ + a02*logZ*logZ;
  a1 = a10 + a11*logZ + a12*logZ*logZ;
  a2 = a20 + a21*logZ + a22*logZ*logZ;
}

double Padova::Lifetime(double dStarMass, double dMetals) 
{
  /* finds stellar lifetime in yr corresponding to stellar 
     mass in solar masses */
  double logStarMass, logLtime, Ltime;
  
  logStarMass = log10(dStarMass);
  CoefInit (dMetals);
  logLtime = a0 + a1 *logStarMass + a2*logStarMass*logStarMass;
  Ltime = pow(10.0, logLtime);
  
  return Ltime;
}

double Padova::StarMass(double dStarLtime, double dMetals) 
{
  /* finds stellar mass in solar masses corresponding to stellar 
     lifetime dStarTime in yr */
  double logStarMass, StarMass;
  double a, b, c;
  
  if(dStarLtime <= 0.0)		/* Time can be zero */
    return 1000;
  
  CoefInit (dMetals);
  c = a0;
  c -= log10(dStarLtime);
  b = a1;
  a = a2;
  if(b*b - 4*a*c < 0.0)		/* time is too small for fitting
				   formula */
    return 1000;
  
  logStarMass = (-b - sqrt(b*b - 4*a*c))/(2*a);
  StarMass = pow(10., logStarMass);
  return StarMass;
}


#if 0
#include "supernova.h"

int
main(int argc, char **argv)
{
    int nsamp;
    int i;
    double dlgm;
    double lgm;
    double Ntot;
    double Mtot;

    PDVAPARAM ppdva;
    SN sn;
    double dStarMass;
    double mass[30], time[30];
    PARTICLE *p;
    double dMetals = 0.0;
    
    assert(argc == 2);
    
    nsamp = atoi(argv[1]);
    dlgm = (2.0 + 1.0)/nsamp;
    snInitialize (&sn);
	//    Ntot = dMSCumNumber(&MSparam, 0.0);
	//    Mtot = dMSCumMass(&MSparam, 0.0);
    
	//    dSTLtimeMStar (ppdva, sn, dStarMass, p);
	//    dSTMStarLtime (ppdva, sn, dStarLtime, p);
    
    PadovaInitialize(&ppdva);
	//    mass = 8;
	//    printf ("first %g %g\n", mass, dSTLtimeMStar (ppdva, mass, p));

    printf("mass(M_sol) time[yr]\n");        
    for(i = 0; i < nsamp; i++) {
	
		//	lgm = -1 + i*dlgm;
		lgm = i*dlgm;

		mass[i] = pow(10.0, lgm);
        time[i] = dSTLtimeMStar (ppdva, mass[i], dMetals);
		printf("%g %g\n", mass[i], time[i]);        

		}

    printf("time[yr] mass(M_sol)\n");        
    for(i = 0; i < nsamp; i++) {

        mass[i] = dSTMStarLtime (ppdva, time[i], dMetals);
		printf("%g %g\n", mass[i], time[i]);        
		}
    
  return 0;
}
#endif
