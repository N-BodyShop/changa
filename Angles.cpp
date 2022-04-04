#ifdef RADIATION
#include "dVector.h"
#include "HydroParticle.h"
#include <math.h>
#include <stdlib.h>


void HydroUtils::initRadNormals() {
  legendre_equal_weight();
}

/*----------------------------------------------------------------------------*/
/*! \fn void init_angles(RadGridS *pRG, const int qmeth, const int NMU_QUAD)
 *  \brief Initialize angles/angle quadratures using Gauss-Legendre
 *  quadrature */

void HydroUtils::legendre_equal_weight()
{
  dVector mutmp[NANG_QUAD], mu[NANGVARS];
  double mutmp1d[2*NMU_QUAD], wtmp1d[2*NMU_QUAD];
  double weightQuad[NANG_QUAD], weights[NANGVARS];
  double dphi, sintheta;
  const double pi = 3.14159265359;

  HydroUtils::gauleg(-1.0, 1.0, mutmp1d, wtmp1d, 2*NMU_QUAD);

  int k = 0;
  for(int i=0; i<NMU_QUAD; i++) {
    dphi = 0.5 * pi / (double) (2*(1+i));
    int l=2*NMU_QUAD-1-i;
    sintheta = sqrt(1.0 - mutmp1d[l] * mutmp1d[l]);
    for(int j=0; j<=i; j++) {
      mutmp[k][0] = sintheta * cos( dphi * (double)(2*j+1));
      mutmp[k][1] = sintheta * sin( dphi * (double)(2*j+1));
      mutmp[k][2] = mutmp1d[l];
      weightQuad[k] = wtmp1d[l] / (double)(i+1);
      printf("wmu %d %d %g %g\n",k,i,weightQuad[k],wtmp1d[l]);
      printf("mu: %d %d %g %g %g %g\n",k,i,mutmp[k][0],mutmp[k][1],mutmp[k][2],dphi * (double)(2*j+1));
      k++;
    }
  }
  
  for (int i=0; i<NANG_QUAD; i++) {
    weightQuad[i] *= 0.125;

    for (int j=0; j<2; j++) {
      for (int k=0; k<2; k++) {
        for (int l=0; l<2; l++) {
          int m=4*j+2*k+l; // this is the quadrant
          int iAngle = m*NANG_QUAD + i;
          if (l == 0)
            mu[iAngle][0] =  mutmp[i][0];
          else
            mu[iAngle][0] = -mutmp[i][0];
          if (k == 0)
            mu[iAngle][1] =  mutmp[i][1];
          else
            mu[iAngle][1] = -mutmp[i][1];
          if (j == 0)
            mu[iAngle][2] =  mutmp[i][2];
          else
            mu[iAngle][2] = -mutmp[i][2];
          weights[iAngle] = weightQuad[i];
        }
      }
    }
  }
    

  for( int iAngle=0; iAngle < NANGVARS; iAngle++) {
    for(int iFreq=0; iFreq < NFREQVARS; iFreq++) {
      int iRad = iAngle*NFREQVARS + iFreq;
      //printf("%d %5.3e\n", iAngle, dotProduct(mu[iAngle], mu[iAngle]));
      radNormals[iRad] = mu[iAngle];
      radWeights[iRad] = weights[iAngle];
    }
  }

  return;
}  

/*----------------------------------------------------------------------------*/
/*! \fn void gauleg(double x1, double x2,  double *x, double *w, int n)
 * gauss-legendre weight routine from numerical recipes */
void HydroUtils::gauleg(double x1, double x2,  double *x, double *w, int n)
{

  double eps = 3.0e-14;
  double xm, xl, z, z1;
  double p1, p2, p3, pp;
  int m, i, j;
  const double pi = 3.14159265359;

  m = (n + 1) / 2;
  xm = 0.5 * (x2 + x1);
  xl = 0.5 * (x2 - x1);

  for (i=1; i<=m; i++) {
    z = cos(pi * ((double)i - 0.25) / ((double)n + 0.5));
    do {
      p1=1.0;
      p2=0.0;
      for(j=1; j<=n; j++) {
        p3 = p2;
        p2 = p1;
        p1 = ((2.0 * (double)j - 1.0) * z * p2 - ((double)j - 1.0) * p3) / (double)j;
      }
      pp = (double)n * (z * p1 - p2) / (z * z - 1.0);
      z1 = z;
      z = z1 - p1 / pp;
    }  while(fabs(z - z1) > eps);
    x[i-1] = xm - xl * z;
    x[n-i] = xm + xl * z;
    w[i-1] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
    w[n-i] = w[i-1];
  }

}

bool HydroUtils::RadiationClean( varVector &prim) {
  bool cleaned = false;
  for( int iFreq = 0; iFreq < NFREQVARS; iFreq++) {
    for( int i = 0; i < NANGVARS; i++) {
      if( prim[iFreq*NANGVARS+i+IRADVARSTART] < 0.) {
        cleaned = true;
        prim[iFreq*NANGVARS+i+IRADVARSTART] = 0.; 
      }
    }
  }
  return cleaned;
}

bool HydroUtils::RadiationCleanCons( varVector &con) {
  bool cleaned = false;

  for( int iFreq = 0; iFreq < NFREQVARS; iFreq++) {
    for( int i = 0; i < NANGVARS; i++) {
      if( con[iFreq*NANGVARS+i+IRADVARSTART] < 0.) {
        cleaned = true;
        con[iFreq*NANGVARS+i+IRADVARSTART] = 0.; 
      }
    }
  }
  cleaned = true;
}

dVector HydroParticle::getRadFlux(){
  const double fourPi = 4*M_PI;
  dVector flux = 0.;
  for( int iFreq = 0; iFreq < NFREQVARS; iFreq++) {
    for( int i = 0; i < NANGVARS; i++) 
      flux += HydroUtils::instance->radNormals[i]*HydroUtils::instance->radWeights[i]*primitives[iFreq*NANGVARS+i+IRADVARSTART];
  }
  return flux*(fourPi*HydroUtils::instance->cInCodeUnits);
}

double  HydroParticle::getRadEnergy(){
  const double fourPi = 4*M_PI;
  double energy = 0.;
  for( int iFreq = 0; iFreq < NFREQVARS; iFreq++) {
    for( int i = 0; i < NANGVARS; i++) 
      energy += HydroUtils::instance->radWeights[i]*primitives[iFreq*NANGVARS+i+IRADVARSTART];
  }
  return fourPi*energy;
} 


#endif /* RADIATION_TRANSFER */
