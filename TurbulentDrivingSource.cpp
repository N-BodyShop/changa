#include "TurbulentDrivingSource.h"
#include "dVector.h"

void TurbulentDrivingSource::initializeParameters( Parameters param) {

  klo = param.dKlo;
  khi = param.dKhi;
  kslope = param.dKslope;
  amplitude = param.dTurbAmplitude;
  seed = param.iTurbSeed;
}

void TurbulentDrivingSource::initialize() { 
  double xp = UserGravity::xPeriod; 
  double yp = UserGravity::yPeriod; 
  double zp = UserGravity::zPeriod; 

  double stir_kx[NSTIR],stir_ky[NSTIR],stir_kz[NSTIR]; 

  for (int i=0; i < NSTIR; i++) { 
     stir_kx[i] = TWOPI*i/xp;
     stir_ky[i] = TWOPI*i/yp;
     stir_kz[i] = TWOPI*i/zp;
  }

  //rescale klo, khi
  srand( seed);
  double k0 = TWOPI/fmax(fmax(xp, yp),zp);
  klo *= k0;
  khi *= k0;
  
  double xAmp[NSTIR][NSTIR][NSTIR], xPhi[NSTIR][NSTIR][NSTIR];
  double yAmp[NSTIR][NSTIR][NSTIR], yPhi[NSTIR][NSTIR][NSTIR];
  double zAmp[NSTIR][NSTIR][NSTIR], zPhi[NSTIR][NSTIR][NSTIR];
  double vScale = cmUnit/secUnit;
  // initialized amplitudes and phase
  for( int i = 0; i < NSTIR; i++) {
    for( int j = 0; j < NSTIR; j++) {
      for( int k = 0; k < NSTIR; k++) {
        //amplitude
        xAmp[i][j][k] = 1.0*rand()/RAND_MAX;
        yAmp[i][j][k] = 1.0*rand()/RAND_MAX;
        zAmp[i][j][k] = 1.0*rand()/RAND_MAX;

        //phases
        xPhi[i][j][k] = TWOPI*rand()/RAND_MAX;
        yPhi[i][j][k] = TWOPI*rand()/RAND_MAX;
        zPhi[i][j][k] = TWOPI*rand()/RAND_MAX;

        double total = sqrt(xAmp[i][j][k]*xAmp[i][j][k] + yAmp[i][j][k]*yAmp[i][j][k] + zAmp[i][j][k]*zAmp[i][j][k]);
        double kvec = sqrt(stir_kx[i]*stir_kx[i] + stir_ky[j]*stir_ky[j] + stir_kz[k]*stir_kz[k]);
        double norm = amplitude*amplitude/(vScale*vScale)*pow(kvec/klo,kslope);
        if( kvec < klo || kvec > khi) norm = 0.;

        xAmp[i][j][k] *= norm/total;
        yAmp[i][j][k] *= norm/total;
        zAmp[i][j][k] *= norm/total;

      }
    }
  }
  // zero out accelerations
  for( int ix = 0; ix < TURB_NX; ix++) 
    for( int iy = 0; iy < TURB_NY; iy++) 
      for( int iz = 0; iz < TURB_NZ; iz++) 
        turbAcc[ix][iy][iz] = 0.;

  // create the acceleration field
  for( int ix = 0; ix < TURB_NX; ix++) {
    double x = xp/(TURB_NX-1)*ix; // this is for the fourier grid zero at the corner
    for( int iy = 0; iy < TURB_NY; iy++) {
      double y = yp/(TURB_NY-1)*iy; // this is for the fourier grid zero at the corner
      for( int iz = 0; iz < TURB_NZ; iz++) {
        double z = zp/(TURB_NZ-1)*iz; // this is for the fourier grid zero at the corner
        for( int i = 0; i < NSTIR; i++) {
          double kx = stir_kx[i]*x;
          for( int j = 0; j < NSTIR; j++) {
            double ky = stir_ky[j]*y;
            for( int k = 0; k < NSTIR; k++) {
              double kz = stir_kz[k]*z;
              double skx = sin(kx+xPhi[i][j][k]);
              double ckx = cos(kx+xPhi[i][j][k]);
              double sky = sin(ky+yPhi[i][j][k]);
              double cky = cos(ky+yPhi[i][j][k]);
              double skz = sin(kz+zPhi[i][j][k]);
              double ckz = cos(kz+zPhi[i][j][k]);

              double imtrigterms = ckx*cky*skz + ckx*sky*ckz + skx*cky*ckz - skx*sky*skz;

              // create divergence free stirring field
              turbAcc[ix][iy][iz][0] += (stir_ky[j]*zAmp[i][j][k] - stir_kz[k]*yAmp[i][j][k])*imtrigterms;
              turbAcc[ix][iy][iz][1] += (stir_kz[k]*xAmp[i][j][k] - stir_kx[i]*zAmp[i][j][k])*imtrigterms;
              turbAcc[ix][iy][iz][2] += (stir_kx[i]*yAmp[i][j][k] - stir_ky[j]*xAmp[i][j][k])*imtrigterms;             
            }
          }
        }
      }
    }
  }
}

void TurbulentDrivingSource::interpolate( dVector &pos, dVector &acc) { 

  double xp = UserGravity::xPeriod; 
  double yp = UserGravity::yPeriod; 
  double zp = UserGravity::zPeriod; 

  int xIndex = int(floor((pos[0]/xp+0.5)*TURB_NX));
  int yIndex = int(floor((pos[1]/yp+0.5)*TURB_NY));
  int zIndex = int(floor((pos[2]/zp+0.5)*TURB_NZ));

  int xIndexP1 = int(ceil((pos[0]/xp+0.5)*TURB_NX));
  int yIndexP1 = int(ceil((pos[1]/yp+0.5)*TURB_NY));
  int zIndexP1 = int(ceil((pos[2]/zp+0.5)*TURB_NZ));

  double dx = xp/TURB_NX;
  double dy = yp/TURB_NY;
  double dz = zp/TURB_NZ;

  //rescaled xd to be dimensionless 
  double xd = (pos[0]+0.5*xp - xp*xIndex/TURB_NX)/dx;
  double yd = (pos[1]+0.5*yp - yp*yIndex/TURB_NY)/dy;
  double zd = (pos[2]+0.5*zp - zp*zIndex/TURB_NZ)/dz;

  // use boundary conditions (periodic)
  xIndex += xIndex < 0 ? TURB_NX : 0;
  yIndex += yIndex < 0 ? TURB_NY : 0;
  zIndex += zIndex < 0 ? TURB_NZ : 0;

  xIndexP1 += xIndexP1 < 0 ? TURB_NX : 0;
  yIndexP1 += yIndexP1 < 0 ? TURB_NY : 0;
  zIndexP1 += zIndexP1 < 0 ? TURB_NZ : 0;

  xIndex -= xIndex >= TURB_NX ? TURB_NX : 0;
  yIndex -= yIndex >= TURB_NY ? TURB_NY : 0;
  zIndex -= zIndex >= TURB_NZ ? TURB_NZ : 0;

  xIndexP1 -= xIndexP1 >= TURB_NX ? TURB_NX : 0;
  yIndexP1 -= yIndexP1 >= TURB_NY ? TURB_NY : 0;
  zIndexP1 -= zIndexP1 >= TURB_NZ ? TURB_NZ : 0;

  // trilinear interpolation -- nomenclature is from Wikipedia
  dVector c00 = turbAcc[xIndex][yIndex][zIndex]*(1.-xd) + turbAcc[xIndexP1][yIndex][zIndex]*xd;
  dVector c01 = turbAcc[xIndex][yIndex][zIndexP1]*(1.-xd) + turbAcc[xIndexP1][yIndex][zIndexP1]*xd;
  dVector c10 = turbAcc[xIndex][yIndexP1][zIndex]*(1.-xd) + turbAcc[xIndexP1][yIndexP1][zIndex]*xd;
  dVector c11 = turbAcc[xIndex][yIndexP1][zIndexP1]*(1.-xd) + turbAcc[xIndexP1][yIndexP1][zIndexP1]*xd;

  dVector c0 = c00*(1.-yd) + c10*yd;
  dVector c1 = c01*(1.-yd) + c11*yd;

  acc = c0*(1.-zd) + c1*zd;
}

void TurbulentDrivingSource::getGravity( double t, double x, double y, double z, double& gx, double& gy, double& gz, double &dtg){

  dVector pos(x,y,z);
  dVector acc;

  interpolate( pos, acc);
  gx = acc[0]; gy = acc[1]; gz = acc[2];
  dtg = sqrt( sqrt(x*x + y*y + z*z)/sqrt(gx*gx + gy*gy + gz*gz));
}
