#ifdef RADIATION
#include "UserRadiationSources.h"
#include "RadiationSource.h"
#include "Opacity.h"
#include "HydroParticle.h"
#include <math.h>

void UserRadiationSources::initialize(){} 

void UserRadiationSources::initializeParameters( Parameters param) {
  position[0][0] = param.dRadPosx1/(3.0857e21*param.dKpcUnit);
  position[0][1] = param.dRadPosy1/(3.0857e21*param.dKpcUnit);
  position[0][2] = param.dRadPosz1/(3.0857e21*param.dKpcUnit);

  position[1][0] = param.dRadPosx2/(3.0857e21*param.dKpcUnit);
  position[1][1] = param.dRadPosy2/(3.0857e21*param.dKpcUnit);
  position[1][2] = param.dRadPosz2/(3.0857e21*param.dKpcUnit);

  normals[0][0] = param.dRadNormx1/(3.0857e21*param.dKpcUnit);
  normals[0][1] = param.dRadNormy1/(3.0857e21*param.dKpcUnit);
  normals[0][2] = param.dRadNormz1/(3.0857e21*param.dKpcUnit);
  double norm = sqrt(dotProduct(normals[0],normals[0]));
  normals[0] /= fmax(norm,1e-30);

  normals[1][0] = param.dRadNormx2/(3.0857e21*param.dKpcUnit);
  normals[1][1] = param.dRadNormy2/(3.0857e21*param.dKpcUnit);
  normals[1][2] = param.dRadNormz2/(3.0857e21*param.dKpcUnit);
  norm = sqrt(dotProduct(normals[1],normals[1]));
  normals[1] /= fmax(norm,1e-30);
  mu0 = param.dRadMu0;
  radiusSource = param.dRadRadius/(3.0857e21*param.dKpcUnit);
  double Teff = param.dRadTeff; double T4 = Teff*Teff*Teff*Teff;
  double cSpeed = HydroUtils::instance->redCInCodeUnits;
  y0 = param.dRadY0;
  deltaY = param.dRadDeltaY;

  pointSource = param.bRadPointSource != 0;
  slabSource = param.bRadSlabSource != 0;
  tau = param.dTauBoundary;
  I0 = ARAD*T4/(param.dErgPerGmUnit*param.dGmPerCcUnit)/(4.*M_PI);
  F0 = ARAD*T4/(param.dErgPerGmUnit*param.dGmPerCcUnit);
}

void UserRadiationSources::applySource( double dt, dVector &pos, double volume, 
                                    varVector &cons, varVector &deltas){
  varVector prim;
  double gamma;
  
  if( !pointSource && !slabSource) return;
  
  HydroUtils::con2Prim(gamma, cons,prim);

  double rho = cons[iRho]*RadiationSource::instance->densityConversion;
  double T = 100.;
  double kappaAbsPl, kappaAbsMean, kappaSca;
  Opacity::computeTotalOpacity( rho, T, prim, kappaAbsPl, kappaAbsMean, kappaSca);

  const double opacityConversion = RadiationSource::instance->opacityConversion; 
  kappaAbsPl *= opacityConversion;
  kappaAbsMean *= opacityConversion;
  kappaSca *= opacityConversion;
  double cSpeed = HydroUtils::instance->redCInCodeUnits;
  double dtau = fmin(dt/tau, 0.5);

  for (int irad = 0; irad < numRadSources; irad++) {
    dVector dp = position[irad] - pos;
    if ((pointSource && (dp[0] * dp[0] + dp[1] * dp[1]) < radiusSource * radiusSource)
      || (slabSource && fabs(pos[1] - y0) < deltaY) && 
         (pos[1]-y0)*normals[irad][1] > 0) { // make sure that the radiation is directional
      double Fy = 0.;
      if( constantFlux ) { //First compute the flux
        for( int iFreq = 0; iFreq < NFREQVARS; iFreq++) {
          double weights = 0.;
          for( int iAng = 0; iAng < NANGVARS; iAng++) {
            int i =  iAng + iFreq*NANGVARS;
            double ynorm = dotProduct(HydroUtils::instance->radNormals[i],normals[irad]);
            Fy += 4.*M_PI*ynorm*cons[i + IRADVARSTART]*HydroUtils::instance->radWeights[i];
          }
        }
      }
      for (int i = 0; i < NRADVARS; i++) {
        double mu = dotProduct(normals[irad], HydroUtils::instance->radNormals[i]);
        //printf("Apply Source %5.3e %5.3e\n", mu, mu0);
        if (mu > mu0) {
          if (constantFlux) {
            deltas[i + IRADVARSTART] += (F0 - Fy)*volume*dtau/(4.*M_PI);
          }
          else {
            deltas[i + IRADVARSTART] += -cons[i + IRADVARSTART] * volume * dtau;
            deltas[i + IRADVARSTART] += I0 * volume * dtau;
          }
          //deltas[i+IRADVARSTART] += dt*volume*Idot*rho*kappaAbsPl;
          //          printf("Apply Source %5.3e\n", dt*volume*Idot);
        }
      }
    }
  }


/*   if(slabSource) {
    for(int irad = 0; irad < numRadSources; irad++) {
      if( fabs(pos[1] - y0) < deltaY) {
        for(int i = 0; i < NRADVARS; i++) {
          double mu = dotProduct(normals[irad],HydroUtils::instance->radNormals[i]);
          //printf("Apply Source %5.3e %5.3e\n", mu, mu0);
          if( mu > mu0) {
            if( constantFlux) {
              deltas[i+IRADVARSTART] += I0*volume*dt;  
            }
            else {
              deltas[i+IRADVARSTART] += -cons[i+IRADVARSTART]*volume*dtau;
              deltas[i+IRADVARSTART] += I0*volume*dtau;
            }
            //deltas[i+IRADVARSTART] += dt*volume*Idot*rho*kappaAbsPl;
//          printf("Apply Source %5.3e\n", dt*volume*Idot);
          }  
        }
      }
    }
  }
 */
}

#endif //RADIATION
