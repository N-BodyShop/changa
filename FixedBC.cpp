#ifdef VORONOI
#include "FixedBC.h"
#include "HydroParticle.h"
#include <math.h>
#include "physconst.h"


void FixedBC::initialize(){} 

void FixedBC::initializeParameters( Parameters param) {
  for( int i = 0; i < 6; i++) {
    boundaries[i][0] = param.dBoundaries[i].x/(3.0857e21*param.dKpcUnit);
    boundaries[i][1] = param.dBoundaries[i].y/(3.0857e21*param.dKpcUnit);
    boundaries[i][2] = param.dBoundaries[i].z/(3.0857e21*param.dKpcUnit);
    rho0[i] = param.dRhoBoundary[i];
    norm2[i] = dotProduct(boundaries[i],boundaries[i]);
    numBoundaries = param.iNumBoundaries <= 6 ? param.iNumBoundaries : 6;
    useRadiation[i] = (param.bUseRadiation[i] == 1);
  }
  tau = param.dTauBoundary;

  // spherical boundary
  center[0] = param.dSphericalBC.x/(3.0857e21*param.dKpcUnit); 
  center[1] = param.dSphericalBC.y/(3.0857e21*param.dKpcUnit); 
  center[2] = param.dSphericalBC.z/(3.0857e21*param.dKpcUnit);
  rhoSph = param.dRhoSpherical;
  radius = param.dSphericalBCRadius/(3.0857e21*param.dKpcUnit);
  tempSph = param.dTempSpherical;

#ifdef RADIATION
  double Teff = param.dRadTeff; double T4 = Teff*Teff*Teff*Teff;
  I0 = ARAD*T4/(param.dErgPerGmUnit*param.dGmPerCcUnit)/(4.*M_PI);

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

#endif
}

void FixedBC::applySource( double dt, dVector &pos, double volume, 
                                    varVector &cons, varVector &deltas){
  if( numBoundaries < 1 && radius < 0) return;

  double dtau = fmin(dt/tau, 0.5);
  // cartesian
  if( numBoundaries >= 1) {
    double gamma;
    varVector prim;
    
    HydroUtils::con2Prim(gamma, cons, prim);

    double rho = cons[iRho];
    
    for( int iBound = 0; iBound < numBoundaries; iBound++) {

      double n2 = dotProduct( boundaries[iBound], pos);
      if( n2 < norm2[iBound]) continue; // do not apply sources  
      double gamma;
      double DeltaRho = (rho0[iBound]-rho)*dtau;
      dVector DeltaP(-cons[ipx]*dtau, -cons[ipy]*dtau, -cons[ipz]*dtau);
      dVector oldP(cons[ipx], cons[ipy], cons[ipz]);
      dVector newP = oldP + DeltaP;

      deltas[iRho] += DeltaRho*volume;
      deltas[ipx] += DeltaP[0]*volume; // zero out momentum
      deltas[ipy] += DeltaP[1]*volume; 
      deltas[ipz] += DeltaP[2]*volume; 
      double oldKE = 0.5*dotProduct(oldP,oldP)/rho;
      double newKE = 0.5*dotProduct(newP,newP)/(rho+DeltaRho);

      double ie = prim[iIntE];
      double s = prim[iS];
      deltas[iE] += (DeltaRho*ie + newKE - oldKE)*volume; // set energy to small amount.
      deltas[iS] += DeltaRho*s*volume;

    //fix radiation
  #ifdef RADIATION 
      for(int i = 0; i < NRADVARS; i++) {
        deltas[i+IRADVARSTART] += -cons[i+IRADVARSTART]*volume*dtau;
        if( useRadiation[iBound]) {
          for( int irad = 0; irad < numRadSources; irad++) { 
            double mu = dotProduct(normals[irad],HydroUtils::instance->radNormals[i]);
            if( mu > mu0) 
              deltas[i+IRADVARSTART] += I0*volume*dtau;
          }
        }
      }
  #endif
    }
  }
  else if( radius > 0) {
    dVector rad = pos-center;
    if( length(rad) > radius) return;

    double gamma;
    varVector prim;
    HydroUtils::con2Prim(gamma, cons, prim);

    double rho = prim[iRho];

    double DeltaRho = (rhoSph-rho)*dtau;
    dVector DeltaP(-cons[ipx]*dtau, -cons[ipy]*dtau, -cons[ipz]*dtau);
    dVector oldP(cons[ipx], cons[ipy], cons[ipz]);
    dVector newP = oldP + DeltaP;

    deltas[iRho] += DeltaRho*volume;
    deltas[ipx] += DeltaP[0]*volume; // zero out momentum
    deltas[ipy] += DeltaP[1]*volume; 
    deltas[ipz] += DeltaP[2]*volume; 
    double oldKE = 0.5*dotProduct(oldP,oldP)/rho;
    double newKE = 0.5*dotProduct(newP,newP)/(rho+DeltaRho);

    double ie = prim[iIntE];
    double s = prim[iS];
    double p = s*pow(rho, gamma), ieNew, cs, temp;
     
    if( tempSph > 0) 
      HydroUtils::EOS( prim, gamma, p, ieNew, cs, tempSph, false);
    else 
      ieNew = ie;

    double sNew = p/pow(rho,gamma);
    deltas[iE] += ((rhoSph*ieNew - rho*ie)*dtau + newKE - oldKE)*volume; // set energy to small amount.
    deltas[iS] += (rhoSph*sNew - rho*s)*dtau*volume;

  }



}
#endif //VORONOI
