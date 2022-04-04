#ifdef RADIATION
#include "RadiationSource.h"
#include "Opacity.h"
#include <math.h>

RadiationSource* RadiationSource::instance = 0;

void RadiationSource::initialize(){ 
  if(instance == 0) 
    instance = this;
} 

void RadiationSource::initializeParameters( Parameters param) {
  opacityConversion = ((3.0857e21*param.dKpcUnit)*(3.0857e21*param.dKpcUnit));
  ergsccConversion = (param.dErgPerGmUnit*param.dGmPerCcUnit); // to check
  densityConversion = param.dGmPerCcUnit;
  noHydro = param.bRadNoHydro != 0;
  noLorentz = param.bRadNoLorentz != 0;
  
}

double RadiationSource::computeGamma( dVector v, dVector norm) {
  const double cSpeed = HydroUtils::instance->cInCodeUnits;
  dVector beta = v/cSpeed;
  if( noLorentz) beta = 0.;

  double lorentz2 = 1./(1.0 - dotProduct(beta, beta));
  double vnc = 1. - dotProduct(beta, norm);
  return vnc*sqrt(lorentz2);
}

void RadiationSource::applySource( double dt, dVector &pos, double volume, 
                                   varVector &cons, varVector &deltas){
  varVector prim;
  double gamma;
  // compute the new conservatives values -- operator split method.
  varVector newCons; newCons=cons + deltas/volume;
  HydroUtils::con2Prim(gamma, newCons,prim);
  double p, ie, cs, T;
  HydroUtils::EOS( prim, gamma, p, ie, cs, T); // get the temperature
  double radWeightsCM[NRADVARS];
  
  double rCSpeed = HydroUtils::instance->redCInCodeUnits;
  double cSpeed = HydroUtils::instance->cInCodeUnits;

  double J;
  varVector primCM; 
  labToCom( prim, primCM, radWeightsCM);

  // implicit update
  double Tnew = T;
  double DeltaE = 0.;
  //double Trad = pow(4.*M_PI*J/(ARAD*ergsccConversion), 0.25);
  computeNewTgas( dt, J, Tnew, DeltaE, primCM, radWeightsCM);

  double rho = prim[iRho]*densityConversion; // get density  
  double dens = prim[iRho];
  double kappaPlanck = 0., kappaMean = 0., kappaSca = 0.; 
  Opacity::computeTotalOpacity( rho, Tnew, prim, kappaPlanck, kappaMean, kappaSca);

  double sigmaAbsPl =   dens*kappaPlanck/opacityConversion;
  double sigmaAbsMean = dens*kappaMean/opacityConversion;
  double sigmaSca = dens*kappaSca/opacityConversion;
 
  dVector fluxCM, oldFluxCM; fluxCM=0.; oldFluxCM = 0.;
  varVector newPrimLab, newPrimCM; newPrimCM=primCM;
  dVector vel(prim[ivx], prim[ivy], prim[ivz]);

  for( int iFreq = 0; iFreq < NFREQVARS; iFreq++) {
    //double emin = emissivity( T);
    double T4 = Tnew*Tnew*Tnew*Tnew;
    //double T4 = T*T*T*T;
    
    const double invFourPi = 1./(4.*M_PI);
    for( int i = 0; i < NANGVARS; i++) {
      int irad = iFreq*NANGVARS + i;
      int ivar = irad + IRADVARSTART;
      dVector beta(prim[ivx], prim[ivy], prim[ivz]); beta *= 1./cSpeed;
      if( noLorentz) beta = 0.;

      double lorentz2 = 1./(1.0 - dotProduct(beta, beta));
      double vnc = 1. - dotProduct(beta, HydroUtils::instance->radNormals[i]);
      const double Gamma = computeGamma( vel, HydroUtils::instance->radNormals[i]);
      
      double cSigmaAbsMeandt = dt*rCSpeed*sigmaAbsMean*Gamma;
      double cSigmaAbsPldt = dt*rCSpeed*sigmaAbsPl*Gamma;
      double cSigmaScatdt = dt*rCSpeed*sigmaSca*Gamma;
      
      double invJfactor = 1./(1. + cSigmaAbsMeandt + cSigmaScatdt);
      double IoldCM = primCM[ivar];
      double InewCM = (cSigmaScatdt*J + cSigmaAbsPldt*ARAD*invFourPi*T4/ergsccConversion + primCM[ivar])*invJfactor;
      oldFluxCM += 4.*M_PI*HydroUtils::instance->radNormals[irad]*radWeightsCM[irad]*IoldCM;
      fluxCM +=    4.*M_PI*HydroUtils::instance->radNormals[irad]*radWeightsCM[irad]*InewCM;
      newPrimCM[ivar] = InewCM;
    }
  }

  comToLab( newPrimCM, newPrimLab);

  dVector flux; flux = 0.; dVector oldFlux; oldFlux = 0.;

  for( int iFreq = 0; iFreq < NFREQVARS; iFreq++) {
    for( int i = 0; i < NANGVARS; i++) {
      int irad = iFreq*NANGVARS + i;
      int ivar = irad + IRADVARSTART;
      deltas[ivar] += volume*(newPrimLab[ivar] - prim[ivar]);

      oldFlux += 4.*M_PI*HydroUtils::instance->radNormals[irad]*HydroUtils::instance->radWeights[irad]*prim[ivar];
      flux +=    4.*M_PI*HydroUtils::instance->radNormals[irad]*HydroUtils::instance->radWeights[irad]*newPrimLab[ivar];
    }
  }

  if(noHydro) return;

  // include changes in matter source terms.
  deltas[iE] += DeltaE*volume;

  dVector dp = (oldFlux-flux)/rCSpeed*volume;
  dVector oldP( cons[ipx], cons[ipy], cons[ipz]);
  oldP *= volume; 
  oldP[0] += deltas[ipx]; oldP[1] += deltas[ipy]; oldP[2] += deltas[ipz];
  double oldKE = 0.5*dotProduct(oldP,oldP)/(dens*volume + deltas[iRho]);
  deltas[ipx] += dp[0]; deltas[ipy] += dp[1]; deltas[ipz] += dp[2];
  dVector newP = oldP + dp;
  double newKE = 0.5*dotProduct(newP,newP)/(dens*volume + deltas[iRho]);
  deltas[iE] += newKE - oldKE;

}


//do the implicit integration
void RadiationSource::computeNewTgas( double dt, double &J, double &Tgas, double &deltaE, varVector &prim, double radWeights[]) {
  // compute the old J from intensities
  const double fourPi = 4*M_PI;
  
#ifdef COOLING_AD
  const double Told = Tgas;
  const double gamma = HydroUtils::defaultGamma, mu = HydroUtils::defaultMu;
  double rho = prim[iRho]*densityConversion;
  double dens = prim[iRho];
  double rCSpeed = HydroUtils::instance->redCInCodeUnits;
  double cSpeed = HydroUtils::instance->cInCodeUnits;

  double RTerm = rho*RGAS/mu/(gamma-1.)/ergsccConversion;

  dVector beta(prim[ivx], prim[ivy], prim[ivz]); beta *= 1./HydroUtils::instance->cInCodeUnits;
  dVector vel(prim[ivx], prim[ivy], prim[ivz]);
  if( noLorentz) beta = 0.;

  double lorentz2 = 1./(1.0 - dotProduct(beta, beta));
  J = 0.; 

  double kappaPlanck = 0., kappaMean = 0., kappaSca = 0.; 
  Opacity::computeTotalOpacity( rho, Tgas, prim, kappaPlanck, kappaMean, kappaSca);

  double sigmaAbsPl =   dens*kappaPlanck/opacityConversion;
  double sigmaAbsMean = dens*kappaMean/opacityConversion;
  double sigmaSca =     dens*kappaSca/opacityConversion;

  double SigmaI = 0., SigmaGamma = 0.;
  for( int iFreq = 0; iFreq < NFREQVARS; iFreq++) {
    for( int i = 0; i < NANGVARS; i++) {
      J += radWeights[i]*prim[iFreq*NANGVARS+i+IRADVARSTART];
      double vnc = 1. - dotProduct(beta, HydroUtils::instance->radNormals[i]);
      //double Gamma = vnc*sqrt(lorentz2);
      const double Gamma = computeGamma( vel, HydroUtils::instance->radNormals[i]);
      // get angular dependent gamma 
      double csigmaGammadt = rCSpeed*(sigmaAbsMean + sigmaSca)*Gamma*dt;
      SigmaI += radWeights[i]*prim[iFreq*NANGVARS+i+IRADVARSTART]/(1.+csigmaGammadt);
      SigmaGamma += radWeights[i]*Gamma/(1. + csigmaGammadt);
    }
  }

  if(!Opacity::isTemperatureDependentOpacity()) {  
 
    double csigmadt = cSpeed*sigmaAbsMean*dt;
    double csigmaPldt = cSpeed*sigmaAbsPl*dt;
    double credsigmadt = rCSpeed*sigmaAbsMean*dt;
    double credsigmaScadt = rCSpeed*sigmaSca*dt;  
    double Jfactor = 1. - (credsigmadt + credsigmaScadt)*SigmaGamma;
    
    // from hydro.py coef4 = csigmaPldt/RTerm*aRad*(1. - credsigmadt*SigmaGamma)

    double coef4 = csigmaPldt/RTerm*ARAD*(1. - csigmadt/csigmaPldt*credsigmadt*SigmaGamma/Jfactor)/ergsccConversion;
    double tconst = -csigmadt/RTerm*SigmaI/Jfactor*fourPi - Told;
    int badCell = FourthPolyRoot( coef4, tconst, Tgas);
    
    if( badCell) {
      Tgas = Told;
      //J = Jold;
      deltaE = 0.;
      return;
    }
    const double T4 = Tgas*Tgas*Tgas*Tgas;
    double Jold = J;
    
    J = (credsigmadt*ARAD/fourPi*T4*SigmaGamma/ergsccConversion + SigmaI)/Jfactor;
    deltaE = RTerm*(Tgas - Told);
  }
  else { // temperature dependent
    // define temperature function
    auto Tfunction = [&]( double temp) -> double
      {
        double kappaPlanck = 0., kappaMean = 0., kappaSca = 0.; 
        Opacity::computeTotalOpacity( rho, temp, prim, kappaPlanck, kappaMean, kappaSca);

        double sigmaAbsPl =   rho*kappaPlanck/opacityConversion;
        double sigmaAbsMean = rho*kappaMean/opacityConversion;

        double csigmadt = cSpeed*sigmaAbsMean*dt;
        double csigmaPldt = cSpeed*sigmaAbsPl*dt;
        double aT4 = ARAD*temp*temp*temp*temp/ergsccConversion;
        double Jfactor = 1. + csigmadt;

        return RTerm*(temp - Told) + csigmaPldt*aT4*(1. - csigmadt*SigmaGamma) - csigmadt*SigmaI*fourPi; 
      };

    // perform Newton-Ralfson iteration
    // first find the bound of the radiation temperature and gas temperature
    const double Jold = J;
    double Trad = pow( Jold*ergsccConversion/ARAD, 0.25);
    double T1 = 0.5*fmin( Trad, Told), T2 = 1.5*fmax( Trad, Told);

    double f1 = Tfunction( T1), f2 = Tfunction( T2);
  
    // check if it is bounded 
    if( f1/f2 > 0.) { // it is not bounded
      printf("not bounded %5.3e %5.3e %5.3e %5.3e %5.3e\n", f1, f2, T1, T2, Trad);
      Tgas = Told;
      //J = Jold;
      deltaE = 0.;
      return;
    } 

    const double h = 0.001, errTol = 1e-8;
    double tol = 1.; 
    const int maxIter = 80;
    int iter = 0;
    double T = fabs(f1) < fabs(f2) ? T1 : T2; // pick the value closer to zero   
  
    while (tol > errTol && iter++ < maxIter) {
      double f = Tfunction( T);
      double dfdT = Tfunction( T*(1. + h)) - f; dfdT /= h*T;
      double Tnew = T - 0.7*f/dfdT;
      tol = fabs(Tnew - T)/Tnew;
      T = Tnew;
    }
    if( iter > maxIter) {
      printf( "exceeded iterations: %d tol = %5.3e T=%5.3e Told=%5.3e\n", iter, tol, T, Told);
    }

    Tgas = T;

    // compute everything again
    double kappaPlanck = 0., kappaMean = 0., kappaSca = 0.; 
    Opacity::computeTotalOpacity( rho, Tgas, prim, kappaPlanck, kappaMean, kappaSca);

    const double T4 = Tgas*Tgas*Tgas*Tgas;
    double sigmaAbsPl =   dens*kappaPlanck/opacityConversion;
    double sigmaAbsMean = dens*kappaMean/opacityConversion;
    double sigmaSca = dens*kappaSca/opacityConversion;
 
    double csigmadt = rCSpeed*sigmaAbsMean*dt;
    double csigmaPldt = rCSpeed*sigmaAbsPl*dt;

    double Jfactor = 1. - rCSpeed*(sigmaAbsMean+sigmaSca)*dt*SigmaGamma;
    
    //from hydro.py J = credsigmadt*aRad/fourPi*T4*SigmaGamma + SigmaI

    J = (SigmaGamma*csigmaPldt*ARAD/fourPi*T4/ergsccConversion + SigmaI)/Jfactor;
    
    deltaE = RTerm*(Tgas - Told);
  }


#endif
}

// transfrom from labToCom
void RadiationSource::labToCom( varVector &prim, varVector &primCM, double radWeightsCM[]) {
  primCM = prim;
  dVector vel(prim[ivx], prim[ivy], prim[ivz]);

  double norm = 0;
  for( int i = 0; i < NANGVARS; i++) {
    const double Gamma = computeGamma( vel, HydroUtils::instance->radNormals[i]);
    radWeightsCM[i] = HydroUtils::instance->radWeights[i]/(Gamma*Gamma);
    norm += radWeightsCM[i];
  }

  // now renormalize
  for( int i = 0; i < NANGVARS; i++)
    radWeightsCM[i] /= norm;

  if( noLorentz) return;


  dVector beta = vel/HydroUtils::instance->cInCodeUnits;
  double lorentz2 = 1./(1.0 - dotProduct(beta, beta));
  for( int iFreq = 0; iFreq < NFREQVARS; iFreq++) {
    for( int i = 0; i < NANGVARS; i++) {
      double vnc = 1. - dotProduct(beta, HydroUtils::instance->radNormals[i]);
      double angCoeff = vnc*vnc*lorentz2; 
      double transCoeff = angCoeff*angCoeff;
      int ivar = iFreq*NANGVARS + i + IRADVARSTART;
      primCM[ivar] = transCoeff*prim[ivar];
    }
  }
}

void RadiationSource::comToLab( varVector &primCM, varVector &prim) {
  prim = primCM;
  if( noLorentz) return;
  dVector beta(prim[ivx], prim[ivy], prim[ivz]); beta *= 1./HydroUtils::instance->cInCodeUnits;
  double lorentz2 = 1./(1.0 - dotProduct(beta, beta));
  for( int iFreq = 0; iFreq < NFREQVARS; iFreq++) {
    for( int i = 0; i < NANGVARS; i++) {
      double vnc = 1. - dotProduct(beta, HydroUtils::instance->radNormals[i]);
      double angCoeff = vnc*vnc*lorentz2; 
      double transCoeff = angCoeff*angCoeff;
      int ivar = iFreq*NANGVARS + i + IRADVARSTART;
      prim[ivar] = primCM[ivar]/transCoeff;
    }
  }
}

// Exact solution for fourth order polynomical with the format
// coef4 * x^4 + x + tconst == 0

int RadiationSource::FourthPolyRoot(const double coef4, const double tconst, double &root)
{

// First, get the double root of
// z^3-4*tconst/coef4 * z - 1/coef4^2==0
  double asquar = coef4 * coef4;
  double acubic = coef4 * asquar;
  double ccubic = tconst * tconst * tconst;
  double delta1 = 0.25 - 64.0 * ccubic * coef4/27.0;
  if(delta1 < 0.0) return -1;
  else delta1 = sqrt(delta1);
  double zroot = 0.0;
  if(delta1 > 1.e11){
    // to avoid small number cancellation
    zroot = pow(delta1,-2.0/3.0)/3.0;
  }else{
    zroot = (pow(0.5 + delta1, 1.0/3.0) -
               pow(-0.5 + delta1, 1.0/3.0));
  }

  zroot *= pow(coef4,-2.0/3.0);
  
  double rcoef = sqrt(zroot);
  double delta2 = -zroot + 2.0/(coef4*rcoef);
  if(delta2 < 0.0) return -1;
  else delta2 = sqrt(delta2);
  root = 0.5 * (delta2 - rcoef);
  if(root < 0.0) return -1;

  return 0;
}
#endif // RADIATION
