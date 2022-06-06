#ifdef VORONOI
#include "HydroParticle.h"
#include "Neighbor.h"
#include <unordered_map>
#include <cmath>
#include "ckio.h"
#include "other-codes/athena/athena_defs.h"
#include "UserSources.h"

// #include<math>
using namespace std;



void Face::computeRotatedVectors()
{
  dVector xaxis; xaxis = 0.; xaxis[0] = 1.;
  dVector yaxis; yaxis = 0.; yaxis[1] = 1.;

  // dVector zaxis(NDIM); zaxis = 0.; zaxis[2] = 1.;

  n0  = normal;
  n0 /= length(n0);

  // make sure not colinear
  crossProduct(n0, xaxis, n1);

  if (length(n1) > normTolerance) n1 /= length(n1);
  else n1 = yaxis;

  // n2.resize(NDIM);
  crossProduct(n0, n1, n2);
  n2 /= length(n2);

  // calculate the inverse rotation

  crossProduct(n1, n2, n0p);
  crossProduct(n2, n0, n1p);
  crossProduct(n0, n1, n2p);

  dVector test;
  crossProduct(n1, n2, test);
  invDetA = 1. / dotProduct(n0, test);
}

void Face::rotate(varVector& v)
{
  varVector rotatedV;

  rotatedV = v;
  rotatedV[ipx] = (v[ipx] * n0[0] + v[ipy] * n0[1] + v[ipz] * n0[2]);
  rotatedV[ipy] = (v[ipx] * n1[0] + v[ipy] * n1[1] + v[ipz] * n1[2]);
  rotatedV[ipz] = (v[ipx] * n2[0] + v[ipy] * n2[1] + v[ipz] * n2[2]);

#ifdef MHD
  rotatedV[iBx] = (v[iBx] * n0[0] + v[iBy] * n0[1] + v[iBz] * n0[2]);
  rotatedV[iBy] = (v[iBx] * n1[0] + v[iBy] * n1[1] + v[iBz] * n1[2]);
  rotatedV[iBz] = (v[iBx] * n2[0] + v[iBy] * n2[1] + v[iBz] * n2[2]);

  rotatedV[iAx] = (v[iAx] * n0[0] + v[iAy] * n0[1] + v[iAz] * n0[2]);
  rotatedV[iAy] = (v[iAx] * n1[0] + v[iAy] * n1[1] + v[iAz] * n1[2]);
  rotatedV[iAz] = (v[iAx] * n2[0] + v[iAy] * n2[1] + v[iAz] * n2[2]);
#endif

  v = rotatedV;
}

void Face::reverse(varVector& v) {
  varVector rotatedV;

  rotatedV = v;
  rotatedV[ipx] = (v[ipx] * n0p[0] + v[ipy] * n1p[0] + v[ipz] * n2p[0]) * invDetA;
  rotatedV[ipy] = (v[ipx] * n0p[1] + v[ipy] * n1p[1] + v[ipz] * n2p[1]) * invDetA;
  rotatedV[ipz] = (v[ipx] * n0p[2] + v[ipy] * n1p[2] + v[ipz] * n2p[2]) * invDetA;

#ifdef MHD

  rotatedV[iBx] = (v[iBx] * n0p[0] + v[iBy] * n1p[0] + v[iBz] * n2p[0]) * invDetA;
  rotatedV[iBy] = (v[iBx] * n0p[1] + v[iBy] * n1p[1] + v[iBz] * n2p[1]) * invDetA;
  rotatedV[iBz] = (v[iBx] * n0p[2] + v[iBy] * n1p[2] + v[iBz] * n2p[2]) * invDetA;

  rotatedV[iAx] = (v[iAx] * n0p[0] + v[iAy] * n1p[0] + v[iAz] * n2p[0]) * invDetA;
  rotatedV[iAy] = (v[iAx] * n0p[1] + v[iAy] * n1p[1] + v[iAz] * n2p[1]) * invDetA;
  rotatedV[iAz] = (v[iAx] * n0p[2] + v[iAy] * n1p[2] + v[iAz] * n2p[2]) * invDetA;

#endif

  v = rotatedV;
}

Face::Face(int            activeRung,
           dVector      & c,
           dVector      & n,
           double         a,
           HydroParticle *center,
           HydroParticle *neighbor,
           bool           debug = false,
           bool           firstOrderReconstruction = false) {
  verbose   = debug;
  centroid  = c;
  normal    = n / length(n);
  validFace = true;
  dVector origin = center->getPosition();

  if (length(n) < tolerance) {
    validFace = false;

    // return;
  }

  iRung = max(center->getHalfRung(), neighbor->getHalfRung());
  if( activeRung != iRung)  //This face is irrelevant
    isActiveFace = false;  //But keep going to prediction steps
  else  
    isActiveFace = true;

  area = a;

  dVector wL             = center->getVelocity();
  dVector wR             = neighbor->getVelocity();
  dVector rL             = center->getPosition();
  dVector rR             = neighbor->getPosition();

  HydroUtils::applyBC( rL, rR);
  dVector rRMinusrL      = rR - rL;
  double  rRMinusrL_norm = length(rRMinusrL);

  dVector wprime;
  dVector wLMinuswR = (wL-wR)/rRMinusrL_norm;
  dVector fx = centroid - (rR - rL) * 0.5;
  debugData1 = rRMinusrL;
  debugData2 = centroid+rL;
  wprime = rRMinusrL/rRMinusrL_norm*dotProduct(wLMinuswR, fx);
  velocity = wprime + (wL + wR) * 0.5;
//  center->reconstruct(this, dt, origin, neighbor, consL, consL1);
//  neighbor->reconstruct(this, dt, origin, center, consR, consR1);

  center->reconstruct(this, origin, neighbor, consL, consL1, firstOrderReconstruction);
  neighbor->reconstruct(this, origin, center, consR, consR1, firstOrderReconstruction);
  
  
  center->getQ( Q);
  consLOrig = consL; consROrig = consR; 
  consL1Orig = consL1; consR1Orig = consR1;

  if (verbose) {
    varVector primL, primR;
    double newGamma;
    HydroUtils::con2Prim( newGamma, consL, primL, false);
    HydroUtils::con2Prim( newGamma, consR, primR, false);
    printf("reconstructed flux: %d %d %5.3e (%5.3e %5.3e %5.3e) %5.3e %5.3e (%5.3e %5.3e %5.3e) %5.3e\n",
           center->iOrder, neighbor->iOrder, primL[0],primL[1], primL[2], primL[3], primL[4],primR[0],primR[1], primR[2], primR[3], primR[4]);
    center->getPrimitives(primL); neighbor->getPrimitives(primR);
    printf("original flux: %d %d %5.3e (%5.3e %5.3e %5.3e) %5.3e %5.3e (%5.3e %5.3e %5.3e) %5.3e\n",
           center->iOrder, neighbor->iOrder, primL[0],primL[1], primL[2], primL[3], primL[4],primR[0],primR[1], primR[2], primR[3], primR[4]);
    center->getConservedPred(consL); neighbor->getConservedPred(consR);
    HydroUtils::con2Prim( newGamma, consL, primL, false);
    HydroUtils::con2Prim( newGamma, consR, primR, false);
    printf("predicted flux: %d %d %5.3e (%5.3e %5.3e %5.3e) %5.3e %5.3e (%5.3e %5.3e %5.3e) %5.3e\n",
           center->iOrder, neighbor->iOrder, primL[0],primL[1], primL[2], primL[3], primL[4],primR[0],primR[1], primR[2], primR[3], primR[4]);
  }

  double rhoL = consL[iRho];
  double oldKEL = 0.5*(consL[ipx]*consL[ipx] + consL[ipy]*consL[ipy] + consL[ipz]*consL[ipz])/rhoL;
  consL[ipx] -= rhoL * velocity[0];
  consL[ipy] -= rhoL * velocity[1];
  consL[ipz] -= rhoL * velocity[2];
  double newKEL = 0.5*(consL[ipx]*consL[ipx] + consL[ipy]*consL[ipy] + consL[ipz]*consL[ipz])/rhoL;
  consL[iE] += newKEL - oldKEL;

  double rhoR = consR[iRho];
  double oldKER = 0.5*(consR[ipx]*consR[ipx] + consR[ipy]*consR[ipy] + consR[ipz]*consR[ipz])/rhoR;
  consR[ipx] -= rhoR * velocity[0];
  consR[ipy] -= rhoR * velocity[1];
  consR[ipz] -= rhoR * velocity[2];
  double newKER = 0.5*(consR[ipx]*consR[ipx] + consR[ipy]*consR[ipy] + consR[ipz]*consR[ipz])/rhoR;
  consR[iE] += newKER - oldKER;

  // First order reconstruction
  double rhoL1 = consL1[iRho];
  double oldKEL1 = 0.5*(consL1[ipx]*consL1[ipx] + consL1[ipy]*consL1[ipy] + consL1[ipz]*consL1[ipz])/rhoL1;
  consL1[ipx] -= rhoL1 * velocity[0];
  consL1[ipy] -= rhoL1 * velocity[1];
  consL1[ipz] -= rhoL1 * velocity[2];
  double newKEL1 = 0.5*(consL1[ipx]*consL1[ipx] + consL1[ipy]*consL1[ipy] + consL1[ipz]*consL1[ipz])/rhoL1;
  consL1[iE] += newKEL1 - oldKEL1;

  double rhoR1 = consR1[iRho];
  double oldKER1 = 0.5*(consR1[ipx]*consR1[ipx] + consR1[ipy]*consR1[ipy] + consR1[ipz]*consR1[ipz])/rhoR1;
  consR1[ipx] -= rhoR1 * velocity[0];
  consR1[ipy] -= rhoR1 * velocity[1];
  consR1[ipz] -= rhoR1 * velocity[2];
  double newKER1 = 0.5*(consR1[ipx]*consR1[ipx] + consR1[ipy]*consR1[ipy] + consR1[ipz]*consR1[ipz])/rhoR1;
  consR1[iE] += newKER1 - oldKER1;

  computeRotatedVectors();

  gammaL = center->getGamma();
  gammaR = neighbor->getGamma();
  rotate(consL); rotate(consR);
  rotate(consL1); rotate(consR1);

  gammaL1 = gammaL; gammaR1 = gammaR;
  HydroUtils::calFlux(gammaL, consL, fluxL);
  HydroUtils::calFlux(gammaR, consR, fluxR);

  HydroUtils::calFlux(gammaL1, consL1, fluxL1);
  HydroUtils::calFlux(gammaR1, consR1, fluxR1);

  if (verbose) {
    varVector primL, primR;
    double newGamma;
    HydroUtils::con2Prim( newGamma, consL, primL, false);
    HydroUtils::con2Prim( newGamma, consR, primR, false);
    printf("flux after trans: %d %d %5.3e (%5.3e %5.3e %5.3e) %5.3e %5.3e (%5.3e %5.3e %5.3e) %5.3e\n",
           center->iOrder, neighbor->iOrder, primL[0],primL[1], primL[2], primL[3], primL[4],primR[0],primR[1], primR[2], primR[3], primR[4]);
  }

}

void Face::solveRiemann(int activeRung, double volume, varVector& deltaQ, double &machNumber, bool predict)
{
  solveRiemann(activeRung, volume, deltaQ, machNumber, predict, NULL, false, false);
}

void Face::solveRiemann(int activeRung, double volume, varVector& deltaQ, double& machNumber, bool predict, bool debug)
{
  solveRiemann(activeRung, volume, deltaQ, machNumber, predict, NULL, debug, false);
}

void Face::solveRiemann(int activeRung, double volume, varVector& deltaQ, double &machNumber, bool predict, varVector *outDelta = NULL, bool debug=false, bool moreDiffusive=false)
{
  if (!validFace) return;
  if (!predict && !isActive()) return;
  double dt = HydroUtils::getTimeStep( iRung-1); // full time steps
  if( predict ) {
    fluxL = fluxL1; fluxR = fluxR1; consL = consL1; consR = consR1;
    dt = HydroUtils::getTimeStep( iRung); // half time steps
  }

  varVector fluxFace, uFace;

  int iRiemannSolver = HydroUtils::riemannSolver;
  iRiemannSolver = min(3, iRiemannSolver + (moreDiffusive ? 1: 0));
  if( iRiemannSolver == 0) {
    bool isOk;
    HydroUtils::HLLCFlux(gammaL, gammaR, fluxL, fluxR, consL, consR, fluxFace, uFace, machNumber, isOk, debug);
    if( !isOk) 
      HydroUtils::HLLFlux(gammaL1, gammaR1, fluxL1, fluxR1, consL1, consR1, fluxFace, uFace, machNumber, debug);
  }
  else if ( iRiemannSolver == 1)
    HydroUtils::HLLFlux(gammaL, gammaR, fluxL, fluxR, consL, consR, fluxFace, uFace, machNumber, debug);
  else if ( iRiemannSolver >= 2)
    HydroUtils::AthenaFlux(gammaL, gammaR, consL, consR, consL1, consR1, fluxFace, machNumber, debug, moreDiffusive);
  varVector delta = 0.;
  double v2 = dotProduct( velocity, velocity);

  // rotate back to lab frame
  reverse(fluxFace);

  // reset conservative variables.
  //consL = consLOrig; consR = consROrig; 
  //consL1 = consL1Orig; consR1 = consR1Orig;

  dVector Q2( fluxFace[ipx], fluxFace[ipy], fluxFace[ipz]);
  fluxFace[ipx] += fluxFace[iRho]*velocity[0];
  fluxFace[ipy] += fluxFace[iRho]*velocity[1];
  fluxFace[ipz] += fluxFace[iRho]*velocity[2];
  fluxFace[iE] += 0.5*fluxFace[iRho]*v2 + dotProduct( Q2, velocity);
#ifdef MHD
  HydroUtils::upwindAFlux(velocity, normal, consL, consR, fluxFace);
  dVector Bfield( 0.5*(consL[iBx]+consR[iBx]),0.5*(consL[iBy]+consR[iBy]),
                0.5*(consL[iBz]+consR[iBz])); 
  double BdotNorm = dotProduct(Bfield,normal);
  fluxFace[iBx] -= BdotNorm*velocity[0];
  fluxFace[iBy] -= BdotNorm*velocity[1];
  fluxFace[iBz] -= BdotNorm*velocity[2];
#endif
#ifdef RADIATION
  HydroUtils::upwindRadFlux( velocity, normal, consL, consR, fluxFace);
  if(debug && HydroUtils::instance->angleTrack >= 0) {
    int iAng = HydroUtils::instance->angleTrack + IRADVARSTART;
    //printf("conL: %5.3e conR: %5.3e\n", consL[iAng], consR[iAng]);
  }
#endif
  
  if( debug) {
    varVector newCons = consL - fluxFace*area/volume*dt, newPrims;
    double newGamma;
    HydroUtils::con2Prim( newGamma, newCons, newPrims, false);
    printf( "newPrims: %5.3e %5.3e %5.3e %5.3e %5.3e\n", newPrims[iRho], newPrims[ivx], newPrims[ivy], newPrims[ivz], newPrims[iIntE]);
    printf( "fluxFace: %5.3e %5.3e %5.3e %5.3e %5.3e\n", fluxFace[iRho], fluxFace[ivx], fluxFace[ivy], fluxFace[ivz], fluxFace[iIntE]);
    printf( "Area: %5.3e %5.3e %5.3e\n", area, volume, dt);
  }
  delta = fluxFace * (-area * dt);
  // make the corrections

  deltaQ = delta;
  if( outDelta != NULL) {
    *outDelta = delta;
  }

}

double HydroParticle::getPressure(bool fast) {
  double p = 0;
  if(fast) {
    double rho = primitives[iRho];
    p = primitives[iS]*pow(rho, gamma);
  } else {
    double cs, ie, T;
    HydroUtils::EOS( primitives, gamma, p, ie, cs, T);
  }
  return p;
}

void HydroParticle::setPrimitives(double rho, double vx, double vy, double vz,
                   double intEnergy) {
  primitives[iRho]  = rho;
  primitives[ivx]   = vx; primitives[ivy] = vy; primitives[ivz] = vz;
  primitives[iIntE] = intEnergy;
  finishPrimitives();
}

void HydroParticle::setPrimitives(int iType, double value) {
  if( iType >= 0 && iType < NTOTALVARS) 
    primitives[iType]  = value;
}

void HydroParticle::setPrimitives( varVector &prim) {
  primitives = prim;
  finishPrimitives();
}

void HydroParticle::finishPrimitives() {

  double p, cs, ie, T;
  // get gammaP
  if( std::isnan( primitives[iE])) {
      printf( "setPrimitives: isnan\n");
    }
  //check("finishPrimitives: before EOS");
  HydroUtils::EOS( primitives, gamma, p, ie, cs, T);
  primitives[iS] = p/pow(primitives[iRho],gamma);
  //check("finishPrimitives: after EOS");
  prim2Con();
  conservedPred = conserved;
}

void HydroParticle::computeGradients( std::unordered_map<int, HydroParticle*>& map,
                                      bool setPred, bool debug) {

  int numNeighbors           = hc->neighborIds.size();
  dVector xaxis(NDIM); xaxis = 0.; xaxis[0] = 1.;
  dVector yaxis(NDIM); yaxis = 0.; yaxis[1] = 1.;
  dVector zaxis(NDIM); zaxis = 0.; zaxis[2] = 1.;
  if( setPred) conservedPred = conserved;

  // for( int j = 0; j < NTOTALVARS; j++) {
  varVector phi_i = conservedPred;
  dVector   matrix[NDIM]; // = {{0.,0.,0.}, {0.,1.,0.}, {0.,0.,1.}}; //
                          // initialize to identity

  matrix[0] = xaxis; matrix[1] = yaxis; matrix[2] = zaxis;

  dVector sumA_ij[NTOTALVARS];

  for (int i = 0; i < NTOTALVARS; i++) {
    sumA_ij[i] = 0.;
  }
  varVector phi_min = conservedPred;
  varVector phi_max = conservedPred;
  double psi_ij[numNeighbors][NTOTALVARS];
  // compute the gradient

  for( int i = 0; i < numNeighbors; i++) {
    int nId = hc->neighborIds[i];
    std::unordered_map<int,HydroParticle*>::const_iterator elem = map.find(nId);

    if( elem == map.end())
      continue;
    HydroParticle *neighbor = map[nId]; //HydroParticleMap::getInstance() ->

    varVector phi_j = 0.;
    if( setPred )
      phi_j = neighbor -> conserved;
    else
      phi_j = neighbor -> conservedPred;

    // find max and min for slope limiter
    for(int j = 0; j < NTOTALVARS; j++) {
      phi_max[j] = max(phi_max[j],phi_j[j]);
      phi_min[j] = min(phi_min[j],phi_j[j]);
    }

    dVector A_ij = hc->faceNormals[i]*hc->faceAreas[i]/volume;
    dVector f_ij = hc->faceCentroids[i];
    dVector p2 = neighbor->getPosition();
    dVector p1 = getPosition();
    HydroUtils::applyBC( p1, p2);
    dVector s_j = neighbor -> getCentroid() + p2 -
        p1;
    dVector c_ij = f_ij - (centroid + s_j)*0.5;
    for(int l = 0; l < NDIM; l++)
      matrix[l] -= A_ij[l] * c_ij;

    for( int j = 0; j < NTOTALVARS; j++)
      sumA_ij[j] += A_ij * 0.5*(phi_i[j] + phi_j[j]);
  }

  HydroUtils::invertMatrix( matrix);

  for( int j = 0; j < NTOTALVARS; j++) {
    dVector gradient = 0.;
    for( int i = 0; i < NDIM; i++) {
      gradient[i] = (matrix[i]*sumA_ij[j]).sum();
      conservedGradients[j][i] = gradient[i];
    }
  }

  // compute the slope limiter
  double theta = 1;//HydroUtils::instance -> theta;
  varVector dPhiMin = (phi_min - conservedPred);
  varVector dPhiMax = (phi_max - conservedPred);
  for( int i = 0; i < numNeighbors; i++) {
    dVector f_ij = hc->faceCentroids[i]-getCentroid();
    for( int j=0; j<NTOTALVARS; j++) {
      double DeltaPhi = dotProduct(conservedGradients[j],f_ij);
      psi_ij[i][j] = 1.;
      if( DeltaPhi < 0.) {
        psi_ij[i][j] = theta*dPhiMin[j]/DeltaPhi;
      }
      else if( DeltaPhi > 0.) {
        psi_ij[i][j] = theta*dPhiMax[j]/DeltaPhi;
      }
    }
  }

  double alpha[NTOTALVARS];
  for( int j=0; j<NTOTALVARS; j++) {
    alpha[j] = 1.;
    for( int i = 0; i < numNeighbors; i++)
       alpha[j] = min(alpha[j], psi_ij[i][j]);
  }

#ifdef MHD
  alpha[iAx] = 1.; alpha[iAy] = 1.; alpha[iAz] = 1.; alpha[iPsi] = 1.;
  alpha[iBx] = 1.; alpha[iBy] = 1.; alpha[iBz] = 1.;
#endif

#ifdef RADIATION
for( int iFreq = 0; iFreq < NFREQVARS; iFreq++) { 
  for( int i = 0; i < NANGVARS; i++)
    //alpha[iFreq*NANGVARS+i+IRADVARSTART] = 0.;
    1;
}
#endif
  double alphaTot = 1.;
  for( int j = 0; j < NHYDROVARS; j++) {
    alphaTot = min(alphaTot, alpha[j]);
  }

  for( int j = 0; j < NHYDROVARS; j++) {
    alpha[j] = alphaTot;
  }  

  for( int j = 0; j < NTOTALVARS; j++)
    //for( int i = 0; i < NDIM; i++)
      conservedGradients[j] *= alpha[j];

}

void HydroParticle::setFacesAndNeighbors(voronoicell_neighbor& cell, bool completeCell) {
  if (!hc) hc = new HydroCell();
  else {
    hc->neighborIds.clear(); hc->faceAreas.clear(); hc->faceNormals.clear();
    hc->faceCentroids.clear(); hc->faceDistance.clear();
  }

  std::vector<double> fNormals;
  std::vector<double> fCentroids;

  cell.neighbors(hc->neighborIds);
  validCell = completeCell;
  maxSearchRadius= sqrt(cell.max_radius_squared());

  iNumNeighbors = hc->neighborIds.size();
  for( int i = 0; i < iNumNeighbors; i++)
    neighbors[i] = hc->neighborIds[i];

  cell.face_areas(hc->faceAreas);

  cell.normals(fNormals);
  cell.face_centroids(fCentroids);
  double cx, cy, cz;
  cell.centroid(cx, cy, cz);
  centroid[0] = cx;   centroid[1] = cy;   centroid[2] = cz;

  int numFaces = hc->faceAreas.size();

  for (int i = 0; i < numFaces; i++) {
    dVector n;
    n[0] = fNormals[3 * i]; n[1] = fNormals[3 * i + 1];
    n[2] = fNormals[3 * i + 2];
    hc->faceNormals.push_back(n);
  }

  int walker = 0;

  for (int i = 0; i < numFaces; i++) {
    dVector fCentroid;
    fCentroid[0] = fCentroids[3 * i]; fCentroid[1] = fCentroids[3 * i + 1];
    fCentroid[2] = fCentroids[3 * i + 2];
    double distanceToPlane = dotProduct(fCentroid, hc->faceNormals[i]);

    hc->faceDistance.push_back(distanceToPlane);
    hc->faceCentroids.push_back(fCentroid);
  }

  volume = cell.volume();

}


bool HydroParticle::check( const char *placeString) {
  if(iOrder == HydroUtils::instance->particleTrack || primitives[iIntE] > 1e30  || abs(primitives[ivz]) > 1e17) {
    printf("track Q at %s = %d %5.3e (%5.3e %5.3e %5.3e) %5.3e %5.3e %5.3e %7.5e\n", placeString, iOrder, primitives[iRho],
                      primitives[ivx], primitives[ivy], primitives[ivz], primitives[iIntE],
                      primitives[iS], gamma, length(particleAcceleration));
    printf("particle at = %d %5.3e %5.3e %5.3e\n",  iOrder, particlePosition[0], particlePosition[1], particlePosition[2]);
  }
  bool isNaNs = false;
  if(HydroUtils::instance->checkNaN) {
    for( int nvar=0; nvar < NTOTALVARS; nvar++) {
      isNaNs = isNaNs || std::isnan(deltaQ[nvar]) || std::isnan(conserved[nvar]);
    }
    if( isNaNs) {
      printf("%s isNaN = %d ", placeString, iOrder);
      for( int nvar=0; nvar < NTOTALVARS; nvar++) {
        if( std::isnan(deltaQ[nvar]) || std::isnan(conserved[nvar])) {
          printf(" %d ", nvar);
        }
      }
      printf("\n");
    }
  }
#ifdef COOLING_AD
  double checkGamma, p, ie, cs, T;
  HydroUtils::EOS( primitives, checkGamma, p, ie, cs, T); // check Temperature
#endif

  if( hc == 0 ) return false;
  return true;
}

void HydroParticle::setTimeStep( std::vector<HydroParticle*> &nearestNeighbors, bool checkNearest) {

  double eff_r = pow(volume / (4. * 3.14151 / 3.), 0.3333); // for 3-d
  double vx    = particleVelocity[0] - primitives[ivx];
  double vy    = particleVelocity[1] - primitives[ivy];
  double vz    = particleVelocity[2] - primitives[ivz];
  double cs    = getSoundSpeed();
  const double minTimeStep = HydroUtils::instance->minTimeStep; 
  double vA = 0., cspeed = 0.;
  double dtSignal = 1e50, dteffCorr = 1.;
#ifdef MHD
  vA = getAlfvenSpeed();
#endif
#ifdef RADIATION
  cspeed = HydroUtils::instance->redCInCodeUnits;
#endif
  double dteff = eff_r / (sqrt(vx * vx + vy * vy + vz * vz) + cs + vA + cspeed);
  if( abs(dteff) < minTimeStep){ 
    printf("Small Time Step dteff: %d %5.3e %5.3e %5.3e v %5.3e vA %5.3e rho %5.3e \n", iOrder, dteff, eff_r, 
		cs, sqrt(vx*vx + vy*vy + vz*vz), vA, primitives[iRho]);
    printf("position: %5.3e %5.3e %5.3e\n", particlePosition[0], particlePosition[1], particlePosition[2]);

  }
  timeStep = dteff;
  numSearch = nearestNeighbors.size();
  double minNeighborTimeStep = 1e30;
  for( int i = 0; i < nearestNeighbors.size(); i++) {
    HydroParticle *hp = nearestNeighbors[i];
    if( checkNearest) 
      minNeighborTimeStep = min( minNeighborTimeStep, hp->getTimeStep());
    
    dVector p2 = hp -> getPosition();
    dVector p1 = getPosition();

    HydroUtils::applyBC( p1, p2);
    dVector dR = p2 - p1;
    dVector nfVelocity(hp->primitives[ivx], hp->primitives[ivy], hp->primitives[ivz]); //neighbor fluid velocity
    dVector myfVelocity(primitives[ivx], primitives[ivy], primitives[ivz]); //neighbor fluid velocity
    dVector dV = nfVelocity-myfVelocity;
    double lengthdR = length(dR);
    double nCs = hp->getSoundSpeed();
    double vA = 0., nvA = 0.;
#ifdef MHD
    vA = getAlfvenSpeed(); nvA = hp->getAlfvenSpeed();
    double vPsi = HydroUtils::instance->vPsi;
    dtSignal = lengthdR/max(nCs + cs + 2.*cspeed + 2.*vPsi - dotProduct(dV, dR)/lengthdR + nvA + vA, TINY_NUMBER);
#else
    dtSignal = lengthdR/max(nCs + cs + 2.*cspeed - dotProduct(dV, dR)/lengthdR, TINY_NUMBER);
#endif
    if( lengthdR < 2.*eff_r) {
      dteffCorr = min(lengthdR/(2.*eff_r), dteffCorr);
    }
    timeStep = min(min(timeStep, abs(dtSignal)), dteffCorr*dteff);
    if( iOrder == HydroUtils::instance->particleTrack && dtSignal < minTimeStep) 
      printf("signalStep: %d %d %.3e, %.3e %.3e %.3e %.3e\n", iOrder, hp->iOrder, timeStep, 
            dtSignal, dteffCorr, dotProduct(dV, dR)/lengthdR, nCs);

#ifdef MHD
    if( abs(dtSignal) < minTimeStep){ 
      printf("Small Time Step: %d %5.3e %5.3e %5.3e %5.3e %5.3e vA: %5.3e %5.3e\n", iOrder, timeStep, lengthdR, nCs, cs, dotProduct(dV,dR)/lengthdR, nvA, vA);
    }
#endif
  }
  if( checkNearest && minNeighborTimeStep > 0.) 
    timeStep = min( timeStep, 1.414*minNeighborTimeStep);
  if( iOrder == HydroUtils::instance->particleTrack) {
    printf("setTimeStep: %d %5.3e %5.3e %5.3e %5.3e %.3e, %.3e, %.3e\n", iOrder, cs, timeStep, dteff, eff_r, dtSignal, dteffCorr,minNeighborTimeStep);
    fflush(stdout);
  }
}

void HydroParticle::prim2Con()
{
  HydroUtils::prim2Con(gamma, primitives, conserved);
}

void HydroParticle::con2Prim()
{
  HydroUtils::con2Prim(gamma, conserved, primitives);
}
void HydroParticle::update() {
#ifdef MHD
  // different Bfield update
  deltaQ[iBx] = 0.; deltaQ[iBy] = 0.; deltaQ[iBz] = 0.;
#endif
  Q           += deltaQ;
  deltaQ       = 0.;
  conserved    = Q / volume;
  conservedPred = conserved;
  //particleMass = Q[iRho];
  con2Prim();
  //check("before EOSClean");
  bool clean = HydroUtils::EOSclean( primitives, gamma);
  //check("after EOSClean");
  HydroUtils::EOSEnergyFix( machNumber, gamma, primitives);
#ifdef RADIATION
  bool radClean = HydroUtils::RadiationClean( primitives);
  if(radClean) {
    printf("RadClean involved at %5.3e %5.3e %5.3e\n", 
      particlePosition[0], particlePosition[1], particlePosition[2]);
  }
#endif
  //check("after EOSEnergyFix");
  // relaxation if necessary
#ifdef DAMPING
  primitives[ivx] *= relaxV; primitives[ivy] *= relaxV; primitives[ivz] *= relaxV;
#endif

  prim2Con();
  //check("after prim2Con");
  setupQ();
  //check("after setupQ");
  //reset();
}

bool HydroParticle::isNeighbor( int nId) {

  if( iNumNeighbors == 0 || !validCell) return true; // assume yes if no idea

  for(int i = 0; i < iNumNeighbors; i++) {
    if( neighbors[i] == nId) return true;
  }

  return false;
}

bool HydroParticle::checkSolution( std::unordered_map<int, HydroParticle*>& map) {
  bool noChange = true;
  for( int i = 0; i < hc->neighborIds.size(); i++) {
     // first get the neighbor
     int nId = hc->neighborIds[i];

     std::unordered_map<int,HydroParticle*>::const_iterator elem = map.find (nId);

     if( elem == map.end()){
       if( nId > 0) {
         printf( "This should never happen %d %d !!\n", iOrder, nId);
         noChange = false;
       }
       continue;
     }

     HydroParticle *neighbor = map[nId]; //HydroParticleMap::getInstance() ->

  }
  return noChange;
}

bool HydroParticle::getChange( int nId, double* change) {
/*  for(int i = 0; i < iNumNeighbors; i++) {
    if( neighbors[i] == nId) {
      for(int j = 0; j < DATAMAX; j++)
        change[j] = changes[i][j];
      return true;
    }
  }

  return false;
*/
  return true;
}


void HydroParticle::solveRiemann(int activeRung, std::unordered_map<int, HydroParticle*>& map, bool predictStep)
{
  double dEgrav = 0.;
  //check();
  if( predictStep && activeRung > getRung()) return;

  bool debug = iOrder == HydroUtils::instance->particleTrack;
  varVector deltaQ1 = 0.;
  for( int i = 0; i < hc->faceAreas.size(); i++) {
     //for(int j = 0; j < DATAMAX; j++)
       //changes[i][j] = 0.;
     // first get the neighbor
     int nId = hc->neighborIds[i];
     std::unordered_map<int,HydroParticle*>::const_iterator elem = map.find (nId);

     if( elem == map.end()){
       if( nId > 0)
         printf( "This should never happen %d %d !!\n", iOrder, nId);
       continue;
     }

     HydroParticle *neighbor = map[nId]; //HydroParticleMap::getInstance() ->

     bool firstOrderReconstruction = predictStep; 
     Face face( activeRung, hc->faceCentroids[i], hc->faceNormals[i], hc->faceAreas[i],
                this, neighbor, false && debug, firstOrderReconstruction);
     if( !face.isValid() || !face.isActive()) continue;
     varVector delta_ij; delta_ij = 0.;
     double faceMachNumber; faceMachNumber = 0.;
     dVector n = hc->faceNormals[i];
     varVector deltaFace; deltaFace = 0.;

     face.solveRiemann( activeRung, volume, deltaFace, faceMachNumber, predictStep, &delta_ij, debug, false);
     // do some checks
     varVector newCons = conserved + deltaFace/volume;
     varVector newPrims;
     double newGamma;
     HydroUtils::con2Prim( newGamma, newCons, newPrims, false);
     if( newPrims[iIntE] > 1e30 || abs(newPrims[ivz]) > 1e17) { 
       if( newPrims[iIntE] > 1e30) 
         printf("excess ie: %d %d %5.3e %5.3e %5.3e\n", iOrder, neighbor->iOrder, primitives[iIntE], newPrims[iIntE], neighbor->primitives[iIntE]);
       if( abs(newPrims[ivz]) > 1e17)
         printf("excess vz: %d %d %5.3e %5.3e %5.3e\n", iOrder, neighbor->iOrder, primitives[ivz], newPrims[ivz], neighbor->primitives[ivz]);
       printf( "Error in solveRiemann: %d\n going to first order", iOrder);

       Face face2( activeRung, hc->faceCentroids[i], hc->faceNormals[i], hc->faceAreas[i],
                this, neighbor, true);
       face2.solveRiemann( activeRung, volume, deltaFace, faceMachNumber, predictStep, &delta_ij, true, false);
     }
     deltaQ1 += deltaFace;
     
     machNumber = max(machNumber, faceMachNumber);
     
     // The change relative in sign compared to Springel is that delta_ij is defined opposite what Springel defines
     dVector avgAcceleration = (neighbor->particleAcceleration + particleAcceleration)*0.5;

     dVector pL             = getPosition();
     dVector pR             = neighbor->getPosition();
     HydroUtils::applyBC( pL, pR);
     dVector rij            = pR-pL;
     dEgrav                 -= 0.5* delta_ij[iRho]*dotProduct( rij,avgAcceleration);
  }

  if( !HydroUtils::instance -> useOldGravMethod)  
    deltaQ1[iE] += dEgrav; 
  dEgrav = 0.;

  if( !predictStep && !isActiveHalf(activeRung)) {
    deltaQ += deltaQ1;
    return; // break don't consolidate until 
  }

  double dt = HydroUtils::getTimeStep( predictStep ? getHalfRung() : getRung()); 

  varVector dQ = deltaQ1;
  deltaQ1[ipx] += (deltaQ1[iRho]*0.5 + Q[iRho])*dt*particleAcceleration[0];
  deltaQ1[ipy] += (deltaQ1[iRho]*0.5 + Q[iRho])*dt*particleAcceleration[1];
  deltaQ1[ipz] += (deltaQ1[iRho]*0.5 + Q[iRho])*dt*particleAcceleration[2];

  varVector dQ2 = deltaQ1;
  dVector dp(dQ[1],dQ[2],dQ[3]);
  dVector dp2(dQ2[1],dQ2[2],dQ2[3]);
  dVector da = dp/(dt*(Q[iRho]+deltaQ1[iRho]));
  dVector da2 = dp2/(dt*(Q[iRho]+deltaQ1[iRho]));
  dEgrav += Q[iRho]*dt*dotProduct(particleVelocity, particleAcceleration);
  double KE1=0.5*(Q[ipx]*Q[ipx]+Q[ipy]*Q[ipy]+Q[ipz]*Q[ipz])/Q[iRho];
  varVector Q2 = Q + deltaQ1;
  double KE2=0.5*(Q2[ipx]*Q2[ipx]+Q2[ipy]*Q2[ipy]+Q2[ipz]*Q2[ipz])/Q2[iRho];
  double deltaKE = KE2-KE1;
  dVector p1(Q[ipx],Q[ipy],Q[ipz]);
  dVector p2(Q2[ipx],Q2[ipy],Q2[ipz]);

  double dEgrav2 = dotProduct(p1,particleAcceleration) + dotProduct(p2,particleAcceleration);
  dEgrav2 *= 0.5*dt;
  if( debug || std::isnan(deltaQ1[iE]) || std::isnan(HydroUtils::instance->useOldGravMethod ? dEgrav2 : dEgrav))  {
    printf("Error: %d Q=%5.3e %5.3e %5.3e %5.3e %5.3e %7.5e\n", iOrder, Q[iRho],deltaQ1[ipx], deltaQ1[ipy],deltaQ1[ipz], deltaQ1[iE], dt);
    printf("Error: %d da=%5.3e %5.3e %5.3e %5.3e\n", iOrder, da[0], da[1], da[2], length(da)-length(particleAcceleration));
  }
  
  deltaQ1[iE] += !HydroUtils::instance->useOldGravMethod ? dEgrav : dEgrav2;
  sources( dt, deltaQ1);

  if( !predictStep) {
    deltaQ += deltaQ1;
  }
  else { 
    conservedPred += deltaQ1/volume;

    HydroUtils::EOScleanCons( conservedPred, gamma);
  }
}

// include the sources terms e.g. energy and MHD terms
void HydroParticle::sources( double dt, varVector &deltas) {

  double edot = HydroUtils::udot( primitives, dt);

  deltas[iE] += Q[iRho]*edot*dt;

#ifdef MHD
  // Evolve the vector potential
  dVector Efield;
  //dVector BPred(conservedPred[iBx],conservedPred[iBy],conservedPred[iBz]);
  dVector BPred(conserved[iBx],conserved[iBy],conserved[iBz]);
  dVector vPredOverC(primitives[ivx],primitives[ipy],primitives[ipz]);
  //dVector vPredOverC(conservedPred[ipx],conservedPred[ipy],conservedPred[ipz]);
  //vPredOverC /= max(conservedPred[iRho], HydroUtils::instance->smallRho);//*HydroUtils::instance->cInCodeUnits;

  crossProduct(vPredOverC, BPred, Efield);
  Efield *= -1.;
  deltas[iAx] += -(Efield[0] - conservedGradients[iPsi][0])*volume*dt;
  deltas[iAy] += -(Efield[1] - conservedGradients[iPsi][1])*volume*dt;
  deltas[iAz] += -(Efield[2] - conservedGradients[iPsi][2])*volume*dt;
  double gradDotA = conservedGradients[iAx][0] + conservedGradients[iAy][1] + conservedGradients[iAz][2];
  deltas[iPsi] += HydroUtils::instance->vPsi*gradDotA*volume*dt;
#endif // MHD

  // include user sources
  for( int iSource = 0; iSource < UserSources::numSources; iSource++)
    UserSources::instances[iSource] -> applySource( dt, particlePosition, volume, conserved, deltas);
}


void HydroParticle::setRegularizedVelocity( double dt, std::unordered_map<int, HydroParticle*>& map)
{
  double eff_r = pow(volume / (4. * 3.14151 / 3.), 0.3333); // for 3-d

  regularizedVelocity = 0.;
  double avg_r = 0.;
  int numNeighbors = 0;

  for( int i = 0; i < hc->faceAreas.size(); i++) {
     // first get the neighbor
     int nId = hc->neighborIds[i];
     std::unordered_map<int,HydroParticle*>::const_iterator elem = map.find (nId);

     if( elem == map.end())
       continue;
     avg_r += hc->faceDistance[i];
     numNeighbors++;
  }
  avg_r /= numNeighbors;

  eff_r = avg_r;
  for( int i = 0; i < hc->faceAreas.size(); i++) {
     // first get the neighbor
     int nId = hc->neighborIds[i];
     std::unordered_map<int,HydroParticle*>::const_iterator elem = map.find (nId);

     if( elem == map.end())
       continue;
     HydroParticle *neighbor = map[nId]; //HydroParticleMap::getInstance() ->

     double nCs = neighbor->getSoundSpeed();
     double cs = getSoundSpeed();
     dVector dV = neighbor->getVelocity() - particleVelocity;
     dVector dR = hc->faceNormals[i];
     double lengthdR = hc->faceDistance[i];
     double signalSpeed = -dotProduct( dV, dR)/lengthdR;
     if( signalSpeed > 0.) {
       double rFactor = (lengthdR)/eff_r;

       double rFactor1 = 0.5, rFactor2 = 0.25;
       if( rFactor < rFactor1) {
         double alpha = 0.0;
         dVector acc = dR*alpha*signalSpeed*min((rFactor1 - rFactor)/(rFactor1 - rFactor2), 1.);
         regularizedVelocity += acc;
       }
     }
  }
}

void HydroParticle::setRelaxVelocity() {
  double dt = HydroUtils::getTimeStep( getRung());
  relaxV = 1.0;
#ifdef DAMPING
  double invRelaxTau = HydroUtils::instance -> invRelaxTau;
  if( invRelaxTau > 0.)
    relaxV = min(exp(-dt*invRelaxTau), 1.);
#endif
}

void HydroParticle::reconstruct(Face          *face,
//                                double         dt,
                                dVector      & origin,
                                HydroParticle *neighbor,
                                varVector    & faceCons,
                                varVector    & faceCons1,
                                bool         firstOrderReconstruction) {


  faceCons = conserved;
  faceCons1 = conserved;

  if( !firstOrderReconstruction) {
    faceCons = conservedPred;
    faceCons1 = conservedPred;
    dVector p1 = origin;
    dVector p2 = particlePosition;
    HydroUtils::applyBC( p1, p2);
    dVector fCentroid = face->getCentroid()+p1;

    fCentroid -= (p2 + centroid); // reset the position
                                                         // relative to local
                                                         // origin
    varVector phi_j = neighbor->conservedPred;
    varVector phi_i = conservedPred;
    double    theta = HydroUtils::instance->theta;
    varVector deltaCons = 0.; deltaCons = 0.;
    varVector alpha = 1.; alpha = 1.;
    bool isNaN = false;
    //return;

    for( int nvar = 0; nvar < NTOTALVARS; nvar++) {
      double DeltaPhi = dotProduct(conservedGradients[nvar],fCentroid);
      double psi_ij = abs(DeltaPhi) >  0. ?  std::max(theta*(phi_j[nvar] -
          phi_i[nvar])/DeltaPhi, 0.) : 1. ;
      alpha[nvar] = std::min(1., psi_ij);
#ifdef MHD
      alpha[iAx] = 1.; alpha[iAy] = 1.; alpha[iAz] = 1.; // don't use limiter
#endif
      deltaCons[nvar] = alpha[nvar]*DeltaPhi;

      isNaN = isNaN || std::isnan( DeltaPhi);
    }
    //if( !isNaN)
    faceCons += deltaCons;
  }

  // clean the rho

  HydroUtils::EOScleanCons( faceCons, gamma);

#ifdef RADIATION
  HydroUtils::RadiationCleanCons( faceCons);
#endif
  return;
}

void HydroParticle::predict( int activeRung) {
  // half time step
  if( activeRung > getRung()) return;
  double hdt = HydroUtils::getTimeStep( getHalfRung());
  dVector cGradients[NTOTALVARS];
  varVector deltaCons; deltaCons = 0.;
  conservedPred = conserved;
  //return;
  for( int nvar = 0; nvar < NTOTALVARS; nvar++) {
    cGradients[nvar] = conservedGradients[nvar];
  }

  dVector w = getVelocity();
  double divRhoV =  cGradients[ipx][0] + cGradients[ipy][1] + cGradients[ipz][2];
  deltaCons[iRho] = -hdt*(divRhoV - dotProduct(cGradients[iRho], w));

  dVector v(primitives[ivx], primitives[ivy], primitives[ivz]);
  dVector vprime = v-w;
  double rhoGamma = pow(conserved[iRho], gamma - 1.);
  dVector pGrad = cGradients[iS] + cGradients[iRho]/conserved[iRho]*conserved[iS]*(gamma-1.);
  pGrad *= rhoGamma;
  dVector rhov = conserved[iRho]*v;

  // The below has no effect if MHD is not turned on
  dVector BdotGradB(0.,0.,0.);
  double BdotGradvdotB = 0.;

#ifdef MHD
  dVector Bfield(conserved[iBx],conserved[iBy],conserved[iBz]);
  dVector pmagGrad;
  pmagGrad[0] = Bfield[0]*cGradients[iBx][0] + Bfield[1]*cGradients[iBy][0] + Bfield[2]*cGradients[iBz][0];
  pmagGrad[1] = Bfield[0]*cGradients[iBx][1] + Bfield[1]*cGradients[iBy][1] + Bfield[2]*cGradients[iBz][1];
  pmagGrad[2] = Bfield[0]*cGradients[iBx][2] + Bfield[1]*cGradients[iBy][2] + Bfield[2]*cGradients[iBz][2];

  pGrad += pmagGrad;

  BdotGradB[0] = dotProduct(Bfield,cGradients[iBx]);
  BdotGradB[1] = dotProduct(Bfield,cGradients[iBy]);
  BdotGradB[2] = dotProduct(Bfield,cGradients[iBz]);

  dVector vdotGradB(0.,0.,0.);
  vdotGradB[0] = dotProduct(v,cGradients[iBx]);
  vdotGradB[1] = dotProduct(v,cGradients[iBy]);
  vdotGradB[2] = dotProduct(v,cGradients[iBz]);

  dVector BdotGradv(0.,0.,0.);
  dVector Gradv[3];
  for( int i = 0; i < 3; i++)
    Gradv[i] = (cGradients[ipx+i] - v[i]*cGradients[iRho])/conserved[iRho];

  BdotGradv[0] = dotProduct(Bfield,Gradv[0]);
  BdotGradv[1] = dotProduct(Bfield,Gradv[1]);
  BdotGradv[2] = dotProduct(Bfield,Gradv[2]);

  double div_vB = (divRhoV - dotProduct(v,cGradients[iRho]))/conserved[iRho];

  deltaCons[iBx] = -hdt*(vdotGradB[0] - BdotGradv[0] + Bfield[0]*div_vB);
  deltaCons[iBy] = -hdt*(vdotGradB[1] - BdotGradv[1] + Bfield[1]*div_vB);
  deltaCons[iBz] = -hdt*(vdotGradB[2] - BdotGradv[2] + Bfield[2]*div_vB);

  dVector tempRHS = vdotGradB + BdotGradv;

  BdotGradvdotB = dotProduct(Bfield, tempRHS);

#endif // MHD

  double p = conserved[iS]*rhoGamma;
  double div_v = (divRhoV - dotProduct(v,cGradients[iRho]))/conserved[iRho];
  deltaCons[ipx] =  -hdt*dotProduct(cGradients[ipx], vprime) - hdt*(conserved[ipx])*div_v +
		    hdt*(conserved[iRho]*particleAcceleration[0] - pGrad[0] - BdotGradB[0]);

  deltaCons[ipy] =  -hdt*dotProduct(cGradients[ipy], vprime) - hdt*(conserved[ipy])*div_v +
		    hdt*(conserved[iRho]*particleAcceleration[1] - pGrad[1] - BdotGradB[1]);

  deltaCons[ipz] =  -hdt*dotProduct(cGradients[ipz], vprime) - hdt*(conserved[ipz])*div_v +
		    hdt*(conserved[iRho]*particleAcceleration[2] - pGrad[2] - BdotGradB[2]);

  deltaCons[iE]  =  -hdt*dotProduct(cGradients[iE], vprime) - hdt*conserved[iE]*div_v -
                     hdt*dotProduct(pGrad, v) - hdt*p*div_v + hdt*BdotGradvdotB
                     + hdt*dotProduct(rhov, particleAcceleration);
  deltaCons[iS]  =  -hdt*dotProduct(cGradients[iS], vprime) - hdt*conserved[iS]*div_v;


  for( int i = NHYDROVARS; i < NVARS; i++) {
    double rhoScalar = conserved[i];
    
    dVector gradRhoScalar = cGradients[i];
    deltaCons[i] = -hdt*( dotProduct( gradRhoScalar, vprime) + rhoScalar*div_v);
  }

#ifdef RADIATION
  const double cspeed = HydroUtils::instance->redCInCodeUnits;
  for( int iFreq = 0; iFreq < NFREQVARS; iFreq++) {
    for( int i = 0; i < NANGVARS; i++) {
      int ivar = iFreq*NANGVARS + i + IRADVARSTART;
      // no half-time predictions for radiation
      //double w = dotProduct( HydroUtils::instance->radNormals[i], cGradients[ivar]);
      //deltaCons[ivar] = -hdt*w*cspeed;//*conserved[ivar]; 
    }
  }
#endif

#ifdef RADIATION
  for( int nvar = 0; nvar < NTOTALVARS; nvar++) {
    if( iOrder == HydroUtils::instance->particleTrack && 
        nvar == (HydroUtils::instance->angleTrack+IRADVARSTART)){
      //printf("predict final: %5.3e %5.3e (%5.3e %5.3e %5.3e) %5.3e\n", conservedPred[nvar], 
      //        deltaCons[nvar], cGradients[nvar][0],cGradients[nvar][1],cGradients[nvar][2], conserved[nvar]);
    }
  }  
#endif

  bool isNaN = false;
  isNaN = isNaN || std::isnan( deltaCons[iRho]);
  isNaN = isNaN || std::isnan( deltaCons[ipx]);
  isNaN = isNaN || std::isnan( deltaCons[ipy]);
  isNaN = isNaN || std::isnan( deltaCons[ipz]);
  isNaN = isNaN || std::isnan( deltaCons[iE]);
  isNaN = isNaN || std::isnan( deltaCons[iS]);

  //if( !isNaN)
  // no prediction for anything
  conservedPred += deltaCons;

  // clean the rho
  //conservedPred[iRho] = max(conservedPred[iRho], HydroUtils::instance->smallRho);
  // clean results
  HydroUtils::EOScleanCons( conservedPred, gamma);
}

void HydroParticle::kick( double fKick) {
  dVector fluidVelocity;

  fluidVelocity[0] = primitives[ivx];
  fluidVelocity[1] = primitives[ivy];
  fluidVelocity[2] = primitives[ivz];
  updateParticleVelocity( fluidVelocity*fKick);
  applyVelocityCorrection( HydroUtils::instance->voroMesh);
  //particleVelocity = 0.;
}

dVector HydroParticle::getVelocityCorrection(bool debug)
{
  dVector vc;
  double  cs = getSoundSpeed();
  const double eta = 0.25, chi = 0.5;
  double r    = pow(volume / (4. * 3.14151 / 3.), 0.3333);
  double d    = length(centroid);
  double test = d / (eta * r);

  if (r < 1e-16) vc = 0.;
  else if (test < 0.9) vc = 0.;
  else if (0.9 <= test <1.1)
    vc = centroid * cs / d * (d - 0.9 * eta * r)/(0.2 * eta * r);
  else vc = centroid * cs / d;

  if(debug) {
    printf( "Correction: %d %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e\n", iOrder, centroid[0], centroid[1], centroid[2], d, r, cs);
  }
  vc *= chi;

  return vc;
}

void HydroParticle::applyVelocityCorrection( const double voroMesh) {
  dVector vc = getVelocityCorrection();

  particleVelocity += voroMesh*vc;
  //particleVelocity += regularizedVelocity;
}

bool HydroParticle::isActive( int activeRung) {
  return activeRung <= getRung();
}

bool HydroParticle::isActiveHalf( int activeRung) {
  return activeRung == getHalfRung();
}

#ifdef MHD

void HydroParticle::mapAtoB() {
  // Map vector potential A to B cell centered values
  conserved[iBx] = conservedGradients[iAz][1] - conservedGradients[iAy][2];
  conserved[iBy] = conservedGradients[iAx][2] - conservedGradients[iAz][0];
  conserved[iBz] = conservedGradients[iAy][0] - conservedGradients[iAx][1];

  conserved[iBx] += globalBfield[0];
  conserved[iBy] += globalBfield[1];
  conserved[iBz] += globalBfield[2];


  primitives[iBx] = conserved[iBx];
  primitives[iBy] = conserved[iBy];
  primitives[iBz] = conserved[iBz];

  // reset the volume integrated B-field
  Q[iBx] = conserved[iBx] * volume;   Q[iBy] = conserved[iBy] * volume;  Q[iBz] = conserved[iBz] * volume;
}
#endif


#ifndef COOLING_NONE
COOL* HydroUtils::cl = NULL;

void  HydroUtils::setupEOS( COOL *cool) {
  cl = cool;
}
#endif
double HydroUtils::defaultGamma = 1.6667;
double HydroUtils::defaultMu    = 1.;
int HydroUtils::riemannSolver = 0;

double HydroUtils::EOSgamma( const double cs2, const double rho, const double p) {
  const double gamma = cs2*rho/p;
  return gamma;
}

bool HydroUtils::EOS( varVector& prim, double &gamma, double &p, double &ie, double &cs, double &temp, bool useEntropy) {

  double rho = prim[iRho];
  double e = prim[iE];
  double s = prim[iS];
  ie = e;
  bool eosSuccess = true;
#if !defined(COOLING_NONE) && !defined(COOLING_AD)
  COOLPARTICLE cp;
  MapToCoolParticle( prim, cp, cl);
  double PoverRho;
  if(!useEntropy) {
    CoolCodePressureOnDensitySoundSpeed( cl, &cp, e, rho, gamma, gamma-1., &PoverRho, &cs);
    p = PoverRho*rho;
    gamma = EOSgamma(cs*cs, rho, p);
    ie = e;
#if defined(COOLING_GRACKLE) || defined(COOLING_HELM) || defined(COOLING_SROEOS) || defined(COOLING_MESA)
    temp = CoolCodeEnergyToTemperature( cl, &cp, ie, rho, 0.);
#else
    temp = CoolCodeEnergyToTemperature( cl, &cp, ie, 0.);
#endif
    if(gamma < 1e-3 )
      printf("gamma is zero %5.3e %5.3e %5.3e %5.3e\n", rho, e, PoverRho, cs);
  }
  else {
    ie = e;
    CoolCodePressureDensity( cl, &cp, p, rho, &ie, &cs, &gamma);
#if defined(COOLING_GRACKLE) || defined(COOLING_HELM) || defined(COOLING_SROEOS) || defined(COOLING_MESA)
    temp = CoolCodeEnergyToTemperature( cl, &cp, ie, rho, 0.);
#else
    temp = CoolCodeEnergyToTemperature( cl, &cp, ie, 0.);
#endif
  }
  if(std::isnan(gamma) )
     printf("gamma is nan %5.3e %5.3e %5.3e %5.3e\n", rho, e, PoverRho, cs);
  //gamma = PoverRho/e + 1.;
#ifdef COOLING_HELM
  eosSuccess = cp.eosSuccess != 0;
#endif
#else
  COOLPARTICLE cp;
  MapToCoolParticle( prim, cp, cl);
//PHIL: NOTE THAT THE BELOW IS ONLY TRUE FOR IDEAL GAS BUT IS APPLIED EVERYWHERE
  //gamma = HydroUtils::defaultGamma;
  if( useEntropy) {
    p = s*pow(rho,gamma);
    double ie = p/(gamma - 1.)/rho;
#if defined(COOLING_GRACKLE) || defined(COOLING_HELM) || defined(COOLING_SROEOS) || defined(COOLING_MESA)
    temp = CoolCodeEnergyToTemperature( cl, &cp, ie, rho, 0.);
#else
    temp = CoolCodeEnergyToTemperature( cl, &cp, ie, 0.);
#endif
  } 
  else {
    p = (gamma - 1.) * rho * e;
#if defined(COOLING_GRACKLE) || defined(COOLING_HELM) || defined(COOLING_SROEOS) || defined(COOLING_MESA)
    temp = CoolCodeEnergyToTemperature( cl, &cp, ie, rho, 0.);
#else
    temp = CoolCodeEnergyToTemperature( cl, &cp, ie, 0.);
#endif
  }
  cs = sqrt(gamma*p/rho);
#endif
  return eosSuccess;
}

void HydroUtils::EOSEnergyFix( double machNumber, double &gamma, varVector& prim) {

  double rho = prim[iRho];
  double e = prim[iE];
  double s = prim[iS];
  if( HydroUtils::instance->useEntropy == 1 &&
      machNumber < HydroUtils::instance->machThreshold) { // use entropy
    double p = s*pow(rho, gamma), cs, ie, temp;
//PHIL; AGAIN ONLY TRUE FOR IDEAL GAS
#if defined(COOLING_AD)
    prim[iE] = p/rho/(gamma - 1.);
#endif
    if( std::isnan( prim[iE])) {
      printf( "EOSEnergyFix: isnan\n");
    }
    HydroUtils::EOS( prim, gamma, p, ie, cs, temp, true); // update gamma
    prim[iE] = ie;
  }
  else {
    double p, cs, ie, temp;
    if( std::isnan( prim[iE])) {
      printf( "EOSEnergyFix: isnan\n");
    }
    HydroUtils::EOS( prim, gamma, p, ie, cs, temp, false); // update gamma
    prim[iS] = p/pow(rho,gamma);
  }
}

bool HydroUtils::EOSclean( varVector& prim, double &gamma) {

  bool eosClean = true;
#ifndef COOLING_NONE
  eosClean = false;
  COOLPARTICLE cp;
  MapToCoolParticle( prim, cp, cl);
  if( prim[iRho] < HydroUtils::instance->smallRho) {


    eosClean = true;
    prim[iRho] = HydroUtils::instance->smallRho;
  }
  double p, cs;
  if( CoolClean( cl, &cp, gamma, &prim[iRho], &prim[iE], &p, &cs) != 0) {
    eosClean = true;
    gamma = EOSgamma(cs*cs, prim[iRho], p);
    prim[iS] = p/pow(prim[iRho],gamma);
  }
#endif
  return eosClean;
}

bool HydroUtils::EOScleanCons( varVector &con, double &gamma) {

  varVector prim;
  HydroUtils::con2Prim( gamma, con, prim, false);

  if( EOSclean( prim, gamma)) {
    prim2Con(gamma, prim, con);
    return true;
  }
  return false;
}

double HydroUtils::udot( varVector& prim, double dt) {

  double rho = prim[iRho];
  double e = prim[iE];
  double eOld = e;
  double edot = 0.;
#ifndef COOLING_NONE
  COOLPARTICLE cp;
  MapToCoolParticle( prim, cp, cl);
  CoolIntegrateEnergyCode( cl, NULL, &cp, &e, 0., rho, 0., NULL, dt*cl->dSecUnit);
  edot = (e-eOld)/dt;
#endif
  return edot;
}

bool HydroUtils::con2Prim(double& gamma, varVector& con, varVector& prim, 
  bool updateGamma = true) {
  double rho   = max( con[iRho], HydroUtils::instance->smallRho);
  double rhovx = con[ipx];
  double rhovy = con[ipy];
  double rhovz = con[ipz];
  double ener  = con[iE];
  double rhos  = con[iS];
  bool success = true;
  prim[iRho] = rho;
  double vx = rhovx / rho, vy = rhovy / rho, vz = rhovz / rho;
  prim[ivx] = vx;
  prim[ivy] = vy;
  prim[ivz] = vz;
#ifdef MHD
  ener -= 0.5*(con[iBx]*con[iBx] + con[iBy]*con[iBy] + con[iBz]*con[iBz]);
#endif
  double e = (ener / rho - 0.5 * (vx * vx + vy * vy + vz * vz));
  prim[iIntE] = e;
  prim[iS] = rhos/rho;
  // scalars
  for( int i = NHYDROVARS; i < NVARS; i++) 
    prim[i] = con[i]/rho;

#ifdef MHD
  prim[iBx]=con[iBx]; prim[iBy]=con[iBy]; prim[iBz]=con[iBz];
  prim[iAx]=con[iAx]; prim[iAy]=con[iAy]; prim[iAz]=con[iAz];
#endif

#ifdef RADIATION
  for( int i = IRADVARSTART; i < IRADVAREND; i++)
    prim[i] = con[i];
#endif

  if(e < 0. || rho < 0.) {
    success = false;
    //printf("internal Energy con2Prim< 0. %5.3e %5.3e\n", e, gamma);
  }
  double p, cs, ie, T;
    if( std::isnan( prim[iE])) {
      printf( "con2Prim: isnan %5.3e %5.3e %5.3e %5.3e %5.3e\n", ener, rho, vx, vy, vz);
    }
  return success && (!updateGamma || HydroUtils::EOS( prim, gamma, p, ie, cs, T));
}

bool HydroUtils::prim2Con(double& gamma, varVector& prim, varVector& con) {
  double rho = max(prim[iRho], HydroUtils::instance->smallRho);
  double vx = prim[ivx], vy = prim[ivy], vz = prim[ivz];
  double e = prim[iIntE];
  double s = prim[iS];

  con[iRho] = rho;
  con[ipx]  = rho * vx;
  con[ipy]  = rho * vy;
  con[ipz]  = rho * vz;
  double ener = rho * e + 0.5 * rho * (vx * vx + vy * vy + vz * vz);
#ifdef MHD
  ener += 0.5*(prim[iBx]*prim[iBx]+prim[iBy]*prim[iBy]+prim[iBz]*prim[iBz]);
#endif
  con[iE] = ener;
  con[iS] = rho*s;

  // do scalars
  for( int i = NHYDROVARS; i < NVARS; i++)
    con[i] = rho*prim[i];

#ifdef MHD
  con[iBx]=prim[iBx]; con[iBy]=prim[iBy]; con[iBz]=prim[iBz];
  con[iAx]=prim[iAx]; con[iAy]=prim[iAy]; con[iAz]=prim[iAz];
#endif

#ifdef RADIATION
  for( int i = IRADVARSTART; i < IRADVAREND; i++)
    con[i] = prim[i];
#endif

  double p, cs, ie, T;
    if( std::isnan( prim[iE])) {
      printf( "prim2Con: isnan\n");
    }
  return HydroUtils::EOS( prim, gamma, p, ie, cs, T);

}

void HydroUtils::calFlux(double &gamma, varVector& con, varVector& flux) {

  double rho = con[iRho];
  double rhovx = con[ipx];
  double rhovy = con[ipy];
  double rhovz = con[ipz];
  double ener = con[iE];
  double vx = rhovx / rho, vy = rhovy / rho, vz = rhovz / rho;
#ifdef MHD
  ener -= 0.5*(con[iBx]*con[iBx] + con[iBy]*con[iBy] + con[iBz]*con[iBz]);
#endif
  double e = ener/rho - 0.5*(vx*vx  + vy*vy + vz*vz);
  double rhos = con[iS];
  double s = rhos/rho;

  varVector prim; prim = 0.;
  prim[iRho] = rho;
  prim[iE] = e;
#ifdef MHD
  prim[iBx] = con[iBx]; prim[iBy] = con[iBy]; prim[iBz] = con[iBz];
  prim[iAx] = con[iAx]; prim[iAy] = con[iAy]; prim[iAz] = con[iAz];
#endif
  // do scalars
  for( int i = NHYDROVARS; i < NVARS; i++) 
    prim[i] = con[i]/rho;

  
  double p = rho*e*(gamma-1);
  // below is correct, but slow with tabulated equations of state
  if( !HydroUtils::instance->useEOSSpeedup) {
    double ie, T, cs;
    HydroUtils::EOS( prim, gamma, p, ie, cs, T);
  }

#ifdef MHD
  dVector Bfield(prim[iBx],prim[iBy],prim[iBz]);
  flux[iRho] = rhovx;
  p += 0.5*dotProduct(Bfield,Bfield);
  double Bx = Bfield[0], By = Bfield[1], Bz = Bfield[2];
  flux[ipx]  = rhovx * vx + p - Bx*prim[iBx];
  flux[ipy]  = rhovy * vx - Bx*prim[iBy]; 
  flux[ipz]  = rhovz * vx - Bx*prim[iBz];
  dVector vel(vx,vy,vz);
  ener = con[iE];
  flux[iE]   = (ener + p) * vx - Bx*dotProduct(vel,Bfield);
  flux[iS]   = rhos * vx;
  flux[iBx]  = 0.;
  flux[iBy]  = By*vx - Bx*vy;
  flux[iBz]  = Bz*vx - Bx*vz;
#else
  flux[iRho] = rhovx;
  flux[ipx]  = rhovx * vx + p;
  flux[ipy]  = rhovy * vx;
  flux[ipz]  = rhovz * vx;
  flux[iE]   = (ener + p) * vx;
  flux[iS]   = rhos * vx;
#endif
  // do the scalars
  for( int i = NHYDROVARS; i < NVARS; i++) {
    flux[i]  = con[i] * vx;
  }
}

// calculate the kinetic energy
double HydroUtils::calKE( varVector &con) {
  return 0.5*(con[ipx]*con[ipx] + con[ipy]*con[ipy] + con[ipz]*con[ipz])/con[iRho];
}

#ifdef RADIATION
void HydroUtils::upwindRadFlux( dVector &vel, dVector &norm, varVector &conL, varVector &conR, 
                                varVector &fluxFace) {
  const double cspeed = HydroUtils::instance->redCInCodeUnits;

  for( int i = 0; i < NRADVARS; i++) { // go over radiation variables
    // include the moving face
    dVector effVel = HydroUtils::instance->radNormals[i]*cspeed - vel;
    double w = dotProduct( effVel, norm);
    int ivar = i+IRADVARSTART;
    if( w > 0.)
      fluxFace[ivar] = conL[ivar]*w; 
    else
      fluxFace[ivar] = conR[ivar]*w;
  }
}
#endif

#ifdef MHD
// calculate the magnetic energy
double HydroUtils::calMagEnergy( varVector &con) {
  return 0.5*(con[iBx]*con[iBx] + con[iBy]*con[iBy] + con[iBz]*con[iBz]);
}
// calculate tue upwind flux
void HydroUtils::upwindAFlux( dVector &vel, dVector &norm,
                              varVector& conL, varVector& conR,
                              varVector& fluxFace) {
  double w = dotProduct( vel, norm);
  if( w > 0.) {
    fluxFace[iAx] = -conR[iAx] * w;
    fluxFace[iAy] = -conR[iAy] * w;
    fluxFace[iAz] = -conR[iAz] * w;
  }
  else {
    fluxFace[iAx] = -conL[iAx] * w;
    fluxFace[iAy] = -conL[iAy] * w;
    fluxFace[iAz] = -conL[iAz] * w;
  }

}
#endif

void HydroUtils::HLLFlux(double gammaL,
                         double gammaR,
                         varVector& FL,
                         varVector& FR,
                         varVector& uL,
                         varVector& uR,
                         varVector& fluxFace,
                         varVector& uFace,
                         double &machNumber,
                         bool debug) {
  double    smallRho = HydroUtils::instance -> smallRho;
  double    rhoL = max(uL[iRho], smallRho), rhoR = max(uR[iRho], smallRho);
  double    vL = uL[ipx] / rhoL, vR = uR[ipx] / rhoR;
  //double    pL = uL[iS]*pow(rhoL,gammaL-1.), pR = uR[iS]*pow(rhoR,gammaR-1.);
  double    pL = FL[ipx] - FL[iRho] * vL, pR = FR[ipx] - FR[iRho] * vR;
  double    csL = sqrt(gammaL * pL / rhoL), csR = sqrt(gammaR * pR / rhoR);

  double sL = min(vL, vR) - max(csL, csR);
  double sR = max(vL, vR) + max(csL, csR);
  machNumber = 0.;
  if (sL >= 0.) {
    fluxFace = FL;
    uFace = uL;
    machNumber = fabs(sL/csL);
  }
  else if (sR <= 0.) {
    fluxFace = FR;
    uFace = uR;
    machNumber = fabs(sR/csR);
  }
  else if ((sL < 0.) && (sR >= 0.)) {
    fluxFace = (FL * sR - FR * sL + (uR - uL) * sL * sR) / (sR - sL);
    uFace = (uR*sR - uL*sL + FL - FR)/(sR - sL);
    machNumber = max(fabs(sL/csL), fabs(sR/csR));
  }

}


void HydroUtils::HLLCFlux( double gammaL, double gammaR, varVector& FL, varVector& FR,
                           varVector& uL, varVector& uR,
                           varVector& fluxFace,
                           varVector& uFace,
                           double &machNumber,
                           bool &isOk,
                           bool debug) {
  varVector flux; flux = 0.;
  double    smallRho = HydroUtils::instance -> smallRho;
  double    rhoL = uL[iRho], rhoR = uR[iRho];
  double    vL = uL[ipx] / rhoL, vR = uR[ipx] / rhoR;
  double    vyL = uL[ipy] / rhoL, vyR = uR[ipy] / rhoR;
  double    vzL = uL[ipz] / rhoL, vzR = uR[ipz] / rhoR;
  //double    pL = uL[iS]*pow(rhoL,gammaL-1.), pR = uR[iS]*pow(rhoR,gammaR-1.);
  double    EL = uL[iE], ER = uR[iE];
  double    iEL =  EL/rhoL - 0.5*(vL*vL + vyL*vyL + vzL*vzL);
  double    iER =  ER/rhoR - 0.5*(vR*vR + vyR*vyR + vzR*vzR);
  double    pL = FL[ipx] - FL[iRho] * vL, pR = FR[ipx] - FR[iRho] * vR;
  double    csL = sqrt(gammaL * pL / rhoL), csR = sqrt(gammaR * pR / rhoR);
  double rhobar = 0.5*(rhoL + rhoR), csbar = 0.5*(csL+csR);
  double pvrs = 0.5*(pL + pR) - 0.5*(vR-vL)*rhobar*csbar;

  isOk = true;
  //if( pvrs < 0) printf("[HLLC Solver] Negative contact pressure\n");
  if( pvrs < 0) isOk = false;

  double pstar = max(0., pvrs);
  double alphaL = pstar > pL ? sqrt(1. + (gammaL+1.)/(2.*gammaL)*(pstar/pL - 1.)) : 1.;
  double alphaR = pstar > pR ? sqrt(1. + (gammaR+1.)/(2.*gammaR)*(pstar/pR - 1.)) : 1.;
  double sL = vL - csL*alphaL, sR = vR + csR*alphaR;

  //double sL = min(vL-csL, vR-csR);
  //double sR = max(vL+csL, vR+csR);
  //double sL = min(vL, vR) - max(csL, csR);
  //double sR = max(vL, vR) + max(csL, csR);
  machNumber = 0.;
  machNumber = max( fabs(sL/csL), fabs(sR/csR));
  double sStar = (pR-pL + rhoL*vL*(sL-vL) - rhoR*vR*(sR-vR))/(rhoL*(sL-vL) - rhoR*(sR-vR));

  if( sL >= 0.) flux = FL;
  else if( sR <= 0) flux = FR;
  else if( sL <= 0. && 0. <= sStar ) {
    machNumber = abs(sStar)/csbar;
    double norm = (sL-vL)/(sL-sStar);
    varVector uStarL = uL*norm;
    uStarL[ipx] = norm*rhoL*sStar;
    uStarL[iE] = norm*uL[iE] + norm*rhoL*(sStar-vL)*(sStar + pL/(rhoL*(sL-vL)));
    flux = FL + (uStarL - uL)*sL;
  }
  else if ( sStar <= 0. && 0. <= sR) {
    machNumber = abs(sStar)/csbar;
    double norm = (sR-vR)/(sR-sStar);
    varVector uStarR = uR*norm;
    uStarR[ipx] = norm*rhoR*sStar;
    uStarR[iE] = norm*uR[iE] + norm*rhoR*(sStar-vR)*(sStar + pR/(rhoR*(sR-vR)));
    flux = FR + (uStarR - uR)*sR;
  }
  fluxFace = flux;
}

void HydroUtils::AthenaFlux( double gammaL, double gammaR,
                             varVector& uL, varVector& uR,
                             varVector& uL1, varVector& uR1,
                             varVector& fluxFace,
                             double &machNumber,
                             bool debug,
                             bool moreDiffusive) {
  Cons1DS AthUl, AthUr, AthUl1, AthUr1, AthFluxFace;
  double Bxi = 0., Bxi1 = 0.;

  AthUl.d = uL[iRho];
  AthUl.Mx = uL[ipx]; AthUl.My = uL[ipy]; AthUl.Mz = uL[ipz];
  AthUl.E = uL[iE];
  AthUl.Gamma = gammaL;

  AthUr.d = uR[iRho];
  AthUr.Mx = uR[ipx]; AthUr.My = uR[ipy]; AthUr.Mz = uR[ipz];
  AthUr.E = uR[iE];
  AthUr.Gamma = gammaR;

  AthUl1.d = uL1[iRho];
  AthUl1.Mx = uL1[ipx]; AthUl1.My = uL1[ipy]; AthUl1.Mz = uL1[ipz];
  AthUl1.E = uL1[iE];
  AthUl1.Gamma = gammaL;

  AthUr1.d = uR1[iRho];
  AthUr1.Mx = uR1[ipx]; AthUr1.My = uR1[ipy]; AthUr1.Mz = uR1[ipz];
  AthUr1.E = uR1[iE];
  AthUr1.Gamma = gammaR;
#ifdef MHD
  Bxi = 0.5*(uL[iBx] + uR[iBx]);
  AthUl.By = uL[iBy]; AthUl.Bz = uL[iBz];
  AthUr.By = uR[iBy]; AthUr.Bz = uR[iBz];
  
  Bxi1 = 0.5*(uL1[iBx] + uR1[iBx]);
  AthUl1.By = uL1[iBy]; AthUl1.Bz = uL1[iBz];
  AthUr1.By = uR1[iBy]; AthUr1.Bz = uR1[iBz];
#endif

  Prim1DS AthWl = Cons1D_to_Prim1D( &AthUl, &Bxi);
  Prim1DS AthWr = Cons1D_to_Prim1D( &AthUr, &Bxi);
  
  Prim1DS AthWl1 = Cons1D_to_Prim1D( &AthUl1, &Bxi1);
  Prim1DS AthWr1 = Cons1D_to_Prim1D( &AthUr1, &Bxi1);

  int isOk = 0;
  if( HydroUtils::riemannSolver == 2 && !debug) 
    isOk = fluxes(AthUl, AthUr, AthWl, AthWr, Bxi, &AthFluxFace, &machNumber);
  else if ( debug || HydroUtils::riemannSolver == 3) { 
    hllFluxes(AthUl, AthUr, AthWl, AthWr, Bxi, &AthFluxFace, &machNumber);
    isOk = 1;
  }

  if( isOk == 0)  
    hllFluxes(AthUl1, AthUr1, AthWl1, AthWr1, Bxi1, &AthFluxFace, &machNumber);

  fluxFace[iRho] = AthFluxFace.d;
  fluxFace[ipx] = AthFluxFace.Mx; fluxFace[ipy] = AthFluxFace.My; fluxFace[ipz] = AthFluxFace.Mz;
  fluxFace[iE] = AthFluxFace.E;

}

double HydroUtils::limiter(double r) {
  double phi;

  phi = max(0., min(1., r)); // minmod limiter
  // phi = max( 0., max(min(2.*r, 1.), min(2.,r))); //superbee limiter
  // phi = max( 0., min(min(2.*r, 0.5*(1.+r)), 2.)); // MC
  // phi = max( 0., min( min(2.*r, (2.+r)/3.), 2.)); // Koren
  // cout << r << endl;
  // phi = 1.;
  // phi = (r + abs(r))/(1. + abs(r)); //van leer
  // phi = min(phi, 2.);
  return phi;
}

void HydroUtils::invertMatrix(dVector matrix[]) {
  dVector n0(NDIM), n1(NDIM), n2(NDIM);
  dVector n0p(NDIM), n1p(NDIM), n2p(NDIM);

  for (int i = 0; i < NDIM; i++) {
    n0[i] = matrix[i][0];
    n1[i] = matrix[i][1];
    n2[i] = matrix[i][2];
  }

  crossProduct(n1, n2, n0p);
  crossProduct(n2, n0, n1p);
  crossProduct(n0, n1, n2p);

  matrix[0] = n0p;
  matrix[1] = n1p;
  matrix[2] = n2p;
}

HydroUtils* HydroUtils::instance = new HydroUtils();

void HydroUtils::setupBC( double x, double y, double z) {
  HydroUtils::instance->xPeriod = x;
  HydroUtils::instance->yPeriod = y;
  HydroUtils::instance->zPeriod = z;
}

void HydroUtils::applyBC( dVector& left, dVector &right) {
  double x1 = left[0];
  double y1 = left[1];
  double z1 = left[2];

  double x2 = right[0];
  double y2 = right[1];
  double z2 = right[2];

  double xp = HydroUtils::instance->xPeriod;
  double yp = HydroUtils::instance->yPeriod;
  double zp = HydroUtils::instance->zPeriod;

  if( x1 > 0. && abs(xp+x2-x1) < abs(x2-x1))
      right[0]+=xp;
  if( x1 < 0. && abs(x2-xp-x1) < abs(x2-x1))
      right[0]-=xp;
  if( y1 > 0. && abs(yp+y2-y1) < abs(y2-y1))
      right[1]+=yp;
  if( y1 < 0. && abs(y2-yp-y1) < abs(y2-y1))
      right[1]-=yp;
  if( z1 > 0. && abs(zp+z2-z1) < abs(z2-z1))
      right[2]+=zp;
  if( z1 < 0. && abs(z2-zp-z1) < abs(z2-z1))
      right[2]-=zp;
}

double HydroUtils::getTimeStep( int _iRung) {
  return HydroUtils::instance -> dtRung[_iRung];
}
#endif //VORONOI
