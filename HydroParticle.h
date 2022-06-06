#ifndef HYDROPARTICLE_HH
#define HYDROPARTICLE_HH

#include <voro++.hh>
using namespace voro;
#include "dVector.h"
#include <vector>
#include <unordered_map>
//#include "UserSources.h"

#define MAX_HYDRO_RUNGS 40

#ifdef MHD
#define NHYDROVARS 13
#else
#define NHYDROVARS 6
#endif

#include "HydroParticleCool.h"

#define NVARS (NHYDROVARS + NSCALARS)
#ifdef RADIATION
// number of angle bins
#define NMU_QUAD 1
#define NANG_QUAD NMU_QUAD*(NMU_QUAD + 1)/2
// BELOW DOES NOT CHANGE
#define QUADRANTS 8
#define NANGVARS (QUADRANTS*NANG_QUAD)
#define NFREQVARS 1
#define NRADVARS NFREQVARS*NANGVARS
#else
#define NRADVARS 0
#endif
#define IRADVARSTART NVARS
#define IRADVAREND (NVARS + NRADVARS)
#define NTOTALVARS (NVARS + NRADVARS)

#define NEIGHBORMAX 128
#define DATAMAX 2

#define iRho 0
#define ivx 1
#define ivy 2
#define ivz 3
#define iIntE 4
#define iS 5
#define ipx 1
#define ipy 2
#define ipz 3
#define iE 4

#ifdef MHD
#define iBx 6
#define iBy 7
#define iBz 8
#define iAx 9
#define iAy 10
#define iAz 11
#define iPsi 12
#endif

#define SPEED_OF_LIGHT 2.99792458e10 // in cgs
#define KPC_CGS 3.0857e21
#ifndef TESTCODE
#include <pup_stl.h>
#endif
#include <iostream>
#include "cooling.h"

class varVector {
  double data[NTOTALVARS];

public:

  varVector() {for (int i = 0; i < NTOTALVARS; i++) data[i]=0.;}

  varVector(int size) {for (int i = 0; i < NTOTALVARS; i++) data[i]=0.;}

  inline int size() {
    return NTOTALVARS;
  }

  inline varVector& operator=(const varVector& v) {
    for (int i = 0; i < NTOTALVARS; i++) data[i] = v.data[i];

    return *this;
  }

  inline varVector& operator=(const double& d) {
    for (int i = 0; i < NTOTALVARS; i++) data[i] = d;

    return *this;
  }

  inline varVector operator-(const varVector& v) {
    varVector a;

    for (int i = 0; i < NTOTALVARS; i++) a.data[i] = data[i] - v.data[i];
    return a;
  }

  inline varVector operator+(const varVector& v) {
    varVector a;

    for (int i = 0; i < NTOTALVARS; i++) a.data[i] = data[i] + v.data[i];
    return a;
  }

  inline varVector operator*(const varVector& v) {
    varVector a;

    for (int i = 0; i < NTOTALVARS; i++) a.data[i] = data[i] * v.data[i];
    return a;
  }

  inline double& operator[](const int& i) {
    return data[i];
  }

  inline void       resize(int size) {}

  inline varVector& operator+=(const varVector& v) {
    for (int i = 0; i < NTOTALVARS; i++) data[i] += v.data[i];

    return *this;
  }

  inline varVector operator*(const double d) {
    varVector a;

    for (int i = 0; i < NTOTALVARS; i++) a.data[i] = data[i] * d;
    return a;
  }

  inline varVector operator/(const double d) {
    varVector a;

    for (int i = 0; i < NTOTALVARS; i++) a.data[i] = data[i] / d;
    return a;
  }

#ifndef TESTCODE
  void pup(PUP::er& p) {
    for (int i = 0; i < NTOTALVARS; i++) p | data[i];
  }
#endif
};


void   crossProduct(dVector& a,
                    dVector& b,
                    dVector& c);
double dotProduct(dVector& a,
                  dVector& b);
double length(dVector& a);

class Face;
class HydroParticle;

class Face {
private:

  dVector normal, centroid;
  dVector n0, n1, n2;
  dVector n0p, n1p, n2p;
  dVector velocity;
  double  invDetA;
  dVector debugData1, debugData2;

  double area, gammaL, gammaR, gammaL1, gammaR1;
  varVector fluxL, fluxR, consL, consR;
  varVector fluxL1, fluxR1, consL1, consR1; //first order reconstruction
  varVector consLOrig, consROrig, consL1Orig, consR1Orig;
  varVector Q;
  bool verbose;
  bool validFace;
  bool isActiveFace;
  int iRung;

  const double tolerance     = 1e-20;
  const double normTolerance = 1e-5;

  void computeRotatedVectors();
  void rotate(varVector& v);
  void reverse(varVector& v);

public:

  Face(int            activeRung,
       dVector      & c,
       dVector      & n,
       double         a,
       HydroParticle *center,
       HydroParticle *neighbor) {
   // Face(dt, c, n, a, center, neighbor, false);
    Face(activeRung, c, n, a, center, neighbor, false, false);
  }

  Face(int            activeRung,
       dVector      & c,
       dVector      & n,
       double         a,
       HydroParticle *center,
       HydroParticle *neighbor,
       bool           debug,
       bool           firstOrderReconstruction);

  dVector getCentroid() {
    return centroid;
  }

  void solveRiemann(int        activeRung,
                    double     volume,
                    varVector& deltaQ,
                    double&    machNumber,
                    bool       predict);

  void solveRiemann(int        activeRung,
                    double     volume,
                    varVector& deltaQ,
                    double&    machNumber,
                    bool       predict,
                    bool       debug);

  void solveRiemann(int        activeRung,
                    double     volume,
                    varVector& deltaQ,
                    double&    machNumber,
                    bool       predict,
		    varVector* outDelta,
                    bool       debug,
                    bool       moreDiffusive);

  inline bool isValid() { return validFace;}
  bool isActive() { return isActiveFace;}

};

class HydroCell {
  friend class HydroParticle;

protected:

  std::vector<double>  faceAreas, faceDistance;
  std::vector<dVector> faceNormals, faceCentroids;
  std::vector<int>     neighborIds;
};

class HydroUtils {
public:
  static double    defaultGamma;
  static double    defaultMu;
  static int       riemannSolver;
  double           cInCodeUnits;
  double           redCInCodeUnits;
  double           smallRho;
  int              checkNaN;
  int              useEntropy;
  int              useOldGravMethod;
  int              useEOSSpeedup;
  bool             useVanLeer;
  double           machThreshold;
  double           xPeriod, yPeriod, zPeriod;
  int              firstOrder;
  int              particleTrack;
  double           theta;
  double           invRelaxTau;
  double           vPsi;
  double           voroMesh;
  double           minTimeStep;
  int              maxRungs;
  double           dtRung[MAX_HYDRO_RUNGS];
#ifdef RADIATION
  int              angleTrack;
  dVector          radNormals[NRADVARS];
  double           radWeights[NRADVARS];
  double           initRadIntensity;
  double           initRadMu;
#endif
  static HydroUtils *instance;
#ifndef COOLING_NONE
  static COOL*     cl;
  static void      setupEOS( COOL *cl);
#endif
  static bool      EOS( varVector& prim, double &gamma, double &p, double &ie, double &cs, double &temp, bool useEntropy=false);
  static double    udot( varVector& prim, double dt);
  static void      EOSEnergyFix( double machNumber, double &gamma, varVector& prim);
  static bool      EOSclean( varVector& prim, double &gamma);
  static bool      EOScleanCons( varVector& con, double &gamma);
  static double    EOSgamma( const double cs2, const double rho, const double p);
  static void      MapToCoolParticle( varVector &prim, COOLPARTICLE &cp, COOL *cl );
  static void      setupBC( double x, double y, double z);
  static void      applyBC( dVector& left, dVector &right);
  static bool      prim2Con(double &gamma,
			                      varVector& prim,
                            varVector& con);
  static bool      con2Prim(double &gamma,
			    varVector& con,
                            varVector& prim) { return con2Prim( gamma, con, prim, true);}
  static bool      con2Prim(double &gamma,
			    varVector& con,
                            varVector& prim,
                            bool updateGamma);
  static void      calFlux(double     &gamma,
                           varVector& con,
                           varVector& flux);
  static double    calKE( varVector &con);
#ifdef MHD
  static double    calMagEnergy( varVector &con);
  static void      upwindAFlux( dVector &vel, dVector &norm,
                                varVector& conL,
                                varVector& conR,
                                varVector& fluxFace);
#endif
#ifdef RADIATION
  static void      upwindRadFlux( dVector &vel, dVector &norm,
                                  varVector &conL,
                                  varVector &conR,
                                  varVector &fluxFace);
#endif
  static void HLLFlux(double     gammaL,
                           double     gammaR,
                           varVector& FL,
                           varVector& FR,
                           varVector& uL,
                           varVector& uR,
                           varVector& fluxFace,
                           varVector& uFace,
                           double &machNumber,
                           bool debug=false);

  static void HLLCFlux(double     gammaL,
                           double     gammaR,
                           varVector& FL,
                           varVector& FR,
                           varVector& uL,
                           varVector& uR,
                           varVector& fluxFace,
                           varVector& uFace,
                           double &machNumber,
                           bool &isOk,
                           bool debug=false);
  static void AthenaFlux(double     gammaL,
                           double     gammaR,
                           varVector& uL,
                           varVector& uR,
                           varVector& uL1,
                           varVector& uR1,
                           varVector& fluxFace,
                           double &machNumber,
                           bool debug=false,
                           bool moreDiffusive=false);
  // static dVector HLLCFlux( double gamma, varVector& FL, varVector& FR,
  // varVector& uL, varVector& uR);
  static double limiter(double r);
  static void invertMatrix(dVector[]);
  static double getTimeStep( int _iRung);

#ifdef RADIATION
  void        initRadNormals();
  void        gauleg( double x1, double x2, double *x, double *w, int n);
  void        legendre_equal_weight();
  static bool RadiationClean( varVector &prim); 
  static bool RadiationCleanCons( varVector &con); 
#endif

};


class HydroParticle {
protected:

  HydroCell *hc;

  double timeStep;

  int iRung;

  //double  particleMass;
  dVector particleVelocity;
  dVector particlePosition;
  dVector particleAcceleration;
  dVector regularizedVelocity;
  double relaxV;
  double volume;
  dVector centroid;
#ifdef MHD
  dVector globalBfield;
#endif
  // const static int iRho=0, ivx=1, ivy=2, ivz=3, iIntE=4;
  // const static int ipx=1, ipy=2, ipz=3, iE=4;

  double gamma;
  varVector primitives, conserved, deltaQ, Q, conservedPred;
  dVector   conservedGradients[NTOTALVARS];

  bool validCell;

  int neighbors[NEIGHBORMAX];
  int iNumNeighbors;

public:
  int iOrder;
  int numSearch;
  double machNumber;
  double maxSearchRadius;

  //double changes[NEIGHBORMAX][DATAMAX];

  HydroParticle() {
    gamma = HydroUtils::defaultGamma;
    hc    = 0;
    iRung = 1; 
  }

  HydroParticle(double x, double y, double z, double vxp, double vyp,
                double vzp) {
    HydroParticle(x, y, z, vxp, vyp, vzp, 0., 0., 0., 0., 0.);
  }

  HydroParticle(double x,
                double y,
                double z,
                double vxp,
                double vyp,
                double vzp,
                double rho,
                double vx,
                double vy,
                double vz,
                double intEnergy) {
    HydroParticle(x,y,z,vxp,vyp,vzp,rho,vx,vy,vz,intEnergy,HydroUtils::defaultGamma);
  }

  HydroParticle(double x,
                double y,
                double z,
                double vxp,
                double vyp,
                double vzp,
                double rho,
                double vx,
                double vy,
                double vz,
                double intEnergy,
		double gammain) {

    particlePosition[0] = x; particlePosition[1] = y; particlePosition[2] = z;

    //particleMass        = 0.;
    particleVelocity[0] = vxp; particleVelocity[1] = vyp;
    particleVelocity[2] = vzp;
    gamma               = gammain;
    setPrimitives(rho, vx, vy, vz, intEnergy);
    hc    = 0;
    iRung = 1; 
    machNumber = 0.;
    iNumNeighbors = 0;
    deltaQ = 0.;

#ifdef MHD
    globalBfield = 0.;
#endif

    maxSearchRadius = 0.;
    numSearch = 0;
  }

  ~HydroParticle() {
    if (hc) delete hc;
  }

#ifndef TESTCODE
  void pup(PUP::er& p)
  {
    // remember to pup your superclass if there is one
    p | volume;
    p | centroid;
    p | iRung;
    p | iOrder;
    p | timeStep;
    p | particlePosition;
    p | particleVelocity;
    p | particleAcceleration;
    p | gamma;
    p | primitives;
    p | conserved;
    p | conservedPred;
    p | Q;
    p | deltaQ;
    p | machNumber;
    for( int i = 0; i < NEIGHBORMAX; i++) {
      p | neighbors[i];
    //  for( int j = 0; j < DATAMAX; j++)
      //  p | changes[i][j];
    }
    p | iNumNeighbors;

    for( int i = 0; i < NTOTALVARS; i++) {
      p | conservedGradients[i];
    }
    p | validCell;
    p | maxSearchRadius;
    p | numSearch;
    p | relaxV;
#ifdef MHD
    p | globalBfield;
#endif //MHD
  }
#endif //TESTCODE
  inline HydroParticle& operator=(const HydroParticle& v) {

    hc       = 0;
    timeStep = v.timeStep;
    iRung    = v.iRung;

    iOrder = v.iOrder;
    particleVelocity = v.particleVelocity;
    particleAcceleration = v.particleAcceleration;
    particlePosition = v.particlePosition;

    gamma      = v.gamma;
    primitives = v.primitives;
    conserved  = v.conserved;
    conservedPred  = v.conservedPred;
    Q          = v.Q;
    deltaQ     = v.deltaQ;
    machNumber = v.machNumber;
    centroid   = v.centroid;
    volume     = v.volume;
    numSearch  = v.numSearch;
    for (int i = 0; i < NTOTALVARS;
         i++) conservedGradients[i] = v.conservedGradients[i];
    iNumNeighbors = v.iNumNeighbors;
    for( int i = 0; i < iNumNeighbors; i++) {
      neighbors[i] = v.neighbors[i];
    //  for( int j = 0; j < DATAMAX; j++)
     //   changes[i][j] = v.changes[i][j];
    }

    validCell = v.validCell;
    maxSearchRadius = v.maxSearchRadius;
#ifdef MHD
    globalBfield = v.globalBfield;
#endif
    return *this;
  }

  double getVolume() {
    return volume;
  }

  void setPrimitives(double rho, double vx, double vy, double vz,
                     double intEnergy);

  void setPrimitives(int iType, double value);

  void setPrimitives( varVector &prim);

  void finishPrimitives();

  void setupQ()
  {
    Q = conserved*volume;
  }

  void QtoCons(){
    conserved = Q/volume;
    con2Prim();
  }

  bool check(){ return check("");}
  bool check( const char *string);
  bool isNeighbor( int nId);
  bool isActive( int activeRung);
  bool isActiveHalf( int activeRung);
  bool getChange( int nId, double *change);
  bool checkSolution( std::unordered_map<int, HydroParticle*>& map);
  int  getNumSearch() {
    return numSearch;
  }

  int  getNumNeighbors() {
    return iNumNeighbors;
  }
  bool isValid() {
    return validCell;
  }

  void getQ( varVector &Qout) {
    Qout = Q;
  }

  void getDeltaQ( varVector &Qout) {
    Qout = deltaQ;
  }

  void getPrimitives(varVector& outVector)
  {
    outVector = primitives;
  }

  void getConserved(varVector& outVector)
  {
    outVector = conserved;
  }

  void getConservedPred(varVector& outVector)
  {
    outVector = conservedPred;
  }

  void getGradients(dVector *outVectors)
  {
    for( int i = 0; i < NTOTALVARS; i++) {
       outVectors[i] = conservedGradients[i];
    }
  }

  inline dVector& getVelocity() {
    return particleVelocity;
  }

  inline dVector& getPosition() {
    return particlePosition;
  }

  inline dVector& getParticleAcceleration() {
    return particleAcceleration;
  }

  inline double& getMass() { return Q[0];}

  inline dVector getCentroid() {
    return centroid;
  }

  // iRung accessors
  inline int getRung() { return iRung-1; } //iRung are guaranteed to be 1 or largers for Voronoi
  inline int getHalfRung() { return iRung; }
  void setRung( int _iRung) { iRung = _iRung; }

  dVector getVelocityCorrection(bool debug=false);

  void    applyVelocityCorrection( const double voroMesh);

  // update the particle position
  void updateParticlePosition(double x, double y, double z) {
    particlePosition[0] = x; particlePosition[1] = y; particlePosition[2] = z;
  }

  void updateParticleVelocity(dVector velocity) {
    particleVelocity = velocity;
  }

  void updateParticleVelocity(double vx, double vy, double vz) {
    particleVelocity[0] = vx; particleVelocity[1] = vy; particleVelocity[2] = vz;
  }

  void updateParticleAcceleration( double ax, double ay, double az) {
    particleAcceleration[0] = ax; particleAcceleration[1] = ay; particleAcceleration[2] = az;
  }

  void drift(double dt) {
    particlePosition += particleVelocity * dt;
  }

  void kick( double fKick = 1.);

  inline double getTimeStep() {
    return timeStep;
  }

  void output() {
    printf("cell: %lu %lu %lu %lu: ", hc->neighborIds.size(),
           hc->faceAreas.size(), hc->faceNormals.size(),
           hc->faceCentroids.size());

    for (int i = 0; i < hc->neighborIds.size(); i++) printf("%d, %5.3e ",
                                                            hc->neighborIds[i],
                                                            hc->faceDistance[i]);
    printf("\n");
  }

  void reset() {
    iNumNeighbors = 0;
    if (hc) delete hc;
    hc = 0;
    validCell = false;
  }

  double getGamma() {
    return gamma;
  }

  double getPressure(bool fast=true);

  double getSoundSpeed() {
    double rho = primitives[iRho];
    double p = getPressure(true);
    double cs = sqrt( gamma*p/rho);
    return cs;
    //return sqrt(gamma*(gamma-1.)*primitives[iIntE]);
  }

#ifdef MHD
  double getAlfvenSpeed() {
    double rho = primitives[iRho];
    double Bx = primitives[iBx], By = primitives[iBy], Bz = primitives[iBz];
    double B2 = Bx*Bx + By*By + Bz*Bz;
    double vAlfven = sqrt( B2/fmax(rho,HydroUtils::instance->smallRho));
    return vAlfven;
  }
#endif

  void   setFacesAndNeighbors(voronoicell_neighbor& cell, bool completeCell);
  void   setTimeStep() {
    std::vector<HydroParticle*> emptyList;
    setTimeStep(emptyList);
  }

  void   setTimeStep( std::vector<HydroParticle*> &nearestNeighbors, bool checkNearest=false);
  void   computeGradients( std::unordered_map<int, HydroParticle*>& map, bool bSetPred=true, bool debug=false);

  void   prim2Con();
  void   con2Prim();
  double getMinDt();

  void   solveRiemann(int activeRung, std::unordered_map<int, HydroParticle*>& map, bool predictStep);
  void   setRegularizedVelocity(double dt, std::unordered_map<int, HydroParticle*>& map);

  void   setRelaxVelocity();

  void   update();

  void   reconstruct(Face          *face,
//                     double         dt,
                     dVector      & origin,
                     HydroParticle *neighbor,
                     varVector    & faceCons,
                     varVector    & faceCons1,
                     bool         firstOrderConstruction = false);
  void   predict( int activeRung);

  void   sources( double dt, varVector &deltas);

#ifdef MHD
  void   mapAtoB();
  void   updateGlobalBfield( const double bx, const double by, const double bz) {
    globalBfield[0] = bx; globalBfield[1] = by; globalBfield[2] = bz;
  }
#endif

#ifdef RADIATION
  dVector getRadFlux();
  double  getRadEnergy(); 
#endif
};

#endif // ifndef HYDROPARTICLE_HH
