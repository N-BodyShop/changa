#ifndef EXTERNALFORCE_HINCLUDED
#define EXTERNALFORCE_HINCLUDED

#include "parameters.h"
/// @brief External force parameters and routines
class ExternalForce {
public:
    int bBodyForce;          ///< Constant acceleration
    double dBodyForceConst;
    int bPatch;              ///< Patch in a disk
    double dCentMass;        ///< Central mass
    double dCentRad;         ///< Radius of central star
    double dOrbDist;         ///< Distance of the patch from the center
    double dOrbFreq;         ///< Orbital frequency
    int bCentralBody;        ///< Mass at the origin
    double dEqRad;           ///< Equatorial radius of central body
    double dJ2;              ///< Oblateness coefficients of central body
    double dJ4;
    double dJ6;
    int bLogarithmicHalo;    ///< Type of halo is logarithmic
    int bDoGasDrag;       /* Apply gas drag force to planetesimals */
    double dRadInfGas;       /* Inflation factor for radius of particles, used for gas drag calculations */
    double dSigma0;       /* Gas surface density at 1 AU */
    int bConstGasProf;    /* Ignore gas surface density and use a constant volume density */
    double dConstGasRho;  /* Value to use for constant gas volume density */
    double dP;            /* Power law slope of gas surface density profile */
    double dQ;            /* Power law slope of gas temperature profile */
    double dT0;           /* Gas temperature at 1 AU */
    double dMu;           /* Mean molecular weight of gas */
    double dCD;           /* Gas drag force coefficient */

    void AddParams(PRM prm);
    void CheckParams(PRM prm, struct parameters &param);
    Vector3D<double> applyGravPotential(GravityParticle *p, int bKepStep) const;
    void applyGasDrag(GravityParticle *p) const;
    inline void pup(PUP::er &p);
    };

inline void ExternalForce::pup(PUP::er &p) {
     p | bBodyForce;
     p | dBodyForceConst;
     p | bPatch;
     p | dCentMass;
     p | dCentRad;
     p | dOrbDist;
     p | dOrbFreq;
     p | bCentralBody;
     p | dEqRad;
     p | dJ2;
     p | dJ4;
     p | dJ6;
     p | bLogarithmicHalo;
     p | bDoGasDrag;
     p | dSigma0;
     p | bConstGasProf;
     p | dConstGasRho;
     p | dRadInfGas;
     p | dP;
     p | dQ;
     p | dT0;
     p | dMu;
     p | dCD;
    }

#endif
