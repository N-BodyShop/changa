#ifndef EXTERNALGRAVITY_HINCLUDED
#define EXTERNALGRAVITY_HINCLUDED

#include "parameters.h"

/// @brief External gravity parameters and routines
class ExternalGravity {
public:
    int bBodyForce;          ///< Constant acceleration
    double dBodyForceConst;
    int bPatch;              ///< Patch in a disk
    double dCentMass;        ///< Central mass
    double dOrbDist;         ///< Distance of the patch from the center
    int bCentralBody;        ///< Mass at the origin
    double dEqRad;           ///< Equatorial radius of central body
    double dJ2;              ///< Oblateness coefficients of central body
    double dJ4;
    double dJ6;

    void AddParams(PRM prm);
    void CheckParams(PRM prm, struct parameters &param);
    Vector3D<double> applyPotential(GravityParticle *p) const;
    inline void pup(PUP::er &p);
    };

inline void ExternalGravity::pup(PUP::er &p) {
     p | bBodyForce;
     p | dBodyForceConst;
     p | bPatch;
     p | dCentMass;
     p | dOrbDist;
     p | bCentralBody;
     p | dEqRad;
     p | dJ2;
     p | dJ4;
     p | dJ6;
    }

#endif
