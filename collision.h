#ifndef COLLISION_HINCLUDED
#define COLLISION_HINCLUDED

#include "parameters.h"
#include <float.h>

class ColliderInfo {
public:
    Vector3D<double> position;
    Vector3D<double> velocity;
    Vector3D<double> w;
    double mass;
    double dtCol;
    double radius;
    int iOrder;
    ColliderInfo() {
        dtCol = DBL_MAX;
        }
    void pup(PUP::er &p) {
        p | position;
        p | velocity;
        p | w;
        p | mass;
        p | dtCol;
        p | radius;
        p | iOrder;
        }
    };

/// @brief Collision parameters and routines
class Collision : public PUP::able {
public:
    int nSmoothCollision; /* number of particles to search for collisions over */

    void AddParams(PRM prm);
    void CheckParams(PRM prm, struct parameters &param);
    void doCollision(GravityParticle* p, ColliderInfo &c);
    void bounce(GravityParticle* p, ColliderInfo &c, double dEpsN, double dEpsT);
    Collision() {}
   
    PUPable_decl(Collision);
    Collision(CkMigrateMessage *m) : PUP::able(m) {}
    inline void pup(PUP::er &p);
    };

inline void Collision::pup(PUP::er &p) {
    p | nSmoothCollision;
    }

#include "smoothparams.h"

/**
 * SmoothParams class for detecting and responding to collisions between particles
 */

class CollisionSmoothParams : public SmoothParams
{
    double dTime, dDelta;
    Collision coll;
    virtual void fcnSmooth(GravityParticle *p, int nSmooth,
               pqSmoothNode *nList);
    virtual int isSmoothActive(GravityParticle *p) {}
    virtual void initSmoothParticle(GravityParticle *p) {}
    virtual void initTreeParticle(GravityParticle *p) {}
    virtual void postTreeParticle(GravityParticle *p) {}
    virtual void initSmoothCache(GravityParticle *p);
    virtual void combSmoothCache(GravityParticle *p1,
                 ExternalSmoothParticle *p2);
public:
    CollisionSmoothParams() {}
    CollisionSmoothParams(int _iType, int am, double _dTime, double _dDelta, Collision *collision) :
        coll (*collision) {
        iType = _iType;
        activeRung = am;
        bUseBallMax = 0;
        dTime = _dTime;
        dDelta = _dDelta;
        }
    PUPable_decl(CollisionSmoothParams);
    CollisionSmoothParams(CkMigrateMessage *m) : SmoothParams(m) {}
    virtual void pup(PUP::er &p) {
        SmoothParams::pup(p);
        p | coll;
        p | dDelta;
        p | dTime;
    }
    };

#endif
