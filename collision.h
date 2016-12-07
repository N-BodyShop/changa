#ifndef COLLISION_HINCLUDED
#define COLLISION_HINCLUDED

#include <float.h>

/// @brief Used to pass information about collision partners between processes
class ColliderInfo {
public:
    Vector3D<double> position;
    Vector3D<double> velocity;
    Vector3D<double> acceleration;
    Vector3D<double> w;
    double mass;
    double dtCol;
    double radius;
    int iOrder;
    int iOrderCol;
    int rung;
    /* If true, this collider is flagged for deletion and will merge with
       its more massive partner */
    int bMergerDelete;
    ColliderInfo() {
        dtCol = DBL_MAX;
        iOrderCol = -1;
        bMergerDelete = 0;
        }
    void pup(PUP::er &p) {
        p | position;
        p | velocity;
        p | acceleration;
        p | w;
        p | mass;
        p | dtCol;
        p | radius;
        p | iOrder;
        p | iOrderCol;
        p | rung;
        p | bMergerDelete;
        }
    };

//#include "parameters.h"

/// @brief Collision parameters and routines
class Collision : public PUP::able {
public:
    int nSmoothCollision; /* number of particles to search for collisions over */
    int bWall;            /* particles will bounce off a wall in the z plane */
    int bAllowMergers;    /* allow particles to merge if they collide at a slow speed */
    int bPerfectAcc;      /* all collisions result in a merger */
    double dWallPos;      /* location of wall along z axis */
    double dEpsN, dEpsT;  /* normal and transverse coefficients of restitution */

    void AddParams(PRM prm);
    void CheckParams(PRM prm, struct parameters &param);
    void doCollision(GravityParticle* p, ColliderInfo &c, int bMerge);
    void checkMerger(ColliderInfo &c1, ColliderInfo &c2);
    double LastKickTime(int rung, double baseTime, double timeNow);
    void setMergerRung(GravityParticle *p, ColliderInfo &c, ColliderInfo &cMerge,
                              double baseStep, double timeNow);
    void doWallCollision(GravityParticle *p);
    void mergeCalc(double r, double m, Vector3D<double> pos,
                   Vector3D<double> vel, Vector3D<double> w,
                   Vector3D<double> acc, Vector3D<double> *posNew,
                   Vector3D<double> *velNew, Vector3D<double> *wNew,
                   Vector3D<double> *aNew, double *radNew, ColliderInfo &c);
    void bounceCalc(double r, double m, Vector3D<double> pos,
                    Vector3D<double> vel, Vector3D<double> w,
                    Vector3D<double> *velNew, Vector3D<double> *wNew,
                    ColliderInfo &c);
    Collision() {
        dEpsN = 1.0;//dEpsN = 0.8;
        dEpsT = 1.0;//dEpsT = 1.0;
        }
   
    PUPable_decl(Collision);
    Collision(CkMigrateMessage *m) : PUP::able(m) {}
    inline void pup(PUP::er &p);
    };

inline void Collision::pup(PUP::er &p) {
    p | nSmoothCollision;
    p | bWall;
    p | dWallPos;
    p | dEpsN;
    p | dEpsT;
    }

#include "smoothparams.h"

/**
 * SmoothParams class for detecting and responding to collisions between particles
 */

class CollisionSmoothParams : public SmoothParams
{
    int bWall;
    double dWallPos;
    double dTime, dDelta;
    Collision coll;
    virtual void fcnSmooth(GravityParticle *p, int nSmooth,
               pqSmoothNode *nList);
    virtual int isSmoothActive(GravityParticle *p) {}
    virtual void initSmoothParticle(GravityParticle *p);
    virtual void initTreeParticle(GravityParticle *p) {}
    virtual void postTreeParticle(GravityParticle *p) {}
    virtual void initSmoothCache(GravityParticle *p);
    virtual void combSmoothCache(GravityParticle *p1,
                 ExternalSmoothParticle *p2);
public:
    CollisionSmoothParams() {}
    CollisionSmoothParams(int _iType, int am, double _dTime, double _dDelta,
                  int _bWall, double _dWallPos, Collision *collision) :
        coll (*collision) {
        iType = _iType;
        activeRung = am;
        bUseBallMax = 0;
        dTime = _dTime;
        dDelta = _dDelta;
        bWall = _bWall;
        dWallPos = _dWallPos;
        }
    PUPable_decl(CollisionSmoothParams);
    CollisionSmoothParams(CkMigrateMessage *m) : SmoothParams(m) {}
    virtual void pup(PUP::er &p) {
        SmoothParams::pup(p);
        p | coll;
        p | dDelta;
        p | dTime;
        p | bWall;
        p | dWallPos;
    }
    };

#endif
