/*
 * Collision detection and response module
 */
#include "ParallelGravity.h"
#include "collision.h"

///
/// @brief initalize parameters for collision detection and handling
///

void Collision::AddParams(PRM prm)
{
    nSmoothCollision = 64;
    prmAddParam(prm, "nSmoothCollision", paramInt, &nSmoothCollision,
        sizeof(int), "nSmoothCollision", "<number of particles to do collision search over> = 64");
    }

void Collision::CheckParams(PRM prm, struct parameters &param)
{
    if (!prmSpecified(prm, "nSmoothCollision"))
        nSmoothCollision = param.nSmooth;
#ifndef COLLISION
    if (param.bCollision)
        CkAbort("COLLISION must be enabled to use collisions\n");
#endif
    }

void Main::doCollisions(double dTime, double dDelta)
{
    CkPrintf("Do collision detection\n");
    CollisionSmoothParams pCS(dTime, dDelta, param.collision);
    double dfBall2OverSoft2 = 4.0*param.dhMinOverSoft*param.dhMinOverSoft;
    treeProxy.startSmooth(&pCS, 0, param.collision->nSmoothCollision,
              dfBall2OverSoft2, CkCallbackResumeThread());
    }

void CollisionSmoothParams::initSmoothCache(GravityParticle *p1)
{

    }

void CollisionSmoothParams::combSmoothCache(GravityParticle *p1,
                              ExternalSmoothParticle *p2)
{

    }

void CollisionSmoothParams::fcnSmooth(GravityParticle *p, int nSmooth, pqSmoothNode *nList)
{
    }
