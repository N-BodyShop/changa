/*
 * Collision detection and response module
 */
#include <float.h>
#include "ParallelGravity.h"
#include "collision.h"
#include "smooth.h"
#include "Reductions.h"
#include "physconst.h"

///
/// @brief initalize parameters for collision detection and handling
///

void Collision::AddParams(PRM prm)
{
    nSmoothCollision = 64;
    prmAddParam(prm, "nSmoothCollision", paramInt, &nSmoothCollision,
        sizeof(int), "nSmoothCollision", "<number of particles to do collision\
                                           search over> = 64");
   
    bWall = 0;
    prmAddParam(prm, "bWall", paramBool, &bWall,
        sizeof(int), "bWall", "<Particles bounce off of a plane pointed in\
                                 the z direction> = 0");

    bAllowMergers = 0;
    prmAddParam(prm, "bAllowMergers", paramBool, &bAllowMergers,
        sizeof(int), "bAllowMergers", "<Particles merge on collision if their\
                                        approach speed is small enough> = 0");

    dWallPos = 0.0;
    prmAddParam(prm, "dWallPos", paramDouble, &dWallPos,
        sizeof(double), "dWallPos", "<Z coordinate of wall> = 0");

    bPerfectAcc = 0;
    prmAddParam(prm, "bPerfectAcc", paramBool, &bPerfectAcc,
        sizeof(int), "bPerfectAcc", "<All collisions result in a merger> = 0");

    iMinBinaryRung = 6;
    prmAddParam(prm, "iMinBinaryRung", paramInt, &iMinBinaryRung,
        sizeof(int), "iMinBinaryRung", "<Dont merge particles in binaries below this rung> = 0");

    dMaxBinaryEcc = 0.9;
    prmAddParam(prm, "dMaxBinaryEcc", paramDouble, &dMaxBinaryEcc,
        sizeof(double), "dMaxBinaryEcc", "<Maximum allowed eccentricty to merge bound particles> = 0.5");

    dBallFac = 2.0;
    prmAddParam(prm, "dBallFac", paramDouble, &dBallFac,
        sizeof(double), "dBallFac", "<Scale factor for collision search radius> = 2.0");

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

/**
 * @brief doCollisions is used to detect and resolve collisions between
 * particles. 
 *
 * A collision search is done using SmoothParams, where the  'dtCol' 
 * field of all particles which will undergo a collision in the time step is
 * set. The function then searches through all of the particles and resolves
 * the soonest collision using the physics specificed in the Collision class.
 * After this, another collision search is done and this process is repeated
 * until there are no more collisions during the current time step.
 *
 * @param dTime The current time in the simulation (in simulation units)
 * @param dDelta The size of the next timestep (in simulation units)
 */
void Main::doCollisions(double dTime, double dDelta)
{
    int bHasCollision;
    int nColl = 0;
    double startTime1, endTime1;
    double tSmooth = 0.0;
    double tRetrieve = 0.0;
    double tResolve = 0.0;
    
    startTime1 = CkWallTimer();
    do {
        bHasCollision = 0;

        CollisionSmoothParams pCS(TYPE_DARK, 0, dTime, dDelta, 
           param.collision.bWall, param.collision.dWallPos,
           param.collision.bAllowMergers, param.collision.dMaxBinaryEcc,
           param.collision.iMinBinaryRung, param.collision);
        double startTime = CkWallTimer();
        treeProxy.startReSmooth(&pCS, CkCallbackResumeThread());
        double endTime = CkWallTimer();
        tSmooth += endTime-startTime;

        startTime = CkWallTimer();
        CkReductionMsg *msgChk;
        treeProxy.getCollInfo(CkCallbackResumeThread((void*&)msgChk));
        ColliderInfo *c = (ColliderInfo *)msgChk->getData();

        // If only one particle participating in the collision detected the imminent collision
        // (due to different search radii), we need to go back and ask for the second collider
        if (c[0].dtCol <= dDelta && c[1].dtCol > dDelta) {
            treeProxy.getCollInfo(c[0].iOrderCol, CkCallbackResumeThread((void*&)msgChk));
            c[1] = *(ColliderInfo *)msgChk->getData();
            c[1].dtCol = c[0].dtCol;
            }
        else if (c[1].dtCol <= dDelta && c[0].dtCol > dDelta) {
            treeProxy.getCollInfo(c[1].iOrderCol, CkCallbackResumeThread((void*&)msgChk));
            c[0] = *(ColliderInfo *)msgChk->getData();
            c[0].dtCol = c[1].dtCol;
            }

        endTime = CkWallTimer();
        tRetrieve += endTime-startTime;

        if (c[0].dtCol <= dDelta) {
            startTime = CkWallTimer();
            bHasCollision = 1;
            // Binary merger
            if (c[0].dtCol < 0.) {
                nColl++;
                CkPrintf("Particles %d and %d stuck in binary...merging\n", c[0].iOrder, c[1].iOrder);
                treeProxy.resolveCollision(param.collision, c[0], c[1], 1,
                                           dDelta, dTime, CkCallbackResumeThread());
                }
            // Collision with wall
            else if (c[0].iOrderCol == -2 && param.collision.bWall) {
                nColl++;
                treeProxy.resolveWallCollision(param.collision, c[0], 
                                               CkCallbackResumeThread());
                }
            // Collision with particle
            else {
                if (c[0].iOrderCol != c[1].iOrder) {
                    CkAbort("Warning: Collider pair mismatch\n");
                    }
                nColl++;
                
                if (param.collision.bAllowMergers) {
                    param.collision.checkMerger(c[0], c[1]);
                    }
                if (c[0].bMergerDelete || c[1].bMergerDelete) {
                    treeProxy.resolveCollision(param.collision, c[0], c[1], 1,
                                               dDelta, dTime, CkCallbackResumeThread());
                    }
                else {
                    treeProxy.resolveCollision(param.collision, c[0], c[1], 0,
                                               dDelta, dTime, CkCallbackResumeThread());
                    }
                }
                endTime = CkWallTimer();
                tResolve += endTime-startTime;
            }

        delete msgChk;
        } while (bHasCollision);
        endTime1 = CkWallTimer();
        CkPrintf("Resolved %d collisions\n", nColl);
        CkPrintf("Collider search took %g seconds\n", tSmooth);
        CkPrintf("Collider retrieval took %g seconds\n", tRetrieve);
        CkPrintf("Collision resolution took %g seconds\n", tResolve);
        CkPrintf("doCollision took %g seconds\n", endTime1-startTime1);
        addDelParticles();

    }

/**
 * @brief Searches for the particle with the soonest dtCol on this TreePiece
 *
 * Contributes two ColliderInfo objects, one for each of the particles
 * participating in the collision. If only one of the particles resides on this
 * tree piece, then the second ColliderInfo object will have a dtCol of DBL_MAX
 * and none of the other fields set.
 */
void TreePiece::getCollInfo(const CkCallback& cb)
{
    double dtMin = DBL_MAX;
    ColliderInfo ci[2];
    ci[0].dtCol = dtMin;
    for (unsigned int i=1; i <= myNumParticles; i++) {
        GravityParticle *p = &myParticles[i];
        if (p->dtCol < dtMin) {
            dtMin = p->dtCol;
            ci[0].position = p->position;
            ci[0].velocity = p->velocity;
            ci[0].acceleration = p->treeAcceleration;
            ci[0].w = p->w;
            ci[0].mass = p->mass;
            ci[0].radius = p->soft*2.;
            ci[0].dtCol = p->dtCol;
            ci[0].iOrder = p->iOrder;
            ci[0].iOrderCol = p->iOrderCol;
            ci[0].rung = p->rung;
            }
        }

    // Check to see if the second collider is here
    int bFoundC1 = 0;
    ci[1].dtCol = DBL_MAX;
    for (unsigned int i=1; i <= myNumParticles; i++) {
        GravityParticle *p = &myParticles[i];
        if (p->dtCol == dtMin && p->dtCol < DBL_MAX) {
            if (bFoundC1) {
                ci[1].position = p->position;
                ci[1].velocity = p->velocity;
                ci[1].acceleration = p->treeAcceleration;
                ci[1].w = p->w;
                ci[1].mass = p->mass;
                ci[1].radius = p->soft*2.;
                ci[1].dtCol = p->dtCol;
                ci[1].iOrder = p->iOrder;
                ci[1].iOrderCol = p->iOrderCol;
                ci[1].rung = p->rung;
                }
            else bFoundC1 = 1;
            }
        }

    contribute(2 * sizeof(ColliderInfo), ci, soonestCollReduction, cb);
    }

/**
 * @brief Searches for a particle with a specific iOrder on this TreePiece
 *
 * Contributes a single collider info object that corresponds to the particle
 * specified by iOrder
 *
 * @param iOrder The iOrder of the particle to retrieve
 */
void TreePiece::getCollInfo(int iOrder, const CkCallback& cb)
{
    double dtMin = DBL_MAX;
    ColliderInfo ci;
    ci.iOrder = -1;
    for (unsigned int i=1; i <= myNumParticles; i++) {
        GravityParticle *p = &myParticles[i];
        if (p->iOrder == iOrder) {
            ci.position = p->position;
            ci.velocity = p->velocity;
            ci.acceleration = p->treeAcceleration;
            ci.w = p->w;
            ci.mass = p->mass;
            ci.radius = p->soft*2.;
            ci.dtCol = p->dtCol;
            ci.iOrder = p->iOrder;
            ci.iOrderCol = p->iOrderCol;
            ci.rung = p->rung;
            }
        }
    contribute(sizeof(ColliderInfo), &ci, findCollReduction, cb);
    }

/**
 * @brief Resolves a collision between a particle and a wall, if the particle
 * resides on this tree piece.
 *
 * @param coll A reference to the collision class that handles collision physics
 * @param c1 Information about the particle that is undergoing a collision
 */
void TreePiece::resolveWallCollision(Collision coll, ColliderInfo &c1, 
                                     const CkCallback& cb) {
    GravityParticle *p;
    for (unsigned int i=1; i <= myNumParticles; i++) {
        p = &myParticles[i];
        if (p->iOrder == c1.iOrder) {
            coll.doWallCollision(p);
            break;
            }
        }

    contribute(cb);
    }

/**
 * @brief Resolves a collision between two particles, if either of them resides
 * on this tree piece.
 *
 * @param coll The collision class object that handles collision physics
 * @param c1 Information about the first particle that is undergoing a collision
 * @param c2 Information about the first particle that is undergoing a collision
 * @param bMerge Whether the collision should result in a merger
 * @param baseStep The size of the current step on this rung
 * @param timeNow The current simulation time
 */
void TreePiece::resolveCollision(Collision coll, ColliderInfo &c1,
                                 ColliderInfo &c2, int bMerge, double baseStep,
                                 double timeNow, const CkCallback& cb) {
        GravityParticle *p;
        // Look for the first collider particle on this tree piece
        int bFoundP1 = 0;
        for (unsigned int i=1; i <= myNumParticles; i++) {
            p = &myParticles[i];
            if (p->iOrder == c1.iOrder) {
                bFoundP1 = 1;
                break;
               }
            }

       if (bFoundP1) {
           if (bMerge) {
               if (p->mass >= c2.mass) {
                   CkPrintf("Merging particle %d into %d\n", c2.iOrder, p->iOrder);
                   coll.doCollision(p, c2, 1);
                   coll.setMergerRung(p, c2, c1, baseStep, timeNow);
                   }
               else {
                   CkPrintf("Delete particle %d\n", p->iOrder);
                   deleteParticle(p);
                   }
               }
           else {
               coll.doCollision(p, c2, 0);
               }
           }

        // Look for the second collider particle on this tree piece
        int bFoundP2 = 0;
        for (unsigned int i=1; i <= myNumParticles; i++) {
            p = &myParticles[i];
            if (p->iOrder == c2.iOrder) {
                bFoundP2 = 1;
                break;
                }
            }

       if (bFoundP2) {
           if (bMerge) {
               if (p->mass > c1.mass) {
                   CkPrintf("Merging particle %d into %d\n", c1.iOrder, p->iOrder);
                   coll.doCollision(p, c1, 1);
                   coll.setMergerRung(p, c1, c2, baseStep, timeNow);
                   }
               else {
                   CkPrintf("Delete particle %d\n", p->iOrder);
                   deleteParticle(p);
                   }
               }
           else {
               coll.doCollision(p, c1, 0);
               }
           }
    
    contribute(cb);
    }

/**
 * @brief Check for a merger between two colliders
 *
 * If they are to merge, set the bMergerDelete flag to True for the less
 * massive collider. A merger will occur if the relative velocity between
 * particles is < vEsc and both particles post collision angular speed is
 * < wMax.
 *
 * @param c1 The first collider participating in the event
 * @param c2 The second collider participating in the event
 */
void Collision::checkMerger(ColliderInfo &c1, ColliderInfo &c2)
{

    // Calculate the post collision spin of the merged particle
    Vector3D<double> posNew, vNew, wNew, aNew;
    // Advance particle positions to moment of collision
    Vector3D<double> pAdjust = c1.velocity*c1.dtCol;
    Vector3D<double> cAdjust = c2.velocity*c2.dtCol;
    double radNew;
    c1.position += pAdjust;
    c2.position += cAdjust;
    mergeCalc(c1.radius, c1.mass, c1.position, c1.velocity, c1.acceleration,
              c1.w, &posNew, &vNew, &wNew, &aNew, &radNew, c2);
    // Revert particle positions back to beginning of step
    c1.position -= pAdjust;
    c2.position -= cAdjust;

    double Mtot = c1.mass + c2.mass;
    double vEsc = sqrt(2.*Mtot/(c1.radius + c2.radius));
    double wMax = sqrt(Mtot/(radNew*radNew*radNew));

    double vRel = (c1.velocity - c2.velocity).length();
    if (vRel > vEsc || wNew.length() > wMax) {
        if (!bPerfectAcc) return;
        }

    // Mark the less massive particle to be consumed and deleted
    if (c1.mass < c2.mass) c1.bMergerDelete = 1;
    else c2.bMergerDelete = 1;
    }

/**
 * @brief Determine the time since a particle on a given rung received a kick
 *
 * @param rung The current rung being considered
 * @param baseTime The size of the current time step
 * @param timeNow The current simulation time
 */
double Collision::LastKickTime(int rung, double baseTime, double timeNow)
{
    double rungTime = baseTime/(1<<rung);
    int nRungSteps = (int)(timeNow/rungTime);
    return timeNow - (rungTime*nRungSteps);
    }

/**
 * @brief Correct the velocity of a particle resulting from a merger if the two
 * particles that merged were on different rungs.
 *
 * @param p A reference to the particle resulting from the merger
 * @param c Information about the less massive particle in the merger
 * @param cMerge Information about the more massive particle
 * @param baseStep The timestep size on the lowest rung
 * @param timeNow The current time in the simulation
 */
void Collision::setMergerRung(GravityParticle *p, ColliderInfo &c, ColliderInfo &cMerge,
                              double baseStep, double timeNow)
{
    double m1 = c.mass;
    double m2 = cMerge.mass;
    double m = m1 + m2;
    if (c.rung != cMerge.rung) {
        double lastkick1 = LastKickTime(c.rung, baseStep, timeNow);
        double lastkick2 = LastKickTime(cMerge.rung, baseStep, timeNow);
        p->rung = c.rung;
        if (lastkick1 > lastkick2) {
            p->velocity += (m1/m)*c.acceleration*(lastkick1-lastkick2);
            }
        else {
            p->velocity += (m1/m)*cMerge.acceleration*(lastkick2-lastkick1);
            }
        }
    }

/**
 * @brief Update the velocity of a particle after it undergoes a collision
 * with a wall.
 *
 * @param p A reference to the particle that is undergoing a collision
 */
void Collision::doWallCollision(GravityParticle *p) {
    // TODO: spin
    // TODO: wall overlap
    p->velocity[2] *= -dEpsN;

    Vector3D<double> vPerp = (0., 0., p->velocity[2]);
    Vector3D<double> vParallel = p->velocity - vPerp;
    p->velocity -= vParallel*(1.-dEpsT);
    }

/**
 * @brief Update the state of a particle after it undergoes a collision
 * with another particle.
 *
 * This is done by advancing the particle by a time interval dtCol and then
 * evaluating the collision at the new position. The particle is then moved
 * back to its original position after its velocity has been updated.
 *
 * @param p A reference to the particle that is undergoing a collision
 * @param c Contains information about the particle that p is colliding with
 * @param bMerge If true, merge the particles together
 */
void Collision::doCollision(GravityParticle *p, ColliderInfo &c, int bMerge)
{
    // If dtCol is negative, this particle is overlapping another one. Move
    // the particles apart along their axis of separation
    if (c.dtCol < 0) {
        Vector3D<double> dx = p->position - c.position;
        double overlap = dx.length()-(p->soft*2.)-c.radius;
        p->position -=  c.mass/(p->mass+c.mass)*overlap*dx.normalize();
        return;
        }

    // Advance particle positions to moment of collision
    Vector3D<double> pAdjust = p->velocity*p->dtCol;
    Vector3D<double> cAdjust = c.velocity*c.dtCol;
    Vector3D<double> pAdv = p->position + pAdjust;
    p->position += pAdjust;
    c.position += cAdjust;    

    Vector3D<double> vNew, wNew;
    double radNew;
    if (bMerge) {
        Vector3D<double> posNew, aNew;
        mergeCalc(p->soft*2., p->mass, pAdv, p->velocity, p->treeAcceleration,
                  p->w, &posNew, &vNew, &wNew, &aNew, &radNew, c);
        p->position = posNew;
        p->treeAcceleration = aNew;
        p->soft = radNew/2.;
        p->mass += c.mass;
        }
    else {
        bounceCalc(p->soft*2., p->mass, pAdv, p->velocity, p->w, &vNew, &wNew, c);
        }
    p->velocity = vNew;
    p->w = wNew;

    // Revert particle positions back to beginning of step
    if (bMerge) {
        p->position -= p->velocity*p->dtCol;
        }
    else {
        p->position -= pAdjust;
        c.position -= cAdjust;
        }
    }


/**
 * @brief Calculates the resulting velocity, spin and acceleration of a particle
 * as it merges with another particle.
 *
 * This function writes the new velocity, spin and acceleration to the location 
 * specified by velNew and wNew.
 * 
 * @param r The radius of the particle
 * @param m The mass of the particle
 * @param pos The position of the particle
 * @param vel The velocity of the particle
 * @param acc The acceleration of the particle
 * @param w The spin of the particle
 * @param posNew Where to store the resulting position of the particle
 * @param velNew Where to store the resulting velocity of the particle
 * @param wNew Where to store the resulting spin of the particle
 * @param aNew Where to store the resulting acceleration of the particle
 * @param radNew Where to store the new radius of the merged particle
 * @param c Contains information about the particle that we are colliding with
 */
void Collision::mergeCalc(double r, double m, Vector3D<double> pos,
                          Vector3D<double> vel, Vector3D<double> w,
                          Vector3D<double> acc, Vector3D<double> *posNew,
                          Vector3D<double> *velNew, Vector3D<double> *wNew,
                          Vector3D<double> *aNew,  double *radNew, ColliderInfo &c)
{
    double r1 = r;
    double r2 = c.radius;
    *radNew = pow(r1*r1*r1 + r2*r2*r2,1.0/3); // Conserves volume

    double m1 = m;
    double m2 = c.mass;
    double M = m1 + m2;
    double i1 = 0.4*m1*r1*r1;
    double i2 = 0.4*m2*r2*r2;
    double i = 0.4*M*r*r;

    Vector3D<double> comPos = (m1*pos + m2*c.position)/M;
    Vector3D<double> comVel = (m1*vel + m2*c.velocity)/M;
    Vector3D<double> comAcc = (m1*acc + m2*c.acceleration)/M;

    Vector3D<double> rc1 = pos - comPos;
    Vector3D<double> rc2 = c.position - comPos;
    Vector3D<double> vc1 = vel - comVel;
    Vector3D<double> vc2 = c.velocity - comVel;

    Vector3D<double> angMom = m1*cross(rc1, vc1) + i1*w + m2*cross(rc2, vc2) + i2*c.w;

    *posNew = comPos;    
    *velNew = comVel;
    *wNew = angMom/i;
    *aNew = comAcc;
    }

/**
 * @brief Calculates the resulting velocity and spin of a particle as it
 * bounces off another particle.
 *
 * This function writes the new velocity and spin to the location specified
 * by velNew and wNew.
 * 
 * @param r The radius of the particle
 * @param m The mass of the particle
 * @param pos The position of the particle
 * @param vel The velocity of the particle
 * @param w The spin of the particle
 * @param velNew Where to store the resulting velocity of the particle
 * @param wNew Where to store the resulting spin of the particle
 * @param c Contains information about the particle that we are colliding with
 */
void Collision::bounceCalc(double r, double m, Vector3D<double> pos,
                           Vector3D<double> vel, Vector3D<double> w,
                           Vector3D<double> *velNew, Vector3D<double> *wNew,
                           ColliderInfo &c)
{
    // Equations come from Richardson 1994
    double r1 = r;
    double m1 = m;
    double m2 = c.mass;
    double M = m1 + m2;
    double mu = m1*m2/M;

    Vector3D<double> N = (c.position - pos).normalize();
    Vector3D<double> v = c.velocity - vel;

    Vector3D<double> R1 = N*r1;
    Vector3D<double> sigma1 = cross(w, R1);
    Vector3D<double> R2 = -N*r1;
    Vector3D<double> sigma2 = cross(c.w, R2);
    Vector3D<double> sigma = sigma2 - sigma1;

    Vector3D<double> u = v + sigma;
    Vector3D<double> uN = dot(u, N)*N;
    Vector3D<double> uT = u - uN;

    double alpha = (5./2.)/mu;
    double beta = 1./(1.+(alpha*mu));

    *velNew = vel + m2/M*((1.+dEpsN)*uN + beta*(1-dEpsT)*uT);
    
    double I1 = 2./5.*m1*r1*r1;
    *wNew = w + beta*mu/I1*(1.-dEpsT)*cross(R1, u);

    }

void CollisionSmoothParams::initSmoothParticle(GravityParticle *p)
{
    double v = p->velocity.length();
    p->fBall = coll.dBallFac*dDelta*v + 4.*p->soft;
    }

void CollisionSmoothParams::initSmoothCache(GravityParticle *p1)
{
    p1->dtCol = DBL_MAX;
    p1->iOrderCol = -1;
    }

void CollisionSmoothParams::combSmoothCache(GravityParticle *p1,
                              ExternalSmoothParticle *p2)
{
    if (p2->dtCol < p1->dtCol) {
        p1->dtCol = p2->dtCol;
        p1->iOrderCol = p2->iOrderCol;
        }
    }

/**
 * @brief Do a neighbor search to look for all possible collisions between
 * particles.
 *
 * This function updates the 'dtCol' field for all particles that will undergo
 * a collision in the next time step. A negative value for dtCol indicates that
 * the field was overloaded and represents the total energy of the two particles.
 * In this case, the two particles are stuck in a binary and should be merged.
 *
 * The 'iOrderCol' is normally set to the iOrder of the particle which this
 * particle is going to collide with. If the collision is with a wall,
 * 'iOrderCol' is set to -2.
 */
void CollisionSmoothParams::fcnSmooth(GravityParticle *p, int nSmooth,
                                      pqSmoothNode *nList)
{
    GravityParticle *q;
    int i;    
    double rq, sr, D, dt, rdotv, vRel2, dr2;
    p->dtCol = DBL_MAX;
    p->iOrderCol = -1;
    if (TYPETest(p, TYPE_DELETED)) return;
    if (bWall) {
        double dt = (dWallPos - p->position[2] - (p->soft*2.))/p->velocity[2];
        if (dt > 0. && dt < dDelta) {
            p->dtCol = dt;
            p->iOrderCol = -2;
            }
    }

    for (i = 0; i < nSmooth; i++) {
        q = nList[i].p;

        // Skip self interactions
        if (p->iOrder == q->iOrder) continue;
        // Ignore deleted particles
        if (TYPETest(q, TYPE_DELETED)) continue;

        Vector3D<double> dx = p->position - q->position;
        Vector3D<double> vRel = p->velocity - q->velocity;

        sr = (p->soft*2.) + (q->soft*2.);
        rdotv = dot(dx, vRel);
        vRel2 = vRel.lengthSquared();
        dr2 = dx.lengthSquared() - sr*sr;
        D = rdotv*rdotv - dr2*vRel2;
        // If D <= 0, continue?
        D = sqrt(D);
        dt = (-rdotv - D)/vRel2;

        if (dt > 0. && dt < dDelta && dt < p->dtCol) {
            p->dtCol = dt;
            p->iOrderCol = q->iOrder;
            }

       // Because there is no gravitational softening, particles can get stuck in
       // binaries. If the particle is not already undergoing a collision, check
       // to see if it is bound to its neighbor. If it is, mark it for a merge
       if (bAllowMergers && p->iOrderCol == -1 && p->rung >= iMinBinaryRung) {
           double pe = -q->mass/dx.length();
           double ke = 0.5*vRel.lengthSquared();
           if ((ke + pe) < 0.) {
               /* This is quick, but not optimal. We need to be sure we don't
               automatically merge any bound particles, as those on highly
               eccentric orbits may be distrupted at apocenter. Therefore
               we assume that these particles will only reach this point of
               the code near pericenter (due to iMinBinaryRung), and we can 
               safely exclude particles that do not meet the following 
               criterion. Some fiddling with iMinBinaryRung and dMaxBinaryEcc
               may be necessary to acheive an appropriate balance of merging.
               */
               if (vRel2 < sqrt((1 + dMaxBinaryEcc)/(1 - dMaxBinaryEcc))*(ke - pe)) {
                   p->dtCol = ke + pe;
                   p->iOrderCol = q->iOrder;
                   }
               }
           }
        }
    }
