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
    
    do {
        bHasCollision = 0;

        CollisionSmoothParams pCS(TYPE_DARK, 0, dTime, dDelta, 
           param.collision->bWall, param.collision->dWallPos, param.collision);
        double dfBall2OverSoft2 = 4.0*param.dhMinOverSoft*param.dhMinOverSoft;
        treeProxy.startSmooth(&pCS, 0, param.collision->nSmoothCollision,
                  dfBall2OverSoft2, CkCallbackResumeThread());

        CkReductionMsg *msgChk;
        treeProxy.getCollInfo(CkCallbackResumeThread((void*&)msgChk));
        ColliderInfo *c = (ColliderInfo *)msgChk->getData();
        if (c[0].dtCol <= dDelta) {
            bHasCollision = 1;
            // Collision with wall
            if (c[0].iOrderCol == -2 && param.collision->bWall) {
                nColl++;
                treeProxy.resolveWallCollision(*(param.collision), c[0], 
                                               CkCallbackResumeThread());
                }
            // Collision with particle
            else {
                if (c[0].dtCol != c[1].dtCol) {
                    CkPrintf("%f %f\n", c[0].dtCol, c[1].dtCol);
                    CkAbort("Warning: Collider pair mismatch\n");
                    }
                nColl++;
                
                if (param.collision->bAllowMergers) {
                    param.collision->checkMerger(c[0], c[1]);
                    }
                if (c[0].bMergerDelete || c[1].bMergerDelete) {
                    treeProxy.resolveMerger(c[0], c[1], CkCallbackResumeThread());
                    }
                else {
                    treeProxy.resolveCollision(*(param.collision), c[0], c[1], 
                                               CkCallbackResumeThread());
                    }
                }
            }

        delete msgChk;
        } while (bHasCollision);
        CkPrintf("Resolved %d collisions\n", nColl);

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
            ci[0].w = p->w;
            ci[0].mass = p->mass;
            ci[0].dtCol = p->dtCol;
            ci[0].iOrder = p->iOrder;
            ci[0].iOrderCol = p->iOrderCol;
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
                ci[1].w = p->w;
                ci[1].mass = p->mass;
                ci[1].dtCol = p->dtCol;
                ci[1].iOrder = p->iOrder;
                ci[1].iOrderCol = p->iOrderCol;
                }
            else bFoundC1 = 1;
            }
        }

    contribute(2 * sizeof(ColliderInfo), ci, soonestCollReduction, cb);
    }

/**
 * @brief Resolves a collision between a particle and a wall, if the particle
 * resides on this tree piece.
 *
 * @param coll A reference to the collision class that handles collision physics
 * @param c1 Information about the particle that is undergoing a collision
 */
void TreePiece::resolveWallCollision(Collision &coll, ColliderInfo &c1, 
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
 */
void TreePiece::resolveCollision(Collision &coll, ColliderInfo &c1,
                                 ColliderInfo &c2, const CkCallback& cb) {
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

       if (bFoundP1) coll.doCollision(p, c2);

        // Look for the second collider particle on this tree piece
        int bFoundP2 = 0;
        for (unsigned int i=1; i <= myNumParticles; i++) {
            p = &myParticles[i];
            if (p->iOrder == c2.iOrder) {
                bFoundP2 = 1;
                break;
                }
            }

       if (bFoundP2) coll.doCollision(p, c1);
    
    contribute(cb);
    }

/**
 * @brief Resolve a merger between two colliders, if either of them exist on
 * this TreePiece.
 *
 * Particle 1 takes on the mass of particle 2. Particle 2 is deleted.
 *
 * @param c1 Information about particle 1
 * @param c2 Information about particle 2
 */
void TreePiece::resolveMerger(ColliderInfo &c1, ColliderInfo &c2, const CkCallback& cb)
{
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
           CkPrintf("Updating mass of particle %d after merger\n", c1.iOrder);
           p->mass += c2.mass;
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
           CkPrintf("Deleting particle %d after merger\n", c2.iOrder);
           deleteParticle(p);
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
    double Mtot = c1.mass + c2.mass;
    double Rtot = c1.radius + c2.radius;
    double vEsc = sqrt(2.*Mtot/Rtot);
    double wMax = sqrt(Mtot/(Rtot*Rtot*Rtot));

    // Calculate the post collision spin of each particle and compare with wMax
    Vector3D<double> vNew;
    Vector3D<double> wNew1, wNew2;
    bounceCalc(c1.radius, c1.mass, c1.position, c1.velocity, c1.w, &vNew, &wNew1, c2);
    bounceCalc(c2.radius, c2.mass, c2.position, c2.velocity, c2.w, &vNew, &wNew2, c1);

    double vRel = (c1.velocity - c2.velocity).length();
    if (vRel > vEsc || wNew1.length() > wMax || wNew2.length() > wMax) return;

    // Mark the less massive particle to be consumed and deleted
    if (c1.mass < c2.mass) c1.bMergerDelete = 1;
    else c2.bMergerDelete = 1;
    }

/**
 * @brief Update the velocity of a particle after it undergoes a collision
 * with a wall.
 *
 * @param p A reference to the particle that is undergoing a collision
 */
void Collision::doWallCollision(GravityParticle *p) {
    // TODO: spin
    p->velocity[2] *= -dEpsN;

    Vector3D<double> vPerp = (0., 0., p->velocity[2]);
    Vector3D<double> vParallel = p->velocity - vPerp;
    p->velocity -= vParallel*(1.-dEpsT);
    }

/**
 * @brief Update the velocity of a particle after it undergoes a collision
 * with another particle.
 *
 * This is done by advancing the particle by a time interval dtCol and then
 * evaluating the collision at the new position. The particle is then moved
 * back to its original position after its velocity has been updated.
 *
 * @param p A reference to the particle that is undergoing a collision
 * @param c Contains information about the particle that p is colliding with
 */
void Collision::doCollision(GravityParticle *p, ColliderInfo &c)
{
    // Advance particle positions to moment of collision
    Vector3D<double> pAdjust = p->velocity*p->dtCol;
    Vector3D<double> cAdjust = c.velocity*c.dtCol;
    p->position += pAdjust;
    c.position += cAdjust;    

    // Bounce
    Vector3D<double> vNew;
    Vector3D<double> wNew;
    bounceCalc(p->soft/2., p->mass, p->position, p->velocity, p->w, &vNew, &wNew, c);
    p->velocity = vNew;
    p->w = wNew;

    // Revert particle positions back to beginning of step
    p->position -= pAdjust;
    c.position -= cAdjust;
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
 * a collision in the next time step.
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
    CkPrintf("Particle smooth %d\n", p->iOrder);
    if (bWall) {
        double dt = (dWallPos - p->position[2] - (p->soft/2.))/p->velocity[2];
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

        Vector3D<double> vRel = p->velocity - q->velocity;
        Vector3D<double> dx = p->position - q->position;

        sr = (p->soft/2.) + (q->soft/2.);
        rdotv = dot(dx, vRel);
        if (rdotv >= 0) continue; // Particles moving apart
        vRel2 = vRel.lengthSquared();
        dr2 = dx.lengthSquared() - sr*sr;
        D = rdotv*rdotv - dr2*vRel2;
        if (D <= 0.0) continue; // Not on a collision trajectory
        D = sqrt(D);
        dt = (-rdotv - D)/vRel2;

        if (dt < 0.) {
            CkAbort("Warning: collision found with dt < 0. Particles overlapping?\n");
            }

        if (dt < dDelta && dt < p->dtCol) {
            p->dtCol = dt;
            p->iOrderCol = q->iOrder;
            }
        }
    }
