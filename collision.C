/*
 * Collision detection and response module
 */
#include <float.h>
#include "ParallelGravity.h"
#include "collision.h"
#include "smooth.h"
#include "Reductions.h"

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
    // Plan: Call initial smooth to find collisions
    //       While there are collisions in the next timestep (treeProxy.getNextColl)
    //           Resolve the soonest collision (treeProxy.resolveNextColl)
    //             See line 738 of feedback.C for an example of how to send back
    //             info from one tree piece
    //           Rebuild tree?
    //           Resmooth

    int bHasCollision;
    
    do {
        bHasCollision = 0;

        CollisionSmoothParams pCS(TYPE_DARK, 0, dTime, dDelta, param.collision);
        double dfBall2OverSoft2 = 4.0*param.dhMinOverSoft*param.dhMinOverSoft;
        treeProxy.startSmooth(&pCS, 0, param.collision->nSmoothCollision,
                  dfBall2OverSoft2, CkCallbackResumeThread());

        CkReductionMsg *msgChk;
        treeProxy.getCollInfo(CkCallbackResumeThread((void*&)msgChk));
        ColliderInfo *c = (ColliderInfo *)msgChk->getData();
        if (c[0].dtCol <= dDelta) {
            bHasCollision = 1;
            CkPrintf("Imminent collision in %f %f\n", c[0].dtCol, c[1].dtCol);
            treeProxy.resolveCollision(*(param.collision), c[0], c[1], CkCallbackResumeThread());
            }

        delete msgChk;
        } while (bHasCollision);

    }

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
                }
            else bFoundC1 = 1;
            }
        }
    contribute(2 * sizeof(ColliderInfo), ci, soonestCollReduction, cb);
    }

void TreePiece::resolveCollision(Collision &coll, ColliderInfo &c1, ColliderInfo &c2, const CkCallback& cb) {
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

void Collision::doCollision(GravityParticle *p, ColliderInfo &c)
{
    CkPrintf("Handling collision between %d and %d\n", p->iOrder, c.iOrder);
    // Advance particle positions to moment of collision
    Vector3D<double> pAdjust = p->velocity*p->dtCol;
    Vector3D<double> cAdjust = c.velocity*c.dtCol;
    p->position += pAdjust;
    c.position += cAdjust;

    CkPrintf("Velocity before bounce: %f, %f, %f\n", p->velocity.x, p->velocity.y, p->velocity.z);
    bounce(p, c, 1.0, 1.0);
    CkPrintf("Velocity after bounce: %f, %f, %f\n", p->velocity.x, p->velocity.y, p->velocity.z);

    // Revert particle positions back to beginning of step
    p->position -= pAdjust;
    c.position -= cAdjust;
    }


void Collision::bounce(GravityParticle *p, ColliderInfo &c, double dEpsN, double dEpsT)
{
    Vector3D<double> p1, q, u, un, ut, n, s1, s2, v, s;
    double a, b, c1, m1, m2, m, mu, alpha, beta, r1, r2, i1, i2;

    m1 = p->mass;
    m2 = c.mass;
    m = m1+m2;
    mu = m1*m2/m;
    r1 = p->soft/2.;
    r2 = c.radius;
    i1 = 0.4*m1*r1*r1;
    i2 = 0.4*m2*r2*r2;
    alpha = 2.5*(1/m1 + 1/m2);
    beta = 1/(1 + alpha*mu);
 
    n = (c.position - p->position);
    a = n.lengthSquared();

	a = 1/sqrt(a);
	n *= a;

    s1 = r1*cross(c.w, n);
    s2 = -r2*cross(c.w, n);

    v = c.velocity - p->velocity;
    s = s2 - s1;
    u = v + s;
    a = dot(u, n);

    if (a >= 0.)
        CkPrintf("%d, %d - near miss? a = %f\n", p->iOrder, c.iOrder, a);

    un = a*n;
    ut = u - un;

    a = (1 + dEpsN);
	b = beta*(1 - dEpsT);
	p1 = a*un + b*ut;

	a = mu*b;
    q = a*cross(n, u);

    a =   m2/m;
	b = -m1/m;
    c1 = r1/i1;

    p->velocity += a*p1;
    p->w += c1*q;

    p->dtCol = DBL_MAX;
    p->iOrderCol = -1;
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
        }
    if (p2->iOrderCol != -1 && p1->iOrderCol == -1) {
        p1->iOrderCol = p2->iOrderCol;
        }
    }

void CollisionSmoothParams::fcnSmooth(GravityParticle *p, int nSmooth, pqSmoothNode *nList)
{
    GravityParticle *q;
    int i;    
    double rq, sr, D, dt, rdotv, vRel2, dr2;

    double rp = p->soft/2.;
    for (i = 0; i < nSmooth; i++) {
        q = nList[i].p;

        // Skip self interactions
        if (p->iOrder == q->iOrder) continue;
        CkAssert(TYPETest(q, TYPE_DARK));

        Vector3D<double> vRel = p->velocity - q->velocity;
        sr = (p->soft/2.) + (q->soft/2.);
        rdotv = dot(nList[i].dx, vRel);
        if (rdotv >= 0) continue; // Particles moving apart
        vRel2 = vRel.lengthSquared();
        dr2 = nList[i].dx.lengthSquared() - sr*sr;
        D = rdotv*rdotv - dr2*vRel2;
        if (D <= 0.0) continue; // Not on a collision trajectory
        D = sqrt(D);
        dt = (-rdotv - D)/vRel2;

        if (dt < dDelta && dt < p->dtCol) {
            p->dtCol = dt;
            p->iOrderCol = q->iOrder;
            }
        }
    }
