/*
 * Collision detection and response module
 */
#include <float.h>
#include "ParallelGravity.h"
#include "collision.h"
#include "smooth.h"

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

    CollisionSmoothParams pCS(TYPE_DARK, 0, dTime, dDelta, param.collision);
    double dfBall2OverSoft2 = 4.0*param.dhMinOverSoft*param.dhMinOverSoft;
    treeProxy.startSmooth(&pCS, 0, param.collision->nSmoothCollision,
              dfBall2OverSoft2, CkCallbackResumeThread());

    CkReductionMsg *msgChk;
    treeProxy.getCollInfo(CkCallbackResumeThread((void*&)msgChk));
    CkReduction::setElement *current = (CkReduction::setElement*) msgChk->getData();

    // Find minimum collision time
    double dtMin = DBL_MAX;
    while (current != NULL) {
        ColliderInfo* c = (ColliderInfo*) &current->data;
        if (c->dtCol < dtMin) dtMin = c->dtCol;
         current = current->next();
        }

    // Get collision info for the two soonest colliders
    if (dtMin < DBL_MAX) {
        ColliderInfo *c1, *c2;
        int bFoundC1 = 0;
        current = (CkReduction::setElement*) msgChk->getData();
        while (current != NULL) {
            ColliderInfo* c = (ColliderInfo*) &current->data;
            if (c->dtCol == dtMin) {
                if (!bFoundC1) {
                    c1 = c;
                    bFoundC1 = 1;
                    }
                else c2 = c;
                }
            current = current->next();
            }
        CkPrintf("Imminent collision in %f\n", c1->dtCol);
        treeProxy.resolveCollision(*(param.collision), *c1, *c2, CkCallbackResumeThread());
        }

    delete msgChk;
    }

void TreePiece::getCollInfo(const CkCallback& cb)
{
    double dtMin = DBL_MAX;
    ColliderInfo ci;
    for (unsigned int i=1; i <= myNumParticles; i++) {
        GravityParticle *p = &myParticles[i];
        if (p->dtCol < dtMin) {
            dtMin = p->dtCol;
            ci.position = p->position;
            ci.velocity = p->velocity;
            ci.w = p->w;
            ci.mass = p->mass;
            ci.dtCol = p->dtCol;
            ci.iOrder = p->iOrder;
            }
        }
    
    contribute(sizeof(ColliderInfo), &ci, CkReduction::set, cb);
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
    //CkPrintf("Handling collision between %d and %d\n", p->iOrder, c.iOrder);
    // Advance particle positions to moment of collision
    Vector3D<double> pAdjust = p->velocity*p->dtCol;
    Vector3D<double> cAdjust = c.velocity*c.dtCol;
    p->position += pAdjust;
    c.position += cAdjust;

    //CkPrintf("Velocity before bounce: %f, %f, %f\n", p->velocity.x, p->velocity.y, p->velocity.z);
    bounce(p, c, 1.0, 1.0);
    //CkPrintf("Velocity after bounce: %f, %f, %f\n", p->velocity.x, p->velocity.y, p->velocity.z);

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
    //CkPrintf("CombSmoothCache called\n");
    if (p2->dtCol < p1->dtCol) {
        //CkPrintf("Setting p1->dtCol to %f\n", p2->dtCol);
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
    CkPrintf("Collision smooth with %d neighbors\n", nSmooth);
    for (i = 0; i < nSmooth; i++) {
        q = nList[i].p;

        CkPrintf("Compare particle %d with particle of mass %g\n", p->iOrder, q->mass);

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
