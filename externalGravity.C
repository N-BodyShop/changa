#include "ParallelGravity.h"
#include "externalGravity.h"

void ExternalGravity::AddParams(PRM prm)
{
    bBodyForce = 0;
    prmAddParam(prm,"bBodyForce",paramBool,&bBodyForce,
                sizeof(int),"bodyforce","use constant body force = -bf");
    dBodyForceConst = 0.0;
    prmAddParam(prm,"dBodyForceConst",paramDouble,&dBodyForceConst,
                sizeof(double),"bodyforceconst",
                "strength of constant bodyforce = 0");

    dCentMass = 1.0;
    prmAddParam(prm,"dCentMass",paramDouble,&dCentMass,
                sizeof(double),
                "fgm","specifies the central mass for Keplerian orbits");
    bPatch = 0;
    prmAddParam(prm,"bPatch",paramBool,&bPatch,
                sizeof(int),
                "patch","enable/disable patch reference frame = -patch");
    dOrbDist = 0.0;
    prmAddParam(prm,"dOrbDist",paramDouble,&dOrbDist,
                sizeof(double),"orbdist","<Patch orbital distance>");

    bCentralBody = 0;
    prmAddParam(prm,"bCentralBody",paramBool,&bCentralBody,
                sizeof(int), "centralbody","<enable/disable central body potential> = 0");
    dEqRad = 1.0;
    prmAddParam(prm,"dEqRad",paramDouble,&dEqRad,
                sizeof(double),"dEqRad","<Equatorial radius of central body> = 1.0");
    dJ2 = 0.0;
    prmAddParam(prm,"dJ2",paramDouble,&dJ2,
                sizeof(double),"dJ2","<Oblateness coefficient J2 of central body> = 0.0");
    dJ4 = 0.0;
    prmAddParam(prm,"dJ4",paramDouble,&dJ4,
                sizeof(double),"dJ4","<Oblateness coefficient J4 of central body> = 0.0");
    dJ6 = 0.0;
    prmAddParam(prm,"dJ6",paramDouble,&dJ6,
                sizeof(double),"dJ6","<Oblateness coefficient J6 of central body> = 0.0");
    }

void ExternalGravity::CheckParams(PRM prm, struct parameters &param)
{
    // Enable external gravity if any of the flags are set
    if (bBodyForce || bPatch || bCentralBody)
        param.bDoExternalGravity = 1;
    }

/*
 * Apply an external gravitational force
 */
void TreePiece::externalGravity(int iKickRung, const ExternalGravity exGrav,
                                const CkCallback& cb)
{
    CkAssert(bBucketsInited);
    for(unsigned int i = 1; i <= myNumParticles; ++i) {
        GravityParticle *p = &myParticles[i];
        if(p->rung >= iKickRung)
            exGrav.applyPotential(p);
        }
    contribute(cb);
    }

void ExternalGravity::applyPotential(GravityParticle *p) const
{
    if (bBodyForce) {
        if(p->position.z > 0.0) {
            p->treeAcceleration.z -= dBodyForceConst;
            p->potential += dBodyForceConst*p->potential;
            double idt2 = dBodyForceConst/p->position.z;
            if(idt2 > p->dtGrav)
                p->dtGrav = idt2;
            }
        else {
            p->treeAcceleration.z += dBodyForceConst;
            p->potential -= dBodyForceConst*p->position.z;
            if (p->position.z != 0.0) {
                double idt2 = -dBodyForceConst/p->position.z;
                if (idt2 > p->dtGrav)
                    p->dtGrav = idt2;
                }
            }
    }

    if (bPatch) {
        double r2 = dOrbDist*dOrbDist + p->position.z*p->position.z;
        double idt2 = dCentMass*pow(r2, -1.5);

        p->treeAcceleration.z -= dCentMass*p->position.z
                                 *pow(r2, -1.5);
        p->potential += dCentMass/sqrt(r2);
        if(idt2 > p->dtGrav)
            p->dtGrav = idt2;
        }

    if(bCentralBody) {
        double px = p->position.x;
        double py = p->position.y;
        double pz = p->position.z;
        double r = p->position.length();
        
        // Legendre polynomials
        double c1 = pz/r;                  // cos(theta)
        double c2 = sqrt(px*px + py*py)/r; // sin(theta)
        double p2 = 0.5*(3.*pow(c1, 2) - 1.);
        double p4 = 1./8.*(35.*pow(c1, 4) - 30.*pow(c1, 2) + 3.);
        double p6 = 1./16.*(231.*pow(c1, 6) - 315.*pow(c1, 4)
                              + 105.*pow(c1, 2) - 5.);
        
        // Theta derivatives of legendre polynomials
        double p2prime = -3.*c1*c2;
        double p4prime = -5./16.*(4.*c1*c2 + 28.*pow(c1, 2)*c2*(2.*pow(c1, 2) - 1.));
        double p6prime = -1./16.*(1386.*c2*pow(c1, 5)
                                 - 1260.*c2*pow(c1, 3)+ 210.*c1*c2);

        double a2 = dJ2*pow(dEqRad/r, 2);
        double a4 = dJ4*pow(dEqRad/r, 4);
        double a6 = dJ6*pow(dEqRad/r, 6);

        p->potential += -dCentMass/r*(1. - a2*p2 - a4*p4 - a6*p6);

        // Acceleration in spherical coordinates
        double ar = -dCentMass/pow(r, 2)*(1. - 3.*a2*p2 -5.*a4*p4 - 7.*a6*p6);
        double atheta = -dCentMass/pow(r, 2)*(a2*p2prime
                          + a4*p4prime + a6*p6prime);

        Vector3D<double> rVec = p->position/r;
        double c = (r*sqrt(px*px+py*py));
        Vector3D<double> thetaVec = (px*pz/c, py*pz/c, -(px*px + py*py)/c);

        p->treeAcceleration += ar*rVec + atheta*thetaVec;

        double idt2 = ar/r;
        if(idt2 > p->dtGrav)
            p->dtGrav = idt2;
        }
    }
