#include "ParallelGravity.h"
#include "externalGravity.h"

///
// @brief initalize parameters for external potential field
//

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
 * @brief This function applies the external potential force to every applicable
 * particle on this TreePiece.
 *
 * Applies an acceleration on every particle on the current tree piece that is on
 * or above the current rung. This function also keeps track of the acceleration on 
 * the potential imparted by the particles.
 *
 * @param iKickRung The current rung that we are on
 * @param exGrav A reference to the ExternalGravity class
 *
 * @return The accumulated acceleration on the potential by the particles on this
 * TreePiece
 */
void TreePiece::externalGravity(int iKickRung, const ExternalGravity exGrav,
                                const CkCallback& cb)
{
    CkAssert(bBucketsInited);

    // Keep track of the forces of the particles on the external potential
    double frameAcc[3];
    frameAcc[0] = 0.0;
    frameAcc[1] = 0.0;
    frameAcc[2] = 0.0;
    Vector3D<double> pFrameAcc;

    for(unsigned int i = 1; i <= myNumParticles; ++i) {
        GravityParticle *p = &myParticles[i];
        if(p->rung >= iKickRung) {
            pFrameAcc = exGrav.applyPotential(p);
            frameAcc[0] += pFrameAcc[0];
            frameAcc[1] += pFrameAcc[1];
            frameAcc[2] += pFrameAcc[2];
            }
        }
    contribute(sizeof(frameAcc), frameAcc, CkReduction::sum_double, cb);
    }

/*
 * @brief This function applies the external potential force to a specific particle
 * by updating its 'treeAcceleration'.
 *
 * @param p The particle to apply the potential to
 *
 * @return The acceleration on the potential by particle p
 */
Vector3D<double> ExternalGravity::applyPotential(GravityParticle *p) const
{
    Vector3D<double> pFrameAcc(0., 0., 0.);
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
        Vector3D<double> thetaVec(px*pz/c, py*pz/c, -(px*px + py*py)/c);
        Vector3D<double> a = ar*rVec + atheta*thetaVec;

        p->treeAcceleration += a;
        pFrameAcc = -a*p->mass/dCentMass;

        double idt2 = fabs(ar/r);
        if(idt2 > p->dtGrav)
            p->dtGrav = idt2;
        }
    return pFrameAcc;
    }

/*
 * @brief For use when transforming into an accelerating frame of reference. Given
 * the acceleration of the frame, this function applies the opposite of this
 * acceleration to all particles on this TreePiece.
 *
 * @param frameAcc The acceleration on the frame at the current time
 */
void TreePiece::applyFrameAcc(int iKickRung, Vector3D<double> frameAcc, const CkCallback& cb)
{
    for (unsigned int i = 1; i <= myNumParticles; ++i) {
        GravityParticle *q = &myParticles[i];
        if(q->rung >= iKickRung)
            q->treeAcceleration -= frameAcc; 
        }

    contribute(cb);
    }
