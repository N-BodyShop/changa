#include "ParallelGravity.h"
#include "externalForce.h"
#include "physconst.h"

///
// @brief initalize parameters for external potential field
//

void ExternalForce::AddParams(PRM prm)
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
                "mcent","specifies the central mass for Keplerian orbits = 1 M_sun");
    dCentRad = 0.0046; // 1 R_sun in AU
    prmAddParam(prm,"dCentRad",paramDouble,&dCentRad,
                sizeof(double),
                "rcent","central radius of the star, delete bodies which will fall inside of this = 1 R_sun");
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

    bLogarithmicHalo = 0;
    prmAddParam(prm,"bLogarithmicHalo",paramBool,&bLogarithmicHalo,
                sizeof(int), "logarithmichalo","<Type of halo is logarithmic> = 0");
     
    bDoGasDrag = 0;
    prmAddParam(prm, "bDoGasDrag", paramBool, &bDoGasDrag,
        sizeof(int), "bDoGasDrag", "<Apply external gas drag force to planetesimals> = 0");
    dRadInfGas = 1.0;
    prmAddParam(prm, "dRadInfGas", paramDouble, &dRadInfGas,
        sizeof(double), "dRadInfGas", "<Inflation factor for radius of particles, used for gas drag calculations> = 1.0");
    dSigma0 = 1700.0;
    prmAddParam(prm, "dSigma0", paramDouble, &dSigma0,
        sizeof(double), "dSigma0", "<Gas surface density at 1 AU (in g/cm^2)> = 1700.0");
    bConstGasProf = 0;
    prmAddParam(prm, "bConstGasProf", paramBool, &bConstGasProf,
        sizeof(int), "bConstGasProf", "<Constant volume density for gas (ignores surface density profile, uses T0 for gas temp)> = 0");
    dConstGasRho = 1e-6;
    prmAddParam(prm, "dConstGasRho", paramDouble, &dConstGasRho,
        sizeof(double), "dConstGasRho", "<Value to use for constant gas volume density (in g/cm^3) = 1e-6");
    dP = 1.5;
    prmAddParam(prm, "dP", paramDouble, &dP,
        sizeof(double), "dP", "<Power law slope of gas surface density profile> = 1.5");
    dQ = 0.5;
    prmAddParam(prm, "dQ", paramDouble, &dQ,
        sizeof(double), "dQ", "<Power law slope of gas temperature profile> = 0.5");
    dT0 = 280.0;
    prmAddParam(prm, "dT0", paramDouble, &dT0,
        sizeof(double), "dT0", "<Gas temperature at 1 AU (in K)> = 280.0");
    dMu = 2.34;
    prmAddParam(prm, "dMu", paramDouble, &dMu,
        sizeof(double), "dMu", "<Mean molecular weight of gas> = 2.34");
    dCD = 2;
    prmAddParam(prm, "dCD", paramDouble, &dCD,
        sizeof(double), "dCD", "<Coefficient of gas drag force> = 2.0");
    }

void ExternalForce::CheckParams(PRM prm, struct parameters &param)
{
    // Enable external force if any of the flags are set
    if (bBodyForce || bPatch || bCentralBody || param.bDoExternalGravity)
        param.bDoExternalForce = 1; 
    // Gas drag requires central point mass
    if (bDoGasDrag && !bCentralBody)
        CkAbort("bDoGasDrag requires bCentralBody to be set\n");
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
 * @param exForce A reference to the ExternalForce class
 *
 * @return The accumulated acceleration on the potential by the particles on this
 * TreePiece
 */
void TreePiece::externalForce(int iKickRung, const ExternalForce& exForce, int bKepStep,
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
            pFrameAcc = exForce.applyGravPotential(p, bKepStep);

            if (exForce.bDoGasDrag) {
                // We don't care about the backreaction on the gas so this
                // force doesn't get applied to the frame
                exForce.applyGasDrag(p);
                }

            frameAcc[0] += pFrameAcc[0];
            frameAcc[1] += pFrameAcc[1];
            frameAcc[2] += pFrameAcc[2];
            }
        }
    contribute(sizeof(frameAcc), frameAcc, CkReduction::sum_double, cb);
    }

/*
 * @brief This function applies the gas drag force to a specific particle
 * by updating its 'treeAcceleration'.
 *
 * @param p The particle to apply the drag force to
 */
void ExternalForce::applyGasDrag(GravityParticle *p) const
{
    // Apply gas drag to dark particles only
    if (!TYPETest(p, TYPE_DARK)) return;

    // Get the cylindrical coordinates of the planetesimal
    // Assume the gas disk lies in the x-y plane
    Vector3D<double> rVec = p->position;
    Vector3D<double> vVec = p->vPred();
    double r = sqrt(rVec[0]*rVec[0] + rVec[1]*rVec[1]);
    double z = rVec[2];

    // Calculate the local sound speed
    double temp = dT0*pow(r, -dQ);

    if (bConstGasProf) temp = dT0;

    double cGas = sqrt(KBOLTZ*temp/(dMu*MHYDR))/CMPERAU*SECONDSPERYEAR/(2*M_PI);

    // Calculate the local volume density of the gas
    double sigmaGas = dSigma0*pow(r, -dP)/MSOLG*CMPERAU*CMPERAU;
    double omega = sqrt(dCentMass/(r*r*r));
    double hGas = cGas/omega;
    // Morishima 2010 eq 3
    double rhoGas = sigmaGas/(sqrt(2*M_PI)*hGas)*exp(-(z*z)/(2*hGas*hGas));

    if (bConstGasProf) rhoGas = dConstGasRho/MSOLG*CMPERAU*CMPERAU*CMPERAU;

    // Force balance due to pressure gradient
    double a1 = -dP + (dQ/2 - 3.0/2.0)*(1 - (z*z/(hGas*hGas))) - dQ;
    // Force balance due to gravity of star
    double a2 = dCentMass*r/pow(r*r + z*z, 3.0/2.0);
    // Morishima 2010 eq 2
    double gasSpeed = sqrt(cGas*cGas*a1 + r*a2);

    // Gas velocity should point in phi direction
    double phi = atan2(rVec[1], rVec[0]);
    double phiX = -sin(phi);
    double phiY = cos(phi);
    Vector3D<double> thetaHat(phiX, phiY, 0);
    Vector3D<double> vGas = gasSpeed*thetaHat;

    double rPl = p->soft*2/dRadInfGas;
    Vector3D<double> vRel = vVec - vGas;
    // Morishima 2010 eq 1
    Vector3D<double> aDrag = -1/(2*p->mass)*dCD*M_PI
                             *rPl*rPl*rhoGas*vRel*vRel.length();
    p->treeAcceleration += aDrag;
    }


/*
 * @brief This function applies the gravitational potential force to a specific particle
 * by updating its 'treeAcceleration'.
 *
 * If the particle is in a Keplerian potential, it gets deleted here if the pericenter distance
 * falls inside the radius of the central body (dCentRad)
 *
 * @param p The particle to apply the potential to
 *
 * @return The acceleration on the potential by particle p
 */
Vector3D<double> ExternalForce::applyGravPotential(GravityParticle *p, int bKepStep) const
{
    Vector3D<double> pFrameAcc(0., 0., 0.);
    if (bBodyForce) {
        if(p->position.z > 0.0) {
            p->treeAcceleration.z -= dBodyForceConst;
            p->potential += dBodyForceConst*p->potential;
            double idt2 = dBodyForceConst/p->position.z;
            if(idt2 > p->dtGrav && !bKepStep)
                p->dtGrav = idt2;
            }
        else {
            p->treeAcceleration.z += dBodyForceConst;
            p->potential -= dBodyForceConst*p->position.z;
            if (p->position.z != 0.0) {
                double idt2 = -dBodyForceConst/p->position.z;
                if (idt2 > p->dtGrav && !bKepStep)
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

#ifdef COLLISION
        double mu = p->mass + dCentMass;
        Vector3D<double> h = cross(p->position, p->vPred());
        Vector3D<double> tmp = cross(p->vPred(), h);
        Vector3D<double> ecc = (tmp/mu) - (p->position/p->position.length());
        double ecc_val = ecc.length();

        double sma = dot(h, h)/(mu*(1 - (ecc_val*ecc_val)));

        double rPeri = sma*(1 - ecc_val);
        double arPeri = -dCentMass/pow(rPeri, 2);
        
        // Only allow dtKep to be set if its uninitalized
	// This should only happen at the beginning of a simulation
        if (p->dtKep == 0) {
          p->dtKep = fabs(arPeri/rPeri);
          }

        if (rPeri < dCentRad) { 
            CkPrintf("Particle %ld will collide with star, deleting\n", p->iOrder);
            deleteParticle(p);
            }
#endif

        }
        
    if(bLogarithmicHalo) {
        double px = p->position.x;
        double py = p->position.y;
        double pz = p->position.z;
        double r = p->position.length();
        
        // Form of Logarithmic Halo Potential
        double v0 = 1.0;
        double qy = 0.9;
        double qz = 0.7;
        double core = sqrt(0.1);
        p->potential += v0*v0*log(px*px + (py/qy)*(py/qy) + (pz/qz)*(pz/qz) + core*core); // Written as double the true value to account for divide-by-two issue in ChaNGa integration
        
        // Acceleration in Cartesian coordinates
        double ax = -v0*v0*px/(px*px + (py/qy)*(py/qy) + (pz/qz)*(pz/qz) + core*core);
        double ay = -v0*v0*py/((qy*qy)*(px*px + (py/qy)*(py/qy) + (pz/qz)*(pz/qz) + core*core));
        double az = -v0*v0*pz/((qz*qz)*(px*px + (py/qy)*(py/qy) + (pz/qz)*(pz/qz) + core*core));
        Vector3D<double> a(ax, ay, az);
        p->treeAcceleration += a;
        
        double idt2 = (v0*v0)/(r*r + core*core);
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
