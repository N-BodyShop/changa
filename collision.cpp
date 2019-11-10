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
    bDelEjected = 0;
    prmAddParam(prm, "bDelEjected", paramBool, &bDelEjected,
        sizeof(int), "bDelEjected", "<Delete particles which travel a distance dDelDist\
                                      from the origin> = 0");
    dDelDist = 1000.0;
    prmAddParam(prm, "dDelDist", paramDouble, &dDelDist,
        sizeof(double), "dDelDist", "<Distance away from the origin before particles are\
                                      deleted> = 1000.0");

    bCollStep = 0;
    prmAddParam(prm, "bCollStep", paramBool, &bCollStep,
        sizeof(int), "bCollStep", "<Place particles on a near-collision trajectory\
                                    on a high rung> = 0");

    nSmoothCollision = 64;
    prmAddParam(prm, "nSmoothCollision", paramInt, &nSmoothCollision,
        sizeof(int), "nSmoothCollision", "<number of particles to do collision\
                                           search over> = 64");

    iCollStepRung = 7;
    prmAddParam(prm, "iCollStepRung", paramInt, &iCollStepRung,
        sizeof(int), "iCollStepRung", "<Rung to place nearly-colliding particles on> = 7");

    dCollStepFac = 2.5;
    prmAddParam(prm, "dCollStepFac", paramDouble, &dCollStepFac,
        sizeof(double), "dCollStepFac", "<Factor by which particle radius is inflated\
                                       when searching for near-collisions> = 2.5");
   
    bAllowMergers = 0;
    prmAddParam(prm, "bAllowMergers", paramBool, &bAllowMergers,
        sizeof(int), "bAllowMergers", "<Particles merge on collision if their\
                                        approach speed is small enough> = 0");
    bWall = 0;
    prmAddParam(prm, "bWall", paramBool, &bWall,
        sizeof(int), "bWall", "<Particles bounce off of a plane pointed in\
                                 the z direction> = 0");

    dWallPos = 0.0;
    prmAddParam(prm, "dWallPos", paramDouble, &dWallPos,
        sizeof(double), "dWallPos", "<Z coordinate of wall> = 0");

    bPerfectAcc = 0;
    prmAddParam(prm, "bPerfectAcc", paramBool, &bPerfectAcc,
        sizeof(int), "bPerfectAcc", "<All collisions result in a merger> = 0");

    dBallFac = 2.0;
    prmAddParam(prm, "dBallFac", paramDouble, &dBallFac,
        sizeof(double), "dBallFac", "<Scale factor for collision search radius> = 2.0");

    dEpsN = 1.0;
    prmAddParam(prm, "dEpsN", paramDouble, &dEpsN,
        sizeof(double), "dEpsN", "<Normal coefficient of restituion for bouncing collisions> = 1.0");

    dEpsT = 1.0;
    prmAddParam(prm, "dEpsT", paramDouble, &dEpsT,
        sizeof(double), "dEpsT", "<Tangential coefficient of restituion for bouncing collisions> = 1.0");

    bSkipP0 = 0;
    prmAddParam(prm, "bSkipP0", paramBool, &bSkipP0,
        sizeof(int), "bSkipP0", "<Don't do collision check for first particle> = 0");

    }

void Collision::CheckParams(PRM prm, struct parameters &param)
{
    if (bPerfectAcc && !bAllowMergers)
        CkAbort("Perfect accretion will not occur if mergers are disabled\n");
    if (bWall && bAllowMergers)
        CkAbort("Mergers cannot occur while bWall is enabled\n");
#ifndef COLLISION
    if (param.bCollision)
        CkAbort("ChaNGa must be compiled with the COLLISION flag in order to use collision detection\n");
#endif
    }

/**
 * @brief Remove particles that are too far from the origin
 *
 * Because we are not using gravitational softening when resolving collisions,
 * scattering events can occasionally be very strong and throw particles to
 * very large distances, which causes ChaNGa to hang. Since these particles
 * are usually no longer important, we just delete them once they get far enough
 * away.
 *
 * @param dDelDist Distance from the origin beyond which a particle is deleted
 */
void TreePiece::delEjected(double dDelDist, const CkCallback& cb)
{
    for (unsigned int i=1; i <= myNumParticles; i++) {
        GravityParticle *p = &myParticles[i];
        double dist = p->position.length();
        if (dist > dDelDist) {
            CkPrintf("%d ejected, deleting\n", p->iOrder);
            deleteParticle(p);
            }
        }
    contribute(cb);
    }

/**
 * @brief Predict near approaches between particles
 *
 * A near-collision search is done using SmoothParams, where the 'rung' 
 * field of all particles which will undergo a near collision in the next
 * time step is set.
 *
 * @param dTime The current time in the simulation (in simulation units)
 * @param dDelta The size of the next timestep (in simulation units)
 * @param activeRung The currently active rung
 */
void Main::doNearCollisions(double dTime, double dDelta, int activeRung)
{
    CollisionSmoothParams pCS(TYPE_DARK, activeRung, dTime, dDelta, 
       param.collision.bWall, param.collision.dWallPos,
       param.collision.bAllowMergers, 1, param.collision);
    double dfBall2OverSoft2 = 4.0*param.dhMinOverSoft*param.dhMinOverSoft;
    treeProxy.startSmooth(&pCS, 0, param.collision.nSmoothCollision,
          dfBall2OverSoft2, CkCallbackResumeThread());

    // Need to ensure that both near-colliding particles end up on the collision rung
    CkReductionMsg *msgChk;
    treeProxy.getNearCollPartners(CkCallbackResumeThread((void*&)msgChk));
    int numIords = msgChk->getSize()/sizeof(int);
    int *data = (int *)msgChk->getData();

    for (int i=0; i < numIords; i++) {
        treeProxy.placeOnCollRung(data[i], param.collision.iCollStepRung, CkCallbackResumeThread());
        }
    }

/**
 * @brief Compile a list of every collision partner that has been identified
 *
 * The 'collision smooth' function only allows you to set properties of one
 * of the two colliding particles. It is possible that particle a will
 * identify an impending collision with particle b, but b will not identify
 * a collision with a. This can cause problems with collision stepping, as
 * only one of the two near-colliders will end up on the proper rung. This
 * function constructs a list of the iorders of all of the collision partners
 * so that we can go back and properly set their rung.
 */
void TreePiece::getNearCollPartners(const CkCallback& cb) {
    std::vector<int> collPartners;
    for (unsigned int i=1; i <= myNumParticles; i++) {
        GravityParticle *p = &myParticles[i];
        if (p->iOrderCol >= 0) {
            collPartners.push_back(p->iOrderCol);
            }
        }

    int* a = &collPartners[0];
    contribute(sizeof(int)*collPartners.size(), a, CkReduction::concat, cb);
    }

/**
 * @brief Detect and resolve collisions between particles. 
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
 * @param activeRung The currently active rung
 */
void Main::doCollisions(double dTime, double dDelta, int activeRung)
{
    int bHasCollision;
    int nColl = 0;
    
    do {
        bHasCollision = 0;

        // Use the smooth framework to search for imminent collisions
        // This sets the 'dtCol' and 'iOrderCol' fields for the particles
        CollisionSmoothParams pCS(TYPE_DARK, activeRung, dTime, dDelta, 
           param.collision.bWall, param.collision.dWallPos,
           param.collision.bAllowMergers, 0, param.collision);
        double dfBall2OverSoft2 = 4.0*param.dhMinOverSoft*param.dhMinOverSoft;
        treeProxy.startSmooth(&pCS, 0, param.collision.nSmoothCollision,
              dfBall2OverSoft2, CkCallbackResumeThread());

        // Once 'dtCol' and 'iOrderCol' are set, we need to determine which
        // collision is going to happen the soonest
        CkReductionMsg *msgChk;
        treeProxy.getCollInfo(CkCallbackResumeThread((void*&)msgChk));
        ColliderInfo *c = (ColliderInfo *)msgChk->getData();

        // If only one particle detected the imminent collision (due to different search
        // radii or velocities), we need to go back and ask for the second collider
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
        delete msgChk;

        // If the collision is going to happen during this timestep, resolve it now
        if (c[0].dtCol <= dDelta) {
            bHasCollision = 1;

            // Collision with wall
            if (c[0].iOrderCol == -2 && param.collision.bWall) {
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

                treeProxy.resolveCollision(param.collision, c[0], c[1], dDelta,
                                           dTime, CkCallbackResumeThread());
                }
            }

        // Velocities and positions may have changed, keep looking for collisions
        } while (bHasCollision);

        // Clean up any merged particles
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
        if (p->dtCol < dtMin && !TYPETest(p, TYPE_DELETED)) {
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
        if (p->dtCol == dtMin && p->dtCol < DBL_MAX && !TYPETest(p, TYPE_DELETED)) {
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
void TreePiece::resolveWallCollision(Collision coll, const ColliderInfo &c1, 
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
void TreePiece::resolveCollision(Collision coll, const ColliderInfo &c1,
                                 const ColliderInfo &c2, double baseStep,
                                 double timeNow, const CkCallback& cb) {
    GravityParticle *p;
    int iCollResult;
    // Colliders need to be handled individually because only one
    // of the two could be on this tree piece
    for (unsigned int i=1; i <= myNumParticles; i++) {
        p = &myParticles[i];
        if (p->iOrder == c1.iOrder) {
            iCollResult = coll.doCollision(p, c2);
            if (iCollResult == MERGE) {
                if (c1.mass >= c2.mass) {
                    CkPrintf("Merge %d into %d\n", p->iOrder, c2.iOrder);
                    coll.setMergerRung(p, c2, c1, baseStep, timeNow);
                }
                else {
                    CkPrintf("Delete %d\n", p->iOrder);
                    deleteParticle(p);
                }
            }
            break;
            }
        }

    for (unsigned int i=1; i <= myNumParticles; i++) {
        p = &myParticles[i];
        if (p->iOrder == c2.iOrder) {
            iCollResult = coll.doCollision(p, c1);
            if (iCollResult == MERGE) {
                if (c2.mass > c1.mass) {
                    CkPrintf("%g %g\n", p->mass, c1.mass);
                    CkPrintf("Merge %d into %d\n", p->iOrder, c1.iOrder);
                    coll.setMergerRung(p, c1, c2, baseStep, timeNow);
                }
                else {
                    CkPrintf("Delete %d\n", p->iOrder);
                    deleteParticle(p);
                }
            }
            break;
            }
        }

    contribute(cb);
    }

/**
 * @brief Undo the kick imparted to all particles
 *
 * This method is meant to be called directly after a opening kick on rung 0. 
 * When using collision stepping, all particles will receive an opening kick
 * on rung 0, but particles eligible for collision stepping will now be on a
 * higher rung and need a smaller kick.
 *
 * @param dDeltaBase The timestep size that was used to calculate the previous kick
 */
void TreePiece::unKickCollStep(int iKickRung, double dDeltaBase, const CkCallback& cb)
{
    for(unsigned int i = 1; i <= myNumParticles; ++i) {
      GravityParticle *p = &myParticles[i];
      if (p->rung >= iKickRung) p->velocity -= dDeltaBase*p->treeAcceleration;
      }
    contribute(cb);
    }

void TreePiece::placeOnCollRung(int iOrder, int collStepRung, const CkCallback& cb) {
    for (unsigned int i = 1; i <= myNumParticles; ++i) {
      GravityParticle*p = &myParticles[i];
      if (p->iOrder == iOrder) p->rung = collStepRung;
      }
    contribute(cb);
    }

/**
 * @brief Place all particles back on rung 0
 *
 * This function is used at the end a collision stepping big step
 */
void TreePiece::resetRungs(const CkCallback& cb)
{
    for(unsigned int i = 1; i <= myNumParticles; ++i) {
      GravityParticle *p = &myParticles[i];
      if(p->rung > 0) p->rung = 0;
      }

    contribute(cb);
}

/**
 * @brief Determine whether there are any particles in the simulation that
 * will undergo a near miss in the next step
 *
 * The way to tell if a particle is undergoing a near miss is to see if it sits
 * on the collision stepping rung.
 *
 * @param collStepRung The rung that particles are moved to after finding a near miss
  */
void TreePiece::getNeedCollStep(int collStepRung, const CkCallback& cb)
{
    int foundNearMiss = 0;
    for(unsigned int i = 1; i <= myNumParticles; ++i) {
      GravityParticle *p = &myParticles[i];
      if (p->rung == collStepRung) {
          foundNearMiss++;
          }
      }

    contribute(sizeof(int), &foundNearMiss, CkReduction::sum_int, cb);
    }


/**
 * @brief Check for a merger between two colliders
 *
 * A merger will occur if the relative velocity between particles
 * is < vEsc and both particles post collision angular speed is
 * < wMax.
 *
 * @param c1 The first collider participating in the event
 * @param c2 The second collider participating in the event
 *
 * @returns 0 if no merger, -1 or 1 if c1 or c2 are to be consumed
 */
int Collision::checkMerger(const ColliderInfo &c1, const ColliderInfo &c2)
{

    // Calculate the post collision spin and velocity of the merged particle
    Vector3D<double> posNew, vNew, wNew, aNew, pAdjust;
    double radNew;
    // Advance particle positions to moment of collision
    pAdjust = c1.position + c1.velocity*c1.dtCol;
    mergeCalc(c1.radius, c1.mass, pAdjust, c1.velocity, c1.acceleration,
              c1.w, &posNew, &vNew, &wNew, &aNew, &radNew, c2);

    double m = c2.mass;
    double M = c1.mass;
    double Mtot = M + m ;



    double r = c2.radius;
    double R = c1.radius;
    Vector3D<double> vRel = (c1.velocity - c2.velocity);
    double V_i =vRel.length();
    Vector3D<double> pRel = (c1.position - c2.position);
    double cos_theta = costheta(vRel, pRel);
    double sintheta = sqrt(1-pow(cos_theta,2));

    double l = (R+r) *(1- sintheta);
    double alpha = (3*r*(pow(l,2))-pow(l,3))/(4*pow(r,3));
    double Mint = alpha * c1.mass;
    double vEsc = sqrt(2.*Mint/(c1.radius + c2.radius));
    double wMax = sqrt(Mtot/(radNew*radNew*radNew));

    double b = sintheta ;
    double bcrit = (R) / (R+r) ;
    if (b > bcrit) {
        CkPrintf ("grazing collision\n");
        } else  {
            CkPrintf ("nongrazing collision\n") ;
        }
    double density = 5.028e-28;
    double RC1 = (Mtot)/(density);

    double cs = 1;
    double dDenFac = 4./3.* M_PI;
    double row1 = Mtot/ (4/3* dDenFac* pow(r,3));
    double QsRD = (cs) * (4/5)* dDenFac * row1 * (pow(RC1,2));
    double Vsy = pow(32*dDenFac,.5)/5 * pow(row1,.5) * RC1;

    //we defined alpha in line 700 upcoming is step c(idk what Mtarget is)
    double Mp = pow(4*r,3);
    double Mtarget = M ;
    double mu = (M + m) / (M*m); 


    // d 
    double Vs = (QsRD) * (Vsy);

    // e for equation 15  i dont understand it's exponet and how to do square root ex on eq16
    double Vbars = sqrt(1/alpha * pow(Vs,2));  
    double mubar = 1; //set by user
    double QsRDa = 1/ alpha * QsRD * pow(Vbars/Vs,2-(3*mubar)) ;
    double Vsa = sqrt((2*(QsRDa)* Mtot)/ mu );

    //step 5 
    // find qr

    double Mlr = M;
    double mtot = M + m;
    double QR = ((QsRD - 1) * ( Mlr/mtot - .5)) / -.5 ;

    // finding Verosion = Vi = V_erosion changed name since we had defined it earlier..

    double V_erosion = sqrt((mtot * QR) / (.5 * mu));
     



    //step 6 (need to check code)
    if (b > bcrit || vEsc < V_i < V_erosion) {
            printf("Hit-and-run collision\n");
            
    }

    if (vRel.length() > vEsc || wNew.length() > wMax) {
        if (!bPerfectAcc) {
            CkPrintf("Merger rejected\n");
            return 0;
            }
        }

// step 7
        double QRn = (((Mlr/.1 *Mtot)-.5)/-.5)*QsRDa;
        double Vsupercat = sqrt(QRn*2 /mtot);

    if (V_i > V_erosion &&V_i<Vsupercat){
        CkPrintf("disruption");
    }
    //Ckprintf("%g",)

    double mass_Mlr = QR/QsRDa;
    double beta = 2.85;
    double Mslr = Mtot*(3-beta)*1-(Mlr/Mtot)/(2*beta);
        //equation 35

    double C = beta*pow((3-beta)*(Mtot-mass_Mlr)/ 4/3 *M_PI*row1*2*beta,beta/3);


    

double gammad= Mint/Mp;
//equation23

double QRDDS = pow(.25*pow(gammad+1,2)/gammad,2/((3*mubar)-1));
 

double deltaV= 2e-3;
double A = -.3*Mlr/Mtot+.3;
double S = (pow(10,A))/ (log(10)*deltaV*((Mtot-Mlr)/Mtot));
Vector3D<double> Vlr = M*c1.velocity + m*c2.velocity;
double norm =((log10(Mslr*deltaV*S*log(10))-A)/-S)-(Vlr.length());


double Vslr=((log10(Mslr*(deltaV/10)*S*log(10))-A)/-S)-norm;
CkPrintf("%g \n ",Vslr);



    // Mark the less massive particle to be consumed and deleted
    if (c1.mass < c2.mass) return -1;
    else return 1;
    
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
void Collision::setMergerRung(GravityParticle *p, const ColliderInfo &c, const ColliderInfo &cMerge,
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
    // TODO: update spin field
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
 *
 * @return The CollType of the collision
 */
int Collision::doCollision(GravityParticle *p, const ColliderInfo &c)
{
    // Determine the collision type
    int iCollType = MERGE;
    // Calculate the post collision spin and velocity of the merged particle
    Vector3D<double> posNew, vNew, wNew, aNew, pAdjust;
    double radNew;
    // Advance particle positions to moment of collision
    pAdjust = p->position + c.velocity*c.dtCol;
    mergeCalc(p->soft*2, p->mass, pAdjust, p->velocity, p->treeAcceleration,
              p->w, &posNew, &vNew, &wNew, &aNew, &radNew, c);
    double Mtot = p->mass + c.mass;
    double vEsc = sqrt(2.*Mtot/(p->soft*2 + c.radius));
    double wMax = sqrt(Mtot/(radNew*radNew*radNew));

    double vRel = (p->velocity - c.velocity).length();
    if (vRel > vEsc || wNew.length() > wMax) {
        if (!bPerfectAcc) {
            CkPrintf("Merger rejected\n");
            iCollType = BOUNCE;
            }
        }

    // If the particles are different sizes, p might not have its dtCol
    // field set. Fix this before going any further.
    p->dtCol = c.dtCol;

    // Advance particle positions to moment of collision
    Vector3D<double> pAdv;
    pAdv = p->position + p->velocity*p->dtCol;

    if (iCollType == MERGE) {
        mergeCalc(p->soft*2., p->mass, pAdv, p->velocity, p->treeAcceleration,
                  p->w, &posNew, &vNew, &wNew, &aNew, &radNew, c);
        CkPrintf("Merger info:\niorder1 iorder2 m1 m2 r1 r2 x1x x1y x1z x2x x2y\
                  x2z xNewx xNewy xNewz v1x v1y v1z v2x v2y v2z vNewx vNewy vNewz\
                   w1x w1y w1z w2x w2y w2z wNewx wNewy wNewz\n");
        CkPrintf("%d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\
                  %g %g %g %g %g %g %g %g %g %g %g %g\n",
                p->iOrder, c.iOrder, p->mass, c.mass, p->soft*2, c.radius,
                p->position.x, p->position.y, p->position.z,
                c.position.x, c.position.y, c.position.z,
                posNew.x, posNew.y, posNew.z,
                p->velocity.x, p->velocity.y, p->velocity.z,
                c.velocity.x, c.velocity.x, c.velocity.z,
                vNew.x, vNew.y, vNew.z,
                p->w.x, p->w.y, p->w.z,
                c.w.x, c.w.y, c.w.z,
                wNew.x, wNew.y, wNew.z);
        p->position = posNew;
        p->treeAcceleration = aNew;
        p->soft = radNew/2.;
        p->mass += c.mass;
                  
        }
    else {
        CkPrintf("Bounce info:\niorder1 iorder2 m1 m2 r1 r2 x1x x1y x1z x2x x2y\
                  x2z xNewx xNewy xNewz v1x v1y v1z v2x v2y v2z vNewx vNewy vNewz\
                   w1x w1y w1z w2x w2y w2z wNewx wNewy wNewz\n");
        CkPrintf("%d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\
                  %g %g %g %g %g %g %g %g %g %g %g %g\n",
                p->iOrder, c.iOrder, p->mass, c.mass, p->soft*2, c.radius,
                p->position.x, p->position.y, p->position.z,
                c.position.x, c.position.y, c.position.z,
                pAdv.x, pAdv.y, pAdv.z,
                p->velocity.x, p->velocity.y, p->velocity.z,
                c.velocity.x, c.velocity.x, c.velocity.z,
                vNew.x, vNew.y, vNew.z,
                p->w.x, p->w.y, p->w.z,
                c.w.x, c.w.y, c.w.z,
                wNew.x, wNew.y, wNew.z);

       
        bounceCalc(p->soft*2., p->mass, pAdv, p->velocity, p->w, &vNew, &wNew, c);
        makeFragment();
        }

    p->velocity = vNew;
    p->w = wNew;

    p->dtCol = DBL_MAX;

    return iCollType;
    }
void TreePiece::makeFragments(Collision coll, const CkCallback& cb)
{ GravityParticle *frag = coll.makeFragment();
    newParticle (frag);
    delete frag;
    int counts[2] ;
    counts[0] = 1; //number of particles formed
    counts[1] = 0; //number of particles deleted

    contribute(2*sizeof(int), counts, CkReduction::sum_int,cb);
    }
GravityParticle* Collision:: makeFragment()
{

    GravityParticle *frag = new GravityParticle();
    frag->mass = 1e-50 ;
    frag->position.x = 1.5; 
    frag->position.y = 1.5;
    frag->position.z = 0;
    frag->velocity.x = 0; 
    frag->velocity.y = 0;
    frag->velocity.z = 0;
    frag->soft = .1/2;
TYPESet(frag, TYPE_DARK);

    CkPrintf("make fragment\n");
    
return frag ;
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
                          Vector3D<double> vel, Vector3D<double> acc,
                          Vector3D<double> w, Vector3D<double> *posNew,
                          Vector3D<double> *velNew, Vector3D<double> *wNew,
                          Vector3D<double> *aNew,  double *radNew, const ColliderInfo &c)
{
    double r1 = r;
    double r2 = c.radius;
    double m1 = m;
    double m2 = c.mass;
    double M = m1 + m2;

    double dDenFac = 4./3.*M_PI;
    double rho1 = m1/(dDenFac*pow(r1, 3.));
    double rho2 = m2/(dDenFac*pow(r2, 3.));

    Vector3D<double> cPosNew = c.position + c.velocity*c.dtCol;

    *radNew = pow(M/(dDenFac*rho1), 1./3.); // Conserves density

    double i1 = 0.4*m1*r1*r1;
    double i2 = 0.4*m2*r2*r2;
    double i = 0.4*M*(* radNew)*(* radNew);

    Vector3D<double> comPos = (m1*pos + m2*cPosNew)/M;
    Vector3D<double> comVel = (m1*vel + m2*c.velocity)/M;
    Vector3D<double> comAcc = (m1*acc + m2*c.acceleration)/M;

    Vector3D<double> rc1 = pos - comPos;
    Vector3D<double> rc2 = cPosNew - comPos;
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
                           const ColliderInfo &c)
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

int CollisionSmoothParams::isSmoothActive(GravityParticle *p)
{
    if(p->rung < activeRung)
	    return 0;

    if (TYPETest(p, TYPE_DARK)) return 1;
    else return 0;
    }

void CollisionSmoothParams::initSmoothParticle(GravityParticle *p)
{
    double v = p->velocity.length();
    // Different particles can have different fBall sizes. This can potentially
    // cause problems where only 1 of 2 colliding particles detects the collision.
    p->fBall = coll.dBallFac*dTimeSub*pow(2., activeRung)*v + 4.*p->soft;
    }

void CollisionSmoothParams::initSmoothCache(GravityParticle *p1)
{
    p1->dtCol = DBL_MAX;
    p1->iOrderCol = -1;
    p1->rung = -1;
    }

void CollisionSmoothParams::combSmoothCache(GravityParticle *p1,
                              ExternalSmoothParticle *p2)
{
    if (p2->dtCol < p1->dtCol) {
        p1->dtCol = p2->dtCol;
        p1->iOrderCol = p2->iOrderCol;
        }

    if (p1->rung < p2->rung) p1->rung = p2->rung;
    }

/**
 * @brief Do a neighbor search to look for all possible collisions between
 * particles.
 *
 * This function updates the 'dtCol' and 'iOrderCol' field for all particles
 * that will undergo a collision in the next time step. If we are in the
 * collision step prediction phase, place any particles that will come close
 * together on the collision stepping rung.
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
    if (coll.bSkipP0) {
        if (p->iOrder == 0) return;
        }

    if (bWall) {
        dt = (dWallPos - p->position[2] - (p->soft*2.))/p->velocity[2];
        if (dt > 0. && dt < dTimeSub) {
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

       // Near-collision search
       if (bNearCollSearch) {
           if (p->iOrderCol == -1) {
               sr = coll.dCollStepFac*2.*(p->soft + q->soft);
               rdotv = dot(dx, vRel);
               vRel2 = vRel.lengthSquared();
               dr2 = dx.lengthSquared() - sr*sr;
               D = rdotv*rdotv - dr2*vRel2;
               if (dx.length() < sr) {
                   p->rung = coll.iCollStepRung;
                   p->iOrderCol = q->iOrder;
                   }
               else if (D > 0.) {
                   D = sqrt(D);
                   dt = (-rdotv - D)/vRel2;
                   if (dt > 0. && dt < dTimeSub) {
                       p->rung = coll.iCollStepRung;
                       p->iOrderCol = q->iOrder;
                       }
                   }
               }
       } else {
           // Collider search
           sr = (p->soft*2.) + (q->soft*2.);
           rdotv = dot(dx, vRel);
           vRel2 = vRel.lengthSquared();
           dr2 = dx.lengthSquared() - sr*sr;
           D = rdotv*rdotv - dr2*vRel2;
           if (D > 0.) {
               D = sqrt(D);
               dt = (-rdotv - D)/vRel2;

               if (dt > 0. && dt < dTimeSub && dt < p->dtCol) {
                   p->dtCol = dt;
                   p->iOrderCol = q->iOrder;
                   }
               }
           }
        }
    }
