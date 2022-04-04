#ifdef VORONOI
#include "ParallelGravity.h"
#include "TreeWalk.h"
#include "DataManager.h"
#include "Opt.h"
#include "smooth.h"
#include "Space.h"
#include "VoronoiSmoothParams.h"
#include "Neighbor.h"
#include "HydroParticle.h"
#include <vector>
#include "UserProblemGenerator.h"

int VoronoiSmoothParams::isSmoothActive(GravityParticle *p)
{
  if( !TYPETest( p, TYPE_GAS)) return 0;
  // for INITIALIZE bActiveOnly is set to VSP_EVERYTHING
  bool activeTest = false, test = false;
  if( hydroStep == BUILD_MESH || hydroStep == MAP_A_TO_B ) {
    activeTest = p->hydroParticle().isActive( activeRung);
  }
  else if (hydroStep == HALF_STEP || hydroStep == FULL_STEP) 
    activeTest = p->hydroParticle().isActiveHalf( activeRung); 
  else 
    bActiveOnly = VSP_EVERYTHING;

  switch( bActiveOnly) {
    case VSP_ACTIVEONLY :
      if( !activeTest) return 0;		// not active
      if( bReSearch && activeTest && p->hydroParticle().isValid()) // no need to research if already valid
        return 0;
      return activeTest; 
      break;
    case VSP_NEIGHBORONLY:
      if( !activeTest && TYPETest(p, TYPE_NbrOfACTIVE)) 
        return bReSearch && p->hydroParticle().isValid() ? 0 : 1; // no need to research if already valid
      return 0;
      break; 
    case VSP_ACTIVE_AND_NEIGHBOR:
      test = activeTest || TYPETest(p, TYPE_NbrOfACTIVE);
      if( bReSearch && test && p->hydroParticle().isValid()) return 0;
      return test;
      break;
    default :
      return !(bReSearch &&  p->hydroParticle().isValid()); // no need to research if already valid
      break;
  }

}

void VoronoiSmoothParams::initSmoothParticle(GravityParticle *p)
{}

void VoronoiSmoothParams::initTreeParticle(GravityParticle *p)
{}

void VoronoiSmoothParams::initSmoothCache(GravityParticle *p)
{}

void VoronoiSmoothParams::combSmoothCache(GravityParticle        *p1,
                                          ExternalSmoothParticle *p2)
{
  p1->iType |= p2->iType;
}

/* Gather only version */
void VoronoiSmoothParams::fcnSmooth(GravityParticle *p, int nSmooth,
                                    pqSmoothNode *nnList)
{

  double rx = HydroUtils::instance->xPeriod;
  double ry = HydroUtils::instance->yPeriod;
  double rz = HydroUtils::instance->zPeriod;

  int firstOrder = HydroUtils::instance->firstOrder;

  rx*=0.1; ry*=0.1;
//  if( !(HydroUtils::instance->xPeriod > HydroUtils::instance->yPeriod) && !(HydroUtils::instance->xPeriod > HydroUtils::instance->zPeriod)) {
    rz*=0.1; // 3-d case else 2-d
    rz = max(max(rx,ry),rz); rx=rz; ry=rz;
 // }

  VoronoiCompute compute(0., 0., 0., rx, ry, rz, p->fBall);

  std::unordered_map<int, HydroParticle*> map;
  std::vector<HydroParticle*> nearestNeighbors;
  nearestNeighbors.reserve( nSmooth);
  
  for (int i = 0; i < nSmooth; ++i) {
    //
    // The smooth compute gives dx in terms of the particle to test vs the local particle is dx is formally reversed.
    //
    double dx = -nnList[i].dx.x;
    double dy = -nnList[i].dx.y;
    double dz = -nnList[i].dx.z;
 
    // skip if it is the same particle
    if( p->iOrder == nnList[i].p->iOrder) continue;

    if( bActiveOnly == VSP_ACTIVEONLY) {
      bool activeTestQ = false;
      bool activeTestP = false;
      GravityParticle *q = nnList[i].p;
      if( hydroStep == BUILD_MESH || hydroStep == MAP_A_TO_B ) {
        activeTestQ = q->hydroParticle().isActive( activeRung);
        activeTestP = p->hydroParticle().isActive( activeRung);
      }
      else if (hydroStep == HALF_STEP || hydroStep == FULL_STEP) {
        activeTestQ = q->hydroParticle().isActiveHalf( activeRung);
        activeTestP = p->hydroParticle().isActiveHalf( activeRung);
      } 
      else break;
      
      if( activeTestP && !activeTestQ)
        TYPESet(q,TYPE_NbrOfACTIVE); // important for multistepping
      else TYPEReset(q,TYPE_NbrOfACTIVE);
    }

    compute.addNeighbor(nnList[i].p->iOrder, dx, dy, dz);
    std::pair<int, HydroParticle*> pair(nnList[i].p->iOrder,
                                        &nnList[i].p->hydroParticle());
    map.insert(pair);
    nearestNeighbors.push_back( &nnList[i].p->hydroParticle());

     if( nnList[i].p->iOrder != nnList[i].p->hydroParticle().iOrder) {
       printf("Error in iOrder consistency, %ld %ld %ld %d\n", nnList[i].p->iOrder,
	     nnList[i].p->hydroParticle().iOrder, p->iOrder, hydroStep);
     }
  }
  
  compute.compute(p->hydroParticle());
  
  // first step or second step
  if( hydroStep == INITIALIZE || hydroStep == BUILD_MESH || hydroStep == MAP_A_TO_B) {
    if ( hydroStep == INITIALIZE) {
      // compute primitives
      double volume        = p->hydroParticle().getVolume();
      double mass          = p->mass;
      double u             = p->u();
      Vector3D<double> v   = p->velocity;
      Vector3D<double> pos = p->position;

      double rho = p-> fDensity;
      //double rho = mass/volume;
      //
      p->hydroParticle().iOrder = p->iOrder;
      p->hydroParticle().updateParticlePosition(pos.x, pos.y, pos.z);
      p->hydroParticle().updateParticleVelocity(v.x, v.y, v.z);
#ifdef COOLING_SROEOS
      p->hydroParticle().setPrimitives(iye, p->fMetals()); // this must be set before
#endif
#ifdef COOLING_MESA
      p->hydroParticle().setPrimitives(iX, p->fMetals()); // this must be set before
#endif

      if( bProblemGenerator) {
        dVector position(pos.x, pos.y, pos.z);
        varVector prim = 0.;
        UserProblemGenerator::instance->applyGenerator( position, prim);
        //p->hydroParticle().updateParticleVelocity(prim[ivx], prim[ivy], prim[ivz]);
    
        p->hydroParticle().setPrimitives( prim);
      } else if (bVoronoiICs) {
        p->hydroParticle().finishPrimitives();
      } else {
#ifdef RADIATION
        dVector xHat(1., 0., 0.);
        for( int iFreq = 0; iFreq < NFREQVARS; iFreq++) { 
          for( int iAng = 0; iAng < NANGVARS; iAng++) {
            int ivar; ivar = IRADVARSTART + iAng + iFreq*NANGVARS;
            if( dotProduct( HydroUtils::instance->radNormals[iAng], xHat) >= HydroUtils::instance->initRadMu) 
	            p->hydroParticle().setPrimitives(ivar, HydroUtils::instance->initRadIntensity); 
// this must be set before
          }
        }

#endif
        p->hydroParticle().setPrimitives(rho, v.x, v.y, v.z, u);
      }
      p->hydroParticle().setupQ();
      //PHIL: CHECK IF THIS IS SUPPOSE TO BE HERE FOR MHD
      p->hydroParticle().setTimeStep( nearestNeighbors);
    } else if( hydroStep == BUILD_MESH) {
      p->hydroParticle().setTimeStep( nearestNeighbors, true);
      if( p->hydroParticle().isValid() && activeRung <= p->hydroParticle().getRung()) {
        if( HydroUtils::instance->useVanLeer) {
          p->hydroParticle().solveRiemann( activeRung, map, true);
        }
        else { 
          p->hydroParticle().computeGradients( map, true);
          p->hydroParticle().predict( activeRung);
        }
      }
    }
#ifdef MHD
    else if( hydroStep == MAP_A_TO_B) {
      if( p->hydroParticle().isValid()) {
        p->hydroParticle().computeGradients( map, true);
        p->hydroParticle().mapAtoB();
      }
    }
#endif //MHD
  }
  else if( hydroStep == HALF_STEP) {
    if( p->hydroParticle().isValid() && activeRung == p->hydroParticle().getHalfRung()) {
      p->hydroParticle().computeGradients( map, false, true);
    }
  }
  else if( hydroStep == FULL_STEP) {
    if( p->hydroParticle().isValid()) 
      p->hydroParticle().solveRiemann( activeRung, map, false);
  }
  else if( hydroStep == CHECK_SOLUTION) {
    p->hydroParticle().checkSolution( map);
  }

}
#endif //VORONOI
