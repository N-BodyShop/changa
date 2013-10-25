#include <charm++.h>
#include "cklists.h"
#include "MultistepOrbLB.h"
#include "TopoManager.h"
#include "ParallelGravity.h"
#include "Vector3D.h"
#include <queue>

extern CProxy_TreePiece treeProxy;
CkpvExtern(int, _lb_obj_index);
using namespace std;
//#define ORB3DLB_NOTOPO_DEBUG CkPrintf

CreateLBFunc_Def(MultistepOrbLB, "Works best with multistepped runs; uses Orb3D_notopo for larger steps, greedy otherwise");

void MultistepOrbLB::init() {
  lbname = "MultistepOrbLB";
  if (CkpvAccess(_lb_obj_index) == -1) {
    CkpvAccess(_lb_obj_index) = LBRegisterObjUserData(sizeof(TaggedVector3D));
  }
}

MultistepOrbLB::MultistepOrbLB(const CkLBOptions &opt): OrbLB(opt, false) {
  init();

  if (CkMyPe() == 0){
    CkPrintf("[%d] MultistepOrbLB created\n",CkMyPe());
  }
}

bool MultistepOrbLB::QueryBalanceNow(int step){
 if(CkMyPe() == 0) CkPrintf("Orb3dOrbLB: Step %d\n", step);
  if(step == 0) return false;
  return true;

}

/// Threshold between ORB (large) and none (small) as fraction of
/// active particles
#define LARGE_PHASE_THRESHOLD 0.0001

/// @brief Implement load balancing: store loads and decide between
/// ORB and none.
void MultistepOrbLB::work(BaseLB::LDStats* stats)
{
#if CMK_LBDB_ON
  // find active objects - mark the inactive ones as non-migratable
  int count;

  int numActiveObjects = 0;
  int numInactiveObjects = 0;

  // to calculate ratio of active particles in phase
  int64_t numActiveParticles = 0;
  int64_t totalNumParticles = 0;
  
  for(int i = 0; i < stats->n_objs; i++){
    stats->to_proc[i] = stats->from_proc[i];
  }

  for(int i = 0; i < stats->n_objs; i++){
    LDObjData &odata = stats->objData[i];
    TaggedVector3D* udata = (TaggedVector3D *)odata.getUserData(CkpvAccess(_lb_obj_index));
    numActiveParticles += udata->numActiveParticles;
    totalNumParticles += udata->myNumParticles;

    if(udata->numActiveParticles == 0){
      numInactiveObjects++;
      if(stats->objData[i].migratable){
        stats->objData[i].migratable = 0;
#ifdef MCLBMSV
        CkPrintf("marking object %d non-migratable (inactive)\n", i);
#endif
        stats->n_migrateobjs--;
      }
    }
    else{
      numActiveObjects++;
    }
  }
  CkPrintf("numActiveObjects: %d, numInactiveObjects: %d\n", numActiveObjects,
	   numInactiveObjects);
  if(numInactiveObjects < 1.0*numActiveObjects) {
    // insignificant number of inactive objects; migrate them anyway
    for(int i = 0; i < stats->n_objs; i++){
      if(!stats->objData[i].migratable){
        stats->objData[i].migratable = 1;
        stats->n_migrateobjs++;
        numActiveObjects++;
        numInactiveObjects--;
      }
    }
    CkPrintf("Migrating all: numActiveObjects: %d, numInactiveObjects: %d\n", numActiveObjects, numInactiveObjects);
  }

  /*
  CkPrintf("**********************************************\n");
  CkPrintf("Object load predictions phase %d\n", phase);
  CkPrintf("**********************************************\n");
  for(int i = 0; i < stats->n_objs; i++){
      int tp = tpCentroids[i].tp;
      int lb = tpCentroids[i].tag;
    CkPrintf("tp %d load %f\n",tp,stats->objData[lb].wallTime);
  }
  CkPrintf("**********************************************\n");
  CkPrintf("Done object load predictions phase %d\n", prevPhase);
  CkPrintf("**********************************************\n");
  */

  // select processors
#ifdef MCLBMSV
  //printData(*stats, phase, NULL);
  CkPrintf("making active processor list\n");
#endif
  count = stats->count;

  // let the strategy take over on this modified instrumented data and processor information
  if((float)numActiveParticles/totalNumParticles > LARGE_PHASE_THRESHOLD){
    if (_lb_args.debug()>=2) {
      CkPrintf("******** BIG STEP *********!\n");
    }
    OrbLB::work(stats);
  }     // end if phase == 0
#endif //CMK_LDB_ON

}

void MultistepOrbLB::pup(PUP::er &p){
  CentralLB::pup(p);
}

#include "MultistepOrbLB.def.h"
