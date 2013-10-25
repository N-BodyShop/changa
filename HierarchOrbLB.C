/**
 * Author: gplkrsh2@illinois.edu (Harshitha Menon)
 * Hierarchical load balancer for ChaNGa. It has 3 level of tree where the top
 * level applies either RefineLB or uses VectorStrategy. At level 1, it applies
 * OrbLB_notopo.
*/

#include "HierarchOrbLB.h"
#include "LBDBManager.h"

#include "GreedyLB.h"
#include "RefineLB.h"
#include "Orb3dLB_notopo.h"

#define  DEBUGF(x)      // CmiPrintf x;

CreateLBFunc_Def(HierarchOrbLB, "Hybrid load balancer")

CkpvExtern(int, _lb_obj_index);
extern BaseLB* AllocateOrb3dLB_notopo();

void HierarchOrbLB::init() {
  lbname = (char *)"HierarchOrbLB";
  thisProxy = CProxy_HierarchOrbLB(thisgroup);
}

HierarchOrbLB::HierarchOrbLB(const CkLBOptions &opt): HybridBaseLB(opt) {
#if CMK_LBDB_ON
  init();

  // decide which load balancer to call
  // IMPORTANT: currently, the greedy LB must allow objects that
  // are not from existing processors.
  refinelb = (CentralLB *)AllocateRefineLB();
  orblb = (CentralLB *)AllocateOrb3dLB_notopo();

  initTree();
#endif
}

HierarchOrbLB::~HierarchOrbLB() {
  delete orblb;
  delete refinelb;
}

void HierarchOrbLB::GetObjsToMigrate(int toPe, double load, LDStats *stats,
    int atlevel, CkVec<LDCommData>& comms, CkVec<LDObjData>& objs) {
  for (int obj=stats->n_objs-1; obj>=0; obj--) {
    LDObjData &objData = stats->objData[obj];
    TaggedVector3D* udata = (TaggedVector3D *)objData.getUserData(CkpvAccess(_lb_obj_index));

    if (!objData.migratable || udata->myNumParticles <= 0) continue;
    if (objData.wallTime <= load) {
      if (_lb_args.debug()>2)
        CkPrintf("[%d] send obj: %d to PE %d (load: %f).\n", CkMyPe(), obj, toPe, objData.wallTime);
      objs.push_back(objData);
      // send comm data
      collectCommData(obj, comms, atlevel);
      load -= objData.wallTime;
      CreateMigrationOutObjs(atlevel, stats, obj);
      stats->removeObject(obj);
      if (load <= 0.0) break;
    }
  }
}

void HierarchOrbLB::work(LDStats* stats) {
#if CMK_LBDB_ON
  LevelData *lData = levelData[currentLevel];

  if (currentLevel == 1) {
    orblb->work(stats);
  }
  else
    refinelb->work(stats);
#endif
}

#include "HierarchOrbLB.def.h"
