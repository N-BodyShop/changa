/**
 * Author: gplkrsh2@illinois.edu (Harshitha Menon)
 * Hierarchical load balancer for ChaNGa. It has 3 level of tree where the top
 * level applies either RefineLB or uses VectorStrategy. At level 1, it applies
 * OrbLB_notopo.
*/

#include "Refiner.h"

#include "HierarchOrbLB.h"

#include "Orb3dLB_notopo.h"

#define  DEBUGF(x)      // CmiPrintf x;

CreateLBFunc_Def(HierarchOrbLB, "Hybrid load balancer")

CkpvExtern(int, _lb_obj_index);
extern BaseLB* AllocateOrb3dLB_notopo();

void HierarchOrbLB::init() {
  lbname = (char *)"HierarchOrbLB";
  thisProxy = CProxy_HierarchOrbLB(thisgroup);
}

HierarchOrbLB::HierarchOrbLB(const CkLBOptions &opt): CBase_HierarchOrbLB(opt) {
#if CMK_LBDB_ON
  init();

  // decide which load balancer to call
  // IMPORTANT: currently, the greedy LB must allow objects that
  // are not from existing processors.
  orblb = (CentralLB *)AllocateOrb3dLB_notopo();

  initTree();
#endif
}

HierarchOrbLB::~HierarchOrbLB() {
  delete orblb;
}

void HierarchOrbLB::GetObjsToMigrate(int toPe, double load, LDStats *stats,
    int atlevel, CkVec<LDCommData>& comms, CkVec<LDObjData>& objs) {
  for (int obj=stats->n_objs-1; obj>=0; obj--) {
    LDObjData &objData = stats->objData[obj];
    if (!objData.migratable) continue;

    TaggedVector3D* udata = (TaggedVector3D *)objData.getUserData(CkpvAccess(_lb_obj_index));

    if (udata->myNumParticles <= 0) continue;
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

// only called on leaves
CLBStatsMsg* HierarchOrbLB::AssembleStats()
{
#if CMK_LBDB_ON
  CLBStatsMsg* msg = HybridBaseLB::AssembleStats();
  // Reset the background time
  msg->bg_walltime = 0.0;
  return msg;
#else
  return NULL;
#endif
}


void HierarchOrbLB::refine(LDStats* stats)
{
  int obj;
  int n_pes = stats->nprocs();

  // get original object mapping
  int* from_procs = Refiner::AllocProcs(n_pes, stats);
  for(obj=0;obj<stats->n_objs;obj++)  {
    int pe = stats->from_proc[obj];
    from_procs[obj] = pe;
  }

  // Get a new buffer to refine into
  int* to_procs = Refiner::AllocProcs(n_pes, stats);

  Refiner refiner(1.05);  // overload tolerance=1.05

  refiner.Refine(n_pes, stats, from_procs, to_procs);

  // Save output
  for(obj=0;obj<stats->n_objs;obj++) {
      int pe = stats->from_proc[obj];
      if (to_procs[obj] != pe) {
        if (_lb_args.debug()>=2)  {
	  CkPrintf("[%d] Obj %d migrating from %d to %d\n",
		 CkMyPe(),obj,pe,to_procs[obj]);
        }
	stats->to_proc[obj] = to_procs[obj];
      }
  }

  if (_lb_args.metaLbOn()) {
    stats->is_prev_lb_refine = 1;
    stats->after_lb_avg = refiner.computeAverageLoad();
    stats->after_lb_max = refiner.computeMax();

    if (_lb_args.debug() > 0)
      CkPrintf("RefineLB> Max load %lf Avg load %lf\n", stats->after_lb_max,
          stats->after_lb_avg);
  }

  // Free the refine buffers
  Refiner::FreeProcs(from_procs);
  Refiner::FreeProcs(to_procs);
}

void HierarchOrbLB::work(LDStats* stats) {
#if CMK_LBDB_ON
  LevelData *lData = levelData[currentLevel];

  if (currentLevel == 1) {
    orblb->work(stats);
  }
  else
    refine(stats);
#endif
}

#include "HierarchOrbLB.def.h"
