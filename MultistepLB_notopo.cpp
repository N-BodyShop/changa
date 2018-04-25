#include <charm++.h>
#include "cklists.h"
#include "MultistepLB_notopo.h"
#include "TopoManager.h"
#include "ParallelGravity.h"
#include "Vector3D.h"
#include <queue>

extern CProxy_TreePiece treeProxy;
CkpvExtern(int, _lb_obj_index);
using namespace std;
//#define ORB3DLB_NOTOPO_DEBUG CkPrintf

CreateLBFunc_Def(MultistepLB_notopo, "Works best with multistepped runs; uses Orb3D_notopo");

void MultistepLB_notopo::init() {
  lbname = "MultistepLB_notopo";
  if (CkpvAccess(_lb_obj_index) == -1)
    CkpvAccess(_lb_obj_index) = LBRegisterObjUserData(sizeof(TaggedVector3D));
}


MultistepLB_notopo::MultistepLB_notopo(const CkLBOptions &opt): CBase_MultistepLB_notopo(opt)
{
  init();
  if (CkMyPe() == 0){
    CkPrintf("[%d] MultistepLB_notopo created\n",CkMyPe());
  }
}

bool MultistepLB_notopo::QueryBalanceNow(int step){
 if(CkMyPe() == 0) CkPrintf("Orb3dLB_notopo: Step %d\n", step);
 //  if(step == 0) return false;
  return true;

}

// helper functions for multistepping
#ifdef MCLBMS

void MultistepLB_notopo::makeActiveProcessorList(BaseLB::LDStats *stats, int numActiveObjs){
  int objsPerProc = 8;
  int expandFactor = 4;
  int procsNeeded;
  procsNeeded = expandFactor*numActiveObjs/objsPerProc > stats->count ? stats->count : expandFactor*numActiveObjs/objsPerProc;

  /* currently, only the first procsNeeded procs are used - could do something more sophisticated here in the future - FIXME */
#ifdef MCLBMSV
  CkPrintf("Processors 0 to %d active\n", procsNeeded-1);
#endif
}
#endif

/// @brief Implement load balancing: store loads and determine active
/// processors and objects, then call ORB3D.
void MultistepLB_notopo::work(BaseLB::LDStats* stats)
{
#if CMK_LBDB_ON
  // find active objects - mark the inactive ones as non-migratable
  int count;

  if(_lb_args.debug() >= 2 && step() > 0) {
      // Write out "particle file" of measured load balance information
      char achFileName[1024];
      sprintf(achFileName, "lb_a.%d.sim", step()-1);
      FILE *fp = fopen(achFileName, "w");
      CkAssert(fp != NULL);

      int num_migratables = stats->n_objs;
      for(int i = 0; i < stats->n_objs; i++) {
        if (!stats->objData[i].migratable) {
          num_migratables--;
        }
      }

      fprintf(fp, "%d %d 0\n", num_migratables, num_migratables);
      for(int i = 0; i < stats->n_objs; i++) {
        if (!stats->objData[i].migratable) continue;

      LDObjData &odata = stats->objData[i];
      TaggedVector3D* udata = (TaggedVector3D *)odata.getUserData(CkpvAccess(_lb_obj_index));
	  fprintf(fp, "%g %g %g %g 0.0 0.0 0.0 %d %d\n",
		  stats->objData[i].wallTime,
		  udata->vec.x,
		  udata->vec.y,
		  udata->vec.z,
		  stats->from_proc[i],
		  udata->tp);
	  }
      fclose(fp);
      }

  int numActiveObjects = 0;
  int numInactiveObjects = 0;
  int minActiveProc = INT_MAX;
  int maxActiveProc = 0;

  for(int i = 0; i < stats->n_objs; i++){
    stats->to_proc[i] = stats->from_proc[i];
  }

  for(int i = 0; i < stats->n_objs; i++){
    if (!stats->objData[i].migratable) continue;

    LDObjData &odata = stats->objData[i];
    TaggedVector3D* udata = (TaggedVector3D *)odata.getUserData(CkpvAccess(_lb_obj_index));

    if(udata->numActiveParticles == 0){
      numInactiveObjects++;
    }
    else{
      numActiveObjects++;
      if(minActiveProc > stats->from_proc[i])
	  minActiveProc = stats->from_proc[i];
      if(maxActiveProc < stats->from_proc[i])
	  maxActiveProc = stats->from_proc[i];
    }
  }
  CkPrintf("numActiveObjects: %d, numInactiveObjects: %d\n", numActiveObjects,
      numInactiveObjects);
  CkPrintf("active PROC range: %d to %d\n", minActiveProc, maxActiveProc);
  if(numActiveObjects < 0.1*numInactiveObjects) {
    // only a small number of active objects, only migrate them
    for(int i = 0; i < stats->n_objs; i++){
      if (!stats->objData[i].migratable) continue;

      LDObjData &odata = stats->objData[i];
      TaggedVector3D* udata =
        (TaggedVector3D *)odata.getUserData(CkpvAccess(_lb_obj_index));
      if(udata->numActiveParticles == 0) {
	  stats->objData[i].migratable = 0;
	  stats->n_migrateobjs--;
      }
    }
  }
  else {
    CkPrintf("Migrating all: numActiveObjects: %d, numInactiveObjects: %d\n", numActiveObjects, numInactiveObjects);
  }

  // select processors
#ifdef MCLBMSV
  //printData(*stats, phase, NULL);
  CkPrintf("making active processor list\n");
#endif
  makeActiveProcessorList(stats, numActiveObjects);
  count = stats->count;

  // let the strategy take over on this modified instrumented data and processor information
  work2(stats,count);
#endif //CMK_LDB_ON

}

//**************************************
// ORB3DLB functions
//**************************************
//

/// @brief ORB3D load balance.
void MultistepLB_notopo::work2(BaseLB::LDStats *stats, int count){
  int numobjs = stats->n_objs;
  int nmig = stats->n_migrateobjs;

  // this data structure is used by the orb3d strategy
  // to balance objects. it is NOT indexed by tree piece index
  // there are as many entries in it as there are
  // migratable (active) tree pieces
 vector<OrbObject> tp_array;
 tp_array.resize(nmig);

  if (_lb_args.debug()>=2) {
    CkPrintf("[work2] ready tp_array data structure\n");
  }

 vector<Event> tpEvents[NDIMS];
  for(int i = 0; i < NDIMS; i++){
    tpEvents[i].reserve(nmig);
  }

  OrientedBox<float> box;

  int numProcessed = 0;

  for(int i = 0; i < numobjs; i++){
    if(!stats->objData[i].migratable) continue;

    float load;
    LDObjData &odata = stats->objData[i];
    TaggedVector3D* udata = (TaggedVector3D *)odata.getUserData(CkpvAccess(_lb_obj_index));
    if(step() == 0){
      load = udata->myNumParticles;
    }
    else{
      load = stats->objData[i].wallTime;
    }

    tpEvents[XDIM].push_back(Event(udata->vec.x,load,numProcessed));
    tpEvents[YDIM].push_back(Event(udata->vec.y,load,numProcessed));
    tpEvents[ZDIM].push_back(Event(udata->vec.z,load,numProcessed));

    tp_array[numProcessed]= OrbObject(i, udata->myNumParticles);
    tp_array[numProcessed].centroid = udata->vec;
    numProcessed++;
  }
  CkAssert(numProcessed==nmig);

  orbPrepare(tpEvents, box, nmig, stats);
  orbPartition(tpEvents,box,stats->count,tp_array, stats);

  refine(stats, numobjs);

  if(_lb_args.debug() >= 2) {
      // Write out "particle file" of load balance information
      char achFileName[1024];
      sprintf(achFileName, "lb.%d.sim", step());
      FILE *fp = fopen(achFileName, "w");
      CkAssert(fp != NULL);

      int num_migratables = numobjs;
      for(int i = 0; i < numobjs; i++) {
        if (!stats->objData[i].migratable) {
          num_migratables--;
        }
      }
      fprintf(fp, "%d %d 0\n", num_migratables, num_migratables);

      for(int i = 0; i < numobjs; i++) {
        if(!stats->objData[i].migratable) continue;

    LDObjData &odata = stats->objData[i];
    TaggedVector3D* udata =
      (TaggedVector3D *)odata.getUserData(CkpvAccess(_lb_obj_index));
	  fprintf(fp, "%g %g %g %g 0.0 0.0 0.0 %d %d\n",
		  stats->objData[i].wallTime,
		  udata->vec.x,
		  udata->vec.y,
		  udata->vec.z,
		  stats->to_proc[i],
		  udata->tp);
	  }
      fclose(fp);
      }
}


void MultistepLB_notopo::pup(PUP::er &p){
  CBase_MultistepLB_notopo::pup(p);
}

#include "MultistepLB_notopo.def.h"
