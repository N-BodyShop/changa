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

CreateLBFunc_Def(MultistepLB_notopo, "Works best with multistepped runs; uses Orb3D_notopo for larger steps, greedy otherwise");

void MultistepLB_notopo::init() {
  lbname = "MultistepLB_notopo";
  if (CkpvAccess(_lb_obj_index) == -1)
    CkpvAccess(_lb_obj_index) = LBRegisterObjUserData(sizeof(TaggedVector3D));
}


MultistepLB_notopo::MultistepLB_notopo(const CkLBOptions &opt): CentralLB(opt)
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

/// Threshold between ORB3D (large) and greedy (small) as fraction of
/// active particles
#define LARGE_PHASE_THRESHOLD 0.0001

/// @brief Implement load balancing: store loads and decide between
/// ORB3D and greedy.
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
  float *ratios = new float[stats->n_objs];
  // save pointers to centroids of treepieces

  int numActiveObjects = 0;
  int numInactiveObjects = 0;

  // to calculate ratio of active particles in phase
  int64_t numActiveParticles = 0;
  int64_t totalNumParticles = 0;

  for(int i = 0; i < stats->n_objs; i++){
    stats->to_proc[i] = stats->from_proc[i];
  }

  for(int i = 0; i < stats->n_objs; i++){
    if (!stats->objData[i].migratable) continue;

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
      if (!stats->objData[i].migratable) continue;

      LDObjData &odata = stats->objData[i];
      TaggedVector3D* udata =
        (TaggedVector3D *)odata.getUserData(CkpvAccess(_lb_obj_index));
      if(!stats->objData[i].migratable && udata->myNumParticles > 0) {
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
  makeActiveProcessorList(stats, numActiveObjects);
  count = stats->count;

  // let the strategy take over on this modified instrumented data and processor information
  if((float)numActiveParticles/totalNumParticles > LARGE_PHASE_THRESHOLD){
    if (_lb_args.debug()>=2) {
      CkPrintf("******** BIG STEP *********!\n");
    }
    work2(stats,count);
  }     // end if phase == 0
  else{
    // greedy(stats,count,phase,prevPhase);
  }
#endif //CMK_LDB_ON

}

//**************************************
// ORB3DLB functions
//**************************************
//
void MultistepLB_notopo::greedy(BaseLB::LDStats *stats, int count){

  int numobjs = stats->n_objs;
  int nmig = stats->n_migrateobjs;
  CkPrintf("[GREEDY] objects total %d active %d\n", numobjs,nmig);

  TPObject *tp_array = new TPObject[nmig];
  int j = 0;
  for(int i = 0; i < stats->n_objs; i++){
    if(!stats->objData[i].migratable) continue;

    tp_array[j].migratable = stats->objData[i].migratable;
    LDObjData &odata = stats->objData[i];
    TaggedVector3D* udata = (TaggedVector3D *)odata.getUserData(CkpvAccess(_lb_obj_index));
    if(step() == 0){
      tp_array[j].load = udata->myNumParticles;
    }
    else{
      tp_array[j].load = stats->objData[i].wallTime;
    }
    tp_array[j].lbindex = i;
    j++;
  }
  mapping = &stats->to_proc;

  CkAssert(j==nmig);

  std::priority_queue<TPObject> objects;
  std::priority_queue<Processor> processors;

  for(int i = 0; i < nmig; i++){
    objects.push(tp_array[i]);
  }

  for(int i = 0; i < count; i++){
    processors.push(Processor(i));
  }

  while(!objects.empty()){
    TPObject obj = objects.top();
    objects.pop();

    Processor p = processors.top();
    processors.pop();

    p.load += obj.load;
    (*mapping)[obj.lbindex] = p.t;

    processors.push(p);
  }

  // diagnostics
  /*
  CkPrintf("**********************************\n");
  CkPrintf("GREEDY CPU LOAD PREDICTIONS phase %d\n", phase);
  CkPrintf("**********************************\n");
  while(!processors.empty()){
    Processor p = processors.top();
    processors.pop();
    CkPrintf("proc %d load %f\n", p.t, p.load);
  }
  */

  CkPrintf("**********************************\n");
  CkPrintf("GREEDY MEASURED CPU LOAD prev %d\n");
  CkPrintf("**********************************\n");
  for(int i = 0; i < stats->count; i++){
    CkPrintf("[pestats] %d %g %g\n",
                               i,
                               stats->procs[i].total_walltime,
                               stats->procs[i].idletime);
  }


  delete []tp_array;
}

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
  CentralLB::pup(p);
}

#include "MultistepLB_notopo.def.h"
