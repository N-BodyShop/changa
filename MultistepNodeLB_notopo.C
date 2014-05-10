#include <charm++.h>
#include "cklists.h"
#include "MultistepNodeLB_notopo.h"
#include "TopoManager.h"
#include "ParallelGravity.h"
#include "Vector3D.h"
#include <queue>

extern CProxy_TreePiece treeProxy;
CkpvExtern(int, _lb_obj_index);
using namespace std;
//#define ORB3DLB_NOTOPO_DEBUG CkPrintf

CreateLBFunc_Def(MultistepNodeLB_notopo, "Works best with multistepped runs; uses Orb3D_notopo for larger steps, greedy otherwise");

//void MultistepLB_notopo::GetMaxNodeInfo(BaseLB::LDStats *stats) {
//
//  int nodeSize = CkNodeSize(0);
//  int nd = 0;
//  int pewithobj;
//
//  for(int i = 0; i < stats->n_objs; i++){
//    LDObjData &odata = stats->objData[i];
//    TaggedVector3D* udata = (TaggedVector3D *)odata.getUserData(CkpvAccess(_lb_obj_index));
//    if (udata->tp == maxPredLdSavedObj) {
//      nd = stats->to_proc[i] / nodeSize;
//      pewithobj = stats->to_proc[i];
//      break;
//    }
//  }
//
//  CkPrintf("NODE %d PE %d contains heavy loaded obj %d with load %f\n", nd, CkMyPe(), maxPredLdSavedObj, maxPredLdSaved);
//
//  for(int i = 0; i < stats->n_objs; i++){
//    LDObjData &odata = stats->objData[i];
//    TaggedVector3D* udata = (TaggedVector3D *)odata.getUserData(CkpvAccess(_lb_obj_index));
//    if (stats->to_proc[i]/nodeSize == nd) {
//      CkPrintf("\t PE [%d] Other Objs there %d load %f\n", stats->to_proc[i], udata->tp, stats->objData[i].wallTime);
//    }
//  }
//}
//
//void MultistepLB_notopo::GetLodInfo(BaseLB::LDStats *stats) {
//  double maxObjLd = 0.0;
//  double load;
//  double curload;
//  int maxLdObj = 0;
//  int prevRung;
//
//  double prevMaxPredLd = maxPredLdSaved;
//  double prevMaxActualLd;
//  int prevMaxPredLdObj = maxPredLdSavedObj;
//
//  maxPredLdSaved = 0.0;
//
//  for(int i = 0; i < stats->n_objs; i++){
//    LDObjData &odata = stats->objData[i];
//    TaggedVector3D* udata = (TaggedVector3D *)odata.getUserData(CkpvAccess(_lb_obj_index));
//    load = udata->busytime;
//    prevRung = udata->prevActiveRung;
//    if (load > maxObjLd) {
//      maxObjLd = load;
//      maxLdObj = udata->tp;
//    }
//
//    curload = stats->objData[i].wallTime;
//    if (curload > maxPredLdSaved) {
//      maxPredLdSaved = curload;
//      maxPredLdSavedObj = udata->tp;
//    }
//
//    if (udata->tp == prevMaxPredLdObj) {
//      prevMaxActualLd = udata->busytime;
//    }
//  }
//  CkPrintf("Previous Load for rung %d is maxloadedobj %d with load %f PRED obj %d with load %f actual load of that is %f\n", prevRung, maxLdObj, maxObjLd, prevMaxPredLdObj, prevMaxPredLd, prevMaxActualLd); 
//}

void MultistepNodeLB_notopo::init() {
  lbname = "MultistepNodeLB_notopo";
  if (CkpvAccess(_lb_obj_index) == -1)
    CkpvAccess(_lb_obj_index) = LBRegisterObjUserData(sizeof(TaggedVector3D));
}


MultistepNodeLB_notopo::MultistepNodeLB_notopo(const CkLBOptions &opt): CentralLB(opt)
{
  init();
  if (CkMyPe() == 0){
    CkPrintf("[%d] MultistepNodeLB_notopo created\n",CkMyPe());
  }
  //existing_tps_on_pe = NULL;
}

bool MultistepNodeLB_notopo::QueryBalanceNow(int step){
 if(CkMyPe() == 0) CkPrintf("Orb3dLB_notopo: Step %d\n", step);
 //  if(step == 0) return false;
  return true;

}

// helper functions for multistepping
#ifdef MCLBMS

void MultistepNodeLB_notopo::makeActiveProcessorList(BaseLB::LDStats *stats, int numActiveObjs){
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
void MultistepNodeLB_notopo::work(BaseLB::LDStats* stats)
{
#if CMK_LBDB_ON
  // find active objects - mark the inactive ones as non-migratable
  int count;
  //GetLodInfo(stats);
  
  if(_lb_args.debug() >= 2 && step() > 0) {
      // Write out "particle file" of measured load balance information
      char achFileName[1024];
      sprintf(achFileName, "lb_a.%d.sim", step()-1);
      FILE *fp = fopen(achFileName, "w");
      CkAssert(fp != NULL);
      fprintf(fp, "%d %d 0\n", stats->n_objs, stats->n_objs);
      for(int i = 0; i < stats->n_objs; i++) {

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

  // to calculate ratio of active particles in phase
  int numActiveParticles = 0;
  int totalNumParticles = 0;
 
//  if (step() == 0) {
//    existing_tps_on_pe = new int[stats->count];
//  } 
//  memset(existing_tps_on_pe, 0, stats->count * sizeof(int));

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
void MultistepNodeLB_notopo::greedy(BaseLB::LDStats *stats, int count){

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
void MultistepNodeLB_notopo::work2(BaseLB::LDStats *stats, int count){
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
  
  //orbPrepare(tpEvents, box, nmig, stats);
  //orbPartition(tpEvents,box,stats->count,tp_array, stats);
  //refine(stats, numobjs);

  orbPrepareMS(tpEvents, box, nmig, stats);
  orbPartitionMS(tpEvents,box,CkNumNodes(),tp_array, stats);
  CkPrintf("MultistepLB> Done OrbPartition Total Nodes %d NodeSize %d\n", CkNumNodes(), CkMyNodeSize());
  balanceTPsNode(stats);
  //// TODO: Cm
  //if (phase < 3) {
  //  //balanceTPsNode(stats);
  //  balanceTPs(stats);
  //  CkPrintf("MultistepLB> Done BalanceTPs\n");
  //} else {
  //  balanceTPsNode(stats);
  //}
  //refineMS(stats, numobjs, &prevMaxPredPe, existing_tps_on_pe);
  //GetMaxNodeInfo(stats);
  refine(stats, numobjs);


  if(_lb_args.debug() >= 2) {
      // Write out "particle file" of load balance information
      char achFileName[1024];
      sprintf(achFileName, "lb.%d.sim", step());
      FILE *fp = fopen(achFileName, "w");
      CkAssert(fp != NULL);
      fprintf(fp, "%d %d 0\n", numobjs, numobjs);
      for(int i = 0; i < numobjs; i++) {
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

class PeLdLesser {
  private:
  double* s;
  public:
  PeLdLesser(double* stats) {
    s = stats;
  }

  bool operator()(int l, int r) {
    return (s[l] < s[r]);
  }
};


class PeLdGreater {
  private:
  double* s;
  public:
  PeLdGreater(double* stats) {
    s = stats;
  }

  bool operator()(int l, int r) {
    return (s[l] < s[r]);
  }
};


void MultistepNodeLB_notopo::balanceTPsNode(BaseLB::LDStats* stats) {
  
  int numNodes = CkNumNodes();
  int nodeSize = CkNodeSize(0);
  double* counts = new double[numNodes];
  memset(counts, 0.0, numNodes * sizeof(double));
  double totalld = 0.0;
  vector<vector<int> > objpemap;
  objpemap.resize(numNodes);
  
  for (int i = 0; i < stats->n_objs; i++) {
    if(!stats->objData[i].migratable) continue;

    int nd = stats->to_proc[i]/nodeSize;
    counts[nd] += stats->objData[i].wallTime;
    totalld += stats->objData[i].wallTime;
    objpemap[nd].push_back(i);
  }
  double avgldperpe = totalld / numNodes;
  vector<int> unldpes;
  vector<int> ovldpes;
  double th = 1.20;
  double uth = 0.9;
  
  int maxcountoftps = 0;
  int pewithmax = -1;
  for (int i = (numNodes-1); i >= 0; i--) {
    if (counts[i] > th*avgldperpe) {
      ovldpes.push_back(i);
    } else if (counts[i] < (uth*avgldperpe)) {
      unldpes.push_back(i);
    }
   // if (counts[i] > maxcountoftps) {
   //   maxcountoftps = counts[i];
   //   pewithmax = i;
   // }
  } 
  if (ovldpes.size() == 0 || unldpes.size() == 0) {
    CkPrintf("No underloaded or overloaded Nodes\n");
    return;
  }
  // make a max heap
  make_heap(ovldpes.begin(), ovldpes.end(), PeLdGreater(counts));
  sort(unldpes.begin(), unldpes.end(), PeLdLesser(counts));

  //CkPrintf("[%d] Maxpewithtps %d thresholds %d and %d\n", pewithmax, maxcountoftps, threshold1, threshold2);
  int undcount = 0;
  CkPrintf("[%d] is the maxLoadedPe with ld %f ovlded %d unldpes %d\n", ovldpes.front(), counts[ovldpes.front()], ovldpes.size(), unldpes.size());
  CkPrintf("[%d] is the minLoadedPe with ld %f\n", unldpes.front(), counts[unldpes.front()]);

  int* tmpcounts = new int[numNodes];
  memcpy(tmpcounts, counts, numNodes * sizeof(int));
  srand(42);
 
  while (undcount < unldpes.size() && ovldpes.size() > 0) {
    int ovlpe = ovldpes.front();
    pop_heap(ovldpes.begin(), ovldpes.end(), PeLdGreater(counts));
    ovldpes.pop_back();
    bool succ = false;
    //for (int i = stats->n_objs-1; i >=0 ; i--) {
    for (int k = 0; k < objpemap[ovlpe].size() ; k++) {
      int i = objpemap[ovlpe][k];
      if (undcount > unldpes.size()) {
        break;
      }
      int to_proc = stats->to_proc[i];
      int to_proc_nd = to_proc/nodeSize;
      int n_proc = -1;
      if (to_proc_nd != ovlpe || !stats->objData[i].migratable) {
        continue;
      }
      if (stats->objData[i].wallTime < 0.2) {
        continue;
      }

      n_proc = unldpes[undcount];
      //if ((counts[n_proc] + stats->objData[i].wallTime) >= (counts[to_proc_nd]-stats->objData[i].wallTime)) {
      if ((counts[n_proc] + stats->objData[i].wallTime) >= (counts[to_proc_nd])) {
        //if(step() == 8) {
        //  CkPrintf("Could not transfer to undloaded pe ld %f obj %f\n",
        //  counts[n_proc], stats->objData[i].wallTime);
        //}
        continue;
      }
      int rand_no = rand() % nodeSize;
      int rand_pe_in_nd = n_proc * nodeSize + rand_no; 
      stats->to_proc[i] = rand_pe_in_nd;
      counts[to_proc_nd] = counts[to_proc_nd] - stats->objData[i].wallTime;
      counts[n_proc] = counts[n_proc] + stats->objData[i].wallTime;
      objpemap[n_proc].push_back(i);
      if (counts[n_proc] > uth*avgldperpe) {
        undcount++;
      }
      if (counts[to_proc_nd] > th*avgldperpe) {
        ovldpes.push_back(ovlpe);
        push_heap(ovldpes.begin(), ovldpes.end(), PeLdGreater(counts));
      }
      if (counts[n_proc] > th*avgldperpe) {
        ovldpes.push_back(n_proc);
        push_heap(ovldpes.begin(), ovldpes.end(), PeLdGreater(counts));
      }
      succ = true;
      //if(step() == 8) {
      //  CkPrintf("transfered to undloaded pe ld %f obj %f new ld for PE %d %f\n",
      //      counts[n_proc], stats->objData[i].wallTime, to_proc, counts[to_proc_nd]);
      //}
      sort(unldpes.begin(), unldpes.end(), PeLdLesser(counts));
      break;

      // if (counts[to_proc] > 70) {
      //   CkPrintf("TP %d is on big proc %d ismigr? %d ld %f\n", i, to_proc,
      //   stats->objData[i].migratable, stats->objData[i].wallTime);
      // }
     // if (counts[to_proc] > 1.5*avgldperpe  && stats->objData[i].wallTime > 1.0) {
     //   n_proc = unldpes[undcount];
     //   stats->to_proc[i] = n_proc;
     //   counts[to_proc] = counts[to_proc] - stats->objData[i].wallTime;
     //   counts[n_proc] = counts[n_proc] + stats->objData[i].wallTime;
     //   if (counts[n_proc] > 0.8*avgldperpe) {
     //     undcount++;
     //   }
     // }
    }
    //if (!succ && step()==8) {

    //  CkPrintf("Couldn't find any obj to migrate from %d pe with ld %f\n",
    //  ovlpe, counts[ovlpe]);
    //}
  }

//  maxcountoftps = 0;
//  pewithmax = -1;
//  for (int i = (stats->count-1); i >= 0; i--) {
//    if (counts[i] > maxcountoftps) {
//      maxcountoftps = counts[i];
//      pewithmax = i;
//    }
//  }
//  CkPrintf("[%d] Afterwards Maxpewithtps %d previously %d\n", pewithmax, maxcountoftps, tmpcounts[pewithmax]);
  CkPrintf("[%d] Afterwards is the maxLoadedPe with ld %f\n", ovldpes.front(), counts[ovldpes.front()]);
  delete[] counts;
  delete[] tmpcounts;
}


void MultistepNodeLB_notopo::balanceTPs(BaseLB::LDStats* stats) {
  
  double* counts = new double[stats->count];
  memset(counts, 0.0, stats->count * sizeof(double));
  double totalld = 0.0;
  vector<vector<int> > objpemap;
  objpemap.resize(stats->count);
  
  for (int i = 0; i < stats->n_objs; i++) {
    if(!stats->objData[i].migratable) continue;

    counts[stats->to_proc[i]] += stats->objData[i].wallTime;
    totalld += stats->objData[i].wallTime;
    objpemap[stats->to_proc[i]].push_back(i);
  }
  double avgldperpe = totalld / stats->count;
  vector<int> unldpes;
  vector<int> ovldpes;
  double th = 1.05;
  double unth = 0.9;
  int total_its = 0;
  
  int maxcountoftps = 0;
  int pewithmax = -1;
  for (int i = (stats->count-1); i >= 0; i--) {
    //if (counts[i] > 4.9) {
    //  CkPrintf("[%d] PE is overloaded with load %f\n", i, counts[i]);
    //}
    if (counts[i] > th*avgldperpe) {
      ovldpes.push_back(i);
    } else if (counts[i] < (unth*avgldperpe)) {
      unldpes.push_back(i);
    }
   // if (counts[i] > maxcountoftps) {
   //   maxcountoftps = counts[i];
   //   pewithmax = i;
   // }
  } 
  if (ovldpes.size() == 0 || unldpes.size() == 0) {
    CkPrintf("No underloaded or overloaded PE\n");
    return;
  }
  // make a max heap
  make_heap(ovldpes.begin(), ovldpes.end(), PeLdGreater(counts));
  sort(unldpes.begin(), unldpes.end(), PeLdLesser(counts));

  //CkPrintf("[%d] Maxpewithtps %d thresholds %d and %d\n", pewithmax, maxcountoftps, threshold1, threshold2);
  int undcount = 0;
  CkPrintf("[%d] is the maxLoadedPe with ld %f ovlded %d unldpes %d\n", ovldpes.front(), counts[ovldpes.front()], ovldpes.size(), unldpes.size());
  CkPrintf("[%d] is the minLoadedPe with ld %f\n", unldpes.front(), counts[unldpes.front()]);

  int* tmpcounts = new int[stats->count];
  memcpy(tmpcounts, counts, stats->count * sizeof(int));
 
  //while (undcount < unldpes.size() && ovldpes.size() > 0 && total_its < stats->count) {
  while (undcount < unldpes.size() && ovldpes.size() > 0) {
    int ovlpe = ovldpes.front();
    pop_heap(ovldpes.begin(), ovldpes.end(), PeLdGreater(counts));
    ovldpes.pop_back();
    bool succ = false;
    if (ovlpe >= stats->count) {
      CkPrintf("ovlpe %d stats count %d\n", ovlpe, stats->count);
      CkAbort("ovle >= count\n");
    }
    //CkPrintf("[%d] ld %f Looking at to make it underloaded contains %d objs\n", ovlpe, counts[ovlpe], objpemap[ovlpe].size());
    //for (int i = stats->n_objs-1; i >=0 ; i--) {
    for (int k = 0; k < objpemap[ovlpe].size() ; k++) {
      int i = objpemap[ovlpe][k];
      if (undcount > unldpes.size()) {
        break;
      }
      int to_proc = stats->to_proc[i];
      int n_proc = -1;
      if (to_proc != ovlpe || !stats->objData[i].migratable) {
        if (counts[ovlpe] > 5.0) {
          //CkPrintf("Working on pe %d ld %f to_proc %d\n", ovlpe, counts[ovlpe], to_proc);
        }
        continue;
      }
      if (stats->objData[i].wallTime < 0.2) {
       if (counts[ovlpe] > 5.0) {
          //CkPrintf("Working on pe %d ld %f walltime %f\n", ovlpe, counts[ovlpe], stats->objData[i].wallTime);
        }
        continue;
      }

      n_proc = unldpes[undcount];
      //if ((counts[n_proc] + stats->objData[i].wallTime) >= (counts[to_proc]-stats->objData[i].wallTime)) {
      if ((counts[n_proc] + stats->objData[i].wallTime) >= (counts[to_proc])) {
        if(step() == 16) {
          //CkPrintf("Could not transfer to undloaded pe ld %f obj %f\n", counts[n_proc], stats->objData[i].wallTime);
        }
        continue;
      }
      stats->to_proc[i] = n_proc;
      counts[to_proc] = counts[to_proc] - stats->objData[i].wallTime;
      counts[n_proc] = counts[n_proc] + stats->objData[i].wallTime;
      objpemap[n_proc].push_back(i);
      if (counts[n_proc] > unth*avgldperpe) {
        undcount++;
      }
      if (counts[to_proc] > th*avgldperpe) {
        ovldpes.push_back(ovlpe);
        push_heap(ovldpes.begin(), ovldpes.end(), PeLdGreater(counts));
      }

      if (counts[n_proc] > th*avgldperpe) {
        ovldpes.push_back(n_proc);
        push_heap(ovldpes.begin(), ovldpes.end(), PeLdGreater(counts));
      }
      succ = true;
      if(step() == 16) {
        //CkPrintf("transfered to undloaded pe %d ld %f obj %f new ld for PE %d %f\n", n_proc, counts[n_proc], stats->objData[i].wallTime, to_proc, counts[to_proc]);
      }
      sort(unldpes.begin(), unldpes.end(), PeLdLesser(counts));
      break;

      // if (counts[to_proc] > 70) {
      //   CkPrintf("TP %d is on big proc %d ismigr? %d ld %f\n", i, to_proc,
      //   stats->objData[i].migratable, stats->objData[i].wallTime);
      // }
     // if (counts[to_proc] > 1.5*avgldperpe  && stats->objData[i].wallTime > 1.0) {
     //   n_proc = unldpes[undcount];
     //   stats->to_proc[i] = n_proc;
     //   counts[to_proc] = counts[to_proc] - stats->objData[i].wallTime;
     //   counts[n_proc] = counts[n_proc] + stats->objData[i].wallTime;
     //   if (counts[n_proc] > 0.8*avgldperpe) {
     //     undcount++;
     //   }
     // }
    }
    if (!succ && step()==16) {

      //CkPrintf("Couldn't find any obj to migrate from %d pe with ld %f\n", ovlpe, counts[ovlpe]);
    }
  }

//  maxcountoftps = 0;
//  pewithmax = -1;
//  for (int i = (stats->count-1); i >= 0; i--) {
//    if (counts[i] > maxcountoftps) {
//      maxcountoftps = counts[i];
//      pewithmax = i;
//    }
//  }
//  CkPrintf("[%d] Afterwards Maxpewithtps %d previously %d\n", pewithmax, maxcountoftps, tmpcounts[pewithmax]);
  CkPrintf("[%d] Afterwards is the maxLoadedPe with ld %f\n", ovldpes.front(), counts[ovldpes.front()]);
  delete[] counts;
  delete[] tmpcounts;
}



void MultistepNodeLB_notopo::pup(PUP::er &p){
  CentralLB::pup(p);
}

#include "MultistepNodeLB_notopo.def.h"
