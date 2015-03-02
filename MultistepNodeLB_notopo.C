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
}

bool MultistepNodeLB_notopo::QueryBalanceNow(int step){
 if(CkMyPe() == 0) CkPrintf("Orb3dLB_notopo: Step %d\n", step);
 return true;
}

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
  
  if(_lb_args.debug() >= 2 && step() > 0) {
      // Write out "particle file" of measured load balance information
      char achFileName[1024];
      sprintf(achFileName, "lb_a.%d.sim", step()-1);
      FILE *fp = fopen(achFileName, "w");
      CkAssert(fp != NULL);
      fprintf(fp, "%d %d 0\n", stats->n_objs, stats->n_objs);
      for(int i = 0; i < stats->n_objs; i++) {
        if(!stats->objData[i].migratable) continue;

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
 
  for(int i = 0; i < stats->n_objs; i++){
    stats->to_proc[i] = stats->from_proc[i];
  }
  
  for(int i = 0; i < stats->n_objs; i++){
    if(!stats->objData[i].migratable) continue;

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
      if(!stats->objData[i].migratable) continue;

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
  count = stats->count;

  // let the strategy take over on this modified instrumented data and processor information
  if((float)numActiveParticles/totalNumParticles > LARGE_PHASE_THRESHOLD){
    if (_lb_args.debug()>=2) {
      CkPrintf("******** BIG STEP *********!\n");
    }
    work2(stats,count);
  }
#endif //CMK_LDB_ON

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
  
  orbPrepare(tpEvents, box, nmig, stats, true);
  orbPartition(tpEvents,box,CkNumNodes(),tp_array, stats, true);
  CkPrintf("MultistepLB> Done OrbPartition Total Nodes %d NodeSize %d\n", CkNumNodes(), CkMyNodeSize());
  balanceTPsNode(stats);
  refine(stats, numobjs);


  if(_lb_args.debug() >= 2) {
      // Write out "particle file" of load balance information
      char achFileName[1024];
      sprintf(achFileName, "lb.%d.sim", step());
      FILE *fp = fopen(achFileName, "w");
      CkAssert(fp != NULL);
      fprintf(fp, "%d %d 0\n", numobjs, numobjs);
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


// Refinement strategy to distribute TreePieces evenly among nodes
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
  } 
  if (ovldpes.size() == 0 || unldpes.size() == 0) {
    CkPrintf("No underloaded or overloaded Nodes\n");
    return;
  }
  // make a max heap
  make_heap(ovldpes.begin(), ovldpes.end(), PeLdGreater(counts));
  sort(unldpes.begin(), unldpes.end(), PeLdLesser(counts));

  int undcount = 0;
  //CkPrintf("[%d] is the maxLoadedPe with ld %f ovlded %d unldpes %d\n", ovldpes.front(), counts[ovldpes.front()], ovldpes.size(), unldpes.size());
  //CkPrintf("[%d] is the minLoadedPe with ld %f\n", unldpes.front(), counts[unldpes.front()]);

  int* tmpcounts = new int[numNodes];
  memcpy(tmpcounts, counts, numNodes * sizeof(int));
  srand(42);
 
  while (undcount < unldpes.size() && ovldpes.size() > 0) {
    int ovlpe = ovldpes.front();
    pop_heap(ovldpes.begin(), ovldpes.end(), PeLdGreater(counts));
    ovldpes.pop_back();
    bool succ = false;
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
      if ((counts[n_proc] + stats->objData[i].wallTime) >= (counts[to_proc_nd])) {
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
      sort(unldpes.begin(), unldpes.end(), PeLdLesser(counts));
      break;
    }
  }

  //CkPrintf("[%d] Afterwards is the maxLoadedPe with ld %f\n", ovldpes.front(), counts[ovldpes.front()]);
  delete[] counts;
  delete[] tmpcounts;
}

// Refinement strategy to distribute the TreePieces across PEs
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
    if (counts[i] > th*avgldperpe) {
      ovldpes.push_back(i);
    } else if (counts[i] < (unth*avgldperpe)) {
      unldpes.push_back(i);
    }
  } 
  if (ovldpes.size() == 0 || unldpes.size() == 0) {
    CkPrintf("No underloaded or overloaded PE\n");
    return;
  }
  // make a max heap
  make_heap(ovldpes.begin(), ovldpes.end(), PeLdGreater(counts));
  sort(unldpes.begin(), unldpes.end(), PeLdLesser(counts));

  int undcount = 0;
  //CkPrintf("[%d] is the maxLoadedPe with ld %f ovlded %d unldpes %d\n", ovldpes.front(), counts[ovldpes.front()], ovldpes.size(), unldpes.size());
  //CkPrintf("[%d] is the minLoadedPe with ld %f\n", unldpes.front(), counts[unldpes.front()]);

  int* tmpcounts = new int[stats->count];
  memcpy(tmpcounts, counts, stats->count * sizeof(int));
 
  while (undcount < unldpes.size() && ovldpes.size() > 0) {
    int ovlpe = ovldpes.front();
    pop_heap(ovldpes.begin(), ovldpes.end(), PeLdGreater(counts));
    ovldpes.pop_back();
    bool succ = false;
    if (ovlpe >= stats->count) {
      CkPrintf("ovlpe %d stats count %d\n", ovlpe, stats->count);
      CkAbort("ovle >= count\n");
    }
    for (int k = 0; k < objpemap[ovlpe].size() ; k++) {
      int i = objpemap[ovlpe][k];
      if (undcount > unldpes.size()) {
        break;
      }
      int to_proc = stats->to_proc[i];
      int n_proc = -1;
      if (to_proc != ovlpe || !stats->objData[i].migratable) {
        continue;
      }
      if (stats->objData[i].wallTime < 0.2) {
        continue;
      }

      n_proc = unldpes[undcount];
      if ((counts[n_proc] + stats->objData[i].wallTime) >= (counts[to_proc])) {
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
      sort(unldpes.begin(), unldpes.end(), PeLdLesser(counts));
      break;
    }
  }

  //CkPrintf("[%d] Afterwards is the maxLoadedPe with ld %f\n", ovldpes.front(), counts[ovldpes.front()]);
  delete[] counts;
  delete[] tmpcounts;
}

void MultistepNodeLB_notopo::pup(PUP::er &p){
  CentralLB::pup(p);
}

#include "MultistepNodeLB_notopo.def.h"
