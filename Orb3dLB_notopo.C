#include <charm++.h>
#include "cklists.h"
#include "Orb3dLB_notopo.h"
#include "ParallelGravity.h"
#include "TopoManager.h"
#include "Vector3D.h"


extern CProxy_TreePiece treeProxy;

extern bool doDumpLB;
extern int lbDumpIteration;
extern bool doSimulateLB;

CkpvExtern(int, _lb_obj_index);

void Orb3dLB_notopo::init() {
  lbname = "Orb3dLB_notopo";
  if (CkpvAccess(_lb_obj_index) == -1)
    CkpvAccess(_lb_obj_index) = LBRegisterObjUserData(sizeof(TaggedVector3D));
}

//#define ORB3DLB_NOTOPO_DEBUG CkPrintf 

using namespace std;

CreateLBFunc_Def(Orb3dLB_notopo, "3d ORB mapping of tree piece space onto 3d processor mesh");

Orb3dLB_notopo::Orb3dLB_notopo(const CkLBOptions &opt): CentralLB(opt)
{
  init();
  if (CkMyPe() == 0){
    CkPrintf("[%d] Orb3dLB_notopo created\n",CkMyPe());
  }
  haveTPCentroids = false;
}

void Orb3dLB_notopo::receiveCentroids(CkReductionMsg *msg){
  if(haveTPCentroids){
    delete tpmsg;
  }
  tpCentroids = (TaggedVector3D *)msg->getData();
  nrecvd = msg->getSize()/sizeof(TaggedVector3D);
  tpmsg = msg;
  haveTPCentroids = true;
  treeProxy.doAtSync();
  CkPrintf("Orb3dLB_notopo: receiveCentroids %d elements, msg length: %d\n", nrecvd, msg->getLength());
}


bool Orb3dLB_notopo::QueryBalanceNow(int step){
  if(step == 0) return false;
  return true;
}

void Orb3dLB_notopo::work(BaseLB::LDStats* stats)
{
  int numobjs = stats->n_objs;
  int nmig = stats->n_migrateobjs;
  double gstarttime = CkWallTimer();

  stats->makeCommHash();

  vector<Event> tpEvents[NDIMS];
  for(int i = 0; i < NDIMS; i++){
    tpEvents[i].reserve(numobjs);
  }
  tps.resize(numobjs);

  OrientedBox<float> box;

  int numProcessed = 0;

  if(doSimulateLB){
    FILE *dumpFile = fopen("lb_dump.dat","r");
    CkAssert(dumpFile != NULL);

    PUP::fromTextFile pff(dumpFile);
    pupDump(pff,stats,tpEvents);

    CkPrintf("read dump\n");
    fclose(dumpFile);
  }
  else{

    for(int i = 0; i < stats->n_objs; i++){
      float load;
      load = stats->objData[i].wallTime;

      LDObjData &odata = stats->objData[i];
      TaggedVector3D* udata = (TaggedVector3D *)odata.getUserData(CkpvAccess(_lb_obj_index));

      tpEvents[XDIM].push_back(Event(udata->vec.x,load,i));
      tpEvents[YDIM].push_back(Event(udata->vec.y,load,i));
      tpEvents[ZDIM].push_back(Event(udata->vec.z,load,i));

      tps[i] = OrbObject(i,udata->myNumParticles);
      tps[i].centroid = udata->vec;
      numProcessed++;
    }
  }

  int flagDump = doDumpLB;
  CkPrintf("doDump %d dumpIter %d step %d\n", flagDump, lbDumpIteration, step());


  if(lbDumpIteration > 0 && (step() == lbDumpIteration)){
    FILE *dumpFile = fopen("lb_dump.dat","w");
    CkAssert(dumpFile != NULL);

    PUP::toTextFile ptf(dumpFile);
    pupDump(ptf,stats,tpEvents);

    CkPrintf("done LB dump, exiting\n");
    fclose(dumpFile);

    CkExit();
    return;
  }

  orbPrepare(tpEvents, box, numobjs, stats);
  orbPartition(tpEvents,box,stats->count, tps, stats);
  
  myrefine(stats);
  refine(stats, numobjs);

  int mcount = 0;
	for(int i = 0; i < numobjs; i++) {
    if (stats->to_proc[i] != stats->from_proc[i]) {
      mcount++;
    }
  }

  CkPrintf("[%d] OrbLB_notopo> migrations count %d started at %f total time %f s\n",
      CkMyPe(), mcount, gstarttime, CkWallTimer() - gstarttime);
  
  if(_lb_args.debug() >= 2) {
	// Write out "particle file" of load balance information
	char achFileName[1024];
	sprintf(achFileName, "lb.%d.sim", step());
	FILE *fp = fopen(achFileName, "w");
	CkAssert(fp != NULL);
	fprintf(fp, "%d %d 0\n", numobjs, numobjs);
	for(int i = 0; i < numobjs; i++) {
	    CkAssert(tps[i].lbindex < stats->n_objs);
	    CkAssert(tps[i].lbindex >= 0);
	    fprintf(fp, "%g %g %g %g 0.0 0.0 0.0 %d 0.0\n",
		stats->objData[tps[i].lbindex].wallTime,
		tps[i].centroid.x,
		tps[i].centroid.y,
		tps[i].centroid.z,
		stats->to_proc[tps[i].lbindex]);
	    }
	fclose(fp);
	}

  if(doSimulateLB){
    CkExit();
    return;
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


void Orb3dLB_notopo::myrefine(BaseLB::LDStats *stats) {
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
    CkPrintf("[%d] No underloaded or overloaded PE\n", CkMyPe());
    return;
  }
  // make a max heap
  make_heap(ovldpes.begin(), ovldpes.end(), PeLdGreater(counts));
  sort(unldpes.begin(), unldpes.end(), PeLdLesser(counts));

  //CkPrintf("[%d] Maxpewithtps %d thresholds %d and %d\n", pewithmax, maxcountoftps, threshold1, threshold2);
  int undcount = 0;
  CkPrintf("[%d] %d is the maxLoadedPe with ld %f ovlded %d unldpes %d\n", CkMyPe(), ovldpes.front(), counts[ovldpes.front()], ovldpes.size(), unldpes.size());
  CkPrintf("[%d] %d is the minLoadedPe with ld %f\n", CkMyPe(), unldpes.front(), counts[unldpes.front()]);

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
  CkPrintf("[%d] %d Afterwards is the maxLoadedPe with ld %f\n", CkMyPe(), ovldpes.front(), counts[ovldpes.front()]);
  delete[] counts;
  delete[] tmpcounts;
}

void Orb3dLB_notopo::pupDump(PUP::er &p, BaseLB::LDStats *stats, vector<Event> *tpEvents){
  stats->pup(p);
  p|stats->count;
  for(int i = XDIM; i <= ZDIM; i++){
    p|tpEvents[i];
  }
  p|tps;
}

#include "Orb3dLB_notopo.def.h"
