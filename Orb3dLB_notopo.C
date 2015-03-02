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
      if (!odata.migrations) continue;

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
  int mcount = 0;
	for(int i = 0; i < numobjs; i++) {
    if (stats->to_proc[i] != stats->from_proc[i]) {
      mcount++;
    }
  }

  refine(stats, numobjs);
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
      if (!stats->objData[i].migratable) continue;

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

void Orb3dLB_notopo::pupDump(PUP::er &p, BaseLB::LDStats *stats, vector<Event> *tpEvents){
  stats->pup(p);
  p|stats->count;
  for(int i = XDIM; i <= ZDIM; i++){
    p|tpEvents[i];
  }
  p|tps;
}

#include "Orb3dLB_notopo.def.h"
