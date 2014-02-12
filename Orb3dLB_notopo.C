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

void Orb3dLB_notopo::setCentroid(CkReductionMsg *msg) {
  if(haveTPCentroids){
    delete tpmsg;
  }
  tpCentroids = (TaggedVector3D *)msg->getData();
  nrecvd = msg->getSize()/sizeof(TaggedVector3D);
  tpmsg = msg;
  haveTPCentroids = true;
//  CkPrintf("[%d] Orb3dLB_notopo: receiveCentroids %d elements, msg length: %d\n", CkMyPe(), nrecvd, msg->getLength()); 
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
  for(int i = 0; i < stats->n_objs; i++){
    LDObjData &odata = stats->objData[i];
    TaggedVector3D* udata = (TaggedVector3D *)odata.getUserData(CkpvAccess(_lb_obj_index));
    if (udata->tag == 21984) {
      CkPrintf("[%d] tag: %d ld: %f, activepart: %d, totalpart: %d\n", CkMyPe(),
      udata->tag, stats->objData[i].wallTime, udata->numActiveParticles,
      udata->myNumParticles);
    }

    if(udata->numActiveParticles == 0){
      if(stats->objData[i].migratable){

       // if (stats->from_proc[i] == 512) {
       //   LDObjData &odata = stats->objData[i];
       //   TaggedVector3D* udata = (TaggedVector3D *)odata.getUserData(CkpvAccess(_lb_obj_index));
       //   CkPrintf("[%d] Making it non-migratable Obj %d tag:%d is from outside but is not migratable activepart %d total part %d ld %f\n", CkMyPe(), i, udata->tag, udata->numActiveParticles, udata->myNumParticles, stats->objData[i].wallTime);
       // }
        stats->objData[i].migratable = 0;
        stats->n_migrateobjs--;
      }
    }
  }

  int numobjs = stats->n_objs;
  int nmig = stats->n_migrateobjs;
//  int nmig = stats->n_objs;
  double gstarttime = CkWallTimer();

  stats->makeCommHash();

  vector<Event> tpEvents[NDIMS];
  for(int i = 0; i < NDIMS; i++){
    tpEvents[i].reserve(nmig);
  }
  tps.resize(nmig);
  CkPrintf("[%d] Total migratable %d and total objs %d\n", CkMyPe(), nmig, numobjs);

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
      if(!stats->objData[i].migratable) {
//        if (stats->from_proc[i] == 512) {
//          LDObjData &odata = stats->objData[i];
//          TaggedVector3D* udata = (TaggedVector3D *)odata.getUserData(CkpvAccess(_lb_obj_index));
//          CkPrintf("[%d] Obj %d tag:%d is from outside but is not migratable activepart %d ld %f\n", CkMyPe(), i, udata->tag, udata->numActiveParticles, stats->objData[i].wallTime);
//        }
        continue;
      }

      float load;
      load = stats->objData[i].wallTime;

      LDObjData &odata = stats->objData[i];
      TaggedVector3D* udata = (TaggedVector3D *)odata.getUserData(CkpvAccess(_lb_obj_index));

      tpEvents[XDIM].push_back(Event(udata->vec.x,load,numProcessed));
      tpEvents[YDIM].push_back(Event(udata->vec.y,load,numProcessed));
      tpEvents[ZDIM].push_back(Event(udata->vec.z,load,numProcessed));

      tps[numProcessed] = OrbObject(i,udata->myNumParticles);
      tps[numProcessed].centroid = udata->vec;
      numProcessed++;
    }
  }
  //CkPrintf("numProcessed %d numMigratable %d\n", numProcessed, nmig);
  if (nmig <= 0) {
    CkPrintf("[%d] No migratable objects for load balancing\n", CkMyPe());
    return;
  }

  int flagDump = doDumpLB;
//  CkPrintf("doDump %d dumpIter %d step %d\n", flagDump, lbDumpIteration, step());
//  if(lbDumpIteration > 0 && (step() == lbDumpIteration)){
//    FILE *dumpFile = fopen("lb_dump.dat","w");
//    CkAssert(dumpFile != NULL);
//
//    PUP::toTextFile ptf(dumpFile);
//    pupDump(ptf,stats,tpEvents);
//
//    CkPrintf("done LB dump, exiting\n");
//    fclose(dumpFile);
//
//    CkExit();
//    return;
//  }

  orbPrepare(tpEvents, box, nmig, stats);
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
