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

//#define ORB3DLB_NOTOPO_DEBUG CkPrintf 

using namespace std;

CreateLBFunc_Def(Orb3dLB_notopo, "3d ORB mapping of tree piece space onto 3d processor mesh");

Orb3dLB_notopo::Orb3dLB_notopo(const CkLBOptions &opt): CentralLB(opt)
{
  lbname = "Orb3dLB_notopo";
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

  stats->makeCommHash();
  CkAssert(nrecvd == nmig);

  vector<Event> tpEvents[NDIMS];
  for(int i = 0; i < NDIMS; i++){
    tpEvents[i].reserve(nrecvd);
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
    for(int i = 0; i < nrecvd; i++){
      TaggedVector3D *data = tpCentroids+i;
      LDObjHandle &handle = data->handle;
      int tag = stats->getHash(handle.id,handle.omhandle.id);

      float load;
      if(step() == 0){
        load = data->myNumParticles;
      }
      else{
        load = stats->objData[tag].wallTime;
      }

      int pe = stats->from_proc[tag];

      //CkPrintf("[mydebug] %d tag %d np %d pe %d load %f vec %f %f %f\n", i, tag, data->myNumParticles, pe, load, data->vec.x, data->vec.y, data->vec.z);

      tpEvents[XDIM].push_back(Event(data->vec.x,load,tag));
      tpEvents[YDIM].push_back(Event(data->vec.y,load,tag));
      tpEvents[ZDIM].push_back(Event(data->vec.z,load,tag));

      tps[tag] = OrbObject(tag,data->myNumParticles);
      tps[tag].centroid = data->vec;

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

  orbPrepare(tpEvents, box, nrecvd, stats);
  orbPartition(tpEvents,box,stats->count, tps, stats);

  refine(stats, numobjs);
  
  if(_lb_args.debug() >= 2) {
	// Write out "particle file" of load balance information
	char achFileName[1024];
	sprintf(achFileName, "lb.%d.sim", step());
	FILE *fp = fopen(achFileName, "w");
	CkAssert(fp != NULL);
	fprintf(fp, "%d %d 0\n", nrecvd, nrecvd);
	for(int i = 0; i < nrecvd; i++) {
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
