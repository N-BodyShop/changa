#include <charm++.h>
#include "cklists.h"
#include "Orb3dLB_notopo.h"
#include "ParallelGravity.h"
#include "TopoManager.h"
#include "Vector3D.h"


extern CProxy_TreePiece treeProxy;

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


CmiBool Orb3dLB_notopo::QueryBalanceNow(int step){
  if(step == 0) return false;
  return true;
}

void Orb3dLB_notopo::work(BaseLB::LDStats* stats)
{
  int numobjs = stats->n_objs;
  int nmig = stats->n_migrateobjs;

  stats->makeCommHash();
  CkAssert(nrecvd == numobjs);

  vector<Event> tpEvents[NDIMS];
  for(int i = 0; i < NDIMS; i++){
    tpEvents[i].reserve(numobjs);
  }
  tps.resize(numobjs);

  OrientedBox<float> box;

  int numProcessed = 0;

  for(int i = 0; i < numobjs; i++){
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


  orbPrepare(tpEvents, box, numobjs, stats);
  orbPartition(tpEvents,box,stats->count, tps, stats);

  refine(stats, numobjs);

}

#include "Orb3dLB_notopo.def.h"

/*@}*/
