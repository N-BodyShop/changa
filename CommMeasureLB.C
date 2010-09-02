#include <charm++.h>
#include "cklists.h"
#include "CommMeasureLB.h"
#include "ParallelGravity.h"
#include "TopoManager.h"
#include "Vector3D.h"

extern CProxy_TreePiece treeProxy;
extern CProxy_DataManager dMProxy;
extern CProxy_CkCacheManager cacheGravPart;
extern CProxy_CkCacheManager cacheGravSmooth;
extern CProxy_CkCacheManager cacheGravNode;

using namespace std;

CreateLBFunc_Def(CommMeasureLB, "Base decisions on recorded communication graph");

CommMeasureLB::CommMeasureLB(const CkLBOptions &opt): CentralLB(opt)
{
  lbname = "CommMeasureLB";
  centroidsAllocated = false;
  if (CkMyPe() == 0){
    CkPrintf("[%d] CommMeasureLB created\n",CkMyPe());
  }
}

void CommMeasureLB::receiveCentroids(CkReductionMsg *msg){
  int i = 0;
   TaggedVector3D * cur = (TaggedVector3D *)msg->getData();
  CkPrintf("CommMeasureLB: receiveCentroids started: %d elements, msg length: %d\n", msg->getGcount(), msg->getLength()); 
  tpCentroids.free();
  
  while(i < msg->getGcount()){
     tpCentroids.push_back(*cur);
     cur = cur + 1;
     i++;
  }
  treeProxy.doAtSync();
  CkPrintf("CommMeasureLB: receiveCentroids done\n");  
  delete msg;
}

//jetley
CmiBool CommMeasureLB::QueryBalanceNow(int step){
  if(step == 0){
    if(CkMyPe() == 0){                          // only one group member need broadcast
      CkPrintf("CommMeasureLB: Step 0, calling treeProxy.receiveProxy(thisgroup)\n");
      treeProxy.receiveProxy(thisgroup);        // broadcast proxy to all treepieces
    }
    firstRound = true;
    return false; 
  }
  if(CkMyPe() == 0)
    CkPrintf("CommMeasureLB: Step %d\n", step);
  return true;

}

void CommMeasureLB::work(BaseLB::LDStats* stats, int count)
{
  int numobjs = stats->n_objs;

  stats->makeCommHash();
  Vector3D<float> *lbcentroids = new Vector3D<float>[numobjs];
  
  //CkPrintf("checking before centroids\n");
  //CmiMemoryCheck();
  CkAssert(tpCentroids.length() == numobjs);
  for(int i = 0; i < numobjs; i++){
    LDObjHandle &handle = tpCentroids[i].handle;
    int tag = stats->getHash(handle.id,handle.omhandle.id);
    lbcentroids[tag] = tpCentroids[i].vec;
  }
  //CkPrintf("checking before comm\n");
  //CmiMemoryCheck();

  const int csz = stats->n_comm;
  for(int i=0; i<csz; i++) {
      LDCommData &cdata = stats->commData[i];
      if(!cdata.from_proc() && cdata.receiver.get_type() == LD_OBJ_MSG)
      {
        int senderID = stats->getHash(cdata.sender);
        int recverID = stats->getHash(cdata.receiver.get_destObj());
        if (stats->complete_flag == 0 && recverID == -1) continue;
        CmiAssert(senderID < numobjs && senderID >= 0);
        CmiAssert(recverID < numobjs && recverID >= 0);
        //comm[senderID][recverID] += cdata.messages;
        //comm[recverID][senderID] += cdata.messages;
        Vector3D<float> dist = lbcentroids[senderID]-lbcentroids[recverID];
        CkPrintf("[edge] %d %d %d %f\n", senderID, recverID, cdata.bytes, dist.length());


      }
      /*
      else if (cdata.receiver.get_type() == LD_OBJLIST_MSG) {
        int nobjs;
        LDObjKey *objs = cdata.receiver.get_destObjs(nobjs);
        int senderID = stats->getHash(cdata.sender);
        for (j=0; j<nobjs; j++) {
           int recverID = stats->getHash(objs[j]);
           if((senderID == -1)||(recverID == -1))
              if (_lb_args.migObjOnly()) continue;
              else CkAbort("Error in search\n");
           comm[senderID][recverID] += cdata.messages;
           comm[recverID][senderID] += cdata.messages;
        }
      }
      */
    }

    delete[] lbcentroids;


}


#include "CommMeasureLB.def.h"

/*@}*/
