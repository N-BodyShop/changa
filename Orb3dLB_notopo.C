#include <charm++.h>
#include "cklists.h"
#include "Orb3dLB_notopo.h"
#include "ParallelGravity.h"
#include "Refiner.h"
#include "TopoManager.h"
#include "Vector3D.h"

extern CProxy_TreePiece treeProxy;

#define ORB3DLB_NOTOPO_DEBUG 

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

  CkAssert(numProcessed == numobjs);
  CkAssert(tpEvents[XDIM].size() == numobjs);
  CkAssert(tpEvents[YDIM].size() == numobjs);
  CkAssert(tpEvents[ZDIM].size() == numobjs);

  mapping = &stats->to_proc;
  from = &stats->from_proc;
  int dim = 0;

  CkPrintf("[Orb3dLB_notopo] sorting\n");
  for(int i = 0; i < NDIMS; i++){
    //tpEvents[i].quickSort();
    sort(tpEvents[i].begin(),tpEvents[i].end());
  }

  box.lesser_corner.x = tpEvents[XDIM][0].position;
  box.lesser_corner.y = tpEvents[YDIM][0].position;
  box.lesser_corner.z = tpEvents[ZDIM][0].position;

  box.greater_corner.x = tpEvents[XDIM][numobjs-1].position;
  box.greater_corner.y = tpEvents[YDIM][numobjs-1].position;
  box.greater_corner.z = tpEvents[ZDIM][numobjs-1].position;

  nextProc = 0;

  procload.resize(stats->count);
  procbox.resize(stats->count);
  for(int i = 0; i < stats->count; i++){
    procload[i] = 0.0;
  }

  orbPartition(tpEvents,box,stats->count);

#ifdef DO_REFINE
  int *from_procs = Refiner::AllocProcs(stats->count, stats);
  int *to_procs = Refiner::AllocProcs(stats->count, stats);
#endif

  int migr = 0;
  for(int i = 0; i < numobjs; i++){
    if(stats->to_proc[i] != stats->from_proc[i]) migr++;
#ifdef DO_REFINE
    int pe = stats->to_proc[i];
    from_procs[i] = pe;
    to_procs[i] = pe;
#endif
  }

  int numRefineMigrated = 0;
#ifdef DO_REFINE
  CkPrintf("[orb3dlb_notopo] refine\n");
  Refiner refiner(1.050);
  refiner.Refine(stats->count,stats,from_procs,to_procs);

  for(int i = 0; i < numobjs; i++){
    if(to_procs[i] != from_procs[i]) numRefineMigrated++;
    stats->to_proc[i] = to_procs[i];
  }
#endif

  double minWall = 0.0;
  double maxWall = 0.0;
  double avgWall = 0.0;

  double minIdle = 0.0;
  double maxIdle = 0.0;
  double avgIdle = 0.0;

  double minBg = 0.0;
  double maxBg = 0.0;
  double avgBg = 0.0;

  CkPrintf("***************************\n");
  CkPrintf("Before LB step %d\n", step());
  CkPrintf("***************************\n");
  CkPrintf("i pe wall idle bg_wall objload\n");
  for(int i = 0; i < stats->count; i++){
    double wallTime = stats->procs[i].total_walltime;
    double idleTime = stats->procs[i].idletime;
    double bgTime = stats->procs[i].bg_walltime;
    /*
    CkPrintf("[pestats] %d %d %f %f %f %f\n", 
        i,
        stats->procs[i].pe, 
        wallTime,
        idleTime,
        bgTime,
        objTime);
        */

    avgWall += wallTime; 
    avgIdle += idleTime; 
    avgBg += bgTime;

    if(i==0 || minWall > wallTime) minWall = wallTime;
    if(i==0 || maxWall < wallTime) maxWall = wallTime;

    if(i==0 || minIdle > idleTime) minIdle = idleTime;
    if(i==0 || maxIdle < idleTime) maxIdle = idleTime;

    if(i==0 || minBg > bgTime) minBg = bgTime;
    if(i==0 || maxBg < bgTime) maxBg = bgTime;

  }

  avgWall /= stats->count;
  avgIdle /= stats->count;
  avgBg /= stats->count;

#if 0
  float minload, maxload, avgload;
  minload = maxload = procload[0];
  avgload = 0.0;
  for(int i = 0; i < stats->count; i++){
    CkPrintf("pe %d load %f box %f %f %f %f %f %f\n", i, procload[i], 
                                procbox[i].lesser_corner.x,
                                procbox[i].lesser_corner.y,
                                procbox[i].lesser_corner.z,
                                procbox[i].greater_corner.x,
                                procbox[i].greater_corner.y,
                                procbox[i].greater_corner.z
                                );
    avgload += procload[i];
    if(minload > procload[i]) minload = procload[i];
    if(maxload < procload[i]) maxload = procload[i];
  }

  avgload /= stats->count;

  CkPrintf("Orb3dLB_notopo stats: min %f max %f avg %f max/avg %f\n", minload, maxload, avgload, maxload/avgload);
#endif

 
  CkPrintf("Orb3dLB_notopo stats: minWall %f maxWall %f avgWall %f maxWall/avgWall %f\n", minWall, maxWall, avgWall, maxWall/avgWall);
  CkPrintf("Orb3dLB_notopo stats: minIdle %f maxIdle %f avgIdle %f maxIdle/avgIdle %f\n", minIdle, maxIdle, avgIdle, maxIdle/avgIdle);
  CkPrintf("Orb3dLB_notopo stats: minBg %f maxBg %f avgBg %f maxBg/avgBg %f\n", minBg, maxBg, avgBg, maxBg/avgBg);
  CkPrintf("Orb3dLB_notopo stats: orb migrated %d refine migrated %d objects\n", migr, numRefineMigrated);

#ifdef DO_REFINE
  // Free the refine buffers
  Refiner::FreeProcs(from_procs);
  Refiner::FreeProcs(to_procs);
#endif

}

void Orb3dLB_notopo::orbPartition(vector<Event> *events, OrientedBox<float> &box, int nprocs){

  ORB3DLB_NOTOPO_DEBUG("partition events %d %d %d nprocs %d\n", 
            events[XDIM].size(),
            events[YDIM].size(),
            events[ZDIM].size(),
            nprocs
            );
  int numEvents = events[XDIM].size();
  CkAssert(numEvents == events[YDIM].size());
  CkAssert(numEvents == events[ZDIM].size());

  if(nprocs == 1){
    ORB3DLB_NOTOPO_DEBUG("base: assign %d tps to proc %d\n", numEvents, nextProc);
    // direct assignment of tree pieces to processors
    //if(numEvents > 0) CkAssert(nprocs != 0);
    float totalLoad = 0.0;
    for(int i = 0; i < events[XDIM].size(); i++){
      Event &ev = events[XDIM][i];
      OrbObject &orb = tps[ev.owner];
      if(orb.numParticles > 0){
        (*mapping)[orb.lbindex] = nextProc;
        totalLoad += ev.load;
        procbox[nextProc].grow(orb.centroid);
      }
      else{
        int fromPE = (*from)[orb.lbindex];
        procload[fromPE] += ev.load;
      }
    }
    procload[nextProc] += totalLoad;

    if(numEvents > 0) nextProc++;
    return;
  }

  // find longest dimension

  int longestDim = XDIM;
  float longestDimLength = box.greater_corner[longestDim] - box.lesser_corner[longestDim];
  for(int i = YDIM; i <= ZDIM; i++){
    float thisDimLength = box.greater_corner[i]-box.lesser_corner[i];
    if(thisDimLength > longestDimLength){
      longestDimLength = thisDimLength;
      longestDim = i;
    }
  }

  ORB3DLB_NOTOPO_DEBUG("dimensions %f %f %f longest %d\n", 
            box.greater_corner[XDIM]-box.lesser_corner[XDIM],
            box.greater_corner[YDIM]-box.lesser_corner[YDIM],
            box.greater_corner[ZDIM]-box.lesser_corner[ZDIM],
            longestDim
          );

  int nlprocs = nprocs/2;
  int nrprocs = nprocs-nlprocs;

  float ratio = (1.0*nlprocs)/(1.0*nrprocs);

  ORB3DLB_NOTOPO_DEBUG("nlprocs %d nrprocs %d ratio %f\n", nlprocs, nrprocs, ratio);

  int splitIndex = partitionRatioLoad(events[longestDim],ratio);
  int nleft = splitIndex;
  int nright = numEvents-nleft;

  OrientedBox<float> leftBox;
  OrientedBox<float> rightBox;

  leftBox = rightBox = box;
  float splitPosition = events[longestDim][splitIndex].position;
  leftBox.greater_corner[longestDim] = splitPosition;
  rightBox.lesser_corner[longestDim] = splitPosition;

  // classify events
  for(int i = 0; i < splitIndex; i++){
    Event &ev = events[longestDim][i];
    CkAssert(ev.owner >= 0);
    CkAssert(tps[ev.owner].partition == INVALID_PARTITION);
    tps[ev.owner].partition = LEFT_PARTITION;
  }
  for(int i = splitIndex; i < numEvents; i++){
    Event &ev = events[longestDim][i];
    CkAssert(ev.owner >= 0);
    CkAssert(tps[ev.owner].partition == INVALID_PARTITION);
    tps[ev.owner].partition = RIGHT_PARTITION;
  }

  vector<Event> leftEvents[NDIMS];
  vector<Event> rightEvents[NDIMS];

  for(int i = 0; i < NDIMS; i++){
    if(i == longestDim){ 
      leftEvents[i].resize(nleft);
      rightEvents[i].resize(nright);
    }
    else{
      leftEvents[i].reserve(nleft);
      rightEvents[i].reserve(nright);
    }
  }

  // copy events of split dimension
  memcpy(&leftEvents[longestDim][0],&events[longestDim][0],sizeof(Event)*nleft);
  memcpy(&rightEvents[longestDim][0],&events[longestDim][splitIndex],sizeof(Event)*nright);
  
  // copy events of other dimensions
  for(int i = XDIM; i <= ZDIM; i++){
    if(i == longestDim) continue;
    for(int j = 0; j < numEvents; j++){
      Event &ev = events[i][j];
      CkAssert(ev.owner >= 0);
      OrbObject &orb = tps[ev.owner];
      CkAssert(orb.partition != INVALID_PARTITION);
      if(orb.partition == LEFT_PARTITION) leftEvents[i].push_back(ev);
      else if(orb.partition == RIGHT_PARTITION) rightEvents[i].push_back(ev);
    }
  }

  // cleanup
  // next, reset the ownership information in the
  // OrbObjects, so that the next invocation may use
  // the same locations for its book-keeping
  vector<Event> &eraseVec = events[longestDim];
  for(int i = 0; i < numEvents; i++){
    Event &ev = eraseVec[i];
    CkAssert(ev.owner >= 0);
    OrbObject &orb = tps[ev.owner];
    CkAssert(orb.partition != INVALID_PARTITION);
    orb.partition = INVALID_PARTITION;
  }

  // free events from parent node,
  // since they are not needed anymore
  // (we have partition all events into the
  // left and right event subsets)
  for(int i = 0; i < NDIMS; i++){
    //events[i].free();
    vector<Event>().swap(events[i]);
  }

  orbPartition(leftEvents,leftBox,nlprocs);
  orbPartition(rightEvents,rightBox,nrprocs);
}

int Orb3dLB_notopo::partitionRatioLoad(vector<Event> &events, float ratio){
  float totalLoad = 0.0;
  for(int i = 0; i < events.size(); i++){
    totalLoad += events[i].load;
  }
  //CkPrintf("************************************************************\n");
  //CkPrintf("partitionEvenLoad start %d end %d total %f\n", tpstart, tpend, totalLoad);
  float lload = 0.0;
  float rload = totalLoad;
  float prevDiff = lload-ratio*rload;
  if(prevDiff < 0.0){
    prevDiff = -prevDiff;
  }

  int consider;
  for(consider = 0; consider < events.size();){
    float newll = lload + events[consider].load;
    float newrl = rload - events[consider].load;

    float newdiff = newll-ratio*newrl;
    if(newdiff < 0.0){
      newdiff = -newdiff;
    }

    ORB3DLB_NOTOPO_DEBUG("consider load %f newdiff %f prevdiff %f\n", events[consider].load, newdiff, prevDiff);

    if(newdiff > prevDiff){
      break;
    }
    else{
      consider++;
      lload = newll;
      rload = newrl;
      prevDiff = newdiff;
    }
  }

  ORB3DLB_NOTOPO_DEBUG("partitionEvenLoad mid %d lload %f rload %f ratio %f\n", consider, lload, rload, lload/rload);
  return consider;
}

#include "Orb3dLB_notopo.def.h"

/*@}*/
