//Abhishek - Extract common code from Orb and mslb_notopo
#ifndef _ORB3DHELPER_H
#define _ORB3DHELPER_H

#include <charm++.h>
#include <charm++.h>
#include "cklists.h"
#include "ParallelGravity.h"
#include "TopoManager.h"

#include "Refiner.h"
#include "MapStructures.h"
#include "TaggedVector3D.h"
#include "Vector3D.h"
#include "CentralLB.h"
#define  ORB3DLB_NOTOPO_DEBUG
// #define  ORB3DLB_NOTOPO_DEBUG CkPrintf

class PeInfo {
  public:
  int idx;
  double load;
  double items;
  PeInfo(int id, double ld, int it) : idx(id), load(ld), items(it) {}
};

class ProcLdGreater {
  public:
  bool operator()(PeInfo& p1, PeInfo& p2) {
    // This can be done based on load or number of tps assigned to a PE
    return (p1.load > p2.load);
  }
};

class Orb3dCommon{
  // pointer to stats->to_proc
  protected:		
    CkVec<int> *mapping;
    CkVec<int> *from;

    CkVec<float> procload;

    /// Take into account memory constraints by limiting the number of pieces
    /// per processor.
    double maxPieceProc;

    /// index of first processor of the group we are considering
    int nextProc;

    // Greedy strategy to assign TreePieces to PEs on a node.
    void orbPePartition(vector<Event> *events, vector<OrbObject> &tp, int node,
        BaseLB::LDStats *stats) {

      std::vector<PeInfo> peinfo;
      float totalLoad = 0.0;
      int firstProc = CkNodeFirst(node);
      int lastProc = firstProc + CkNodeSize(node) - 1;
      for (int i = firstProc; i <= lastProc; i++) {
        peinfo.push_back(PeInfo(i, 0.0, 0));
      }
      // Make a heap of processors belonging to this node
      std::make_heap(peinfo.begin(), peinfo.end(), ProcLdGreater());

      int nextProc;
      for(int i = 0; i < events[XDIM].size(); i++){
        Event &ev = events[XDIM][i];
        OrbObject &orb = tp[ev.owner];

        // Pop the least loaded PE from the heap and assign TreePiece to it
        PeInfo p = peinfo.front();
        pop_heap(peinfo.begin(), peinfo.end(), ProcLdGreater());
        peinfo.pop_back();

        nextProc = p.idx;

        if(orb.numParticles > 0){
          (*mapping)[orb.lbindex] = nextProc;
          procload[nextProc] += ev.load;
          p.load += ev.load;
          p.items += 1;
          totalLoad += ev.load;
        } else{
          int fromPE = (*from)[orb.lbindex];
          procload[fromPE] += ev.load;
        }

        peinfo.push_back(p);
        push_heap(peinfo.begin(), peinfo.end(), ProcLdGreater());
      }
    }

    void orbPartition(vector<Event> *events, OrientedBox<float> &box, int nprocs,
        vector<OrbObject> & tp, BaseLB::LDStats *stats,
        bool node_partition=false){

      ORB3DLB_NOTOPO_DEBUG("partition events %d %d %d nprocs %d\n", 
          events[XDIM].size(),
          events[YDIM].size(),
          events[ZDIM].size(),
          nprocs
          );
      int numEvents = events[XDIM].size();
      CkAssert(numEvents == events[YDIM].size());
      CkAssert(numEvents == events[ZDIM].size());

      if(numEvents == 0)
	return;

      if(nprocs == 1){
        ORB3DLB_NOTOPO_DEBUG("base: assign %d tps to proc %d\n", numEvents, nextProc);
        if (!stats->procs[nextProc].available) {
          nextProc++;
          return;
        }

        // If we are doing orb partition at the node level, then call
        // orbPePartition to assign the treepieces to the PEs belonging to the node.
        if (node_partition) {
          orbPePartition(events, tp, nextProc, stats);
        } else {
          // direct assignment of tree pieces to processors
          //if(numEvents > 0) CkAssert(nprocs != 0);
          float totalLoad = 0.0;
          for(int i = 0; i < events[XDIM].size(); i++){
            Event &ev = events[XDIM][i];
            OrbObject &orb = tp[ev.owner];
            if(orb.numParticles > 0){
              (*mapping)[orb.lbindex] = nextProc;
              totalLoad += ev.load;
            }
            else{
              int fromPE = (*from)[orb.lbindex];
              if (fromPE < 0 || fromPE >= procload.size()) {
                CkPrintf("[%d] trying to access fromPe %d nprocs %d\n", CkMyPe(), fromPE, procload.size());
                CkAbort("Trying to access a PE which is outside the range\n");
              }
              procload[fromPE] += ev.load;
            }
          }
          procload[nextProc] += totalLoad;
        }

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

      float ratio = (1.0*nlprocs)/(1.0*(nlprocs+nrprocs));

      // sum background load on each side of the processor split
      float bglprocs = 0.0;
      for(int np = nextProc; np < nextProc + nlprocs; np++)
        bglprocs += stats->procs[np].bg_walltime;
      float bgrprocs = 0.0;
      for(int np = nextProc + nlprocs; np < nextProc + nlprocs + nrprocs; np++)
        bgrprocs += stats->procs[np].bg_walltime;

      ORB3DLB_NOTOPO_DEBUG("nlprocs %d nrprocs %d ratio %f\n", nlprocs, nrprocs, ratio);

      int splitIndex = partitionRatioLoad(events[longestDim],ratio,bglprocs,
                                          bgrprocs);
      if(splitIndex == numEvents) {
        ORB3DLB_NOTOPO_DEBUG("evenly split 0 load\n");
        splitIndex = splitIndex/2;
      }
      int nleft = splitIndex;
      int nright = numEvents-nleft;

#if 0
      if(nright < nrprocs) {  // at least one piece per processor
        nright = nrprocs;
        nleft = splitIndex = numEvents-nright;
        CkAssert(nleft >= nlprocs);
      }
      else if(nleft < nlprocs) {
        nleft = splitIndex = nlprocs;
        nright = numEvents-nleft;
        CkAssert(nright >= nrprocs);
      }
#endif

      if(nleft > nlprocs*maxPieceProc) {
	  nleft = splitIndex = (int) (nlprocs*maxPieceProc);
	  nright = numEvents-nleft;
	  }
      else if (nright > nrprocs*maxPieceProc) {
	  nright = (int) (nrprocs*maxPieceProc);
	  nleft = splitIndex = numEvents-nright;
	  }
      CkAssert(splitIndex >= 0);
      CkAssert(splitIndex < numEvents);

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
        CkAssert(tp[ev.owner].partition == INVALID_PARTITION);
        tp[ev.owner].partition = LEFT_PARTITION;
      }
      for(int i = splitIndex; i < numEvents; i++){
        Event &ev = events[longestDim][i];
        CkAssert(ev.owner >= 0);
        CkAssert(tp[ev.owner].partition == INVALID_PARTITION);
        tp[ev.owner].partition = RIGHT_PARTITION;
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
          OrbObject &orb = tp[ev.owner];
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
        OrbObject &orb = tp[ev.owner];
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
      orbPartition(leftEvents,leftBox,nlprocs,tp, stats);
      orbPartition(rightEvents,rightBox,nrprocs,tp, stats);
    }

    void orbPrepare(vector<Event> *tpEvents, OrientedBox<float> &box, int
    numobjs, BaseLB::LDStats * stats, bool node_partition=false){

      int nmig = stats->n_migrateobjs;
      if(dMaxBalance < 1.0)
        dMaxBalance = 1.0;

      // If using node based orb partition, then the maxPieceProc is total
      // migratable objs / total number of node.
      if (node_partition) {
        maxPieceProc = dMaxBalance * nmig / CkNumNodes();
      } else {
        maxPieceProc = dMaxBalance*nmig/stats->count;
      }

      if(maxPieceProc < 1.0)
        maxPieceProc = 1.01;

      CkAssert(tpEvents[XDIM].size() == numobjs);
      CkAssert(tpEvents[YDIM].size() == numobjs);
      CkAssert(tpEvents[ZDIM].size() == numobjs);

      mapping = &stats->to_proc;
      from = &stats->from_proc;

      CkPrintf("[Orb3dLB_notopo] sorting\n");
      for(int i = 0; i < NDIMS; i++){
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
      for(int i = 0; i < stats->count; i++){
        procload[i] = stats->procs[i].bg_walltime;
      }

    }

    void refine(BaseLB::LDStats *stats, int numobjs){
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

      double *predLoad = new double[stats->count];
      double *predCount = new double[stats->count];
      for(int i = 0; i < stats->count; i++){
        predLoad[i] = 0.0;
        predCount[i] = 0.0;
      }

      double maxObjLoad = 0.0;
      
      for(int i = 0; i < numobjs; i++){
        double ld = stats->objData[i].wallTime;
        int proc = stats->to_proc[i];
        predLoad[proc] += ld; 
        predCount[proc] += 1.0; 
        if(ld > maxObjLoad)
            maxObjLoad = ld;
      }

      double minWall = 0.0;
      double maxWall = 0.0;
      double avgWall = 0.0;

      double minIdle = 0.0;
      double maxIdle = 0.0;
      double avgIdle = 0.0;

      double minBg = 0.0;
      double maxBg = 0.0;
      double avgBg = 0.0;

      double avgPred = 0.0;
      double minPred = 0.0;
      double maxPred = 0.0;

      double avgPiece = 0.0;
      double minPiece = 0.0;
      double maxPiece = 0.0;

      CkPrintf("***************************\n");
      //  CkPrintf("Before LB step %d\n", step());
      for(int i = 0; i < stats->count; i++){
        double wallTime = stats->procs[i].total_walltime;
        double idleTime = stats->procs[i].idletime;
        double bgTime = stats->procs[i].bg_walltime;
        double pred = predLoad[i];
        double npiece = predCount[i];
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
        avgPred += pred;
        avgPiece += npiece;

        if(i==0 || minWall > wallTime) minWall = wallTime;
        if(i==0 || maxWall < wallTime) maxWall = wallTime;

        if(i==0 || minIdle > idleTime) minIdle = idleTime;
        if(i==0 || maxIdle < idleTime) maxIdle = idleTime;

        if(i==0 || minBg > bgTime) minBg = bgTime;
        if(i==0 || maxBg < bgTime) maxBg = bgTime;

        if(i==0 || minPred > pred) minPred = pred;
        if(i==0 || maxPred < pred) maxPred = pred;

        if(i==0 || minPiece > npiece) minPiece = npiece;
        if(i==0 || maxPiece < npiece) maxPiece = npiece;

      }

      avgWall /= stats->count;
      avgIdle /= stats->count;
      avgBg /= stats->count;
      avgPred /= stats->count;
      avgPiece /= stats->count;

#ifdef PRINT_LOAD_PERCENTILES
      double accumVar = 0;
      vector<double> objectWallTimes;
      for(int i = 0; i < stats->count; i++){
        double wallTime = stats->procs[i].total_walltime;
        objectWallTimes.push_back(wallTime);
        accumVar += (wallTime - avgWall) * (wallTime - avgWall);
      }
      double stdDev = sqrt(accumVar / stats->count);
      CkPrintf("Average load: %.3f\n", avgWall);
      CkPrintf("Standard deviation: %.3f\n", stdDev);

      std::sort(objectWallTimes.begin(), objectWallTimes.end());
      CkPrintf("Object load percentiles: \n");
      double increment = (double) objectWallTimes.size() / 10;
      int j = 0;
      double index = 0;
      for (int j = 0; j < 100; j += 10) {
        index += increment;
        CkPrintf("%d: %.3f\n", j, objectWallTimes[(int) index]);
      }
      CkPrintf("100: %.3f\n", objectWallTimes.back());
#endif

      delete[] predLoad;
      delete[] predCount;

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


      CkPrintf("Orb3dLB_notopo stats: maxObjLoad %f\n", maxObjLoad);
      CkPrintf("Orb3dLB_notopo stats: minWall %f maxWall %f avgWall %f maxWall/avgWall %f\n", minWall, maxWall, avgWall, maxWall/avgWall);
      CkPrintf("Orb3dLB_notopo stats: minIdle %f maxIdle %f avgIdle %f minIdle/avgIdle %f\n", minIdle, maxIdle, avgIdle, minIdle/avgIdle);
      CkPrintf("Orb3dLB_notopo stats: minPred %f maxPred %f avgPred %f maxPred/avgPred %f\n", minPred, maxPred, avgPred, maxPred/avgPred);
      CkPrintf("Orb3dLB_notopo stats: minPiece %f maxPiece %f avgPiece %f maxPiece/avgPiece %f\n", minPiece, maxPiece, avgPiece, maxPiece/avgPiece);

      CkPrintf("Orb3dLB_notopo stats: minBg %f maxBg %f avgBg %f maxBg/avgBg %f\n", minBg, maxBg, avgBg, maxBg/avgBg);
      CkPrintf("Orb3dLB_notopo stats: orb migrated %d refine migrated %d objects\n", migr, numRefineMigrated);

#ifdef DO_REFINE
      // Free the refine buffers
      Refiner::FreeProcs(from_procs);
      Refiner::FreeProcs(to_procs);
#endif

    }

#define LOAD_EQUAL_TOLERANCE 1.02
    int partitionRatioLoad(vector<Event> &events, float ratio, float bglp, float bgrp){

      float approxBgPerEvent = (bglp + bgrp) / events.size();
      float totalLoad = bglp + bgrp;
      for(int i = 0; i < events.size(); i++){
        totalLoad += events[i].load;
      }
      //CkPrintf("************************************************************\n");
      //CkPrintf("partitionEvenLoad start %d end %d total %f\n", tpstart, tpend, totalLoad);
      float perfectLoad = ratio * totalLoad;
      ORB3DLB_NOTOPO_DEBUG("partitionRatioLoad bgl %f bgr %f\n",
                           bglp, bgrp);
      int splitIndex = 0;
      float prevLoad = 0.0;
      float leftLoadAtSplit = 0.0;
      for(splitIndex = 0; splitIndex < events.size(); splitIndex++){

        leftLoadAtSplit += events[splitIndex].load + approxBgPerEvent;

        if (leftLoadAtSplit > perfectLoad) {
          if ( fabs(leftLoadAtSplit - perfectLoad) < fabs(prevLoad - perfectLoad) ) {
            splitIndex++;
          }
          else {
            leftLoadAtSplit = prevLoad;
          }
          break;
        }
        prevLoad = leftLoadAtSplit;
      }

      ORB3DLB_NOTOPO_DEBUG("partitionEvenLoad mid %d lload %f rload %f ratio %f\n", splitIndex, leftLoadAtSplit, totalLoad - leftLoadAtSplit, leftLoadAtSplit / totalLoad);
      return splitIndex;
    }

}; //end class

#endif
