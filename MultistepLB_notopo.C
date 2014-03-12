#include <charm++.h>
#include "cklists.h"
#include "MultistepLB_notopo.h"
#include "TopoManager.h"
#include "ParallelGravity.h"
#include "Vector3D.h"
#include <queue>
#include "CkLoopAPI.h"

extern CProxy_TreePiece treeProxy;
using namespace std;
//#define ORB3DLB_NOTOPO_DEBUG CkPrintf

CreateLBFunc_Def(MultistepLB_notopo, "Works best with multistepped runs; uses Orb3D_notopo for larger steps, greedy otherwise");


MultistepLB_notopo::MultistepLB_notopo(const CkLBOptions &opt): CentralLB(opt)
{
  lbname = "MultistepLB_notopo";

  if (CkMyPe() == 0){
    CkPrintf("[%d] MultistepLB_notopo created\n",CkMyPe());
  }

  
  haveTPCentroids = false;
  clearPeLoad();
}

/// @brief Get position centroids of all TreePieces
/// @param msg Reduction message with a concatenation of all centroids.
void MultistepLB_notopo::receiveCentroids(CkReductionMsg *msg){
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



bool MultistepLB_notopo::QueryBalanceNow(int step){
 if(CkMyPe() == 0) CkPrintf("Orb3dLB_notopo: Step %d\n", step);
 //  if(step == 0) return false;
  return true;

}

// helper functions for multistepping
#ifdef MCLBMS

// determine phase based on lastActiveRung as saved in map.tpCentroids
unsigned int MultistepLB_notopo::determinePhase(unsigned int lastActiveRung){
  return lastActiveRung;
}

// merge data instrumented in previous iteration of computation with data from earlier iterations
// this data must be stored according to the tree piece 
// index of the treepiece (obtained from treepiece[.].tag)
// since the ldb index of an object changes from iteration
// to iteration.
void MultistepLB_notopo::mergeInstrumentedData(int phase, BaseLB::LDStats *stats){

  int i, len;
  int whichPos;
  int numAdditional;

  if(phase == -1){
#ifdef MCLBMSV
    CkPrintf("phase = -1, discarding\n");
#endif
    //return;
    phase = 0;
  }

  if (_lb_args.debug()>=2) {
    CkPrintf("**********************************************\n");
    CkPrintf("Actual object loads phase %d\n", phase);
    CkPrintf("**********************************************\n");
    for(int i = 0; i < stats->n_objs; i++){
      int tp = tpCentroids[i].tp;
      int lb = tpCentroids[i].tag;
      CkPrintf("tp %d load %f\n",tp,stats->objData[lb].wallTime);
    }

    CkPrintf("Processor stats\n");
    for(int i = 0; i < stats->count; i++){
      CkPrintf("Proc %d with %d tps: wall %g idle %g busy %g bg_wall %g\n",
          stats->procs[i].pe,
          stats->procs[i].n_objs, stats->procs[i].total_walltime,
          stats->procs[i].idletime,
          stats->procs[i].total_walltime-stats->procs[i].idletime,
          stats->procs[i].bg_walltime);
    }
  }
  int whichPeMaxLd;
  double peMaxLd = 0.0;
  if (step() > 0) {
    for(int i = 0; i < stats->count; i++){
      //if (step() == 9) {
      //  CkPrintf("Proc %d with %d tps: wall %g idle %g busy %g bg_wall %g\n",
      //      stats->procs[i].pe,
      //      stats->procs[i].n_objs, stats->procs[i].total_walltime,
      //      stats->procs[i].idletime,
      //      stats->procs[i].total_walltime-stats->procs[i].idletime,
      //      stats->procs[i].bg_walltime);
      //}
      if (peMaxLd < 
          stats->procs[i].total_walltime-stats->procs[i].idletime) {
        peMaxLd = 
          stats->procs[i].total_walltime-stats->procs[i].idletime;
        whichPeMaxLd = i;
      }
    }


    CkPrintf("MaxLoadedPe for prevstep is %d with ld %f\n", whichPeMaxLd,
        peMaxLd);
    for(int i = 0; i < stats->n_objs; i++){
      int tp = tpCentroids[i].tp;
      int lb = tpCentroids[i].tag;
      double ld = tpCentroids[i].busytime;
      if (stats->from_proc[lb] == whichPeMaxLd) {
        CkPrintf("\t TP %d on %d PrevPredLd %f ActualLD %f Particles %d\n",
            tp, whichPeMaxLd, objDataPrevPred[tp], ld, objDataPrevPredPart[tp]);
            //stats->objData[lb].wallTime, objDataPrevPredPart[tp]);
      }

      if (tp == 239226) {
        CkPrintf("\t Sidelined!!! TP %d PrevPredLd %f ActualLD %f Particles %d\n\n",
            tp, objDataPrevPred[tp], ld, objDataPrevPredPart[tp]);
      }

    }

    for(int i = 0; i < stats->n_objs; i++){
      int tp = tpCentroids[i].tp;
      int lb = tpCentroids[i].tag;
      double ld = tpCentroids[i].busytime;

      if (stats->from_proc[lb] == prevMaxPredPe) {
        CkPrintf("\t TP %d on PrevMaxPred %d PrevPredLd %f ActualLD %f Particles %d\n",
            tp, prevMaxPredPe, objDataPrevPred[tp], ld, objDataPrevPredPart[tp]);
      }
    }

  }


  // if (step() == 6) {
  //   objDataPrevPred.resize(stats->n_objs);
  //   objDataPrevPredPart.resize(stats->n_objs);
  // }

  len = savedPhaseStats.length();

  if(phase > len-1){
    numAdditional = phase-len+1;
    while(numAdditional > 0){
      savedPhaseStats.push_back(LightweightLDStats());
#ifdef MCLBMSV
      CkPrintf("Making new entry for phase %d (%d)\n", savedPhaseStats.length()-1, phase);
#endif
      savedPhaseStats[savedPhaseStats.length()-1].n_objs = 0;  // so that we can identify this phase's entry as being unfilled when we get the data for it
      numAdditional--;
    }
    len = savedPhaseStats.length();
    savedPhaseStats[len-1].objData.resize(stats->n_objs);
    savedPhaseStats[len-1].n_activeparticles.resize(stats->n_objs);
    for(i = 0; i < stats->n_objs; i++){
      int tp = tpCentroids[i].tp;
      int lb = tpCentroids[i].tag;
      savedPhaseStats[len-1].objData[tp] = stats->objData[lb];
      if (prevNumActiveParts.length() > 0) {
        savedPhaseStats[len-1].n_activeparticles[tp] = prevNumActiveParts[tp];
      } else {
        savedPhaseStats[len-1].n_activeparticles[tp] =
          tpCentroids[i].myNumParticles;
      }
    }
    whichPos = len-1;
  }
  else{ 
    whichPos = phase;
    if(savedPhaseStats[phase].n_objs == 0){       // haven't yet populated this phase
#ifdef MCLBMSV
      CkPrintf("Found unpopulated entry for phase %d\n", phase);
#endif
      savedPhaseStats[phase].objData.resize(stats->n_objs);
      savedPhaseStats[phase].n_activeparticles.resize(stats->n_objs);
      for(i = 0; i < stats->n_objs; i++){
        int tp = tpCentroids[i].tp;
        int lb = tpCentroids[i].tag;
        savedPhaseStats[phase].objData[tp] = stats->objData[lb];
        if (prevNumActiveParts.length() > 0) {
          savedPhaseStats[phase].n_activeparticles[tp] = prevNumActiveParts[tp];
        } else {
          savedPhaseStats[phase].n_activeparticles[tp] =
            tpCentroids[i].myNumParticles;
        }

      }
    }
    else{        // filled this phase out some time in the past - merge current with previous data
#ifdef MCLBMSV
      CkPrintf("Found previous entry for phase %d - merging\n", phase);
#endif
      for(i = 0; i < stats->n_objs; i++){
        int tp = tpCentroids[i].tp;
        int lb = tpCentroids[i].tag;
        savedPhaseStats[phase].objData[tp] = stats->objData[lb];
        if (prevNumActiveParts.length() > 0) {
          savedPhaseStats[phase].n_activeparticles[tp] = prevNumActiveParts[tp];
        } else {
          savedPhaseStats[phase].n_activeparticles[tp] =
            tpCentroids[i].myNumParticles;
        }

      }
    }
  }
  savedPhaseStats[whichPos].n_objs = stats->n_objs;
  savedPhaseStats[whichPos].n_migrateobjs = stats->n_migrateobjs;
#ifdef MCLBMSV
  //printData(savedPhaseStats[whichPos], phase, NULL);
#endif
}

void MultistepLB_notopo::saveStatsData(int phase, BaseLB::LDStats *stats){

  int i, len;
  int whichPos;
  int numAdditional;

  if(phase == -1){
#ifdef MCLBMSV
    CkPrintf("phase = -1, discarding\n");
#endif
    //return;
    phase = 0;
  }

  int whichPeMaxLd;
  double peMaxLd = 0.0;
  if (step() > 0) {
    for(int i = 0; i < stats->count; i++){
      //if (step() == 9) {
      //  CkPrintf("Proc %d with %d tps: wall %g idle %g busy %g bg_wall %g\n",
      //      stats->procs[i].pe,
      //      stats->procs[i].n_objs, stats->procs[i].total_walltime,
      //      stats->procs[i].idletime,
      //      stats->procs[i].total_walltime-stats->procs[i].idletime,
      //      stats->procs[i].bg_walltime);
      //}
      if (peMaxLd < 
          stats->procs[i].total_walltime-stats->procs[i].idletime) {
        peMaxLd = 
          stats->procs[i].total_walltime-stats->procs[i].idletime;
        whichPeMaxLd = i;
      }
    }


    CkPrintf("MaxLoadedPe for prevstep is %d with ld %f\n", whichPeMaxLd,
        peMaxLd);
    for(int i = 0; i < stats->n_objs; i++){
      int tp = tpCentroids[i].tp;
      int lb = tpCentroids[i].tag;
      double ld = tpCentroids[i].busytime;
      if (stats->from_proc[lb] == whichPeMaxLd) {
        CkPrintf("\t TP %d on %d PrevPredLd %f ActualLD %f Particles %d\n",
            tp, whichPeMaxLd, objDataPrevPred[tp], ld, objDataPrevPredPart[tp]);
      }

      if (tp == 239226) {
        CkPrintf("\t Sidelined!!! TP %d PrevPredLd %f ActualLD %f Particles %d\n\n",
            tp, objDataPrevPred[tp], ld, objDataPrevPredPart[tp]);
      }

    }

    for(int i = 0; i < stats->n_objs; i++){
      int tp = tpCentroids[i].tp;
      int lb = tpCentroids[i].tag;
      double ld = tpCentroids[i].busytime;

      if (stats->from_proc[lb] == prevMaxPredPe) {
        CkPrintf("\t TP %d on PrevMaxPred %d PrevPredLd %f ActualLD %f Particles %d\n",
            tp, prevMaxPredPe, objDataPrevPred[tp], ld, objDataPrevPredPart[tp]);
      }
    }
  }
}

// whether we have instrumented data for this phase
bool MultistepLB_notopo::havePhaseData(int phase){
  return (savedPhaseStats.length() > phase && savedPhaseStats[phase].n_objs > 0);
}

void MultistepLB_notopo::printData(BaseLB::LDStats &stats, int phase, int *revObjMap){
  int i;
  
  CkPrintf("---- data (%d): %d objects ----\n", phase, stats.n_objs);
  for(i = 0; i < stats.n_objs; i++){
     CkPrintf("%d: %g\n", i, 
	       stats.objData[i].wallTime);
  }
  CkPrintf("---- end data (%d) ----\n", phase);
}

void MultistepLB_notopo::makeActiveProcessorList(BaseLB::LDStats *stats, int numActiveObjs){
  int objsPerProc = 8;
  int expandFactor = 4;
  int procsNeeded;
  procsNeeded = expandFactor*numActiveObjs/objsPerProc > stats->count ? stats->count : expandFactor*numActiveObjs/objsPerProc;
  CkPrintf("Processors 0 to %d active\n", procsNeeded-1);

  /* currently, only the first procsNeeded procs are used - could do something more sophisticated here in the future - FIXME */
#ifdef MCLBMSV
#endif
}
#endif

/// Threshold between ORB3D (large) and greedy (small) as fraction of
/// active particles
#define LARGE_PHASE_THRESHOLD 0.0001

/// @brief Implement load balancing: store loads and decide between
/// ORB3D and greedy.
void MultistepLB_notopo::work(BaseLB::LDStats* stats)
{
#if CMK_LBDB_ON
  // find active objects - mark the inactive ones as non-migratable
  int count;
  CkPrintf("[%d]  ^^^^^^^^^^^^^ CkMyNodeSize %d ^^^^^^^^^^^^\n", CkMyPe(), CkMyNodeSize());
  
  stats->makeCommHash();
  for(int i = 0; i < stats->n_objs; i++){
    LDObjHandle &handle = tpCentroids[i].handle;
    tpCentroids[i].tag = stats->getHash(handle.id, handle.omhandle.id);
  }
  int phase = determinePhase(tpCentroids[0].activeRung);
  int prevPhase = tpCentroids[0].prevActiveRung;
  float *ratios = new float[stats->n_objs];
  // save pointers to centroids of treepieces

  int numActiveObjects = 0;
  int numInactiveObjects = 0;

  // to calculate ratio of active particles in phase
  int numActiveParticles = 0;
  int totalNumParticles = 0;

  // This won't work if stats->count changes in every step like in hierarchical
  // probably?
  //if (existing_tps_on_pe == NULL) {
  //  existing_tps_on_pe = new int[stats->count];
  //}
  //memcpy(existing_tps_on_pe, 0, stats->count * sizeof(int));

  
  for(int i = 0; i < stats->n_objs; i++){
    stats->to_proc[i] = stats->from_proc[i];
  }
  // update phase data 
  if (_lb_args.debug()>=2) {
    CkPrintf("merging previous phase %d data; current phase: %d\n", prevPhase, phase);
  }

  saveStatsData(prevPhase, stats); 
//  if (step() == 0) {
//    return;
//  }

  if (prevNumActiveParts.length() == 0) {
    prevNumActiveParts.resize(stats->n_objs);
    objDataPrevPred.resize(stats->n_objs);
    objDataPrevPredPart.resize(stats->n_objs);
  }
  peddloads = new double[stats->count];
  for (int i = 0; i < stats->count; i++) {
    peddloads[i] = 0.0;
  }
 
  unsigned int maxMovedWithin = 0;
  unsigned int maxMoved = 0;
  int maxMovedPe, maxMovedWithinPe;
  unsigned int totalMoved = 0;
  unsigned int totalMovedWithin = 0;
  for(int i = 0; i < stats->n_objs; i++){
    int tp = tpCentroids[i].tp;
    int lb = tpCentroids[i].tag;
    peddloads[stats->from_proc[lb]] += tpCentroids[i].ddtime;
    
    if(tpCentroids[i].myNumParticles != 0){
      ratios[tp] = tpCentroids[i].numActiveParticles/(float)tpCentroids[i].myNumParticles;
    }
    else{
      ratios[tp] = 1.0;
    }
    numActiveParticles += tpCentroids[i].numActiveParticles;
    totalNumParticles += tpCentroids[i].myNumParticles;
    prevNumActiveParts[tp] = tpCentroids[i].numActiveParticles;

    if (tp == 239226) {
      CkPrintf("In LB TP %d has numActiveParticles %d\n", tp,
      tpCentroids[i].numActiveParticles);
    }
    if (maxMovedWithin < tpCentroids[i].numMovedParticlesWithin) {
      maxMovedWithin = tpCentroids[i].numMovedParticlesWithin;
      maxMovedWithinPe = stats->from_proc[lb];
    }
    if (maxMoved< tpCentroids[i].numMovedParticles) {
      maxMoved= tpCentroids[i].numMovedParticles;
      maxMovedPe = stats->from_proc[lb];
    }
    totalMoved += tpCentroids[i].numMovedParticles;
    totalMovedWithin += tpCentroids[i].numMovedParticlesWithin;

    if(tpCentroids[i].numActiveParticles == 0){
      numInactiveObjects++;

     // existing_tps_on_pe[stats->from_proc[i]] =
     //   existing_tps_on_pe[stats->from_proc[i]] + 1;

      if(stats->objData[lb].migratable){
        stats->objData[lb].migratable = 0;
     #ifdef MCLBMSV
        CkPrintf("marking object %d non-migratable (inactive)\n", tpCentroids[i].tag);
     #endif
        stats->n_migrateobjs--;
      }
    }
    else{
      numActiveObjects++;
    }
  }
  CkPrintf("maxParticlesMovedWithin %u (PE %d) maxParticlesMoved %u (PE %d) totalMovedWithin %u totalMoved %u\n", maxMovedWithin, maxMovedWithinPe,
  maxMoved, maxMovedPe, totalMovedWithin, totalMoved);

  // Do a prefix sum
//  for (int i = 1; i < stats->count; i++) {
//    existing_tps_on_pe[i] = existing_tps_on_pe[i] + existing_tps_on_pe[i-1];
//  }

  CkPrintf("numActiveObjects: %d, numInactiveObjects: %d\n", numActiveObjects,
	   numInactiveObjects);
/*  if(numInactiveObjects < 1.0*numActiveObjects) {
	// insignificant number of inactive objects; migrate them anyway
  	for(int i = 0; i < stats->n_objs; i++){
    	    int lb = tpCentroids[i].tag;
            if(!stats->objData[lb].migratable
		&& tpCentroids[i].myNumParticles > 0){
        	stats->objData[lb].migratable = 1;
        	stats->n_migrateobjs++;
		numActiveObjects++;
		numInactiveObjects--;
		}
	    }
  	CkPrintf("Migrating all: numActiveObjects: %d, numInactiveObjects: %d\n", numActiveObjects, numInactiveObjects);
	}
*/


  // TODO: Remove
  for(int i = 0; i < stats->n_objs; i++){
    int tp = tpCentroids[i].tp;
    int lb = tpCentroids[i].tag;

    objDataPrevPred[tp] = stats->objData[lb].wallTime;
    objDataPrevPredPart[tp] = tpCentroids[i].numActiveParticles;
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

  delete []ratios;

  // let the strategy take over on this modified instrumented data and processor information
  if((float)numActiveParticles/totalNumParticles > LARGE_PHASE_THRESHOLD){
    //if (_lb_args.debug()>=2) {
      CkPrintf("******** FULL WORK STEP *********!\n");
    //}
    work2(stats,count,phase,prevPhase);
  }     // end if phase == 0
  else{
      CkPrintf("******** NO LB AT ALL STEP *********!\n");
    // greedy(stats,count,phase,prevPhase);
  }
  //delete[] existing_tps_on_pe;
#endif //CMK_LDB_ON

}

//**************************************
// ORB3DLB functions
//**************************************
//
void MultistepLB_notopo::greedy(BaseLB::LDStats *stats, int count, int phase, int prevPhase){

  int numobjs = stats->n_objs;
  int nmig = stats->n_migrateobjs;
  CkPrintf("[GREEDY] objects total %d active %d\n", numobjs,nmig);

  TPObject *tp_array = new TPObject[nmig];
  int j = 0;
  for(int i = 0; i < stats->n_objs; i++){
    int lb = tpCentroids[i].tag;

    if(!stats->objData[lb].migratable) continue;
    tp_array[j].migratable = stats->objData[lb].migratable;
    if(step() == 0){
      tp_array[j].load = tpCentroids[i].myNumParticles; 
    }
    else{
      tp_array[j].load = stats->objData[lb].wallTime;
    }
    tp_array[j].lbindex = lb;
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
  CkPrintf("GREEDY MEASURED CPU LOAD prev %d\n", prevPhase);
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
void MultistepLB_notopo::work2(BaseLB::LDStats *stats, int count, int phase, int prevPhase){

  int pewithmaxddload = -1;
  unsigned int maxmig;
  double maxddload = 0.0;


  for (int i = 0; i < stats->count; i++) {
    if (peddloads[i] > maxddload) {
      maxddload = peddloads[i];
      pewithmaxddload = i;
    }
  }

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

  int maxldtp, maxldidx, maxldpart;
  float maxldtpld = 0.0;
  int fxdtpid = -1;

  for(int i = 0; i < numobjs; i++){
    int lb = tpCentroids[i].tag;


    if (tpCentroids[i].tp == 142747) {
      fxdtpid = lb;
    }

    if (stats->from_proc[lb] == pewithmaxddload) {
      CkPrintf("** [%d] DD load of %f for tp %d migs %u mig within %u total parts %u\n", pewithmaxddload,
        tpCentroids[i].ddtime, tpCentroids[i].tp,
        tpCentroids[i].numMovedParticles,
        tpCentroids[i].numMovedParticlesWithin, tpCentroids[i].myNumParticles);
    }
    if(!stats->objData[lb].migratable) continue;
 
    float load;
    if(step() == 0){
      load = tpCentroids[i].myNumParticles;
    }
    else{
      load = stats->objData[lb].wallTime;
    }

    // CkPrintf("Before calling Orb %d %f \n",lb, load);
    if (maxldtpld < load) {
      maxldtpld = load;
      maxldtp = tpCentroids[i].tp;
      maxldidx = lb;
      maxldpart = tpCentroids[i].numActiveParticles;
    }
    //if (tpCentroids[i].numActiveParticles > 1.3*7560) {
    //  CkPrintf("[%d] overloaded TP on PE %d with ld %f and activeParticles %u\n",
    //  tpCentroids[i].tp, stats->from_proc[lb], stats->objData[lb].wallTime,
    //  tpCentroids[i].numActiveParticles);
    //}
    

    tpEvents[XDIM].push_back(Event(tpCentroids[i].vec.x,load,numProcessed));
    tpEvents[YDIM].push_back(Event(tpCentroids[i].vec.y,load,numProcessed));
    tpEvents[ZDIM].push_back(Event(tpCentroids[i].vec.z,load,numProcessed));

    tp_array[numProcessed]= OrbObject(lb,tpCentroids[i].myNumParticles);
    tp_array[numProcessed].centroid = tpCentroids[i].vec;

   
    tp_array[numProcessed].lbindex = lb;
    numProcessed++;

  }
  CkPrintf("** DD max load %f for pe %d numobjs %d\n", maxddload,
  pewithmaxddload, numobjs);
  delete[] peddloads;
  CkAssert(numProcessed==nmig);
  
  orbPrepare(tpEvents, box, nmig, stats);
  orbPartition(tpEvents,box,CkNumNodes(),tp_array, stats);
  CkPrintf("MultistepLB> Done OrbPartition\n");
  //orbPartition(tpEvents,box,stats->count,tp_array, stats);
  balanceTPs(stats);

  CkPrintf("MultistepLB> stats: maxldtptag: %d maxldidx: %d ld: %f particles: %d PE: %d\n", maxldtp,
  maxldidx, maxldtpld, maxldpart, stats->to_proc[maxldidx]);
  CkPrintf("MultistepLB> PE %d from_proc %d to_proc %d\n", 142747,
    stats->from_proc[fxdtpid], stats->to_proc[fxdtpid]);

  refine(stats, numobjs, &prevMaxPredPe);

  for(int i = 0; i < numobjs; i++){
    int lb = tpCentroids[i].tag;
    int tp = tpCentroids[i].tp;

    if (stats->to_proc[lb] == prevMaxPredPe) {
      CkPrintf("\t In this LB TP %d is in %d PE with LD %f and cur active particles %d\n",
          tpCentroids[i].tp, prevMaxPredPe, stats->objData[lb].wallTime,
          tpCentroids[i].numActiveParticles);
    }
  }


  if(_lb_args.debug() >= 2) {
      // Write out "particle file" of load balance information
      char achFileName[1024];
      sprintf(achFileName, "lb.%d.sim", step());
      FILE *fp = fopen(achFileName, "w");
      CkAssert(fp != NULL);
      fprintf(fp, "%d %d 0\n", nrecvd, nrecvd);
      for(int i = 0; i < nrecvd; i++) {
	  CkAssert(tpCentroids[i].tag < numobjs);
	  CkAssert(tpCentroids[i].tag >= 0);
	  fprintf(fp, "%g %g %g %g 0.0 0.0 0.0 %d %d\n",
		  stats->objData[tpCentroids[i].tag].wallTime,
		  tpCentroids[i].vec.x,
		  tpCentroids[i].vec.y,
		  tpCentroids[i].vec.z,
		  stats->to_proc[tpCentroids[i].tag],
		  tpCentroids[i].tp);
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

void MultistepLB_notopo::balanceTPs(BaseLB::LDStats* stats) {
  
  double* counts = new double[stats->count];
  memset(counts, 0.0, stats->count * sizeof(double));
  double totalld = 0.0;
  vector<vector<int> > objpemap;
  objpemap.resize(stats->count);
  
  for (int i = 0; i < stats->n_objs; i++) {
    counts[stats->to_proc[i]] += stats->objData[i].wallTime;
    totalld += stats->objData[i].wallTime;
    objpemap[stats->to_proc[i]].push_back(i);
  }
  double avgldperpe = totalld / stats->count;
  vector<int> unldpes;
  vector<int> ovldpes;
  double th = 1.2;
  
  int maxcountoftps = 0;
  int pewithmax = -1;
  for (int i = (stats->count-1); i >= 0; i--) {
    if (counts[i] > th*avgldperpe) {
      ovldpes.push_back(i);
    } else if (counts[i] < (0.7*avgldperpe)) {
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
      int n_proc = -1;
      if (to_proc != ovlpe || !stats->objData[i].migratable) {
        continue;
      }
      if (stats->objData[i].wallTime < 0.6) {
        continue;
      }

      n_proc = unldpes[undcount];
      if ((counts[n_proc] + stats->objData[i].wallTime) >=
        (counts[to_proc]-stats->objData[i].wallTime)) {
        if(step() == 8) {
          CkPrintf("Could not transfer to undloaded pe ld %f obj %f\n",
          counts[n_proc], stats->objData[i].wallTime);
        }
        continue;
      }
      stats->to_proc[i] = n_proc;
      counts[to_proc] = counts[to_proc] - stats->objData[i].wallTime;
      counts[n_proc] = counts[n_proc] + stats->objData[i].wallTime;
      if (counts[n_proc] > 0.8*avgldperpe) {
        undcount++;
      }
      if (counts[to_proc] > th*avgldperpe) {
        ovldpes.push_back(ovlpe);
        push_heap(ovldpes.begin(), ovldpes.end(), PeLdGreater(counts));
      }
      succ = true;
      if(step() == 8) {
        CkPrintf("transfered to undloaded pe ld %f obj %f new ld for PE %d %f\n",
            counts[n_proc], stats->objData[i].wallTime, to_proc, counts[to_proc]);
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
    if (!succ && step()==8) {

      CkPrintf("Couldn't find any obj to migrate from %d pe with ld %f\n",
      ovlpe, counts[ovlpe]);
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
/*
void MultistepLB_notopo::balanceTPs(BaseLB::LDStats* stats) {
  int* counts = new int[stats->count];
  memset(counts, 0, stats->count * sizeof(int));
  for (int i = 0; i < stats->n_objs; i++) {
    counts[stats->to_proc[i]] = counts[stats->to_proc[i]] + 1;
  }
  int avgtpsperpe = stats->n_objs / stats->count;
  vector<int> unldpes;
  vector<int> ovldpes;
  
  int threshold = 3*avgtpsperpe;
  int threshold1 = 5*avgtpsperpe;
  int threshold2 = 7*avgtpsperpe;
  int maxcountoftps = 0;
  int pewithmax = -1;
  for (int i = (stats->count-1); i >= 0; i--) {
    if (counts[i] > 50) {
      ovldpes.push_back(i);
    } else if (counts[i] < (threshold1)) {
      unldpes.push_back(i);
    }
    if (counts[i] > maxcountoftps) {
      maxcountoftps = counts[i];
      pewithmax = i;
    }
  }
  CkPrintf("[%d] Maxpewithtps %d thresholds %d and %d\n", pewithmax, maxcountoftps, threshold1, threshold2);
  int undcount = 0;

  int* tmpcounts = new int[stats->count];
  memcpy(tmpcounts, counts, stats->count * sizeof(int));

  for (int i = stats->n_objs-1; i >=0 ; i--) {
    if (undcount > unldpes.size()) {
      break;
    }
    int to_proc = stats->to_proc[i];
    int n_proc = -1;
    if (counts[to_proc] > 70) {
      CkPrintf("TP %d is on big proc %d ismigr? %d ld %f\n", i, to_proc,
      stats->objData[i].migratable, stats->objData[i].wallTime);
    }
    if (counts[to_proc] > threshold1 && stats->objData[i].wallTime < 0.02) {
      n_proc = unldpes[undcount];
      stats->to_proc[i] = n_proc;
      counts[to_proc] = counts[to_proc] - 1;
      counts[n_proc] = counts[n_proc] + 1;
      if (counts[n_proc] > threshold1) {
        undcount++;
      }
    }
  }

  maxcountoftps = 0;
  pewithmax = -1;
  for (int i = (stats->count-1); i >= 0; i--) {
    if (counts[i] > maxcountoftps) {
      maxcountoftps = counts[i];
      pewithmax = i;
    }
  }
  CkPrintf("[%d] Afterwards Maxpewithtps %d previously %d\n", pewithmax, maxcountoftps, tmpcounts[pewithmax]);
  delete[] counts;
  delete[] tmpcounts;
}
*/
void MultistepLB_notopo::receiveAvgLoad(double avg_load) {
  avg_load_after_lb = avg_load;
}

void MultistepLB_notopo::getLoadInfo(double& avg_load, double& my_load) {
  avg_load = avg_load_after_lb;
  my_load = my_load_after_lb;
}

void MultistepLB_notopo::clearPeLoad() {
  my_load_after_lb = 0.0;
  tpscount = 0;
  tpsregisteredfordd = 0;
  tpsregisteredforacc = 0;
  tpsonpe.clear();
  tpsonpeforacc.clear();
}

void MultistepLB_notopo::addToPeLoad(double tpload) {
  //__sync_add_and_fetch(&my_load_after_lb, tpload);
  my_load_after_lb += tpload;
}

void MultistepLB_notopo::addTpCount() {
  tpscount++;
}

void tpPar(int start, int end, void *result, int pnum, void * param) {
  std::vector<TreePiece*> *tps = (std::vector<TreePiece*> *) param;

  for (int i = start; i <= end; i++) {
    (*tps)[i]->unshuffleParticlesWoDDCb();
  }
}

void tpParForAcc(int start, int end, void *result, int pnum, void * param) {
  std::vector<TreePiece*> *tps = (std::vector<TreePiece*> *) param;

  for (int i = start; i <= end; i++) {
    (*tps)[i]->shuffleAfterQDSpecificOpt();
    //if ((*tps)[i]->myPeIs() == 24645) {
   // if ((*tps)[i]->myPeIs() != CkMyPe()) {
   //   CkPrintf("[%d] I am from %d but running on %d\n", (*tps)[i]->thisIndex,
   //   (*tps)[i]->myPeIs(), CkMyPe());
   // }
  }
    //if (CkMyPe() == 0 || (CkMyPe() >= 24645 && CkMyPe() < 24676)) {
    //  CkPrintf("[%d] ^^^^ Got %d TPs ^^^^\n", CkMyPe(), (end-start+1));
    //}
}

void MultistepLB_notopo::addTpForDD(TreePiece *tp) {
//  if (tpsregisteredfordd == 0) {
//    timestart = CkWallTimer();
//  }
  tpsregisteredfordd++;
  int avgtpspe = 262144/CkNumPes();
  if (tpsregisteredfordd <= avgtpspe) {
    tp->unshuffleParticlesWoDDCb();
    return;
  }
  tpsonpe.push_back(tp);
  if (tpsregisteredfordd == tpscount) {
    // CkLoopstuff
   // if (CkMyPe() == 0) {
   //   CkPrintf("**Time to do CkLoop count %d  ckloop done on %d on PE %d******\n", tpscount, tpsonpe.size(), CkMyPe());
   // }
//    double endtime = CkWallTimer();
    //if (CkMyPe()%1024 == 0) {
    //  CkPrintf("time taken in before starting loop %f\n", endtime - timestart);
    //}
    int numchunks = (tpsonpe.size() >= 64) ? 64 : tpsonpe.size();
    CkLoop_Parallelize(tpPar, 1, &tpsonpe, numchunks, 0, tpsonpe.size()-1, 1, NULL, CKLOOP_NONE);
    tpsregisteredfordd = 0;
  }
}

void MultistepLB_notopo::addTpForAcceptSorted(TreePiece *tp) {
  if (tpsregisteredforacc == 0) {
    timestart = CkWallTimer();

   // if (CkMyPe()%1111 == 0) {
   //   CkPrintf("[%d] has Memory %f MB tpsregisteredforacc %d\n", CkMyPe(),
   //       CmiMemoryUsage()/(1024.0*1024.0), tpsregisteredforacc);
   // }
  }
  int avgtpspe = 262144/CkNumPes();
  tpsregisteredforacc++;
  if (tpsregisteredforacc <= avgtpspe) {
    tp->shuffleAfterQDSpecificOpt();
    return;
  }
  tpsonpeforacc.push_back(tp);
  if (tpsregisteredforacc == tpscount) {
    // CkLoopstuff
    if (CkMyPe() == 0) {
      CkPrintf("**Time to do CkLoop count %d  ckloop done on %d on PE %d******\n", tpscount, tpsonpeforacc.size(), CkMyPe());
    }
    double endtime = CkWallTimer();
    int numchunks = (tpsonpeforacc.size() >= 64) ? 64 : tpsonpeforacc.size();
    CkLoop_Parallelize(tpParForAcc, 1, &tpsonpeforacc, numchunks, 0, tpsonpeforacc.size()-1, 1, NULL, CKLOOP_NONE);
   // for (int i = 0; i < tpsonpeforacc.size(); i++) {
   //   tpsonpeforacc[i]->shuffleAfterQDSpecificOpt();
   // }
    tpsregisteredforacc = 0;
  }
}

void MultistepLB_notopo::tpWorkLocalContinue(TpWorkMsg *msg) {
  //CkPrintf("[%d] Multistep calling work on tp %d nodesize %d\n", CkMyPe(),
  //(msg->tp)->getIndex());
  //(msg->tp)->calculateGravityRemoteForeign(msg->buckets, msg->chunkNum);
  //(msg->tp)->doForeignBuckets(msg);
  //CkPrintf("[%d] Multistep tp callign continueForeignBucketLocal ptr %d\n", CkMyPe(), msg);
  (msg->tp)->doForeignBucketsLocal(msg);
}

void MultistepLB_notopo::tpWorkEwaldContinue(TpWorkMsg *msg) {
  //CkPrintf("[%d] Multistep calling work on tp %d nodesize %d\n", CkMyPe(),
  //(msg->tp)->getIndex());
  //(msg->tp)->calculateGravityRemoteForeign(msg->buckets, msg->chunkNum);
  //(msg->tp)->doForeignBuckets(msg);
  //CkPrintf("[%d] Multistep tp callign continueForeignBucketEwald ptr %d\n", CkMyPe(), msg);
  (msg->tp)->doForeignBucketsEwald(msg);
}

void MultistepLB_notopo::tpWorkRemoteContinue(TpWorkMsg *msg) {
  //CkPrintf("[%d] Multistep calling work on tp %d nodesize %d\n", CkMyPe(),
  //(msg->tp)->getIndex());
  //(msg->tp)->calculateGravityRemoteForeign(msg->buckets, msg->chunkNum);
  //CkPrintf("[%d] Multistep tp callign continueForeignBucketRemote ptr %d\n", CkMyPe(), msg);
  (msg->tp)->doForeignBucketsRemote(msg);
}


void MultistepLB_notopo::tpWorkLocal(TpWorkMsg *msg) {
  //CkPrintf("[%d] Multistep calling work on tp %d nodesize %d\n", CkMyPe(),
  //(msg->tp)->getIndex());
  //(msg->tp)->calculateGravityRemoteForeign(msg->buckets, msg->chunkNum);
  //(msg->tp)->doForeignBuckets(msg);
  //CkPrintf("[%d] Multistep tp callign startForeignBucketLocal ptr %d\n", CkMyPe(), msg);
  (msg->tp)->startForeignBucketLocal(msg);
}

void MultistepLB_notopo::tpWorkEwald(TpWorkMsg *msg) {
  //CkPrintf("[%d] Multistep calling work on tp %d nodesize %d\n", CkMyPe(),
  //(msg->tp)->getIndex());
  //(msg->tp)->calculateGravityRemoteForeign(msg->buckets, msg->chunkNum);
  //(msg->tp)->doForeignBuckets(msg);
  //CkPrintf("[%d] Multistep tp callign startForeignBucketEwald ptr %d\n", CkMyPe(), msg);
  (msg->tp)->startForeignBucketEwald(msg);
}

void MultistepLB_notopo::tpWorkRemote(TpWorkMsg *msg) {
  //CkPrintf("[%d] Multistep calling work on tp %d nodesize %d\n", CkMyPe(),
  //(msg->tp)->getIndex());
  //(msg->tp)->calculateGravityRemoteForeign(msg->buckets, msg->chunkNum);
  //CkPrintf("[%d] Multistep tp callign startForeignBucketRemote ptr %d\n", CkMyPe(), msg);
  (msg->tp)->startForeignBucketRemote(msg);
}

void MultistepLB_notopo::pup(PUP::er &p){
  CentralLB::pup(p);
  if(p.isPacking() && haveTPCentroids) {
    // if checkpointing, no need to 
    // keep around the centroid message
    delete tpmsg;
    haveTPCentroids = false;
  }
  p | haveTPCentroids; 
  p | procsPerNode;
  p | savedPhaseStats;
}

void LightweightLDStats::pup(PUP::er &p){
  p|n_objs;
  p|n_migrateobjs;
  p|n_activeparticles;
  p|objData;
}

#include "MultistepLB_notopo.def.h"
