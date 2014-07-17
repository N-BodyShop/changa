#include <charm++.h>
#include "cklists.h"
#include "MultistepOrbLB.h"
#include "TopoManager.h"
#include "ParallelGravity.h"
#include "Vector3D.h"
#include <queue>

extern CProxy_TreePiece treeProxy;
using namespace std;
//#define ORB3DLB_NOTOPO_DEBUG CkPrintf

CreateLBFunc_Def(MultistepOrbLB, "Works best with multistepped runs; uses Orb3D_notopo for larger steps, greedy otherwise");


MultistepOrbLB::MultistepOrbLB(const CkLBOptions &opt): OrbLB(opt)
{
  lbname = "MultistepOrbLB";

  if (CkMyPe() == 0){
    CkPrintf("[%d] MultistepOrbLB created\n",CkMyPe());
  }
  
  haveTPCentroids = false;
}

/// @brief Get position centroids of all TreePieces
/// @param msg Reduction message with a concatenation of all centroids.
void MultistepOrbLB::receiveCentroids(CkReductionMsg *msg){
  if(haveTPCentroids){
    delete tpmsg;
  }
  tpCentroids = (TaggedVector3D *)msg->getData();
  nrecvd = msg->getSize()/sizeof(TaggedVector3D);
  tpmsg = msg;
  haveTPCentroids = true;
  treeProxy.doAtSync();
  CkPrintf("MultistepOrbLB: receiveCentroids %d elements, msg length: %d\n", nrecvd, msg->getLength()); 
}

bool MultistepOrbLB::QueryBalanceNow(int step){
 if(CkMyPe() == 0) CkPrintf("Orb3dOrbLB: Step %d\n", step);
  if(step == 0) return false;
  return true;

}

// helper functions for multistepping
#ifdef MCLBMS

// determine phase based on lastActiveRung as saved in map.tpCentroids
unsigned int MultistepOrbLB::determinePhase(unsigned int lastActiveRung){
  return lastActiveRung;
}

// merge data instrumented in previous iteration of computation with data from earlier iterations
// this data must be stored according to the tree piece 
// index of the treepiece (obtained from treepiece[.].tag)
// since the ldb index of an object changes from iteration
// to iteration.
void MultistepOrbLB::mergeInstrumentedData(int phase, BaseLB::LDStats *stats){

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

  /*
  CkPrintf("**********************************************\n");
  CkPrintf("Actual object loads phase %d\n", phase);
  CkPrintf("**********************************************\n");
  for(int i = 0; i < stats->n_objs; i++){
    int tp = tpCentroids[i].tp;
    int lb = tpCentroids[i].tag;
    CkPrintf("tp %d load %f\n",tp,stats->objData[lb].wallTime);
  }
  CkPrintf("**********************************************\n");
  CkPrintf("Done actual object loads phase %d\n", phase);
  CkPrintf("**********************************************\n");
  */
  
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
    for(i = 0; i < stats->n_objs; i++){
      int tp = tpCentroids[i].tp;
      int lb = tpCentroids[i].tag;
      savedPhaseStats[len-1].objData[tp] = stats->objData[lb];
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
      for(i = 0; i < stats->n_objs; i++){
        int tp = tpCentroids[i].tp;
        int lb = tpCentroids[i].tag;
        savedPhaseStats[phase].objData[tp] = stats->objData[lb];
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
      }
    }
  }
  savedPhaseStats[whichPos].n_objs = stats->n_objs;
  savedPhaseStats[whichPos].n_migrateobjs = stats->n_migrateobjs;
#ifdef MCLBMSV
  //printData(savedPhaseStats[whichPos], phase, NULL);
#endif
}

// whether we have instrumented data for this phase
bool MultistepOrbLB::havePhaseData(int phase){
  return (savedPhaseStats.length() > phase && savedPhaseStats[phase].n_objs > 0);
}

void MultistepOrbLB::printData(BaseLB::LDStats &stats, int phase, int *revObjMap){
  int i;
  
  CkPrintf("---- data (%d): %d objects ----\n", phase, stats.n_objs);
  for(i = 0; i < stats.n_objs; i++){
     CkPrintf("%d: %g\n", i, 
	       stats.objData[i].wallTime);
  }
  CkPrintf("---- end data (%d) ----\n", phase);
}

#endif

/// Threshold between ORB (large) and none (small) as fraction of
/// active particles
#define LARGE_PHASE_THRESHOLD 0.0001

/// @brief Implement load balancing: store loads and decide between
/// ORB and none.
void MultistepOrbLB::work(BaseLB::LDStats* stats)
{
#if CMK_LBDB_ON
  // find active objects - mark the inactive ones as non-migratable
  int count;
  
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
  int64_t numActiveParticles = 0;
  int64_t totalNumParticles = 0;
  
  for(int i = 0; i < stats->n_objs; i++){
    stats->to_proc[i] = stats->from_proc[i];
  }
  // update phase data 
  if (_lb_args.debug()>=2) {
    CkPrintf("merging previous phase %d data; current phase: %d\n", prevPhase, phase);
  }
  mergeInstrumentedData(prevPhase, stats); 
  
  for(int i = 0; i < stats->n_objs; i++){
    int tp = tpCentroids[i].tp;
    int lb = tpCentroids[i].tag;
    if(tpCentroids[i].myNumParticles != 0){
      ratios[tp] = tpCentroids[i].numActiveParticles/(float)tpCentroids[i].myNumParticles;
    }
    else{
      ratios[tp] = 1.0;
    }
    numActiveParticles += tpCentroids[i].numActiveParticles;
    totalNumParticles += tpCentroids[i].myNumParticles;

    if(tpCentroids[i].numActiveParticles == 0){
      numInactiveObjects++;
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
  CkPrintf("numActiveObjects: %d, numInactiveObjects: %d\n", numActiveObjects,
	   numInactiveObjects);
  if(numInactiveObjects < 1.0*numActiveObjects) {
	// insignificant number of inactive objects; migrate them anyway
  	for(int i = 0; i < stats->n_objs; i++){
    	    int lb = tpCentroids[i].tag;
            if(!stats->objData[lb].migratable){
        	stats->objData[lb].migratable = 1;
        	stats->n_migrateobjs++;
		numActiveObjects++;
		numInactiveObjects--;
		}
	    }
  	CkPrintf("Migrating all: numActiveObjects: %d, numInactiveObjects: %d\n", numActiveObjects, numInactiveObjects);
	}

  // get load information for this phase, if possible
  // after this, stats->objData[] is indexed by tree piece
  // index since we are copying data from savedPhaseStats,
  // which was written into using tree piece indices
  if(havePhaseData(phase)){
#ifdef MCLBMSV
    CkPrintf("phase %d data available\n", phase);
#endif
    CkPrintf("phase %d data available\n", phase);
    for(int i = 0; i < stats->n_objs; i++){
      int tp = tpCentroids[i].tp;
      int lb = tpCentroids[i].tag;
      stats->objData[lb].wallTime = savedPhaseStats[phase].objData[tp].wallTime;
    }
  }
  else if(havePhaseData(0)){
#ifdef MCLBMSV
    CkPrintf("phase %d data unavailable, using phase 0 loads\n", phase);
#endif
    CkPrintf("phase %d data unavailable, using phase 0 loads\n", phase);
    //CkPrintf("using phase 0 loads\n", phase);
    for(int i = 0; i < stats->n_objs; i++){
      int tp = tpCentroids[i].tp;
      int lb = tpCentroids[i].tag;
      stats->objData[lb].wallTime = ratios[tp]*savedPhaseStats[0].objData[tp].wallTime;
    }
  }
  else{
#ifdef MCLBMSV
    CkPrintf("phase %d data unavailable\n", phase);
#endif
    CkPrintf("phase %d data unavailable\n", phase);
    delete[] ratios;
    return;
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

  delete []ratios;

  // let the strategy take over on this modified instrumented data and processor information
  if((float)numActiveParticles/totalNumParticles > LARGE_PHASE_THRESHOLD){
    if (_lb_args.debug()>=2) {
      CkPrintf("******** BIG STEP *********!\n");
    }
    OrbLB::work(stats);
  }     // end if phase == 0
#endif //CMK_LDB_ON

}

void MultistepOrbLB::pup(PUP::er &p){
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

#include "MultistepOrbLB.def.h"
