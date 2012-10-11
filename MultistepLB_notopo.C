#include <charm++.h>
#include "cklists.h"
#include "MultistepLB_notopo.h"
#include "TopoManager.h"
#include "ParallelGravity.h"
#include "Vector3D.h"
#include <queue>

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



CmiBool MultistepLB_notopo::QueryBalanceNow(int step){
 if(CkMyPe() == 0) CkPrintf("Orb3dLB_notopo: Step %d\n", step);
  if(step == 0) return false;
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

  // tune alpha as needed - this is the merge parameter
  double alpha = 0.0;
  double savedWall;
  
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

  /* currently, only the first procsNeeded procs are used - could do something more sophisticated here in the future - FIXME */
#ifdef MCLBMSV
  CkPrintf("Processors 0 to %d active\n", procsNeeded-1);
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
  makeActiveProcessorList(stats, numActiveObjects);
  count = stats->count;

  delete []ratios;

  // let the strategy take over on this modified instrumented data and processor information
  if((float)numActiveParticles/totalNumParticles > LARGE_PHASE_THRESHOLD){
    if (_lb_args.debug()>=2) {
      CkPrintf("******** BIG STEP *********!\n");
    }
    work2(stats,count,phase,prevPhase);
  }     // end if phase == 0
  else{
    // greedy(stats,count,phase,prevPhase);
  }
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
    int tp = tpCentroids[i].tag;
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
  int numobjs = stats->n_objs;
  int nmig = stats->n_migrateobjs;

  // this data structure is used by the orb3d strategy
  // to balance objects. it is NOT indexed by tree piece index
  // there are as many entries in it as there are
  // migratable (active) tree pieces
 vector<OrbObject> tp_array;
 tp_array.reserve(nmig);

  if (_lb_args.debug()>=2) {
    CkPrintf("[work2] ready tp_array data structure\n");
  }

 vector<Event> tpEvents[NDIMS];
  for(int i = 0; i < NDIMS; i++){
    tpEvents[i].reserve(nmig);
  }

  OrientedBox<float> box;

  int numProcessed = 0;

  for(int i = 0; i < numobjs; i++){
    int tp = tpCentroids[i].tp;
    int lb = tpCentroids[i].tag;


    if(!stats->objData[lb].migratable) continue;
 
    float load;
    if(step() == 0){
      load = tpCentroids[i].myNumParticles;
    }
    else{
      load = stats->objData[lb].wallTime;
    }

    // CkPrintf("Before calling Orb %d %f \n",lb, load);

    tpEvents[XDIM].push_back(Event(tpCentroids[i].vec.x,load,numProcessed));
    tpEvents[YDIM].push_back(Event(tpCentroids[i].vec.y,load,numProcessed));
    tpEvents[ZDIM].push_back(Event(tpCentroids[i].vec.z,load,numProcessed));

    tp_array[numProcessed]= OrbObject(lb,tpCentroids[i].myNumParticles);
    tp_array[numProcessed].centroid = tpCentroids[i].vec;

   
    tp_array[numProcessed].lbindex = lb;
    numProcessed++;
  }
  CkAssert(numProcessed==nmig);

  
  orbPrepare(tpEvents, box, nmig, stats);
  orbPartition(tpEvents,box,stats->count,tp_array, stats);

  refine(stats, numobjs);
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
  p|objData;
}

#include "MultistepLB_notopo.def.h"
