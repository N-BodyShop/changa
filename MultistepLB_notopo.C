#include <charm++.h>
#include "cklists.h"
#include "MultistepLB_notopo.h"
#include "TopoManager.h"
#include "ParallelGravity.h"
#include "Vector3D.h"
#include <queue>

extern CProxy_TreePiece treeProxy;
using namespace std;
#define ORB3DLB_NOTOPO_DEBUG

CreateLBFunc_Def(MultistepLB_notopo, "Works best with multistepped runs; uses Orb3D_notopo for larger steps, greedy otherwise");


MultistepLB_notopo::MultistepLB_notopo(const CkLBOptions &opt): CentralLB(opt)
{
  lbname = "MultistepLB_notopo";

  if (CkMyPe() == 0){
    CkPrintf("[%d] MultistepLB_notopo created\n",CkMyPe());
  }

  
  haveTPCentroids = false;

}

void MultistepLB_notopo::receiveCentroids(CkReductionMsg *msg){

  if(haveTPCentroids){
    haveTPCentroids = false;
    delete tpmsg; 
  }

  tpCentroids = (TaggedVector3D *)msg->getData();
  nrecvd = msg->getGcount();
  if (_lb_args.debug()>=2) {
    CkPrintf("MultistepLB_notopo: receiveCentroids started: %d elements, msg length: %d\n", msg->getGcount(), msg->getLength()); 
  }
  haveTPCentroids = true;
  tpmsg = msg;
  treeProxy.doAtSync();
  if (_lb_args.debug()>=2) {
    CkPrintf("MultistepLB_notopo: receiveCentroids done\n");  
  }
}

CmiBool MultistepLB_notopo::QueryBalanceNow(int step){
 if(CkMyPe() == 0) CkPrintf("Orb3dLB_notopo: Step %d\n", step);
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

#define LARGE_PHASE_THRESHOLD 0.0

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

    if(false && tpCentroids[i].numActiveParticles == 0){
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
#ifdef MCLBMSV
  CkPrintf("numActiveObjects: %d, numInactiveObjects: %d\n", numActiveObjects, numInactiveObjects);
#endif



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
    greedy(stats,count,phase,prevPhase);
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

void MultistepLB_notopo::work2(BaseLB::LDStats *stats, int count, int phase, int prevPhase){
  int numobjs = stats->n_objs;
  int nmig = stats->n_migrateobjs;

  // this data structure is used by the orb3d strategy
  // to balance objects. it is NOT indexed by tree piece index
  // there are as many entries in it as there are
  // migratable (active) tree pieces
  OrbObject *tp_array = new OrbObject[nmig];

  if (_lb_args.debug()>=2) {
    CkPrintf("[work2] ready tp_array data structure\n");
  }

 CkVec<Event> tpEvents[NDIMS];
  for(int i = 0; i < NDIMS; i++){
    tpEvents[i].reserve(nmig);
  }

  OrientedBox<float> box;

  int numProcessed = 0;

  int j = 0;
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
 //   int lb = tpCentroids[i].tag;

    // CkPrintf("Before calling Orb %d %f \n",lb, load);

    tpEvents[XDIM].push_back(Event(tpCentroids[i].vec.x,load,j));
    tpEvents[YDIM].push_back(Event(tpCentroids[i].vec.y,load,j));
    tpEvents[ZDIM].push_back(Event(tpCentroids[i].vec.z,load,j));

    tp_array[j]= OrbObject(lb);
    tp_array[j].centroid = tpCentroids[i].vec;

   
    tp_array[j].lbindex = lb;
    j++;
    numProcessed++;

  }
  CkAssert(j==nmig);
  CkAssert(numProcessed == nmig);
  CkAssert(tpEvents[XDIM].length() == nmig);
  CkAssert(tpEvents[YDIM].length() == nmig);
  CkAssert(tpEvents[ZDIM].length() == nmig);


  mapping = &stats->to_proc;
  int dim = 0;

  CkPrintf("[Orb3dLB_notopo] sorting\n");
  for(int i = 0; i < NDIMS; i++){
    tpEvents[i].quickSort();
/*  for(int j = 0; j < tpEvents[i].length(); j++){
	int own = tpEvents[i][j].owner;
	CkPrintf("After sorting %d %d \n", own, tp_array[own].lbindex);

  }
  */ 
  }
  box.lesser_corner.x = tpEvents[XDIM][0].position;
  box.lesser_corner.y = tpEvents[YDIM][0].position;
  box.lesser_corner.z = tpEvents[ZDIM][0].position;

  box.greater_corner.x = tpEvents[XDIM][nmig-1].position;
  box.greater_corner.y = tpEvents[YDIM][nmig-1].position;
  box.greater_corner.z = tpEvents[ZDIM][nmig-1].position;

  nextProc = 0;

  procload.resize(stats->count);
  procbox.resize(stats->count);
  for(int i = 0; i < stats->count; i++){
    procload[i] = 0.0;
  }
 


  
  orbPartition(tpEvents,box,stats->count,tp_array);
 /* 
  float *objload = new float[stats->count];
  for(int i = 0; i < stats->count; i++){
    objload[i] = 0.0;
  }
  for(j = 0; j < nmig; j++){
    float load = tp_array[j].load;
    int lb = tp_array[j].lbindex;
    int pe = stats->to_proc[lb];

    objload[pe] += load; 
  }
*/
  /*
  CkPrintf("******************************\n");
  CkPrintf("CPU LOAD PREDICTIONS phase %d\n", phase);
  CkPrintf("******************************\n");
  for(int i = 0; i < stats->count; i++){
    CkPrintf("[pestats] %d %g \n", 
                               i,
                               objload[i]);
  }
  */


  CkPrintf("******************************\n");
  CkPrintf("MEASURED CPU LOAD prev %d\n", prevPhase);
  CkPrintf("******************************\n");
  for(int i = 0; i < stats->count; i++){
    CkPrintf("[pestats] %d %g %g\n", 
                               i,
                               stats->procs[i].total_walltime,
                               stats->procs[i].idletime
                               ); 
  }


 // delete[] objload;
  delete[] tp_array;
}

void MultistepLB_notopo::orbPartition(CkVec<Event> *events, OrientedBox<float> &box, int nprocs, OrbObject * tp){

  ORB3DLB_NOTOPO_DEBUG("partition events %d %d %d nprocs %d\n", 
            events[XDIM].length(),
            events[YDIM].length(),
            events[ZDIM].length(),
            nprocs
            );
  int numEvents = events[XDIM].length();
  CkAssert(numEvents == events[YDIM].length());
  CkAssert(numEvents == events[ZDIM].length());

  if(nprocs <= 1){
    ORB3DLB_NOTOPO_DEBUG("base: assign %d tps to proc %d\n", numEvents, nextProc);
    // direct assignment of tree pieces to processors
    if(numEvents > 0) CkAssert(nprocs != 0);
    float totalLoad = 0.0;
    for(int i = 0; i < events[XDIM].length(); i++){
      Event &ev = events[XDIM][i];
      OrbObject &orb = tp[ev.owner];
      (*mapping)[orb.lbindex] = nextProc;
      totalLoad += ev.load;
      procbox[nextProc].grow(orb.centroid);
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
  if(splitIndex == numEvents) {
      ORB3DLB_NOTOPO_DEBUG("evenly split 0 load\n");
      splitIndex = splitIndex/2;
      }
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
    CkAssert(tp[ev.owner].partition == INVALID_PARTITION);
    tp[ev.owner].partition = LEFT_PARTITION;
  }
  for(int i = splitIndex; i < numEvents; i++){
    Event &ev = events[longestDim][i];
    CkAssert(ev.owner >= 0);
    CkAssert(tp[ev.owner].partition == INVALID_PARTITION);
    tp[ev.owner].partition = RIGHT_PARTITION;
  }

  CkVec<Event> leftEvents[NDIMS];
  CkVec<Event> rightEvents[NDIMS];

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
  memcpy(leftEvents[longestDim].getVec(),events[longestDim].getVec(),sizeof(Event)*nleft);
  memcpy(rightEvents[longestDim].getVec(),events[longestDim].getVec()+splitIndex,sizeof(Event)*nright);
  
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
  CkVec<Event> &eraseVec = events[longestDim];
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
    events[i].free();
  }

  orbPartition(leftEvents,leftBox,nlprocs,tp);
  orbPartition(rightEvents,rightBox,nrprocs,tp);
}

#define ZERO_THRESHOLD 0.00001

#if 0
void Orb3dLB_notopo::directMap(int tpstart, int tpend, int nodestart, int nodeend){
  //CkPrintf("[Orb3dLB_notopo] mapping %d objects to Node (%d,%d,%d)\n", ntp, nodes[0].x, nodes[0].y, nodes[0].z);

  std::priority_queue<TPObject> pq_obj;
  std::priority_queue<Processor> pq_proc;

  float load = 0.0;
  CkAssert(nodestart==(nodeend-1));
  for(int i = tpstart; i < tpend; i++){
    //CkPrintf("obj %d thisindex %d %d %f %f %f %f to node %d %d %d\n", tp[i].lbindex, tp[i].index, tp[i].nparticles, tp[i].load, tp[i].centroid.x, tp[i].centroid.y, tp[i].centroid.z, nodes[0].x, nodes[0].y, nodes[0].z);
    load += tps[i].load;
    pq_obj.push(tps[i]);
  }
  //CkPrintf("node %d %d %d total load %f\n", nodes[0].x, nodes[0].y, nodes[0].z, load);
  
  for(int i = 0; i < procsPerNode; i++){
    Processor p;
    p.load = 0.0;
    p.t = i;
    pq_proc.push(p);
  }

  int currentZeroProc = 0;
  while(!pq_obj.empty()){
    TPObject tp = pq_obj.top();
    pq_obj.pop();

    // spread around zero-load objects
    // disabled to reduce the number of migrations, and
    // check whether this might solve the BG/P crash
    if(tp.load < ZERO_THRESHOLD){
      /*
      (*mapping)[tp.lbindex] = nodes[0].procRanks[currentZeroProc];
      currentZeroProc = currentZeroProc+1;
      if(currentZeroProc == procsPerNode){
        currentZeroProc = 0;
      }
      */
    }
    else{
      // if object has some non-zero load, assign it to a proc greedily
      Processor p = pq_proc.top();
      pq_proc.pop();

      //CkPrintf("proc %d load %f gets obj %d load %f\n", p.t, p.load, tp.lbindex, tp.load);

      p.load += tp.load;
      (*mapping)[tp.lbindex] = nodes[nodestart].procRanks[p.t];
#ifdef PRINT_BOUNDING_BOXES
      nodes[nodestart].box.grow(tp.centroid);
#endif

      pq_proc.push(p);
    }
  }
}
#endif

#define LOAD_EQUAL_TOLERANCE 1.02

int MultistepLB_notopo::partitionRatioLoad(CkVec<Event> &events, float ratio){
  float totalLoad = 0.0;
  for(int i = 0; i < events.length(); i++){
    totalLoad += events[i].load;
  }
  float lload = 0.0;
  float rload = totalLoad;
  float prevDiff = lload-ratio*rload;
  if(prevDiff < 0.0){
    prevDiff = -prevDiff;
  }

  int consider;
  for(consider = 0; consider < events.length();){
    float newll = lload + events[consider].load;
    float newrl = rload - events[consider].load;
    
    float newdiff = newll-ratio*newrl;
    if(newdiff < 0.0){
      newdiff = -newdiff;
    }

    //CkPrintf("consider load %f newdiff %f prevdiff %f\n", tp[consider].load, newdiff, prevDiff);

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

  //CkPrintf("partitionEvenLoad lload %f rload %f\n", lload, rload);

  return consider;
}

void MultistepLB_notopo::pup(PUP::er &p){
  CentralLB::pup(p);
  if(p.isPacking()){
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

/*@}*/
