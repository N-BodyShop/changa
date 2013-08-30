#include <charm++.h>
#include "cklists.h"
#include "MultistepLB.h"
#include "TopoManager.h"
#include "ParallelGravity.h"
#include "Vector3D.h"
#include <queue>

extern CProxy_TreePiece treeProxy;
using namespace std;

CreateLBFunc_Def(MultistepLB, "Works best with multistepped runs; uses Orb3D for larger steps, greedy otherwise");


MultistepLB::MultistepLB(const CkLBOptions &opt): CBase_MultistepLB(opt)
{
  lbname = "MultistepLB";

  if (CkMyPe() == 0){
    CkPrintf("[%d] MultistepLB created\n",CkMyPe());
  }

  TopoManager tmgr;

  int ppn = tmgr.getDimNT();

  int nx = tmgr.getDimNX();
  int ny = tmgr.getDimNY();
  int nz = tmgr.getDimNZ();
  int numnodes = nx*ny*nz; 

  if (CkMyPe() == 0){
    CkPrintf("[%d] Multistep Topo %d %d %d %d %d \n",CkMyPe(), nx, ny, nz, numnodes, ppn);
  }

  compares[0] = comparx;
  compares[1] = compary;
  compares[2] = comparz;

  pc[0] = pcx;
  pc[1] = pcy;
  pc[2] = pcz;

  haveTPCentroids = false;

}

void MultistepLB::receiveCentroids(CkReductionMsg *msg){

  if(haveTPCentroids){
    haveTPCentroids = false;
    delete tpmsg; 
  }

  tpCentroids = (TaggedVector3D *)msg->getData();
  nrecvd = msg->getGcount();
  if (_lb_args.debug()>=2) {
    CkPrintf("MultistepLB: receiveCentroids started: %d elements, msg length: %d\n", msg->getGcount(), msg->getLength()); 
  }
  haveTPCentroids = true;
  tpmsg = msg;
  treeProxy.doAtSync();
  if (_lb_args.debug()>=2) {
    CkPrintf("MultistepLB: receiveCentroids done\n");  
  }
}

bool MultistepLB::QueryBalanceNow(int step){
  if(step == 0){
    return false; 
  }
  if (_lb_args.debug()>=1) {
    if(CkMyPe() == 0){
      CkPrintf("MultistepLB: Step %d\n", step);
    }
  }
  return true;

}

// helper functions for multistepping
#ifdef MCLBMS

// determine phase based on lastActiveRung as saved in map.tpCentroids
unsigned int MultistepLB::determinePhase(unsigned int lastActiveRung){
  return lastActiveRung;
}

// merge data instrumented in previous iteration of computation with data from earlier iterations
// this data must be stored according to the tree piece 
// index of the treepiece (obtained from treepiece[.].tag)
// since the ldb index of an object changes from iteration
// to iteration.
void MultistepLB::mergeInstrumentedData(int phase, BaseLB::LDStats *stats){

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
      savedPhaseStats.push_back(LightweightLDStats1());
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
bool MultistepLB::havePhaseData(int phase){
  return (savedPhaseStats.length() > phase && savedPhaseStats[phase].n_objs > 0);
}

void MultistepLB::printData(BaseLB::LDStats &stats, int phase, int *revObjMap){
  int i;
  
  CkPrintf("---- data (%d): %d objects ----\n", phase, stats.n_objs);
  for(i = 0; i < stats.n_objs; i++){
     CkPrintf("%d: %g\n", i, 
	       stats.objData[i].wallTime);
  }
  CkPrintf("---- end data (%d) ----\n", phase);
}

void MultistepLB::makeActiveProcessorList(BaseLB::LDStats *stats, int numActiveObjs){
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

#define LARGE_PHASE_THRESHOLD 0.10

void MultistepLB::work(BaseLB::LDStats* stats)
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
  //if(true){
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
void MultistepLB::greedy(BaseLB::LDStats *stats, int count, int phase, int prevPhase){

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

void MultistepLB::work2(BaseLB::LDStats *stats, int count, int phase, int prevPhase){
  int numobjs = stats->n_objs;
  int nmig = stats->n_migrateobjs;

  if (_lb_args.debug()>=2) {
    CkPrintf("[work2] %d objects allocating %d bytes for tp\n", nmig, nmig*sizeof(TPObject));
  }
  CkPrintf("[ORB3D] objects total %d active %d\n", numobjs,nmig);

  // this data structure is used by the orb3d strategy
  // to balance objects. it is NOT indexed by tree piece index
  // there are as many entries in it as there are
  // migratable (active) tree pieces
  TPObject *tp_array = new TPObject[nmig];

  if (_lb_args.debug()>=2) {
    CkPrintf("[work2] ready tp_array data structure\n");
  }

  int j = 0;
  for(int i = 0; i < numobjs; i++){
    int lb = tpCentroids[i].tag;

    if(!stats->objData[lb].migratable) continue;
    tp_array[j].centroid.x = tpCentroids[i].vec.x;
    tp_array[j].centroid.y = tpCentroids[i].vec.y;
    tp_array[j].centroid.z = tpCentroids[i].vec.z;
    tp_array[j].migratable = true; 
    if(step() == 0){
      tp_array[j].load = tpCentroids[i].myNumParticles;
    }
    else{
      tp_array[j].load = stats->objData[lb].wallTime;
    }
    tp_array[j].lbindex = lb;
    j++;
  }
  CkAssert(j==nmig);

  mapping = &stats->to_proc;
  int dim = 0;
  TopoManager tmgr;

  procsPerNode = tmgr.getDimNT();

  int nx = tmgr.getDimNX();
  int ny = tmgr.getDimNY();
  int nz = tmgr.getDimNZ();
  int numnodes = nx*ny*nz; 

  if (_lb_args.debug()>=2) {
    CkPrintf("[work2] %d numnodes allocating %d bytes for nodes\n", numnodes, numnodes*sizeof(Node));
  }
  Node *nodes = new Node[numnodes];

  for(int i = 0; i < stats->count; i++){
    int t;
    int x,y,z;
    int node;
    tmgr.rankToCoordinates(i,x,y,z,t);
    
    node = z*nx*ny + y*nx + x; 
    nodes[node].x = x;
    nodes[node].y = y;
    nodes[node].z = z;
    nodes[node].procRanks.push_back(i);
    //CkPrintf("node %d,%d,%d (%d) gets t %d\n", nodes[node].x, nodes[node].y, nodes[node].z, node, t);
  }

  if (_lb_args.debug()>=2) {
    CkPrintf("[work2] map\n");
  }
  map(tp_array,nmig,numnodes,nodes,nx,ny,nz,dim);

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


  delete[] objload;
  delete[] tp_array;
  delete[] nodes;


}

void MultistepLB::map(TPObject *tp, int ntp, int nn, Node *nodes, int xs, int ys, int zs, int dim){
  //CkPrintf("ntp: %d np: %d dim: %d path: 0x%x\n",ntp,np,dim,path);
  if(nn == 1){
    directMap(tp,ntp,nodes);
  }
  else{
    int totalTp = ntp;
    
    qsort(tp,ntp,sizeof(TPObject),compares[dim]);
    qsort(nodes,nn,sizeof(Node),pc[dim]);
    // tp and ntp are modified to hold the particulars of
    // the left/dn/near partitions
    // tp2 and totalTp-ntp hold the objects in the 
    // right/up/far partitions
    TPObject *tp2 = partitionEvenLoad(tp,ntp);
    Node *nodes2 = halveNodes(nodes,nn);
    int d = nextDim(dim,xs,ys,zs); 
    if(d == 0){
      xs >>= 1;
    }
    else if(d == 1){
      ys >>= 1;
    }
    else{
      zs >>= 1;
    }
    map(tp,ntp,nn/2,nodes,xs,ys,zs,d);
    map(tp2,totalTp-ntp,nn/2,nodes2,xs,ys,zs,d);
  }
}

#define ZERO_THRESHOLD 0.00001

void MultistepLB::directMap(TPObject *tp, int ntp, Node *nodes){

  float load = 0.0;
  for(int i = 0; i < ntp; i++){
    //CkPrintf("obj %d thisindex %d %d %f %f %f %f to node %d %d %d\n", tp[i].lbindex, tp[i].index, tp[i].nparticles, tp[i].load, tp[i].centroid.x, tp[i].centroid.y, tp[i].centroid.z, nodes[0].x, nodes[0].y, nodes[0].z);
    load += tp[i].load;
  }
  //CkPrintf("node %d %d %d total load %f\n", nodes[0].x, nodes[0].y, nodes[0].z, load);
  
  std::priority_queue<TPObject> pq_obj;
  std::priority_queue<Processor> pq_proc;

  for(int i = 0; i < ntp; i++){
    pq_obj.push(tp[i]);
  }

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
    if(tp.load < ZERO_THRESHOLD){
      (*mapping)[tp.lbindex] = nodes[0].procRanks[currentZeroProc];
      currentZeroProc = currentZeroProc+1;
      if(currentZeroProc == procsPerNode){
        currentZeroProc = 0;
      }
    }
    else{
      // if object has some non-zero load, assign it to a proc greedily
      Processor p = pq_proc.top();
      pq_proc.pop();

      //CkPrintf("proc %d load %f gets obj %d load %f\n", p.t, p.load, tp.lbindex, tp.load);

      p.load += tp.load;
      (*mapping)[tp.lbindex] = nodes[0].procRanks[p.t];

      pq_proc.push(p);
    }
  }
}

int MultistepLB::nextDim(int dim_, int xs, int ys, int zs){
  int max = xs;
  int dim = 0;
  if(max < ys){
    max = ys;
    dim = 1;
  }
  if(max < zs){
    max = zs;
    dim = 2;
  }
  return dim; 
}




































































































#define LOAD_EQUAL_TOLERANCE 1.02

TPObject *MultistepLB::partitionEvenLoad(TPObject *tp, int &ntp){
  float totalLoad = 0.0;
  for(int i = 0; i < ntp; i++){
    totalLoad += tp[i].load;
  }
  float lload = 0.0;
  float rload = totalLoad;
  float prevDiff = lload-rload;
  if(prevDiff < 0.0){
    prevDiff = -prevDiff;
  }

  int consider;
  for(consider = 0; consider < ntp;){
    float newll = lload + tp[consider].load;
    float newrl = rload - tp[consider].load;
    
    float newdiff = newll-newrl;
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

  ntp = consider;
  return (tp+consider);
}

Node *MultistepLB::halveNodes(Node *start, int np){
  Node *ret = start;
  ret = start+np/2;
  return ret;
}

void MultistepLB::pup(PUP::er &p){
  CentralLB::pup(p);
  if(p.isPacking() && haveTPCentroids){
    // if checkpointing, no need to 
    // keep around the centroid message
    delete tpmsg;
    haveTPCentroids = false;
  }
  p | haveTPCentroids; 
  p | procsPerNode;
  p | savedPhaseStats;
}

void LightweightLDStats1::pup(PUP::er &p){
  p|n_objs;
  p|n_migrateobjs;
  p|objData;
}

#include "MultistepLB.def.h"
