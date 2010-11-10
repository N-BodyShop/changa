#include <charm++.h>
#include "cklists.h"
#include "MultistepLB.h"
#include "TopoManager.h"
#include "ParallelGravity.h"
#include "Vector3D.h"
#include <queue>

extern CProxy_TreePiece treeProxy;
using namespace std;

extern double msLargePhaseThreshold;
extern double msExpandFactor;

CreateLBFunc_Def(MultistepLB, "Works best with multistepped runs; uses OrbSmooth for larger steps, RoundRobin otherwise");


MultistepLB::MultistepLB(const CkLBOptions &opt): CentralLB(opt)
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
    delete tpmsg; 
  }

  tpCentroids = (TaggedVector3D *)msg->getData();
  if(_lb_args.debug() >= 1){
    CkPrintf("MultistepLB: receiveCentroids started: %d elements, msg length: %d\n", msg->getGcount(), msg->getLength()); 
  }
  haveTPCentroids = true;
  tpmsg = msg;
  /*
  map.tpCentroids.free();
  
  while(i < msg->getGcount()){
     map.tpCentroids.push_back(*cur);
     cur = cur + 1;
     i++;
  }
  */
  treeProxy.doAtSync();
  if(_lb_args.debug() >= 1){
    CkPrintf("MultistepLB: receiveCentroids done\n");  
  }
  //delete msg;
}

//jetley
CmiBool MultistepLB::QueryBalanceNow(int step){
  if(CkMyPe() == 0){
    CkPrintf("MultistepLB: Step %d\n", step);
  }
  if(step == 0){
    return false;
  }
  else{
    return true;
  }
}






// helper functions for multistepping
#ifdef MCLBMS

// determine phase based on lastActiveRung as saved in map.tpCentroids
unsigned int MultistepLB::determinePhase(unsigned int lastActiveRung){
  return lastActiveRung;
}

// merge data instrumented in previous iteration of computation with data from earlier iterations
void MultistepLB::mergeInstrumentedData(int phase, BaseLB::LDStats *stats, int *lbToTp){

  int len;
  int whichPos;

  // tune alpha as needed - this is the merge parameter
  double alpha = 0.0;
  
  if(phase == -1){
#ifdef MCLBMSV
    CkPrintf("phase = -1, discarding\n");
#endif
    return;
  }
  
  len = savedPhaseStats.length();
  
  // must extend array
  if(phase > len-1){
    int numAdditional = phase-len+1;
    while(numAdditional > 0){
      CondensedLDStats clds;
      savedPhaseStats.push_back(clds);
#ifdef MCLBMSV
      CkPrintf("Making new entry for phase %d (%d)\n", savedPhaseStats.length()-1, phase);
#endif
      savedPhaseStats[savedPhaseStats.length()-1].n_objs = 0;  // so that we can identify this phase's entry as being unfilled when we get the data for it
      numAdditional--;
    }
    len = savedPhaseStats.length();
    whichPos = len-1;
    CondensedLDStats &condensed = savedPhaseStats[whichPos];
    // initialize 
    condensed.init(stats, lbToTp);
  }
  // enough space for this phase 
  else{ 
    whichPos = phase;
    CondensedLDStats &condensed = savedPhaseStats[whichPos];
    // haven't yet populated this phase
    if(condensed.n_objs == 0){       
      // initialize 
      condensed.init(stats, lbToTp);
    }
    else{        
      // merge current with previous data
      condensed.merge(stats, lbToTp, alpha, true);
    }
  }
}

// whether we have instrumented data for this phase
bool MultistepLB::havePhaseData(int phase){
  return (savedPhaseStats.length() > phase && savedPhaseStats[phase].n_objs > 0);
}

void MultistepLB::printData(CondensedLDStats *stats, int phase){
  int i;
  
  CkPrintf("---- data phase %d %d objects ----\n", phase, stats->n_objs);
  for(i = 0; i < stats->n_objs; i++){
     CkPrintf("%d %d %f\n", i, stats->migratable[i], stats->wallTime[i]);
  }
  CkPrintf("---- end data ----\n", phase);
}

void MultistepLB::makeActiveProcessorList(CondensedLDStats *stats, int na, int n, bool largePhase){
  double expandFactor = msExpandFactor;
  if(!largePhase){
    int procsNeeded = expandFactor*(((float)na*stats->count)/n);
    if(procsNeeded > stats->count){
      procsNeeded = stats->count;
    }
    else if(procsNeeded < 1){
      procsNeeded = 1;
    }

    if(_lb_args.debug() >= 0){
      CkPrintf("small: Processors 0 to %d active %d total %d\n", procsNeeded-1, na, n);
    }
    stats->count = procsNeeded;
  }
  else{
    CkPrintf("large: Processors 0 to %d active %d total %d\n", stats->count-1, na, n);
    // don't modify stats->count
  }
}
#endif


void MultistepLB::work(BaseLB::LDStats* stats, int count)
{
#if CMK_LBDB_ON
  // find active objects - mark the inactive ones as non-migratable
  int i;
  
  stats->makeCommHash();
  int *lbToTp = new int[stats->n_objs];
  for(int i = 0; i < stats->n_objs; i++){
    LDObjHandle &handle = tpCentroids[i].handle;
    tpCentroids[i].tag = stats->getHash(handle.id, handle.omhandle.id);
    lbToTp[tpCentroids[i].tag] = tpCentroids[i].tpindex;
  }

  int phase = determinePhase(tpCentroids[0].activeRung);
  int prevPhase = tpCentroids[0].prevActiveRung;

  int numActiveObjects = 0;
  int numInactiveObjects = 0;

  // to calculate ratio of active particles in phase
  int numActiveParticles = 0;
  int totalNumParticles = 0;
  
  for(int i = 0; i < stats->n_objs; i++){
    stats->to_proc[i] = stats->from_proc[i];
  }
  // update phase data 
  if(_lb_args.debug() >= 1){
    CkPrintf("merging previous phase %d data; current phase: %d\n", prevPhase, phase);
  }
  mergeInstrumentedData(prevPhase, stats, lbToTp); 

  bool useRatios;
  CondensedLDStats *statsToUse = NULL;
  if(havePhaseData(phase)){
    statsToUse = &(savedPhaseStats[phase]); 
    useRatios = false;
  }
  else if(havePhaseData(0)){
    statsToUse = &(savedPhaseStats[0]);
    useRatios = true;
  }
  else{
    return;
  }

  
  float *ratios = new float[statsToUse->n_objs];

  // update migratable info in current stats 
  // also record the ratios of active to total
  // particles for each object
  for(int i = 0; i < statsToUse->n_objs; i++){
    int np = tpCentroids[i].myNumParticles;
    int nap = tpCentroids[i].numActiveParticles;
    int tpindex = tpCentroids[i].tpindex;
    if(np == 0 || nap == 0){
      numInactiveObjects++;
      ratios[tpindex] = 0.0;
      if(statsToUse->migratable[tpindex]){
        statsToUse->migratable[tpindex] = 0;
        statsToUse->n_migrateobjs--;
      }
    }
    else{
      ratios[tpindex] = nap/(float)np;
      numActiveObjects++;
    }
    numActiveParticles += nap; 
    totalNumParticles += np; 
  }

  if(_lb_args.debug() >= 2){
    CkPrintf("LDSTATS:\n");
    for(int i = 0; i < statsToUse->n_objs; i++){
      int lbindex = tpCentroids[i].tag;
      int tpindex = tpCentroids[i].tpindex;
      CkPrintf("work1 %d %d %d %f\n", tpindex, lbindex, statsToUse->migratable[tpindex], statsToUse->wallTime[tpindex]);
    }
  }


#ifdef MCLBMSV
  CkPrintf("numActiveObjects: %d, numInactiveObjects: %d\n", numActiveObjects, numInactiveObjects);
#endif

/*
  if(havePhaseData(phase)){
#ifdef MCLBMSV
    CkPrintf("phase %d data available\n", phase);
#endif
    for(i = 0; i < stats->n_objs; i++){
      stats->objData[i].cpuTime = savedPhaseStats[phase].cpuTime[i];
      stats->objData[i].wallTime = savedPhaseStats[phase].wallTime[i];
    }
  }
  else if(havePhaseData(0)){
#ifdef MCLBMSV
    CkPrintf("phase %d data unavailable, using phase 0 loads\n", phase);
#endif
    for(i = 0; i < stats->n_objs; i++){
      stats->objData[i].cpuTime = ratios[i]*savedPhaseStats[0].objData[i].cpuTime;
      stats->objData[i].wallTime = ratios[i]*savedPhaseStats[0].objData[i].wallTime;
    }
  }
  else{
#ifdef MCLBMSV
    CkPrintf("phase %d data unavailable\n", phase);
#endif
    delete[] ratios;
    return;
  }
  */


  // select processors
  float activeRatio = (float)numActiveParticles/totalNumParticles;
  bool largePhase = activeRatio > msLargePhaseThreshold;
  makeActiveProcessorList(statsToUse, numActiveParticles, totalNumParticles, largePhase);
  int savedCount = count;
  count = statsToUse->count;

  // let the strategy take over on merged data and processor information
  if(largePhase){
    if(_lb_args.debug() >= 1){
      CkPrintf("******** BIG STEP activeRatio %f activeObjects %d processors %d *********!\n", activeRatio, numActiveObjects, statsToUse->count);
    }
    work2(statsToUse, ratios, count, &stats->to_proc, useRatios, phase);
  }     // end if large phase
  else{ 
    //fprintf(stderr, "before alloc\n");
    CmiMemoryCheck();
    SmallPhaseObject *tp = new SmallPhaseObject[statsToUse->n_migrateobjs];
    int j = 0;
    // iterate through BaseLB::LDStats. therefore,
    // order to get the load information from savedPhaseStats,
    // we will need to translate from lbindex to tpindex
    if(_lb_args.debug() >= 1){
      fprintf(stderr, "smallPhaseGreedy activeRatio %f activeObjects %d processors %d\n", activeRatio, numActiveObjects, statsToUse->count);
      fflush(stderr);
    }
    for(int i = 0; i < stats->n_objs; i++){
      int tpindex = lbToTp[i]; 
      if(statsToUse->migratable[tpindex] == 0){
        continue;
      }
      else{
        tp[j].load = statsToUse->wallTime[tpindex];
        tp[j].lbindex = i;
        tp[j].tpindex = tpindex;
        if(_lb_args.debug() >= 2){
          fprintf(stderr, "smallPhaseGreedy tpindex %d lbindex %d load %f\n", tpindex, i, tp[j].load);
          fflush(stderr);
        }
        j++;
      }
    }
    //fprintf(stderr, "after copies\n");
    CmiMemoryCheck();
    greedySmallPhase(tp, statsToUse->n_migrateobjs, statsToUse->count, &(stats->to_proc));
    //fprintf(stderr, "after greedySmallPhase\n");
    CmiMemoryCheck();
    /*
    for(int i = 0; i < stats->n_objs; i++){
      if(stats->objData[i].migratable == 0) continue;
      stats->to_proc[i] = i%stats->count;
    }
    */
    delete[] tp;
    //fprintf(stderr, "after delete\n");
    CmiMemoryCheck();
  }
  if(_lb_args.debug() >= 2){
    fprintf(stderr, "\n\nLB decisions:\n\n");
    for(int i = 0; i < stats->count; i++){
      fprintf(stderr, "tpindex %d lbindex %d from %d to %d\n", lbToTp[i], i, stats->from_proc[i], stats->to_proc[i]);
    }
  }

  for(int i = 0; i < stats->count; i++){
    fprintf(stderr, "procload step %d proc %d nobjs %d walltime %f cputime %f idletime %f\n", 
                                    step()-1, 
                                    i, 
                                    stats->procs[i].n_objs, 
                                    stats->procs[i].total_walltime,
                                    stats->procs[i].total_cputime, 
                                    stats->procs[i].idletime);
  }
  delete []ratios;
  delete []lbToTp;
#endif //CMK_LDB_ON

}

void MultistepLB::greedySmallPhase(SmallPhaseObject *tp, int ntp, int np, CkVec<int> *mapping){
  std::priority_queue<SmallPhaseObject> pq_obj;
  std::priority_queue<Processor> pq_proc;

  for(int i = 0; i < ntp; i++){
    pq_obj.push(tp[i]);
  }

  for(int i = 0; i < np; i++){
    Processor p;
    p.load = 0.0;
    p.t = i;
    pq_proc.push(p);
  }

  while(!pq_obj.empty()){
    SmallPhaseObject obj = pq_obj.top();
    pq_obj.pop();

    Processor p = pq_proc.top();
    pq_proc.pop();

    p.load += obj.load;
    if(_lb_args.debug() >= 2){
      fprintf(stderr,"smallPhaseGreedy tpindex %d lbindex %d load %f to rank %d weight %f\n", obj.tpindex, obj.lbindex, obj.load, p.t, p.load);
      fflush(stderr);
    }
    (*mapping)[obj.lbindex] = p.t;
    pq_proc.push(p);
  }

  if(_lb_args.debug() >= 1){
    fprintf(stderr, "smallPhaseGreedy PROCESSOR LOADS\n");
    while(!pq_proc.empty()){
      Processor p = pq_proc.top();
      fprintf(stderr, "rank %d load %f\n", p.t, p.load);
      pq_proc.pop();
    }
  }

}

//**************************************
// ORB3DLB functions
//**************************************

// use ratios only when data isn't available 
// for a particular phase.
void MultistepLB::work2(CondensedLDStats *stats, float *ratios, int count, CkVec<int> *m, bool useRatios, int phase){
  int numobjs = stats->n_objs;
  int nmig = stats->n_migrateobjs;

  if(_lb_args.debug() >= 1){
    CkPrintf("[work2] %d objects allocating %d bytes for tp\n", nmig, nmig*sizeof(TPObject));
  }
  TPObject *tp = new TPObject[nmig];

  if(_lb_args.debug() >= 1){
    CkPrintf("[work2] ready tp data structure phase %d\n", phase);
  }

  int j = 0;
  for(int i = 0; i < stats->n_objs; i++){
    int tpindex = tpCentroids[i].tpindex;
    int lbindex = tpCentroids[i].tag;

    if(!stats->migratable[tpindex]) continue;
    tp[j].centroid.x = tpCentroids[i].vec.x;
    tp[j].centroid.y = tpCentroids[i].vec.y;
    tp[j].centroid.z = tpCentroids[i].vec.z;
    tp[j].migratable = stats->migratable[tpindex];

    if(step() == 0){
      tp[j].load = tpCentroids[i].myNumParticles;
    }
    else if(useRatios){
      tp[j].load = ratios[tpindex]*stats->wallTime[tpindex];
    }
    else{
      tp[j].load = stats->wallTime[tpindex];
    }
    tp[j].lbindex = lbindex;
    tp[j].tpindex = tpindex;
    if(_lb_args.debug() >= 2){
      CkPrintf("work2 %d %d %d %f\n", tpindex, lbindex, stats->migratable[tpindex], tp[j].load);
    }
    j++;
  }


  mapping = m;
  int dim = 0;
  TopoManager tmgr;

  procsPerNode = tmgr.getDimNT();

  int nx = tmgr.getDimNX();
  int ny = tmgr.getDimNY();
  int nz = tmgr.getDimNZ();
  int numnodes = nx*ny*nz; 

  if(_lb_args.debug() >= 1){
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
    nodes[node].num = i+1;
  }

  if(_lb_args.debug() >= 1){
    CkPrintf("[work2] map\n");
  }
  map(tp,nmig,numnodes,nodes,nx,ny,nz,dim);

  /*
  int migr = 0;
  float *objload = new float[stats->count];
  for(int i = 0; i < stats->count; i++){
    objload[i] = 0.0;
  }
  for(int i = 0; i < numobjs; i++){
    objload[stats->from_proc[i]] += stats->objData[i].wallTime;
    if(stats->to_proc[i] != stats->from_proc[i]) migr++;
  }
  
  CkPrintf("***************************\n");
  CkPrintf("Before LB step %d\n", step());
  CkPrintf("***************************\n");
  CkPrintf("i pe wall cpu idle bg_wall bg_cpu objload\n");
  for(int i = 0; i < stats->count; i++){
    CkPrintf("[pestats] %d %d %f %f %f %f %f %f\n", 
                               i,
                               stats->procs[i].pe, 
                               stats->procs[i].total_walltime, 
                               stats->procs[i].total_cputime, 
                               stats->procs[i].idletime,
                               stats->procs[i].bg_walltime,
                               stats->procs[i].bg_cputime,
                               objload[i]);
  }
  CkPrintf("%d objects migrating\n", migr);
  */

  //delete[] objload;
  delete[] tp;
  delete[] nodes;


}

void MultistepLB::map(TPObject *tp, int ntp, int nn, Node *nodes, int xs, int ys, int zs, int dim){

  if(_lb_args.debug() >= 2){
    CkPrintf("ntp %d nn %d dim %d\n",ntp,nn,dim);
  }
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
  if(_lb_args.debug() >= 1){
    CkPrintf("directMap %d objects\n", ntp);
    if(_lb_args.debug() >= 2){
      CkPrintf("nodes 0x%x nodes.procRanks vec 0x%x size %d init %d\n", nodes, nodes[0].procRanks.getVec(), nodes[0].procRanks.length(), nodes[0].num);
      if(nodes == NULL || nodes[0].procRanks.getVec() == NULL || nodes[0].procRanks.length() == 0){
        CkAbort("bad nodes\n");
      }
      for(int i = 0; i < nodes[0].procRanks.length(); i++){
        CkPrintf("%d\n", nodes[0].procRanks[i]);
      }
    }
  }
  for(int i = 0; i < ntp; i++){
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
      if(_lb_args.debug() >= 2){
        CkPrintf("directMap tpindex %d lbindex %d to rank %d\n", tp.tpindex, tp.lbindex, nodes[0].procRanks[p.t]);
      }
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
  float partition1Load = 0.0;
  float totalLoad = 0.0;
  for(int i = 0; i < ntp; i++){
    totalLoad += tp[i].load;
  }
  float halfLoad = 0.5*totalLoad;
  //CkPrintf("partitionEvenLoad total load %f half load %f\n", totalLoad, halfLoad);
  int split = -1;

  for(int i = 0; i < ntp; i++){
    // if including this element in partition1 brings us closer
    // to halfLoad, do it
    //if((partition1Load+tp[i].load-halfLoad) < (halfLoad-partition1Load)){
    if((partition1Load+tp[i].load <= halfLoad) ||
       (partition1Load < halfLoad && partition1Load+tp[i].load > halfLoad)){
      partition1Load += tp[i].load;
      split++;
    }
    else{
      break;
    }
  }

  //float lload = 0.0;
  //for(int i = 0; i < split+1; i++){
  //  lload += tp[i].load;
  //}
  //float rload = totalLoad-lload;

  //CkPrintf("partitionEvenLoad partition1Load %f lsplit %f rsplit %f\n", partition1Load, lload, rload);

  ntp = split+1;
  return (tp+split+1);
}

Node *MultistepLB::halveNodes(Node *start, int np){
  Node *ret = start;
  ret = start+np/2;
  return ret;
}



#include "MultistepLB.def.h"

/*@}*/
