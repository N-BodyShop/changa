#include <charm++.h>
#include "cklists.h"
#include "MultistepLB.h"
#include "TopoManager.h"
#include "ParallelGravity.h"
#include "Vector3D.h"
#include <queue>

extern CProxy_TreePiece treeProxy;
using namespace std;

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
  CkPrintf("MultistepLB: receiveCentroids started: %d elements, msg length: %d\n", msg->getGcount(), msg->getLength()); 
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
  CkPrintf("MultistepLB: receiveCentroids done\n");  
  //delete msg;
}

//jetley
CmiBool MultistepLB::QueryBalanceNow(int step){
  if(CkMyPe() == 0)
    CkPrintf("MultistepLB: Step %d\n", step);
  if(step == 0){
    return false;
  }
  else{
    return true;
  }
}

void MultistepLB::rec_divide(int n, Partition &p)
{
  int i;
  int midpos;
  int n1, n2;
  double load1, currentload;
  int maxdir, count;
  Partition p1, p2;

  if (_lb_args.debug()>=2) {
    CmiPrintf("rec_divide starts: partition n:%d count:%d load:%f (%d %d %d, %d %d %d)\n", n, p.count, p.load, p.origin[0], p.origin[1], p.origin[2], p.corner[0], p.corner[1], p.corner[2]);
  }

  if (n==1) {		// we are done in this branch
    partitions[currentp++] = p;
    return;
  }
/*
  if (p.origin.x==p.corner.x && p.origin.y==p.corner.y && p.origin.z==p.corner.z) 
     CmiAbort("AlgRecBisection failed in recursion.\n"); 
*/
  if (_lb_args.debug()>=2) {
    CmiPrintf("{\n");
  }

  // divide into n1 and n2 two subpartitions
  n2 = n/2;
  n1 = n-n2;

  // subpartition n1 should have this load
  load1 = (1.0*n1/n) * p.load;
  if (_lb_args.debug()>=2)
    CmiPrintf("goal: n1: %d with load1: %f; n2: %d load2: %f\n", n1, load1, n2, p.load-load1);

  p1 = p;
  p1.refno = ++refno;
  p1.bkpes.resize(0);

  p2 = p;
  p2.refno = ++refno;
  p2.bkpes.resize(0);

  // determine the best division direction
  int maxSpan=-1;
  maxdir = XDIR;
  for (i=XDIR; i<=ZDIR; i++) {
    int myspan = p.corner[i] - p.origin[i];
    if (myspan > maxSpan) {
      maxdir = i;
      maxSpan = myspan;
    }
  }

  // other two dimensions
  int dir2 = (maxdir+1)%3;
  int dir3 = (maxdir+2)%3;

  currentload = 0.0;
  // counting background load
  if (!_lb_args.ignoreBgLoad()) {
    CmiAssert(p.bkpes.size() == n);
    // first n1 processors
    for (i=0; i<n1; i++) currentload += statsData->procs[p.bkpes[i]].bg_walltime;
  }

  count = 0;
  midpos = p.origin[maxdir];
  for (i=0; i<nObjs; i++) {
    // not belong to this partition
    if (computeLoad[vArray[maxdir][i].id].refno != p.refno) continue;
    if (vArray[maxdir][i].v<p.origin[maxdir]) continue;
    if (vArray[maxdir][i].v>p.corner[maxdir]) break;

    int cid = vArray[maxdir][i].id;	// this compute ID
    // check if this compute is within the partition
    if ( computeLoad[cid].v[dir2] >= p.origin[dir2] &&
	 computeLoad[cid].v[dir2] <= p.corner[dir2] &&
	 computeLoad[cid].v[dir3] >= p.origin[dir3] &&
	 computeLoad[cid].v[dir3] <= p.corner[dir3]  ) {
      // this compute is set to the first partition
      if (currentload <= load1) {
	computeLoad[cid].refno = p1.refno;
        currentload += computeLoad[cid].load;
        count ++;
	midpos = computeLoad[cid].v[maxdir];
      }
      else {	// or the next partition
	computeLoad[cid].refno = p2.refno;
      }
    }
  }
#ifdef DEBUG
//  CmiPrintf("X:cur:%d, prev:%d load:%f %f\n", cur, prev, currentload, prevload);
  CmiPrintf("DIR:%d %d load:%f\n", maxdir, midpos, currentload);
#endif

  p1.corner[maxdir] = midpos;
  p2.origin[maxdir] = midpos;

  p1.load = currentload;
  p1.count = count;
  p2.load = p.load - p1.load;
  p2.count = p.count - p1.count;

  // assign first n1 copy of background to p1, and rest to p2
  if (!_lb_args.ignoreBgLoad()) {
    for (i=0; i<n; i++)
      if (i<n1) p1.bkpes.push_back(p.bkpes[i]);
      else p2.bkpes.push_back(p.bkpes[i]);
  }

  if (_lb_args.debug()>=2) {
    CmiPrintf("p1: n:%d count:%d load:%f\n", n1, p1.count, p1.load);
    CmiPrintf("p2: n:%d count:%d load:%f\n", n2, p2.count, p2.load);
    CmiPrintf("}\n");
  }

  rec_divide(n1, p1);
  rec_divide(n2, p2);
}

void MultistepLB::setVal(int x, int y, int z)
{
  int i;
  for (i=0; i<nObjs; i++) {
    computeLoad[i].tv = 1000000.0*computeLoad[i].v[x]+
			1000.0*computeLoad[i].v[y]+
			computeLoad[i].v[z];
  }
#if 0
  CmiPrintf("original:%d\n", x);
  for (i=0; i<numComputes; i++) 
    CmiPrintf("%d ", computeLoad[i].tv);
  CmiPrintf("\n");
#endif
}

int MultistepLB::sort_partition(int x, int p, int r)
{
  double mid = computeLoad[vArray[x][p].id].tv;
  int i= p;
  int j= r;
  while (1) {
    while (computeLoad[vArray[x][j].id].tv > mid && j>i) j--;
    while (computeLoad[vArray[x][i].id].tv < mid && i<j) i++;
    if (i<j) {
      if (computeLoad[vArray[x][i].id].tv == computeLoad[vArray[x][j].id].tv)
      {
	if (computeLoad[vArray[x][i].id].tv != mid) CmiAbort("my god!\n");
	if (i-p < r-j) i++;
	else j--;
	continue;
      }
      VecArray tmp = vArray[x][i];
      vArray[x][i] = vArray[x][j];
      vArray[x][j] = tmp;
    }
    else
      return j;
  }
}

void MultistepLB::qsort1(int x, int p, int r)
{
  if (p<r) {
    int q = sort_partition(x, p, r);
//CmiPrintf("midpoint: %d %d %d\n", p,q,r);
    qsort1(x, p, q-1);
    qsort1(x, q+1, r);
  }
}

void MultistepLB::quicksort(int x)
{
  int y = (x+1)%3;
  int z = (x+2)%3;
  setVal(x, y, z);
  qsort1(x, 0, nObjs-1);

#if 0
  CmiPrintf("result for :%d\n", x);
  for (int i=0; i<nObjs; i++) 
    CmiPrintf("%d ", computeLoad[vArray[x][i].id].tv);
  CmiPrintf("\n");
#endif
}

void MultistepLB::mapPartitionsToNodes()
{
  int i,j;
#if 1
  if (!_lb_args.ignoreBgLoad()) {
      // processor mapping has already been determined by the background load pe
    for (i=0; i<npartition; i++) partitions[i].node = partitions[i].bkpes[0];
  }
  else {
    for (i=0; i<P; i++) partitions[i].node = i;
  }
#else
  PatchMap *patchMap = PatchMap::Object();

  int **pool = new int *[P];
  for (i=0; i<P; i++) pool[i] = new int[P];
  for (i=0; i<P; i++) for (j=0; j<P; j++) pool[i][j] = 0;

  // sum up the number of nodes that patches of computes are on
  for (i=0; i<numComputes; i++)
  {
    for (j=0; j<P; j++)
      if (computeLoad[i].refno == partitions[j].refno) 
      {
	int node1 = patchMap->node(computes[i].patch1);
	int node2 = patchMap->node(computes[i].patch2);
	pool[j][node1]++;
	pool[j][node2]++;
      }
  }
#ifdef DEBUG
  for (i=0; i<P; i++) {
    for (j=0; j<P; j++) CmiPrintf("%d ", pool[i][j]);
    CmiPrintf("\n");
  }
#endif
  while (1)
  {
    int index=-1, node=0, eager=-1;
    for (j=0; j<npartition; j++) {
      if (partitions[j].node != -1) continue;
      int wantmost=-1, maxnodes=-1;
      for (k=0; k<P; k++) if (pool[j][k] > maxnodes && !partitions[k].mapped) {wantmost=k; maxnodes = pool[j][k];}
      if (maxnodes > eager) {
	index = j; node = wantmost; eager = maxnodes;
      }
    }
    if (eager == -1) break;
    partitions[index].node = node;
    partitions[node].mapped = 1;
  }

  for (i=0; i<P; i++) delete [] pool[i];
  delete [] pool;
#endif

/*
  if (_lb_args.debug()) {
    CmiPrintf("partition load: ");
    for (i=0; i<npartition; i++) CmiPrintf("%f ", partitions[i].load);
    CmiPrintf("\n");
    CmiPrintf("partitions to nodes mapping: ");
    for (i=0; i<npartition; i++) CmiPrintf("%d ", partitions[i].node);
    CmiPrintf("\n");
  }
*/
  if (_lb_args.debug()) {
    CmiPrintf("After partitioning: \n");
    for (i=0; i<npartition; i++) {
      double bgload = 0.0;
      if (!_lb_args.ignoreBgLoad())
        bgload = statsData->procs[partitions[i].bkpes[0]].bg_walltime;
      CmiPrintf("[%d=>%d] (%d,%d,%d) (%d,%d,%d) load:%f count:%d objload:%f\n", i, partitions[i].node, partitions[i].origin[0], partitions[i].origin[1], partitions[i].origin[2], partitions[i].corner[0], partitions[i].corner[1], partitions[i].corner[2], partitions[i].load, partitions[i].count, partitions[i].load-bgload);
    }
    for (i=npartition; i<P; i++) CmiPrintf("[%d] --------- \n", i);
  }

}

// helper functions for multistepping
#ifdef MCLBMS

// determine phase based on lastActiveRung as saved in map.tpCentroids
unsigned int MultistepLB::determinePhase(unsigned int lastActiveRung){
  return lastActiveRung;
}

// merge data instrumented in previous iteration of computation with data from earlier iterations
void MultistepLB::mergeInstrumentedData(int phase, BaseLB::LDStats *stats){

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
    condensed.init(stats);
  }
  // enough space for this phase 
  else{ 
    whichPos = phase;
    CondensedLDStats &condensed = savedPhaseStats[whichPos];
    // haven't yet populated this phase
    if(condensed.n_objs == 0){       
      // initialize 
      condensed.init(stats);
    }
    else{        
      // merge current with previous data
      condensed.merge(stats, alpha, true);
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
     CkPrintf("%d %d %f %f \n", i, stats->tpindex[i], stats->cpuTime[i],
	       stats->wallTime[i]);
  }
  CkPrintf("---- end data ----\n", phase);
}

void MultistepLB::makeActiveProcessorList(CondensedLDStats *stats, int numActiveObjs, int numTotalObjects){
  int objsPerProc = 8;
  double expandFactor = 4.0;
  int procsNeeded = expandFactor*(((float)numActiveObjs)/numTotalObjects);
  if(procsNeeded > stats->count){
    procsNeeded = stats->count;
  }

#ifdef MCLBMSV
  CkPrintf("Processors 0 to %d active\n", procsNeeded-1);
#endif

  stats->count = procsNeeded;
}
#endif

#define LARGE_PHASE_THRESHOLD 0.1

void MultistepLB::work(BaseLB::LDStats* stats, int count)
{
#if CMK_LBDB_ON
  // find active objects - mark the inactive ones as non-migratable
  int i;
  
  stats->makeCommHash();
  for(i = 0; i < stats->n_objs; i++){
    LDObjHandle &handle = tpCentroids[i].handle;
    tpCentroids[i].tag = stats->getHash(handle.id, handle.omhandle.id);
  }

  //CkPrintf("LDSTATS:\n");
  //for(i = 0; i < stats->n_objs; i++){
  //  int lbindex = tpCentroids[i].tag;
  //  CkPrintf("%d %d %f\n", lbindex, tpCentroids[i].tpindex, stats->objData[lbindex].wallTime);
  //}

  int phase = determinePhase(tpCentroids[0].activeRung);
  int prevPhase = tpCentroids[0].prevActiveRung;

  int numActiveObjects = 0;
  int numInactiveObjects = 0;

  // to calculate ratio of active particles in phase
  int numActiveParticles = 0;
  int totalNumParticles = 0;
  
  for(i = 0; i < stats->n_objs; i++){
    stats->to_proc[i] = stats->from_proc[i];
  }
  // update phase data 
  CkPrintf("merging previous phase %d data; current phase: %d\n", prevPhase, phase);
  mergeInstrumentedData(prevPhase, stats); 

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
  for(i = 0; i < statsToUse->n_objs; i++){
    int np = tpCentroids[i].myNumParticles;
    int nap = tpCentroids[i].numActiveParticles;
    int lbindex = tpCentroids[i].tag;
    if(np == 0){
      numInactiveObjects++;
      ratios[lbindex] = 0.0;
      if(statsToUse->migratable[lbindex]){
        statsToUse->migratable[lbindex] = 0;
        statsToUse->n_migrateobjs--;
      }
    }
    else{
      ratios[lbindex] = nap/(float)np;
      numActiveObjects++;
    }
    numActiveParticles += nap; 
    totalNumParticles += np; 
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
  makeActiveProcessorList(statsToUse, numActiveObjects, numActiveObjects+numInactiveObjects);
  count = statsToUse->count;

  // let the strategy take over on merged data and processor information
  if((float)numActiveParticles/totalNumParticles > LARGE_PHASE_THRESHOLD){
    CkPrintf("******** BIG STEP *********!\n");
    work2(statsToUse, ratios, count, &stats->to_proc, useRatios, phase);
  }     // end if phase == 0
  else{ // not even greedy; round-robin
    for(i = 0; i < stats->n_objs; i++){
      if(stats->objData[i].migratable == 0) continue;
      stats->to_proc[i] = i%stats->count;
    }
  }
  delete []ratios;
#endif //CMK_LDB_ON

}

//**************************************
// ORB3DLB functions
//**************************************

// use ratios only when data isn't available 
// for a particular phase.
void MultistepLB::work2(CondensedLDStats *stats, float *ratios, int count, CkVec<int> *m, bool useRatios, int phase){
  int numobjs = stats->n_objs;
  int nmig = stats->n_migrateobjs;

  CkPrintf("[work2] %d objects allocating %d bytes for tp\n", nmig, nmig*sizeof(TPObject));
  TPObject *tp = new TPObject[nmig];

  CkPrintf("[work2] ready tp data structure\n");

  int j = 0;
  for(int i = 0; i < stats->n_objs; i++){
    int lbindex = tpCentroids[i].tag;
    if(!stats->migratable[lbindex]) continue;
    tp[j].centroid.x = tpCentroids[i].vec.x;
    tp[j].centroid.y = tpCentroids[i].vec.y;
    tp[j].centroid.z = tpCentroids[i].vec.z;
    tp[j].migratable = stats->migratable[lbindex];
    if(step() == 0){
      tp[j].load = tpCentroids[i].myNumParticles;
    }
    else if(useRatios){
      tp[j].load = ratios[i]*stats->wallTime[lbindex];
    }
    else{
      tp[j].load = stats->wallTime[lbindex];
    }
    tp[j].lbindex = lbindex;
    stats->tpindex[lbindex] = tpCentroids[i].tpindex;
    j++;
  }

  //printData(stats,phase);

  mapping = m;
  int dim = 0;
  TopoManager tmgr;

  procsPerNode = tmgr.getDimNT();

  int nx = tmgr.getDimNX();
  int ny = tmgr.getDimNY();
  int nz = tmgr.getDimNZ();
  int numnodes = nx*ny*nz; 

  CkPrintf("[work2] %d numnodes allocating %d bytes for nodes\n", numnodes, numnodes*sizeof(Node));
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

  CkPrintf("[work2] map\n");
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
