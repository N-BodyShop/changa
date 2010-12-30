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

//**************************************
// ORB3DLB functions
//**************************************
static int comparx(const void *a, const void *b){
  TPObject *ta = (TPObject *)a;
  TPObject *tb = (TPObject *)b;
  return (int)(ta->centroid.x-tb->centroid.x);
}
static int compary(const void *a, const void *b){
  TPObject *ta = (TPObject *)a;
  TPObject *tb = (TPObject *)b;
  return (int)(ta->centroid.y-tb->centroid.y);
}
static int comparz(const void *a, const void *b){
  TPObject *ta = (TPObject *)a;
  TPObject *tb = (TPObject *)b;
  return (int)(ta->centroid.z-tb->centroid.z);
}

static int pcx(const void *a, const void *b){
  Node *ta = (Node *)a;
  Node *tb = (Node *)b;
  return (int)(ta->x-tb->x);
}
static int pcy(const void *a, const void *b){
  Node *ta = (Node *)a;
  Node *tb = (Node *)b;
  return (int)(ta->y-tb->y);
}
static int pcz(const void *a, const void *b){
  Node *ta = (Node *)a;
  Node *tb = (Node *)b;
  return (int)(ta->z-tb->z);
}
//**************************************


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
    haveTPCentroids = false;
    delete tpmsg; 
  }

  tpCentroids = (TaggedVector3D *)msg->getData();
  nrecvd = msg->getGcount();
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
  if(step == 0){
    if(CkMyPe() == 0){                          // only one group member need broadcast
      CkPrintf("MultistepLB: Step 0, calling treeProxy.receiveProxy(thisgroup)\n");
      treeProxy.receiveProxy(thisgroup);        // broadcast proxy to all treepieces
    }
    firstRound = true;
    return false; 
  }
  if(CkMyPe() == 0)
    CkPrintf("MultistepLB: Step %d\n", step);
  return true;

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

  int i, len;
  int whichPos;
  int numAdditional;

  // tune alpha as needed - this is the merge parameter
  double alpha = 0.0;
  double savedCpu, savedWall;
  
  if(phase == -1){
#ifdef MCLBMSV
    CkPrintf("phase = -1, discarding\n");
#endif
    return;
  }
  
  len = savedPhaseStats.length();
  
  if(phase > len-1){
    numAdditional = phase-len+1;
    while(numAdditional > 0){
      savedPhaseStats.push_back(BaseLB::LDStats());
#ifdef MCLBMSV
      CkPrintf("Making new entry for phase %d (%d)\n", savedPhaseStats.length()-1, phase);
#endif
      savedPhaseStats[savedPhaseStats.length()-1].n_objs = 0;  // so that we can identify this phase's entry as being unfilled when we get the data for it
      numAdditional--;
    }
    len = savedPhaseStats.length();
    for(i = 0; i < stats->n_objs; i++)
      savedPhaseStats[len-1].objData.push_back(stats->objData[i]);
    whichPos = len-1;
  }
  else{ 
    whichPos = phase;
    if(savedPhaseStats[phase].n_objs == 0){       // haven't yet populated this phase
#ifdef MCLBMSV
      CkPrintf("Found unpopulated entry for phase %d\n", phase);
#endif
      for(i = 0; i < stats->n_objs; i++)
        savedPhaseStats[phase].objData.push_back(stats->objData[i]);
    }
    else{        // filled this phase out some time in the past - merge current with previous data
#ifdef MCLBMSV
      CkPrintf("Found previous entry for phase %d - merging\n", phase);
#endif
      for(i = 0; i < stats->n_objs; i++){
        savedCpu =  savedPhaseStats[phase].objData[i].cpuTime;
        savedWall =  savedPhaseStats[phase].objData[i].wallTime;
        
        savedPhaseStats[phase].objData[i]= stats->objData[i];
        
        savedPhaseStats[phase].objData[i].cpuTime = alpha*savedCpu + (1.0-alpha)*stats->objData[i].cpuTime;
        savedPhaseStats[phase].objData[i].wallTime = alpha*savedWall + (1.0-alpha)*stats->objData[i].wallTime;
        
      }
    }
  }
  savedPhaseStats[whichPos].n_objs=  stats->n_objs;
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
     CkPrintf("%d: %g %g\n", i, stats.objData[i].cpuTime,
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
  int phase = determinePhase(tpCentroids[0].activeRung);
  int prevPhase = tpCentroids[0].prevActiveRung;
  float *ratios = new float[stats->n_objs];
  // save pointers to centroids of treepieces
  Vector3D<float> **pCentroids = new Vector3D<float> *[stats->n_objs];

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
  
  for(i = 0; i < stats->n_objs; i++){
    ratios[tpCentroids[i].tag] = tpCentroids[i].numActiveParticles/(float)tpCentroids[i].myNumParticles;
    numActiveParticles += tpCentroids[i].numActiveParticles;
    totalNumParticles += tpCentroids[i].myNumParticles;
    pCentroids[i] = &tpCentroids[i].vec;

    if(tpCentroids[i].numActiveParticles == 0){
      numInactiveObjects++;
      if(stats->objData[tpCentroids[i].tag].migratable){
        stats->objData[tpCentroids[i].tag].migratable = 0;
#ifdef MCLBMSV
        CkPrintf("marking object %d non-migratable (inactive)\n", tpCentroids[i].tag);
#endif
        stats->n_migrateobjs--;
      }
    }
    else{
      numActiveObjects++;
      //CkPrintf("object %d (proc %d) active, ratio: %f\n", map.tpCentroids[i].tag, 
                                              //stats->from_proc[map.tpCentroids[i].tag], ratios[map.tpCentroids[i].tag]);
    }
  }
#ifdef MCLBMSV
  CkPrintf("numActiveObjects: %d, numInactiveObjects: %d\n", numActiveObjects, numInactiveObjects);
#endif
  if(havePhaseData(phase)){
#ifdef MCLBMSV
    CkPrintf("phase %d data available\n", phase);
#endif
    for(i = 0; i < stats->n_objs; i++){
      stats->objData[i].cpuTime = savedPhaseStats[phase].objData[i].cpuTime;
      stats->objData[i].wallTime = savedPhaseStats[phase].objData[i].wallTime;
    }
  }
  else if(havePhaseData(0)){
#ifdef MCLBMSV
    CkPrintf("phase %d data unavailable, using phase 0 loads\n", phase);
#endif
    //CkPrintf("using phase 0 loads\n", phase);
    for(i = 0; i < stats->n_objs; i++){
      stats->objData[i].cpuTime = ratios[i]*savedPhaseStats[0].objData[i].cpuTime;
      stats->objData[i].wallTime = ratios[i]*savedPhaseStats[0].objData[i].wallTime;
    }
  }
  else{
#ifdef MCLBMSV
    CkPrintf("phase %d data unavailable\n", phase);
#endif
    delete[] pCentroids;
    delete[] ratios;
    return;
  }
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
#ifdef NOTDF
  statsData = stats;

  P = count;

  // calculate total number of migratable objects
  nObjs = stats->n_migrateobjs;
#ifdef MCLBMSV
  CkPrintf("OrbLB: num objects: %d\n", nObjs);
#endif

  // create computeLoad and calculate tentative computes coordinates
  computeLoad = new ComputeLoad[nObjs];
  for (i=XDIR; i<=ZDIR; i++) vArray[i] = new VecArray[nObjs];

  // v[0] = XDIR  v[1] = YDIR v[2] = ZDIR
  // vArray[XDIR] is an array holding the x vector for all computes
  int objIdx = 0;
  for (i=0; i<stats->n_objs; i++) {
    LDObjData &odata = stats->objData[i];
    if (odata.migratable == 0) continue;
#ifdef MCLBMSV
    CkPrintf("OrbLB: considering object %d\n", i);
#endif
    computeLoad[objIdx].id = i;
    /*
    computeLoad[objIdx].v[XDIR] = odata.objID().id[0];
    computeLoad[objIdx].v[YDIR] = odata.objID().id[1];
    computeLoad[objIdx].v[ZDIR] = odata.objID().id[2];
    */
    Vector3D<float> *pvec = pCentroids[i];
    computeLoad[objIdx].v[XDIR] = pvec->x;
    computeLoad[objIdx].v[YDIR] = pvec->y;
    computeLoad[objIdx].v[ZDIR] = pvec->z;
    computeLoad[objIdx].load = _lb_args.useCpuTime()?odata.cpuTime:odata.wallTime;
    computeLoad[objIdx].refno = 0;
    computeLoad[objIdx].partition = NULL;
    for (int k=XDIR; k<=ZDIR; k++) {
        vArray[k][objIdx].id = objIdx;
        vArray[k][objIdx].v = computeLoad[objIdx].v[k];
    }
#ifdef DEBUG
    CmiPrintf("Object %d: %d %d %d load:%f\n", objIdx, computeLoad[objIdx].v[XDIR], computeLoad[objIdx].v[YDIR], computeLoad[objIdx].v[ZDIR], computeLoad[objIdx].load);
#endif
    objIdx ++;
  }
  CmiAssert(nObjs == objIdx);

  double t = CkWallTimer();

  quicksort(XDIR);
  quicksort(YDIR);
  quicksort(ZDIR);
#ifdef DEBUG
  CmiPrintf("qsort1 time: %f\n", CkWallTimer() - t);
#endif

  npartition = P;
  partitions = new Partition[npartition];

  double totalLoad = 0.0;
  int minx, miny, minz, maxx, maxy, maxz;
  minx = maxx= computeLoad[0].v[XDIR];
  miny = maxy= computeLoad[0].v[YDIR];
  minz = maxz= computeLoad[0].v[ZDIR];
  for (i=1; i<nObjs; i++) {
    totalLoad += computeLoad[i].load;
    if (computeLoad[i].v[XDIR] < minx) minx = computeLoad[i].v[XDIR];
    else if (computeLoad[i].v[XDIR] > maxx) maxx = computeLoad[i].v[XDIR];
    if (computeLoad[i].v[YDIR] < miny) miny = computeLoad[i].v[YDIR];
    else if (computeLoad[i].v[YDIR] > maxy) maxy = computeLoad[i].v[YDIR];
    if (computeLoad[i].v[ZDIR] < minz) minz = computeLoad[i].v[ZDIR];
    else if (computeLoad[i].v[ZDIR] > maxz) maxz = computeLoad[i].v[ZDIR];
  }

  top_partition.origin[XDIR] = minx;
  top_partition.origin[YDIR] = miny;
  top_partition.origin[ZDIR] = minz;
  top_partition.corner[XDIR] = maxx;
  top_partition.corner[YDIR] = maxy; 
  top_partition.corner[ZDIR] = maxz;

  top_partition.refno = 0;
  top_partition.load = 0.0;
  top_partition.count = nObjs;

  // if we take background load into account
  if (!_lb_args.ignoreBgLoad()) {
    top_partition.bkpes.resize(0);
    double total = totalLoad;
    for (i=0; i<P; i++) {
      double bkload = stats->procs[i].bg_walltime;
      total += bkload;
    }
    double averageLoad = total / P;
    for (i=0; i<P; i++) {
      double bkload = stats->procs[i].bg_walltime;
      if (bkload < averageLoad) top_partition.bkpes.push_back(i);
      else CkPrintf("MultistepLB Info> PE %d with %f background load will have 0 object.\n", i, bkload);
    }
    npartition = top_partition.bkpes.size();
    // formally add these bg load to total load
    for (i=0; i<npartition; i++) 
      totalLoad += stats->procs[top_partition.bkpes[i]].bg_walltime; 
    if (_lb_args.debug()>=2) {
      CkPrintf("BG load: ");
      for (i=0; i<P; i++)  CkPrintf(" %f", stats->procs[i].bg_walltime);
      CkPrintf("\n");
      CkPrintf("Partition BG load: ");
      for (i=0; i<npartition; i++)  CkPrintf(" %f", stats->procs[top_partition.bkpes[i]].bg_walltime);
      CkPrintf("\n");
    }
  }

  top_partition.load = totalLoad;

  currentp = 0;
  refno = 0;

  // recursively divide
  rec_divide(npartition, top_partition);

  // mapping partitions to nodes
  mapPartitionsToNodes();

  // this is for sanity check
  int *num = new int[P];
  for (i=0; i<P; i++) num[i] = 0;

  for (i=0; i<nObjs; i++)
  {
    for (j=0; j<npartition; j++)
      if (computeLoad[i].refno == partitions[j].refno)   {
        computeLoad[i].partition = partitions+j;
        num[j] ++;
    }
    CmiAssert(computeLoad[i].partition != NULL);
  }

  for (i=0; i<npartition; i++)
    if (num[i] != partitions[i].count) 
      CmiAbort("MultistepLB: Compute counts don't agree!\n");

  delete [] num;

  // Save output
  for(int obj = 0; obj < stats->n_migrateobjs; obj++){
    int frompe = stats->from_proc[computeLoad[obj].id];
    int tope = computeLoad[obj].partition->node;
    if(frompe != tope){
      stats->to_proc[computeLoad[obj].id] = tope;
    }
    CkPrintf("%d(%d): %d -> %d\n", computeLoad[obj].id, obj, frompe, tope);
  }
  /*
  objIdx = 0;
  for(int obj=0;obj<stats->n_objs;obj++) {
      stats->to_proc[obj] = stats->from_proc[obj];
      LDObjData &odata = stats->objData[obj];
      if (odata.migratable == 0) { continue; }
      int frompe = stats->from_proc[obj];
      int tope = computeLoad[objIdx].partition->node;
      if (frompe != tope) {
        if (_lb_args.debug() >= 3) {
              CkPrintf("[%d] Obj %d migrating from %d to %d\n",
                     CkMyPe(),obj,frompe,tope);
        }
	stats->to_proc[obj] = tope;
      }
      objIdx ++;
  }
*/
#define min(a,b) (a <= b ? a : b)
#define ceil(a,b) (a%b == 0? a/b : a/b+1)
  // jetley - smooth out the load profile
  int stride = 4, step = 4;             // if step = stride, groups of processors are disjoint and productive
  int k;
  int group;
  double th = 0.15;
  
  double *procWeights = new double[count];
  CkVec<WeightObject> *objectList = new CkVec<WeightObject>[count];        // one for each processor
  Vector3D<float> *procCentroids = new Vector3D<float>[count];             // one for each processor - COM of all objects on processor
  Vector3D<float> *objCentroids = new Vector3D<float>[stats->n_objs];      // one for each object
  
  // initialize COM's to zero
  for(i = 0; i < count; i++){
    procCentroids[i].x = 0;
    procCentroids[i].y = 0;
    procCentroids[i].z = 0;
    procWeights[i] = 0.0;
  }

  // set up centroids of individual objects
  for(i = 0; i < stats->n_objs; i++)
    objCentroids[tpCentroids[i].tag] = tpCentroids[i].vec;      // now objCentroids[i] contains centroid of object i
    
  // populate objectList, calc. COM's and weights 
  double tempWeight;
  for(i = 0; i < stats->n_objs; i++){
    LDObjData &odata = stats->objData[i];
    if(odata.migratable == 0) continue;
    tempWeight = odata.wallTime*stats->procs[stats->to_proc[i]].pe_speed;
    (objectList[stats->to_proc[i]]).push_back(WeightObject(i,tempWeight));      // membership
    procCentroids[stats->to_proc[i]] += tempWeight*objCentroids[i];             // COM
    procWeights[stats->to_proc[i]] += tempWeight;                               // weight
  }
  for(i = 0; i < count; i++)
    procCentroids[i] /= procWeights[i];
    
  // begin balancing
  // for each group g
#ifdef MCLBMSV
  CkPrintf("MCLB: Start balancing\n");
#endif
  for(i = 0; i+stride <= count; i += step){
    CkVec<int> heavy, light;
    double idealAvg = 0.0; 
    bool done; 
#ifdef MCLBMSV
    CkPrintf("MCLB: balance [%d,%d]\n", i, i+stride-1);
#endif
    // calc. idealAvg for this group of processors
    for(j = i; j < i + stride; j++)
      for(k = 0; k < (objectList[j]).size(); k++){
        CkAssert((objectList[j])[k].weight == stats->objData[(objectList[j])[k].idx].wallTime*stats->procs[j].pe_speed);
        idealAvg += (objectList[j])[k].weight;
      }
    idealAvg /= stride;
#ifdef MCLBMSV
    CkPrintf("MCLB: idealAvg: %f\n", idealAvg);
#endif
    
    // compute heavy and light sets - H(g), L(g)
    for(j = i; j < i + stride; j++){
      if(procWeights[j] >= (1+th)*idealAvg){
        heavy.push_back(j);
        // CkPrintf("MCLB: %d(%f) heavy\n", j, procWeights[j]);
      }
      else if(procWeights[j] < (1-th)*idealAvg){
        light.push_back(j);
        // CkPrintf("MCLB: %d(%f) light\n", j, procWeights[j]);
      }
    }
    // for each heavy processor p
    for(j = 0; j < heavy.size(); j++){
      int donor = heavy[j];
      done = false;
      // CkPrintf("MCLB: Picked heavy proc. %d(%f)\n", donor, procWeights[donor]);
      // pick an o in O(p)
      // CkPrintf("MCLB: Pick an object\n");
      for(int r = 0; r < objectList[donor].size(); r++){ 
        int donatedObj = (objectList[donor])[r].idx;
        double donatedWeight = (objectList[donor])[r].weight;
        // CkPrintf("MCLB: Trying object %d(%f)\n", donatedObj, donatedWeight);
        // pick a p' in L(g) such that a swap is viable (weight) and advisable (distance)
        int acceptor = -1;
        float closestDist = 3.4e+38;    // about the largest fp value possible 
        int savedLightIndex = -1;
        for(k = 0; k < light.size(); k++){
          if(donatedWeight+procWeights[light[k]] <= (1+th)*idealAvg && 
            procWeights[donor]-donatedWeight >= (1-th)*idealAvg &&
            closestDist > (procCentroids[light[k]]-objCentroids[donatedObj]).length()){
            closestDist = (procCentroids[light[k]]-objCentroids[donatedObj]).length();
            acceptor = light[k];
            savedLightIndex = k;
          }
        }
        // found an acceptor p' for object o
        if(acceptor > 0){
          // CkPrintf("MCLB: Picked light proc. %d(%f)\n", acceptor, procWeights[acceptor]);
          procCentroids[donor] *= procWeights[donor];
          procCentroids[donor] -= objCentroids[donatedObj]*donatedWeight;
          procWeights[donor] -= donatedWeight;
          procCentroids[donor] /= procWeights[donor];
          
          procCentroids[acceptor] *= procWeights[acceptor];
          procCentroids[acceptor] += objCentroids[donatedObj]*donatedWeight;
          procWeights[acceptor] += donatedWeight;
          procCentroids[acceptor] /= procWeights[acceptor];
          
          objectList[donor].remove(r); // expensive - FIXME
          objectList[acceptor].push_back(WeightObject(donatedObj, donatedWeight));
          //actually shift the donated object
          stats->to_proc[donatedObj] = acceptor;
          
          // if heavy processor isn't so anymore
          if(procWeights[donor] < (1+th)*idealAvg){
            // CkPrintf("MCLB: Proc. %d(%d) not heavy anymore\n", donor, j);
            heavy.remove(j);
            done = true;
          }
          if(procWeights[acceptor] >= (1-th)*idealAvg){
            // CkPrintf("MCLB: Proc. %d(%d) not light anymore\n", acceptor, savedLightIndex);
            light.remove(savedLightIndex);  
          }
          if(done)
            break;
        }// end if acceptor found
      }// end pick an object
      //if(!done)
        //CkPrintf("MultistepLB: Warning - couldn't flatten processor %d\n", donor);
    }// end for each heavy processor
  }// end for each group
  
  // free memory
  delete [] procCentroids;
  delete [] objectList;
  delete [] objCentroids;
  delete [] procWeights;
  
  delete [] computeLoad;
  for (i=0; i<3; i++) delete [] vArray[i];
  delete [] partitions;

  if (_lb_args.debug() >= 1)
    CkPrintf("MultistepLB finished time: %fs\n", CkWallTimer() - t);
#else // Orb3dLB
  CkPrintf("******** BIG STEP *********!\n");
  work2(stats,count);
#endif  // MCLBMS_ORBSMOOTH
    
  }     // end if phase == 0
  else{ // not even greedy; round-robin
    for(i = 0; i < stats->n_objs; i++){
      if(stats->objData[i].migratable == 0) continue;
      stats->to_proc[i] = i%stats->count;
    }
  }
  delete[] pCentroids;
#endif //CMK_LDB_ON

}

//**************************************
// ORB3DLB functions
//**************************************

void MultistepLB::work2(BaseLB::LDStats *stats, int count){
  int numobjs = stats->n_objs;
  int nmig = stats->n_migrateobjs;

  CkPrintf("[work2] %d objects allocating %d bytes for tp\n", nmig, nmig*sizeof(TPObject));
  TPObject *tp = new TPObject[nmig];

  CkPrintf("[work2] ready tp data structure\n");

  int j = 0;
  for(int i = 0; i < stats->n_objs; i++){
    int tag = tpCentroids[i].tag;
    if(!stats->objData[tag].migratable) continue;
    tp[j].centroid.x = tpCentroids[i].vec.x;
    tp[j].centroid.y = tpCentroids[i].vec.y;
    tp[j].centroid.z = tpCentroids[i].vec.z;
    tp[j].migratable = stats->objData[tag].migratable;
    if(step() == 0){
      tp[j].load = tpCentroids[i].myNumParticles;
    }
    else{
      tp[j].load = stats->objData[tag].wallTime;
    }
    tp[j].lbindex = tag;
    j++;
  }

  mapping = &stats->to_proc;
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
