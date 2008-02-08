/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

/**
 * \addtogroup CkLdb
*/
/*@{*/

#include <charm++.h>
#include "cklists.h"
#include "MetisCosmoLB.h"
#include "ParallelGravity.h"

#include "clustering.cpp"

extern CProxy_TreePiece treeProxy;

CreateLBFunc_Def(MetisCosmoLB, "Use Metis(tm) to partition object graph; assignment of objects based on scaling heuristic");

MetisCosmoLB::MetisCosmoLB(const CkLBOptions &opt): CentralLB(opt)
{
  lbname = "MetisCosmoLB";
  centroidsAllocated = false;
  if (CkMyPe() == 0)
    CkPrintf("[%d] MetisCosmoLB created\n",CkMyPe());
}

void MetisCosmoLB::receiveCentroids(CkReductionMsg *msg){
  int i = 0;
  // CkReduction::setElement *cur = (CkReduction::setElement *)msg->getData();
   TaggedVector3D * cur = (TaggedVector3D *)msg->getData();
  //CmiMemoryCheck();
  CkPrintf("MetisCosmoLB: receiveCentroids started: %d elements, msg length: %d\n", msg->getGcount(), msg->getLength()); 
  map.tpCentroids.free();
  
  while(i < msg->getGcount()){
    /*
     TaggedVector3D *tv = (TaggedVector3D *)cur->data; 
     CkAssert(tv != NULL);
     */
     //map.tpCentroids.push_back(*tv);
     map.tpCentroids.push_back(*cur);
     // cur = cur->next();
     cur = cur + 1;
     //CkPrintf("%d done\n", i);
     i++;
  }
  //CkPrintf("MetisCosmoLB: loop done\n");  
  
  /*
  if(!centroidsAllocated){
    CkPrintf("centroids not allocated\n");
    centroidsAllocated = true;
    while(cur != NULL){
      TaggedVector3D *tv = (TaggedVector3D *)cur->data; 
      map.tpCentroids.push_back(*tv);
      cur = cur->next();
    }
  }
  else{
    i = 0;
    CkPrintf("centroids allocated\n");
    while(cur != NULL){
      TaggedVector3D *tv = (TaggedVector3D *)cur->data; 
      map.tpCentroids[i] = *tv;
      cur = cur->next();
      i++;
    }
  }
  */
  treeProxy.doAtSync();
  //CmiMemoryCheck();
  CkPrintf("MetisCosmoLB: receiveCentroids done\n");  
  delete msg;
}

//jetley
CmiBool MetisCosmoLB::QueryBalanceNow(int step){
  /*
  if(step == 0){
#if COSMO_MCLB > 1
    CkPrintf("MetisCosmoLB: Step 0, calling treeProxy.receiveProxy(thisgroup)\n");
#endif
    if(CkMyPe() == 0)                          // only one group member need broadcast
      treeProxy.receiveProxy(thisgroup);        // broadcast proxy to all treepieces
    return false;
  }
#if COSMO_MCLB > 1
  CkPrintf("MetisCosmoLB: QueryBalanceNow returning true\n");
#endif
  return true;
  */
  if(step == 0){
    if(CkMyPe() == 0)                          // only one group member need broadcast
      CkPrintf("MetisCosmoLB: Step 0, calling treeProxy.receiveProxy(thisgroup)\n");
      treeProxy.receiveProxy(thisgroup);        // broadcast proxy to all treepieces
    firstRound = true;
    return false; 
  }
//#if COSMO_MCLB > 1
  if(CkMyPe() == 0)
    CkPrintf("MetisCosmoLB: Step %d\n", step);
//#endif
  return true;

}

void MetisCosmoLB::rec_divide(int n, Partition &p)
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

void MetisCosmoLB::setVal(int x, int y, int z)
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

int MetisCosmoLB::sort_partition(int x, int p, int r)
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

void MetisCosmoLB::qsort(int x, int p, int r)
{
  if (p<r) {
    int q = sort_partition(x, p, r);
//CmiPrintf("midpoint: %d %d %d\n", p,q,r);
    qsort(x, p, q-1);
    qsort(x, q+1, r);
  }
}

void MetisCosmoLB::quicksort(int x)
{
  int y = (x+1)%3;
  int z = (x+2)%3;
  setVal(x, y, z);
  qsort(x, 0, nObjs-1);

#if 0
  CmiPrintf("result for :%d\n", x);
  for (int i=0; i<nObjs; i++) 
    CmiPrintf("%d ", computeLoad[vArray[x][i].id].tv);
  CmiPrintf("\n");
#endif
}

void MetisCosmoLB::mapPartitionsToNodes()
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
unsigned int MetisCosmoLB::determinePhase(unsigned int lastActiveRung){
  return lastActiveRung;
}

// merge data instrumented in previous iteration of computation with data from earlier iterations
void MetisCosmoLB::mergeInstrumentedData(int phase, BaseLB::LDStats *stats){

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
  printData(savedPhaseStats[whichPos], phase, NULL);
#endif
}

// whether we have instrumented data for this phase
bool MetisCosmoLB::havePhaseData(int phase){
  return (savedPhaseStats.length() > phase && savedPhaseStats[phase].n_objs > 0);
}

void MetisCosmoLB::printData(BaseLB::LDStats &stats, int phase, int *revObjMap){
  int i,j;
  
  CkPrintf("---- data (%d): %d objects ----\n", phase, stats.n_objs);
  for(i = 0; i < stats.n_objs; i++){
    CkPrintf("%d: %f\n", i, stats.objData[i].cpuTime);
  }
  CkPrintf("---- end data (%d) ----\n", phase);
}

void MetisCosmoLB::makeActiveProcessorList(BaseLB::LDStats *stats, int numActiveObjs){
  int i;
  int objsPerProc = 8;
  int expandFactor = 4;
  int procsNeeded;
  procsNeeded = expandFactor*numActiveObjs/objsPerProc > stats->count ? stats->count : expandFactor*numActiveObjs/objsPerProc;

  /* currently, only the first procsNeeded procs are used - could do something more sophisticated here in the future - FIXME */
  /*
  newStats->procs = new BaseLB::ProcStats [procsNeeded];
  for(i = 0; i < procsNeeded; i++){
    newStats->procs[i] = stats->procs[i];
  }
  */
  stats->count = procsNeeded;
#ifdef MCLBMSV
  CkPrintf("Processors 0 to %d active\n", procsNeeded-1);
#endif
}
#endif

void MetisCosmoLB::work(BaseLB::LDStats* stats, int count)
{
  if (_lb_args.debug() >= 2) {
    CkPrintf("[%d] In MetisCosmoLB Strategy...\n", CkMyPe());
  }
  int i, j, m;
  int option = 0;

  // jetley - statistics
  int communicating = 0;
  int wellPlaced = 0;

  CkPrintf("In work()...");
  stats->makeCommHash();
  for(i = 0; i < stats->n_objs; i++){
    LDObjHandle &handle = map.tpCentroids[i].handle;
    map.tpCentroids[i].tag = stats->getHash(handle.id, handle.omhandle.id);
    CkPrintf("tpCentroids[%d].tag = %d\n", i, map.tpCentroids[i].tag);
  }

  removeNonMigratable(stats, count);

  int numobjs = stats->n_objs;

  // allocate space for the computing data
  double *objtime = new double[numobjs];
  int *origmap = new int[numobjs];
  LDObjHandle *handles = new LDObjHandle[numobjs];
  for(i=0;i<numobjs;i++) {
    objtime[i] = 0.0;
    origmap[i] = 0;
  }

  for (i=0; i<stats->n_objs; i++) {
      LDObjData &odata = stats->objData[i];
      if (!odata.migratable) 
        CmiAbort("MetisCosmoLB does not support nonmigratable objects.\n");
      /*
      origmap[odata[i].id.id[0]] = j;
      cputime[odata[i].id.id[0]] = odata[i].cpuTime;
      handles[odata[i].id.id[0]] = odata[i].handle;
      */
      int frompe = stats->from_proc[i];
      origmap[i] = frompe;
      objtime[i] = odata.wallTime*stats->procs[frompe].pe_speed;
      handles[i] = odata.handle;
  }
  CkPrintf("initialized, ");

  // to convert the weights on vertices to integers
  
  int *newmap;
  int sameMapFlag = 1;

  if (count < 1) {
    CkPrintf("error: Number of Pe less than 1!");
  }
  else if (count == 1) {
    newmap = origmap;
    sameMapFlag = 1;
  }
  else {
    sameMapFlag = 0;
    newmap = new int[numobjs];
	// call clustering.cpp function here
    CkPrintf("calling clustering\n");
    clustering(numobjs, objtime, map.tpCentroids, count, newmap);
    CkPrintf("clustering done\n");
    
  }

	
  /* jetley - at this point, 'newmap' gives us clusters of objects. each cluster goes to some processor.
   we now need to map these clusters to particular processors. base this final mapping on some heuristic.
   for cosmology on BG, 'nearby' clusters are mapped to physically proximal processors.
  */
  
  int *clusterWeights = NULL; // done on purpose, no need for clusterWeights
  
  CkPrintf("calling assign\n");
  map.assign(newmap, clusterWeights, stats->to_proc, numobjs, count);
  if(_lb_args.debug() >= 3)
	  for(i=0; i<numobjs; i++)
		  CkPrintf("[%d] Obj %d migrating from %d to %d\n", CkMyPe(),i,stats->from_proc[i],stats->to_proc[i]);
  
  // efficacy check
/*
  for(i = 0; i < numobjs; i++)
    for(j = 0; j < numobjs; j++)
      if(comm[i][j] != 0 && map.isNeighborOf(stats->to_proc[i], stats->to_proc[j]))
        wellPlaced++;

  if(communicating != 0) // gives FP error otherwise; FIXME - why should communicating equal 0 ?
  CkPrintf("metis-cosmo-eff: Comm. graph is %f-complete, %f fraction of communicating objects well-placed\n", (float)(communicating/(numobjs*(numobjs-1))), (float)(wellPlaced/communicating));
*/

  delete[] objtime;
  delete[] origmap;
  if(newmap != origmap)
    delete[] newmap;
  if(handles) delete []handles;
  if (_lb_args.debug() >= 1) {
   CkPrintf("[%d] MetisCosmoLB done! \n", CkMyPe());
  }
}

#include "MetisCosmoLB.def.h"

/*@}*/
