#include <charm++.h>
#include "cklists.h"
#include "MultistepLB.h"
#include "TopoManager.h"
#include "ParallelGravity.h"
#include "Vector3D.h"
#include <queue>

extern CProxy_TreePiece treeProxy;
CkpvExtern(int, _lb_obj_index);
using namespace std;

CreateLBFunc_Def(MultistepLB, "Works best with multistepped runs; uses Orb3D for larger steps, greedy otherwise");

void MultistepLB::init() {
  lbname = "MultistepLB";

  if (CkpvAccess(_lb_obj_index) == -1)
    CkpvAccess(_lb_obj_index) = LBRegisterObjUserData(sizeof(TaggedVector3D));

  compares[0] = comparx;
  compares[1] = compary;
  compares[2] = comparz;

  pc[0] = pcx;
  pc[1] = pcy;
  pc[2] = pcz;
}

MultistepLB::MultistepLB(const CkLBOptions &opt): CBase_MultistepLB(opt)
{
  init();

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

  int numActiveObjects = 0;
  int numInactiveObjects = 0;

  // to calculate ratio of active particles in phase
  int64_t numActiveParticles = 0;
  int64_t totalNumParticles = 0;

  for(int i = 0; i < stats->n_objs; i++){
    stats->to_proc[i] = stats->from_proc[i];
  }

  for(int i = 0; i < stats->n_objs; i++){
    if (!stats->objData[i].migratable) continue;

    LDObjData &odata = stats->objData[i];
    TaggedVector3D* udata = (TaggedVector3D *)odata.getUserData(CkpvAccess(_lb_obj_index));

    numActiveParticles += udata->numActiveParticles;
    totalNumParticles += udata->myNumParticles;

    if(udata->numActiveParticles == 0){
      numInactiveObjects++;
      if(stats->objData[i].migratable){
        stats->objData[i].migratable = 0;
#ifdef MCLBMSV
        CkPrintf("marking object %d non-migratable (inactive)\n", i);
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

  // let the strategy take over on this modified instrumented data and processor information
  if((float)numActiveParticles/totalNumParticles > LARGE_PHASE_THRESHOLD){
  //if(true){
    if (_lb_args.debug()>=2) {
      CkPrintf("******** BIG STEP *********!\n");
    }
    work2(stats,count);
  }     // end if phase == 0
  else{
    greedy(stats,count);
  }
#endif //CMK_LDB_ON

}

//**************************************
// ORB3DLB functions
//**************************************
//
void MultistepLB::greedy(BaseLB::LDStats *stats, int count){

  int numobjs = stats->n_objs;
  int nmig = stats->n_migrateobjs;
  CkPrintf("[GREEDY] objects total %d active %d\n", numobjs,nmig);

  TPObject *tp_array = new TPObject[nmig];
  int j = 0;
  for(int i = 0; i < stats->n_objs; i++){
    if(!stats->objData[i].migratable) continue;
    tp_array[j].migratable = stats->objData[i].migratable;

    LDObjData &odata = stats->objData[i];
    TaggedVector3D* udata = (TaggedVector3D *)odata.getUserData(CkpvAccess(_lb_obj_index));
    if(step() == 0){
      tp_array[j].load = udata->myNumParticles;
    }
    else{
      tp_array[j].load = stats->objData[i].wallTime;
    }
    tp_array[j].lbindex = i;
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
  CkPrintf("GREEDY MEASURED CPU LOAD\n");
  CkPrintf("**********************************\n");
  for(int i = 0; i < stats->count; i++){
    CkPrintf("[pestats] %d %g %g\n",
                               i,
                               stats->procs[i].total_walltime,
                               stats->procs[i].idletime);
  }



  delete []tp_array;
}

void MultistepLB::work2(BaseLB::LDStats *stats, int count){
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
    if(!stats->objData[i].migratable) continue;

    LDObjData &odata = stats->objData[i];
    TaggedVector3D* udata = (TaggedVector3D *)odata.getUserData(CkpvAccess(_lb_obj_index));

    tp_array[j].centroid = udata->vec;
    tp_array[j].migratable = true;
    if(step() == 0){
      tp_array[j].load = udata->myNumParticles;
    }
    else{
      tp_array[j].load = stats->objData[i].wallTime;
    }
    tp_array[j].lbindex = i;
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
  CkPrintf("MEASURED CPU LOAD\n");
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
  CBase_MultistepLB::pup(p);
  p | procsPerNode;
}

#include "MultistepLB.def.h"
