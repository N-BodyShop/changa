#include <charm++.h>
#include "cklists.h"
#include "Orb3dLB.h"
#include "ParallelGravity.h"
#include "TopoManager.h"
#include "Vector3D.h"

extern CProxy_TreePiece treeProxy;

using namespace std;

CreateLBFunc_Def(Orb3dLB, "3d ORB mapping of tree piece space onto 3d processor mesh");

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

Orb3dLB::Orb3dLB(const CkLBOptions &opt): CentralLB(opt)
{
  lbname = "Orb3dLB";
  centroidsAllocated = false;
  if (CkMyPe() == 0){
    CkPrintf("[%d] Orb3dLB created\n",CkMyPe());
    TopoManager tmgr;

    int ppn = tmgr.getDimNT();

    int nx = tmgr.getDimNX();
    int ny = tmgr.getDimNY();
    int nz = tmgr.getDimNZ();
    int numnodes = nx*ny*nz; 

    CkPrintf("[%d] Orb3dLB Topo %d %d %d %d %d \n",CkMyPe(), nx, ny, nz, numnodes, ppn);
  }
  compares[0] = comparx;
  compares[1] = compary;
  compares[2] = comparz;

  pc[0] = pcx;
  pc[1] = pcy;
  pc[2] = pcz;

  haveTPCentroids = false;

}

void Orb3dLB::receiveCentroids(CkReductionMsg *msg){
  if(haveTPCentroids){
    delete tpmsg;
  }
  tpCentroids = (TaggedVector3D *)msg->getData();
  nrecvd = msg->getGcount();
  tpmsg = msg;
  haveTPCentroids = true;
  // TaggedVector3D * cur = (TaggedVector3D *)msg->getData();
  CkPrintf("Orb3dLB: receiveCentroids started: %d elements, msg length: %d\n", msg->getGcount(), msg->getLength()); 
  //tpCentroids.free();
  
  /*
  while(i < msg->getGcount()){
     tpCentroids.push_back(*cur);
     cur = cur + 1;
     i++;
  }
  */
  treeProxy.doAtSync();
  CkPrintf("Orb3dLB: receiveCentroids done\n");  
  // delete msg later
  //delete msg;
}

//jetley
CmiBool Orb3dLB::QueryBalanceNow(int step){
  // Now, by the first step, we already have
  // centroids. of course, since no gravity has
  // been calculated yet, load can only be 
  // estimated via num particles in treepiece
  if(CkMyPe() == 0)
    CkPrintf("Orb3dLB: Step %d\n", step);
  //return true;
  if(step == 0){
    return false;
  }
  else{
    return true;
  }

}

void Orb3dLB::work(BaseLB::LDStats* stats, int count)
{
  int numobjs = stats->n_objs;
  int nmig = stats->n_migrateobjs;

  CkPrintf("[orb3dlb] %d objects allocating %d bytes for tp\n", numobjs, numobjs*sizeof(TPObject));
  TPObject *tp = new TPObject[numobjs];

  stats->makeCommHash();
  CkPrintf("[orb3dlb] ready tp data structure\n");
  if(nrecvd != numobjs){
    CkAbort("wrong tpCentroids length\n");
  }
  for(int i = 0; i < stats->n_objs; i++){
    LDObjHandle &handle = tpCentroids[i].handle;
    int tag = stats->getHash(handle.id,handle.omhandle.id);
    tp[tag].centroid.x = tpCentroids[i].vec.x;
    tp[tag].centroid.y = tpCentroids[i].vec.y;
    tp[tag].centroid.z = tpCentroids[i].vec.z;
    tp[tag].migratable = stats->objData[tag].migratable;
    if(step() == 0){
      tp[tag].load = tpCentroids[i].myNumParticles;
    }
    else{
      tp[tag].load = stats->objData[tag].wallTime;
    }
    tp[tag].lbindex = tag;
  }

  mapping = &stats->to_proc;
  int dim = 0;
  TopoManager tmgr;

  procsPerNode = tmgr.getDimNT();

  int nx = tmgr.getDimNX();
  int ny = tmgr.getDimNY();
  int nz = tmgr.getDimNZ();
  int numnodes = nx*ny*nz; 

  CkPrintf("[orb3dlb] %d numnodes allocating %d bytes for nodes\n", numnodes, numnodes*sizeof(Node));
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

  CkPrintf("[orb3dlb] map\n");
  map(tp,numobjs,numnodes,nodes,nx,ny,nz,dim);

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

void Orb3dLB::map(TPObject *tp, int ntp, int nn, Node *nodes, int xs, int ys, int zs, int dim){
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

void Orb3dLB::directMap(TPObject *tp, int ntp, Node *nodes){
  //CkPrintf("[Orb3dLB] mapping %d objects to Node (%d,%d,%d)\n", ntp, nodes[0].x, nodes[0].y, nodes[0].z);

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

int Orb3dLB::nextDim(int dim_, int xs, int ys, int zs){
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

TPObject *Orb3dLB::partitionEvenLoad(TPObject *tp, int &ntp){
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

Node *Orb3dLB::halveNodes(Node *start, int np){
  Node *ret = start;
  ret = start+np/2;
  return ret;
}

#include "Orb3dLB.def.h"

/*@}*/
