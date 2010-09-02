#include <charm++.h>
#include "cklists.h"
#include "Orb3dLB.h"
#include "ParallelGravity.h"
#include "TopoManager.h"
#include "Vector3D.h"

extern CProxy_TreePiece treeProxy;
extern CProxy_DataManager dMProxy;
extern CProxy_CkCacheManager cacheGravPart;
extern CProxy_CkCacheManager cacheGravSmooth;
extern CProxy_CkCacheManager cacheGravNode;

using namespace std;

CreateLBFunc_Def(Orb3dLB, "3d ORB mapping of tree piece space onto 3d processor mesh");

int comparx(const void *a, const void *b){
  TPObject *ta = (TPObject *)a;
  TPObject *tb = (TPObject *)b;
  return (int)(ta->centroid.x-tb->centroid.x);
}
int compary(const void *a, const void *b){
  TPObject *ta = (TPObject *)a;
  TPObject *tb = (TPObject *)b;
  return (int)(ta->centroid.y-tb->centroid.y);
}
int comparz(const void *a, const void *b){
  TPObject *ta = (TPObject *)a;
  TPObject *tb = (TPObject *)b;
  return (int)(ta->centroid.z-tb->centroid.z);
}

int pcx(const void *a, const void *b){
  Node *ta = (Node *)a;
  Node *tb = (Node *)b;
  return (int)(ta->x-tb->x);
}
int pcy(const void *a, const void *b){
  Node *ta = (Node *)a;
  Node *tb = (Node *)b;
  return (int)(ta->y-tb->y);
}
int pcz(const void *a, const void *b){
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

  //thisgroup[CkMyPe()].waitForGraphParts();
}

void Orb3dLB::receiveCentroids(CkReductionMsg *msg){
  int i = 0;
   TaggedVector3D * cur = (TaggedVector3D *)msg->getData();
  CkPrintf("Orb3dLB: receiveCentroids started: %d elements, msg length: %d\n", msg->getGcount(), msg->getLength()); 
  tpCentroids.free();
  
  while(i < msg->getGcount()){
     tpCentroids.push_back(*cur);
     cur = cur + 1;
     i++;
  }
  treeProxy.doAtSync();
  CkPrintf("Orb3dLB: receiveCentroids done\n");  
  delete msg;
}

//jetley
CmiBool Orb3dLB::QueryBalanceNow(int step){
  if(step == 0){
    if(CkMyPe() == 0){                          // only one group member need broadcast
      CkPrintf("Orb3dLB: Step 0, calling treeProxy.receiveProxy(thisgroup)\n");
      treeProxy.receiveProxy(thisgroup);        // broadcast proxy to all treepieces
    }
    firstRound = true;
    return false; 
  }
  if(CkMyPe() == 0)
    CkPrintf("Orb3dLB: Step %d\n", step);
  return true;

}

void Orb3dLB::work(BaseLB::LDStats* stats, int count)
{
  int numobjs = stats->n_objs;
  TPObject *tp = new TPObject[numobjs];

  stats->makeCommHash();
  CkAssert(tpCentroids.length() == stats->n_objs);
  for(int i = 0; i < stats->n_objs; i++){
    LDObjHandle &handle = tpCentroids[i].handle;
    int tag = stats->getHash(handle.id,handle.omhandle.id);
    tp[tag].centroid.x = tpCentroids[i].vec.x;
    tp[tag].centroid.y = tpCentroids[i].vec.y;
    tp[tag].centroid.z = tpCentroids[i].vec.z;
    tp[tag].load = stats->objData[i].wallTime;
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

  map(tp,numobjs,numnodes,nodes,nx,ny,nz,dim);

  float *procload = new float[stats->count];
  for(int i = 0; i < stats->count; i++){
    procload[i] = 0.0;
  }

  for(int i = 0; i < numobjs; i++){
    //CkPrintf("obj %d to %d (%f %f %f)\n", stats->objData[i].id().id[0], stats->to_proc[i], tp[i].centroid.x, tp[i].centroid.y, tp[i].centroid.z);
    procload[stats->to_proc[i]] += tp[i].load;
  }

  for(int i = 0; i < stats->count; i++){
    CkPrintf("proc %d load %f\n", i, procload[i]);
  }

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

void Orb3dLB::directMap(TPObject *tp, int ntp, Node *nodes){
  //CkPrintf("[Orb3dLB] mapping %d objects to Node (%d,%d,%d)\n", ntp, nodes[0].x, nodes[0].y, nodes[0].z);

  for(int i = 0; i < ntp; i++){
    CkPrintf("obj %d %f %f %f %f to node %d %d %d\n", tp[i].lbindex, tp[i].load, tp[i].centroid.x, tp[i].centroid.y, tp[i].centroid.z, nodes[0].x, nodes[0].y, nodes[0].z);
  }
  
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

  while(!pq_obj.empty()){
    TPObject tp = pq_obj.top();
    pq_obj.pop();

    Processor p = pq_proc.top();
    pq_proc.pop();

    p.load += tp.load;
    (*mapping)[tp.lbindex] = nodes[0].procRanks[p.t];

    pq_proc.push(p);
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

TPObject *Orb3dLB::partitionEvenLoad(TPObject *tp, int &ntp){
  float partition1Load = 0.0;
  float totalLoad = 0.0;
  for(int i = 0; i < ntp; i++){
    totalLoad += tp[i].load;
  }
  float halfLoad = 0.5*totalLoad;
  int split = -1;

  for(int i = 0; i < ntp; i++){
    // if including this element in partition1 brings us closer
    // to halfLoad, do it
    if((partition1Load+tp[i].load-halfLoad) < (halfLoad-partition1Load)){
      partition1Load += tp[i].load;
      split++;
    }
  }

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
