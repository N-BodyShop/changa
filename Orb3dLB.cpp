#include <charm++.h>
#include "cklists.h"
#include "Orb3dLB.h"
#include "ParallelGravity.h"
#include "TopoManager.h"
#include "Vector3D.h"

extern CProxy_TreePiece treeProxy;
CkpvExtern(int, _lb_obj_index);

using namespace std;

CreateLBFunc_Def(Orb3dLB, "3d ORB mapping of tree piece space onto 3d processor mesh");

static int comparx(const void *a, const void *b){
  const TPObject *ta = (const TPObject*)(a);
  const TPObject *tb = (const TPObject*)(b);
  return ta->centroid.x < tb->centroid.x ? -1 : ta->centroid.x > tb->centroid.x ? 1 : 0;
}
static int compary(const void *a, const void *b){
  const TPObject *ta = (const TPObject*)(a);
  const TPObject *tb = (const TPObject*)(b);
  return ta->centroid.y < tb->centroid.y ? -1 : ta->centroid.y > tb->centroid.y ? 1 : 0;
}
static int comparz(const void *a, const void *b){
  const TPObject *ta = (const TPObject*)(a);
  const TPObject *tb = (const TPObject*)(b);
  return ta->centroid.z < tb->centroid.z ? -1 : ta->centroid.z > tb->centroid.z ? 1 : 0;
}

static int pcx(const void *a, const void *b){
  const Node *ta = (const Node*)(a);
  const Node *tb = (const Node*)(b);
  return ta->x < tb->x ? -1 : ta->x > tb->x ? 1 : 0;
}
static int pcy(const void *a, const void *b){
  const Node *ta = (const Node*)(a);
  const Node *tb = (const Node*)(b);
  return ta->y < tb->y ? -1 : ta->y > tb->y ? 1 : 0;
}
static int pcz(const void *a, const void *b){
  const Node *ta = (const Node*)(a);
  const Node *tb = (const Node*)(b);
  return ta->z < tb->z ? -1 : ta->z > tb->z ? 1 : 0;
}

void Orb3dLB::init() {
  lbname = "Orb3dLB";
  if (CkpvAccess(_lb_obj_index) == -1)
    CkpvAccess(_lb_obj_index) = LBRegisterObjUserData(sizeof(TaggedVector3D));

  compares[0] = comparx;
  compares[1] = compary;
  compares[2] = comparz;

  pc[0] = pcx;
  pc[1] = pcy;
  pc[2] = pcz;
}

Orb3dLB::Orb3dLB(const CkLBOptions &opt): CBase_Orb3dLB(opt) {
  init();
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
}

//jetley
bool Orb3dLB::QueryBalanceNow(int step){
  // Now, by the first step, we already have
  // centroids. of course, since no gravity has
  // been calculated yet, load can only be
  // estimated via num particles in treepiece
  if(CkMyPe() == 0)
    CkPrintf("Orb3dLB: Step %d\n", step);
  if(step == 0){
    return false;
  }
  else{
    return true;
  }

}

void Orb3dLB::work(BaseLB::LDStats* stats)
{
  int numobjs = stats->n_objs;

  CkPrintf("[orb3dlb] %d objects allocating %d bytes for tp\n", numobjs, numobjs*sizeof(TPObject));
  tps.resize(numobjs);

  for(int i = 0; i < numobjs; i++){
    if(!stats->objData[i].migratable) continue;

    LDObjData &odata = stats->objData[i];
    TaggedVector3D* udata = (TaggedVector3D *)odata.getUserData(CkpvAccess(_lb_obj_index));

    tps[i].centroid = udata->vec;
    tps[i].migratable = stats->objData[i].migratable;
    if(step() == 0){
      tps[i].load = udata->myNumParticles;
    }
    else{
      tps[i].load = stats->objData[i].wallTime;
    }
    tps[i].lbindex = i;
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
  nodes.resize(numnodes);

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
  }

  CkPrintf("[orb3dlb] map\n");
  int tpstart = 0;
  int tpend = numobjs;
  int nodestart = 0;
  int nodeend = numnodes;
  map(tpstart,tpend,nodestart,nodeend,nx,ny,nz,dim);

#ifdef PRINT_BOUNDING_BOXES
  for(int i = 0; i < numnodes; i++){
    CkPrintf("bb of node %d %f %f %f %f %f %f\n", i,
                    nodes[i].box.lesser_corner.x,
                    nodes[i].box.lesser_corner.y,
                    nodes[i].box.lesser_corner.z,
                    nodes[i].box.greater_corner.x,
                    nodes[i].box.greater_corner.y,
                    nodes[i].box.greater_corner.z
                    );
  }
#endif
}

void Orb3dLB::map(int tpstart, int tpend, int nodestart, int nodeend, int xs, int ys, int zs, int dim){
  int nn = nodeend-nodestart;
  CkAssert(nn > 0);
  if(nn == 1){
    directMap(tpstart,tpend,nodestart,nodeend);
  }
  else{
    int totalTp = tpend-tpstart;

    float tpload = 0.0;
    for(int i = tpstart; i < tpend; i++){
      tpload += tps[i].load;
    }

    int numProcRanks = nodes[nodestart].procRanks.length();
    for(int i = nodestart+1; i < nodeend; i++){
      CkAssert(nodes[i].procRanks.length() == numProcRanks);
    }

    qsort(tps.getVec()+tpstart,totalTp,sizeof(TPObject),compares[dim]);
    qsort(nodes.getVec()+nodestart,nn,sizeof(Node),pc[dim]);
    // tp and ntp are modified to hold the particulars of
    // the left/dn/near partitions
    // tp2 and totalTp-ntp hold the objects in the
    // right/up/far partitions
    int tpmid;
    int nodemid;
    partitionEvenLoad(tpstart, tpend, tpmid);
    halveNodes(nodestart, nodeend, nodemid);

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
    map(tpstart,tpmid,nodestart,nodemid,xs,ys,zs,d);
    map(tpmid,tpend,nodemid,nodeend,xs,ys,zs,d);
  }
}

#define ZERO_THRESHOLD 0.00001

void Orb3dLB::directMap(int tpstart, int tpend, int nodestart, int nodeend){
  std::priority_queue<TPObject> pq_obj;
  std::priority_queue<Processor> pq_proc;

  float load = 0.0;
  CkAssert(nodestart==(nodeend-1));
  for(int i = tpstart; i < tpend; i++){
    load += tps[i].load;
    pq_obj.push(tps[i]);
  }

  for(int i = 0; i < procsPerNode; i++){
    Processor p;
    p.load = 0.0;
    p.t = i;
    pq_proc.push(p);
  }

  // int currentZeroProc = 0;
  while(!pq_obj.empty()){
    TPObject tp = pq_obj.top();
    pq_obj.pop();

    // spread around zero-load objects
    // disabled to reduce the number of migrations, and
    // check whether this might solve the BG/P crash
    if(tp.load < ZERO_THRESHOLD){
    }
    else{
      // if object has some non-zero load, assign it to a proc greedily
      Processor p = pq_proc.top();
      pq_proc.pop();

      p.load += tp.load;
      (*mapping)[tp.lbindex] = nodes[nodestart].procRanks[p.t];
#ifdef PRINT_BOUNDING_BOXES
      nodes[nodestart].box.grow(tp.centroid);
#endif

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

void Orb3dLB::partitionEvenLoad(int tpstart, int tpend, int &tpmid){
  float totalLoad = 0.0;
  for(int i = tpstart; i < tpend; i++){
    totalLoad += tps[i].load;
  }
  float lload = 0.0;
  float rload = totalLoad;
  float prevDiff = lload-rload;
  if(prevDiff < 0.0){
    prevDiff = -prevDiff;
  }

  int consider;
  for(consider = tpstart; consider < tpend;){
    float newll = lload + tps[consider].load;
    float newrl = rload - tps[consider].load;

    float newdiff = newll-newrl;
    if(newdiff < 0.0){
      newdiff = -newdiff;
    }

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

  tpmid = consider;
}

void Orb3dLB::halveNodes(int nodestart, int nodeend, int &nodemid){
  int np = nodeend-nodestart;
  nodemid = nodestart+np/2;
}

#include "Orb3dLB.def.h"
