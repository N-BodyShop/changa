#include <charm++.h>
#include "cklists.h"
#include "Orb3dLB.h"
#include "ParallelGravity.h"
#include "TopoManager.h"
#include "Vector3D.h"

extern CProxy_TreePiece treeProxy;

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

Orb3dLB::Orb3dLB(const CkLBOptions &opt): CentralLB(opt)
{
  lbname = "Orb3dLB";
  centroidsAllocated = false;
  if (CkMyPe() == 0) CkPrintf("[%d] Orb3dLB created\n",CkMyPe());
  compares[0] = comparx;
  compares[1] = compary;
  compares[2] = comparz;
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
  CmiUInt4 path = 0x1;
  int dim = 0;
  TopoManager tmgr;

  procsPerNode = tmgr.getDimNT();

  int nx = tmgr.getDimNX();
  int ny = tmgr.getDimNY();
  int nz = tmgr.getDimNZ();
  int numnodes = nx*ny*nz; 

  int x1[NDIMS] = {0,0,0};
  int x2[NDIMS] = {nx-1,ny-1,nz-1};

  map(tp,numobjs,numnodes,x1,x2,dim);

  for(int i = 0; i < numobjs; i++){
    CkPrintf("%f %f %f %d\n", tp[i].centroid.x, tp[i].centroid.y, tp[i].centroid.z, stats->to_proc[i]);
  }

#ifdef ORB3DLB_VISUALIZE
  CkVec<Vector3D<float> > procCentroids;
  CkVec<int> procNumTPs;
  procCentroids.resize(stats->count);
  procNumTPs.resize(stats->count);
  procCentroids.length() = stats->count;
  procNumTPs.length() = stats->count;
  for(int i = 0; i < stats->count; i++){
    procCentroids[i].x = 0.0;
    procCentroids[i].y = 0.0;
    procCentroids[i].z = 0.0;
    procNumTPs[i] = 0;
  }

  for(int i = 0; i < numobjs; i++){
    CkAssert(i == tp[i].lbindex);
    int proc = stats->to_proc[i];
    procCentroids[proc] += tp[i].centroid;
    procNumTPs[proc]++;
  }

  for(int i = 0; i < stats->count; i++){
    if(procNumTPs[i] > 0){
      procCentroids[i] = procCentroids[i]/(1.0*procNumTPs[i]);
    }
  }

  float minx = procCentroids[0].x;
  float miny = procCentroids[0].y;
  float minz = procCentroids[0].z;

  float maxx = procCentroids[0].x;
  float maxy = procCentroids[0].y;
  float maxz = procCentroids[0].z;

  for(int i = 1; i < stats->count; i++){
    if(procNumTPs[i] > 0){
      if(procCentroids[i].x > maxx){
        maxx = procCentroids[i].x;
      }
      if(procCentroids[i].y > maxy){
        maxy = procCentroids[i].y;
      }
      if(procCentroids[i].z > maxz){
        maxz = procCentroids[i].z;
      }

      if(procCentroids[i].x < minx){
        minx = procCentroids[i].x;
      }
      if(procCentroids[i].y < miny){
        miny = procCentroids[i].y;
      }
      if(procCentroids[i].z < minz){
        minz = procCentroids[i].z;
      }
    }
  }

  maxx -= minx;
  maxy -= miny;
  maxz -= minz;

  Vector3D<float> minv;
  minv.x = minx;
  minv.y = miny;
  minv.z = minz;
  
  Vector3D<float> r;
  r.x = maxx;
  r.y = maxy;
  r.z = maxz;

  for(int i = 0; i < stats->count; i++){
    if(procNumTPs[i] > 0){
      procCentroids[i] -= minv;
      procCentroids[i] /= r;
      CkPrintf("%f %f %f %d\n", procCentroids[i].x, procCentroids[i].y, procCentroids[i].z, i);
    }
  }
#endif

}

void Orb3dLB::map(TPObject *tp, int ntp, int nn, int *x1, int *x2, int dim){
  //CkPrintf("ntp: %d np: %d dim: %d path: 0x%x\n",ntp,np,dim,path);
  if(nn == 1){
    CkAssert(x1[0]==x2[0]);
    CkAssert(x1[1]==x2[1]);
    CkAssert(x1[2]==x2[2]);
    directMap(tp,ntp,x1);
  }
  else{
    int totalTp = ntp;
    // xx0 is the greater corner of the left/dn/near partition
    int xx0[NDIMS];
    // xx1 is the lesser corner of the right/up/far partition
    int xx1[NDIMS];
    
    qsort(tp,ntp,sizeof(TPObject),compares[dim]);
    // tp and ntp are modified to hold the particulars of
    // the left/dn/near partitions
    // tp2 and totalTp-ntp hold the objects in the 
    // right/up/far partitions
    TPObject *tp2 = partitionEvenLoad(tp,ntp);
    int d = nextDim(dim); 
    for(int i = 0; i < NDIMS; i++){
      if(i == dim){
        xx0[dim] = x2[dim]-(x2[dim]-x1[dim]+1)/2; 
        xx1[dim] = x1[dim]+(x2[dim]-x1[dim]+1)/2;
      }
      else{
        xx0[i] = x2[i];
        xx1[i] = x1[i];
      }
    }
    map(tp,ntp,nn/2,x1,xx0,d);
    map(tp2,totalTp-ntp,nn/2,xx1,x2,d);
  }
}
void Orb3dLB::directMap(TPObject *tp, int ntp, int *x){
  /*
  int numshifts = 0;
  CmiUInt4 correctbits = 0x0; 
  CmiUInt4 bit;
  while(path > 0x1){
    bit = path & 0x1;
    path = path >> 1;
    correctbits = correctbits << 1;
    correctbits |= bit;
    numshifts++;
  }
  */

/*
#ifdef USE_TOPOMGR
  TopoManager tm;
#endif
  int min[NDIMS] = {0,0,0};
#ifdef USE_TOPOMGR
  int max[NDIMS] = {tm.getDimNX()-1,tm.getDimNY()-1,tm.getDimNZ()-1};
  int coresPerNode = tm.getDimNT();
#else
  int max[NDIMS] = {XMAX-1,YMAX-1,ZMAX-1};
#endif
  int dim = 0;
  for(int i = 0; i < numshifts; i++){
    bit = correctbits & 0x1;
    if(bit){
      min[dim] = min[dim]+(max[dim]-min[dim]+1)/2;
    }
    else{
      max[dim] = max[dim]-(max[dim]-min[dim]+1)/2;
    }
    dim = (dim+1)%NDIMS;
    correctbits = correctbits >> 1;
  }

  CkAssert(min[0] == max[0]);
  CkAssert(min[1] == max[1]);
  CkAssert(min[2] == max[2]);
*/

  CkPrintf("[Orb3dLB] mapping %d objects to node (%d,%d,%d)\n", ntp, x[0],x[1],x[2]);
  int t = 0;
  TopoManager tmgr;
  for(int i = 0; i < ntp; i++){
    (*mapping)[tp[i].lbindex] = tmgr.coordinatesToRank(x[0],x[1],x[2],t);
    t = (t+1)%procsPerNode;
  }
}

int Orb3dLB::nextDim(int dim){
  return (dim+1)%NDIMS; 
}

TPObject *Orb3dLB::partitionEvenLoad(TPObject *tp, int &ntp){
  float partition1Load = 0.0;
  float totalLoad = 0.0;
  for(int i = 0; i < ntp; i++){
    totalLoad += tp[i].load;
  }
  float halfLoad = totalLoad/2.0;
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

#include "Orb3dLB.def.h"

/*@}*/
