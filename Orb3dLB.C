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

int pcx(const void *a, const void *b){
  Proc *ta = (Proc *)a;
  Proc *tb = (Proc *)b;
  return (int)(ta->x-tb->x);
}
int pcy(const void *a, const void *b){
  Proc *ta = (Proc *)a;
  Proc *tb = (Proc *)b;
  return (int)(ta->y-tb->y);
}
int pcz(const void *a, const void *b){
  Proc *ta = (Proc *)a;
  Proc *tb = (Proc *)b;
  return (int)(ta->z-tb->z);
}

Orb3dLB::Orb3dLB(const CkLBOptions &opt): CentralLB(opt)
{
  lbname = "Orb3dLB";
  centroidsAllocated = false;
  if (CkMyPe() == 0) CkPrintf("[%d] Orb3dLB created\n",CkMyPe());
  compares[0] = comparx;
  compares[1] = compary;
  compares[2] = comparz;

  pc[0] = pcx;
  pc[1] = pcy;
  pc[2] = pcz;
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

  Proc *procs = new Proc[stats->count];

  for(int i = 0; i < stats->count; i++){
    procs[i].rank = i;
    tmgr.rankToCoordinates(i,procs[i].x,procs[i].y,procs[i].z,procs[i].t);
  }

  map(tp,numobjs,numnodes,procs,dim);

  float *procload = new float[stats->count];
  for(int i = 0; i < stats->count; i++){
    procload[i] = 0.0;
  }

  for(int i = 0; i < numobjs; i++){
    CkPrintf("obj %d to %d (%f %f %f)\n", stats->objData[i].id().id[0], stats->to_proc[i], tp[i].centroid.x, tp[i].centroid.y, tp[i].centroid.z);
    procload[stats->to_proc[i]] += tp[i].load;
  }

  for(int i = 0; i < stats->count; i++){
    CkPrintf("proc %d load %f\n", i, procload[i]);
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

void Orb3dLB::map(TPObject *tp, int ntp, int np, Proc *procs, int dim){
  //CkPrintf("ntp: %d np: %d dim: %d path: 0x%x\n",ntp,np,dim,path);
  if(np == 1){
    directMap(tp,ntp,procs);
  }
  else{
    int totalTp = ntp;
    // xx0 is the greater corner of the left/dn/near partition
    int xx0[NDIMS];
    // xx1 is the lesser corner of the right/up/far partition
    int xx1[NDIMS];
    
    qsort(tp,ntp,sizeof(TPObject),compares[dim]);
    qsort(procs,np,sizeof(Proc),pc[dim]);
    // tp and ntp are modified to hold the particulars of
    // the left/dn/near partitions
    // tp2 and totalTp-ntp hold the objects in the 
    // right/up/far partitions
    TPObject *tp2 = partitionEvenLoad(tp,ntp);
    Proc *procs2 = halveProcessors(procs,np);
    int d = nextDim(dim); 
    map(tp,ntp,np/2,procs,d);
    map(tp2,totalTp-ntp,np/2,procs2,d);
  }
}
void Orb3dLB::directMap(TPObject *tp, int ntp, Proc *procs){
  CkPrintf("[Orb3dLB] mapping %d objects to proc (%d,%d,%d,%d)\n", ntp, procs[0].x, procs[0].y, procs[0].z, procs[0].t);
  TopoManager tmgr;
  for(int i = 0; i < ntp; i++){
    (*mapping)[tp[i].lbindex] = procs[0].rank;
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

Proc *Orb3dLB::halveProcessors(Proc *start, int np){
  Proc *ret = start;
  for(int i = 0; i < np/2; i++,ret=ret+1);
  return ret;
}

#include "Orb3dLB.def.h"

/*@}*/
