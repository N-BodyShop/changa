#include <charm++.h>
#include "cklists.h"
#include "HilbertLB.h"
#include "ParallelGravity.h"
#include "TopoManager.h"
#include "Vector3D.h"
#include "hilbert.h"

#define ALPHA 4
#define HILBERTLB_DEBUG(...) CkPrintf(__VA_ARGS__);

//#define HILBERTLB_DEBUG(...) fprintf(stderr,__VA_ARGS__);

using namespace std;

CreateLBFunc_Def(HilbertLB, "Hilbert Order load balancing");

extern CProxy_TreePiece treeProxy;

int comparekey(const void* a, const void* b){
	const TPObject *ta = (const TPObject*)a;
	const TPObject *tb = (const TPObject*)b;
	if((*ta).key == (*tb).key){return 0;}
	else if((*ta).key < (*tb).key){return -1;}
	else{return 1;}
}

HilbertLB::HilbertLB(const CkLBOptions &opt): CentralLB(opt)
{
	lbname = "HilbertLB";
/* Add other necessary information */
}

void HilbertLB::receiveCentroids(CkReductionMsg *msg){
  if(haveTPCentroids){
    delete tpmsg;
  }
  tpCentroids = (CkReduction::setElement *)msg->getData();
  CkReduction::setElement *cur = tpCentroids;
  nrecvd = 0;
  while(cur != NULL){
    CkAssert(cur->dataSize == sizeof(TaggedVector3D));
    nrecvd++;
    cur = cur->next();
  }
  tpmsg = msg;
  haveTPCentroids = true;
  CkPrintf("HilbertLB: receiveCentroids %d elements, msg length %d\n", nrecvd, msg->getLength()); 
  treeProxy.doAtSync();
}

void HilbertLB::work(BaseLB::LDStats* stats){

  HILBERTLB_DEBUG("[%d] HilbertLB work\n", CkMyPe());
  int numobjs = stats->n_objs;
  int nmig = stats->n_migrateobjs;
  mapping = &stats->to_proc;

  tps.resize(numobjs);
  stats->makeCommHash();

  CkAssert(nrecvd == numobjs);

  /* Find the bounding box of all centroids. */
  OrientedBox<float> univBB;
  getBoundingBox(univBB);
  HILBERTLB_DEBUG("univ bb %f %f %f %f %f %f\n", 
      univBB.lesser_corner.x,
      univBB.lesser_corner.y,
      univBB.lesser_corner.z,
      univBB.greater_corner.x,
      univBB.greater_corner.y,
      univBB.greater_corner.z
      );

  float xsz, ysz, zsz;	

  /* Each tp is put into a bin in each dimension */
  numxbins = numybins = numzbins = (1<<21);

  xsz = (univBB.greater_corner.x-univBB.lesser_corner.x);
  ysz = (univBB.greater_corner.y-univBB.lesser_corner.y);
  zsz = (univBB.greater_corner.z-univBB.lesser_corner.z);

  HILBERTLB_DEBUG("[%d] HilbertLB start setup \n", CkMyPe());

  CkReduction::setElement *cur = tpCentroids;
  while(cur != NULL){
    TaggedVector3D *data = (TaggedVector3D *) cur->data;
    LDObjHandle &handle = data->handle;
    int tag = stats->getHash(handle.id,handle.omhandle.id);
    tps[tag].centroid = data->vec;

    // skip zero-load tree pieces: keep them where they are

    /*
    CkPrintf("tree piece %d box %f %f %f %f %f %f\n",
              data->tp, 
              data->box.lesser_corner.x,
              data->box.lesser_corner.y,
              data->box.lesser_corner.z,
              data->box.greater_corner.x,
              data->box.greater_corner.y,
              data->box.greater_corner.z
            );
    */
    tps[tag].migratable = stats->objData[tag].migratable;
    if(step() == 0){
      tps[tag].load = data->myNumParticles;
    }
    else{
      tps[tag].load = stats->objData[tag].wallTime;
    }

    totalLoad += tps[tag].load;
    tps[tag].lbindex = tag;
    tps[tag].tpindex = data->tp;

    int xbin = ((tps[tag].centroid.x-univBB.lesser_corner.x)/xsz)*(1.0*numxbins);
    int ybin = ((tps[tag].centroid.y-univBB.lesser_corner.y)/ysz)*(1.0*numxbins);
    int zbin = ((tps[tag].centroid.z-univBB.lesser_corner.z)/zsz)*(1.0*numxbins);

    /*
    HILBERTLB_DEBUG("tp %d tag %d centroid %f %f %f bin %d %d %d max %d\n",
                    data->tp, tag,
                    tps[tag].centroid.x,
                    tps[tag].centroid.y,
                    tps[tag].centroid.z,
                    xbin,
                    ybin,
                    zbin,
                    numxbins
                    );
    */
    tps[tag].key = generateKey(xbin,ybin,zbin);

    //HILBERTLB_DEBUG("key %lu\n", tps[tag].key);

    cur = cur->next();
  }

  HILBERTLB_DEBUG("[%d] HilbertLB done setup \n", CkMyPe());

  // sort tp's by keys; the key for each tp
  // is obtained from the bin it was placed in
  qsort(tps.getVec(),numobjs,sizeof(TPObject),comparekey);

  TPObject *pointer = tps.getVec();
  /*
  for(int i = 0; i < numobjs; i++){
    if(pointer->load > ZERO_LOAD_THRESHOLD){
      HILBERTLB_DEBUG("TP %d id %d centroid %f %f %f key %lx load %f\n",
          i,
          pointer->tpindex,
          pointer->centroid.x,
          pointer->centroid.y,
          pointer->centroid.z,
          pointer->key,
          pointer->load
          );
    }
    pointer++;
  }
  */


  loadThreshold = totalLoad/(1.0*ALPHA*stats->count);

  // construct tree and obtain buckets thereof
  // start by placing all particles (indices [0..numobjs))
  // under root node (key 1, depth 0)

  
  root = new HilbertNode(HKey(1),0,numobjs,tps.getVec());
  buildBuckets(CmiUInt8(1),0,0,numobjs,root);

  HILBERTLB_DEBUG("LoadThreshold %f totalLoad %f buckets %d\n",loadThreshold, totalLoad, bucketList.length());

  // sort by Hilbert keys
  bucketList.quickSort();

  for(int i = 0; i < bucketList.length(); i++){
    /*
    HILBERTLB_DEBUG("BUCKET bb %d %f %f %f %f %f %f %d %f\n",
                     i, 
                     bucketList[i].box.lesser_corner.x,
                     bucketList[i].box.lesser_corner.y,
                     bucketList[i].box.lesser_corner.z,
                     bucketList[i].box.greater_corner.x,
                     bucketList[i].box.greater_corner.y,
                     bucketList[i].box.greater_corner.z,
                     bucketList[i].numobjs,
                     bucketList[i].load
                     );
    */

    //CkPrintf("bucket %d hilbertid %lx\n", i, bucketList[i].hilbertID);

    TPObject *pointer = bucketList[i].tp;
    for(int j = 0; j < bucketList[i].numobjs; j++){
      if(pointer->load > ZERO_LOAD_THRESHOLD){
        /*
        HILBERTLB_DEBUG("TP %d bucket %d centroid %f %f %f\n",
            pointer->tpindex,
            i,
            pointer->centroid.x,
            pointer->centroid.y,
            pointer->centroid.z
            );
        */
      }
      pointer++;
    }
  }

  //CkPrintf("[%d] HilbertLB built tree Bucket Len= %d\n", CkMyPe(), bucketList.length());

  int index;
  int numProcs = stats->count;
  int toProc = 0;
  int length;
  float currLoad = 0.0;
  float avgLoad = totalLoad/numProcs;
  int mapCount=0;
  CkPrintf("Average Load: %f numProcs: %d\n",avgLoad,numProcs);

  /* assign buckets to procs based on an average */
  int numPes = stats->count;
  peload.resize(numPes);
  pebb.resize(numPes);

  for(int i = 0; i < numPes; i++){
    peload[i] = 0.0;
    pebb[i].reset();
  }

  for(int i = 0; i < bucketList.length();){
    if(currLoad < avgLoad){
      currLoad += bucketList[i].load;
      mapCount++;
      mapBucket(i,toProc);
      i++;
      //CkPrintf("HilbertLB: Assigned Load:%f CurrentLoad: %f  Avg Load:%f toProc:%d\n",bucketList[i].load,currLoad,avgLoad,toProc);
    }
    else{
      mapCount = 0;
      toProc++;
      totalLoad -= currLoad;
      numProcs--;
      avgLoad = totalLoad/numProcs;
      currLoad = 0;
    }
  }

  for(int i = 0; i < numPes; i++){
    CkPrintf("pe %d load %f bb %f %f %f %f %f %f\n", 
             i, peload[i],
             pebb[i].lesser_corner.x,
             pebb[i].lesser_corner.y,
             pebb[i].lesser_corner.z,
             pebb[i].greater_corner.x,
             pebb[i].greater_corner.y,
             pebb[i].greater_corner.z
             );
  }

  /*
  for(int i = 0; i < numobjs; i++){
    CkPrintf("obj %d goes to proc %d\n", i, stats->to_proc[i]);
  }

  CkPrintf("[%d] HilbertLB done assignment\n", CkMyPe());
  */

  /* Clean-up */

  totalLoad = 0;
  bucketList.length() = 0;
}

void HilbertLB::mapBucket(int bucketIndex, int toProc){
  int length = bucketList[bucketIndex].numobjs;
  for(int i = 0; i < length; i++){
    TPObject &tp = bucketList[bucketIndex].tp[i];
    (*mapping)[tp.lbindex] = toProc;

    CkPrintf("proc %d gets tp %d\n",toProc,tp.tpindex);

    peload[toProc] += tp.load;
    pebb[toProc].grow(tp.centroid);
  }
}

void HilbertLB::getBoundingBox(OrientedBox<float> &univBB){
  CkReduction::setElement *cur = tpCentroids;
  while(cur != NULL){
    TaggedVector3D *data = (TaggedVector3D *)cur->data;
    univBB.grow(data->vec);
    cur = cur->next();
  }

  float pad = 0.001;
  univBB.greater_corner = univBB.greater_corner*pad+univBB.greater_corner;
  univBB.lesser_corner = univBB.lesser_corner-pad*univBB.lesser_corner;

}

CmiUInt8 HilbertLB::generateKey(int xbin, int ybin, int zbin){
  CmiUInt8 prepend = CmiUInt8(1);
  prepend <<= 63;
  CmiUInt8 mask = CmiUInt8(0x1);
  CmiUInt8 k = CmiUInt8(0x0);
  int shiftBy = 0;

  CkAssert(xbin < (1 << 21));
  CkAssert(ybin < (1 << 21));
  CkAssert(zbin < (1 << 21));

  for(int j = 0; j < 21; j++){
    k |= ((zbin & mask) <<  shiftBy);
    k |= ((ybin & mask) << (shiftBy+1));
    k |= ((xbin & mask) << (shiftBy+2));
    mask <<= 1;
    // minus 1 because mask itself has shifted
    // left by one position
    shiftBy += 2;
  }
  k |= prepend;
  return k;
}


void HilbertLB::buildBuckets(CmiUInt8 key, int depth, int start, int end, HilbertNode *parent){
  int numTps = end-start;

  if(numTps == 0) return;

  float currLoad = 0.0;
  for(int i = start; i < end; i++){
    currLoad += tps[i].load;
  }

  CkAssert(depth < 64);

  //HILBERTLB_DEBUG("key %ld depth %d start %d end %d load %f\n", key, depth, start, end, currLoad);

  if(currLoad < loadThreshold || numTps == 1){
    LBBucket b;
    b.load = currLoad;
    b.tp = &(tps[start]);
    b.numobjs = numTps;

    for(int i = start; i < end; i++){
      // sum not weighted
      if(tps[i].load < ZERO_LOAD_THRESHOLD) continue;
      b.centroid += tps[i].centroid;
      b.box.grow(tps[i].box);
    }
    b.centroid /= numTps;

    b.hilbertID = hilbert3d(b.centroid.x, b.centroid.y, b.centroid.z);
    b.tpStartIndex = start;
    bucketList.push_back(b);
    /*
    HILBERTLB_DEBUG("MAKE BUCKET: key %ld currLoad: %f loadThreshold: %f\n", 
                     key, depth, start, end,
                     currLoad,loadThreshold);
    */

    return;
  }

  // otherwise, too heavy, need to split into two children
  // 1. find split point: first tp whose key falls under second child
  CmiUInt8 firstChildKey = (key<<1);
  CmiUInt8 secondChildKey = firstChildKey+1;
  int childDepth = depth+1;
  CmiUInt8 testKey = (secondChildKey << (64-(1*childDepth+1)));

  int splitIndex = binary_search_ge(testKey,tps.getVec(),start,end);
  parent->split(splitIndex,tps.getVec());

  buildBuckets(firstChildKey,childDepth,start,splitIndex,parent->getChild(0));
  buildBuckets(secondChildKey,childDepth,splitIndex,end,parent->getChild(1));
}

int HilbertLB::binary_search_ge(CmiUInt8 check, TPObject *treePieces, int start, int end){
  int lo = start;
  int hi = end;
  int mid;
  while(lo < hi){
    mid = lo+((hi-lo)>>1);
    if(treePieces[mid].key >= check){
      hi = mid;
    }
    else{
      lo = mid+1;
    }
  }
  return lo;
}

bool intersection(OrientedBox<float> &b1, OrientedBox<float> &b2){
  Vector3D<float> displacement = b1.greater_corner-b1.lesser_corner;

  bool ret = false;

  Vector3D<float> test = b1.lesser_corner;
  ret |= b2.contains(test);

  test.x += displacement.x;
  ret |= b2.contains(test);

  test.y += displacement.y;
  ret |= b2.contains(test);

  test.x -= displacement.x; 
  ret |= b2.contains(test);

  test.z += displacement.z;
  ret |= b2.contains(test);

  test.x += displacement.x;
  ret |= b2.contains(test);

  test.y -= displacement.y;
  ret |= b2.contains(test);

  test.x -= displacement.x;
  ret |= b2.contains(test);

  return ret;
}


#include "HilbertLB.def.h"
