#include <charm++.h>
#include "cklists.h"
#include "HilbertLB.h"
#include "ParallelGravity.h"
#include "TopoManager.h"
#include "Vector3D.h"
#include "hilbert.h"

#define UNIVERSAL_BIT 0x8000000000000000
#define ALPHA 2

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

  CkPrintf("[%d] HilbertLB work\n", CkMyPe());
  int numobjs = stats->n_objs;
  int nmig = stats->n_migrateobjs;
  mapping = &stats->to_proc;

  tps = new TPObject[numobjs];
  stats->makeCommHash();

  CkAssert(nrecvd == numobjs);

  /* Find the bounding box of all centroids. */
  OrientedBox<double> univBB;
  getBoundingBox(univBB);
  CkPrintf("[%d] HilbertLB univ bb %g %g %g %g %g %g\n", 
                CkMyPe(),
                univBB.lesser_corner.x,
                univBB.lesser_corner.y,
                univBB.lesser_corner.z,
                univBB.greater_corner.x,
                univBB.greater_corner.y,
                univBB.greater_corner.z
                  );

  float xres, yres, zres;	

  /* Each tp is put into a bin in each dimension */
  xbin = new CmiUInt8[numobjs];
  ybin = new CmiUInt8[numobjs];
  zbin = new CmiUInt8[numobjs];

  numxbins = numybins = numzbins = (1<<21);

  xres = (univBB.greater_corner.x - univBB.lesser_corner.x) / numxbins;
  yres = (univBB.greater_corner.y - univBB.lesser_corner.y) / numybins;
  zres = (univBB.greater_corner.z - univBB.lesser_corner.z) / numzbins;

  bit_mask = 1;
  int x = numxbins; /* Assuming same number of bins for each dimension */
  numshifts = -1;
  while(x != 0){
    x >>= 1;
    numshifts++;
    bit_mask <<= 1;
  }

  CkPrintf("[%d] HilbertLB start setup \n", CkMyPe());

  CkReduction::setElement *cur = tpCentroids;
  while(cur != NULL){
    TaggedVector3D *data = (TaggedVector3D *) cur->data;
    LDObjHandle &handle = data->handle;
    int tag = stats->getHash(handle.id,handle.omhandle.id);
    tps[tag].centroid.x = data->vec.x;
    tps[tag].centroid.y = data->vec.y;
    tps[tag].centroid.z = data->vec.z;
    tps[tag].migratable = stats->objData[tag].migratable;
    if(step() == 0){
      tps[tag].load = data->myNumParticles;
    }
    else{
      tps[tag].load = stats->objData[tag].wallTime;
    }

    totalLoad += tps[tag].load;
    tps[tag].lbindex = tag;
    tps[tag].key = generateKey(tag);

    xbin[tag] = (tps[tag].centroid.x - univBB.lesser_corner.x)/xres;
    ybin[tag] = (tps[tag].centroid.y - univBB.lesser_corner.y)/yres;
    zbin[tag] = (tps[tag].centroid.z - univBB.lesser_corner.z)/zres;

    cur = cur->next();
  }

  CkPrintf("[%d] HilbertLB done setup \n", CkMyPe());

  loadThreshold = totalLoad / (ALPHA * stats->count);
  qsort(tps,numobjs,sizeof(TPObject),comparekey);
  buildBuckets(0,numobjs);
  bucketList.quickSort();

  CkPrintf("[%d] HilbertLB built tree \n", CkMyPe());

  int index;
  int numProcs = stats->count;
  int toProc = 0;
  int length;
  float currLoad = 0.0;
  float avgLoad = totalLoad / numProcs;

  for(int i = 0; i < bucketList.length(); i++){

    index = bucketList[i].tpStartIndex;
    length = bucketList[i].numobjs;

    for(int j = index; j < length; j++){
      if(currLoad < avgLoad){
        (*mapping)[bucketList[i].tp[j].lbindex] = toProc;
        currLoad += bucketList[i].tp[j].load;
      }
      else{
        toProc++;
        totalLoad -= currLoad;
        numProcs--;
        avgLoad = totalLoad/numProcs;
        currLoad = 0;
        (*mapping)[bucketList[i].tp[j].lbindex] = toProc;
        currLoad += bucketList[i].tp[j].load;
      }
    }
  }

  CkPrintf("[%d] HilbertLB done assignment\n", CkMyPe());
  /* assign buckets to procs based on an average */


  /* Clean-up */

  delete[] tps;
  delete[] xbin;
  delete[] ybin;
  delete[] zbin;

  bucketList.length() = 0;
}

void HilbertLB::getBoundingBox(OrientedBox<double> &univBB){
  CkReduction::setElement *cur = tpCentroids;
  while(cur != NULL){
    TaggedVector3D *data = (TaggedVector3D *)cur;
    univBB.grow(data->vec);
    cur = cur->next();
  }
}

CmiUInt8 HilbertLB::generateKey(int i){
	CmiUInt8 key;
	int mask = bit_mask;
	
	/* For each centroid: use the <x,y,z> bin numbers to generate the 64bit key by interleaving bits.
	   prepend the key with a '1' in order to distinguish between root and children.
	*/
		key = 0;
		while(mask != 0){
			key |= (xbin[i] & mask);
			key <<= 1;
			key |= (ybin[i] & mask);
			key <<= 1;
			key |= (zbin[i] & mask);
			mask >>= 1;			
		}
		key |= UNIVERSAL_BIT;
	return key;
}

void HilbertLB::buildBuckets(int index, int numobjs){
	float currLoad = 0.0;
	for(int i = index; i < index + numobjs; i++){
		currLoad += tps[i].load;
	}
	if(currLoad < loadThreshold){
		LBBucket b;
		b.setLoad(currLoad);
		b.setTP(&(tps[index]),numobjs);
		newCentroid(currLoad, b, numobjs);
		b.hilbertID = hilbert3d(b.getx(),b.gety(),b.getz());
		b.setIndex(index);
		bucketList.push_back(b);
		return;
	}
	else{
		buildBuckets(index,numobjs/2);
		buildBuckets(index+(numobjs/2),numobjs-(numobjs/2));
	}
}

void HilbertLB::newCentroid(float totalLoad, LBBucket b, int numobjs){
	float x,y,z;
	x = 0.0;
	y = 0.0;
	z = 0.0;
	for(int i = 0; i < numobjs; i++){
		x += (b.getx() * ((float) b.getLoad() / totalLoad));
		y += (b.gety() * ((float) b.getLoad() / totalLoad));
		z += (b.getz() * ((float) b.getLoad() / totalLoad));
	}
	b.setCentroid(x,y,z);
}

#include "HilbertLB.def.h"
