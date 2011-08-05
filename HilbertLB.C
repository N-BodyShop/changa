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
	centroidsAllocated = false;
/* Add other necessary information */
}

void HilbertLB::receiveCentroids(CkReductionMsg *msg){
	if(haveTPCentroids){
		delete tpmsg;
	}
	tpCentroids = (TaggedVector3D *)msg->getData();
	nrecvd = msg->getGcount();
	tpmsg = msg;
	haveTPCentroids = true;
	CkPrintf("HilbertLB: receiveCentroids started: %d elements, msg length: %d\n", msg->getGcount(), msg->getLength()); 
	treeProxy.doAtSync();
	CkPrintf("HilbertLB: receiveCentroids done\n");
}

void HilbertLB::work(BaseLB::LDStats* stats){
	
	int numobjs = stats->n_objs;
	int nmig = stats->n_migrateobjs;
	mapping = &stats->to_proc;

	
	tp = new TPObject[numobjs];
	stats->makeCommHash();

	if(nrecvd != numobjs)
		CkAbort("wrong tpCentroids length\n");

	for(int i = 0; i < numobjs; i++){
		LDObjHandle &handle = tpCentroids[i].handle;
		int tag = stats->getHash(handle.id, handle.omhandle.id);
		tp[tag].centroid.x = tpCentroids[i].vec.x;
		tp[tag].centroid.y = tpCentroids[i].vec.y;
		tp[tag].centroid.z = tpCentroids[i].vec.z;
		tp[tag].migratable = stats->objData[tag].migratable;
		tp[tag].load = tpCentroids[i].myNumParticles;
		totalLoad += tp[tag].load;
		tp[tag].lbindex = tag;
		tp[tag].key = generateKey(tag);
	}
	
	normalizeCoordinates(numobjs);
	
	bit_mask = 1;
	int x = numxbins; /* Assuming same number of bins for each dimension */
	numshifts = -1;
	while(x != 0){
		x >>= 1;
		numshifts++;
		bit_mask <<= 1;
	}

	
	loadThreshold = totalLoad / (ALPHA * stats->count);
	qsort(tp,numobjs,sizeof(TPObject),comparekey);
	buildBuckets(0,numobjs);
	bucketList.quickSort();
	
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

	/* assign buckets to procs based on an average */
	

	/* Clean-up */

	delete[] tp;
	delete[] xbin;
	delete[] ybin;
	delete[] zbin;
}

void HilbertLB::normalizeCoordinates(int numobjs){
	
	float xmin = tpCentroids[0].vec.x;
	float ymin = tpCentroids[0].vec.y;
	float zmin = tpCentroids[0].vec.z;
	float xmax = tpCentroids[0].vec.x;
	float ymax = tpCentroids[0].vec.y;
	float zmax = tpCentroids[0].vec.z;
	
	float xres, yres, zres;	

	xbin = new CmiUInt8[numobjs];
	ybin = new CmiUInt8[numobjs];
	zbin = new CmiUInt8[numobjs];
	
	/* Find the coordinates of the bounding box(mins and max) for the centroids. */
	for(int i = 1; i < numobjs; i++){
		xmin = (tpCentroids[i].vec.x < xmin) ? tpCentroids[i].vec.x : xmin;
		xmax = (tpCentroids[i].vec.x > xmax) ? tpCentroids[i].vec.x : xmax;
		ymin = (tpCentroids[i].vec.y < ymin) ? tpCentroids[i].vec.y : ymin;
		ymax = (tpCentroids[i].vec.y > ymax) ? tpCentroids[i].vec.y : ymax;
		zmin = (tpCentroids[i].vec.z < zmin) ? tpCentroids[i].vec.z : zmin;
		zmax = (tpCentroids[i].vec.z > zmax) ? tpCentroids[i].vec.z : zmax;
	}
	
	/* Calculate the number of cells for each dimension */
	/*	numxbins = (xmax - xmin)/resolution;
		numybins = (ymax - ymin)/resolution;
		numzbins = (zmax - zmin)/resolution;
	*/
		numxbins = 2097151;
		numybins = 2097151;
		numzbins = 2097151;
		
		xres = (xmax - xmin) / numxbins;
		yres = (ymax - ymin) / numybins;
		zres = (zmax - zmin) / numzbins;
	
	/* Calculate which cell the centroid belongs to for each dimension */
	for(int i = 0; i < numobjs; i++){
		xbin[i] = (long)((tpCentroids[i].vec.x - xmin) / xres);
		ybin[i] = (long)((tpCentroids[i].vec.y - ymin) / yres);
		zbin[i] = (long)((tpCentroids[i].vec.z - zmin) / zres);
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
		currLoad += tp[i].load;
	}
	if(currLoad < loadThreshold){
		LBBucket b;
		b.setLoad(currLoad);
		b.setTP(&(tp[index]),numobjs);
		newCentroid(currLoad, b, numobjs);
		b.hilbertID = hilbert3d(b.getx(),b.gety(),b.getz());
		b.setIndex(index);
		bucketList.insertAtEnd(b);
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
