#ifndef _HILBERTLB_H_
#define _HILBERTLB_H_

#include "CentralLB.h"
#include "MapStructures.h"
#include "HilbertLB.decl.h"
#include "TaggedVector3D.h"
#include <queue>
#include "Bucket.h"

class HilbertLB : public CentralLB{

private:
	TaggedVector3D *tpCentroids;
	bool haveTPCentroids;
	int nrecvd, numxbins,numybins, numzbins, numshifts, bit_mask;
	long *xbin, *ybin, *zbin;
	CkReductionMsg *tpmsg;	
	TPObject *tp;
	CmiBool centroidsAllocated;
	float loadThreshold;
	float totalLoad;	

	CkVec<LBBucket> bucketList; 
	CkVec<int> *mapping;

public:
	HilbertLB(const CkLBOptions &);
	HilbertLB(CkMigrateMessage *m):CentralLB(m) {lbname = "HilbertLB";}
	void work(BaseLB::LDStats*);
	void receiveCentroids(CkReductionMsg *msg);
	void normalizeCoordinates(int numobjs);
	long generateKey(int i);
	void buildBuckets(int index, int numobjs);
	void newCentroid(long totalLoad, LBBucket b, int numobjs);
};

#endif
