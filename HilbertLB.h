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
	bool haveTPCentroids;
	int nrecvd, numxbins,numybins, numzbins, numshifts, bit_mask;
	CmiUInt8 *xbin, *ybin, *zbin;
        CkReduction::setElement *tpCentroids;
	CkReductionMsg *tpmsg;	
	TPObject *tps;
	float loadThreshold;
	float totalLoad;	

	CkVec<LBBucket> bucketList; 
	CkVec<int> *mapping;

        void getBoundingBox(OrientedBox<double> &univBB);

public:
	HilbertLB(const CkLBOptions &);
	HilbertLB(CkMigrateMessage *m):CentralLB(m) {lbname = "HilbertLB";}
	void work(BaseLB::LDStats*);
	void receiveCentroids(CkReductionMsg *msg);
	void normalizeCoordinates(int numobjs);
	CmiUInt8 generateKey(int i);
	void buildBuckets(int index, int numobjs);
	void mapBuckets(int bucketIndex, int toProc);
	void newCentroid(float totalLoad, LBBucket b, int numobjs);
};

#endif
