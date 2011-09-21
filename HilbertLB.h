#ifndef _HILBERTLB_H_
#define _HILBERTLB_H_

#include "CentralLB.h"
#include "MapStructures.h"
#include "HilbertLB.decl.h"
#include "TaggedVector3D.h"
#include <queue>
#include "Bucket.h"

#define ZERO_LOAD_THRESHOLD 1e-4

typedef CmiUInt8 HKey;

bool intersection(OrientedBox<float> &b1, OrientedBox<float> &b2);

struct HilbertNode {
  HKey key;
  CkVec<HilbertNode*> children;

  int tpStart;
  int tpEnd;
  OrientedBox<float> box;

  HilbertNode(HKey k, int tpstart, int tpend, TPObject *tps){
    tpStart = tpstart;
    tpEnd = tpend;
    key = k;
    for(int i = tpStart; i < tpEnd; i++){
      TPObject &tp = tps[i];

      if(tp.load < ZERO_LOAD_THRESHOLD) continue;
      box.grow(tp.centroid);
    }
  }

  void split(int i, TPObject *tps){
    children.push_back(new HilbertNode(key<<1,tpStart,i,tps));
    children.push_back(new HilbertNode((key<<1)+1,i,tpEnd,tps));

    if(intersection(children[0]->box,children[1]->box)){
      CkPrintf("hkey %ld and %ld intersect bbs %f %f %f %f %f %f %f %f %f %f %f %f\n",
              children[0]->key,
              children[1]->key,
              children[0]->box.lesser_corner.x,
              children[0]->box.lesser_corner.y,
              children[0]->box.lesser_corner.z,
              children[0]->box.greater_corner.x,
              children[0]->box.greater_corner.y,
              children[0]->box.greater_corner.z,
              children[1]->box.lesser_corner.x,
              children[1]->box.lesser_corner.y,
              children[1]->box.lesser_corner.z,
              children[1]->box.greater_corner.x,
              children[1]->box.greater_corner.y,
              children[1]->box.greater_corner.z
              );
    }
  }

  HilbertNode *getChild(int i){
    return children[i];
  }
};

class HilbertLB : public CentralLB{

private:
	bool haveTPCentroids;
	int nrecvd, numxbins,numybins, numzbins, numshifts;
        CkReduction::setElement *tpCentroids;
	CkReductionMsg *tpmsg;	
	CkVec<TPObject> tps;
	float loadThreshold;
	float totalLoad;	

	CkVec<LBBucket> bucketList; 
	CkVec<int> *mapping;
        CkVec<float> peload;
        CkVec<OrientedBox<float> > pebb;

        void getBoundingBox(OrientedBox<float> &univBB);
        int binary_search_ge(CmiUInt8 check, TPObject *treePieces, int start, int end);

        HilbertNode *root;

public:
	HilbertLB(const CkLBOptions &);
	HilbertLB(CkMigrateMessage *m):CentralLB(m) {lbname = "HilbertLB";}
	void work(BaseLB::LDStats*);
	void receiveCentroids(CkReductionMsg *msg);
	void normalizeCoordinates(int numobjs);
	CmiUInt8 generateKey(int xbin, int ybin, int zbin);
        void buildBuckets(CmiUInt8 key, int depth, int start, int end, HilbertNode *node);
	void mapBucket(int bucketIndex, int toProc);
};

#endif
