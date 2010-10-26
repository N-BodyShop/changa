/*****************************************************************************
 * $Source$
 * $Author$
 * $Date$
 * $Revision$
 *****************************************************************************/

/**
 * \addtogroup CkLdb
*/
/*@{*/

#ifndef _MULTISTEPLB_H_
#define _MULTISTEPLB_H_

#define MCLBMS          // multistepping enabled
#define MCLB_ORBSMOOTH  // orbsmooth for large steps
#define MCLB_RR         // round robin otherwise

#include "CentralLB.h"

#include "MapStructures.h"
#include "ScaleTranMapBG.h"
#include "ScaledORBMapBG.h"

#include "MultistepLB.decl.h"

//**************************************
// ORB3DLB functions
//**************************************
static int comparx(const void *a, const void *b){
  TPObject *ta = (TPObject *)a;
  TPObject *tb = (TPObject *)b;
  return (int)(ta->centroid.x-tb->centroid.x);
}
static int compary(const void *a, const void *b){
  TPObject *ta = (TPObject *)a;
  TPObject *tb = (TPObject *)b;
  return (int)(ta->centroid.y-tb->centroid.y);
}
static int comparz(const void *a, const void *b){
  TPObject *ta = (TPObject *)a;
  TPObject *tb = (TPObject *)b;
  return (int)(ta->centroid.z-tb->centroid.z);
}

static int pcx(const void *a, const void *b){
  Node *ta = (Node *)a;
  Node *tb = (Node *)b;
  return (int)(ta->x-tb->x);
}
static int pcy(const void *a, const void *b){
  Node *ta = (Node *)a;
  Node *tb = (Node *)b;
  return (int)(ta->y-tb->y);
}
static int pcz(const void *a, const void *b){
  Node *ta = (Node *)a;
  Node *tb = (Node *)b;
  return (int)(ta->z-tb->z);
}
//**************************************


void CreateMultistepLB();
BaseLB * AllocateMultistepLB();

class WeightObject{
  public:

  int idx;
  double weight;

  bool operator<= (const WeightObject& rhs){
    return weight > rhs.weight;
  }
  bool operator>= (const WeightObject& rhs){
    return weight < rhs.weight;
  }
  WeightObject(int _idx, double _weight) : idx(_idx), weight(_weight){};
  WeightObject() : idx(0), weight(0.0){};
};


// struct to hold information such as 
// cputime, walltime etc. 
struct CondensedLDStats {
  CkVec<double> cpuTime;
  CkVec<double> wallTime;
  CkVec<bool> migratable;
  CkVec<int> tpindex;

  int n_objs;
  int n_migrateobjs; 
  int count;

  public: 
  void init(BaseLB::LDStats *stats){
    int n = stats->n_objs;

    cpuTime.reserve(n);
    wallTime.reserve(n);
    migratable.reserve(n);
    tpindex.reserve(n);

    cpuTime.length() = n;
    wallTime.length() = n;
    migratable.length() = n;
    tpindex.length() = n;

    copy(stats);
  }

  void copy(BaseLB::LDStats *stats){
    merge(stats,0.0,false);
  }

  void merge(BaseLB::LDStats *stats, double alpha, bool _m){
    for(int i = 0; i < stats->n_objs; i++){
      if(_m){
        cpuTime[i] = alpha*cpuTime[i] + (1.0-alpha)*stats->objData[i].cpuTime;
        wallTime[i] = alpha*wallTime[i] + (1.0-alpha)*stats->objData[i].cpuTime;
      }
      else{
        cpuTime[i] = stats->objData[i].cpuTime;
        wallTime[i] = stats->objData[i].cpuTime;
      }
      migratable[i] = stats->objData[i].migratable;
    }
    n_objs = stats->n_objs;
    count = stats->count;
    n_migrateobjs = stats->n_migrateobjs;
  }

  void pup(PUP::er &p){
    p | n_objs;
    p | n_migrateobjs;
    p | count;

    p | cpuTime;
    p | wallTime;
    p | migratable;
    p | tpindex;
  }
  
};

class MultistepLB : public CentralLB {
private:
  bool haveTPCentroids;
  ComparatorFn compares[NDIMS];
  ComparatorFn pc[NDIMS];

  TaggedVector3D *tpCentroids;
  CkReductionMsg *tpmsg;
  CkVec<int> *mapping;

  int procsPerNode;

  CkVec<CondensedLDStats> savedPhaseStats;       // stats saved from previous phases
  
  CmiBool QueryBalanceNow(int step);
  //int prevPhase;

  unsigned int determinePhase(unsigned int activeRung);
  void makeActiveProcessorList(CondensedLDStats *stats, int numActiveObjs, int numTotalObjs);
  void mergeInstrumentedData(int phase, BaseLB::LDStats *phaseStats);
  bool havePhaseData(int phase); 
  void printData(CondensedLDStats *stats, int phase);
public:
  MultistepLB(const CkLBOptions &);
  MultistepLB(CkMigrateMessage *m):CentralLB(m) { 
    lbname = "MultistepLB"; 
  }

  ~MultistepLB(){
    if(haveTPCentroids){
      haveTPCentroids = false;
      delete tpmsg;
      tpmsg = NULL;
    }
  }

  void work(BaseLB::LDStats* stats, int count);
  void receiveCentroids(CkReductionMsg *msg);
  //ScaleTranMapBG map;
  //
  void pup(PUP::er &p){
    CentralLB::pup(p);

    p | haveTPCentroids;
    p | savedPhaseStats;

    if(p.isUnpacking()){
      compares[0] = comparx;
      compares[1] = compary;
      compares[2] = comparz;

      pc[0] = pcx;
      pc[1] = pcy;
      pc[2] = pcz;

      lbname = "MultistepLB"; 
    }


  }

public:/* <- Sun CC demands Partition be public for ComputeLoad to access it. */

  class Partition {
  public:
    int refno;
    double load;			// total load in this set
    int origin[3];			// box coordinates
    int corner[3];
    int  count;				// number of objects in this partition
    int node, mapped;
    CkVec<int>   bkpes;			// background processors
  public:
    Partition(): refno(0), load(0.0), node(-1), mapped(0) {};
  };

private:  
  struct ComputeLoad {
    int id;
    int v[3];
    double load;
    int  refno;
    double  tv;
    Partition * partition;
  };
  
  
  struct VecArray {
    int v;
    int id;
  };
  
  enum {XDIR=0, YDIR, ZDIR};
  
  LDStats* statsData;
  int P;
  ComputeLoad *computeLoad;
  int nObjs;
  VecArray  *(vArray[3]);
  Partition *partitions;
  Partition top_partition;
  int npartition;
  int currentp, refno;
  
  void strategy();
  void rec_divide(int, Partition&);
  void setVal(int x, int y, int z);
  int sort_partition(int x, int p, int r);
  void qsort1(int x, int p, int r);
  void quicksort(int x);
  void mapPartitionsToNodes();

public:
  double overLoad;
  
//**************************************
// ORB3DLB functions
//**************************************
//
  void work2(CondensedLDStats *stats, float *ratios, int count, CkVec<int> *mapping, bool useRatios, int phase);
  void directMap(TPObject *tp, int ntp, Node *nodes);
  void map(TPObject *tp, int ntp, int nn, Node *procs, int xs, int ys, int zs, int dim);
  int nextDim(int dim, int xs, int ys, int zs);
  TPObject *partitionEvenLoad(TPObject *tp, int &ntp);
  Node *halveNodes(Node *start, int np);
};

#endif /* _MultistepLB */

/*@}*/
