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

  int n_objs;
  int n_migrateobjs; 
  int count;

  public: 
  void init(BaseLB::LDStats *stats, int *lbToTp){
    int n = stats->n_objs;

    cpuTime.reserve(n);
    wallTime.reserve(n);
    migratable.reserve(n);

    cpuTime.length() = n;
    wallTime.length() = n;
    migratable.length() = n;

    copy(stats, lbToTp);
  }

  void copy(BaseLB::LDStats *stats, int *lbToTp){
    merge(stats,lbToTp,0.0,false);
  }

  void merge(BaseLB::LDStats *stats, int *lbToTp, double alpha, bool _m){
    for(int i = 0; i < stats->n_objs; i++){
      int tpindex = lbToTp[i];
      if(_m){
        cpuTime[tpindex] = alpha*cpuTime[tpindex] 
                           + (1.0-alpha)*stats->objData[i].cpuTime;
        wallTime[tpindex] = alpha*wallTime[tpindex] 
                           + (1.0-alpha)*stats->objData[i].cpuTime;
      }
      else{
        cpuTime[tpindex] = stats->objData[i].cpuTime;
        wallTime[tpindex] = stats->objData[i].cpuTime;
      }
      //migratable[tpindex] = stats->objData[i].migratable;
    }
    n_objs = stats->n_objs;
    count = stats->count;
    //n_migrateobjs = stats->n_migrateobjs;
  }

  void pup(PUP::er &p){
    p | n_objs;
    p | n_migrateobjs;
    p | count;

    p | cpuTime;
    p | wallTime;
    p | migratable;
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
  void makeActiveProcessorList(CondensedLDStats *stats, float ratio, bool largePhase);
  void mergeInstrumentedData(int phase, BaseLB::LDStats *phaseStats, int *lbToTp);
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

public:
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

  void greedySmallPhase(SmallPhaseObject *tp, int ntp, int np, CkVec<int> *mapping);
};

#endif /* _MultistepLB */

/*@}*/
