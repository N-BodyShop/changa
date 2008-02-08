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

#ifndef _METISCOSMOLB_H_
#define _METISCOSMOLB_H_

#include "CentralLB.h"

#include "MapStructures.h"
#include "ScaleTranMapBG.h"
#include "ScaledORBMapBG.h"

#include "MetisCosmoLB.decl.h"

void CreateMetisCosmoLB();
BaseLB * AllocateMetisCosmoLB();

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


class MetisCosmoLB : public CentralLB {
private:
  CmiBool firstRound; 
  CmiBool centroidsAllocated;

  CkVec<BaseLB::LDStats> savedPhaseStats;       // stats saved from previous phases
  
  CmiBool QueryBalanceNow(int step);
  //int prevPhase;

  unsigned int determinePhase(unsigned int activeRung);
  void makeActiveProcessorList(BaseLB::LDStats *stats, int numActiveObjs);
  void mergeInstrumentedData(int phase, BaseLB::LDStats *phaseStats);
  bool havePhaseData(int phase); 
  void printData(BaseLB::LDStats &stats, int phase, int *revObjMap);
public:
  MetisCosmoLB(const CkLBOptions &);
  MetisCosmoLB(CkMigrateMessage *m):CentralLB(m) { lbname = "MetisCosmoLB"; }
  void work(BaseLB::LDStats* stats, int count);
  void receiveCentroids(CkReductionMsg *msg);
  //ScaleTranMapBG map;
  ScaledORBMapBG map;

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
  void qsort(int x, int p, int r);
  void quicksort(int x);
  void mapPartitionsToNodes();

public:
  double overLoad;
  
};

#endif /* _METISCOSMOLB */

/*@}*/
