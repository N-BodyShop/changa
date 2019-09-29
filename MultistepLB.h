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

#include "MultistepLB.decl.h"

void CreateMultistepLB();
BaseLB * AllocateMultistepLB();

//**************************************

class MultistepLB : public CBase_MultistepLB {
private:
  ComparatorFn compares[NDIMS];
  ComparatorFn pc[NDIMS];

  CkVec<int> *mapping;

  int procsPerNode;

  void init();
  bool QueryBalanceNow(int step);

  void makeActiveProcessorList(BaseLB::LDStats *stats, int numActiveObjs);
public:
  MultistepLB(const CkLBOptions &);
  MultistepLB(CkMigrateMessage *m) : CBase_MultistepLB(m) {
    init();
  }

  void work(BaseLB::LDStats* stats);

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

public:
  double overLoad;

//**************************************
// ORB3DLB functions
//**************************************
//
  void work2(BaseLB::LDStats* stats, int count);
  void greedy(BaseLB::LDStats* stats, int count);
  void directMap(TPObject *tp, int ntp, Node *nodes);
  void map(TPObject *tp, int ntp, int nn, Node *procs, int xs, int ys, int zs, int dim);
  int nextDim(int dim, int xs, int ys, int zs);
  TPObject *partitionEvenLoad(TPObject *tp, int &ntp);
  Node *halveNodes(Node *start, int np);

  void pup(PUP::er &p);
};

#endif /* _MultistepLB */

/*@}*/
