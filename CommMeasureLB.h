/**
 * \addtogroup CkLdb
*/
/*@{*/

#ifndef _CommMeasureLB_H_
#define _CommMeasureLB_H_

#include "CentralLB.h"
#include "CommMeasureLB.decl.h"
#include "TaggedVector3D.h"
#include <queue>

void CreateCommMeasureLB();
BaseLB * AllocateCommMeasureLB();

#define NDIMS 3
typedef int (*ComparatorFn) (const void *, const void *);

class CommMeasureLB : public CentralLB {
private:
  CmiBool firstRound; 
  CmiBool centroidsAllocated;
  // things are stored in here before work
  // is ever called.
  CkVec<TaggedVector3D> tpCentroids;

  CmiBool QueryBalanceNow(int step);

  void printData(BaseLB::LDStats &stats, int phase, int *revObjMap);

public:
  CommMeasureLB(const CkLBOptions &);
  CommMeasureLB(CkMigrateMessage *m):CentralLB(m) { lbname = "CommMeasureLB"; }
  void work(BaseLB::LDStats* stats, int count);
  void receiveCentroids(CkReductionMsg *msg);
};

#endif /* _CommMeasureLB_H_ */

/*@}*/
