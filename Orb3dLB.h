/**
 * \addtogroup CkLdb
*/
/*@{*/

#ifndef _ORB3DLB_H_
#define _ORB3DLB_H_

#include "CentralLB.h"
#include "MapStructures.h"
#include "Orb3dLB.decl.h"
#include "TaggedVector3D.h"
#include <queue>

void CreateOrb3dLB();
BaseLB * AllocateOrb3dLB();

class Orb3dLB : public CBase_Orb3dLB {
  friend class MultistepLB;
private:
  bool firstRound;
  bool centroidsAllocated;
  ComparatorFn compares[NDIMS];
  ComparatorFn pc[NDIMS];
  // pointer to stats->to_proc
  CkVec<int> *mapping;

  int procsPerNode;

  CkVec<TPObject> tps;
  CkVec<Node> nodes;

  void init();
  bool QueryBalanceNow(int step);

  void printData(BaseLB::LDStats &stats, int phase, int *revObjMap);

public:
  Orb3dLB(const CkLBOptions &);
  Orb3dLB(CkMigrateMessage *m) : CBase_Orb3dLB(m) {init();}
  void work(BaseLB::LDStats* stats);
  void directMap(int tpstart, int tpend, int nodestart, int nodeend);
  void map(int tpstart, int tpend, int nodestart, int nodeend, int xs, int ys, int zs, int dim);
  int nextDim(int dim, int xs, int ys, int zs);
  void partitionEvenLoad(int tpstart, int tpend, int &tpmid);
  void halveNodes(int nodestart, int nodeend, int &nodemid);


};

#endif /* _ORB3DLB_H_ */

/*@}*/
