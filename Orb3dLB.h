/**
 * \addtogroup CkLdb
*/
/*@{*/

#ifndef _ORB3DLB_H_
#define _ORB3DLB_H_

#include "CentralLB.h"
#include "Orb3dLB.decl.h"
#include "TaggedVector3D.h"

void CreateOrb3dLB();
BaseLB * AllocateOrb3dLB();

#define NDIMS 3
class Centroid3d{
  public:
  float x;
  float y;
  float z;

  float *pointers[NDIMS];

  Centroid3d(){
    pointers[0] = &x;
    pointers[1] = &y;
    pointers[2] = &z;
  }

  float& operator[](int i){
    return *(pointers[i]);
  }

};

class TPObject{
  public:

  Vector3D<float> centroid;
  float load;
  int index;
  int lbindex;
};

class Proc {
  public:
  int x, y, z, t;
  int rank;
};

typedef int (*ComparatorFn) (const void *, const void *);

class Orb3dLB : public CentralLB {
private:
  CmiBool firstRound; 
  CmiBool centroidsAllocated;
  ComparatorFn compares[NDIMS];
  ComparatorFn pc[NDIMS];
  // pointer to stats->to_proc
  CkVec<int> *mapping;

  int procsPerNode;

  // things are stored in here before work
  // is ever called.
  CkVec<TaggedVector3D> tpCentroids;

  CmiBool QueryBalanceNow(int step);

  void printData(BaseLB::LDStats &stats, int phase, int *revObjMap);

public:
  Orb3dLB(const CkLBOptions &);
  Orb3dLB(CkMigrateMessage *m):CentralLB(m) { lbname = "Orb3dLB"; }
  void work(BaseLB::LDStats* stats, int count);
  void receiveCentroids(CkReductionMsg *msg);
  void directMap(TPObject *tp, int ntp, Proc *proc);
  void map(TPObject *tp, int ntp, int np, Proc *procs, int dim);
  int nextDim(int dim);
  TPObject *partitionEvenLoad(TPObject *tp, int &ntp);
  Proc *halveProcessors(Proc *start, int np);

};

#endif /* _ORB3DLB_H_ */

/*@}*/
