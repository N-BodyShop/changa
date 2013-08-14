/**
 * \addtogroup CkLdb
*/
/*@{*/

#ifndef _ORB3DLB_NOTOPO_H_
#define _ORB3DLB_NOTOPO_H_

#include "Orb3dLB_notopo.decl.h"
#include <queue>
#include "Orb3dLBCommon.h"

#include "pup_stl.h"

void CreateOrb3dLB_notopo();
BaseLB * AllocateOrb3dLB_notopo();

class Orb3dLB_notopo : public CentralLB, public Orb3dCommon {
private:

  vector<OrbObject> tps;
  // things are stored in here before work
  // is ever called.
  TaggedVector3D *tpCentroids;
  CkReductionMsg *tpmsg;
 
  bool QueryBalanceNow(int step);
  void printData(BaseLB::LDStats &stats, int phase, int *revObjMap);

  void pupDump(PUP::er &, BaseLB::LDStats *, vector<Event> *);

//  void orbPartition(vector<Event> *events, OrientedBox<float> &box, int procs);
//  int partitionRatioLoad(vector<Event> &events, float ratio);

public:
  Orb3dLB_notopo(const CkLBOptions &);
  Orb3dLB_notopo(CkMigrateMessage *m):CentralLB(m) { lbname = "Orb3dLB_notopo"; }
  void work(BaseLB::LDStats* stats);
  void receiveCentroids(CkReductionMsg *msg);

};

#endif /* _Orb3dLB_notopo_H_ */

/*@}*/
