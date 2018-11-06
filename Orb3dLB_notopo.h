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

/// @brief Load balancer that divides work according to a 3D spatial
/// ORB.
class Orb3dLB_notopo : public CBase_Orb3dLB_notopo, public Orb3dCommon {
private:

  vector<OrbObject> tps;
  // things are stored in here before work
  // is ever called.

  void init();
  bool QueryBalanceNow(int step);
  void printData(BaseLB::LDStats &stats, int phase, int *revObjMap);

  void pupDump(PUP::er &, BaseLB::LDStats *, vector<Event> *);

public:
  Orb3dLB_notopo(const CkLBOptions &);
  Orb3dLB_notopo(CkMigrateMessage *m): CBase_Orb3dLB_notopo(m) {init();}
  void work(BaseLB::LDStats* stats);

};

#endif /* _Orb3dLB_notopo_H_ */

/*@}*/
