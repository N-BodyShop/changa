/**
 * Author: gplkrsh2@illinois.edu (Harshitha Menon)
*/

#ifndef HIERARCH_ORBLB_H
#define HIERARCH_ORBLB_H

#include "HierarchOrbLB.decl.h"
#include "HybridBaseLB.h"

#include "CentralLB.h"

class HierarchOrbLB : public CBase_HierarchOrbLB
{
public:
  HierarchOrbLB(const CkLBOptions &);
  HierarchOrbLB(CkMigrateMessage *m): CBase_HierarchOrbLB(m) {init();}
  ~HierarchOrbLB();

protected:
  CentralLB *orblb;

  void init();
  virtual bool QueryBalanceNow(int step) {
    if(step == 0) return false;
    return true;
  };
  virtual bool
  QueryMigrateStep(int) { return true; };
  virtual void work(LDStats* stats);
  virtual CLBStatsMsg* AssembleStats();

private:
  CProxy_HierarchOrbLB  thisProxy;
  void refine(LDStats* stats);
};

#endif /* HIERARCH_ORBLB_H */
