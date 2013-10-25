/**
 * Author: gplkrsh2@illinois.edu (Harshitha Menon)
*/

#ifndef HIERARCH_ORBLB_H
#define HIERARCH_ORBLB_H

#include "HierarchOrbLB.decl.h"
#include "HybridBaseLB.h"

#include "CentralLB.h"

void CreateHierarchOrbLB();

class HierarchOrbLB : public HybridBaseLB
{
public:
  HierarchOrbLB(const CkLBOptions &);
  HierarchOrbLB(CkMigrateMessage *m): HybridBaseLB(m) {init();}
  ~HierarchOrbLB();

protected:
  CentralLB *orblb;
  CentralLB *refinelb;

  void init();
  virtual bool QueryBalanceNow(int step) {
    if(step == 0) return false;
    return true;
  };  
  virtual bool QueryMigrateStep(int) { return true; };  
  virtual void work(LDStats* stats);
  virtual void GetObjsToMigrate(int toPe, double load, LDStats *stats, 
      int atlevel, CkVec<LDCommData>& comms, CkVec<LDObjData>& objs);

private:
  CProxy_HierarchOrbLB  thisProxy;

};

#endif /* HIERARCH_ORBLB_H */
