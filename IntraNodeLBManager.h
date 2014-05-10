#ifndef INTRA_NODE_LB_MANAGER_H
#define INTRA_NODE_LB_MANAGER_H

#include <vector>

//#include "ParallelGravity.decl.h"

using std::vector;

class IntraNodeLBManager : public CBase_IntraNodeLBManager {

 public:
  IntraNodeLBManager(int dummy, CkGroupID gid);

  IntraNodeLBManager(CkMigrateMessage *m);

  void pup(PUP::er &p);
 
  void registerTP();
  void finishedTPWork();

  vector<int> getOtherIdlePes();

 private:
  int total_tps_;
  int tps_done_;
  int num_loc_mgr_;
  CkGroupID *loc_mgr_;
};
#endif
