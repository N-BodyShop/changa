#include "ParallelGravity.h"

#include "IntraNodeLBManager.h"

CpvDeclare(int, tpGravDone);

IntraNodeLBManager::IntraNodeLBManager(CkMigrateMessage *m) :
  CBase_IntraNodeLBManager(m) {
  total_tps_ = 0;
  tps_done_ = 0;
  CpvInitialize(int, tpGravDone );
  CpvAccess(tpGravDone) = 1;
}

IntraNodeLBManager::IntraNodeLBManager(int dummy, CkGroupID gid) : 
    total_tps_(0), tps_done_(0), num_loc_mgr_(1) {
  loc_mgr_ = new CkGroupID[1];
  loc_mgr_[0] = gid;
  CpvInitialize(int, tpGravDone );
  CpvAccess(tpGravDone) = 1;
}


void IntraNodeLBManager::registerTP() {
  if (total_tps_ > 0) {
    return;
  }

  tps_done_ = 0;
  CkCacheArrayCounter local_chares;
  for (int i = 0; i < num_loc_mgr_; ++i) {
    CkLocMgr *mgr = (CkLocMgr *)CkLocalBranch(loc_mgr_[i]);
    mgr->iterate(local_chares);
  }
  total_tps_ = local_chares.count;
  CpvAccess(tpGravDone) = 0;
}

void IntraNodeLBManager::finishedTPWork() {
  tps_done_++;

  // All the work for the TPs in this PE is done so reset the counters for the
  // next iteration
  if (tps_done_ == total_tps_) {
    total_tps_ = 0;
    tps_done_ = 0;
    CpvAccess(tpGravDone) = 1;
  }
}

vector<int> IntraNodeLBManager::getOtherIdlePes() {
  vector<int> idlepes;
  int nsize = CkNodeSize(CkMyNode());
  for (int i = 0; i < nsize; i++) {
    if (CpvAccessOther(tpGravDone, i) == 1 &&
        (CkNodeFirst(CkMyNode()) + i)!=CkMyPe()) {
      idlepes.push_back(CkNodeFirst(CkMyNode()) + i);
    }
  }

  return idlepes;
}

void IntraNodeLBManager::pup(PUP::er &p) {
  CBase_IntraNodeLBManager::pup(p);
  p | num_loc_mgr_;
  if (p.isUnpacking()) {
    loc_mgr_ = new CkGroupID[num_loc_mgr_];
  }
  PUParray(p, loc_mgr_, num_loc_mgr_);
}
