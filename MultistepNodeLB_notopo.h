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

#ifndef _MULTISTEPNODELB_NOTOPO_H_
#define _MULTISTEPNODELB_NOTOPO_H_

#define MCLBMS          // multistepping enabled
#define MCLB_ORBSMOOTH  // orbsmooth for large steps
#define MCLB_RR         // round robin otherwise

#include "Orb3dLBCommon.h"

#include "MultistepNodeLB_notopo.decl.h"

void CreateMultistepLB_notopo();
BaseLB * AllocateMultistepNodeLB_notopo();


/// @brief Multistep load balancer where no processor topology
/// information is used.
///
/// This balancer recognizes different "phases" (called rungs in other
/// parts of the code), and uses loads based on measurements of the
/// previous calculation at the same phase.  For large phases, (i.e.,
/// when many particles are active, the TreePieces are divided among
/// by 3 dimensional ORB based on the centroids of the TreePieces.
/// For small phases, a greedy algorithm is used.
///
class MultistepNodeLB_notopo : public CentralLB, public Orb3dCommon {
private:

  int prevMaxPredPe;
  //int* existing_tps_on_pe;

  void init();
  bool QueryBalanceNow(int step);
  void makeActiveProcessorList(BaseLB::LDStats *stats, int numActiveObjs);


public:
  MultistepNodeLB_notopo(const CkLBOptions &);
  MultistepNodeLB_notopo(CkMigrateMessage *m):CentralLB(m) { 
    init();
  }
    
  void work(BaseLB::LDStats* stats);
  void balanceTPs(BaseLB::LDStats* stats);
  void balanceTPsNode(BaseLB::LDStats* stats);

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
  
//**************************************
// ORB3DLB functions
//**************************************
//
  void work2(BaseLB::LDStats* stats, int count);
  void greedy(BaseLB::LDStats* stats, int count);
  //void GetLodInfo(BaseLB::LDStats* stats);
  //void GetMaxNodeInfo(BaseLB::LDStats* stats);
 
  void pup(PUP::er &p);

  int maxPredLdSavedObj;
  double maxPredLdSaved;
};

#endif /* _MultistepNodeLB_notopo */

/*@}*/
