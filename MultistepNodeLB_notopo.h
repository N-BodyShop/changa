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
/// information is used. This first performs orb partition at the node level and
/// assigns the partitioned TreePieces to the PEs belonging to the node. Finally
/// after this assignment, a refinement is done.
class MultistepNodeLB_notopo : public CBase_MultistepNodeLB_notopo, public Orb3dCommon {
private:

  int prevMaxPredPe;

  void init();
  bool QueryBalanceNow(int step);

public:
  MultistepNodeLB_notopo(const CkLBOptions &);
  MultistepNodeLB_notopo(CkMigrateMessage *m) : CBase_MultistepNodeLB_notopo(m) {
    init();
  }
    
  void work(BaseLB::LDStats* stats);
  void balanceTPs(BaseLB::LDStats* stats);
  void balanceTPsNode(BaseLB::LDStats* stats);

private:  
  
  enum {XDIR=0, YDIR, ZDIR};
  
public:
  
//**************************************
// ORB3DLB functions
//**************************************
//
  void work2(BaseLB::LDStats* stats, int count);
  void greedy(BaseLB::LDStats* stats, int count);
 
  void pup(PUP::er &p);
};

#endif /* _MultistepNodeLB_notopo */

/*@}*/
