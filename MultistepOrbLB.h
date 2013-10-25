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

#ifndef _MULTISTEPORBLB_H_
#define _MULTISTEPORBLB_H_

#define MCLBMS          // multistepping enabled

#include "MultistepOrbLB.decl.h"

#include "TaggedVector3D.h"
#include "OrbLB.h"

void CreateMultistepOrbLB();
BaseLB * AllocateMultistepOrbLB();

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
class MultistepOrbLB : public OrbLB {
private:
  bool QueryBalanceNow(int step);
  void init();

  unsigned int determinePhase(unsigned int activeRung);
  void makeActiveProcessorList(BaseLB::LDStats *stats, int numActiveObjs);
  void mergeInstrumentedData(int phase, BaseLB::LDStats *phaseStats);
  bool havePhaseData(int phase);
  void printData(BaseLB::LDStats &stats, int phase, int *revObjMap);

public:
  MultistepOrbLB(const CkLBOptions &);
  MultistepOrbLB(CkMigrateMessage *m):OrbLB(m) {
    init();
  }

  void work(BaseLB::LDStats* stats);

public:

  void pup(PUP::er &p);
};

#endif /* _MultistepOrbLB */

/*@}*/
