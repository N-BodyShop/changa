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

//**************************************
// ORB3DLB functions
//**************************************


class LightweightLDStats {
  public:
  int n_objs;
  int n_migrateobjs;
  CkVec<LDObjData> objData;

  void pup(PUP::er &p);
};


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
  bool firstRound;
  // things are stored in here before work
  // is ever called.
  bool haveTPCentroids;
  int nrecvd;
  TaggedVector3D *tpCentroids;
  CkReductionMsg *tpmsg;

 // CkVec<OrbObject> tps;
  int procsPerNode;

  CkVec<LightweightLDStats> savedPhaseStats;      /// stats saved from previous phases
  
  bool QueryBalanceNow(int step);
  //int prevPhase;

  unsigned int determinePhase(unsigned int activeRung);
  void makeActiveProcessorList(BaseLB::LDStats *stats, int numActiveObjs);
  void mergeInstrumentedData(int phase, BaseLB::LDStats *phaseStats);
  bool havePhaseData(int phase); 
  void printData(BaseLB::LDStats &stats, int phase, int *revObjMap);


//  void orbPartition(CkVec<Event> *events, OrientedBox<float> &box, int procs, OrbObject * tp);
//  int partitionRatioLoad(CkVec<Event> &events, float ratio);


public:
  MultistepOrbLB(const CkLBOptions &);
  MultistepOrbLB(CkMigrateMessage *m):OrbLB(m) { 
    lbname = "MultistepOrbLB"; 
     }
    
  void work(BaseLB::LDStats* stats);
  void receiveCentroids(CkReductionMsg *msg);

public:
  
  void pup(PUP::er &p);
};

#endif /* _MultistepOrbLB */

/*@}*/
