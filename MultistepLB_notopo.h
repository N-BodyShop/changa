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

#ifndef _MULTISTEPLB_NOTOPO_H_
#define _MULTISTEPLB_NOTOPO_H_

#define MCLBMS          // multistepping enabled
#define MCLB_ORBSMOOTH  // orbsmooth for large steps
#define MCLB_RR         // round robin otherwise

#include "Orb3dLBCommon.h"

#include "MultistepLB_notopo.decl.h"

void CreateMultistepLB_notopo();
BaseLB * AllocateMultistepLB_notopo();

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
class MultistepLB_notopo : public CBase_MultistepLB_notopo, public Orb3dCommon {
private:
  bool firstRound;
  // things are stored in here before work
  // is ever called.
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
  MultistepLB_notopo(const CkLBOptions &);
  MultistepLB_notopo(CkMigrateMessage *m):CBase_MultistepLB_notopo(m) { 
    lbname = "MultistepLB_notopo"; 
     }
    
  void work(BaseLB::LDStats* stats);
  void receiveCentroids(CkReductionMsg *msg);
  //ScaleTranMapBG map;

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
  void work2(BaseLB::LDStats* stats, int count, int phase, int prevPhase);
  void greedy(BaseLB::LDStats* stats, int count, int phase, int prevPhase);
 
  void pup(PUP::er &p);
};

#endif /* _MultistepLB_notopo */

/*@}*/
