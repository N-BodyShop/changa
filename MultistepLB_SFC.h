#ifndef _MULTISTEPLB_SFC_H_
#define _MULTISTEPLB_SFC_H_

#include "MultistepLB_SFC.decl.h"
#include "Vector3D.h"
#include "cosmoType.h"
#include "SFC.h"
#include "CentralLB.h"

void Orb_PrintLBStats(BaseLB::LDStats *stats, int numobjs);

/// @brief Multistep load balancer using Space Filling Curve
///
/// This balancer recognizes different "phases" (called rungs in other
/// parts of the code), and uses loads based on measurements of the
/// previous calculation at the same phase.  For large phases, (i.e.,
/// when many particles are active, the TreePieces are divided among
/// the processors using a Space Filling Curve based on the centroids
/// of the TreePieces.
///

class MultistepLB_SFC : public CBase_MultistepLB_SFC {
private:
    void init();
    bool QueryBalanceNow(int step);

    decltype(BaseLB::LDStats::to_proc) *mapping;
    decltype(BaseLB::LDStats::from_proc) *from;
    /// total computational cost to be balanced
    double dTotalLoad;
    /// Maximum number of pieces per processor
    double maxPieceProc;

public:
    MultistepLB_SFC(const CkLBOptions &);
    MultistepLB_SFC(CkMigrateMessage *m) : CBase_MultistepLB_SFC(m) {
        init();
    }

    class SFCObject
    {
    public:
        /// index into LB stats->objData
        int lbindex;
        /// Spacial location of TreePiece
        Vector3D<cosmoType> centroid;
        SFC::Key key;
        /// computational cost of this object
        double load;

        SFCObject() : lbindex(-1), load(0) {}
        SFCObject(int _lbindex, double _load) :
            lbindex(_lbindex),
            load(_load)
        {
        }
        bool operator<(const SFCObject &o) const{
            return key < o.key;
        }
    };

    void work(BaseLB::LDStats* stats);
    void work2(BaseLB::LDStats* stats);
    void sfcPrepare(std::vector<SFCObject> &tp_array,
                    OrientedBox<float> &boundingBox,
                    int nObjs, BaseLB::LDStats * stats,
                    bool node_partition=false);
    void sfcPartition(int nProcs, std::vector<SFCObject> & tp,
                      BaseLB::LDStats *stats, bool node_partition=false);
    void pup(PUP::er &p);
};


#endif /* _MultistepLB_notopo */
