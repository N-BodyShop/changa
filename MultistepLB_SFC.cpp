#include <charm++.h>
#include "MultistepLB_SFC.h"
#include "ParallelGravity.h"
#include "Vector3D.h"
#include "formatted_string.h"

CkpvExtern(int, _lb_obj_index);
using namespace std;

#if CHARM_VERSION > 61002
static void lbinit()
{
    LBRegisterBalancer<MultistepLB_SFC>("MultistepLB_SFC",
      "Works best with multistepped runs; uses SFC distribution");
}
#else
CreateLBFunc_Def(MultistepLB_SFC,
                 "Works best with multistepped runs; uses SFC distribution");
#endif

void MultistepLB_SFC::init() {
    lbname = "MultistepLB_SFC";
    if (CkpvAccess(_lb_obj_index) == -1)
        CkpvAccess(_lb_obj_index) = LBRegisterObjUserData(sizeof(TaggedVector3D));
}


MultistepLB_SFC::MultistepLB_SFC(const CkLBOptions &opt): CBase_MultistepLB_SFC(opt)
{
    init();
    if (CkMyPe() == 0){
        CkPrintf("[%d] MultistepLB_SFC created\n",CkMyPe());
    }
}

bool MultistepLB_SFC::QueryBalanceNow(int step){
    if(CkMyPe() == 0) CkPrintf("LB_SFC: Step %d\n", step);
    return true;
}

/// @brief Implement load balancing: store loads and determine active
/// processors and objects, sort by SFC, then divide up among processors.
/// @param stats The Load Balancer statistics object.
void MultistepLB_SFC::work(BaseLB::LDStats* stats)
{
#if CMK_LBDB_ON
    // find active objects - mark the inactive ones as non-migratable
    const auto num_objs = stats->objData.size();

    if(_lb_args.debug() >= 2 && step() > 0) {
        // Write out "particle file" of measured load balance information
        auto achFileName = make_formatted_string("lb_a.%d.sim", step()-1);
        write_LB_particles(stats, achFileName.c_str(), true);
    }

    int numActiveObjects = 0;
    int numInactiveObjects = 0;
    int minActiveProc = INT_MAX;
    int maxActiveProc = 0;

    for(int i = 0; i < num_objs; i++){
        stats->to_proc[i] = stats->from_proc[i];
    }

    for(int i = 0; i < num_objs; i++){
        if (!stats->objData[i].migratable) continue;

        LDObjData &odata = stats->objData[i];
        TaggedVector3D* udata = (TaggedVector3D *)odata.getUserData(CkpvAccess(_lb_obj_index));

        if(udata->myNumParticles == 0){ // ignore pieces with no particles
            stats->objData[i].migratable = 0;
            stats->n_migrateobjs--;
            continue;
        }
        if(udata->numActiveParticles == 0){
            numInactiveObjects++;
        }
        else{
            numActiveObjects++;
            if(minActiveProc > stats->from_proc[i])
                minActiveProc = stats->from_proc[i];
            if(maxActiveProc < stats->from_proc[i])
                maxActiveProc = stats->from_proc[i];
        }
    }
    CkPrintf("numActiveObjects: %d, numInactiveObjects: %d\n", numActiveObjects,
             numInactiveObjects);
    CkPrintf("active PROC range: %d to %d\n", minActiveProc, maxActiveProc);
    if(numActiveObjects < 0.1*numInactiveObjects) {
        // only a small number of active objects, only migrate them
        for(int i = 0; i < stats->objData.size(); i++){
            if (!stats->objData[i].migratable) continue;

            LDObjData &odata = stats->objData[i];
            TaggedVector3D* udata =
              (TaggedVector3D *)odata.getUserData(CkpvAccess(_lb_obj_index));
            if(udata->numActiveParticles == 0) {
                stats->objData[i].migratable = 0;
                stats->n_migrateobjs--;
            }
        }
    }
    else {
        CkPrintf("Migrating all: numActiveObjects: %d, numInactiveObjects: %d\n",
                 numActiveObjects, numInactiveObjects);
    }

    // let the strategy take over on this modified instrumented data and processor information
    work2(stats);
#endif //CMK_LDB_ON
}

/// @brief SFC load balance.
void MultistepLB_SFC::work2(BaseLB::LDStats *stats){
    const int numobjs = stats->objData.size();
    const int nmig = stats->n_migrateobjs;

    // this data structure is used by the SFC strategy
    // to balance objects. it is NOT indexed by tree piece index
    // there are as many entries in it as there are
    // migratable (active) tree pieces
    vector<SFCObject> tp_array;
    tp_array.resize(nmig);

    if (_lb_args.debug()>=2) {
        CkPrintf("[work2] ready tp_array data structure\n");
    }

    int numProcessed = 0;

    dTotalLoad = 0.0;
    for(int i = 0; i < numobjs; i++){
        if(!stats->objData[i].migratable) continue;

        float load;
        LDObjData &odata = stats->objData[i];
        TaggedVector3D* udata = (TaggedVector3D *)odata.getUserData(CkpvAccess(_lb_obj_index));
        if(step() == 0){ // no load information, balance by particle numbers
            load = udata->myNumParticles;
        }
        else{
            load = stats->objData[i].wallTime;
        }

        tp_array[numProcessed] = SFCObject(i, load);
        tp_array[numProcessed].centroid = udata->vec;
        numProcessed++;
        dTotalLoad += load;
    }
    CkAssert(numProcessed==nmig);

    sfcPrepare(tp_array, nmig, stats);
    sfcPartition(stats->nprocs(),tp_array, stats);

    // refine(stats, numobjs);
    Orb_PrintLBStats(stats, numobjs);

    if(_lb_args.debug() >= 2) {
        // Write out "particle file" of load balance information
        auto achFileName = make_formatted_string("lb.%d.sim", step());
        write_LB_particles(stats, achFileName.c_str(), false);
    }
}

/// @brief Prepare structures for the ORB partition.
/// @param tp_array Reference to Vector of Objects representing TreePieces.
/// @param nObjs Number of tree pieces to partition.
/// @param stats Data from the load balancing framework.
/// @param node_partition Are we partitioning on nodes.
void MultistepLB_SFC::sfcPrepare(vector<SFCObject> &tp_array,
                                 int nObjs,
                                 BaseLB::LDStats *stats,
                                 bool node_partition){

    OrientedBox<float> boundingBox;
    int nmig = stats->n_migrateobjs;
    if(dMaxBalance < 1.0)
        dMaxBalance = 1.0;

    // If using node based orb partition, then the maxPieceProc is total
    // migratable objs / total number of node.
    if (node_partition) {
        maxPieceProc = dMaxBalance * nmig / CkNumNodes();
    } else {
        maxPieceProc = dMaxBalance*nmig/stats->nprocs();
    }

    if(maxPieceProc < 1.0)
        maxPieceProc = 1.01;

    CkAssert(tp_array.size() == nObjs);

    mapping = &stats->to_proc;
    from = &stats->from_proc;

    CkPrintf("[LB_SFC] sorting\n");
    for(int i = 0; i < nObjs; i++)
        boundingBox.grow(tp_array[i].centroid);

    // N.B. code below from TreePiece::assignKeys().
    // Refactoring is a possibility.
    // get longest axis
    Vector3D<float> bsize = boundingBox.size();
    float max = (bsize.x > bsize.y) ? bsize.x : bsize.y;
    max = (max > bsize.z) ? max : bsize.z;
    //
    // Make the bounding box cubical.
    //
    Vector3D<float> bcenter = boundingBox.center();
    // The magic number below is approximately 2^(-19)
    const float fEps = 1.0 + 1.91e-6;  // slop to ensure keys fall
                                       // between 0 and 1.
    bsize = Vector3D<float>(fEps*0.5*max);
    boundingBox = OrientedBox<float>(bcenter-bsize, bcenter+bsize);
    if(verbosity > 1)
        ckout << "TreePiece: Bounding box now: " << boundingBox << endl;

    for(unsigned int i = 0; i < nObjs; ++i) {
        tp_array[i].key = SFC::generateKey(tp_array[i].centroid, boundingBox);
    }
    sort(tp_array.begin(),tp_array.end());
}

/// @brief Partition treepieces among processors by
/// dividing the SFC as evenly as possible.
/// @param nprocs Number of processors over which to partition the
/// pieces. N.B. if node_partition is true, then this is the number of nodes.
/// @param tp Vector of TreePiece data.
/// @param stats Load balance data
void MultistepLB_SFC::sfcPartition(int nProcs, vector<SFCObject> & tp,
                                   BaseLB::LDStats *stats,
                                   bool node_partition){

    double loadPrev = 0.0; // load on all lower processors
    int iCurrPiece = 0;    // Piece under consideration
    for (int iProc = 0; iProc < nProcs; iProc++) {
        if (!stats->procs[iProc].available)
            continue;
        // always assign one piece to a processor
        SFCObject &oSFC = tp[iCurrPiece];
        (*mapping)[oSFC.lbindex] = iProc;
        double loadCurr = oSFC.load;
        iCurrPiece++;
        int nCurrPiece = 1;     // number of pieces on this processor
        double loadTarget = (iProc+1)*dTotalLoad/nProcs;
        double dLoadError = fabs(loadPrev + loadCurr - loadTarget);

        while ((nCurrPiece < maxPieceProc)
               && fabs(tp[iCurrPiece].load + loadPrev + loadCurr - loadTarget)
               < dLoadError) {  // add pieces to this processor to get
                                // the closest to the target load
            oSFC = tp[iCurrPiece];
            loadCurr += oSFC.load;
            (*mapping)[oSFC.lbindex] = iProc;
            dLoadError = fabs(loadPrev + loadCurr - loadTarget);
            iCurrPiece++;
            nCurrPiece++;
        }
        loadPrev += loadCurr;
    }
    CkAssert(iCurrPiece == tp.size());
}

void MultistepLB_SFC::pup(PUP::er &p){
    CBase_MultistepLB_SFC::pup(p);
}

#include "MultistepLB_SFC.def.h"
