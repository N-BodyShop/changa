#include <charm++.h>
#include "cklists.h"
#include "MultistepLB_notopo.h"
#include "TopoManager.h"
#include "ParallelGravity.h"
#include "Vector3D.h"
#include <queue>
#include "formatted_string.h"

extern CProxy_TreePiece treeProxy;
CkpvExtern(int, _lb_obj_index);
using namespace std;
//#define ORB3DLB_NOTOPO_DEBUG CkPrintf

#if CHARM_VERSION > 61002
static void lbinit()
{
    LBRegisterBalancer<MultistepLB_notopo>("MultistepLB_notopo",
      "Works best with multistepped runs; uses Orb3D_notopo");
}
#else
CreateLBFunc_Def(MultistepLB_notopo, "Works best with multistepped runs; uses Orb3D_notopo");
#endif

void MultistepLB_notopo::init() {
  lbname = "MultistepLB_notopo";
  if (CkpvAccess(_lb_obj_index) == -1)
    CkpvAccess(_lb_obj_index) = LBRegisterObjUserData(sizeof(TaggedVector3D));
}


MultistepLB_notopo::MultistepLB_notopo(const CkLBOptions &opt): CBase_MultistepLB_notopo(opt)
{
  init();
  if (CkMyPe() == 0){
    CkPrintf("[%d] MultistepLB_notopo created\n",CkMyPe());
  }
}

bool MultistepLB_notopo::QueryBalanceNow(int step){
 if(CkMyPe() == 0) CkPrintf("Orb3dLB_notopo: Step %d\n", step);
 //  if(step == 0) return false;
  return true;

}

// helper functions for multistepping
#ifdef MCLBMSV

void MultistepLB_notopo::makeActiveProcessorList(BaseLB::LDStats *stats, int numActiveObjs){
  int objsPerProc = 8;
  int expandFactor = 4;
  int procsNeeded;
  procsNeeded = expandFactor*numActiveObjs/objsPerProc > stats->nprocs() ? stats->nprocs() : expandFactor*numActiveObjs/objsPerProc;

  /* currently, only the first procsNeeded procs are used - could do something more sophisticated here in the future - FIXME */
  CkPrintf("Processors 0 to %d active\n", procsNeeded-1);
}
#endif

/// @brief Implement load balancing: store loads and determine active
/// processors and objects, then call ORB3D.
void MultistepLB_notopo::work(BaseLB::LDStats* stats)
{
#if CMK_LBDB_ON
  // find active objects - mark the inactive ones as non-migratable
  int count;
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
    CkPrintf("Migrating all: numActiveObjects: %d, numInactiveObjects: %d\n", numActiveObjects, numInactiveObjects);
  }

  // select processors
#ifdef MCLBMSV
  CkPrintf("making active processor list\n");
  makeActiveProcessorList(stats, numActiveObjects);
#endif
  count = stats->nprocs();

  // let the strategy take over on this modified instrumented data and processor information
  work2(stats,count);
#endif //CMK_LDB_ON

}

//**************************************
// ORB3DLB functions
//**************************************
//

/// @brief ORB3D load balance.
void MultistepLB_notopo::work2(BaseLB::LDStats *stats, int count){
  const int numobjs = stats->objData.size();
  int nmig = stats->n_migrateobjs;

  // this data structure is used by the orb3d strategy
  // to balance objects. it is NOT indexed by tree piece index
  // there are as many entries in it as there are
  // migratable (active) tree pieces
 vector<OrbObject> tp_array;
 tp_array.resize(nmig);

  if (_lb_args.debug()>=2) {
    CkPrintf("[work2] ready tp_array data structure\n");
  }

 vector<Event> tpEvents[NDIMS];
  for(int i = 0; i < NDIMS; i++){
    tpEvents[i].reserve(nmig);
  }

  OrientedBox<float> box;

  int numProcessed = 0;

  for(int i = 0; i < numobjs; i++){
    if(!stats->objData[i].migratable) continue;

    float load;
    LDObjData &odata = stats->objData[i];
    TaggedVector3D* udata = (TaggedVector3D *)odata.getUserData(CkpvAccess(_lb_obj_index));
    if(step() == 0){
      load = udata->myNumParticles;
    }
    else{
      load = stats->objData[i].wallTime;
    }

    tpEvents[XDIM].push_back(Event(udata->vec.x,load,numProcessed));
    tpEvents[YDIM].push_back(Event(udata->vec.y,load,numProcessed));
    tpEvents[ZDIM].push_back(Event(udata->vec.z,load,numProcessed));

    tp_array[numProcessed]= OrbObject(i, udata->myNumParticles);
    tp_array[numProcessed].centroid = udata->vec;
    numProcessed++;
  }
  CkAssert(numProcessed==nmig);

  orbPrepare(tpEvents, box, nmig, stats);
  orbPartition(tpEvents,box,stats->nprocs(),tp_array, stats);

  refine(stats, numobjs);

  if(_lb_args.debug() >= 2) {
      // Write out "particle file" of load balance information
      auto achFileName = make_formatted_string("lb.%d.sim", step());
      write_LB_particles(stats, achFileName.c_str(), false);
  }
}

void Orb_PrintLBStats(BaseLB::LDStats *stats, int numobjs)
{
    std::vector<double> predLoad(stats->nprocs(), 0.0);
    std::vector<int> predCount(stats->nprocs(), 0);
    double maxObjLoad = 0.0;

    int migr = 0;
    for(int i = 0; i < numobjs; i++){
        LDObjData &odata = stats->objData[i];
        TaggedVector3D* udata = (TaggedVector3D *)odata.getUserData(CkpvAccess(_lb_obj_index));
        if(udata->myNumParticles == 0) // ignore empty TreePieces
            continue;
        if(stats->to_proc[i] != stats->from_proc[i])
            migr++;
        double ld = stats->objData[i].wallTime;
        int proc = stats->to_proc[i];
        predLoad[proc] += ld;
        predCount[proc] ++;
        if(ld > maxObjLoad)
            maxObjLoad = ld;
    }

    double minWall = 0.0;
    double maxWall = 0.0;
    double avgWall = 0.0;

    double minIdle = 0.0;
    double maxIdle = 0.0;
    double avgIdle = 0.0;

    double minBg = 0.0;
    double maxBg = 0.0;
    double avgBg = 0.0;

    double avgPred = 0.0;
    double minPred = 0.0;
    double maxPred = 0.0;

    double avgPiece = 0.0;
    double minPiece = 0.0;
    double maxPiece = 0.0;

    CkPrintf("***************************\n");
    for(int i = 0; i < stats->nprocs(); i++){
        double wallTime = stats->procs[i].total_walltime;
        double idleTime = stats->procs[i].idletime;
        double bgTime = stats->procs[i].bg_walltime;
        double pred = predLoad[i];
        double npiece = predCount[i];

        avgWall += wallTime;
        avgIdle += idleTime;
        avgBg += bgTime;
        avgPred += pred;
        avgPiece += npiece;

        if(i==0 || minWall > wallTime) minWall = wallTime;
        if(i==0 || maxWall < wallTime) maxWall = wallTime;

        if(i==0 || minIdle > idleTime) minIdle = idleTime;
        if(i==0 || maxIdle < idleTime) maxIdle = idleTime;

        if(i==0 || minBg > bgTime) minBg = bgTime;
        if(i==0 || maxBg < bgTime) maxBg = bgTime;

        if(i==0 || minPred > pred) minPred = pred;
        if(i==0 || maxPred < pred) maxPred = pred;

        if(i==0 || minPiece > npiece) minPiece = npiece;
        if(i==0 || maxPiece < npiece) maxPiece = npiece;

    }

    avgWall /= stats->nprocs();
    avgIdle /= stats->nprocs();
    avgBg /= stats->nprocs();
    avgPred /= stats->nprocs();
    avgPiece /= stats->nprocs();

#ifdef PRINT_LOAD_PERCENTILES
    double accumVar = 0;
    vector<double> objectWallTimes;
    for(int i = 0; i < stats->nprocs(); i++){
        double wallTime = stats->procs[i].total_walltime;
        objectWallTimes.push_back(wallTime);
        accumVar += (wallTime - avgWall) * (wallTime - avgWall);
    }
    double stdDev = sqrt(accumVar / stats->nprocs());
    CkPrintf("Average load: %.3f\n", avgWall);
    CkPrintf("Standard deviation: %.3f\n", stdDev);

    std::sort(objectWallTimes.begin(), objectWallTimes.end());
    CkPrintf("Object load percentiles: \n");
    double increment = (double) objectWallTimes.size() / 10;
    int j = 0;
    double index = 0;
    for (int j = 0; j < 100; j += 10) {
        index += increment;
        CkPrintf("%d: %.3f\n", j, objectWallTimes[(int) index]);
    }
    CkPrintf("100: %.3f\n", objectWallTimes.back());
#endif

    CkPrintf("LB stats: maxObjLoad %f\n", maxObjLoad);
    CkPrintf("LB stats: minWall %f maxWall %f avgWall %f maxWall/avgWall %f\n", minWall, maxWall, avgWall, maxWall/avgWall);
    CkPrintf("LB stats: minIdle %f maxIdle %f avgIdle %f minIdle/avgIdle %f\n", minIdle, maxIdle, avgIdle, minIdle/avgIdle);
    CkPrintf("LB stats: minPred %f maxPred %f avgPred %f maxPred/avgPred %f\n", minPred, maxPred, avgPred, maxPred/avgPred);
    CkPrintf("LB stats: minPiece %f maxPiece %f avgPiece %f maxPiece/avgPiece %f\n", minPiece, maxPiece, avgPiece, maxPiece/avgPiece);
    CkPrintf("LB stats: minBg %f maxBg %f avgBg %f maxBg/avgBg %f\n", minBg, maxBg, avgBg, maxBg/avgBg);
    CkPrintf("LB stats: orb migrated %d objects\n", migr);
}

/// @brief  Write out TreePieces as "particles" into a simple file
/// that can be converted into a tipsy file for visualization and
/// analysis.
/// @param stats LB structure
/// @param achFileName file to write.
/// @param bFrom use "from" processor if true, otherwise, use "to" processor
void write_LB_particles(BaseLB::LDStats* stats, const char *achFileName, bool bFrom)
{
    const auto num_objs = stats->objData.size();
    FILE *fp = fopen(achFileName, "w");
    CkAssert(fp != NULL);

    int num_migratables = num_objs;
    for(int i = 0; i < num_objs; i++) {
        if (!stats->objData[i].migratable) {
            num_migratables--;
        }
    }

    fprintf(fp, "%d %d 0\n", num_migratables, num_migratables);
    for(int i = 0; i < num_objs; i++) {
        if (!stats->objData[i].migratable) continue;

        LDObjData &odata = stats->objData[i];
        TaggedVector3D* udata = (TaggedVector3D *)odata.getUserData(CkpvAccess(_lb_obj_index));
        int proc;
        if(bFrom) 
            proc = stats->from_proc[i];
        else
            proc = stats->to_proc[i];
        fprintf(fp, "%g %g %g %g 0.0 0.0 0.0 %d %d\n",
                stats->objData[i].wallTime,
                udata->vec.x, udata->vec.y, udata->vec.z,
                proc, udata->tp);
    }
    fclose(fp);
}

void MultistepLB_notopo::pup(PUP::er &p){
  CBase_MultistepLB_notopo::pup(p);
}

#include "MultistepLB_notopo.def.h"
