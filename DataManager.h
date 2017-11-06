/** @file DataManager.h
 The DataManager holds data that TreePieces need access to, but not
 their own copy.
 @author Graeme Lufkin (gwl@u.washington.edu)
*/

#ifndef DATAMANAGER_H
#define DATAMANAGER_H

#include <vector>
#include <map>
#include <string>
#include "GenericTreeNode.h"
#include "ParallelGravity.decl.h"

#if CHARM_VERSION > 60401 && CMK_BALANCED_INJECTION_API
#include "ckBIconfig.h"
#endif


struct TreePieceDescriptor{
	TreePiece *treePiece;
        Tree::GenericTreeNode *root;

	TreePieceDescriptor(){}
	TreePieceDescriptor(TreePiece *tp_, GenericTreeNode *r){
		treePiece = tp_;
                root = r;
	}
};

#ifdef CUDA
struct UpdateParticlesStruct{
  CkCallback *cb;
  DataManager *dm;
  VariablePartData *buf;
  int size;
};


// store pointers to flattened buffers if !gpuFree 
struct PendingBuffers {
  CkVec<CudaMultipoleMoments> *moments;
  CkVec<CompactPartData> *particles;
  int chunk;
  //CudaMultipoleMoments *moments;
  //int nMoments;
  //CompactPartData *particles;
  //int nParticles;
};

#endif

/** The DataManager is used to store information that all TreePieces will need,
 but will not modify.  The first example is the list of splitter keys and the
 responsible chare for each interval.  This data is given to the DataManager by
 the Sorter.  The DataManager then instructs all the TreePieces on its node
 to evaluate the boundary keys.
 */
class DataManager : public CBase_DataManager {
	friend class TreePiece;
        friend class OctTreeBuildPhaseIWorker;

	/// The array of TreePieces I hold data for.
	CProxy_TreePiece treePieces;

protected:

	/// The array of splitter keys for the sort.
	std::vector<SFC::Key> boundaryKeys;

	/// An array identifying which chare is responsible for each interval of keys.
	std::vector<int> responsibleIndex;

	/// An array with how many particles are held by each TreePiece when sorted.
	std::vector<int> particleCounts;

	/// A list of roots of the TreePieces in this node
	// holds chare array indices of registered treepieces
	CkVec<TreePieceDescriptor> registeredTreePieces;
#ifdef CUDA
	//CkVec<int> registeredTreePieceIndices;
        int cumNumReplicatedNodes;
        int treePiecesDone;
        int savedChunk;
        int treePiecesDonePrefetch;
        int treePiecesDoneLocalComputation;
        // XXX - assumes that only one chunk can be on the gpu
        // at a given time
        int treePiecesDoneRemoteChunkComputation;
        int treePiecesWantParticlesBack;
        int savedNumTotalParticles;
        int savedNumTotalNodes;
        // keeps track of buckets of particles that were
        // received during the prefetch and which were subsequently
        // shipped off to the gpu - XXX
        // not including cached particles in postfrefetch entities shipped to gpu
        // since it is hard to count their number given just the pointer to the cache entry
        // * either do not concern yourself with cached particles
        // * or for each entry, get key, find bucket node in CM, DM or TPs and get number
        // for now, the former

        std::map<NodeKey, int> cachedPartsOnGpu;
        // local particles that have been copied to the gpu
        //std::map<NodeKey, int> localPartsOnGpu;

        //std::map<Tree::NodeKey, GenericTreeNode *> missedNodesOnGpu;
        //std::map<Tree::NodeKey, ExternalGravityParticle *> missedPartsOnGpu;

        // can the gpu accept a chunk of remote particles/nodes?
        bool gpuFree;
        // queue that stores all pending chunk transfers
        CkQ<PendingBuffers *> pendingChunkTransferQ;

        // last remote chunk's size in moments and particles
        int lastChunkMoments;
        int lastChunkParticles;

#ifdef CUDA_INSTRUMENT_WRS
        int activeRung;
        int treePiecesDoneInitInstrumentation;
#endif
#endif
	/// The root of the combined trees
	Tree::GenericTreeNode * root;
	/// Table for nodes created by the DataManager to combine trees.
	/// Kept track of here so we can delete them when done.
	CkVec<GenericTreeNode *> nodeTable;

	/// Number of chunks in which the tree was splitted during last combine operation
	int oldNumChunks;
	/// Nodes currently used as roots for remote computation
	Tree::NodeKey *chunkRoots;
        /// Lookup table for the chunkRoots
        Tree::NodeLookupType chunkRootTable;

public:

	/* 
	 ** Cooling 
	 */
	COOL *Cool;
	/// @brief log of star formation events.
	///
	/// Star formation events are stored on the data manager since there
	/// is no need to migrate them with the TreePiece.
	StarLog *starLog;
    HMStarLog *hmStarLog;
	/// @brief Lock for accessing starlog from TreePieces
	CmiNodeLock lockStarLog;
	/// @brief Lock for accessing hmstarlog from TreePieces
    CmiNodeLock lockHMStarLog;

	DataManager(const CkArrayID& treePieceID);
	DataManager(CkMigrateMessage *);

#ifdef CUDA
        //void serializeNodes(GenericTreeNode *start, CudaMultipoleMoments *&postPrefetchMoments, CompactPartData *&postPrefetchParticles);
		//void serializeNodes(GenericTreeNode *start);
        void donePrefetch(int chunk); // serialize remote chunk wrapper
        void serializeLocalTree();

        // actual serialization methods
        PendingBuffers *serializeRemoteChunk(GenericTreeNode *);
	void serializeLocal(GenericTreeNode *);
        void freeLocalTreeMemory();
        void freeRemoteChunkMemory(int chunk);
        void transferParticleVarsBack();
        void updateParticles(UpdateParticlesStruct *data);
        void initiateNextChunkTransfer();
#ifdef CUDA_INSTRUMENT_WRS
        int initInstrumentation();
#endif
        DataManager(){} 
#endif
        void clearInstrument(CkCallback &cb);

private:
        void init();

public:

	~DataManager() {
	    for (unsigned int i = 0; i < nodeTable.length(); i++) {
      		delete nodeTable[i];
    		}
    	    nodeTable.clear();

	    CoolFinalize(Cool);
	    delete starLog;
        delete hmStarLog;
	    CmiDestroyLock(lockStarLog);
        CmiDestroyLock(lockHMStarLog);
	    }

	/// Called by ORB Sorter, save the list of which TreePiece is
	/// responsible for which interval.
	void acceptResponsibleIndex(const int* responsible, const int n,
				    const CkCallback& cb);
	/// Called by the Sorter, I save these final keys and the list
	/// of which TreePiece is responsible for which interval.
	/// This routine then calls TreePiece::unshuffleParticles to
	/// move the particles around.
	/// @param keys vector of boundary keys
	/// @param responsible vector of which piece is responsible
	/// for which interval
	/// @param bins number of particles in each interval.
	void acceptFinalKeys(const SFC::Key* keys, const int* responsible, uint64_t* bins, const int n, const CkCallback& cb);
	void pup(PUP::er& p);

#ifdef CUDA
        /*
        std::map<NodeKey, int> &getLocalPartsOnGpuTable(){
          return localPartsOnGpu;
        }
        */
        std::map<NodeKey, int> &getCachedPartsOnGpuTable(){
          return cachedPartsOnGpu;
        }
#endif
	// Functions used to create a tree inside the DataManager comprising
	// all the trees in the TreePieces in the local node
private:
	Tree::GenericTreeNode *buildProcessorTree(int n, Tree::GenericTreeNode **gtn);
	int createLookupRoots(Tree::GenericTreeNode *node, Tree::NodeKey *keys);
public:

/// \brief Collect roots of treepieces on this node.
///
/// The roots are stored in registeredChares to be used by TreePiece
/// combineLocalTrees.
    void notifyPresence(Tree::GenericTreeNode *root, TreePiece *treePiece);
    void combineLocalTrees(CkReductionMsg *msg);
    void getChunks(int &num, Tree::NodeKey *&roots);
    inline Tree::GenericTreeNode *chunkRootToNode(const Tree::NodeKey k) {
      NodeLookupType::iterator iter = chunkRootTable.find(k);
      if (iter != chunkRootTable.end()) return iter->second;
      else return NULL;
    }
    inline Tree::GenericTreeNode *getRoot() { return root; }
    void initCooling(double dGmPerCcUnit, double dComovingGmPerCcUnit,
		     double dErgPerGmUnit, double dSecUnit, double dKpcUnit,
		     COOLPARAM inParam, const CkCallback& cb);
    void initStarLog(std::string _fileName, const CkCallback &cb);
    void initHMStarLog(std::string _fileName, const CkCallback &cb);
    void dmCoolTableRead(double *dTableData, int nData, const CkCallback& cb);
    void CoolingSetTime(double z, // redshift
			double dTime, // Time
			const CkCallback& cb);
    void SetStarCM(double dCenterOfMass[4], const CkCallback& cb);
    void memoryStats(const CkCallback& cb);
    void resetReadOnly(Parameters param, const CkCallback &cb);

  public:
  static Tree::GenericTreeNode *pickNodeFromMergeList(int n, GenericTreeNode **gtn, int &nUnresolved, int &pickedIndex);
};

inline static void setBIconfig()
{
#if CHARM_VERSION > 60401 && CMK_BALANCED_INJECTION_API
    if (CkMyRank()==0) {
#define GNI_BI_DEFAULT    64
      uint16_t cur_bi = ck_get_GNI_BIConfig();
      if (cur_bi > GNI_BI_DEFAULT) {
        ck_set_GNI_BIConfig(GNI_BI_DEFAULT);
      }
    }
    if (CkMyPe() == 0)
      CkPrintf("Balanced injection is set to %d.\n", ck_get_GNI_BIConfig());
#endif
}

/** @brief Control recording of Charm++ projections logs
 *
 *  The constructors for this class are also used to set default
 *  node-wide communication parameters.
 */
class ProjectionsControl : public CBase_ProjectionsControl { 
  public: 
  ProjectionsControl() {
    setBIconfig();
    LBTurnCommOff();
    LBSetPeriod(0.0); // no need for LB interval: we are using Sync Mode
  } 
  ProjectionsControl(CkMigrateMessage *m) : CBase_ProjectionsControl(m) {
    setBIconfig();
    LBTurnCommOff();
    LBSetPeriod(0.0); // no need for LB interval: we are using Sync Mode
  } 
 
  void on(CkCallback cb) { 
    if(CkMyPe() == 0){ 
      CkPrintf("\n\n**** PROJECTIONS ON *****\n\n"); 
    } 
    traceBegin();  
    contribute(cb); 
  } 
 
  void off(CkCallback cb) { 
    if(CkMyPe() == 0){ 
      CkPrintf("\n\n**** PROJECTIONS OFF *****\n\n"); 
    } 
    traceEnd();  
    contribute(cb); 
  } 

  void pup(PUP::er &p){
    CBase_ProjectionsControl::pup(p);
  }
}; 


#endif //DATAMANAGER_H
