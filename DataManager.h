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

#ifdef COOLING_MOLECULARH
#define COOL_NV 5
#elif COOLING_METAL
#define COOL_NV 4
#elif COOLING_COSMO || COOLING_BOLEY
#define COOL_NV 1
#endif

/// @brief Information about TreePieces on an SMP node.
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
    /// Pointer to callback so it can be freed.
    CkCallback *cb;
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
        /// @brief counter for the number of tree nodes that are
        /// replicated by TreePieces that share the same address space.
        int cumNumReplicatedNodes;
        int treePiecesDone;
        int treePiecesDoneUdot;
        int savedChunk;
        int treePiecesDonePrefetch;
        int treePiecesDoneLocalComputation;
        // XXX - assumes that only one chunk can be on the gpu
        // at a given time
        int treePiecesDoneRemoteChunkComputation;
        int treePiecesWantParticlesBack;
        /// Reference count for Pieces that have finished updating
        /// their acclerations.
        int treePiecesParticlesUpdated;
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

        // TreePiece counter for multi-threaded GPU host buffer copy
	int treePiecesBufferFilled;

        // can the gpu accept a chunk of remote particles/nodes?
        bool gpuFree;

        /// Callback pointer to pass to HAPI.
        CkCallback *localTransferCallback;

        PendingBuffers *currentChunkBuffers;
        // queue that stores all pending chunk transfers
        CkQ<PendingBuffers *> pendingChunkTransferQ;

        // last remote chunk's size in moments and particles
        int lastChunkMoments;
        int lastChunkParticles;
        /// host buffer to transfer remote moments to GPU
        CudaMultipoleMoments *bufRemoteMoments;
        /// host buffer to transfer remote particles to GPU
        CompactPartData *bufRemoteParts;

        /// Vector to accumulate localMoments for transfering to GPU
        CkVec<CudaMultipoleMoments> localMoments;
        /// host buffer to transfer local moments to GPU
        CudaMultipoleMoments *bufLocalMoments;
        /// host buffer to transfer local particles to GPU
        CompactPartData *bufLocalParts;
        /// host buffer to transfer initial accelerations to GPU
        VariablePartData *bufLocalVars;

	// Pointers to particle and tree data on GPU
	CudaMultipoleMoments *d_localMoments;
        CudaMultipoleMoments *d_remoteMoments;
        CompactPartData *d_localParts;
	CompactPartData *d_remoteParts;
        VariablePartData *d_localVars;
        size_t sMoments;
        size_t sCompactParts;
        size_t sVarParts;

#ifndef COOLING_NONE
  clDerivsData *h_CoolData;
  STIFF *h_Stiff;
	double *h_ymin;
	double *h_y0;
	double *h_y1;
	double *h_q;
	double *h_d;
	double *h_rtau;
	double *h_ys;
	double *h_qs;
	double *h_rtaus;
	double *h_scrarray;
#ifdef CUDA
  // Pointers to cooling data on GPU
  clDerivsData *d_CoolData;
  STIFF *d_Stiff;
  double *d_y;
  double *d_dtg;
	double *d_ymin;
	double *d_y0;
	double *d_y1;
	double *d_q;
	double *d_d;
	double *d_rtau;
	double *d_ys;
	double *d_qs;
	double *d_rtaus;
	double *d_scrarray;

  // Used to determine if memory needs to be reallocated
  // Total number of gas particles on this node
  int numTotalGasParts;
  // Total number of gas particles allocated
  int numGasParts;
#endif // CUDA
#endif // COOLING_NONE

	int numStreams;
	cudaStream_t *streams;

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
	int Cool_nv; /// Number of variables for ODE solver
#ifdef CUDA
	COOL *d_Cool;
#if !defined(COOLING_NONE) && !defined(COOLING_BOLEY)
	RATES_T *d_Rates_T;
	UVSPECTRUM *d_Uvspectrum;
#if defined(COOLING_MOLECULARH) || defined(COOLING_METAL)
	float *d_MetalCoolln;
	float *d_MetalHeatln;
#endif
#endif
#endif
	/// @brief log of star formation events.
	///
	/// Star formation events are stored on the data manager since there
	/// is no need to migrate them with the TreePiece.
	StarLog *starLog;
	/// @brief Lock for accessing starlog from TreePieces
	CmiNodeLock lockStarLog;

	DataManager(const CkArrayID& treePieceID);
	DataManager(CkMigrateMessage *);

        void startLocalWalk();
        void resumeRemoteChunk();
#ifdef CUDA
	void createStreams(int _numStreams, const CkCallback& cb);
        void donePrefetch(int chunk); // serialize remote chunk wrapper
        void serializeLocalTree();
        void assignCUDAStreams(const CkCallback& cb);

        void setupuDot(int activeRung, int bAll, const CkCallback& cb);
        void setupuDotDone(const CkCallback& cb);

#ifdef GPU_LOCAL_TREE_WALK
        void transformLocalTreeRecursive(GenericTreeNode *node, CkVec<CudaMultipoleMoments>& localMoments);
#endif //GPU_LOCAL_TREE_WALK
        
        // actual serialization methods
        PendingBuffers *serializeRemoteChunk(GenericTreeNode *);
	void serializeLocal(GenericTreeNode *);
	void transferLocalToGPU(int nParts, GenericTreeNode *node);
        void freeLocalTreeMemory();
        void freeRemoteChunkMemory(int chunk);
        void transferParticleVarsBack();
        void updateParticles(UpdateParticlesStruct *data);
        void updateParticlesFreeMemory(UpdateParticlesStruct *data);
        void initiateNextChunkTransfer();
        DataManager(){ }

#endif

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
	    CmiDestroyLock(lockStarLog);
#ifdef CUDA
#if !defined(COOLING_NONE) && !defined(COOLING_BOLEY)
	    cudaFree(d_Rates_T);
	    cudaFree(d_MetalCoolln);
	    cudaFree(d_MetalHeatln);
	    cudaFree(d_Uvspectrum);
#endif

	    cudaFree(d_Cool);

	    cudaFree(d_CoolData);
	    cudaFree(d_Stiff);
	    cudaFree(d_dtg);

	    cudaFree(d_y);
	    cudaFree(d_ymin);
	    cudaFree(d_y0);
	    cudaFree(d_q);
	    cudaFree(d_d);
	    cudaFree(d_rtau);
	    cudaFree(d_ys);
	    cudaFree(d_qs);
	    cudaFree(d_rtaus);
	    cudaFree(d_scrarray);
	    cudaFree(d_y1);

            for (int i = 0; i < numStreams; i++) {
                cudaStreamDestroy(streams[i]);
	    }
	    delete[] streams;
#endif
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
    void clearRegisteredPieces(const CkCallback& cb);
    void combineLocalTrees(CkReductionMsg *msg);
    void getChunks(int &num, Tree::NodeKey *&roots);
    inline Tree::GenericTreeNode *chunkRootToNode(const Tree::NodeKey k) {
      NodeLookupType::iterator iter = chunkRootTable.find(k);
      if (iter != chunkRootTable.end()) return iter->second;
      else return NULL;
    }
    inline Tree::GenericTreeNode *getRoot() { return root; }
#ifdef CUDA
    void allocCoolParticleBlock(int numParts, int bFree);
#endif
    void initCooling(double dGmPerCcUnit, double dComovingGmPerCcUnit,
		     double dErgPerGmUnit, double dSecUnit, double dKpcUnit,
		     COOLPARAM inParam, const CkCallback& cb);
    void initStarLog(std::string _fileName, const CkCallback &cb);
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
#ifdef CUDA
    // GPUs are assigned to nodes in a round-robin fashion. This allows the user to define
    // one virtual node per device and utilize multiple GPUs on a single node
    // Beacuse devices are assigned per-PE, this is a convenient place to call setDevice
    // Note that this code has nothing to do with initalizing projections
    int numGpus;
    cudaGetDeviceCount(&numGpus);
    cudaSetDevice(CmiMyNode() % numGpus);
#endif
    setBIconfig();
    LBTurnCommOff();
#ifndef LB_MANAGER_VERSION
    // Older Charm++ requires this to avoid excessive delays between successive LBs even
    // when using AtSync mode
    LBSetPeriod(0.0);
#endif
  } 
  ProjectionsControl(CkMigrateMessage *m) : CBase_ProjectionsControl(m) {
    setBIconfig();
    LBTurnCommOff();
#ifndef LB_MANAGER_VERSION
    // Older Charm++ requires this to avoid excessive delays between successive LBs even
    // when using AtSync mode
    LBSetPeriod(0.0);
#endif
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
