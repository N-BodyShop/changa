#ifndef _HOST_CUDA_H_
#define _HOST_CUDA_H_


#include <cuda_runtime.h>
#include "cuda_typedef.h"

#define THREADS_PER_BLOCK 128

#ifdef GPU_LOCAL_TREE_WALK
#define THREADS_PER_WARP 32
#define WARPS_PER_BLOCK (THREADS_PER_BLOCK / THREADS_PER_WARP)
#define WARP_INDEX (threadIdx.x >> 5)
#endif //GPU_LOCAL_TREE_WALK

#ifdef CUDA_2D_TB_KERNEL
#define PARTS_PER_BLOCK 16
#define NODES_PER_BLOCK (THREADS_PER_BLOCK/PARTS_PER_BLOCK)

#define THREADS_PER_BLOCK_PART 128
#define PARTS_PER_BLOCK_PART 16
#define NODES_PER_BLOCK_PART (THREADS_PER_BLOCK_PART/PARTS_PER_BLOCK_PART)
#endif

// FIXME - find appropriate values
#define NUM_INIT_MOMENT_INTERACTIONS_PER_BUCKET 100
#define NUM_INIT_PARTICLE_INTERACTIONS_PER_BUCKET 100


// Usetable fns only
#define CLRATES( _cl, _Rate, _T, _rho, _ZMetal, _columnL, _Rate_Phot_H2_stellar) cudaClRates_Table( _cl, _Rate, _T, _rho, _ZMetal, _columnL, _Rate_Phot_H2_stellar)
#define CLEDOTINSTANT( _cl, _Y, _Rate, _rho, _ZMetal, _Heat, _Cool ) cudaClEdotInstant_Table( _cl, _Y, _Rate, _rho, _ZMetal, _Heat, _Cool)


#define EPS 1e-5
#define M_H      1.672e-24

#ifdef CUBICTABLEINTERP
#define TABLEFACTOR 2
#else 
#define TABLEFACTOR 1
#endif

#ifdef CUBICTABLEINTERP
#define TABLEINTERP( _rname ) (wTln0*RT0->_rname+wTln1*RT1->_rname+wTln0d*RT0d->_rname+wTln1d*RT1d->_rname)
#else
#define TABLEINTERP( _rname ) (wTln0*RT0->_rname+wTln1*RT1->_rname)
#endif

#define TABLEINTERPLIN( _rname ) (wTln0*RT0->_rname+wTln1*RT1->_rname)

#define CL_Rgascode         8.2494e7
#define CL_Eerg_gm_degK     CL_Rgascode
#define CL_ev_degK          1.0/1.1604e4
#define CL_Eerg_gm_ev       CL_Eerg_gm_degK/CL_ev_degK
#define CL_Eerg_gm_degK3_2  1.5*CL_Eerg_gm_degK /* 1.23741e8*/

/** @brief Data and parameters for requesting gravity calculations on
 * the GPU. */
typedef struct _CudaRequest{
        /// CUDA stream to handle memory operations and kernel launches for this request
        /// Allocation and deallocation of the stream is handled by the DataManager
	cudaStream_t stream;

	/// for accessing device memory
	CudaMultipoleMoments *d_localMoments;
	CudaMultipoleMoments *d_remoteMoments;
	CompactPartData *d_localParts;
	CompactPartData *d_remoteParts;
	VariablePartData *d_localVars;
	size_t sMoments;
	size_t sCompactParts;
	size_t sVarParts;

        /// can either be a ILCell* or an ILPart*
	void *list;
	int *bucketMarkers;     /**< index in the cell or particle
                                 * list for the cells/particles for a
                                 * given bucket */
	int *bucketStarts;      /**< index in the particle array on
                                 * the GPU of the bucket */
	int *bucketSizes;       /**< number of particles in a bucket
                                 * on the GPU  */
	int numInteractions;    /**< Total number of interactions in
                                 * this request  */
	int numBucketsPlusOne;  /**< Number of buckets affected by
                                 *  this request */
        void *tp;               /**< Pointer to TreePiece that made
                                 * this request */
        /// pointer to off-processor Node/Particle buffer.
        void *missedNodes;
        void *missedParts;
        /// Size of the off-processor data buffer.
        size_t sMissed;

	/// these buckets were finished in this work request
	int *affectedBuckets;
        void *cb;               /**< Callback  */
        void *state;            /**< Pointer to state of walk that
                                 * created this request  */
        cudatype fperiod;       /**< Size of periodic volume  */

        // TODO: remove these later if we don't use COSMO_PRINT_BK.
        /// is this a node or particle computation request?
        bool node;
        /// is this a remote or local computation?
        bool remote;
#ifdef GPU_LOCAL_TREE_WALK
  int firstParticle;
  int lastParticle;
  int rootIdx;
  cosmoType theta;
  cosmoType thetaMono;
  int nReplicas;
  cudatype fperiodY;  // Support periodic boundary condition in more dimensions
  cudatype fperiodZ;  // Support periodic boundary condition in more dimensions
#endif //GPU_LOCAL_TREE_WALK
}CudaRequest;

/// Device memory pointers used by most functions in HostCUDA
typedef struct _CudaDevPtr{
    void *d_list;
    int *d_bucketMarkers;
    int *d_bucketStarts;
    int *d_bucketSizes;
}CudaDevPtr;

void allocatePinnedHostMemory(void **, size_t);
void freePinnedHostMemory(void *);

void DataManagerTransferLocalTree(void *moments, size_t sMoments,
                                  void *compactParts, size_t sCompactParts,
                                  void *varParts, size_t sVarParts,
				  void **d_localMoments, void **d_compactParts, void **d_varParts,
				  cudaStream_t stream, int numParticles,
                                  void *callback);
void DataManagerTransferRemoteChunk(void *moments, size_t sMoments,
                                  void *compactParts, size_t sCompactParts,
				  void **d_remoteMoments, void **d_remoteParts,
				  cudaStream_t stream,
                                  void *callback);

void CudaCoolSetTime( CudaCOOL *cl, double dTime, double z, cudaStream_t stream );
void TransferParticleVarsBack(VariablePartData *hostBuffer, size_t size, void *d_varParts, cudaStream_t stream, void *cb);

void TreePieceCellListDataTransferLocal(CudaRequest *data);
void TreePieceCellListDataTransferRemote(CudaRequest *data);
void TreePieceCellListDataTransferRemoteResume(CudaRequest *data);


void TreePiecePartListDataTransferLocal(CudaRequest *data);
void TreePiecePartListDataTransferLocalSmallPhase(CudaRequest *data, CompactPartData *parts, int len);
void TreePiecePartListDataTransferRemote(CudaRequest *data);
void TreePiecePartListDataTransferRemoteResume(CudaRequest *data);

void TreePieceODESolver(CudaSTIFF *d_CudaStiff, double *d_y, double *d_dtg, double  **y_in, double tstart, std::vector<double> dtg, int numParts, cudaStream_t stream);

void DummyKernel(void *cb);

#endif
