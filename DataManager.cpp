/** \file DataManager.cpp
 \brief Implementation of the DataManager
 \author Graeme Lufkin (gwl@u.washington.edu)
*/

#include "ParallelGravity.h"
#include "DataManager.h"
#include "Reductions.h"
#include "formatted_string.h"

#ifdef CUDA
#include "hapi.h"
#include "cuda_typedef.h"
#include "SFC.h"
#endif

#include "Compute.h"
#include "TreeWalk.h"

void printTreeGraphViz(GenericTreeNode *node, ostream &out, const string &name);

DataManager::DataManager(const CkArrayID& treePieceID) {
  init();
  treePieces = CProxy_TreePiece(treePieceID);
}

DataManager::DataManager(CkMigrateMessage *m) : CBase_DataManager(m) {
  init();
}

void DataManager::init() {
  root = NULL;
  oldNumChunks = 0;
  chunkRoots = NULL;
#ifdef CUDA
  treePiecesDone = 0;
  treePiecesDonePrefetch = 0;
  treePiecesDoneLocalComputation = 0;
  treePiecesDoneRemoteChunkComputation = 0;
  treePiecesWantParticlesBack = 0;
  treePiecesParticlesUpdated = 0;
  gpuFree = true;

#endif
  Cool = CoolInit();
  LWData = LymanWernerTableInit();
  starLog = new StarLog();
  hmStarLog = new HMStarLog();
  lockStarLog = CmiCreateLock();
  lockHMStarLog = CmiCreateLock();
}

#ifdef CUDA
/// @brief Initialize CUDA streams
/// @param _numStreams Total number of streams to create
void DataManager::createStreams(int _numStreams, const CkCallback& cb) {
  numStreams = _numStreams;
  streams = new cudaStream_t[numStreams];
  for (int i = 0; i < numStreams; i++) {
      hapiCheck(cudaStreamCreate(&streams[i]));
  }
  contribute(cb);
}
#endif

/**
 * Fill in responsibleIndex after ORB decomposition
 */
void DataManager::acceptResponsibleIndex(const int* responsible, const int n,
					 const CkCallback& cb) {
    responsibleIndex.resize(n);
    copy(responsible, responsible + n, responsibleIndex.begin());
    contribute(cb);
    }

void DataManager::acceptFinalKeys(const SFC::Key* keys, const int* responsible, uint64_t* bins, const int n, const CkCallback& cb) {

  //should not assign responsibility or place to a treepiece that will get no particles
  int ignored = 0;
  for (int i = 0; i < n - 1; i ++){
    if (bins[i] == 0)
      ignored++;
  }
  boundaryKeys.resize(n - ignored);
  responsibleIndex.resize(n - 1 - ignored);
  particleCounts.resize(n - 1 - ignored);

  //if all treepieces receiving particles, copy everything
  if (ignored == 0){
    copy(keys, keys + n, boundaryKeys.begin());
    copy(responsible, responsible + n - 1, responsibleIndex.begin());
    copy(bins, bins + n - 1, particleCounts.begin());
  } else {
    boundaryKeys[0] = keys[0];
    int idx = 0;
    int last = 0;
      
    for (int i = 0; i < n; i ++){
      //skip empty treepieces by copying in chunks
      if (bins[i] == 0 || i == n - 1){
	copy(keys + last + 1, keys + i + 1, boundaryKeys.begin() + idx + 1);
	copy(responsible + last, responsible + i, responsibleIndex.begin() + idx);
	copy(bins + last, bins + i, particleCounts.begin() + idx);
	idx += i - last;
	last = i+1;

      } 
    }
  }


  if(verbosity >= 3 && CkMyPe()==0){
    std::vector<int>::iterator iter1;
    std::vector<int>::iterator iter2;
    CkPrintf("responsible,particleCounts:");
    for(iter1=responsibleIndex.begin(),iter2=particleCounts.begin();iter1!=responsibleIndex.end();iter1++,iter2++){
      CkPrintf("(%d,%d),",*iter1,*iter2);
    }
    CkPrintf("\n");
    std::vector<SFC::Key>::iterator iter3;
    CkPrintf("Keys:");
    for(iter3=boundaryKeys.begin();iter3!=boundaryKeys.end();iter3++){
      CkPrintf("%s,", make_formatted_string(*iter3).c_str());
    }
    CkPrintf("\n");
  }

  CkCallback unshuffleCallback(CkIndex_TreePiece::unshuffleParticles(0),treePieces);

  contribute(sizeof(CkCallback), &cb, callbackReduction, unshuffleCallback);
}

///
/// @brief Class to represent key pairs that are to be sorted on the
/// first key of the pair.
///
class KeyDouble {
  SFC::Key first;
  SFC::Key second;
public:
  inline bool operator<(const KeyDouble& k) const {
    return first < k.first;
  }
};

void DataManager::pup(PUP::er& p) {
    CBase_DataManager::pup(p);
    p | treePieces;
}

void DataManager::notifyPresence(Tree::GenericTreeNode *root, TreePiece *tp) {
  CmiLock(__nodelock);
  registeredTreePieces.push_back(TreePieceDescriptor(tp, root));
#ifdef CUDA
  //gpuFree = true;
  //registeredTreePieceIndices.push_back(index);
#if COSMO_PRINT_BK > 1
  CkPrintf("(%d) notifyPresence called by %d, length: %d\n", CkMyPe(), tp->getIndex(), registeredTreePieces.length());
#endif
#endif
  CmiUnlock(__nodelock);
}

/// \brief Clear registeredTreePieces on this node.
void DataManager::clearRegisteredPieces(const CkCallback& cb) {
    registeredTreePieces.removeAll();
    contribute(cb);
}


/// \brief Build a local tree inside the node.
///
/// This will be an exact superset of all the trees in this
/// processor. Only the minimum number of nodes is duplicated.  The
/// actual treebuilding is in a call to DataManager::buildProcessorTree.
/// @param msg contains the callback for when we are done.
///
void DataManager::combineLocalTrees(CkReductionMsg *msg) {
  int totalChares = registeredTreePieces.size();
  if (totalChares > 0) {
#ifdef PRINT_MERGED_TREE
    auto fout = make_formatted_string("cache.%d.dot",CkMyPe());
    ofstream ofs(fout.c_str());
#endif

#ifdef PRINT_MERGED_TREE
    for(int i = 0; i < registeredTreePieces.length(); i++){
      ostringstream oss;
      ostringstream name;

      oss << "tree." << registeredTreePieces[i].treePiece->getIndex() << "." << CkMyPe() << ".dot";
      name << "tree_" << registeredTreePieces[i].treePiece->getIndex();

      ofstream ofs1;
      ofs1.open(oss.str().c_str());
      printTreeGraphViz(registeredTreePieces[i].root,ofs1,name.str());
      ofs1.close();
    }
#endif

    // delete old tree
    for (int i = 0; i < nodeTable.length(); i++) {
      delete nodeTable[i];
    }
    nodeTable.clear();

#ifdef CUDA
    cumNumReplicatedNodes = 0;
#endif

    CkVec<Tree::GenericTreeNode*> gtn;
    for(int i = 0; i < registeredTreePieces.length(); i++){
      gtn.push_back(registeredTreePieces[i].root);
    }
    root = buildProcessorTree(totalChares, &gtn[0]);

#ifdef PRINT_MERGED_TREE
    ostringstream dmName;
    dmName << "dm_" << CkMyNode();

    printTreeGraphViz(root,ofs,dmName.str());
    ofs.close();
#endif

    // select the roots for the chunks in which computation will be divided
    if (_numChunks != oldNumChunks) {
      delete[] chunkRoots;
      oldNumChunks = _numChunks;
      chunkRoots = new Tree::NodeKey[_numChunks];
      root->getChunks(_numChunks, chunkRoots);
    } else {
      // TODO: update the chunks roots accordingly to some criteria

      // For now, simply riorder the keys if we have randomization on
      // (otherwise createLookupRoots fails if they are not sorted)
      root->getChunks(_numChunks, chunkRoots);
    }
    int numMappedRoots = createLookupRoots(root, chunkRoots);
    CkAssert(numMappedRoots = _numChunks);

    //Ramdomize the prefetchRoots
    if(_randChunks){
      for (int i=_numChunks; i>1; --i) {
        int r = rand();
        int k = (int) ((((float)r) * i) / (((float)RAND_MAX) + 1));
        Tree::NodeKey tmp = chunkRoots[i-1];
        chunkRoots[i-1] = chunkRoots[k];
        chunkRoots[k] = tmp;
      }
    }

  }
#ifdef CUDA
  gpuFree = true;
#endif
  contribute(*(CkCallback*)msg->getData());
  delete msg;
}

/// @brief Pick a node out of equivalent nodes on different
/// TreePieces.
/// If one of the nodes is internal to a TreePiece, return that one.
/// Otherwise pick from among the others.
/// @param n Number of equivalent nodes.
/// @param gtn Array of equivalent nodes.
/// @param nUnresolved Count of boundary nodes in the array (returned).
/// @param pickedIndex Index of picked Node.
Tree::GenericTreeNode *DataManager::pickNodeFromMergeList(int n, GenericTreeNode **gtn, int &nUnresolved, int &pickedIndex){
  int pick = -1;
  nUnresolved = 0;

  for (int i=0; i<n; ++i) {
    Tree::NodeType nt = gtn[i]->getType();
    if (nt == Tree::Internal || nt == Tree::Bucket) {
      // we can use this directly, noone else can have it other than NL
      CkAssert(nUnresolved == 0);
      pickedIndex = i;
      return gtn[i];
      // no change to the count of replicated nodes
    } else if (nt == Tree::Boundary) {
      // let's count up how many boundaries we find
      pick = i;
      nUnresolved++;
    } else {
      // here it can be NonLocal, NonLocalBucket or Empty. In all cases nothing to do.
    }
  }

  if(nUnresolved == 0){
    CkAssert(pick < 0);
    // only NonLocal (or Empty). any is good
    pickedIndex = 0;
    return gtn[0];
  }
  else{
    // multiple boundary nodes: return anyone of them
    pickedIndex = pick;
    return gtn[pick];
  } 
}

const char *typeString(NodeType type);
/**
 * \brief Build common tree for all pieces in a node.
 *
 * Given an array of pointers to an identical treenode in multiple
 * treepieces, return a node whose decendents will contain the union
 * of all those trees.  This is done recursively by calling this
 * function on each of the children of the treenode.  The recursion
 * stops if we hit an node that is totally contained in a single
 * processor, or there is only one copy of the node, or we have a node
 * that is non-local to all the treepieces.
 *
 * @param n number of nodes to process.
 * @param gtn array of nodes to process.  This contains the pointers
 * to the copies in each treepiece of an identical node.
 */
Tree::GenericTreeNode *DataManager::buildProcessorTree(int n, Tree::GenericTreeNode **gtn) {
#ifdef CUDA
  cumNumReplicatedNodes += (n-1);
#endif
  int nUnresolved;
  int pickedIndex;
  GenericTreeNode *pickedNode = pickNodeFromMergeList(n,gtn,nUnresolved,pickedIndex);

  /*
  ostringstream oss;
  for(int i = 0; i < n; i++){
    oss << "(" << gtn[i]->getKey() << "," << typeString(gtn[i]->getType()) << ")";
    if(gtn[i] == pickedNode) oss << " * ";
    oss << "; ";
  }
  CkPrintf("[%d] %s\n", CkMyPe(), oss.str().c_str());
  */

  if(nUnresolved <= 1){
    return pickedNode;
  }
  else{
   // more than one boundary, need recursion
    Tree::GenericTreeNode *newNode = pickedNode->clone();
#if INTERLIST_VER > 0
    //newNode->particlePointer = (GravityParticle *)0;
    //newNode->firstParticle = -1;
    //newNode->lastParticle = -1;
    newNode->startBucket = -1;
#endif
    // keep track if all the children are internal, in which case we have to
    // change this node type too from boundary to internal
    bool isInternal = true;
    nodeTable.push_back(newNode);
    CkVec<Tree::GenericTreeNode*> newgtn;
    // Recurse into common children.
    for (int child=0; child<gtn[0]->numChildren(); ++child) {
      for (int i=0; i<n; ++i) {
        if (gtn[i]->getType() == Tree::Boundary){
          GenericTreeNode *childNode = gtn[i]->getChildren(child);
          newgtn.push_back(childNode);
          //CkPrintf("[%d] (%llu,%s) add child (%llu,%s)\n", CkMyPe(), pickedNode->getKey(), typeString(pickedNode->getType()), childNode->getKey(), typeString(childNode->getType()));
        }
      }
      Tree::GenericTreeNode *ch = buildProcessorTree(newgtn.length(), &newgtn[0]);
      newgtn.length() = 0;
      newNode->setChildren(child, ch);
      if (ch->getType() == Tree::Boundary || ch->getType() == Tree::NonLocal || ch->getType() == Tree::NonLocalBucket) isInternal = false;
    }
    if (isInternal) {
      newNode->setType(Internal);
    }
    return newNode;
  }
}

/// \brief Fill in chunkRootTable, mapping keys to nodes.
int DataManager::createLookupRoots(Tree::GenericTreeNode *node, Tree::NodeKey *keys) {
  // assumes that the keys are ordered in tree depth first!
  if (node->getKey() == *keys) {
    // ok, found a chunk root, we can end the recursion
    //CkPrintf("mapping key %s\n",keyBits(*keys,KeyBits).c_str());
    chunkRootTable[*keys] = node;
    return 1;
  }
  // need to continue the recursion on the children
  int count = 0;
  for (int i=0; i<node->numChildren(); ++i) {
    Tree::GenericTreeNode *child = node->getChildren(i);
    int partial;
    if (child != NULL) {
      // child present, get the recursion going
      partial = createLookupRoots(node->getChildren(i), keys);
      keys += partial;
    } else {
      // the child does not exist, count the keys falling under it
      Tree::NodeKey childKey = node->getChildKey(i);
      for (partial=0; ; ++partial, ++keys) {
        int k;
        for (k=0; k<NodeKeyBits-1; ++k) {
          if (childKey == ((*keys)>>k)) break;
        }
        if (((*keys)|(~0 << k)) == ~0) break;
      }
      // add the last key found to the count
      ++partial;
      ++keys;
      //CkPrintf("missed keys of %s, %d\n",keyBits(childKey,63).c_str(),partial);
    }
    count += partial;
  }
  return count;
}

/// @brief return the number of chunks and the roots of the remote
/// walk subtrees.
/// @param num number of chunks (returned)
/// @param roots roots of the chunks (returned)
/// The remote walk is broken up into "chunks", which are subtrees.
/// This method returns the number of these chunks and the roots of
/// the corresponding subtrees.

void DataManager::getChunks(int &num, Tree::NodeKey *&roots) {
  num = oldNumChunks;
  roots = chunkRoots;
}

/*
 * obtain memory utilization statistics
 */
void DataManager::memoryStats(const CkCallback& cb)
{
    int mem = CmiMemoryUsage()/(1024*1024);
    contribute(sizeof(int), &mem, CkReduction::max_int, cb);
    }

/*
 * reset readonly variables after a restart
 */
void DataManager::resetReadOnly(Parameters param, const CkCallback &cb) 
{
    /*
     * Insert any variables that can change due to a restart.
     */
    _cacheLineDepth = param.cacheLineDepth;
    verbosity = param.iVerbosity;
    dExtraStore = param.dExtraStore;
    dMaxBalance = param.dMaxBalance;
    dFracLoadBalance = param.dFracLoadBalance;
    nIOProcessor = param.nIOProcessor;
    theta = param.dTheta;
    thetaMono = theta*theta*theta*theta;
#if CMK_SMP
    bUseCkLoopPar = param.bUseCkLoopPar;
#else
    bUseCkLoopPar = 0;
#endif
    contribute(cb);
    // parameter structure requires some cleanup
    delete param.stfm;
    free(param.csm);
    delete param.feedback;
}
  
	 
const char *typeString(NodeType type);

#ifdef CUDA
#include "HostCUDA.h"

/// @brief After all treepieces have registered, tranfer tree of entire node
/// to GPU.
void DataManager::serializeLocalTree(){
  CmiLock(__nodelock);
  treePiecesDone++;
#if COSMO_PRINT_BK > 1
  CkPrintf("(%d) serializeLocalTree treePiecesDone: %d, registered: %d\n", CkMyPe(), treePiecesDone, registeredTreePieces.length());
#endif
  if(treePiecesDone == registeredTreePieces.length()){
    treePiecesDone = 0;
    CmiUnlock(__nodelock);

    if(verbosity > 1)
        CkPrintf("[%d] Registered tree pieces length: %lu\n", CkMyPe(), registeredTreePieces.length());
    serializeLocal(root);
    if(verbosity > 1)
        CkPrintf("[%d] Registered tree pieces length after serialize local: %lu\n", CkMyPe(), registeredTreePieces.length());
  }
  else
      CmiUnlock(__nodelock);
}

/// @brief Callback from local data transfer to GPU
/// Indicate the transfer is done, and start the local gravity walks
/// on the treepieces on this node.
void DataManager::startLocalWalk() {
    delete localTransferCallback;

    for(int i = 0; i < registeredTreePieces.length(); i++){
      if(verbosity > 1) CkPrintf("[%d] GravityLocal %d\n", CkMyPe(), i);
      int in = registeredTreePieces[i].treePiece->getIndex();
      treePieces[in].commenceCalculateGravityLocal((intptr_t)d_localMoments, 
		                                   (intptr_t)d_localParts, 
						   (intptr_t)d_localVars,
						   (intptr_t)streams, numStreams,
		                                   sMoments, sCompactParts, sVarParts);
      if(registeredTreePieces[0].treePiece->bEwald) {
          EwaldMsg *msg = new (8*sizeof(int)) EwaldMsg;
          msg->fromInit = false;
          // Make priority lower than gravity or smooth.
          *((int *)CkPriorityPtr(msg)) = 3*numTreePieces + in + 1;
          CkSetQueueing(msg,CK_QUEUEING_IFIFO);
          treePieces[in].calculateEwald(msg);
      }
    }

    freePinnedHostMemory(bufLocalMoments);
    freePinnedHostMemory(bufLocalParts);
    freePinnedHostMemory(bufLocalVars);
}

/// @brief Callback from remote data transfer to GPU.
/// The data for remote interactions is on the GPU, so continue the
/// remote walk.
void DataManager::resumeRemoteChunk() {
  if(verbosity > 1) CkPrintf("[%d] resumeRemoteChunk registered: %lu\n", CkMyPe(), registeredTreePieces.length());
  int chunk = 0;
  chunk = currentChunkBuffers->chunk;
  delete currentChunkBuffers->moments;
  delete currentChunkBuffers->particles;
  delete currentChunkBuffers->cb;
  delete currentChunkBuffers;

  if(bufRemoteMoments != NULL)
      freePinnedHostMemory(bufRemoteMoments);
  if(bufRemoteParts != NULL)
      freePinnedHostMemory(bufRemoteParts);

    // resume each treepiece's startRemoteChunk, now that the nodes
    // are properly labeled and the particles accounted for
    for(int i = 0; i < registeredTreePieces.length(); i++){
      if(verbosity > 1) CkPrintf("[%d] resumeRemoteChunk %d\n", CkMyPe(), i);
      int in = registeredTreePieces[i].treePiece->getIndex();
#if COSMO_PRINT_BK > 1
      CkPrintf("(%d) dm->%d\n", CkMyPe(), in);
#endif
      treePieces[in].continueStartRemoteChunk(chunk, (intptr_t)d_remoteMoments, (intptr_t)d_remoteParts);
    }
}

/// @brief record when all TreePieces have finished their prefetch.
void DataManager::donePrefetch(int chunk){
  CmiLock(__nodelock);

  savedChunk = chunk;

  treePiecesDonePrefetch++;
  if(treePiecesDonePrefetch == registeredTreePieces.length()){
    treePiecesDonePrefetch = 0;
#ifdef HAPI_TRACE
    double starttime = CmiWallTimer();
#endif
    currentChunkBuffers = serializeRemoteChunk(root);
#ifdef HAPI_TRACE
    traceUserBracketEvent(CUDA_SER_TREE, starttime, CmiWallTimer());
#endif
    PendingBuffers *buffers = currentChunkBuffers;
    if(gpuFree){
      gpuFree = false;
      lastChunkMoments = buffers->moments->length();
      lastChunkParticles = buffers->particles->length();

      CkCallback *remoteChunkTransferCallback
          = new CkCallback(CkIndex_DataManager::resumeRemoteChunk(), CkMyNode(),
                           dMProxy);
      buffers->cb = remoteChunkTransferCallback;
      // XXX copies can be saved here.
      size_t sRemMoments = lastChunkMoments*sizeof(CudaMultipoleMoments);
      if(sRemMoments > 0) {
          allocatePinnedHostMemory((void **)&bufRemoteMoments, sRemMoments);
          memcpy(bufRemoteMoments, buffers->moments->getVec(), sRemMoments);
          }
      else
          bufRemoteMoments = NULL;
      size_t sRemParts = lastChunkParticles*sizeof(CompactPartData);
      if(sRemParts > 0) {
          allocatePinnedHostMemory((void **)&bufRemoteParts, sRemParts);
          memcpy(bufRemoteParts, buffers->particles->getVec(), sRemParts);
          }
      else
          bufRemoteParts = NULL;

      // Transfer moments and particle cores to gpu
      DataManagerTransferRemoteChunk(bufRemoteMoments, sRemMoments,
                                     bufRemoteParts, sRemParts,
                                     (void **)&d_remoteMoments,  (void **)&d_remoteParts,
				     streams[0],
                                     remoteChunkTransferCallback);

    }
    else{
      // enqueue pendingbuffers
      pendingChunkTransferQ.enq(buffers);
    }
    
  }
  CmiUnlock(__nodelock);
}

typedef std::map<KeyType, CkCacheEntry<KeyType>*> cacheType;

/// @brief add a node to the moment list and record its index
static inline void addNodeToList(GenericTreeNode *nd,
                                 CkVec<CudaMultipoleMoments> &list,
                                 int &index) {
    nd->nodeArrayIndex = index;
    list.push_back(CudaMultipoleMoments(nd->moments));
    index++;
}

/// @brief add a node the moment list (using pointer) and record its index
static inline void addNodeToListPtr(GenericTreeNode *nd,
                                    CkVec<CudaMultipoleMoments> *list,
                                    int &index) {
    nd->nodeArrayIndex = index;
    list->push_back(CudaMultipoleMoments(nd->moments));
    index++;
}

/// @brief add a node with walk information to the moment list and
/// record its index
#ifdef GPU_LOCAL_TREE_WALK
static inline void addTreeNodeToList(GenericTreeNode *nd,
                                     CkVec<CudaMultipoleMoments> &list,
                                     int &index) {
    nd->nodeArrayIndex = index;
    CudaMultipoleMoments cmm(nd->moments);
    cmm.lesser_corner = nd->boundingBox.lesser_corner;
    cmm.greater_corner = nd->boundingBox.greater_corner;
    list.push_back(cmm);
    index++;
    }
#endif //GPU_LOCAL_TREE_WALK

const char *typeString(NodeType type);

/// @brief Create a vector of remote nodes from a remote prefetch.
PendingBuffers *DataManager::serializeRemoteChunk(GenericTreeNode *node){
  CkQ<GenericTreeNode *> queue;
  int chunk = savedChunk;

  int numTreePieces = registeredTreePieces.length();
  int numNodes = 0;
  int numParticles = 0;
  int totalNumBuckets = 0;

  cacheType *wholeNodeCache = cacheNode.ckLocalBranch()->getCache();
  cacheType *ctNode = &wholeNodeCache[chunk];
  cacheType *wholePartCache = cacheGravPart.ckLocalBranch()->getCache();
  cacheType *ctPart = &wholePartCache[chunk];

  // find out number of particles and nodes cached
  // get them from cache - iterate and count each type

  CkVec<CudaMultipoleMoments> *postPrefetchMoments = new CkVec<CudaMultipoleMoments>;
  CkVec<CompactPartData> *postPrefetchParticles = new CkVec<CompactPartData>;
  PendingBuffers *pendingBuffers = new PendingBuffers;

  // XXX - better way to estimate NL, NLB, C, CB nodes/particles? 
  // thse are just guessed initial sizes for CkVecs
  numNodes = ctNode->size();
  numParticles = ctPart->size();

  postPrefetchMoments->reserve(numNodes);
  postPrefetchParticles->reserve(numParticles);

  postPrefetchMoments->length() = 0;
  postPrefetchParticles->length() = 0;

  // fill up postPrefetchMoments with node moments
  int nodeIndex = 0;
  int partIndex = 0;

#ifdef CUDA_DM_PRINT_TREES
  CkPrintf("*************\n");
  CkPrintf("[%d] DM remote chunk %d\n", CkMyPe(), chunk);
  CkPrintf("*************\n");
#endif
  queue.enq(node);
  while(!queue.isEmpty()){
    GenericTreeNode *node = queue.deq();
    NodeType type = node->getType();

    if(type == Empty || type == CachedEmpty || type == Internal || type == Bucket){ // skip
      continue;
    }// B, NL, NLBu, CBu, C 
    else if(type == Boundary){
      // enqueue children
      for(int i = 0; i < node->numChildren(); i++){
	GenericTreeNode *child = node->getChildren(i);
	queue.enq(child);
      }
    }
    else if(type == NonLocal){
      // need node moments; also, must enqueue children so that complete list of 
      // used nodes can be obtained
        addNodeToListPtr(node,postPrefetchMoments,nodeIndex);
    }
    else if(type == NonLocalBucket || type == CachedBucket){
      if(type == CachedBucket){
          addNodeToListPtr(node,postPrefetchMoments,nodeIndex);
      }
      // if this is a NonLocalBucket, don't need node itself, just its particles
      ExternalGravityParticle *parts;
      int nParticles = node->lastParticle-node->firstParticle+1;
      NodeKey key = node->getKey();
      // N.B. Key for particles is shifted to distinguish it from the Key
      // for the node.
      key <<= 1;

      cacheType::iterator p = ctPart->find(key);
      if (p != ctPart->end() && p->second->replyRecvd) {
        // found particles
        // mark presence and add to data to ship
        parts = (ExternalGravityParticle *)p->second->data;
        cachedPartsOnGpu[key] = partIndex;
#ifdef CUDA_DM_PRINT_TREES
        CkPrintf("(%d) type %s parts (key %ld) start: %d\n", CkMyPe(), 
                                                            typeString(type), key, partIndex);
#endif
        // put particles in array:
        for(int i = 0; i < nParticles; i++){
          postPrefetchParticles->push_back(CompactPartData(parts[i]));
          partIndex++;
        }
      }
    }
    else if(type == Cached){
      addNodeToListPtr(node,postPrefetchMoments,nodeIndex);
      // put children into queue, if available
      for(int i = 0 ; i < node->numChildren(); i++){
	GenericTreeNode *child = node->getChildren(i);
        if(child){// available to dm
    	  queue.enq(child);
        }
        else{ // look in cache
    	  NodeKey childKey = node->getChildKey(i);
          cacheType::iterator p = ctNode->find(childKey);
          if (p != ctNode->end() && p->second->replyRecvd) {
            // found node, enqueue
    	    queue.enq((GenericTreeNode *)p->second->data);
    	  }
        }
      }
    }
  }// end while queue not empty

#ifdef CUDA_DM_PRINT_TREES
  CkPrintf("*************\n");
#endif

  pendingBuffers->moments = postPrefetchMoments;
  pendingBuffers->particles = postPrefetchParticles;
  pendingBuffers->chunk = chunk;

  return pendingBuffers;

}// end serializeNodes

/// @brief gather local nodes and particles and send to GPU
/// @param nodeRoot Root of tree to walk.
void DataManager::serializeLocal(GenericTreeNode *nodeRoot){
  /// queue for breadth first treewalk.
  CkQ<GenericTreeNode *> queue;

  int numTreePieces = registeredTreePieces.length();
  int numNodes = 0;
  int numParticles = 0;

  for(int i = 0; i < numTreePieces; i++){
    TreePiece *tp = registeredTreePieces[i].treePiece;
    numNodes += tp->getNumNodes();
    numParticles += tp->getDMNumParticles();
  }
  numNodes -= cumNumReplicatedNodes;

  localMoments.reserve(numNodes);

  localMoments.length() = 0;

  // fill up postPrefetchMoments with node moments
  int nodeIndex = 0;

#ifdef CUDA_DM_PRINT_TREES
  CkPrintf("*************\n");
  CkPrintf("[%d] DM local tree\n", CkMyPe());
  CkPrintf("*************\n");
#endif
  double  starttime = CmiWallTimer();
  // Walk local tree
  queue.enq(nodeRoot);
  while(!queue.isEmpty()){
    GenericTreeNode *node = queue.deq();
    NodeType type = node->getType();

#ifdef CUDA_DM_PRINT_TREES
    //CkPrintf("Process [%d] %ld (%s)\n", CkMyPe(), node->getKey(), typeString(type));
#endif

    if(type == Empty || type == CachedEmpty){ // skip
      continue;
    }
    else if(type == Bucket || type == NonLocalBucket){ // NLB
      // don't need the particles, only the moments
#ifdef GPU_LOCAL_TREE_WALK
        addTreeNodeToList(node,localMoments,nodeIndex);
#else
        addNodeToList(node,localMoments,nodeIndex);
#endif //GPU_LOCAL_TREE_WALK
    }
    else if(type == Boundary || type == Internal){ // B,I 
#ifdef GPU_LOCAL_TREE_WALK
        addTreeNodeToList(node,localMoments,nodeIndex);
#else
        addNodeToList(node,localMoments,nodeIndex);
#endif //GPU_LOCAL_TREE_WALK
      for(int i = 0; i < node->numChildren(); i++){
        GenericTreeNode *child = node->getChildren(i);
        queue.enq(child);
      }
    }
  }// end while queue not empty

#ifdef HAPI_TRACE
  traceUserBracketEvent(SER_LOCAL_WALK, starttime, CmiWallTimer());
#endif

  // used later, when copying particle vars back to the host
  savedNumTotalParticles = numParticles;
  savedNumTotalNodes = localMoments.length();

#ifdef CUDA_DM_PRINT_TREES
  CkPrintf("*************\n");
#endif

#if COSMO_PRINT_BK > 1
  CkPrintf("(%d): DM->GPU local tree\n", CkMyPe());
#endif
  size_t sLocalParts = numParticles*sizeof(CompactPartData);
  size_t sLocalMoments = localMoments.length()*sizeof(CudaMultipoleMoments);
  allocatePinnedHostMemory((void **)&bufLocalParts, sLocalParts);
  allocatePinnedHostMemory((void **)&bufLocalMoments, sLocalMoments);

  int pTPindex = 0;
  treePiecesBufferFilled = 0;
  for(int i = 0; i < numTreePieces; i++){
      treePieces[registeredTreePieces[i].treePiece->getIndex()].fillGPUBuffer((intptr_t) bufLocalParts,
		      (intptr_t) bufLocalMoments, (intptr_t) localMoments.getVec(), pTPindex,
		      numParticles, (intptr_t) nodeRoot);
      pTPindex += registeredTreePieces[i].treePiece->getDMNumParticles();
      }
}

///
/// @brief After all pieces have filled the buffer, initiate the transfer.
/// @param numParticles total number of particles on this node
/// @param node root of tree
///
void DataManager::transferLocalToGPU(int numParticles, GenericTreeNode *node)
{
    CmiLock(__nodelock);
    treePiecesBufferFilled++;
    if(treePiecesBufferFilled == registeredTreePieces.length()){
        treePiecesBufferFilled = 0;
        CmiUnlock(__nodelock);
    }
    else {
        CmiUnlock(__nodelock);
        return;
    }

  double starttime = CmiWallTimer();
#ifdef GPU_LOCAL_TREE_WALK
  transformLocalTreeRecursive(node, localMoments);
#endif //GPU_LOCAL_TREE_WALK
#ifdef HAPI_TRACE
  traceUserBracketEvent(SER_LOCAL_TRANSFORM, starttime, CmiWallTimer());
#endif

  localTransferCallback
      = new CkCallback(CkIndex_DataManager::startLocalWalk(), CkMyNode(), dMProxy);

  // XXX copies can be saved here.
  starttime = CmiWallTimer();
  size_t sLocalVars = numParticles*sizeof(VariablePartData);
  size_t sLocalParts = numParticles*sizeof(CompactPartData);
  size_t sLocalMoments = localMoments.length()*sizeof(CudaMultipoleMoments);

  memcpy(bufLocalMoments, localMoments.getVec(), sLocalMoments);
#ifdef HAPI_TRACE
  traceUserBracketEvent(SER_LOCAL_MEMCPY, starttime, CmiWallTimer());
#endif

  allocatePinnedHostMemory((void **)&bufLocalVars, sLocalVars);

  // Transfer moments and particle cores to gpu
  DataManagerTransferLocalTree(bufLocalMoments, sLocalMoments, bufLocalParts,
                               sLocalParts, bufLocalVars, sLocalVars,
			       (void **)&d_localMoments, (void **)&d_localParts, (void **)&d_localVars,
			       streams[0], numParticles,
                               localTransferCallback);
}

#ifdef GPU_LOCAL_TREE_WALK
// Add more information to each Moment, basically transform moment to a computable tree node
void DataManager::transformLocalTreeRecursive(GenericTreeNode *node, CkVec<CudaMultipoleMoments>& localMoments) {
  NodeType type = node->getType();
  int node_index = node->nodeArrayIndex;

  if(type == Empty || type == CachedEmpty){ // skip
    return;
  } else if(type == Bucket || type == NonLocalBucket) {
    localMoments[node_index].type = (int)type;
    localMoments[node_index].nodeArrayIndex = node_index;
    localMoments[node_index].particleCount = node->particleCount;
    if(type == NonLocalBucket) // NonLocalBucket has no particles on
                               // the GPU.
        localMoments[node_index].bucketSize = 0;
    for (int i = 0; i < 2; i ++) {
      localMoments[node_index].children[i] = -1;
    }
  } else if(type == Boundary || type == Internal){ // B,I
    localMoments[node_index].type = (int)type;
    localMoments[node_index].nodeArrayIndex = node_index;
    localMoments[node_index].particleCount = node->particleCount;

    localMoments[node_index].bucketStart = INT_MAX;
    localMoments[node_index].bucketSize = 0;

    for (int i = 0; i < 2; i ++) {
      localMoments[node_index].children[i] = -1;
    }
    for(int i = 0; i < node->numChildren(); i++){
      GenericTreeNode *child = node->getChildren(i);
      int child_index = child->nodeArrayIndex;
      localMoments[node_index].children[i] = child_index;
      transformLocalTreeRecursive(child, localMoments);

      // child_index == -1 can indicate an empty node or a non-local node.
      if (child_index != -1 && localMoments[child_index].bucketSize > 0) {
        localMoments[node_index].bucketStart = std::min(localMoments[node_index].bucketStart, localMoments[child_index].bucketStart);
        localMoments[node_index].bucketSize += localMoments[child_index].bucketSize;
      }
    }
  }
}
#endif //GPU_LOCAL_TREE_WALK

void updateParticlesCallback(void *, void *);

/// @brief Copy particle accelerations back from GPU to host memory and
///        deallocate the device memory
/// This is triggered when all TreePieces call finishBucket
void DataManager::transferParticleVarsBack(){
  UpdateParticlesStruct *data;
  CmiLock(__nodelock);
  treePiecesWantParticlesBack++;
  if(treePiecesWantParticlesBack == registeredTreePieces.length()){
    treePiecesWantParticlesBack = 0; 
    VariablePartData *buf;
    
    if(savedNumTotalParticles > 0){
#ifdef PINNED_HOST_MEMORY
      allocatePinnedHostMemory((void **)&buf, savedNumTotalParticles*sizeof(VariablePartData));
#else
      buf = (VariablePartData *) malloc(savedNumTotalParticles*sizeof(VariablePartData));
#endif
    }
    else{
      buf = NULL;
    }

    data = new UpdateParticlesStruct;
    data->cb = new CkCallback(updateParticlesCallback, data);
    data->dm = this;
    data->buf = buf;
    data->size = savedNumTotalParticles;

    if(verbosity > 1) CkPrintf("[%d] transferParticleVarsBack\n", CkMyPe());
    TransferParticleVarsBack(buf, 
                             savedNumTotalParticles*sizeof(VariablePartData),
			     d_localVars,
			     streams[0],
                             data->cb);
    
    cudaFree(d_localMoments);
    cudaFree(d_localParts);
    cudaFree(d_localVars);
    cudaFree(d_remoteMoments);
    cudaFree(d_remoteParts); 

#ifdef CUDA_PRINT_ERRORS
    printf("transferParticleVarsBack: %s\n", cudaGetErrorString( cudaGetLastError() ) );
#endif
  }
  CmiUnlock(__nodelock);
}

void DataManager::updateParticles(UpdateParticlesStruct *data){
  int partIndex = 0;

  VariablePartData *deviceParticles = data->buf;

#ifdef CUDA_PRINT_TRANSFER_BACK_PARTICLES
  CkPrintf("(%d) In DM::updateParticles %d tps\n", CkMyPe(), registeredTreePieces.length());
#endif

  for(int i = 0; i < registeredTreePieces.length(); i++){
    TreePiece *tp = registeredTreePieces[i].treePiece;
    int numParticles = tp->getNumParticles();
#ifdef CUDA_PRINT_TRANSFER_BACK_PARTICLES
    CkPrintf("(%d) tp %d, numParticles: %d\n", CkMyPe(), tp->getIndex(), numParticles);
#endif
#ifdef CHANGA_REFACTOR_MEMCHECK 
    CkPrintf("(%d) memcheck before updating tp %d particles\n", CkMyPe(), tp->getIndex());
    CmiMemoryCheck();
#endif

    // N.B. passing a pointer to an entry method is not normally done.
    // It is OK here because we know that the treePiece is on our own
    // SMP node, but we need a kludgey cast to get it to work.
    treePieces[registeredTreePieces[i].treePiece->getIndex()].updateParticles((intptr_t) data, partIndex);
    partIndex += numParticles;

  }

}

void updateParticlesCallback(void *param, void *msg){  
  UpdateParticlesStruct *data = (UpdateParticlesStruct *)param;
  data->dm->updateParticles(data);
}

/// @brief clean up buffer for GPU transfer back.
void DataManager::updateParticlesFreeMemory(UpdateParticlesStruct *data)
{
    CmiLock(__nodelock);
    treePiecesParticlesUpdated++;
    if(treePiecesParticlesUpdated == registeredTreePieces.length()){
        treePiecesParticlesUpdated = 0;

        if(data->size > 0){
#ifdef PINNED_HOST_MEMORY
            freePinnedHostMemory(data->buf);
#else
            free(data->buf);
#endif
        }
        delete (data->cb);
        delete data;
    }
    CmiUnlock(__nodelock);
}

#endif // CUDA
