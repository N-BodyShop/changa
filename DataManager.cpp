/** \file DataManager.cpp
 \brief Implementation of the DataManager
 \author Graeme Lufkin (gwl@u.washington.edu)
*/

#include "ParallelGravity.h"
#include "DataManager.h"
#include "Reductions.h"

#ifdef CUDA
#include "cuda_typedef.h"
#include "SFC.h"

#ifdef CUDA_INSTRUMENT_WRS
//#define GPU_INSTRUMENT_WRS
//#include "wr.h"
void hapi_initInstrument(int, char);
void hapi_clearInstrument();
#endif

#endif

DataManager::DataManager(const CkArrayID& treePieceID) {
  init();
  treePieces = CProxy_TreePiece(treePieceID);
}

DataManager::DataManager(CkMigrateMessage *m) : CBase_DataManager(m) {
  init();
}

void DataManager::init() {
  splitters = NULL;
  numSplitters = 0;
  root = NULL;
  oldNumChunks = 0;
  chunkRoots = NULL;
#ifdef CUDA
  treePiecesDone = 0;
  treePiecesDonePrefetch = 0;
  treePiecesDoneLocalComputation = 0;
  treePiecesDoneRemoteChunkComputation = 0;
  treePiecesWantParticlesBack = 0;
#ifdef CUDA_INSTRUMENT_WRS
  treePiecesDoneInitInstrumentation = 0;
#endif

  gpuFree = true;
#endif
  Cool = CoolInit();
}

/**
 * Fill in responsibleIndex after ORB decomposition
 */
void DataManager::acceptResponsibleIndex(const int* responsible, const int n,
					 const CkCallback& cb) {
    responsibleIndex.resize(n);
    copy(responsible, responsible + n, responsibleIndex.begin());
    contribute(cb);
    }

void DataManager::acceptFinalKeys(const SFC::Key* keys, const int* responsible, unsigned int* bins, const int n, const CkCallback& cb) {

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
      CkPrintf("%016llx,",*iter3);
    }
    CkPrintf("\n");
  }

	contribute(sizeof(CkCallback), &cb, callbackReduction,
		   CkCallback(CkIndex_TreePiece::unshuffleParticles(0),
			      treePieces));
}

class KeyDouble {
  SFC::Key first;
  SFC::Key second;
public:
  inline bool operator<(const KeyDouble& k) const {
    return first < k.first;
  }
};

/**
 * Collect and sort the boundaries of all treepieces so that they are
 * available to each treepiece as they do their tree build.  The
 * treebuild continues with a call to TreePiece::startOctTreeBuild.
 */
void DataManager::collectSplitters(CkReductionMsg *m) {
  numSplitters = m->getSize() / sizeof(SFC::Key);
  CkAssert(! (numSplitters&1)); // must be even
  CkAssert(numSplitters > 0);
  delete[] splitters;
  splitters = new SFC::Key[numSplitters];
  SFC::Key* splits = static_cast<SFC::Key *>(m->getData());
  std::copy(splits, splits + numSplitters, splitters);
  // The splitters come in pairs (1st and last key of each treepiece).
  // Sort them as pairs.
  KeyDouble* splitters2 = (KeyDouble *)splitters;
  std::sort(splitters2, splitters2 + (numSplitters>>1));
  for (unsigned int i=1; i<numSplitters; ++i) {
    if (splitters[i] < splitters[i-1]) {
      if(CkMyNode()==0)
        CkAbort("Keys not ordered");
    }
  }
  splitters[0] = SFC::firstPossibleKey;
  contribute(0, 0, CkReduction::concat, CkCallback(CkIndex_TreePiece::startOctTreeBuild(0), treePieces));
  delete m;
  if(verbosity > 3)
    ckerr << "DataManager " << CkMyNode() << ": Collected splitters" << endl;

}

void DataManager::pup(PUP::er& p) {
    CBase_DataManager::pup(p);
    p | treePieces;
}

#ifdef CUDA
void DataManager::notifyPresence(Tree::GenericTreeNode *root, TreePiece *tp, int index) {
#else
void DataManager::notifyPresence(Tree::GenericTreeNode *root) {
#endif

  CmiLock(__nodelock);
  registeredChares.push_back(root);

#ifdef CUDA
  registeredTreePieces.push_back(TreePieceDescriptor(tp, index));
  //gpuFree = true;
  //registeredTreePieceIndices.push_back(index);
#if COSMO_PRINT_BK > 1
  CkPrintf("(%d) notifyPresence called by %d, length: %d\n", CkMyPe(), index, registeredTreePieces.length());
#endif
#endif
  CmiUnlock(__nodelock);
}

/// \brief Build a local tree inside the node.
///
/// This will be an exact superset of all the trees in this
/// processor. Only the minimum number of nodes is duplicated.  The
/// actual treebuilding is in a call to DataManager::buildProcessorTree.
/// @param msg contains the callback for when we are done.
///
void DataManager::combineLocalTrees(CkReductionMsg *msg) {
  int totalChares = registeredChares.size();
  if (totalChares > 0) {
#if COSMO_DEBUG > 0
    char fout[100];
    sprintf(fout,"cache.%d.%d",CkMyPe(),iterationNo);
    ofs = new ofstream(fout);
#endif
    Tree::GenericTreeNode **gtn = registeredChares.getVec();
    // delete old tree
    Tree::NodeLookupType::iterator nodeIter;
    for (nodeIter = nodeLookupTable.begin(); nodeIter != nodeLookupTable.end(); nodeIter++) {
      delete nodeIter->second;
    }
    nodeLookupTable.clear();
#ifdef CUDA
    cumNumReplicatedNodes = 0;
#endif
    root = buildProcessorTree(totalChares, gtn);
    registeredChares.removeAll();
#if COSMO_DEBUG > 0
    printGenericTree(root,*ofs);
    ofs->close();
    delete ofs;
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
      srand((CkMyNode()+1)*1000);
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
  contribute(0, 0, CkReduction::random, *(CkCallback*)msg->getData());
  delete msg;
}

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
  //CkAssert(n > 1); // the recursion should have stopped before!
  int pick = -1;
  int count = 0;
#ifdef CUDA
  cumNumReplicatedNodes += (n-1);
#endif
  for (int i=0; i<n; ++i) {
    Tree::NodeType nt = gtn[i]->getType();
    if (nt == Tree::Internal || nt == Tree::Bucket) {
      // we can use this directly, noone else can have it other than NL
#if COSMO_DEBUG > 0
      (*ofs) << "cache "<<CkMyPe()<<": "<<keyBits(gtn[i]->getKey(),63)<<" using Internal node"<<endl;
#endif
      return gtn[i];
      // no change to the count of replicated nodes
    } else if (nt == Tree::Boundary) {
      // let's count up how many boundaries we find
      pick = i;
      count++;
    } else {
      // here it can be NonLocal, NonLocalBucket or Empty. In all cases nothing to do.
    }
  }
  if (count == 0) {
    // only NonLocal (or Empty). any is good
#if COSMO_DEBUG > 0
    (*ofs) << "cache "<<CkMyPe()<<": "<<keyBits(gtn[0]->getKey(),63)<<" using NonLocal node"<<endl;
#endif
    return gtn[0];
  } else if (count == 1) {
    // only one single Boundary, all others are NonLocal, use this directly
#if COSMO_DEBUG > 0
    (*ofs) << "cache "<<CkMyPe()<<": "<<keyBits(gtn[pick]->getKey(),63)<<" using Boundary node"<<endl;
#endif
    return gtn[pick];
  } else {
    // more than one boundary, need recursion
    Tree::GenericTreeNode *newNode = gtn[pick]->clone();
#if INTERLIST_VER > 0
    //newNode->particlePointer = (GravityParticle *)0;
    //newNode->firstParticle = -1;
    //newNode->lastParticle = -1;
    newNode->startBucket = -1;
#endif
    // keep track if all the children are internal, in which case we have to
    // change this node type too from boundary to internal
    bool isInternal = true;
#if COSMO_DEBUG > 0
    (*ofs) << "cache "<<CkMyPe()<<": "<<keyBits(newNode->getKey(),63)<<" duplicating node"<<endl;
#endif
    nodeLookupTable[newNode->getKey()] = newNode;
    Tree::GenericTreeNode **newgtn = new Tree::GenericTreeNode*[count];
    // Recurse into common children.
    for (int child=0; child<gtn[0]->numChildren(); ++child) {
      for (int i=0, j=0; i<n; ++i) {
        if (gtn[i]->getType() == Tree::Boundary) newgtn[j++]=gtn[i]->getChildren(child);
      }
      Tree::GenericTreeNode *ch = buildProcessorTree(count, newgtn);
      newNode->setChildren(child, ch);
      if (ch->getType() == Tree::Boundary || ch->getType() == Tree::NonLocal || ch->getType() == Tree::NonLocalBucket) isInternal = false;
    }
    delete[] newgtn;
    if (isInternal) {
      newNode->setType(Internal);
#if COSMO_DEBUG > 0
      (*ofs) << "cache "<<CkMyPe()<<": "<<keyBits(newNode->getKey(),63)<<" converting to Internal"<<endl;
#endif
    }
    return newNode;
  }
}

/// \brief Fill in chunkRootTable, mapping keys to nodes.
int DataManager::createLookupRoots(Tree::GenericTreeNode *node, Tree::NodeKey *keys) {
  // assumes that the keys are ordered in tree depth first!
  if (node->getKey() == *keys) {
    // ok, found a chunk root, we can end the recursion
    //CkPrintf("mapping key %s\n",keyBits(*keys,63).c_str());
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
        for (k=0; k<63; ++k) {
          if (childKey == ((*keys)>>k)) break;
        }
        if (((*keys)|(~0 << k)) == ~0) break;
        //if (k==63) break;
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
    dExtraStore = param.dExtraStore;
    contribute(cb);
    }
  
	 
const char *typeString(NodeType type);

#ifdef CUDA
#include "HostCUDA.h"

#ifdef CUDA_INSTRUMENT_WRS
int DataManager::initInstrumentation(){
  CmiLock(__nodelock);
  int saved = treePiecesDoneInitInstrumentation;
  treePiecesDoneInitInstrumentation++;
  if(treePiecesDoneInitInstrumentation == registeredTreePieces.length()){
    activeRung = registeredTreePieces[0].tp->getActiveRung();
    hapi_initInstrument(registeredTreePieces.length(), BOTTOM_EWALD_KERNEL+1);
    treePiecesDoneInitInstrumentation = 0;
  }
  CmiUnlock(__nodelock);
  return saved;
}
#endif

void DataManager::serializeLocalTree(){
  CmiLock(__nodelock);
  treePiecesDone++;
#if COSMO_PRINT_BK > 1
  CkPrintf("(%d) serializeLocalTree treePiecesDone: %d, registered: %d\n", CkMyPe(), treePiecesDone, registeredTreePieces.length());
#endif
  if(treePiecesDone == registeredTreePieces.length()){
    treePiecesDone = 0;

#ifdef CUDA_TRACE
    double starttime = CmiWallTimer();
#endif
    serializeLocal(root);
#ifdef CUDA_TRACE
    traceUserBracketEvent(CUDA_SER_TREE, starttime, CmiWallTimer());
#endif
    // resume each treepiece's startRemoteChunk, now that the nodes
    // are properly labeled and the particles accounted for
    for(int i = 0; i < registeredTreePieces.length(); i++){
      int in = registeredTreePieces[i].index;
#if COSMO_PRINT_BK > 1
      CkPrintf("(%d) dm->%d\n", CkMyPe(), in);
#endif
      treePieces[in].commenceCalculateGravityLocal();
    }
  }
  CmiUnlock(__nodelock);
}

void DataManager::donePrefetch(int chunk){
  CmiLock(__nodelock);

  //if(savedChunk < 0){
  //  savedChunk = chunk;
  //}
  //CkAssert(savedChunk == chunk);
  savedChunk = chunk;

  treePiecesDonePrefetch++;
  if(treePiecesDonePrefetch == registeredTreePieces.length()){
    treePiecesDonePrefetch = 0;
#ifdef CUDA_TRACE
    double starttime = CmiWallTimer();
#endif
    PendingBuffers *buffers = serializeRemoteChunk(root);
#ifdef CUDA_TRACE
    traceUserBracketEvent(CUDA_SER_TREE, starttime, CmiWallTimer());
#endif
    if(gpuFree){
      gpuFree = false;
      lastChunkMoments = buffers->moments->length();
      lastChunkParticles = buffers->particles->length();
      //CkPrintf("(%d) DM donePrefetch gpuFree, transferring 0x%x (%d); 0x%x (%d) \n", CkMyPe(), buffers->moments->getVec(), lastChunkMoments, buffers->particles->getVec(), lastChunkParticles);

      // Transfer moments and particle cores to gpu
#ifdef CUDA_INSTRUMENT_WRS
      DataManagerTransferRemoteChunk(buffers->moments->getVec(), lastChunkMoments, buffers->particles->getVec(), lastChunkParticles, 0, activeRung);
#else
      DataManagerTransferRemoteChunk(buffers->moments->getVec(), lastChunkMoments, buffers->particles->getVec(), lastChunkParticles);
#endif

      delete buffers->moments;
      delete buffers->particles;
      delete buffers;
      // resume each treepiece's startRemoteChunk, now that the nodes
      // are properly labeled and the particles accounted for
      for(int i = 0; i < registeredTreePieces.length(); i++){
        int in = registeredTreePieces[i].index;
        //CkPrintf("(%d) dm->%d chunk %d\n", CkMyPe(), in, chunk);
        treePieces[in].continueStartRemoteChunk(chunk);
      }
    }
    else{
      // enqueue pendingbuffers
      //CkPrintf("(%d) DM donePrefetch gpu not free, enqueuing\n", CkMyPe());
      pendingChunkTransferQ.enq(buffers);
    }
    
  }
  CmiUnlock(__nodelock);
}

typedef std::map<CkCacheKey, CkCacheEntry*> cacheType;

#ifdef CUDA_DM_PRINT_TREES 
#define addNodeToList(nd, list, index) \
      { \
        nd->nodeArrayIndex = index; \
        nd->wasNeg = false; \
        list.push_back(CudaMultipoleMoments(nd->moments));\
        CkPrintf("(%d) node %d: %ld (%s)\n", CkMyPe(), index, nd->getKey(), typeString(type));\
        index++;\
      }
#define addNodeToListPtr(nd, list, index) \
      { \
        nd->nodeArrayIndex = index; \
        nd->wasNeg = false; \
        list->push_back(CudaMultipoleMoments(nd->moments));\
        CkPrintf("(%d) node %d: %ld (%s)\n", CkMyPe(), index, nd->getKey(), typeString(type));\
        index++;\
      }

#else
#define addNodeToList(nd, list, index) \
      { \
        nd->nodeArrayIndex = index; \
        nd->wasNeg = false; \
        list.push_back(CudaMultipoleMoments(nd->moments));\
        index++;\
      }
#define addNodeToListPtr(nd, list, index) \
      { \
        nd->nodeArrayIndex = index; \
        nd->wasNeg = false; \
        list->push_back(CudaMultipoleMoments(nd->moments));\
        index++;\
      }

#endif


const char *typeString(NodeType type);

PendingBuffers *DataManager::serializeRemoteChunk(GenericTreeNode *node){
  CkQ<GenericTreeNode *> queue;
  int chunk = savedChunk;

  int numTreePieces = registeredTreePieces.length();
  int numNodes = 0;
  int numParticles = 0;
  int numCachedNodes = 0;
  int numCachedParticles = 0;
  int totalNumBuckets = 0;

  cacheType *wholeNodeCache = cacheNode[CkMyPe()].getCache();
  cacheType *ctNode = &wholeNodeCache[chunk];
  cacheType *wholePartCache = cacheGravPart[CkMyPe()].getCache();
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

  //postPrefetchMoments = new CudaMultipoleMoments[numNodes];
  //postPrefetchParticles = new CompactPartData[numParticles];

  // needed so we know how many particles there are in each bucket
  //int *bmarks = new int[totalNumBuckets+1];

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
      addNodeToListPtr(node,postPrefetchMoments,nodeIndex)
    }
    else if(type == NonLocalBucket || type == CachedBucket){
      if(type == CachedBucket){
        addNodeToListPtr(node,postPrefetchMoments,nodeIndex)
      }
      // if this is a NonLocalBucket, don't need node itself, just its particles
      ExternalGravityParticle *parts;
      int nParticles = node->lastParticle-node->firstParticle+1;
      NodeKey key = node->getKey();
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
      addNodeToListPtr(node,postPrefetchMoments,nodeIndex)
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


void DataManager::serializeLocal(GenericTreeNode *node){
  CkQ<GenericTreeNode *> queue;

  int numTreePieces = registeredTreePieces.length();
  int numNodes = 0;
  int numParticles = 0;
  int numCachedNodes = 0;
  int numCachedParticles = 0;

  for(int i = 0; i < numTreePieces; i++){
    TreePiece *tp = registeredTreePieces[i].tp;
    numNodes += tp->getNumNodes();
    numParticles += tp->getDMNumParticles();
  }
  numNodes -= cumNumReplicatedNodes;

  CkVec<CudaMultipoleMoments> localMoments;
  CkVec<CompactPartData> localParticles;

  localMoments.reserve(numNodes);
  localParticles.resize(numParticles);

  localMoments.length() = 0;

  // fill up postPrefetchMoments with node moments
  int nodeIndex = 0;
  int partIndex = 0;

#ifdef CUDA_DM_PRINT_TREES
  CkPrintf("*************\n");
  CkPrintf("[%d] DM local tree\n", CkMyPe());
  CkPrintf("*************\n");
#endif
  queue.enq(node);
  while(!queue.isEmpty()){
    GenericTreeNode *node = queue.deq();
    NodeType type = node->getType();

#ifdef CUDA_DM_PRINT_TREES
    //CkPrintf("Process [%d] %ld (%s)\n", CkMyPe(), node->getKey(), typeString(type));
#endif

    if(type == Empty || type == CachedEmpty){ // skip
      continue;
    }
    else if(type == Bucket){ // Bu, follow pointer
      // copy particles
      // need both the node moments and the particles
      NodeKey bucketKey = node->getKey();

      // nodes
      addNodeToList(node,localMoments,nodeIndex)
    }
    else if(type == NonLocalBucket){ // NLB
      // don't need the particles, only the moments
      addNodeToList(node,localMoments,nodeIndex)
    }
    else if(type == Boundary || type == Internal){ // B,I 
      addNodeToList(node,localMoments,nodeIndex)
      for(int i = 0; i < node->numChildren(); i++){
        GenericTreeNode *child = node->getChildren(i);
        queue.enq(child);
      }
    }
  }// end while queue not empty

  // used later, when copying particle vars back to the host
  savedNumTotalParticles = numParticles;
  savedNumTotalNodes = localMoments.length();

  for(int i = 0; i < registeredTreePieces.length(); i++){
    TreePiece *tp = registeredTreePieces[i].tp;
    tp->getDMParticles(localParticles.getVec(), partIndex);
  }
#ifdef CUDA_DM_PRINT_TREES
  CkPrintf("*************\n");
#endif
  CkAssert(partIndex == numParticles);
#if COSMO_PRINT_BK > 1
  CkPrintf("(%d): DM->GPU local tree\n", CkMyPe());
#endif
  // Transfer moments and particle cores to gpu
#ifdef CUDA_INSTRUMENT_WRS
  DataManagerTransferLocalTree(localMoments.getVec(), localMoments.length(), localParticles.getVec(), partIndex, 0, activeRung);
#else
  DataManagerTransferLocalTree(localMoments.getVec(), localMoments.length(), localParticles.getVec(), partIndex, CkMyPe());
#endif

}// end serializeLocal

#if 0
TreePieceDescriptor *DataManager::findKeyInDescriptors(SFC::Key particleKey){
  for(int i = 0; i < registeredTreePieces.length(); i++){
    if(registeredTreePieces[i].firstParticleKey == particleKey){
      return &registeredTreePieces[i];
    }
}

return 0;
}
#endif


void DataManager::freeLocalTreeMemory(){
  CmiLock(__nodelock);
  treePiecesDoneLocalComputation++;
  if(treePiecesDoneLocalComputation == registeredTreePieces.length()){
    treePiecesDoneLocalComputation = 0; 
#ifdef CUDA_INSTRUMENT_WRS
    FreeDataManagerLocalTreeMemory(savedNumTotalNodes > 0, savedNumTotalParticles > 0, 0, activeRung);
#else
    FreeDataManagerLocalTreeMemory(savedNumTotalNodes > 0, savedNumTotalParticles > 0);
#endif
  }
  CmiUnlock(__nodelock);
}

void DataManager::freeRemoteChunkMemory(int chunk){
  CmiLock(__nodelock);
  treePiecesDoneRemoteChunkComputation++;
  if(treePiecesDoneRemoteChunkComputation == registeredTreePieces.length()){
    treePiecesDoneRemoteChunkComputation = 0; 
    //CkPrintf("(%d) DM freeRemoteChunkMemory chunk %d called\n", CkMyPe(), chunk);
#ifdef CUDA_INSTRUMENT_WRS
    FreeDataManagerRemoteChunkMemory(chunk, (void *)this, lastChunkMoments != 0, lastChunkParticles != 0, 0, activeRung);
#else
    FreeDataManagerRemoteChunkMemory(chunk, (void *)this, lastChunkMoments != 0, lastChunkParticles != 0);
#endif
  }
  CmiUnlock(__nodelock);
}

void initiateNextChunkTransfer(void *dm_){
  DataManager *dm = (DataManager *)dm_;
  dm->initiateNextChunkTransfer();
}

void DataManager::initiateNextChunkTransfer(){
  PendingBuffers *next = 0;
  if(next = pendingChunkTransferQ.deq()){
    // Transfer moments and particle cores to gpu
    int chunk = next->chunk;
    //CkPrintf("(%d) DM initiateNextChunkTransfer chunk %d (%d moments, %d particles) called\n", CkMyPe(), chunk, next->moments->length(), next->particles->length());
    lastChunkMoments = next->moments->length();
    lastChunkParticles = next->particles->length();

    CkPrintf("(%d) DM initiateNextChunkTransfer chunk %d, 0x%x (%d); 0x%x (%d) \n", CkMyPe(), next->moments->getVec(), lastChunkMoments, next->particles->getVec(), lastChunkParticles);
#ifdef CUDA_INSTRUMENT_WRS
    DataManagerTransferRemoteChunk(next->moments->getVec(), next->moments->length(), next->particles->getVec(), next->particles->length(), 0, activeRung);
#else
    DataManagerTransferRemoteChunk(next->moments->getVec(), next->moments->length(), next->particles->getVec(), next->particles->length());
#endif

    delete next->moments;
    delete next->particles;
    delete next;
    // resume each treepiece's startRemoteChunk, now that the nodes
    // are properly labeled and the particles accounted for
    for(int i = 0; i < registeredTreePieces.length(); i++){
      int in = registeredTreePieces[i].index;
      //CkPrintf("(%d) dm->%d\n", CkMyPe(), in);
      treePieces[in].continueStartRemoteChunk(chunk);
    }
  }
  else{
    //CkPrintf("(%d) DM initiateNextChunkTransfer no chunks found, gpu is free\n", CkMyPe());
    gpuFree = true;
  }
}

void updateParticlesCallback(void *, void *);

void allocatePinnedHostMemory(void **ptr, int size);
void freePinnedHostMemory(void *ptr);

void DataManager::transferParticleVarsBack(){
  UpdateParticlesStruct *data;
  CmiLock(__nodelock);
  treePiecesWantParticlesBack++;
  if(treePiecesWantParticlesBack == registeredTreePieces.length()){
    treePiecesWantParticlesBack = 0; 
    VariablePartData *buf;
    
    if(savedNumTotalParticles > 0){
#ifdef CUDA_USE_CUDAMALLOCHOST
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

#ifdef CUDA_INSTRUMENT_WRS
    TransferParticleVarsBack(buf, savedNumTotalParticles*sizeof(VariablePartData), data->cb, savedNumTotalNodes > 0, savedNumTotalParticles > 0, 0, activeRung);
#else
    TransferParticleVarsBack(buf, savedNumTotalParticles*sizeof(VariablePartData), data->cb, savedNumTotalNodes > 0, savedNumTotalParticles > 0);
#endif
  }
  CmiUnlock(__nodelock);
}

void DataManager::updateParticles(UpdateParticlesStruct *data){
  int partIndex = 0;

  // FIXME - not required?
  CmiLock(__nodelock);

  VariablePartData *deviceParticles = data->buf;

#ifdef CUDA_PRINT_TRANSFER_BACK_PARTICLES
  CkPrintf("(%d) In DM::updateParticles %d tps\n", CkMyPe(), registeredTreePieces.length());
#endif

  for(int i = 0; i < registeredTreePieces.length(); i++){
    TreePiece *tp = registeredTreePieces[i].tp;
    int numParticles = tp->getNumParticles();
#ifdef CUDA_PRINT_TRANSFER_BACK_PARTICLES
    CkPrintf("(%d) tp %d, numParticles: %d\n", CkMyPe(), tp->getIndex(), numParticles);
#endif
#ifdef CHANGA_REFACTOR_MEMCHECK 
    CkPrintf("(%d) memcheck before updating tp %d particles\n", CkMyPe(), tp->getIndex());
    CmiMemoryCheck();
#endif

    if(tp->largePhase()){
      for(int j = 1; j <= numParticles; j++){
        if(tp->isActive(j)){
#ifndef CUDA_NO_ACC_UPDATES
          // FIXME - used to be +=
          tp->myParticles[j].treeAcceleration.x = deviceParticles[partIndex].a.x; //+ tp->myParticles[j].treeAcceleration;
          tp->myParticles[j].treeAcceleration.y = deviceParticles[partIndex].a.y; //+ tp->myParticles[j].treeAcceleration;
          tp->myParticles[j].treeAcceleration.z = deviceParticles[partIndex].a.z; //+ tp->myParticles[j].treeAcceleration;
          tp->myParticles[j].potential = deviceParticles[partIndex].potential;
#endif
#ifdef CUDA_PRINT_TRANSFER_BACK_PARTICLES
          CkPrintf("particle %d device: (%f,%f,%f) host: (%f,%f,%f)\n",
              j, 
              deviceParticles[partIndex].a.x,
              deviceParticles[partIndex].a.y,
              deviceParticles[partIndex].a.z,
              tp->myParticles[j].treeAcceleration.x,
              tp->myParticles[j].treeAcceleration.y,
              tp->myParticles[j].treeAcceleration.z);
#endif
        }
        partIndex++;
      }
    }
    else{
      for(int j = 1; j <= numParticles; j++){
        if(tp->isActive(j)){
#ifndef CUDA_NO_ACC_UPDATES
          // FIXME - used to be +=
          tp->myParticles[j].treeAcceleration.x = deviceParticles[partIndex].a.x; // + tp->myParticles[j].treeAcceleration;
          tp->myParticles[j].treeAcceleration.y = deviceParticles[partIndex].a.y; // + tp->myParticles[j].treeAcceleration;
          tp->myParticles[j].treeAcceleration.z = deviceParticles[partIndex].a.z; // + tp->myParticles[j].treeAcceleration;
          tp->myParticles[j].potential = deviceParticles[partIndex].potential;
#endif
#ifdef CUDA_PRINT_TRANSFER_BACK_PARTICLES
          CkPrintf("particle %d device: (%f,%f,%f) host: (%f,%f,%f)\n",
              j, 
              deviceParticles[partIndex].a.x,
              deviceParticles[partIndex].a.y,
              deviceParticles[partIndex].a.z,
              tp->myParticles[j].treeAcceleration.x,
              tp->myParticles[j].treeAcceleration.y,
              tp->myParticles[j].treeAcceleration.z);
#endif
          partIndex++;
        }     
      }
    }

#ifdef CHANGA_REFACTOR_MEMCHECK 
    CkPrintf("(%d) memcheck after updating tp %d particles\n", CkMyPe(), tp->getIndex());
    CmiMemoryCheck();
#endif

    // tell treepiece to go ahead with 
    // iteration wrap-up
    treePieces[registeredTreePieces[i].index].continueWrapUp();
  }

  registeredTreePieces.length() = 0;
  // FIXME - not required?
  CmiUnlock(__nodelock); 
}

void updateParticlesCallback(void *param, void *msg){  
  UpdateParticlesStruct *data = (UpdateParticlesStruct *)param;
  data->dm->updateParticles(data);
  if(data->size > 0){
#ifdef CUDA_USE_CUDAMALLOCHOST
    freePinnedHostMemory(data->buf);
#else
    free(data->buf);
#endif
  }
  delete (data->cb);
  delete data;
}

#endif

void DataManager::clearInstrument(CkCallback &cb){
#ifdef CUDA_INSTRUMENT_WRS
  hapi_clearInstrument();
  contribute(0,0,CkReduction::concat,cb);
#endif
}


