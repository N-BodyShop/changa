/** \file DataManager.cpp
 Implementation of the DataManager
 \author Graeme Lufkin (gwl@u.washington.edu)
*/

#include "ParallelGravity.h"
#include "DataManager.h"
#include "Reductions.h"

#ifdef CUDA
#include "cuda_typedef.h"
#include "SFC.h"
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
}

/*
void DataManager::acceptCandidateKeys(const SFC::Key* keys, const int n, int isRefine, const CkCallback& cb) {
  if (localCache == NULL) {
    localCache = cacheManagerProxy.ckLocalBranch();
  }
	boundaryKeys.resize(n);
	copy(keys, keys + n, boundaryKeys.begin());
	//tell the TreePieces on this node to evaluate the splitter keys
	map<int,GenericTreeNode*> *myTreePieces = localCache->getRegisteredChares();
	for(map<int,GenericTreeNode*>::iterator iter = myTreePieces->begin(); iter != myTreePieces->end(); ++iter)
		treePieces[iter->first].evaluateBoundaries(isRefine, cb);
}
*/

void DataManager::acceptFinalKeys(const SFC::Key* keys, const int* responsible, unsigned int* bins, const int n, const CkCallback& cb) {
	boundaryKeys.resize(n);
	copy(keys, keys + n, boundaryKeys.begin());
	responsibleIndex.resize(n - 1);
	copy(responsible, responsible + n - 1, responsibleIndex.begin());
	particleCounts.resize(n - 1);
	copy(bins, bins + n - 1, particleCounts.begin());

  if(verbosity >= 3 && CkMyPe()==0){
    std::vector<int>::iterator iter1;
    std::vector<int>::iterator iter2;
    CkPrintf("responsible,particleCounts:");
    for(iter1=responsibleIndex.begin(),iter2=particleCounts.begin();iter1!=responsibleIndex.end();iter1++,iter2++){
      CkPrintf("(%d,%d),",*iter1,*iter2);
    }
    CkPrintf("\n");
    /*CkPrintf("Particle Counts:");
    for(iter=particleCounts.begin();iter!=particleCounts.end();iter++){
      if(*iter==0)
        CkPrintf("%d,",*iter);
    }
    CkPrintf("\n");*/
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

	//tell my TreePieces to move the particle data to the responsible chare
	//for(vector<int>::iterator iter = myTreePieces.begin(); iter != myTreePieces.end(); ++iter)
	//	treePieces[*iter].unshuffleParticles(cb);
}

class KeyDouble {
  SFC::Key first;
  SFC::Key second;
public:
  inline bool operator<(const KeyDouble& k) const {
    return first < k.first;
  }
};

void DataManager::collectSplitters(CkReductionMsg *m) {
  numSplitters = m->getSize() / sizeof(SFC::Key);
//  numSplitters = 2 * numTreePieces;
  delete[] splitters;
  splitters = new SFC::Key[numSplitters];
  SFC::Key* splits = static_cast<SFC::Key *>(m->getData());
  std::copy(splits, splits + numSplitters, splitters);
  KeyDouble* splitters2 = (KeyDouble *)splitters;
  //sort(splitters, splitters + numSplitters);
  std::sort(splitters2, splitters2 + numTreePieces);
  for (unsigned int i=1; i<numSplitters; ++i) {
    if (splitters[i] < splitters[i-1]) {
      //for (unsigned int j=0; j<numSplitters; ++j)
      //  CkPrintf("%d: Key %d = %016llx\n",thisIndex,j,splitters[j]);
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
void DataManager::notifyPresence(Tree::GenericTreeNode *root, TreePiece *tp, int index, SFC::Key firstParticleKey, int numParticles) {
#else
void DataManager::notifyPresence(Tree::GenericTreeNode *root) {
#endif

  CmiLock(__nodelock);
  registeredChares.push_back(root);

#ifdef CUDA
  registeredTreePieces.push_back(TreePieceDescriptor(tp, index, firstParticleKey, numParticles));
  //registeredTreePieceIndices.push_back(index);
  CkPrintf("(%d) notifyPresence called by %d, length: %d\n", CkMyPe(), index, registeredTreePieces.length());
#endif
  CmiUnlock(__nodelock);
}

void DataManager::combineLocalTrees(CkReductionMsg *msg) {
  // build a local tree inside the node. This will be an exact superset of all
  // the trees in this processor. Only the minimum number of nodes is duplicated
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
  contribute(0, 0, CkReduction::random, *(CkCallback*)msg->getData());
  delete msg;
}

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

const char *typeString(NodeType type);

#ifdef CUDA
#include "HostCUDA.h"

void DataManager::serializeLocalTree(){
  CmiLock(__nodelock);
  treePiecesDone++;
  if(treePiecesDone == registeredTreePieces.length()){
    treePiecesDone = 0;
    serializeLocal(root);
    // resume each treepiece's startRemoteChunk, now that the nodes
    // are properly labeled and the particles accounted for
    for(int i = 0; i < registeredTreePieces.length(); i++){
      int in = registeredTreePieces[i].index;
      //CkPrintf("(%d) dm->%d\n", CkMyPe(), in);
      treePieces[in].commenceCalculateGravityLocal();
    }
  }
  CmiUnlock(__nodelock);
}

void DataManager::donePrefetch(int chunk){
  CmiLock(__nodelock);
  treePiecesDone++;
  if(treePiecesDonePrefetch == registeredTreePieces.length()){
    treePiecesDonePrefetch = 0;
    serializeRemoteChunk(root);
    // resume each treepiece's startRemoteChunk, now that the nodes
    // are properly labeled and the particles accounted for
    for(int i = 0; i < registeredTreePieces.length(); i++){
      int in = registeredTreePieces[i].index;
      //CkPrintf("(%d) dm->%d\n", CkMyPe(), in);
      treePieces[in].continueStartRemoteChunk();
    }
  }
  CmiUnlock(__nodelock);
}

#if 0
void DataManager::donePrefetch(){

	treePiecesDone++;
	if(treePiecesDone == registeredTreePieces.length()){
		treePiecesDone = 0;
		CkPrintf("(%d) registered: %d\n", CkMyPe(), registeredTreePieces.length());
		serializeNodes(root);

		// resume each treepiece's startRemoteChunk, now that the nodes
		// are properly labeled and the particles accounted for
		for(int i = 0; i < registeredTreePieces.length(); i++){
			int in = registeredTreePieces[i].index;
			CkPrintf("(%d) dm->%d\n", CkMyPe(), in);
			treePieces[in].continueStartRemoteChunk();
		}
		// don't delete list of registered treepieces yet; need the descriptors
		// to update particle acceleration, etc. on the treepieces when callbacks
		// received.
		//registeredTreePieces.length() = 0;
		//registeredTreePieceIndices.length() = 0;
	}
}

void DataManager::serializeNodes(GenericTreeNode *node){
  CkQ<GenericTreeNode *> queue;
  queue.enq(node);

  int numTreePieces = registeredTreePieces.length();
  int numNodes = 0;
  int numParticles = 0;
  int numCachedNodes = 0;
  int numCachedParticles = 0;
  int totalNumBuckets = 0;

  std::map<CkCacheKey, CkCacheEntry*> *cache = cacheManagerProxy[CkMyPe()].getCache();

  for(int i = 0; i < numTreePieces; i++){
	TreePiece *tp = registeredTreePieces[i].tp;
    numNodes += tp->getNumNodes();
    numParticles += tp->getNumParticles();
    totalNumBuckets += tp->getNumBuckets();
  }
  numNodes -= cumNumReplicatedNodes;

  // find out number of particles and nodes cached
  // get them from cache - iterate and count each type
  std::map<CkCacheKey, CkCacheEntry*>::iterator it;

  for(it = cache->begin();it != cache->end(); it++){
    CkCacheEntryType *type = it->second->type;

    /*
    if((EntryTypeGravityParticle * e = dynamic_cast<EntryTypeGravityParticle *>(type))){
      numCachedParticles += e->;
    }
    else
    */
    if(dynamic_cast<EntryTypeGravityNode *>(type)){
      numCachedNodes++;
    }
    // don't need to count smooth particles yet
  }

  CudaMultipoleMoments *postPrefetchMoments;
  CompactPartData *postPrefetchParticles;

  postPrefetchMoments = new CudaMultipoleMoments[numNodes+numCachedNodes];
  postPrefetchParticles = new CompactPartData[numParticles+numCachedParticles];

  // needed so we know how many particles there are in each bucket
  //int *bmarks = new int[totalNumBuckets+1];

  // fill up postPrefetchMoments with node moments
  int nodeIndex = 0;
  int partIndex = 0;
  int bucketIndex = 0;

  while(!queue.isEmpty()){
    GenericTreeNode *node = queue.deq();
    NodeType type = node->getType();

    if(type == Empty || type == CachedEmpty){ // skip
      continue;
    }
    else if(type == Bucket){ // follow pointer
      // copy particles
      int start = node->firstParticle;
      int end = node->lastParticle;
      GravityParticle *gravParts = node->particlePointer;
      NodeKey bucketKey = node->getKey();

      node->nodeArrayIndex = nodeIndex;
      postPrefetchMoments[nodeIndex] = node->moments;
#ifdef CUDA_PRINT_POST_PREFETCH_LIST
    CkPrintf("node %d: %ld (%s)\n", nodeIndex, node->getKey(), typeString(type));
#endif
      nodeIndex++;
      localPartsOnGpu[bucketKey << 1] = partIndex;
      //bmarks[bucketIndex] = partIndex;
#ifdef CUDA_PRINT_POST_PREFETCH_LIST
      CkPrintf("local bucket %d: %ld (parts: %d-%d)\n", bucketIndex, bucketKey, partIndex, partIndex+(end-start));
#endif
      for(int i = start; i <= end; i++){
        postPrefetchParticles[partIndex] = gravParts[i-start];

        SFC::Key particleKey = gravParts[i-start].key;
        // if the particle is the first of some registered treepiece, record this fact
        TreePieceDescriptor *descriptor = 0;
        descriptor = descriptor = findKeyInDescriptors(particleKey);
        if(descriptor){
			// the particleArrayStartIndex is used to update the particles of treepieces
        	// once all computation on the gpu has finished and the accelerations, velocities,
        	// potentials have been copied back.
        	descriptor->particleArrayStartIndex = partIndex;
        }
        partIndex++;
      }
      bucketIndex++;
    }
    else if(type == CachedBucket || type == NonLocalBucket){ // request cache for particles
      int numParticles = node->lastParticle-node->firstParticle+1;
      // ask cache for particles:
      NodeKey key = node->getKey();
      node->nodeArrayIndex = nodeIndex;
      postPrefetchMoments[nodeIndex] = node->moments;
#ifdef CUDA_PRINT_POST_PREFETCH_LIST
      CkPrintf("node %d: %ld (%s)\n", nodeIndex, node->getKey(), typeString(type));
#endif
      nodeIndex++;
      /*
      NodeKey partKey = key << 1;
      ExternalGravityParticle *parts = (ExternalGravityParticle *)  cacheManagerProxy[CkMyPe()].requestDataNoFetch(partKey, 0);
      if(parts){
        cachedPartsOnGpu[partKey] = partIndex;
#ifdef CUDA_PRINT_POST_PREFETCH_LIST
        CkPrintf("cached bucket %d: %ld (parts: %d-%d)\n", bucketIndex, key, partIndex, partIndex+numParticles-1);
#endif
        for(int i = 0; i < numParticles; i++){
          postPrefetchParticles[partIndex] = parts[i];
          partIndex++;
        }
        bucketIndex++;
      }
      */
    }
    else if(type == NonLocal || type == Cached){
    	NodeKey key = node->getKey();
    	node->nodeArrayIndex = nodeIndex;
    	postPrefetchMoments[nodeIndex] = node->moments;
#ifdef CUDA_PRINT_POST_PREFETCH_LIST
    	CkPrintf("node %d: %ld (%s)\n", nodeIndex, node->getKey(), typeString(type));
#endif
    	nodeIndex++;

    	// enqueue children if available
    	for(int i = 0; i < node->numChildren(); i++){
    		GenericTreeNode *child = node->getChildren(i);
    		if(child){// available to dm
    			queue.enq(child);
    		}
    		else{ // look in cache
    			NodeKey childKey = node->getChildKey(i);
    			CkCacheEntry *cacheEntry = cacheManagerProxy[CkMyPe()].requestCacheEntryNoFetch(childKey, 0);
    			if(cacheEntry && cacheEntry->replyRecvd){ // available in cache
    				GenericTreeNode *child = (GenericTreeNode *)cacheEntry->data;
    				queue.enq(child);
    			}
    		}
    	}
    }
    else{ // Boundary or Internal
    	node->nodeArrayIndex = nodeIndex;
    	postPrefetchMoments[nodeIndex] = node->moments;
#ifdef CUDA_PRINT_POST_PREFETCH_LIST
    	CkPrintf("node %d: %ld (%s)\n", nodeIndex, node->getKey(), typeString(type));
#endif
    	nodeIndex++;
    	for(int i = 0; i < node->numChildren(); i++){
    		GenericTreeNode *child = node->getChildren(i);
    		queue.enq(child);
    	}
    }
  }// end while queue not empty
  CkAssert(bucketIndex == totalNumBuckets);
  CkAssert(partIndex == numParticles+numCachedParticles);

  // Transfer moments and particle cores to gpu
  int bufferID = CkMyPe();
  DataManagerTransfer(postPrefetchMoments, numNodes+numCachedNodes, postPrefetchParticles, numParticles+numCachedParticles, bufferID);

  /*
  for(int i = 0; i < numParticles+numCachedParticles; i++){
	  if(postPrefetchParticles[i].)
  }
  */

}// end serializeNodes
#endif

void DataManager::serializeRemoteChunk(GenericTreeNode *node){
  CkQ<GenericTreeNode *> queue;

  int numTreePieces = registeredTreePieces.length();
  int numNodes = 0;
  int numParticles = 0;
  int numCachedNodes = 0;
  int numCachedParticles = 0;
  int totalNumBuckets = 0;

  std::map<CkCacheKey, CkCacheEntry*> *cache = cacheManagerProxy[CkMyPe()].getCache();

  // find out number of particles and nodes cached
  // get them from cache - iterate and count each type
  std::map<CkCacheKey, CkCacheEntry*>::iterator it;
  int cacheSize = cache->size();

  CkVec<CudaMultipoleMoments> postPrefetchMoments;
  CkVec<CompactPartData> postPrefetchParticles;

  // FIXME - better way to estimate NL, NLB, C, CB nodes/particles? 
  // thse are just estimates, initial sizes for CkVec's
  numNodes = size;
  numParticles = size;

  postPrefetchMoments.reserve(numNodes);
  postPrefetchParticles.reserve(numParticles);

  postPrefetchMoments.length() = 0;
  postPrefetchParticles.length() = 0;

  //postPrefetchMoments = new CudaMultipoleMoments[numNodes];
  //postPrefetchParticles = new CompactPartData[numParticles];

  // needed so we know how many particles there are in each bucket
  //int *bmarks = new int[totalNumBuckets+1];

  // fill up postPrefetchMoments with node moments
  int nodeIndex = 0;
  int partIndex = 0;
  int bucketIndex = 0;

  queue.enq(node);
  while(!queue.isEmpty()){
    GenericTreeNode *node = queue.deq();
    NodeType type = node->getType();

    // get all cached/cachedbucket nodes later, by iterating over the cache's contents.
    // look only for NL, NLB nodes right now
    if(type == Empty || type == CachedEmpty || type == Internal || type == Bucket || type == CachedBucket || type == Cached){ // skip
      continue;
    }
    else if(type == Boundary){
	for(int i = 0; i < node->numChildren(); i++){
    		GenericTreeNode *child = node->getChildren(i);
    		queue.enq(child);
    	}
    }
    else if(type == NonLocal){
      node->nodeArrayIndex = nodeIndex;
      postPrefetchMoments.push_back(node->moments);
#ifdef CUDA_PRINT_POST_PREFETCH_LIST
      CkPrintf("node %d: %ld (%s)\n", nodeIndex, node->getKey(), typeString(type));
#endif
      nodeIndex++;
    }
    else if(type == NonLocalBucket){
      node->nodeArrayIndex = nodeIndex;
      postPrefetchMoments.push_back(node->moments);
#ifdef CUDA_PRINT_POST_PREFETCH_LIST
      CkPrintf("node %d: %ld (%s)\n", nodeIndex, node->getKey(), typeString(type));
#endif
      nodeIndex++;

      // FIXME also need particles  
    }
  }// end while queue not empty

  // FIXME - start here
  for(it = cache->begin(); it != cache->end(); it++){
    
  }

  // Transfer moments and particle cores to gpu
  // FIXME - fix bufferid
  int bufferID = CkMyPe();
  DataManagerTransfer(postPrefetchMoments, numNodes, postPrefetchParticles, numParticles, bufferID);

  /*
  for(int i = 0; i < numParticles+numCachedParticles; i++){
	  if(postPrefetchParticles[i].)
  }
  */

}// end serializeNodes


void DataManager::serializeLocal(GenericTreeNode *node){
  CkQ<GenericTreeNode *> queue;

  int numTreePieces = registeredTreePieces.length();
  int numNodes = 0;
  int numParticles = 0;
  int numCachedNodes = 0;
  int numCachedParticles = 0;
  int totalNumBuckets = 0;

  std::map<CkCacheKey, CkCacheEntry*> *cache = cacheManagerProxy[CkMyPe()].getCache();

  for(int i = 0; i < numTreePieces; i++){
	TreePiece *tp = registeredTreePieces[i].tp;
    numNodes += tp->getNumNodes();
    numParticles += tp->getNumParticles();
    totalNumBuckets += tp->getNumBuckets();
  }
  numNodes -= cumNumReplicatedNodes;

  CudaMultipoleMoments *localMoments;
  CompactPartData *localParticles;

  localMoments = new CudaMultipoleMoments[numNodes];
  localParticles = new CompactPartData[numParticles];

  // needed so we know how many particles there are in each bucket
  //int *bmarks = new int[totalNumBuckets+1];

  // fill up postPrefetchMoments with node moments
  int nodeIndex = 0;
  int partIndex = 0;
  int bucketIndex = 0;

  queue.enq(node);
  while(!queue.isEmpty()){
    GenericTreeNode *node = queue.deq();
    NodeType type = node->getType();

    if(type == Empty || type == CachedEmpty){ // skip
      continue;
    }
    else if(type == Bucket){ // Bu, follow pointer
      // copy particles
      int start = node->firstParticle;
      int end = node->lastParticle;
      GravityParticle *gravParts = node->particlePointer;
      NodeKey bucketKey = node->getKey();

      node->nodeArrayIndex = nodeIndex;
      localMoments[nodeIndex] = node->moments;
#ifdef CUDA_PRINT_POST_PREFETCH_LIST
    CkPrintf("node %d: %ld (%s)\n", nodeIndex, node->getKey(), typeString(type));
#endif
      nodeIndex++;
      localPartsOnGpu[bucketKey << 1] = partIndex;
      //bmarks[bucketIndex] = partIndex;
#ifdef CUDA_PRINT_POST_PREFETCH_LIST
      CkPrintf("local bucket %d: %ld (parts: %d-%d)\n", bucketIndex, bucketKey, partIndex, partIndex+(end-start));
#endif
      for(int i = start; i <= end; i++){
        localParticles[partIndex] = gravParts[i-start];
        partIndex++;
      }
      bucketIndex++;
    }
    else if(type == NonLocalBucket){ // NLB
    	node->nodeArrayIndex = nodeIndex;
    	localMoments[nodeIndex] = node->moments;
#ifdef CUDA_PRINT_POST_PREFETCH_LIST
    	CkPrintf("node %d: %ld (%s)\n", nodeIndex, node->getKey(), typeString(type));
#endif
    	nodeIndex++;
    }
    else if(type == Boundary || type == Internal){ // B,I 
    	node->nodeArrayIndex = nodeIndex;
    	localMoments[nodeIndex] = node->moments;
#ifdef CUDA_PRINT_POST_PREFETCH_LIST
    	CkPrintf("node %d: %ld (%s)\n", nodeIndex, node->getKey(), typeString(type));
#endif
    	nodeIndex++;
    	for(int i = 0; i < node->numChildren(); i++){
    		GenericTreeNode *child = node->getChildren(i);
    		queue.enq(child);
    	}
    }
  }// end while queue not empty
  CkAssert(bucketIndex == totalNumBuckets);
  CkAssert(partIndex == numParticles);

  // Transfer moments and particle cores to gpu
  // FIXME - bufferID
  int bufferID = CkMyPe();
  DataManagerTransfer(localMoments, numNodes, localParticles, numParticles, bufferID);

}// end serializeLocal

TreePieceDescriptor *DataManager::findKeyInDescriptors(SFC::Key particleKey){
	for(int i = 0; i < registeredTreePieces.length(); i++){
		if(registeredTreePieces[i].firstParticleKey == particleKey){
			return &registeredTreePieces[i];
		}
	}
	return 0;
}

void DeleteHostMoments(CudaMultipoleMoments *array){
	delete [] array;
}

void DeleteHostParticles(CompactPartData *array){
	delete [] array;
}
#endif

