/** \file DataManager.cpp
 Implementation of the DataManager
 \author Graeme Lufkin (gwl@u.washington.edu)
*/

#include "ParallelGravity.h"
#include "DataManager.h"
#include "Reductions.h"

#ifdef CUDA
#include "cuda_typedef.h"
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

void DataManager::notifyPresence(Tree::GenericTreeNode *root, TreePiece *tp) {
  CmiLock(__nodelock);
  registeredChares.push_back(root);
#ifdef CUDA
  registeredTreePieces.push_back(tp);
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

#ifdef CUDA
void DataManager::donePrefetch(){
  static int done = 0;
  done++;
  if(done == registeredTreePieces.length()){
    done = 0;
    //serializeNodes();
    serializeNodes(root);
    registeredTreePieces.length() = 0;
  }
}

void DataManager::serializeNodes(GenericTreeNode *node){
//void serializeNodes(){
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
    TreePiece *tp = registeredTreePieces[i];
    numNodes += tp->getNodeLookupTable().size();
    numParticles += tp->getNumParticles();
    totalNumBuckets += tp->getNumBuckets();
  }
  numNodes -= cumNumReplicatedNodes;   

  // find out number of particles and nodes cached
  // get them from cache - iterate and count each type
  std::map<CkCacheKey, CkCacheEntry*>::iterator it;

  for(it = cache->begin();it != cache->end(); it++){
    CkCacheEntryType *type = it->second->type;

    if(dynamic_cast<EntryTypeGravityParticle *>(type)){
      numCachedParticles++;
    }
    else if(dynamic_cast<EntryTypeGravityNode *>(type)){
      numCachedNodes++;
    }
    // don't need to count smooth particles yet
  }
  
  CudaMultipoleMoments *postPrefetchMoments = new CudaMultipoleMoments[numNodes+numCachedNodes];
  CompactPartData *postPrefetchParticles = new CompactPartData[numParticles+numCachedParticles];
  // needed so we know how many particles there are in each bucket
  int *bmarks = new int[totalNumBuckets+1];
  
  // fill up postPrefetchMoments with node moments
  int nodeIndex = 0;
  int partIndex = 0;
  int bucketIndex = 0;

  while(!queue.isEmpty()){
    GenericTreeNode *node = queue.deq();
    node->nodeArrayIndex = nodeIndex;
    // put node moments in postPrefetchMoments
    postPrefetchMoments[nodeIndex] = node->moments;
    nodeIndex++;
    
    NodeType type = node->getType();
    if(type == Empty || type == CachedEmpty){ // skip
      continue;
    }
    else if(type == Bucket){ // follow pointer
      // copy particles
      int start = node->firstParticle;
      int end = node->lastParticle;
      GravityParticle *gravParts = node->particlePointer;
      // give particles with 'particleKey' the index 'bucketIndex' 
      NodeKey particleKey = node->getKey();

      localPartsOnGpu[particleKey << 1] = bucketIndex;
      bmarks[bucketIndex] = partIndex;
      for(int i = start; i <= end; i++){
        postPrefetchParticles[partIndex] = gravParts[i-start];
        partIndex++;
      }
      bucketIndex++;

      continue;
    }
    else if(type == CachedBucket || type == NonLocalBucket){ // request cache
      int numParticles = node->lastParticle-node->firstParticle+1;
      // ask cache for particles:
      // FIXME - works only for one chunk 
      NodeKey key = node->getKey();
      ExternalGravityParticle *parts = (ExternalGravityParticle *)
                                            cacheManagerProxy[CkMyPe()].requestDataNoFetch(key, 0);
      if(parts){

        cachedPartsOnGpu[key << 1] = bucketIndex;
        bmarks[bucketIndex] = partIndex;
        cachedPartsOnGpu[key] = bucketIndex;
        for(int i = 0; i < numParticles; i++){
          postPrefetchParticles[partIndex] = parts[i];
          partIndex++;
        }
        bucketIndex++;
      }
      continue;
    }
    // enqueue children if available
    for(int i = 0; i < node->numChildren(); i++){
      GenericTreeNode *child = node->getChildren(i);
      if(child){// available to dm
        queue.enq(child);
      }
      else{ // look in cache
        // FIXME - works only for one chunk 
        CkCacheEntry *cacheEntry = cacheManagerProxy[CkMyPe()].requestCacheEntryNoFetch(node->getChildKey(i), 0);
        if(cacheEntry && cacheEntry->replyRecvd){ // available in cache
          GenericTreeNode *child = (GenericTreeNode *)cacheEntry->data;
          queue.enq(child);
        }
      }
    }// end for each child of node
  }// end while queue not empty  

}// end serializeNodes
#endif

