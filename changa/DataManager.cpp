/** \file DataManager.cpp
 Implementation of the DataManager
 \author Graeme Lufkin (gwl@u.washington.edu)
*/

#include "ParallelGravity.h"
#include "DataManager.h"
#include "Reductions.h"

using std::string;
using std::vector;
using std::map;

DataManager::DataManager(const CkArrayID& treePieceID) {
	treePieces = CProxy_TreePiece(treePieceID);
	localCache = NULL;
	//interp = 0;
}
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
