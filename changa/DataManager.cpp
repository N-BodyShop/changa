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
	//interp = 0;
}
void DataManager::acceptCandidateKeys(const SFC::Key* keys, const int n, const CkCallback& cb) {
	boundaryKeys.resize(n);
	copy(keys, keys + n, boundaryKeys.begin());
	//tell the TreePieces on this node to evaluate the splitter keys
	for(vector<int>::iterator iter = myTreePieces.begin(); iter != myTreePieces.end(); ++iter)
		treePieces[*iter].evaluateBoundaries(cb);
}

void DataManager::acceptFinalKeys(const SFC::Key* keys, const int* responsible, const int* bins, const int n, const CkCallback& cb) {
	boundaryKeys.resize(n);
	copy(keys, keys + n, boundaryKeys.begin());
	responsibleIndex.resize(n - 1);
	copy(responsible, responsible + n - 1, responsibleIndex.begin());
	particleCounts.resize(n - 1);
	copy(bins, bins + n - 1, particleCounts.begin());

  if(verbosity >= 3 && CkMyPe()==0){
    std::vector<int>::iterator iter;
    CkPrintf("responsible:");
    for(iter=responsibleIndex.begin();iter!=responsibleIndex.end();iter++){
      CkPrintf("%d,",*iter);
    }
    CkPrintf("\n");
    CkPrintf("Particle Counts:");
    for(iter=particleCounts.begin();iter!=particleCounts.end();iter++){
      CkPrintf("%d,",*iter);
    }
    CkPrintf("\n");
    std::vector<SFC::Key>::iterator iter2;
    CkPrintf("Keys:");
    for(iter2=boundaryKeys.begin();iter2!=boundaryKeys.end();iter2++){
      CkPrintf("%016llx,",*iter2);
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
