/** \file DataManager.cpp
 Implementation of the DataManager
 \author Graeme Lufkin (gwl@u.washington.edu)
*/

#include "DataManager.h"

using std::string;
using std::vector;
using std::map;

DataManager::DataManager(const CkArrayID& treePieceID) {
	treePieces = CProxy_TreePiece(treePieceID);
	//interp = 0;
}
/*
void DataManager::loadParticles(const string& filename, const int nbodies, const CkCallback& cb) {
	//figure out how many particles this node should read in
	int segmentSize = nbodies / CkNumPes();
	int beginParticleNum = CkNodeFirst(CkMyNode()) * segmentSize;
	int endParticleNum;
	if(CkMyNode() == CkNumNodes() - 1)
		endParticleNum = nbodies;
	else
		endParticleNum = beginParticleNum + CkNodeSize(CkMyNode()) * segmentSize;
	
	int myNumParticles = endParticleNum - beginParticleNum;

	if(verbosity >= 3)
		cout << "N" << CkMyNode() << "P" << CkMyPe() << ": DataManager: I'll take particles " << beginParticleNum << " to " << (endParticleNum - 1) << endl;

	datafilename = filename;
	//open our part of the file
	PartialTipsyFile ptf(datafilename, beginParticleNum, endParticleNum);
	if(!ptf.loadedSuccessfully()) {
		cerr << "N" << CkMyNode() << ": DataManager: Fatal: tipsy file not loaded!" << endl;
		cb.send(0);
		return;
	}
	
	int myNumTreePieces = myTreePieces.size();
	
	if(myNumParticles < myNumTreePieces) {
		cerr << "N" << CkMyNode() << ": DataManager: Fatal: I need at least one particle per TreePiece!" << endl;
		cb.send(0);
		return;
	}
	
	if(verbosity >= 3)
		cout << "N" << CkMyNode() << "P" << CkMyPe() << ": DataManager: I've got " << myNumTreePieces << " TreePieces and " << myNumParticles << " particles" << endl;
	
	int baseNumParticles = myNumParticles / myNumTreePieces;
	int x = myNumParticles % myNumTreePieces;
	
	FullParticle* myParticles = new FullParticle[baseNumParticles + 1];
	int diskOrder = beginParticleNum;
	int particle = 0;
	for(vector<int>::iterator iter = myTreePieces.begin(); iter != myTreePieces.end(); ++iter) {
		int n;
		if((iter - myTreePieces.begin()) < x)
			n = baseNumParticles + 1;
		else
			n = baseNumParticles;
		
		for(int i = 0; i < n; ++i, ++particle) {
			if(particle < ptf.h.nsph)
				myParticles[i] = ptf.gas[particle];
			else if(particle < ptf.h.nsph + ptf.h.ndark)
				myParticles[i] = ptf.darks[particle - ptf.h.nsph];
			else
				myParticles[i] = ptf.stars[particle - ptf.h.nsph - ptf.h.ndark];
			
			myParticles[i].diskOrder = diskOrder++;
		}

		treePieces[*iter].receiveParticles(myParticles, n, cb);
	}
	
	if(diskOrder != endParticleNum)
		cerr << "DataManager: diskOrder not right!" << endl;
	if(particle != myNumParticles)
		cerr << "DataManager: number of particles missed!" << endl;
	
	if(verbosity >= 3)
		cout << "N" << CkMyNode() << "P" << CkMyPe() << ": DataManager: Particles loaded and distributed to TreePieces" << endl;
	
	delete[] myParticles;	
}
*/
void DataManager::acceptCandidateKeys(const Key* keys, const int n, const CkCallback& cb) {
	boundaryKeys.resize(n);
	copy(keys, keys + n, boundaryKeys.begin());
	//tell the TreePieces on this node to evaluate the splitter keys
	for(vector<int>::iterator iter = myTreePieces.begin(); iter != myTreePieces.end(); ++iter)
		treePieces[*iter].evaluateBoundaries(cb);
}

void DataManager::acceptFinalKeys(const Key* keys, const int* responsible, const int* bins, const int n, const CkCallback& cb) {
	boundaryKeys.resize(n);
	copy(keys, keys + n, boundaryKeys.begin());
	responsibleIndex.resize(n - 1);
	copy(responsible, responsible + n - 1, responsibleIndex.begin());
	particleCounts.resize(n - 1);
	copy(bins, bins + n - 1, particleCounts.begin());
	
	contribute(sizeof(CkCallback), &cb, callbackReduction, CkCallback(CkIndex_TreePiece::unshuffleParticles(0), treePieces));
	
	//tell my TreePieces to move the particle data to the responsible chare
	//for(vector<int>::iterator iter = myTreePieces.begin(); iter != myTreePieces.end(); ++iter)
	//	treePieces[*iter].unshuffleParticles(cb);
}
/*
void DataManager::compilePerParticleFunction(const string& code, const CkCallback& cb) {
	delete interp;
	string prefix = "#include \"Particle.h\"\nextern \"C\" void perParticle(FullParticle& p) {\n";
	string suffix = "}\n";
	string includePath = "/host/hermes/home/gwl/projects/structures";
	cout << "DataManager: Compiling code:\n" << (prefix + code + suffix) << endl;
	interp = new CkCppInterpreter((prefix + code + suffix).c_str(), includePath.c_str());
	if(interp->valid()) {
		cout << "DataManager: Compiled successfully" << endl;
		userPerParticleFunction = (perParticleFunction) interp->lookup("perParticle");
		for(vector<int>::iterator iter = myTreePieces.begin(); iter != myTreePieces.end(); ++iter)
			treePieces[*iter].applyPerParticleFunction(cb);
	} else {
		int i = 1;
		//contribute(sizeof(int), &i, CkReduction::sum_int, cb);
		cout << CkMyNode() << " DataManager: Could not compile correctly." << endl;
	}
}

void DataManager::compileTreeRequestFunction(const string& code, const CkCallback& cb) {
}
*/
