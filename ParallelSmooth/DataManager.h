/** \file DataManager.h
 The DataManager holds data that TreePieces need access to, but not
 their own copy.
 \author Graeme Lufkin (gwl@u.washington.edu)
*/

#ifndef DATAMANAGER_H
#define DATAMANAGER_H

#include <vector>
#include <map>
#include <string>

//#include "ckdll.h"

#include "NChilada.h"

#include "Particle.h"

#include "RequestResponse.h"

//typedef void (*perParticleFunction)(FullParticle&);
//typedef TreeRequest* (*treeRequestCreationFunction)(const FullParticle&);

/** The DataManager is used to store information that all TreePieces will need,
 but will not modify.  The first example is the list of splitter keys and the
 responsible chare for each interval.  This data is given to the DataManager by
 the Sorter.  The DataManager then instructs all the TreePieces on its node
 to evaluate the boundary keys.
 The DataManager also reads in all the particles for its node, and then hands
 them out to the TreePieces who will handle them.
 */
class DataManager : public NodeGroup {
	friend class TreePiece;
	/// The array of TreePieces I hold data for.
	CProxy_TreePiece treePieces;

protected:
	
	/// The filename of the tipsy file.
	std::string datafilename;
	
	/// An array telling me which TreePieces are actually on my node.
	std::vector<int> myTreePieces;
	
	/// The array of splitter keys for the sort.
	std::vector<Key> boundaryKeys;
	
	/// An array identifying which chare is responsible for each interval of keys.
	std::vector<int> responsibleIndex;
	
	/// An array with how many particles are held by each TreePiece when sorted.
	std::vector<int> particleCounts;
	
	/// A pointer to a user defined per-particle function.
	//perParticleFunction userPerParticleFunction;
	
	/// A pointer to a used defined operation that creates a TreeRequest.
	//treeRequestCreationFunction createUserTreeRequest;
	
	/// The code interpreter on this node.
	//CkCppInterpreter* interp;
	
public:
	
	DataManager(const CkArrayID& treePieceID);
	
	~DataManager() { }
	
	/** Read my node's share of particles, and distribute them to my TreePieces.
	 To use, you must have created a DataManager and array of TreePieces and tied them together
	 with treePieces.registerWithDataManager()
	 You need to know the total number of particles in the simulation.
	 Example:
		PartialTipsyFile ptf(filename, 0, 1);
		dataManager.loadParticles(filename, ptf.fullHeader.nbodies, CkCallback(...))
	 When the callback returns, the particles will be loaded into the TreePieces, given
	 keys relative to the bounding box, and locally sorted.  Global sorting takes place
	 with a Sorter.  The callback will receive a CkReductionMsg containing no data.
	 */
	void loadParticles(const string& filename, const int nbodies, const CkCallback& cb);
	
	/// Called by the Sorter, I ask my TreePieces to evaluate these splitter keys
	void acceptCandidateKeys(const Key* keys, const int n, const CkCallback& cb);
	/// Called by the Sorter, I save these final keys and the list of which TreePiece is responsible for which interval
	void acceptFinalKeys(const Key* keys, const int* responsible, const int* bins, const int n, const CkCallback& cb);
	
	//void compilePerParticleFunction(const string& code, const CkCallback& cb);
	
	//void compileTreeRequestFunction(const string& code, const CkCallback& cb);
};

#endif //DATAMANAGER_H
