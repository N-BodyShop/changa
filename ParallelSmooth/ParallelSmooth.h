/** @file ParallelSmooth.h
 @author Graeme Lufkin (gwl@u.washington.edu)
 */
 
#ifndef PARALLELSMOOTH_H
#define PARALLELSMOOTH_H

#include <string>

#include "Particle.h"
#include "PeriodicBoundaryConditions.h"

#include "RequestResponse.h"

#include "ParallelSmooth.decl.h"

class Main : public Chare {
	/// The filename of the simulation to process.
	string filename;
	/// The array of TreePieces
	CProxy_Smooth_TreePiece smoothTreePieces;
	CProxy_DataManager dataManager;
	CProxy_Sorter sorter;
	/// The total number of particles in the simulation.
	int nbodies;
	/// The number of TreePieces that get created.
	int numTreePieces;
	/// The percent tolerance used in the particle sorting.
	double tolerance;
	/// The filename prefix to use for output files.
	string outputPrefix;
	
	/// The number of neighbors to find the radius of.
	int numNeighbors;
	
	double startTime;
	
	int numDensityBins;
	double minDensity, maxDensity;
	
public:
	
	Main(CkArgMsg* m);

	void loadSortBuild();
	void doCalculations();
	void makeHistogram(CkReductionMsg* m);
	void doStressTest(CkReductionMsg* m);
};

#endif //PARALLELSMOOTH_H
