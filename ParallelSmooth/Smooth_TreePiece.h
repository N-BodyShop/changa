/** @file Smooth_TreePiece.h
 Smooth_TreePiece extends tree behaviors to allow smoothing operations
 as used in SPH.
 @author Graeme Lufkin (gwl@u.washington.edu)
*/

#ifndef SMOOTH_TREEPIECE_H
#define SMOOTH_TREEPIECE_H

#include <string>

#include "PeriodicBoundaryConditions.h"
#include "SPH_Kernel.h"

#include "TreePiece.h"
#include "RequestResponse.h"


enum SmoothOperation { 
	Density,
	VelocityDivCurl,
	BallCount,
	VelocityDispersion,
	User
};

#include "ParallelSmooth.decl.h"

class Smooth_TreePiece : public TreePiece {
	/// The periodic boundary conditions for the simulation I'm in
	PeriodicBoundaryConditions<double> pbc;
	
	/// The array of \c Neighbor_TreePieces I am a member of.
	CProxy_Smooth_TreePiece smoothTreePieces;
	
	/// The number of nearest neighbors to look for.
	int numNeighbors;
	
	/// A table of pending responses.
	std::map<int, PendingResponse> pendingResponses;
	
	/// The serial number of requests I get.
	int mySerialNumber;
	
	/// The number of my particles whose operation has completed.
	int numComplete;
	/// The number of particles that are being smoothed in this operation
	int totalNumPending;
	
	/// A callback to hold on to.
	CkCallback callback;
	
	void preOrderTraverse(TreeNode* node, TreeRequest* req, Response* resp);
	void localFirstTraverse(TreeNode* node, TreeRequest* req, Response* resp);
	void localOnlyTraverse(TreeNode* node, TreeRequest* req, Response* resp);
	void remoteOnlyTraverse(TreeNode* node, TreeRequest* req, Response* resp);
	
	int numDensityBins;
	
public:
	
	Smooth_TreePiece(const PeriodicBoundaryConditions<double> & periodic) : pbc(periodic), smoothTreePieces(thisArrayID), mySerialNumber(0) { }
	
	Smooth_TreePiece(CkMigrateMessage* m) { }
	
	~Smooth_TreePiece() { }
	
	/** Perform the n-nearest neighbor search.
	 The callback will receive a \c CkReductionMsg containing no data.
	 */
	void findSmoothingRadius(const int n, const CkCallback& cb);
	
	/** Perform a smooth operation on the particles.
	 The operation to perform is specified with an integer.
	 The callback will receive a \c CkReductionMsg containing no data.
	 */
	void performSmoothOperation(const SmoothOperation op, const CkCallback& cb);
	
	void handleTreeRequest(TreeRequest* req);
	
	void receiveResponse(Response* resp);
	
	/// Write my calculated values to disk
	void saveInformation(const std::string& prefix, const CkCallback& cb);
	
	/// Find the minimum and maximum particle densities
	void minmaxDensity(const CkCallback& cb);
	void makeDensityHistogram(const int numDensityBins, const double minDensity, const double maxDensity, const CkCallback& cb);

	/// Do a smoothing operation with a density criterion
	void densityCutOperation(const SmoothOperation op, double minDensity, const CkCallback& cb);
};

#endif //SMOOTH_TREEPIECE_H
