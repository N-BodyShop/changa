/** \file Smooth_TreePiece.h
 Smooth_TreePiece extends tree behaviors to allow smoothing operations
 as used in SPH.
 \author Graeme Lufkin (gwl@u.washington.edu)
*/

#ifndef SMOOTH_TREEPIECE_H
#define SMOOTH_TREEPIECE_H

#include <string>

#include "PeriodicBoundaryConditions.h"
#include "SPH_Kernel.h"

#include "TreePiece.h"
#include "RequestResponse.h"

#include "ParallelSmooth.decl.h"

class Smooth_TreePiece : public TreePiece {
	/// The periodic boundary conditions for the simulation I'm in
	PeriodicBoundaryConditions<double> pbc;
	
	/// The array of Neighbor_TreePieces I am a member of.
	CProxy_Smooth_TreePiece smoothTreePieces;
	
	/// The number of nearest neighbors to look for.
	int numNeighbors;
	
	/// A table of pending responses.
	std::map<int, PendingResponse> pendingResponses;
	
	/// The serial number of requests I get.
	int mySerialNumber;
	
	/// The number of my particles whose operation has completed.
	int numComplete;
	
	/// A callback to hold on to.
	CkCallback callback;
	
	/// The kernel to use when smoothing
	SPH::Kernel* kernel;
	
	void preOrderTraverse(TreeNode* node, TreeRequest* req, Response* resp);
	void localFirstTraverse(TreeNode* node, TreeRequest* req, Response* resp);
		
public:
	
	Smooth_TreePiece(const PeriodicBoundaryConditions<double> & periodic) : pbc(periodic), smoothTreePieces(thisArrayID), mySerialNumber(0) {
		kernel = new SPH::SplineKernel;
	}
	
	Smooth_TreePiece(CkMigrateMessage* m) { }
	
	~Smooth_TreePiece() {
		delete kernel;
	}
	
	/** Perform the n-nearest neighbor search.
	 The callback will receive a CkReductionMsg containing no data.
	 */
	void findSmoothingRadius(const int n, const CkCallback& cb);
	
	/** Perform a smooth operation on the particles.
	 The operation to perform is specified with an integer.
	 The callback will receive a CkReductionMsg containing no data.
	 */
	void performSmoothOperation(const int opID, const CkCallback& cb);
	
	void handleTreeRequest(TreeRequest* req);
	
	void receiveResponse(Response* resp);
	
	/// Write my calculated values to disk
	void saveInformation(const std::string& prefix, const CkCallback& cb);
};

#endif //SMOOTH_TREEPIECE_H
