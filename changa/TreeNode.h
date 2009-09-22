/** @file TreeNode.h
 This header defines the structures used to make up trees.
 @author Graeme Lufkin (gwl@u.washington.edu)
 @date 01 Aug 2005 - Deleted the class SFCTreeNode in favor of GenericTreeNode, Filippo
 */

#ifndef TREENODE_H
#define TREENODE_H

#include <sstream>

#include "SFC.h"

namespace TreeStuff {

using namespace SFC;

inline std::string keyBits(const Key k, const int numBits) {
	std::ostringstream oss;
	//oss << "N";
	bool ready = false;
	for(int i = 0; i < numBits; i++) {
	  Key k2 = k & (static_cast<Key>(1) << (62 - i));
	  if (ready) oss << (k2 ? 1 : 0);
	  else if (k2 != 0) ready = true;
	}
	return oss.str();
}

// This converts a "radius to the furthest particle" to the standard
// Barnes-Hut cell size.  We use it to keep the definition of theta
// consistent.

const double opening_geometry_factor = 2 / sqrt(3.0);

/// The maximum number of particles in a bucket.
extern int maxBucketSize;
	

} //close namespace TreeStuff

#endif //TREENODE_H
