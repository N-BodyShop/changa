/** \file GravityTreeNode.h
 This file defines a tree node suitable for doing gravity calculations on.
 @author Graeme Lufkin (gwl@u.washington.edu)
 @version 1.0
 */

#ifndef GRAVITYTREENODE_H
#define GRAVITYTREENODE_H

#include "GenericTreeNode.h"
#include "MultipoleMoments.h"

namespace Tree {

/// A gravity tree node has multipole moments of the enclosed mass
/* Note that GravityTreeNode does not need to be a BinaryTreeNode.  This
 was chosen for purely binary-biased reasons.  We could have gravity
 tree nodes of all types that simply add the multipoles to the node. */
class GravityTreeNode : public BinaryTreeNode {
public:
		
	MultipoleMoments moments;
};

} //close namespace Tree

#endif //GRAVITYTREENODE_H
