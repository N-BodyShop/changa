/*
 * FMM implementation based on Joachim's implementation in PKDGRAV2M
 *
 * Proposed strategy:
 * 
 * Start walk on the "owner" of the root cell.  For the first "shifts"
 * call a remote entry to start checking on another TreePiece.
 * 
 */

/*
 * Summary of PKDGRAV2M walk:
 *
 * Given Check list (includes periodic replicas) and current cell:
 *
 * Process the check list:
 * 	For each cell on the list:
 *		if local expansion is valid
 *			if softening is valid
 *				Add cell moments to current local expansion
 *			else if cell satisfies monopole
 *				Add to particle list
 *			else undecided
 *		if undecided
 *			if check cell is larger
 *				if not a bucket
 *				     add children to checklist
 *				else
 *					if node-bucket expansion is valid
 *						if softening is valid
 *						     add to cell interact list
 *						else if monopole is valid
 *							add to	particle list
 *						else
 *							add particles
 *					else
 *					        add particles
 *			else
 *				if current cell is bucket
 *					if node-bucket expansion is valid
 *						if softening is valid
 *						     add to cell interact list
 *						else if monopole is valid
 *							add to	particle list
 *						else
 *							open current cell
 *					else
 *					        open current cell
 *		  if undecided
 *			keep on checklist
 *		  else if should be opened
 *			if not a bucket
 *				add children to checklist
 *			else place particles on interact list
 *		  else
 *			add to either local, cell, or particle list
 *       if current cell is not a bucket
 *		current cell is first child
 *		add sibling to checklist
 *		push particle list, cell list, check list, local
 * 							expansion on to stack
 *		shift local expansion
 *	 else
 *		calculate forces on all particles
 *		current cell is "next cell"
 *		pop stack for particle, cell, check, and local expansion
 */

/*
 * Below is an implementation of openCriterion() that checks all the
 * above conditions.
 * There are now several possible "returns" from an opencriterion
 * test.   Here is what doWork() should do based on these returns:
 *
 * ACCEPT_LOCAL: accept source (i.e. "check") cell as local expansion
 * 	for destination cell: use momLocrAddMomr5() to contribute to
 * 	local expansion.
 * ACCEPT_MONO: accept source cell as softened monopole for
 * 	destination cell, place on particle interaction list.
 * ACCEPT_MULTIPOLE: accept as remote expansion for destination cell;
 * 	place on node interaction list.
 * OPEN_CHECK: replace check cell with children on check list;
 * 		continue processing check list.
 * OPEN_MYNODE: if myNode is not a bucket then
 *		  1) push particle interaction list, node interaction
 * list, check list (with first child), and local expansion (shifted
 * to second child) onto stack.  Use momShiftLocr()
 *		  2) add second child to check list
 *	          3) shift local expansion to center of first child
 *		  4) start examining check list with first child as "myNode"
 *              else (a bucket)
 *                1) Calculate the gravity interaction! Use momEvalLocr().
 *		  2) pop the next cell, with lists and expansion off
 * stack, and start examining checklist.
 */

// Below "myNode" is the node for which we are examining the check
// list, accumulating the local expansion and interaction terms.
// "node" is the node under consideration for acceptance.

int openCriterionNodeFMM(Tree::GenericTreeNode *node, // "source" node
			 Tree::GenericTreeNode *myNode,	// "destination" node
                    Vector3D<double> offset,
                    int localIndex // requesting TreePiece
                    ) {
  // Note that some of this could be pre-calculated into an "opening radius"
  double radius = TreeStuff::opening_geometry_factor * node->moments.radius / theta;
  if(radius < node->moments.radius)
      radius = node->moments.radius;

  double myRadius = TreeStuff::opening_geometry_factor * myNode->moments.radius / theta;
  if(myRadius < myNode->moments.radius)
      myRadius = myNode->moments.radius;

  Sphere<double> s(node->moments.cm + offset, radius);
  Sphere<double> myS(myNode->moments.cm, myRadius);

  if(!Space::intersect(s, myS)) {
      if(!openSoftening(node, myNode, offset))
	  return ACCEPT_LOCAL;
      }
  else {
      double monoRadius = TreeStuff::opening_geometry_factor*node->moments.radius/thetaMono;
      Sphere<double> sM(node->moments.cm + offset, monoRadius);
      if(!Space::intersect(myNode->boundingBox, sM))
	  return ACCEPT_MONO;
      }
  
  if(radius > myRadius) {	// Check Cell is larger
      if(node->getType()==Tree::Bucket || node->getType()==Tree::CachedBucket || node->getType()==Tree::NonLocalBucket) {
	  if(Space::intersect(myNode->boundingBox, s)) {
	      return ACCEPT_PARTICLES;
	      }
	  else {
	      if(!openSoftening(node, myNode, offset)) {
		  return ACCEPT_MULTIPOLE;
		  }
	      else {
		  monoRadius = TreeStuff::opening_geometry_factor*node->moments.radius/thetaMono;
		  Sphere<double> sM(node->moments.cm + offset, monoRadius);
		  if(Space::intersect(myNode->boundingBox, sM))
		      return ACCEPT_PARTICLES;
		  else
		      return ACCEPT_MONO;
		  }
	      }
	  }
      else { // not a bucket
	  return OPEN_CHECK;
	  }
      }
  else { // myNode is larger
      if(myNode->getType()==Tree::Bucket || myNode->getType()==Tree::CachedBucket || myNode->getType()==Tree::NonLocalBucket){
	  if(Space::intersect(myNode->boundingBox, s))
	      return OPEN_MYNODE;
	  else {
	      // Well separated, now check softening
	      if(!openSoftening(node, myNode, offset)) {
		  return ACCEPT_MULTIPOLE;
		  }
	      else {      // Open as monopole?
		  monoRadius = TreeStuff::opening_geometry_factor*node->moments.radius/thetaMono;
		  Sphere<double> sM(node->moments.cm + offset, monoRadius);
		  if(Space::intersect(myNode->boundingBox, sM))
		      return OPEN_MYNODE;
		  else
		      return ACCEPT_MONO;
		  }
	      }
	  }
      else
	  return OPEN_MYNODE;
    }
}

int FMMCompute::doWork(GenericTreeNode *node, // Node to test
			  TreeWalk *tw,
			  State *state,
			  int chunk,
			  int reqID, // bucket we are working with
			  bool isRoot, 
			  bool &didcomp, int awi)
{
    // so that we have a quick return in case of empty nodes
    if(node->getType() == Empty || node->getType() == CachedEmpty){
	return DUMP;
	}
    bool open = openCriterion(tp, node, reqID, state);
    }

