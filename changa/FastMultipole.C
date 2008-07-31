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
 * reqID should now refer to a cell instead of a bucket.
 */
/*
 * There are now several possible "returns" from an opencriterion
 * test:
 * 0: Unknown: destination cell needs to be opened?
 * -1: accept source cell as local expansion for destination cell.
 * -3: accept souce cell as softened monopole for destination cell.
 * 1: open source cell.
 * -2: accept as remote expansion for destination cell.
 */

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

