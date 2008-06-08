/* values returned by compute object to tree walk object */
#ifndef CODES_H
#define CODES_H

#define KEEP 0  // descend further down this path
#define KEEP_REMOTE_BUCKET 6    // request particles of a remote bucket
#define KEEP_LOCAL_BUCKET 7     // directly access particles of a local bucket
#define DUMP 1  // no need to descend down this path  
#define ERROR 2

/* values, in addition to the ones above 
 * returned by optimzation object to compute object */
#define NOP 3 
#define COMPUTE 4
#define DEFER 5 // opt. object doesn't know enough about the context of the computation
		// this happens when we have to check whether a cached node's ancestors 
		// haven't already been involved in a computation - an ancestor check must
		// be performed, based on which the compute obj. makes a decision 
           
enum WalkType {TopDown, Double, BottomUp, BucketIterator, InvalidWalk};
enum ComputeType {Gravity, Prefetch, List, BucketEwald, Smooth, ReSmooth,
		  InvalidCompute};

enum OptType {Local, Remote, Pref, InvalidOpt};


#ifdef CHANGA_REFACTOR_WALKCHECK
#define CHECK_INDEX 1
#define CHECK_BUCKET 16
#endif

#endif // CODES_H
