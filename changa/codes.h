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
// opt. object doesn't know enough about the context of the computation
// this happens when we have to check whether a cached node's ancestors 
// haven't already been involved in a computation - an ancestor check must
// be performed, based on which the compute obj. makes a decision 
#define DEFER 5 
           
#define INTERSECT -1
#define CONTAIN 1 
#define NO_INTERSECT 0

enum WalkType {TopDown, LocalTarget, BottomUp, BucketIterator, InvalidWalk};
enum ComputeType {Gravity, Prefetch, List, BucketEwald, Smooth, ReSmooth,
		  InvalidCompute};

enum OptType {Local, Remote, Pref, Double, InvalidOpt};

#define INTERLIST_LEVELS 64
#define NUM_NODE_TYPES 11

// debug
#ifdef CHANGA_REFACTOR_WALKCHECK
#define CHECK_INDEX 1
#define CHECK_BUCKET 16
#endif

#define TEST_BUCKET 0
#endif // CODES_H
