/** @file Reductions.h
 */

#include "Reductions.decl.h"
 
extern CkReduction::reducerType growOrientedBox_float;
extern CkReduction::reducerType growOrientedBox_double;

extern CkReduction::reducerType minmax_int;
extern CkReduction::reducerType minmax_float;
extern CkReduction::reducerType minmax_double;

/// My reduction that returns one copy of the contributed CkCallback
extern CkReduction::reducerType callbackReduction;
