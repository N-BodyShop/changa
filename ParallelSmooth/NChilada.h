/** \file NChilada.h
 \author Graeme Lufkin (gwl@u.washington.edu)
*/

#ifndef NCHILADA_H
#define NCHILADA_H

#include <string>

#include "liveViz.h"

#include "Particle.h"

#include "NChilada.decl.h"

/** The verbosity level, used to control how much diagnostic information gets output.
 Rule of thumb usage is: Level 1 explains program flow.  Level 2 gives timing information.
 Level 3 gives other useful information.  Level 4 and up gives probably useless stats. */
extern int verbosity;

/// My reduction that grows axis-aligned boxes.
extern CkReduction::reducerType boxGrowthReduction;

/// The chare ID of the array of TreePieces you create.  You must set this!
extern CkArrayID treePieceID;

/// The group ID of your DataManager.  You must set this!
extern CkGroupID dataManagerID;

/// Register the CCS handlers, this must be called in your mainchare's constructor.
//void initializeNChiladaCcs();

#endif //NCHILADA_H
