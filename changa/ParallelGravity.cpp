/** @file ParallelGravity.cpp
 */
 
#include <iostream>

#include <popt.h>

#include "StreamingStrategy.h"

#include "ParallelGravity.h"

using namespace std;

int verbosity;

Main::Main(CkArgMsg* m) {
	verbosity = 0;
	theta = 0.7;
	numTreePieces = CkNumPes();
	bucketSize = 12;
	
	poptOption optionsTable[] = {
		{"verbose", 'v', POPT_ARG_NONE | POPT_ARGFLAG_ONEDASH | POPT_ARGFLAG_SHOW_DEFAULT, 0, 1, "be verbose about what's going on", "verbosity"},
		{"theta", 't', POPT_ARG_DOUBLE | POPT_ARGFLAG_ONEDASH | POPT_ARGFLAG_SHOW_DEFAULT, &theta, 0, "the opening angle to particle-cell interaction", "opening angle"},
		{"pieces", 'p', POPT_ARG_INT | POPT_ARGFLAG_ONEDASH | POPT_ARGFLAG_SHOW_DEFAULT, &numTreePieces, 0, "the number of TreePieces to create", "num TreePieces"},
		{"bucketSize", 'b', POPT_ARG_INT | POPT_ARGFLAG_ONEDASH | POPT_ARGFLAG_SHOW_DEFAULT, &bucketSize, 0, "the maxiumum number of particles in a bucket", "bucketSize"},
		POPT_AUTOHELP
		POPT_TABLEEND
	};
	
	poptContext context = poptGetContext("ParallelGravity", m->argc, const_cast<const char **>(m->argv), optionsTable, 0);
	
	poptSetOtherOptionHelp(context, " [OPTION ...] basefile");
	
	int rc;
	while((rc = poptGetNextOpt(context)) >= 0) {
		switch(rc) {
			case 1: //increase verbosity
				verbosity++;
				break;
		}
	}
	
	if(rc < -1) {
		cerr << "Argument error: " << poptBadOption(context, POPT_BADOPTION_NOALIAS) << " : " << poptStrerror(rc) << endl;
		poptPrintUsage(context, stderr, 0);
		CkExit();
		return;
	}

	const char* fname = poptGetArg(context);
	
	if(fname == 0) {
		cerr << "You must provide a base filename to calculate gravity on" << endl;
		poptPrintUsage(context, stderr, 0);
		CkExit();
		return;
	} else
		basefilename = fname;
	
	poptFreeContext(context);
	
	if(verbosity)
		cerr << "Verbosity level " << verbosity << endl;
	
	//ComlibInstanceHandle cinst = CkGetComlibInstance();
	//StreamingStrategy* strategy = new StreamingStrategy;
	//strategy->enableShortArrayMessagePacking();
	//cinst.setStrategy(strategy);
	
    pieces = CProxy_TreePiece::ckNew(numTreePieces, numTreePieces);
	if(verbosity)
		cerr << "Created " << numTreePieces << " pieces of tree" << endl;
	
	CProxy_Main(thishandle).nextStage();
}

void Main::nextStage() {
	double startTime;
	if(verbosity)
		cerr << "Loading particles ..." << endl;
	startTime = CkWallTimer();
	pieces.load(basefilename, CkCallbackResumeThread());
	if(verbosity > 1)
		cerr << "Main: Loading particles took " << (CkWallTimer() - startTime) << " seconds." << endl;

	if(verbosity)
		cerr << "Building trees ..." << endl;
	startTime = CkWallTimer();
	pieces.buildTree(bucketSize, CkCallbackResumeThread());
	if(verbosity > 1)
		cerr << "Main: Building trees took " << (CkWallTimer() - startTime) << " seconds." << endl;
	
	//pieces.report(CkCallbackResumeThread());
	//cerr << "Trees written" << endl;
	//CkExit();
	
	if(verbosity)
		cerr << "Calculating gravity (direct) ..." << endl;
	startTime = CkWallTimer();
	pieces.calculateGravityDirect(CkCallbackResumeThread());
	if(verbosity > 1)
		cerr << "Main: Calculating gravity took " << (CkWallTimer() - startTime) << " seconds." << endl;
	/*
	if(verbosity)
		cerr << "Outputting accelerations ..." << endl;
	startTime = CkWallTimer();
	pieces[0].outputAccelerations(OrientedBox<double>(), "acc", CkCallbackResumeThread());
	if(verbosity > 1)
		cerr << "Main: Outputting took " << (CkWallTimer() - startTime) << " seconds." << endl;
	*/
	/*
	if(verbosity)
		cerr << "Calculating gravity (tree, theta = " << theta << ") ..." << endl;
	startTime = CkWallTimer();
	pieces.calculateGravityTree(theta, CkCallbackResumeThread());
	if(verbosity > 1)
		cerr << "Main: Calculating gravity took " << (CkWallTimer() - startTime) << " seconds." << endl;
	*/
	
	if(verbosity)
		cerr << "Calculating gravity (tree bucket, theta = " << theta << ") ..." << endl;
	startTime = CkWallTimer();
	pieces.calculateGravityBucketTree(theta, CkCallbackResumeThread());
	if(verbosity > 1)
		cerr << "Main: Calculating gravity took " << (CkWallTimer() - startTime) << " seconds." << endl;
	
	/*
	if(verbosity)
		cerr << "Outputting accelerations ..." << endl;
	startTime = CkWallTimer();
	pieces[0].outputAccelerations(OrientedBox<double>(), "acc2", CkCallbackResumeThread());
	if(verbosity > 1)
		cerr << "Main: Outputting took " << (CkWallTimer() - startTime) << " seconds." << endl;
	*/
	/*
	if(verbosity)
		cerr << "Outputting statistics ..." << endl;
	startTime = CkWallTimer();
	Interval<unsigned int> dummy;
	pieces[0].outputStatistics(dummy, dummy, dummy, dummy, CkCallbackResumeThread());
	if(verbosity > 1)
		cerr << "Main: Outputting took " << (CkWallTimer() - startTime) << " seconds." << endl;
	*/
	
	if(verbosity)
		cerr << "Outputting relative errors ..." << endl;
	startTime = CkWallTimer();
	pieces[0].outputRelativeErrors(Interval<double>(), CkCallbackResumeThread());
	if(verbosity > 1)
		cerr << "Main: Outputting took " << (CkWallTimer() - startTime) << " seconds." << endl;
	
	cerr << "Done." << endl;
	CkExit();
}


#include "ParallelGravity.def.h"

/*
#define CK_TEMPLATES_ONLY
#include "CacheManager.def.h"
#undef CK_TEMPLATES_ONLY
#include "CacheManager.def.h"
*/
