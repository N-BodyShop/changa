/** @file ParallelGravity.cpp
 */
 
#include <iostream>

#include <popt.h>

#include "StreamingStrategy.h"

#include "ParallelGravity.h"
#include "CacheManager.h"

using namespace std;

int verbosity;
CProxy_TreePiece treeProxy;
bool _cache;
int _cacheLineDepth;
int numIterations=1;

Main::Main(CkArgMsg* m) {
	verbosity = 0;
	theta = 0.7;
	numTreePieces = CkNumPes();
	bucketSize = 12;
	_cache = false;
	_cacheLineDepth=1;

	poptOption optionsTable[] = {
		{"verbose", 'v', POPT_ARG_NONE | POPT_ARGFLAG_ONEDASH | POPT_ARGFLAG_SHOW_DEFAULT, 0, 1, "be verbose about what's going on", "verbosity"},
		{"theta", 't', POPT_ARG_DOUBLE | POPT_ARGFLAG_ONEDASH | POPT_ARGFLAG_SHOW_DEFAULT, &theta, 0, "the opening angle to particle-cell interaction", "opening angle"},
		{"pieces", 'p', POPT_ARG_INT | POPT_ARGFLAG_ONEDASH | POPT_ARGFLAG_SHOW_DEFAULT, &numTreePieces, 0, "the number of TreePieces to create", "num TreePieces"},
		{"bucketSize", 'b', POPT_ARG_INT | POPT_ARGFLAG_ONEDASH | POPT_ARGFLAG_SHOW_DEFAULT, &bucketSize, 0, "the maxiumum number of particles in a bucket", "bucketSize"},
		{"cache", 'c', POPT_ARG_NONE | POPT_ARGFLAG_ONEDASH | POPT_ARGFLAG_SHOW_DEFAULT, &_cache, 1,"should the cache manager be used"},
		{"cacheLineDepth",'d', POPT_ARG_INT | POPT_ARGFLAG_ONEDASH | POPT_ARGFLAG_SHOW_DEFAULT,&_cacheLineDepth,0,"The depth of the cacheLine tree to be fetched ","cache Line Depth"},
		{"numIterations",'n', POPT_ARG_INT | POPT_ARGFLAG_ONEDASH | POPT_ARGFLAG_SHOW_DEFAULT,&numIterations,0,"Number of Iterations of the walk bucket tree algorithm ","numIterations"},
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
	cerr<<"cache "<<_cache<<endl;
	cerr<<"cacheLineDepth "<<_cacheLineDepth<<endl;
	
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
	
	ComlibInstanceHandle cinst = CkGetComlibInstance();
	StreamingStrategy* strategy = new StreamingStrategy(10,50);
	cinst.setStrategy(strategy);

	cacheManagerProxy = CProxy_CacheManager::ckNew();
    pieces = CProxy_TreePiece::ckNew(numTreePieces, numTreePieces);
	treeProxy = pieces;
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
	/*
	for(int i =0; i<numIterations; i++){
		if(verbosity)
			cerr << "Calculating gravity (direct) ..." << endl;
		startTime = CkWallTimer();
		pieces.calculateGravityDirect(CkCallbackResumeThread());
		if(verbosity > 1)
			cerr << "Main: Calculating gravity took " << (CkWallTimer() - startTime) << " seconds." << endl;
	}
	*/
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
	
	// the cached walk
	for(int i =0; i<numIterations; i++){
		if(verbosity)
			cerr << "Calculating gravity (tree bucket, theta = " << theta << ") ..." << endl;
		startTime = CkWallTimer();
		pieces.calculateGravityBucketTree(theta, CkCallbackResumeThread());
		if(verbosity)
			cerr << "Main: Calculating gravity took " << (CkWallTimer() - startTime) << " seconds." << endl;
		startTime = CkWallTimer();
		if(i > 2){
			pieces.startlb(CkCallbackResumeThread());
		}
		if(verbosity){
			cerr<< "Main: Load Balancing step took "<<(CkWallTimer() - startTime) << " seconds." << endl;
		}
	}
	
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
	/*
	if(verbosity)
		cerr << "Outputting relative errors ..." << endl;
	startTime = CkWallTimer();
	pieces[0].outputRelativeErrors(Interval<double>(), CkCallbackResumeThread());
	if(verbosity > 1)
		cerr << "Main: Outputting took " << (CkWallTimer() - startTime) << " seconds." << endl;
	
	cerr << "Done." << endl;
	*/
	CkExit();
}


#include "ParallelGravity.def.h"
#include "CacheManager.def.h"

/*
#define CK_TEMPLATES_ONLY
#include "CacheManager.def.h"
#undef CK_TEMPLATES_ONLY
#include "CacheManager.def.h"
*/
