/** @file ParallelGravity.cpp
 */
 
#include <iostream>

//#include <popt.h>
#include <unistd.h>

#include "ParallelGravity.h"
#include "CacheManager.h"

extern char *optarg;
extern int optind, opterr, optopt;

using namespace std;

int verbosity;
CProxy_TreePiece treeProxy;
bool _cache;
int _cacheLineDepth;
unsigned int _yieldPeriod;
int numIterations=1;
DomainsDec domainDecomposition;
GenericTrees useTree;
CProxy_TreePiece streamingProxy;

Main::Main(CkArgMsg* m) {
	verbosity = 0;
	theta = 0.7;
	numTreePieces = CkNumPes();
	bucketSize = 12;
	_cache = false;
	_cacheLineDepth=1;
	printBinaryAcc=1;
	_yieldPeriod=100000000;
/*
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
		*/
	
	const char *optstring = "vt:p:b:cd:n:z:y:";
	int c;
	while((c=getopt(m->argc,m->argv,optstring))>0){
		if(c == -1){
			break;
		}
		switch(c){
			case 'v':
				verbosity++;
				break;
			case 't':
				theta = atof(optarg);
				break;
			case 'p':
				numTreePieces = atoi(optarg);
				break;
			case 'b':
				bucketSize = atoi(optarg);
				break;
			case 'c':
				_cache = true;
				break;
			case 'd':
				_cacheLineDepth = atoi(optarg);
				break;
			case 'n':
				numIterations = atoi(optarg);
				break;
			case 'z':
				printBinaryAcc = atoi(optarg);
				break;
			case 'y':
				_yieldPeriod = atoi(optarg);
				break;
		};
	}

	// hardcoding some parameters, later may be full options
	domainDecomposition = SFC_dec;
	useTree = Binary_Oct;

	if (verbosity) 
	  ckerr << "yieldPeriod set to " << _yieldPeriod << endl;
	if(_cacheLineDepth <= 0)
		CkAbort("Cache Line depth must be greater than 0");

	const char *fname;
	if(optind  < m->argc){
		fname = m->argv[optind];
	}else{
		ckerr<<"Base file name missing\n";
		CkExit();
	}
	basefilename = fname;
	
	if (verbosity) {
	  ckerr<<"cache "<<_cache<<endl;
	  ckerr<<"cacheLineDepth "<<_cacheLineDepth<<endl;
	  if(printBinaryAcc==1)
	    ckerr<<"particle accelerations to be printed in binary format..."<<endl;
	  else
	    ckerr<<"particles accelerations to be printed in ASCII format..."<<endl;

	  ckerr << "Verbosity level " << verbosity << endl;
	}

	cacheManagerProxy = CProxy_CacheManager::ckNew();

	CProxy_BlockMap myMap=CProxy_BlockMap::ckNew(); 
	CkArrayOptions opts(numTreePieces); 
	opts.setMap(myMap);

	pieces = CProxy_TreePiece::ckNew(numTreePieces,opts);
	treeProxy = pieces;
	
	streamingProxy = pieces;
	//StreamingStrategy* strategy = new StreamingStrategy(10,50);
	//ComlibAssociateProxy(strategy, streamingProxy);

	if(verbosity)
	  ckerr << "Created " << numTreePieces << " pieces of tree" << endl;
	
	CProxy_Main(thishandle).nextStage();
}

void Main::nextStage() {
	double startTime;
	
	piecedata *totaldata = new piecedata;
	
	ckerr << "Loading particles ...";
	startTime = CkWallTimer();
	pieces.load(basefilename, CkCallbackResumeThread());
	ckerr << " took " << (CkWallTimer() - startTime) << " seconds." << endl;

	ckerr << "Building trees ...";
	startTime = CkWallTimer();
	pieces.buildTree(bucketSize, CkCallbackResumeThread());
	ckerr << " took " << (CkWallTimer() - startTime) << " seconds." << endl;
	
	//pieces.report(CkCallbackResumeThread());
	//cerr << "Trees written" << endl;
	//CkExit();
	
	/*for(int i =0; i<numIterations; i++){
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
	
	// DEBUGGING
	//CkStartQD(CkCallback(CkIndex_TreePiece::quiescence(),pieces));

	// the cached walk
	for(int i =0; i<numIterations; i++){
	  ckerr << "Calculating gravity (tree bucket, theta = " << theta << ") ...";
	  startTime = CkWallTimer();
	  pieces.calculateGravityBucketTree(theta, CkCallbackResumeThread());
	  ckerr << " took " << (CkWallTimer() - startTime) << " seconds." << endl;
	  if(i >= 1){
	    startTime = CkWallTimer();
	    pieces.startlb(CkCallbackResumeThread());
	    ckerr<< "Load Balancing step took "<<(CkWallTimer() - startTime) << " seconds." << endl;
	  }
	  
#if COSMO_STATS > 0
	  CkPrintf("Iteration # %d stats\n",i);
	  CkCallback cb(CkCallback::resumeThread);
	    
	  totaldata->setcallback(cb);
	  pieces[0].getPieceValues(totaldata);
	  totaldata = (piecedata *) cb.thread_delay();
	  ckerr << "Total statistics iteration " << i << ":\n Number of MAC checks: " << totaldata->MACChecks << endl;
	  ckerr << " Number of particle-cell interactions: " << totaldata->CellInteractions << endl;
	  ckerr << " Number of particle-particle interactions: " << totaldata->ParticleInteractions << endl;
	  ckerr << " Total mass of the tree : " << totaldata->totalmass << endl;
	  totaldata->reset();
#endif
	}

#if COSMO_DEBUG > 1
	ckerr << "Outputting accelerations ...";
	startTime = CkWallTimer();
	if(printBinaryAcc)
	  pieces[0].outputAccelerations(OrientedBox<double>(), "acc2", CkCallbackResumeThread());
	else
	  pieces[0].outputAccASCII(OrientedBox<double>(), "acc2", CkCallbackResumeThread());
	
	ckerr << " took " << (CkWallTimer() - startTime) << " seconds." << endl;
#endif

#if COSMO_STATS > 0
	ckerr << "Outputting statistics ...";
	startTime = CkWallTimer();
	Interval<unsigned int> dummy;
	
	pieces[0].outputStatistics(dummy, dummy, dummy, dummy, totaldata->totalmass, CkCallbackResumeThread());

	ckerr << " took " << (CkWallTimer() - startTime) << " seconds." << endl;
#endif

#if COSMO_DEBUG > 1
	ckerr << "Outputting relative errors ...";
	startTime = CkWallTimer();
	pieces[0].outputRelativeErrors(Interval<double>(), CkCallbackResumeThread());
	ckerr << " took " << (CkWallTimer() - startTime) << " seconds." << endl;
#endif

	ckerr << "Done." << endl;
	
	ckerr << endl << "******************" << endl << endl; 
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
