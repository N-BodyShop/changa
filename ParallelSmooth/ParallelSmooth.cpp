/** \file ParallelSmooth.cpp
 \author Graeme Lufkin (gwl@u.washington.edu)
*/

#include <iostream>
#include <popt.h>

#include "pup_stl.h"
//#include "liveViz.h"

#include "SPH_Kernel.h"

#include "NChilada.h"
#include "ParallelSmooth.h"
#include "Smooth_TreePiece.h"
#include "SPH_Requests.h"

using std::string;
using std::ofstream;

Main::Main(CkArgMsg* m) {
	
	double allPeriod = 0;
	verbosity = 0;
	char* prefix = 0;
	tolerance = 0.01;
	numNeighbors = 32;
	
	poptOption optionsTable[] = {
		{"neighbors", 'n', POPT_ARG_INT | POPT_ARGFLAG_ONEDASH | POPT_ARGFLAG_SHOW_DEFAULT, &numNeighbors, 0, "number of nearest neighbors to find radius of", "num neighbors"},
		{"xp", 'x', POPT_ARG_DOUBLE | POPT_ARGFLAG_ONEDASH | POPT_ARGFLAG_SHOW_DEFAULT, &pbc.xPeriod, 0, "x-axis period for periodic boundary conditions", "xPeriod"},
		{"yp", 'y', POPT_ARG_DOUBLE | POPT_ARGFLAG_ONEDASH | POPT_ARGFLAG_SHOW_DEFAULT, &pbc.yPeriod, 0, "y-axis period for periodic boundary conditions", "yPeriod"},
		{"zp", 'z', POPT_ARG_DOUBLE | POPT_ARGFLAG_ONEDASH | POPT_ARGFLAG_SHOW_DEFAULT, &pbc.zPeriod, 0, "z-axis period for periodic boundary conditions", "zPeriod"},
		{"period", 'p', POPT_ARG_DOUBLE | POPT_ARGFLAG_ONEDASH | POPT_ARGFLAG_SHOW_DEFAULT, &allPeriod, 2, "period for periodic boundary conditions (all dimensions)", "period"},
		{"verbose", 'v', POPT_ARG_NONE | POPT_ARGFLAG_ONEDASH | POPT_ARGFLAG_SHOW_DEFAULT, 0, 1, "be verbose about what's going on", "verbosity"},
		{"output", 'o', POPT_ARG_STRING | POPT_ARGFLAG_ONEDASH, &prefix, 0, "prefix for output filenames", "string"},
		{"tolerance", 't', POPT_ARG_DOUBLE | POPT_ARGFLAG_ONEDASH | POPT_ARGFLAG_SHOW_DEFAULT, &tolerance, 0, "percent tolerance for sorting into TreePieces", "percent tolerance"},
		POPT_AUTOHELP
		POPT_TABLEEND
	};
	
	poptContext context = poptGetContext("Neighbors", m->argc, const_cast<const char **>(m->argv), optionsTable, 0);
	
	poptSetOtherOptionHelp(context, " [OPTION ...] filename");
	
	int rc;
	while((rc = poptGetNextOpt(context)) >= 0) {
		switch(rc) {
			case 1: //increase verbosity
				verbosity++;
				break;
			case 2: //allPeriod sets x-, y-, and z-period
				pbc.xPeriod = pbc.yPeriod = pbc.zPeriod = allPeriod;
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
	
	poptFreeContext(context);
	
	if(fname == 0) {
		cerr << "You must provide a filename to find groups in!" << endl;
		poptPrintUsage(context, stderr, 0);
		CkExit();
		return;
	} else
		filename = fname;
		
	if(pbc.xPeriod < 0 || pbc.yPeriod < 0 || pbc.zPeriod < 0) {
		cerr << "Periods for periodic boundary conditions must be positive!" << endl;
		CkExit();
		return;
	}

	if(verbosity)
		cout << "Arguments processed, starting main." << endl;
	
	//figure out how many particles in the whole file
	PartialTipsyFile ptf(filename, 0, 1);
	if(!ptf.loadedSuccessfully()) {
		cerr << "Couldn't load the tipsy file \"" << filename << "\". Maybe it's not a tipsy file?" << endl;
		CkExit();
		return;
	}
	
	nbodies = ptf.fullHeader.nbodies;
	
	if(numNeighbors < 2 || numNeighbors > nbodies) {
		cerr << "Number of neighbors must be positive, and smaller than the total number of particles!" << endl;
		CkExit();
		return;
	}


	if(prefix == 0) {
		string s(filename);
		unsigned int pos = s.find_last_of("/");
		if(pos != string::npos)
			s.erase(0, pos + 1);
		outputPrefix = s.c_str();
	} else
		outputPrefix = prefix;
	
	cout << "The full tipsy file has " << nbodies << " particles." << endl;
	cout << "The number of nearest neighbors to look at is " << numNeighbors << endl;
	cout << "Periodic boundary conditions: x: " << pbc.xPeriod << " y: " << pbc.yPeriod << " z: " << pbc.zPeriod << endl;
	cout << "Output files will be made with prefix \"" << outputPrefix << "\"" << endl;
	if(verbosity)
		cout << "Verbosity level " << verbosity << endl;
	
	//initializeNChiladaCcs();
	
	//create an array of TreePieces (XXX: how many should we make?)
	numTreePieces = min(CkNumPes() * 6, ptf.fullHeader.nbodies / 200);
	smoothTreePieces = CProxy_Smooth_TreePiece::ckNew(pbc, numTreePieces);
	treePieceID = smoothTreePieces;
	
	//create the DataManager
	dataManager = CProxy_DataManager::ckNew(smoothTreePieces);
	dataManagerID = dataManager;
	
	/*
	CkBbox3d box;
	box.min = CkVector3d(-0.5, -0.5, -0.5);
	box.max = CkVector3d(0.5, 0.5, 0.5);
	liveVizConfig cfg(false, false, box);
	liveVizInit(cfg, smoothTreePieces, CkCallback(CkIndex_TreePiece::generateImage(0), smoothTreePieces));
	*/
	CProxy_Main(thishandle).loadSortBuild();
}

void Main::loadSortBuild() {
	if(verbosity)
		cout << "Main: Registering Neighbor_TreePieces" << endl;
	smoothTreePieces.registerWithDataManager(dataManager, CkCallbackResumeThread());
	if(verbosity)
		cout << "Main: Loading particles" << endl;
	dataManager.loadParticles(filename, nbodies, CkCallbackResumeThread());
	sorter = CProxy_Sorter::ckNew();
	if(verbosity)
		cout << "Main: Sorting particles" << endl;
	startTime = CkTimer();
	sorter.startSorting(dataManager, numTreePieces, tolerance, CkCallbackResumeThread());
	if(verbosity > 1)
		cout << "Main: Sorting took " << (CkTimer() - startTime) << " seconds." << endl;
	if(verbosity)
		cout << "Main: Building trees" << endl;
	startTime = CkTimer();
	smoothTreePieces.startTreeBuild(CkCallbackResumeThread());
	if(verbosity > 1)
		cout << "Main: Tree build took " << (CkTimer() - startTime) << " seconds." << endl;
	CProxy_Main(thishandle).doCalculations();
}

void Main::doCalculations() {
	if(verbosity)
		cout << "Main: Starting neighbor search" << endl;
	startTime = CkTimer();
	smoothTreePieces.findSmoothingRadius(numNeighbors, CkCallbackResumeThread());
	if(verbosity > 1)
		cout << "Main: Radius finding took " << (CkTimer() - startTime) << " seconds." << endl;
	if(verbosity)
		cout << "Main: Starting density calculation" << endl;
	startTime = CkTimer();
	smoothTreePieces.performSmoothOperation(DensityOperation, CkCallbackResumeThread());
	if(verbosity > 1)
		cout << "Main: Density calculation took " << (CkTimer() - startTime) << " seconds." << endl;
	if(verbosity)
		cout << "Main: Starting velocity calculations" << endl;
	startTime = CkTimer();
	smoothTreePieces.performSmoothOperation(VelDivCurlOperation, CkCallbackResumeThread());
	if(verbosity > 1)
		cout << "Main: Velocity div/curl/mean calculation took " << (CkTimer() - startTime) << " seconds." << endl;
	if(verbosity)
		cout << "Main: Starting velocity dispersion calculation" << endl;
	startTime = CkTimer();
	smoothTreePieces.performSmoothOperation(VelDispOperation, CkCallbackResumeThread());
	if(verbosity > 1)
		cout << "Main: Velocity dispersion calculation took " << (CkTimer() - startTime) << " seconds." << endl;
	
	ofstream outgroup((string(outputPrefix) + ".radius").c_str());
	outgroup.write(&nbodies, sizeof(int));
	outgroup.close();
	outgroup.open((string(outputPrefix) + ".density").c_str());
	outgroup.write(&nbodies, sizeof(int));
	outgroup.close();
	outgroup.open((string(outputPrefix) + ".divv").c_str());
	outgroup.write(&nbodies, sizeof(int));
	outgroup.close();
	outgroup.open((string(outputPrefix) + ".curlv").c_str());
	outgroup.write(&nbodies, sizeof(int));
	outgroup.close();
	outgroup.open((string(outputPrefix) + ".meanv").c_str());
	outgroup.write(&nbodies, sizeof(int));
	outgroup.close();
	outgroup.open((string(outputPrefix) + ".veldisp").c_str());
	outgroup.write(&nbodies, sizeof(int));
	outgroup.close();
	outgroup.open((string(outputPrefix) + ".phase").c_str());
	outgroup.write(&nbodies, sizeof(int));
	outgroup.close();
	outgroup.open((string(outputPrefix) + ".mach").c_str());
	outgroup.write(&nbodies, sizeof(int));
	outgroup.close();
	outgroup.open((string(outputPrefix) + ".speed").c_str());
	outgroup.write(&nbodies, sizeof(int));
	outgroup.close();
	if(verbosity)
		cout << "Main: Writing values to file" << endl;
	smoothTreePieces.saveInformation(outputPrefix, CkCallbackResumeThread());
	
	//cout << "Waiting for CCS control." << endl;
	CkExit();
}

extern void initRequestResponsePUP();
extern void initSPHPUP();

#include "ParallelSmooth.def.h"
