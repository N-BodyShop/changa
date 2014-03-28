#include <deque>
#include <iostream>
#include <vector>

#include "charm++.h"

using namespace std;


class ScaledORBMapBG{
#ifdef CMK_VERSION_BLUEGENE
	BGLTorusManager *manager;
#else
	bool isVnodeMode;
#endif
	int xsize,ysize,zsize;
	bool *avail;			// array indicating processor availability

        // cluster information
        Cluster *clusterArray;
        int nClusters;
        
	void map(/*Volume <float> &clusterVol,*/ ClusterSection &section, Volume <int> &procVol, int clusters);
        void q_sort(int axis, int left, int right);	// FIXME - improve to bubblesort below threshold number of entries
        void enqueueNeighbors(int x, int y, int z, deque<int> &q, int dist);
	
	public:
        CkVec <TaggedVector3D> tpCentroids;
	
        void assign(int *from, int *clusterWeights, CkVec<int> &to, int numobjs, int numprocs);
        void sortOnAxis(int axis, ClusterSection section);
        
        bool isNeighborOf(int pe1, int pe2);
  	
};


