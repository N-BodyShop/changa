#include <deque>
#include <iostream>
#include <vector>

#include "charm++.h"
#include "BGLTorus.h"

using namespace std;

class ScaleTranMapBG{
#ifdef CMK_VERSION_BLUEGENE
	BGLTorusManager *manager;
#else
	bool isVnodeMode;
#endif
	int xsize,ysize,zsize;
        int collisions;
	Cluster *clusterArray;
	int nClusters;
	
	inline void translate(Vector3D<float> &vec, float xdist, float ydist, float zdist);
	inline void scale(Vector3D<float> &vec, int xlen, int ylen, int zlen);
	void roundOff(Vector3D<float> &vec);

	int map(Vector3D<float> &vec, int numprocs, bool *avail, int *ringRadius);
	void enqueueNeighbors(int x, int y, int z, deque<int> &q, int dist);
	
	public:
	void assign(int *from, double *clusterWeights, CkVec<int> &to, int numobjs, int numprocs);
        CkVec <TaggedVector3D> tpCentroids;
        bool isNeighborOf(int pe1, int pe2);
};
