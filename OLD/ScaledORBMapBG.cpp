#include "MapStructures.h"
#include "ScaledORBMapBG.h"

#include "ParallelGravity.h"

void ScaledORBMapBG::assign(int *from, int *clusterWeights, CkVec<int> &to, int numobjs, int count){
        int i;

	// Get number of clusters
        nClusters = count;
	clusterArray = new Cluster[nClusters];
        
        float *clusterRadius = new float[nClusters];    // radii
        for(i = 0; i < count; i++){
          clusterRadius[i] = 0.0;
        }
	// Get cluster centroids and cluster counts by going through 'from'
	// from[i] gives the cluster of treepiece i
	Vector3D<float> vec;
        int j;
	for(i = 0; i < numobjs; i++)	{
          j = (int)tpCentroids[i].tag; 
          clusterArray[from[j]].count++;
          clusterArray[from[j]].centroid[0] += tpCentroids[j].vec.x;
          clusterArray[from[j]].centroid[1] += tpCentroids[j].vec.y;
          clusterArray[from[j]].centroid[2] += tpCentroids[j].vec.z;
	  clusterArray[from[j]].tps.push_back(i);
	}

    	for(i = 0; i < nClusters; i++)
            for(int k = 0; k < 3; k++)
	      clusterArray[i].centroid[k] /= clusterArray[i].count;
        
	// we have the centers of mass of the clusters at this point
		
#define abs(x) (x > 0 ? x : -1*x)
 
#ifdef CMK_VERSION_BLUEGENE
	manager = new BGLTorusManager;
	xsize = manager->getDimNX();
	ysize = manager->getDimNY();
	zsize = manager->getDimNZ();
#else
	isVnodeMode = false;
	xsize = 2; 
	ysize = 2; 
	zsize = 2; 
#endif
	// we now have the torus dimensions
        int l[3] = {0,0,0};
        int g[3] = {xsize-1,ysize-1,zsize-1};
        Volume<int> entireProcSpace(l,g);       
        // all clusters are assigned to the entire 
        // processor space to begin with
        // the idea is to recursively refine this 
        // assignment to successively smaller 
        // partitions of processors, till each 
        // cluster has just one processor

        // construct initial simulation space
        float lesser[3] = {-.5,-.5,-.5};
        float greater[3] = {.5,.5,.5};
        Volume<float> clusterVol(lesser, greater);
        
        // create initial array section - 
        // entire array of clusters
        ClusterSection section(0,nClusters-1);
        
#if COSMO_MCLB > 1 
        CkPrintf("Map \n{");
        for(i = 0; i < nClusters; i++)
          CkPrintf("(%f,%f,%f), ", clusterArray[i].centroid[0],clusterArray[i].centroid[1],clusterArray[i].centroid[2]);
        CkPrintf("}\n");
#endif
        map(section, entireProcSpace, nClusters);
                
	// write out the assignments in the 'to' array
        for(i = 0; i < numobjs; i++){
          to[i] = clusterArray[from[i]].proc;
#if COSMO_MCLB > 1 
          CkPrintf("ScaledORBMapBG.C: TP %d to proc. %d\n", i, to[i]);
#endif  
          if (to[i] == -1)
            CmiAbort("ScaledORBMapBG.C: assignment out of range\n");
        }
        
        delete [] clusterRadius;
        delete [] clusterArray;
}

void ScaledORBMapBG::map(ClusterSection &section, Volume <int> &procVol, int totalClusters){
 
  int i;

#if COSMO_MCLB > 1 
  CkPrintf("map: (%d %d) %d %d %d %d %d %d\n", section.beg, section.end, procVol.lesser[0],procVol.lesser[1],procVol.lesser[2],procVol.greater[0],procVol.greater[1],procVol.greater[2]);
#endif 
  if(totalClusters == 1){      
    CmiAssert(section.end == section.beg);
    CmiAssert(procVol.lesser[0] == procVol.greater[0]);
    CmiAssert(procVol.lesser[1] == procVol.greater[1]);
    CmiAssert(procVol.lesser[2] == procVol.greater[2]);
#ifdef CMK_VERSION_BLUEGENE
    clusterArray[section.beg].proc = manager->coordinatesToRank(procVol.lesser[0], procVol.lesser[1], procVol.lesser[2]);
#else
    clusterArray[section.beg].proc = procVol.lesser[0] + xsize*procVol.lesser[1]+ xsize*ysize*procVol.lesser[2];
#endif
    return;
  }
  // decide which axis to split on based on partition dimensions
  int axis = procVol.decideSplitAxis();
  
  /* split the processor space along the same dimension
    FIXME - tacit assumption that number of processors (and therefore the number of clusters) is power of two. This might not happen, 
    and so, to avoid load imbalance, we must divide the processors in the ratio determined by the loads
    induced by the corresponding partitions of clusters, and not down the middle. Need to have sheared partitions - do later.
    
    Also, currently, we split the simulation space into partitions based on the number of clusters (half to each,) and not the
    actual load induced by the two sets of clusters  
  */
  
  Volume <int> newProcVol;
  procVol.halve(axis, newProcVol);

  // sort cluster section on centroids along appropriate axis
#if COSMO_MCLB > 1 
  CkPrintf("(ScaledORBMapBG::map) Calling sortOnAxis(%d, %d, %d)\n", axis, section.beg, section.end);
#endif 
  sortOnAxis(axis, section);

  // this step decides the processor partition to which each set of clusters is assigned
  int clusters = totalClusters/2;

  ClusterSection cs1 = ClusterSection(section.beg, section.beg+clusters-1);
  ClusterSection cs2 = ClusterSection(section.beg+clusters, section.end);

  CmiAssert(2*clusters == totalClusters);

  map(cs1, procVol, clusters);
  map(cs2, newProcVol, totalClusters-clusters);
  
}

void ScaledORBMapBG::sortOnAxis(int axis, ClusterSection section){
  q_sort(axis, section.beg, section.end);
#if COSMO_MCLB > 1 
  CkPrintf("(%d:%d) - \n", section.beg, section.end);
  for(int i = section.beg; i <= section.end; i++)
    CkPrintf("%f, ", clusterArray[i].centroid[axis]);
  CkPrintf("\n");
#endif 
}

void ScaledORBMapBG::q_sort(int axis, int left, int right)
{
  Cluster pivot;
  int pivotPosition;
  int l_hold, r_hold;

  l_hold = left;
  r_hold = right;
  pivot = clusterArray[left];
  while (left < right)
  {
    while ((clusterArray[right].centroid[axis] >= pivot.centroid[axis]) && (left < right))
      right--;
    if (left != right)
    {
      clusterArray[left] = clusterArray[right];
      left++;
    }
    while ((clusterArray[left].centroid[axis] <= pivot.centroid[axis]) && (left < right))
      left++;
    if (left != right)
    {
      clusterArray[right] = clusterArray[left];
      right--;
    }
  }
  clusterArray[left] = pivot;
  pivotPosition = left;
  left = l_hold;
  right = r_hold;
  if (left < pivotPosition)
    q_sort(axis, left, pivotPosition-1);
  if (right > pivotPosition)
    q_sort(axis, pivotPosition+1, right);
}

#ifdef CMK_VERSION_BLUEGENE
#define push(a,b,c) q.push_back(manager->coordinatesToRank(a,b,c))
#else
#define push(a,b,c) q.push_back(a + xsize*b+ xsize*ysize*c)
#endif

#define prev(a,x) (a ==0? x##size-1:a -1)

void ScaledORBMapBG::enqueueNeighbors(int x, int y, int z, deque<int> &q, int dist){
	
        // In Coprocessor mode, each proc. has 6 neighbors, with 
	// node coordinates matching processor coordinates
#ifdef CMK_VERSION_BLUEGENE
	//if(!manager->isVnodeMode()){
          if((manager->getProcsPerNode() == 1)){

#else
		if(!isVnodeMode){
#endif
		push(prev(x,x), y, z);
		push(x, prev(y,y), z);
		push(x, y, prev(z,z));
		push((x+1)%xsize, y, z);
		push(x, (y+1)%ysize, z);
		push(x, y, (z+1)%zsize);
	}
	// VN mode - each processor has 13 neighbors,
	// enqueue them all. Each node holds two processors,
	// and processor coordinates needn't correspond to node
	// coordinates
	else{
		push(x, (y+1)%ysize, z);
		push(x, (y+2)%ysize, z);
		push(x, prev(y,y), z);
		push(x, prev(prev(y,y),y), z);
		
		push(x, y, (z+1)%zsize);
		push(x, y, (z+2)%zsize);
		push(x, y, prev(z,z));
		push(x, y, prev(prev(z,z),z));
		
		push((x+1)%xsize, y, z);
		push((x+2)%xsize, y, z);
		push(prev(x,x), y, z);
		push(prev(prev(x, x),x), y, z);
		
#ifdef CMK_VERSION_BLUEGENE
		if(manager->coordinatesToRank(x,y,z)%2 == 0){
#else
		if((x + y * xsize + z * xsize * ysize)%2 == 0){
#endif
			push((x+3)%xsize, y, z);
		}
		else{
			push(prev(prev(prev(x,x),x),x), y, z);
		}
	}
}


bool ScaledORBMapBG::isNeighborOf(int p1, int p2){
  deque<int> neighbors;
  int x,y,z;

  int i;

#ifdef CMK_VERSION_BLUEGENE
  manager->rankToCoordinates(p1, x, y, z);
#else
  x = p1 % xsize;
  y = (p1 % (xsize * ysize)) / xsize;
  z = p1 / (xsize * ysize);    
#endif
  
  enqueueNeighbors(x, y, z, neighbors, 1);
  for(i = 0; i < neighbors.size(); i++)
    if(neighbors[i] == p2)
      return true;
  
  return false;

}

