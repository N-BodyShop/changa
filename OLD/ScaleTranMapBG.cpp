#include "MapStructures.h"
#include "ScaleTranMapBG.h"

#include "ParallelGravity.h"
#include <deque>

class Object{
  double weight;// object weight
  int idx;      // which object

  inline bool operator<= (const Object &rhs){
    return weight >= rhs.weight;
  }
  inline bool operator>= (const Object &rhs){
    return weight <= rhs.weight;
  }
};

#define prev(a,x) (a ==0? x##size-1:a -1)
#define abs(x) (x < 0 ? (-1)*x : x)
#define wrap(a,x) (a <0 ? (x##size-(((-1)*(a))%x##size))%x##size: (a)%x##size)
void ScaleTranMapBG::assign(int *from, double *objWeights, CkVec<int> &to, int numobjs, int count){
        int i,j;
        double th = 0.15;

        // tpCentroids gives the centroids of each object, objweights, their weights
        
        double procWeights[count];
        int procCounts[count];
/*
        Vector3D<float> procCentroids[count];
        CkVec<Object> procLists[count];
*/
        double idealAvg = 0.0;
        
/*
        for(i = 0; i < count; i++){
          procCentroids[i].x = 0;
          procCentroids[i].y = 0;
          procCentroids[i].z = 0;
          procWeights[i] = 0;
          procCounts[i] = 0;
        }

        for (i = 0; i < numobjs; i++){
          procLists[from[i]].push_back(Object(i,objWeight[i]));
          procCounts[from[i]]++;
          idealAvg += objWeights[i];
          procWeights[from[i]] += objWeights[i];
          procCentroids[from[i]] += (objWeights[i]*tpCentroids[i]);     //center of mass
        }
      
        idealAvg /= count;
        for(i = 0; i < count; i++){
          procCentroids[i] /= procWeights[i];
          procLists[i].quickSort();     // descending order, because of operator definition
        }

        // have processor centers of mass now
        
        deque q<float *>;
        for(i = 0; i < count; i++){
          float myCoords[0];
          manager->coordinatesToRank(i,myCoords[0], myCoords[1], myCoords[2]);
          
          if(procWeight[i] > (1+th)){
          // try to shift weight out to neighbours
          // enqueue neighbours
          // first, 1-away neighbours: 6 of them
          q.push_back(wrap(myCoords[0]+1,x), myCoords[1], myCoords[2]);
          q.push_back(myCoords[0], wrap(myCoords[1]+1,y), myCoords[2]);
          q.push_back(myCoords[0], myCoords[1], wrap(myCoords[2]+1,z));
          q.push_back(wrap(myCoords[0]-1,x), myCoords[1], myCoords[2]);
          q.push_back(myCoords[0], wrap(myCoords[1]-1, y), myCoords[2]);
          q.push_back(myCoords[0], myCoords[1], wrap(myCoords[2]-1, z));
          

          // 2-away neighbours: 4+4+4 = 12
          // xy plane
          q.push_back(wrap(myCoords[0]+1,x), wrap(myCoords[1]+1,y), myCoords[2]); 
          q.push_back(wrap(myCoords[0]+1,x), wrap(myCoords[1]-1,y), myCoords[2]); 
          q.push_back(wrap(myCoords[0]-1,x), wrap(myCoords[1]+1,y), myCoords[2]); 
          q.push_back(wrap(myCoords[0]-1,x), wrap(myCoords[1]-1,y), myCoords[2]); 
          
          // zy plane
          q.push_back(myCoords[0], wrap(myCoords[1]+1,y), wrap(myCoords[2]+1,z)); 
          q.push_back(myCoords[0], wrap(myCoords[1]-1,y), wrap(myCoords[2]+1,z)); 
          q.push_back(myCoords[0], wrap(myCoords[1]+1,y), wrap(myCoords[2]-1,z)); 
          q.push_back(myCoords[0], wrap(myCoords[1]-1,y), wrap(myCoords[2]-1,z)); 
          
          // xz plane
          q.push_back(wrap(myCoords[0]+1,x), myCoords[1], wrap(myCoords[2]+1,z)); 
          q.push_back(wrap(myCoords[0]+1,x), myCoords[1], wrap(myCoords[2]-1,z)); 
          q.push_back(wrap(myCoords[0]-1,x), myCoords[1], wrap(myCoords[2]+1,z)); 
          q.push_back(wrap(myCoords[0]-1,x), myCoords[1], wrap(myCoords[2]-1,z)); 

          // 3-away neighbours: 4+4
          q.push_back(wrap(myCoords[0]+1,x), wrap(myCoords[1]+1,y), wrap(myCoords[2]+1,z)); 
          q.push_back(wrap(myCoords[0]+1,x), wrap(myCoords[1]-1,y), wrap(myCoords[2]+1,z)); 
          q.push_back(wrap(myCoords[0]-1,x), wrap(myCoords[1]+1,y), wrap(myCoords[2]+1,z)); 
          q.push_back(wrap(myCoords[0]-1,x), wrap(myCoords[1]-1,y), wrap(myCoords[2]+1,z)); 
          
          q.push_back(wrap(myCoords[0]+1,x), wrap(myCoords[1]+1,y), wrap(myCoords[2]-1,z)); 
          q.push_back(wrap(myCoords[0]+1,x), wrap(myCoords[1]-1,y), wrap(myCoords[2]-1,z)); 
          q.push_back(wrap(myCoords[0]-1,x), wrap(myCoords[1]+1,y), wrap(myCoords[2]-1,z)); 
          q.push_back(wrap(myCoords[0]-1,x), wrap(myCoords[1]-1,y), wrap(myCoords[2]-1,z)); 

          // check which objects to shift 
          while(!q.empty()){
          }
          }
        }
        */
	
}

void ScaleTranMapBG::roundOff(Vector3D<float> &vec){
//	int largestSmaller;
  int x,y,z;
  x = (int)((vec.x+0.5)/2.0);
  y = (int)((vec.y+0.5)/2.0);
  z = (int)((vec.z+0.5)/2.0);

  vec.x = x;
  vec.y = y;
  vec.z = z;
}

inline void ScaleTranMapBG::translate(Vector3D<float> &vec, float xdist, float ydist, float zdist){
	vec.x += xdist; 
	vec.y += ydist; 
	vec.z += zdist; 
}

inline void ScaleTranMapBG::scale(Vector3D<float> &vec, int xtl, int ytl, int ztl){
	vec.x *= xtl;
	vec.y *= ytl;
	vec.z *= ztl;
}

int ScaleTranMapBG::map(Vector3D<float> &vec, int numprocs, bool *avail, int *ringRadius){
	deque<int> q;
#ifdef CMK_VERSION_BLUEGENE
        int rank = manager->coordinatesToRank((int)vec.x, (int)vec.y, (int)vec.z);
	q.push_back(manager->coordinatesToRank((int)vec.x, (int)vec.y, (int)vec.z));
        //CkPrintf("Insert %d into queue\n", manager->coordinatesToRank((int)vec.x,(int)vec.y,(int)vec.z));
#else
        int rank = (int)vec.x + ((int)vec.y) * xsize + ((int)vec.z) * xsize * ysize;
	q.push_back((int)vec.x + ((int)vec.y) * xsize + ((int)vec.z) * xsize * ysize);
        //CkPrintf("Insert %d into queue\n", (((int)vec.x) + xsize*((int)vec.y)+ xsize*ysize*((int)vec.z)));
#endif
	int numChecked = 0;
        int cur;
	bool found = false;
		
        while(!found){
          //CkPrintf("Enq. processor %d %d-neighbours:\n", manager->coordinatesToRank((int)vec.x, (int)vec.y, (int)vec.z), ringRadius[rank]);
          enqueueNeighbors((int)vec.x, (int)vec.y, (int)vec.z, q, ringRadius[rank]);
          while(!q.empty()){
            cur = q.front();
            q.pop_front();
            //CkPrintf("Deq'd processor %d - ", cur);
            if(avail[cur]){
              found = true;
              //CkPrintf("free\n");
              avail[cur] = false;
              return cur;
            }
            else{
              //CkPrintf("not free\n");
              collisions++;
              CmiAssert(collisions < numprocs*numprocs);        // A tighter bound is numprocs*(numprocs+1)/2
            }
          }
          ringRadius[rank]++;
        }
        /*
	if(q.empty() && !found)
		CmiAbort("ScaleTranMapBG: map failed\n");
        */
	return -1;
}
#ifdef CMK_VERSION_BLUEGENE
#define push(a,b,c) q.push_back(manager->coordinatesToRank(a,b,c));//CkPrintf("Insert %d into queue\n", manager->coordinatesToRank(a,b,c));
#else
#define push(a,b,c) q.push_back(a + xsize*b+ xsize*ysize*c); //CkPrintf("Insert %d into queue\n", (a + xsize*b+ xsize*ysize*c));
#endif


void ScaleTranMapBG::enqueueNeighbors(int x, int y, int z, deque<int> &q, int d){
	
	int dx, dy, dz;		// distance from x,y,z in each dimension
	int xn, yn, zn;		// neigbour coordinates
        int rank;
	// In Coprocessor mode, each proc. has 6 neighbors, with 
	// node coordinates matching processor coordinates
	// currently, support only for co-processor mode
#ifdef CMK_VERSION_BLUEGENE
	//if(!manager->isVnodeMode()){
          if(manager->getProcsPerNode() == 1){
#else
	if(!isVnodeMode){
#endif
                for(dz = (-1)*d; dz <= d; dz++){
                  for(dx=0, dy = d-abs(dz); dx <= d-abs(dz) && dy >= 0;){
                    //CkPrintf("1: %d,%d,%d (%d,%d,%d)\n", wrap(x+dx,x), wrap(y+dy,y), wrap(z+dz,z), x+dx, y+dy, z+dz);
                    /*
                    CmiAssert(wrap(x+dx, x) >= 0);
                    CmiAssert(wrap(y+dy, y) >= 0);
                    CmiAssert(wrap(z+dz, z) >= 0);
                    CmiAssert(wrap(x-dx, x) >= 0);
                    CmiAssert(wrap(y-dy, y) >= 0);
                    */
                    push(wrap(x+dx,x), wrap(y+dy,y), wrap(z+dz,z));
                    if(dz != (-1)*d && dz != d){         // so that one point isn't enqueued twice when dx and dy are both zero
                      //CkPrintf("1: %d,%d,%d (%d,%d,%d)\n", wrap(x-dx,x), wrap(y-dy,y), wrap(z+dz,z), x-dx, y-dy, z+dz);
                      push(wrap(x-dx,x), wrap(y-dy,y), wrap(z+dz,z));
                    }
                    dx++; dy--;
                  }
                  for(dx=1, dy=abs(dz)-d+1; dx < d-abs(dz) && dy < 0;){
                    //CkPrintf("2: %d,%d,%d (%d,%d,%d)\n", wrap(x+dx,x), wrap(y+dy,y), wrap(z+dz,z), x+dx, y+dy, z+dz);
                    /*
                    CmiAssert(wrap(x+dx, x) >= 0);
                    CmiAssert(wrap(y+dy, y) >= 0);
                    CmiAssert(wrap(z+dz, z) >= 0);
                    CmiAssert(wrap(x-dx, x) >= 0);
                    CmiAssert(wrap(y-dy, y) >= 0);
                    */
                    push(wrap(x+dx,x), wrap(y+dy,y), wrap(z+dz,z));
                    if(dz != (-1)*d && dz != d){
                      //CkPrintf("2: %d,%d,%d (%d,%d,%d)\n", wrap(x-dx,x), wrap(y-dy,y), wrap(z+dz,z), x-dx, y-dy, z+dz);
                      push(wrap(x-dx,x), wrap(y-dy,y), wrap(z+dz,z));
                    }
                    dx++; dy++;
                  }
                }
                                                                                        
			
        }
	// VN mode - each processor has 13 neighbors,
	// enqueue them all. Each node holds two processors,
	// and processor coordinates needn't correspond to node
	// coordinates
	else{
		// FIXME - for now, CmiAbort if we get here - too little time! 
		CmiAbort("VN mode not supported by scaling heuristic\n");
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

bool ScaleTranMapBG::isNeighborOf(int p1, int p2){
  deque<int> neighbors;
  int x,y,z;

  int i;

#ifdef CMK_VERSION_BLUEGENE
  manager->coordinatesToRank(p1, x, y, z);
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
