#ifndef _MAP_STRUCTURES_H_
#define _MAP_STRUCTURES_H_
#include <iostream>

#include "Vector3D.h"
#include "TaggedVector3D.h"
#include "OrientedBox.h"
using namespace std;

#define XDIM 0
#define YDIM 1
#define ZDIM 2
#define NDIMS 3

#define LEFT_PARTITION 0
#define RIGHT_PARTITION 1
#define INVALID_PARTITION -1

template <class T>
class Volume{
	
  public:
  T lesser[3];          // lesser corner
  T greater[3];         // greater corner 
  
  void halve(int whichAxis, Volume<T> &secondhalf);
  void split(int whichAxis, T *where, Volume<T> &secondhalf);
  int decideSplitAxis();

  Volume(){}
  Volume(T *l, T *g){
    for (int i= 0; i < 3; i++){
     lesser[i] = l[i];
     greater[i] = g[i];
    }
  }
};

class Cluster{
  public:
  int load;     // total load induced by this cluster
  int proc;     // which processor this cluster will be assigned to - one cluster per processor
  Volume <int>* partition;      // the lowest level process partition that this cluster belongs to
  float centroid[3];    // cluster centroid - calculated like center of mass, with unit mass for each participating TP. FIXME - refine
  int count;    // this cluster comprises 'count' tree pieces

  CkVec<int> tps;   // to check metis efficacy only. remove later on FIXME

  Cluster& operator= (Cluster& in){
    load = in.load;
    proc = in.proc;
    partition = in.partition;
    count = in.count;
    for(int i=0; i < 2; i++)
      centroid[i] = in.centroid[i];
    return *this;
  }

  Cluster(){
    for(int i = 0; i < 3; i++)
      centroid[i] = 0.0;
      count = 0;
      proc = -1;
  }
};

class ClusterSection{   // indicates which section of the clusterArray is to be worked upon
  public:
  int beg, end;
  ClusterSection(int b, int e): beg(b), end(e){}
};

template <class T>
int Volume<T>::decideSplitAxis(){

  int which = 0;
  T max = abs(greater[0] - lesser[0]);
  
  for(int i = 1; i < 3; i++)
    if(abs(greater[i] - lesser[i]) >= max){
      max = abs(greater[i] - lesser[i]);
      which = i;
    }
  return which;
}

template <class T>
void Volume<T>::halve(int axis, Volume<T> &secondhalf){

  secondhalf = *this;
  greater[axis] = (lesser[axis]+greater[axis])/2;
  secondhalf.lesser[axis] = greater[axis]+1;
}

template <class T>
void Volume<T>::split(int axis, T* where, Volume<T> &secondhalf){
  secondhalf = *this;
  greater[axis] = where[axis];
  secondhalf.lesser[axis] = where[axis];
}

#define NDIMS 3
class Centroid3d{
  public:
  float x;
  float y;
  float z;

  float *pointers[NDIMS];

  Centroid3d(){
    pointers[0] = &x;
    pointers[1] = &y;
    pointers[2] = &z;
  }

  float& operator[](int i){
    return *(pointers[i]);
  }

};

// Each tpobject has three events, one for 
// each component of its centroid. This helps
// us avoid sorting after each recursive partitioning
// step.

class TPObject;

struct Event {
  int owner;
  float load;
  float position;

  Event(float pos, float ld, int o) : 
    position(pos),
    load(ld),
    owner(o)
  {
  }

  Event() : 
    owner(-1),
    load(0.0),
    position(0.0)
  {
  }

  bool operator<=(const Event &e){
    return position <= e.position;
  }

  bool operator>=(const Event &e){
    return position >= e.position;
  }
};

struct OrbObject {
  int partition;
  // XXX do we really need this field?
  int lbindex;

  OrbObject() : 
    partition(-1),
    lbindex(-1)
  {
  }

  OrbObject(int tag) : 
    partition(-1),
    lbindex(tag)
  {
  }
};

class TPObject{
  public:

  Vector3D<float> centroid;
  float load;
  //int index;
  int lbindex;
  bool migratable;
  //int nparticles;

  int whichPartition;

  bool operator<(const TPObject &t) const{
    return load < t.load;
  }

  TPObject() : 
    load(0.0),
    lbindex(-1),
    whichPartition(-1)
  {
  }

};

class Processor{
  public:

  float load;
  int t;

  bool operator<(const Processor &p) const{
    return load > p.load;
  }

};

class Node {
  public:
  int x, y, z;
  CkVec<int> procRanks;
#ifdef PRINT_BOUNDING_BOXES
  OrientedBox<float> box;
#endif
};

typedef int (*ComparatorFn) (const void *, const void *);
#endif // _MAP_STRUCTURES_H_
