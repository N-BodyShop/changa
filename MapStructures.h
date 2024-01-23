#ifndef _MAP_STRUCTURES_H_
#define _MAP_STRUCTURES_H_
///
/// @file MapStructures.h
/// Data structures for spacially aware load balancing.
///
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

#define NDIMS 3

// Each tpobject has three events, one for 
// each component of its centroid. This helps
// us avoid sorting after each recursive partitioning
// step.

class TPObject;

///
/// @brief structure for sorting TreePieces in one dimension for load
/// balancing.
///
struct Event {
    /// Index within TreePiece array or stats->objData array of this TreePiece.
  int owner;
    /// load of this TreePiece
  float load;
    /// Centroid of this TreePiece in a particular dimension.
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

  bool operator<(const Event &e) const {
    return position < e.position;
  }

  bool operator<=(const Event &e) const {
    return position <= e.position;
  }

  bool operator>=(const Event &e) const {
    return position >= e.position;
  }

  void pup(PUP::er &p){
    p|owner;
    p|load;
    p|position;
  }
};

///
/// @brief data for Orb3D load balancing.
///
struct OrbObject {
  int partition;
  /// index into LB stats->objData
  int lbindex;
  /// Spacial location of TreePiece
  Vector3D<float> centroid;
  /// Particles in the TreePiece
  int numParticles;

  OrbObject(int tag) : 
    partition(-1),
    lbindex(tag),
    numParticles(0)
   {
   }

  OrbObject() : 
    partition(-1),
    lbindex(-1),
    numParticles(0)
  {
  }

  OrbObject(int tag, int np) : 
    partition(-1),
    lbindex(tag),
    numParticles(np)
  {
  }

  void pup(PUP::er &p){
    p|partition;
    p|lbindex;
    p|centroid;
    p|numParticles;
  }
};

///
/// @brief hold information needed by the load balancer for each TreePiece.
///
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

///
/// @brief Hold per processor information for load balancing.
///
class Processor{
  public:

  float load;
  int t;

  bool operator<(const Processor &p) const{
    return load > p.load;
  }
  Processor(){
    load = 0.0;
    t = -1;
  }

  Processor(int pe){
    load = 0.0;
    t = pe;
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
