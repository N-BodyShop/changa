/** \file Particle.h
 This file defines the basic particle data structures.
 \author Graeme Lufkin (gwl@u.washington.edu)
 */
 
#ifndef PARTICLE_H
#define PARTICLE_H

#include <iostream>
#include <bitset>

#include "TipsyFile.h"
#include "Vector3D.h"
#include "OrientedBox.h"

using namespace Tipsy;

/** A 64-bit key identifying a position.
 This is used to keep track of particles' locations and tree nodes.
 It is generated from a space-filling curve.
 */
typedef unsigned long long Key;

typedef unsigned char byte;

/** A Key that signals that a tree node is actually a bucket of particles. */
const Key placeholderMask = ((Key) 1 << 63);

/** The very first possible key a particle can take on. */
const Key firstPossibleKey = static_cast<Key>(0);
/** The very last possible key a particle can take on. */
const Key lastPossibleKey = ~(static_cast<Key>(1) << 63);

/** The base class for particles. 
 This type of particle contains only a key (made from a position).
 */
class BareParticle {
protected:
	/** Given the floating point numbers for the location, construct the key. 
	 The key uses 21 of the 23 bits for the floats of the x, y, and z coordinates
	 of the position vector.  This process will only make sense if the position
	 coordinates are in the range [1,2).  The matissa bits are taken, and interleaved
	 in xyz order to form the key.  This makes the key a position on the z-ordering
	 space-filling curve. */
	inline Key makeKey(Vector3D<float> v) const {
		int* ix = reinterpret_cast<int *>(&v.x);
		int* iy = reinterpret_cast<int *>(&v.y);
		int* iz = reinterpret_cast<int *>(&v.z);
		Key key = 0;
		for(int mask = (1 << 22); mask > 2; mask >>= 1) {
			key <<= 3;
			if(*ix & mask)
				key += 4;
			if(*iy & mask)
				key += 2;
			if(*iz & mask)
				key += 1;
		}
		return key;
	}

public:
	Key key;
	
	/** Default constructor, does nothing. */
	BareParticle() { }
	
	/** Make a particle with just a key. */
	BareParticle(Key k) : key(k) { }
	
	/** Given a position, make a key for this particle. */
	BareParticle(const Vector3D<float>& pos) {
		key = makeKey(pos);
	}
	
	/// Copy constructor
	BareParticle(const BareParticle& bp) : key(bp.key) {
	}
	
	~BareParticle() { }
	
	/// Assignment operator
	BareParticle& operator=(const BareParticle& bp) {
		key = bp.key;
		return *this;
	}
	
	/** Comparison operator, used for a sort based on key values. */
	bool operator< (const BareParticle& bp) const {
		return key < bp.key;
	}
	
	/// Output operator, used for formatted display
	friend std::ostream& operator<< (std::ostream& os, const BareParticle& bp) {
		std::bitset<32> firstHalf(bp.key >> 32);
		std::bitset<32> secondHalf(bp.key);
		os << firstHalf << secondHalf;
		return os;
	}
	
};

/** This particle carries fields along with it.
 */
class FullParticle : public BareParticle {
public:
		
	double mass;
	Vector3D<double> position;
	Vector3D<double> velocity;
	double rho;
	double temperature;
	double eps;
	double h;
	double metals;
	double phi;
	double tform;
	
	double smoothingRadius;
	double divv;
	Vector3D<double> curlv;
	Vector3D<double> meanVelocity;
	double density;
	double velDispSq;
	int neighborCount;
	double colorValue;
	int diskOrder;
	
	/// The array index of this particle's parent's chare (-1 means no parent).
	int groupParentChare;
	/// The array index of this particle's parent (on the parent's chare).
	int groupParentParticle;
	
	FullParticle() { 
		key = 0;
		groupParentChare = -1;
		groupParentParticle = 0;
		colorValue = 0;
	}

	FullParticle(const gas_particle& p) {
		mass = p.mass;
		position = p.pos;
		velocity = p.vel;
		rho = p.rho;
		temperature = p.temp;
		eps = p.hsmooth;
		h = 0;
		metals = p.metals;
		phi = p.phi;
		tform = 0;
		groupParentChare = -1;
		groupParentParticle = 0;
		colorValue = 0;
	}

	FullParticle(const dark_particle& p) {
		mass = p.mass;
		position = p.pos;
		velocity = p.vel;
		rho = 0;
		temperature = 0;
		eps = p.eps;
		h = 0;
		metals = 0;
		phi = p.phi;
		tform = 0;
		groupParentChare = -1;
		groupParentParticle = 0;
		colorValue = 0;
	}

	FullParticle(const star_particle& p) {
		mass = p.mass;
		position = p.pos;
		velocity = p.vel;
		rho = 0;
		temperature = 0;
		eps = p.eps;
		h = 0;
		metals = p.metals;
		phi = p.phi;
		tform = p.tform;
		groupParentChare = -1;
		groupParentParticle = 0;
		colorValue = 0;
	}
	
	template<class T>
	inline Key generateKey(const OrientedBox<T>& boundingBox) {
		return key = makeKey((position - boundingBox.lesser_corner) / (boundingBox.greater_corner - boundingBox.lesser_corner) + Vector3D<float>(1, 1, 1));
	}
	
	/// Output operator, used for formatted display
	friend std::ostream& operator<< (std::ostream& os, const FullParticle& bp) {
		std::bitset<32> firstHalf(bp.key >> 32);
		std::bitset<32> secondHalf(bp.key);
		os << firstHalf << secondHalf;
		return os;
	}
};



#endif //PARTICLE_H
