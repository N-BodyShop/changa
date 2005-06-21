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
#include "SFC.h"

using namespace Tipsy;
using namespace SFC;

typedef unsigned char byte;

/** A Key that signals that a tree node is actually a bucket of particles. */
const Key placeholderMask = ((Key) 1 << 63);

/** The base class for particles. 
 This type of particle contains only a key (made from a position).
 */
class BareParticle {
protected:

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
	void initializeExtraFields() {
		key = 0;
		smoothingRadius = 0;
		density = 0;
		velDispSq = 0;
		neighborCount = 0;
		colorValue = 0;
		diskOrder = 0;
		groupParentChare = -1;
		groupParentParticle = 0;
	}
	
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
		initializeExtraFields();
	}

	FullParticle(Key k) : BareParticle(k) { }
	
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
		initializeExtraFields();
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
		initializeExtraFields();
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
		initializeExtraFields();
	}
		
	/// Output operator, used for formatted display
	friend std::ostream& operator<< (std::ostream& os, const FullParticle& bp) {
		std::bitset<32> firstHalf(bp.key >> 32);
		std::bitset<32> secondHalf(bp.key);
		os << firstHalf << secondHalf;
		return os;
	}
	void pup(PUP::er &p)
            {
		p | key;
                p | mass;
                p | position;
                p | velocity;
		p | rho;
		p | temperature;
                p | eps;
		p | h;
		p | metals;
		p | phi;
		p | tform;
		p | smoothingRadius;
		p | divv;
		p | curlv;
		p | meanVelocity;
		p | density;
		p | velDispSq;
		p | neighborCount;
		p | colorValue;
		p | diskOrder;
		p | groupParentChare;
		p | groupParentParticle;
                }

};



#endif //PARTICLE_H
