/** \file MultipoleMoments.h
 This file defines the representation of a multipole expansion and 
 operations between expansions.
 @author Graeme Lufkin (gwl@u.washington.edu)
 @version 1.0
 */

#ifndef MULTIPOLEMOMENTS_H
#define MULTIPOLEMOMENTS_H

#include <cmath>

#include "Vector3D.h"

/// A representation of a multipole expansion.
class MultipoleMoments {
public:
	/// A physical size for this multipole expansion, calculated by an external function using some other information
	double radius;

	/// The total mass represented by this expansion
	double totalMass;
	/// The center of mass (zeroth order multipole), not divided by the total mass
	Vector3D<double> cm;
	//Tensor for higher order moments goes here (not yet implemented)
	
	MultipoleMoments() : radius(0), totalMass(0) { }
	
	/// Add two expansions together, using parallel axis theorem
	MultipoleMoments& operator+=(const MultipoleMoments& m) {
		//radius gets set by external function
		totalMass += m.totalMass;
		cm += m.cm;
		//add higher order components here
		return *this;
	}
	
	/// Add the contribution of a particle to this multipole expansion
	template <typename ParticleType>
	MultipoleMoments& operator+=(const ParticleType& p) {
		totalMass += p.mass;
		cm += p.mass * p.position;
		//add higher order components here
		return *this;
	}
	
	/// Subtract an expansion from this larger one, yielding the leftover
	MultipoleMoments operator-(const MultipoleMoments& m) {
		MultipoleMoments newMoments;
		newMoments.totalMass = totalMass - m.totalMass;
		newMoments.cm = cm - m.cm;
		//subtract off higher order components here
		return newMoments;
	}
	
	/// Reset this expansion to nothing
	void clear() {
		radius = 0;
		totalMass = 0;
		cm.x = cm.y = cm.z = 0;
		//clear higher order components here
	}
	
};

#ifdef __CHARMC__

#include "pup.h"

inline void operator|(PUP::er& p, MultipoleMoments& m) {
	p | m.radius;
	p | m.totalMass;
	p | m.cm;
}

#endif //__CHARMC__

//What follows are criteria for deciding the size of a multipole

/// Given an enclosing box, set the multipole expansion size to the distance from the center of mass to the farthest corner of the box
inline void calculateRadiusFarthestCorner(MultipoleMoments& m, const OrientedBox<double>& box) {
	Vector3D<double> delta1 = m.cm / m.totalMass - box.lesser_corner;	
	Vector3D<double> delta2 = box.greater_corner - m.cm / m.totalMass;
	delta1.x = (delta1.x > delta2.x ? delta1.x : delta2.x);
	delta1.y = (delta1.y > delta2.y ? delta1.y : delta2.y);
	delta1.z = (delta1.z > delta2.z ? delta1.z : delta2.z);
	m.radius = delta1.length();
}

/// Given the positions that make up a multipole expansion, set the distance to the farthest particle from the center of mass
template <typename ParticleType>
inline void calculateRadiusFarthestParticle(MultipoleMoments& m, const ParticleType* begin, const ParticleType* end) {
	Vector3D<double> cm = m.cm / m.totalMass;
	double d;
	m.radius = 0;
	for(ParticleType* iter = begin; iter != end; ++iter) {
		d = (cm - iter->position).lengthSquared();
		if(d > m.radius)
			m.radius = d;
	}
	m.radius = sqrt(m.radius);
}

#endif //MULTIPOLEMOMENTS_H
