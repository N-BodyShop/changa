/** \file MultipoleMoments.h
 This file defines the representation of a multipole expansion and 
 operations between expansions.
 @author Graeme Lufkin (gwl@u.washington.edu)
 @version 1.0
 */

#ifndef MULTIPOLEMOMENTS_H
#define MULTIPOLEMOMENTS_H

#include <cmath>
#include <assert.h>

#include "Vector3D.h"

/// A representation of a multipole expansion.
class MultipoleMoments {
public:
	/// A physical size for this multipole expansion, calculated by an external function using some other information
	double radius;
	double soft;

	/// The total mass represented by this expansion
	double totalMass;
	/// The center of mass (zeroth order multipole)
	Vector3D<double> cm;
	//Tensor for higher order moments goes here
	double xx, xy, xz, yy, yz, zz;
	
	MultipoleMoments() : radius(0), totalMass(0) { 
	    soft = 0;
		cm.x = cm.y = cm.z = 0;
		//clear higher order components here
		xx = xy = xz = yy = yz = zz = 0;
	    }
	
	/// Add two expansions together, using parallel axis theorem
	MultipoleMoments& operator+=(const MultipoleMoments& m) {
		//radius gets set by external function
	    double m1 = totalMass;
		totalMass += m.totalMass;
		if(totalMass == 0.0) {
		    soft = 0.5*(soft + m.soft);
		    cm = 0.5*(cm + m.cm);
		    return *this;
		    }
		soft = (m1*soft + m.totalMass*m.soft)/totalMass;
		Vector3D<double> cm1 = cm;
		cm = (m1*cm + m.totalMass*m.cm)/totalMass;
		//add higher order components here
		Vector3D<double> dr = cm1 - cm;
		xx += m1*dr[0]*dr[0];
		yy += m1*dr[1]*dr[1];
		zz += m1*dr[2]*dr[2];
		xy += m1*dr[0]*dr[1];
		xz += m1*dr[0]*dr[2];
		yz += m1*dr[1]*dr[2];
		dr = m.cm - cm;
		xx += m.xx + m.totalMass*dr[0]*dr[0];
		yy += m.yy + m.totalMass*dr[1]*dr[1];
		zz += m.zz + m.totalMass*dr[2]*dr[2];
		xy += m.xy + m.totalMass*dr[0]*dr[1];
		xz += m.xz + m.totalMass*dr[0]*dr[2];
		yz += m.yz + m.totalMass*dr[1]*dr[2];
		return *this;
	}
	
	/// Add the contribution of a particle to this multipole expansion
	template <typename ParticleType>
	MultipoleMoments& operator+=(const ParticleType& p) {
	    double m1 = totalMass;
		totalMass += p.mass;
		Vector3D<double> cm1 = cm;
		cm = (m1*cm + p.mass * p.position)/totalMass;
		//add higher order components here
		Vector3D<double> dr = cm1 - cm;

		xx += m1*dr[0]*dr[0];
		yy += m1*dr[1]*dr[1];
		zz += m1*dr[2]*dr[2];
		xy += m1*dr[0]*dr[1];
		xz += m1*dr[0]*dr[2];
		yz += m1*dr[1]*dr[2];
		
		dr = p.position - cm;
		
		xx += p.mass*dr[0]*dr[0];
		yy += p.mass*dr[1]*dr[1];
		zz += p.mass*dr[2]*dr[2];
		xy += p.mass*dr[0]*dr[1];
		xz += p.mass*dr[0]*dr[2];
		yz += p.mass*dr[1]*dr[2];
		
		return *this;
	}
	
	/// Subtract an expansion from this larger one, yielding the leftover
	MultipoleMoments operator-(const MultipoleMoments& m) {
		MultipoleMoments newMoments;
		newMoments.totalMass = totalMass - m.totalMass;
		newMoments.soft = (totalMass*soft - m.totalMass*m.soft)
		    /newMoments.totalMass;
		newMoments.cm = (totalMass*cm - m.totalMass*m.cm)
		    /newMoments.totalMass;
		//subtract off higher order components here
		Vector3D<double> dr = cm - newMoments.cm;
		newMoments.xx = xx + totalMass*dr[0]*dr[0];
		newMoments.yy = yy + totalMass*dr[1]*dr[1];
		newMoments.zz = zz + totalMass*dr[2]*dr[2];
		newMoments.xy = xy + totalMass*dr[0]*dr[1];
		newMoments.xz = xz + totalMass*dr[0]*dr[2];
		newMoments.yz = yz + totalMass*dr[1]*dr[2];
		dr = m.cm - newMoments.cm;
		newMoments.xx -= m.xx + m.totalMass*dr[0]*dr[0];
		newMoments.yy -= m.yy + m.totalMass*dr[1]*dr[1];
		newMoments.zz -= m.zz + m.totalMass*dr[2]*dr[2];
		newMoments.xy -= m.xy + m.totalMass*dr[0]*dr[1];
		newMoments.xz -= m.xz + m.totalMass*dr[0]*dr[2];
		newMoments.yz -= m.yz + m.totalMass*dr[1]*dr[2];
		return newMoments;
	}
	
	/// Reset this expansion to nothing
	void clear() {
	    soft = 0;
		radius = 0;
		totalMass = 0;
		cm.x = cm.y = cm.z = 0;
		//clear higher order components here
		xx = xy = xz = yy = yz = zz = 0;
	}
	
};

#ifdef __CHARMC__

#include "pup.h"

inline void operator|(PUP::er& p, MultipoleMoments& m) {
	p | m.radius;
	p | m.totalMass;
	p | m.soft;
	p | m.cm;
	p | m.xx;
	p | m.xy;
	p | m.xz;
	p | m.yy;
	p | m.yz;
	p | m.zz;
}

#endif //__CHARMC__

//What follows are criteria for deciding the size of a multipole

/// Given an enclosing box, set the multipole expansion size to the distance from the center of mass to the farthest corner of the box
inline void calculateRadiusFarthestCorner(MultipoleMoments& m, const OrientedBox<double>& box) {
	Vector3D<double> delta1 = m.cm - box.lesser_corner;	
	Vector3D<double> delta2 = box.greater_corner - m.cm;
	delta1.x = (delta1.x > delta2.x ? delta1.x : delta2.x);
	delta1.y = (delta1.y > delta2.y ? delta1.y : delta2.y);
	delta1.z = (delta1.z > delta2.z ? delta1.z : delta2.z);
	m.radius = delta1.length();
}

/// Given the positions that make up a multipole expansion, set the distance to the farthest particle from the center of mass
template <typename ParticleType>
inline void calculateRadiusFarthestParticle(MultipoleMoments& m, ParticleType* begin, const ParticleType* end) {
	Vector3D<double> cm = m.cm;
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
