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
#ifdef HEXADECAPOLE
#include "moments.h"
#endif

/// A representation of a multipole expansion.
class MultipoleMoments {
public:
	/// A physical size for this multipole expansion, calculated by an external function using some other information
	double radius;
	double soft;		/* Effective softening */
	double softBound;	/* Radius to which softening extends */

	/// The total mass represented by this expansion
	double totalMass;
	/// The center of mass (zeroth order multipole)
	Vector3D<double> cm;
#ifdef HEXADECAPOLE
	MOMR mom;
#else						\
	//Tensor for higher order moments goes here
	double xx, xy, xz, yy, yz, zz;
#endif
	
	MultipoleMoments() : radius(0), totalMass(0) { 
	    soft = 0;
	    softBound = 0;
		cm.x = cm.y = cm.z = 0;
		//clear higher order components here
#ifdef HEXADECAPOLE
		momClearMomr(&mom);
#else
		xx = xy = xz = yy = yz = zz = 0;
#endif
	    }
	
	/// Add two expansions together, using parallel axis theorem
	MultipoleMoments& operator+=(const MultipoleMoments& m) {
		//radius gets set by external function
	    double m1 = totalMass;
		totalMass += m.totalMass;
		if(totalMass == 0.0) {
		    soft = 0.5*(soft + m.soft);
		    cm = 0.5*(cm + m.cm);
#ifdef HEXADECAPOLE
		    Vector3D<double> dr = cm - m.cm;
		    if(softBound >= m.softBound)
			softBound = softBound + 0.5*dr.length();
		    else
			softBound = m.softBound + 0.5*dr.length();
#endif
		    return *this;
		    }
		soft = (m1*soft + m.totalMass*m.soft)/totalMass;
		Vector3D<double> cm1 = cm;
		cm = (m1*cm + m.totalMass*m.cm)/totalMass;
#ifdef HEXADECAPOLE
		Vector3D<double> dr = cm1 - cm;
		double sMax1 = dr.length() + softBound;
		momShiftMomr(&mom, dr.x, dr.y, dr.z);
		MOMR mom2 = m.mom;
		dr = m.cm - cm;
		double sMax2 = dr.length() + m.softBound;
		if(sMax1 > sMax2)
		    softBound = sMax1;
		else
		    softBound = sMax2;
		momShiftMomr(&mom2, dr.x, dr.y, dr.z);
		momAddMomr(&mom, &mom2);
#else
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
#endif
		return *this;
	}
	
	/// Add the contribution of a particle to this multipole expansion
	template <typename ParticleType>
	MultipoleMoments& operator+=(const ParticleType& p) {
	    double m1 = totalMass;
		totalMass += p.mass;
		soft = (m1*soft + p.mass*p.soft)/totalMass;
		Vector3D<double> cm1 = cm;
		cm = (m1*cm + p.mass * p.position)/totalMass;
#ifdef HEXADECAPOLE
		// XXX this isn't the most efficient way, but it
		// retains the semantics of this function.  It would
		// be better to do this many particles at a time, then
		// you could first determine the center of mass, then
		// do a momMakeMomr(); momAddMomr() for each particle.
		Vector3D<double> dr = cm1 - cm;
		double sMax1 = dr.length() + softBound;
		if(m1 == 0.0)	// first particle is being added.
		    sMax1 = 0.0;
		momShiftMomr(&mom, dr.x, dr.y, dr.z);
		dr = p.position - cm;
		// spline softening: extends to 2*soft
		double sMax2 = dr.length() + 2.0*p.soft;
		if(sMax1 > sMax2)
		    softBound = sMax1;
		else
		    softBound = sMax2;
		MOMR momPart;
		momMakeMomr(&momPart, p.mass, dr.x, dr.y, dr.z);
		momAddMomr(&mom, &momPart);
#else
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
#endif
		
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
#ifdef HEXADECAPOLE
		Vector3D<double> dr = cm - newMoments.cm;
		newMoments.mom = mom;
		momShiftMomr(&mom, dr.x, dr.y, dr.z);
		MOMR mom2 = m.mom;
		dr = m.cm - newMoments.cm;
		momShiftMomr(&mom2, dr.x, dr.y, dr.z);
		momSubMomr(&newMoments.mom, &mom2);
		// I don't think we can safely reduce the softening bound
		newMoments.softBound = softBound;
#else
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
#endif
		return newMoments;
	}
	
	/// Reset this expansion to nothing
	void clear() {
	    soft = 0;
	    softBound = 0;
		radius = 0;
		totalMass = 0;
		cm.x = cm.y = cm.z = 0;
		//clear higher order components here
#ifdef HEXADECAPOLE
		momClearMomr(&mom);
#else
		xx = xy = xz = yy = yz = zz = 0;
#endif
	}
	
};

#ifdef __CHARMC__

#include "pup.h"

inline void operator|(PUP::er& p, MultipoleMoments& m) {
	p | m.radius;
	p | m.totalMass;
	p | m.soft;
	p | m.softBound;
	p | m.cm;
#ifdef HEXADECAPOLE
	p((char *) &m.mom, sizeof(m.mom)); /* PUPs as bytes */
#else
	p | m.xx;
	p | m.xy;
	p | m.xz;
	p | m.yy;
	p | m.yz;
	p | m.zz;
#endif
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
