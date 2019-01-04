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
#include <pup.h>

#include <OrientedBox.h>
#include <Vector3D.h>
#ifdef HEXADECAPOLE
#include "moments.h"
#endif

#include "SSEdefs.h"

#if CMK_SSE && defined(HEXADECAPOLE) 
/*
 ** This is a new fast version of QEVAL which evaluates
 ** the interaction due to the reduced moment 'm'.
 ** This version is nearly two times as fast as a naive
 ** implementation.
 **
 ** OpCount = (*,+) = (103,72) = 175 - 8 = 167
 */
inline 
void momEvalMomr(MOMR *m,SSEcosmoType dir0,SSEcosmoType x,SSEcosmoType y,
		 SSEcosmoType z,SSEcosmoType *fPot,SSEcosmoType *ax,
		 SSEcosmoType *ay,SSEcosmoType *az)
{
	const SSEcosmoType onethird = 1.0/3.0;
	SSEcosmoType xx,xy,xz,yy,yz,zz;
	SSEcosmoType xxx,xxy,xxz,xyy,yyy,yyz,xyz;
	SSEcosmoType tx,ty,tz,dir2,g2,g3,g4;
	SSEcosmoType dir;

	dir = -dir0;
	dir2 = dir*dir;
	g2 = 3.0*dir*dir2*dir2;
	g3 = -5.0*g2*dir2;
	g4 = -7.0*g3*dir2;
	/*
	 ** Calculate the funky distance terms.
	 */
	xx = 0.5*x*x;
	xy = x*y;
	xz = x*z;
	yy = 0.5*y*y;
	yz = y*z;
	zz = 0.5*z*z;
	xxx = x*(onethird*xx - zz);
	xxz = z*(xx - onethird*zz);
	yyy = y*(onethird*yy - zz);
	yyz = z*(yy - onethird*zz);
	xx -= zz;
	yy -= zz;
	xxy = y*xx;
	xyy = x*yy;
	xyz = xy*z;
	/*
	 ** Now calculate the interaction up to Hexadecapole order.
	 */
	tx = g4*(m->xxxx*xxx + m->xyyy*yyy + m->xxxy*xxy + m->xxxz*xxz + m->xxyy*xyy + m->xxyz*xyz + m->xyyz*yyz);
	ty = g4*(m->xyyy*xyy + m->xxxy*xxx + m->yyyy*yyy + m->yyyz*yyz + m->xxyy*xxy + m->xxyz*xxz + m->xyyz*xyz);
	tz = g4*(-m->xxxx*xxz - (m->xyyy + m->xxxy)*xyz - m->yyyy*yyz + m->xxxz*xxx + m->yyyz*yyy - m->xxyy*(xxz + yyz) + m->xxyz*xxy + m->xyyz*xyy);
	g4 = 0.25*(tx*x + ty*y + tz*z);
	xxx = g3*(m->xxx*xx + m->xyy*yy + m->xxy*xy + m->xxz*xz + m->xyz*yz);
	xxy = g3*(m->xyy*xy + m->xxy*xx + m->yyy*yy + m->yyz*yz + m->xyz*xz);
	xxz = g3*(-(m->xxx + m->xyy)*xz - (m->xxy + m->yyy)*yz + m->xxz*xx + m->yyz*yy + m->xyz*xy);
	g3 = onethird*(xxx*x + xxy*y + xxz*z);
	xx = g2*(m->xx*x + m->xy*y + m->xz*z);
	xy = g2*(m->yy*y + m->xy*x + m->yz*z);
	xz = g2*(-(m->xx + m->yy)*z + m->xz*x + m->yz*y);
	g2 = 0.5*(xx*x + xy*y + xz*z);
	dir *= m->m;
	dir2 *= -(dir + 5.0*g2 + 7.0*g3 + 9.0*g4);
	*fPot += dir + g2 + g3 + g4;
	*ax += xx + xxx + tx + x*dir2;
	*ay += xy + xxy + ty + y*dir2;
	*az += xz + xxz + tz + z*dir2;
	}
/*
 ** This is a new fast version of QEVAL which evaluates
 ** the interaction due to the reduced moment 'm'.
 ** This version is nearly two times as fast as a naive
 ** implementation.
 **
 ** March 23, 2007: This function now uses unit vectors
 ** which reduces the required precision in the exponent
 ** since the highest power of r is now 5 (g4 ~ r^(-5)).
 **
 ** OpCount = (*,+) = (106,72) = 178 - 8 = 170
 **
 */
inline
void momEvalFmomrcm(FMOMR *m,SSEcosmoType u,SSEcosmoType dir,SSEcosmoType x,
		    SSEcosmoType y,SSEcosmoType z,
		    SSEcosmoType *fPot,SSEcosmoType *ax,SSEcosmoType *ay,
		    SSEcosmoType *az,SSEcosmoType *magai) {
    const SSEcosmoType onethird = 1.0f/3.0f;
    SSEcosmoType xx,xy,xz,yy,yz,zz;
    SSEcosmoType xxx,xxy,xxz,xyy,yyy,yyz,xyz;
    SSEcosmoType tx,ty,tz,g0,g2,g3,g4;

    u *= dir;
    g0 = dir;
    g2 = 3*dir*u*u;
    g3 = 5*g2*u;
    g4 = 7*g3*u;
    /*
     ** Calculate the trace-free distance terms.
     */
    x *= dir;
    y *= dir;
    z *= dir;
    xx = 0.5*x*x;
    xy = x*y;
    xz = x*z;
    yy = 0.5*y*y;
    yz = y*z;
    zz = 0.5*z*z;
    xxx = x*(onethird*xx - zz);
    xxz = z*(xx - onethird*zz);
    yyy = y*(onethird*yy - zz);
    yyz = z*(yy - onethird*zz);
    xx -= zz;
    yy -= zz;
    xxy = y*xx;
    xyy = x*yy;
    xyz = xy*z;
    /*
     ** Now calculate the interaction up to Hexadecapole order.
     */
    tx = g4*(m->xxxx*xxx + m->xyyy*yyy + m->xxxy*xxy + m->xxxz*xxz + m->xxyy*xyy + m->xxyz*xyz + m->xyyz*yyz);
    ty = g4*(m->xyyy*xyy + m->xxxy*xxx + m->yyyy*yyy + m->yyyz*yyz + m->xxyy*xxy + m->xxyz*xxz + m->xyyz*xyz);
    tz = g4*(-m->xxxx*xxz - (m->xyyy + m->xxxy)*xyz - m->yyyy*yyz + m->xxxz*xxx + m->yyyz*yyy - m->xxyy*(xxz + yyz) + m->xxyz*xxy + m->xyyz*xyy);
    g4 = 0.25*(tx*x + ty*y + tz*z);
    xxx = g3*(m->xxx*xx + m->xyy*yy + m->xxy*xy + m->xxz*xz + m->xyz*yz);
    xxy = g3*(m->xyy*xy + m->xxy*xx + m->yyy*yy + m->yyz*yz + m->xyz*xz);
    xxz = g3*(-(m->xxx + m->xyy)*xz - (m->xxy + m->yyy)*yz + m->xxz*xx + m->yyz*yy + m->xyz*xy);
    g3 = onethird*(xxx*x + xxy*y + xxz*z);
    xx = g2*(m->xx*x + m->xy*y + m->xz*z);
    xy = g2*(m->yy*y + m->xy*x + m->yz*z);
    xz = g2*(-(m->xx + m->yy)*z + m->xz*x + m->yz*y);
    g2 = 0.5*(xx*x + xy*y + xz*z);
    g0 *= m->m;
    *fPot += -(g0 + g2 + g3 + g4);
    g0 += 5*g2 + 7*g3 + 9*g4;
    *ax += dir*(xx + xxx + tx - x*g0);
    *ay += dir*(xy + xxy + ty - y*g0);
    *az += dir*(xz + xxz + tz - z*g0);
    *magai = g0*dir;
    }
#endif

/// A representation of a multipole expansion.
class MultipoleMoments {
	friend class CudaMultipoleMoments;
	/// A physical size for this multipole expansion, calculated
	/// by an external function using some other information
	cosmoType radius;
public:
	cosmoType soft;		/* Effective softening */

	/// The total mass represented by this expansion
	cosmoType totalMass;
	/// The center of mass (zeroth order multipole)
	Vector3D<cosmoType> cm;
#ifdef HEXADECAPOLE
	FMOMR mom;
#else						\
	//Tensor for higher order moments goes here
	double xx, xy, xz, yy, yz, zz;
#endif
	
	MultipoleMoments() : radius(0), totalMass(0) { 
	    soft = 0;
		cm.x = cm.y = cm.z = 0;
		//clear higher order components here
#ifdef HEXADECAPOLE
		momClearFmomr(&mom);
#else
		xx = xy = xz = yy = yz = zz = 0;
#endif
	    }
	
	/// Add two expansions together, using parallel axis theorem
	MultipoleMoments& operator+=(const MultipoleMoments& m) {
		//radius gets set by external function
	    cosmoType m1 = totalMass;
		totalMass += m.totalMass;
		if(totalMass == 0.0) {
		    soft = 0.5*(soft + m.soft);
		    cm = 0.5*(cm + m.cm);
		    return *this;
		    }
		if(m1 == 0.0) {	/* just copy over argument */
		    *this = m;
		    return *this;
		    }
		if(m.totalMass == 0.0)
		    return *this;
		soft = (m1*soft + m.totalMass*m.soft)/totalMass;
		Vector3D<cosmoType> cm1 = cm;
		cm = (m1*cm + m.totalMass*m.cm)/totalMass;
#ifdef HEXADECAPOLE
		Vector3D<cosmoType> dr = cm1 - cm;
		momShiftFmomr(&mom, radius, dr.x, dr.y, dr.z);
		FMOMR mom2 = m.mom;
		dr = m.cm - cm;
		momShiftFmomr(&mom2, m.radius, dr.x, dr.y, dr.z);
		momScaledAddFmomr(&mom, radius, &mom2, m.radius);
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
	    cosmoType m1 = totalMass;
		totalMass += p.mass;
               if(totalMass == 0.0) {
                   soft = 0.5*(soft + p.soft);
                   cm = 0.5*(cm + p.position);
                   return *this;
                   }
		soft = (m1*soft + p.mass*p.soft)/totalMass;
		Vector3D<cosmoType> cm1 = cm;
		cm = (m1*cm + p.mass * p.position)/totalMass;
#ifdef HEXADECAPOLE
		// XXX this isn't the most efficient way, but it
		// retains the semantics of this function.  It would
		// be better to do this many particles at a time, then
		// you could first determine the center of mass, then
		// do a momMakeMomr(); momAddMomr() for each particle.
		Vector3D<cosmoType> dr = cm1 - cm;
		momShiftFmomr(&mom, radius, dr.x, dr.y, dr.z);
		dr = p.position - cm;
		FMOMR momPart;
		momMakeFmomr(&momPart, p.mass, radius, dr.x, dr.y, dr.z);
		momAddFmomr(&mom, &momPart);
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
               if(newMoments.totalMass == 0.0) {
                   soft = 0.5*(soft - m.soft);
                   cm = 0.5*(cm - m.cm);
                   return *this;
                   }
		newMoments.soft = (totalMass*soft - m.totalMass*m.soft)
		    /newMoments.totalMass;
		newMoments.cm = (totalMass*cm - m.totalMass*m.cm)
		    /newMoments.totalMass;
#ifdef HEXADECAPOLE
		Vector3D<cosmoType> dr = cm - newMoments.cm;
		newMoments.mom = mom;
		momShiftFmomr(&mom, radius, dr.x, dr.y, dr.z);
		FMOMR mom2 = m.mom;
		dr = m.cm - newMoments.cm;
		momShiftFmomr(&mom2, m.radius, dr.x, dr.y, dr.z);
		momScaledSubFmomr(&newMoments.mom, radius, &mom2, m.radius);
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
		radius = 0;
		totalMass = 0;
		cm.x = cm.y = cm.z = 0;
		//clear higher order components here
#ifdef HEXADECAPOLE
		momClearFmomr(&mom);
#else
		xx = xy = xz = yy = yz = zz = 0;
#endif
	}
	inline cosmoType getRadius() const {return radius;}
	inline cosmoType& getRadius() {return radius;}
	friend void operator|(PUP::er& p, MultipoleMoments& m);
	friend void calculateRadiusFarthestCorner(MultipoleMoments& m,
					      const OrientedBox<double>& box);
	template<typename ParticleType>
	    friend void calculateRadiusFarthestParticle(MultipoleMoments& m,
						const ParticleType * begin,
						const ParticleType * end);
	friend void calculateRadiusBox(MultipoleMoments& m,
				       const OrientedBox<double>& box);
};

#ifdef __CHARMC__

#include "pup.h"

inline void operator|(PUP::er& p, MultipoleMoments& m) {
	p | m.radius;
	p | m.totalMass;
	p | m.soft;
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
	Vector3D<cosmoType> delta1 = m.cm - box.lesser_corner;	
	Vector3D<cosmoType> delta2 = box.greater_corner - m.cm;
	delta1.x = (delta1.x > delta2.x ? delta1.x : delta2.x);
	delta1.y = (delta1.y > delta2.y ? delta1.y : delta2.y);
	delta1.z = (delta1.z > delta2.z ? delta1.z : delta2.z);
	cosmoType newradius = delta1.length();
#ifdef HEXADECAPOLE
	momRescaleFmomr(&m.mom, newradius, m.radius);
#endif
	m.radius = newradius;
}

/// Given an enclosing box, set the multipole expansion size to the
/// distance from the center of the box to the farthest corner of the box
inline void calculateRadiusBox(MultipoleMoments& m,
			       const OrientedBox<double>& box) {
	Vector3D<cosmoType> delta = box.greater_corner - box.lesser_corner;
	cosmoType newradius = 0.5*delta.length();
#ifdef HEXADECAPOLE
	if(m.totalMass > 0.0)
	    momRescaleFmomr(&m.mom, newradius, m.radius);
#endif
	m.radius = newradius;
}

/// Given the positions that make up a multipole expansion, set the distance to the farthest particle from the center of mass
template <typename ParticleType>
inline void calculateRadiusFarthestParticle(MultipoleMoments& m, const ParticleType* begin, const ParticleType* end) {
	Vector3D<cosmoType> cm = m.cm;
	cosmoType d;
	cosmoType newradius = 0;
	for(const ParticleType* iter = begin; iter != end; ++iter) {
		d = (cm - iter->position).lengthSquared();
		if(d > newradius)
			newradius = d;
            }
        if(newradius > 0.0) {
            newradius = sqrt(newradius);
#ifdef HEXADECAPOLE
            momRescaleFmomr(&m.mom, newradius, m.radius);
#endif
            m.radius = newradius;
            }
}

#endif //MULTIPOLEMOMENTS_H
