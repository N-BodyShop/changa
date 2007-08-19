#include <stdio.h>
#include <math.h>
#include "moments.h"

void momClearMomr(MOMR *mr)
{
	mr->m = 0.0;
	mr->xx = 0.0;
	mr->yy = 0.0;
	mr->xy = 0.0;
	mr->xz = 0.0;
	mr->yz = 0.0;
	mr->xxx = 0.0;
	mr->xyy = 0.0;
	mr->xxy = 0.0;
	mr->yyy = 0.0;
	mr->xxz = 0.0;
	mr->yyz = 0.0;
	mr->xyz = 0.0;
	mr->xxxx = 0.0;
	mr->xyyy = 0.0;
	mr->xxxy = 0.0;
	mr->yyyy = 0.0;
	mr->xxxz = 0.0;
	mr->yyyz = 0.0;
	mr->xxyy = 0.0;
	mr->xxyz = 0.0;
	mr->xyyz = 0.0;
    }


/*
 ** This function adds the complete moment ma to the complete moment mc
 */
void momAddMomc(MOMC *mc,MOMC *ma)
{
	mc->m += ma->m;
	mc->xx += ma->xx;
	mc->yy += ma->yy;
	mc->xy += ma->xy;
	mc->xz += ma->xz;
	mc->yz += ma->yz;
	mc->xxx += ma->xxx;
	mc->xyy += ma->xyy;
	mc->xxy += ma->xxy;
	mc->yyy += ma->yyy;
	mc->xxz += ma->xxz;
	mc->yyz += ma->yyz;
	mc->xyz += ma->xyz;
	mc->xxxx += ma->xxxx;
	mc->xyyy += ma->xyyy;
	mc->xxxy += ma->xxxy;
	mc->yyyy += ma->yyyy;
	mc->xxxz += ma->xxxz;
	mc->yyyz += ma->yyyz;
	mc->xxyy += ma->xxyy;
	mc->xxyz += ma->xxyz;
	mc->xyyz += ma->xyyz;
	mc->zz += ma->zz;
	mc->xzz += ma->xzz;
	mc->yzz += ma->yzz;
	mc->zzz += ma->zzz;
	mc->xxzz += ma->xxzz;
	mc->xyzz += ma->xyzz;
	mc->xzzz += ma->xzzz;
	mc->yyzz += ma->yyzz;
	mc->yzzz += ma->yzzz;
	mc->zzzz += ma->zzzz;
	}


/*
 ** This function adds the reduced moment ma to the reduced moment mc
 */
void momAddMomr(MOMR *mr,MOMR *ma)
{
	mr->m += ma->m;
	mr->xx += ma->xx;
	mr->yy += ma->yy;
	mr->xy += ma->xy;
	mr->xz += ma->xz;
	mr->yz += ma->yz;
	mr->xxx += ma->xxx;
	mr->xyy += ma->xyy;
	mr->xxy += ma->xxy;
	mr->yyy += ma->yyy;
	mr->xxz += ma->xxz;
	mr->yyz += ma->yyz;
	mr->xyz += ma->xyz;
	mr->xxxx += ma->xxxx;
	mr->xyyy += ma->xyyy;
	mr->xxxy += ma->xxxy;
	mr->yyyy += ma->yyyy;
	mr->xxxz += ma->xxxz;
	mr->yyyz += ma->yyyz;
	mr->xxyy += ma->xxyy;
	mr->xxyz += ma->xxyz;
	mr->xyyz += ma->xyyz;
	}


/*
 ** This function multiply-adds the complete moment ma
 */
void momMulAddMomc(MOMC *mc,momFloat m,MOMC *ma)
{
	mc->m += m*ma->m;
	mc->xx += m*ma->xx;
	mc->yy += m*ma->yy;
	mc->xy += m*ma->xy;
	mc->xz += m*ma->xz;
	mc->yz += m*ma->yz;
	mc->xxx += m*ma->xxx;
	mc->xyy += m*ma->xyy;
	mc->xxy += m*ma->xxy;
	mc->yyy += m*ma->yyy;
	mc->xxz += m*ma->xxz;
	mc->yyz += m*ma->yyz;
	mc->xyz += m*ma->xyz;
	mc->xxxx += m*ma->xxxx;
	mc->xyyy += m*ma->xyyy;
	mc->xxxy += m*ma->xxxy;
	mc->yyyy += m*ma->yyyy;
	mc->xxxz += m*ma->xxxz;
	mc->yyyz += m*ma->yyyz;
	mc->xxyy += m*ma->xxyy;
	mc->xxyz += m*ma->xxyz;
	mc->xyyz += m*ma->xyyz;
	mc->zz += m*ma->zz;
	mc->xzz += m*ma->xzz;
	mc->yzz += m*ma->yzz;
	mc->zzz += m*ma->zzz;
	mc->xxzz += m*ma->xxzz;
	mc->xyzz += m*ma->xyzz;
	mc->xzzz += m*ma->xzzz;
	mc->yyzz += m*ma->yyzz;
	mc->yzzz += m*ma->yzzz;
	mc->zzzz += m*ma->zzzz;
	}


/*
 ** This function multiply-adds the reduced moment ma
 */
void momMulAddMomr(MOMR *mr,momFloat m,MOMR *ma)
{
	mr->m += m*ma->m;
	mr->xx += m*ma->xx;
	mr->yy += m*ma->yy;
	mr->xy += m*ma->xy;
	mr->xz += m*ma->xz;
	mr->yz += m*ma->yz;
	mr->xxx += m*ma->xxx;
	mr->xyy += m*ma->xyy;
	mr->xxy += m*ma->xxy;
	mr->yyy += m*ma->yyy;
	mr->xxz += m*ma->xxz;
	mr->yyz += m*ma->yyz;
	mr->xyz += m*ma->xyz;
	mr->xxxx += m*ma->xxxx;
	mr->xyyy += m*ma->xyyy;
	mr->xxxy += m*ma->xxxy;
	mr->yyyy += m*ma->yyyy;
	mr->xxxz += m*ma->xxxz;
	mr->yyyz += m*ma->yyyz;
	mr->xxyy += m*ma->xxyy;
	mr->xxyz += m*ma->xxyz;
	mr->xyyz += m*ma->xyyz;
	}


/*
 ** This function subtracts the complete moment ma from the complete moment mc
 */
void momSubMomc(MOMC *mc,MOMC *ma)
{
	mc->m -= ma->m;
	mc->xx -= ma->xx;
	mc->yy -= ma->yy;
	mc->xy -= ma->xy;
	mc->xz -= ma->xz;
	mc->yz -= ma->yz;
	mc->xxx -= ma->xxx;
	mc->xyy -= ma->xyy;
	mc->xxy -= ma->xxy;
	mc->yyy -= ma->yyy;
	mc->xxz -= ma->xxz;
	mc->yyz -= ma->yyz;
	mc->xyz -= ma->xyz;
	mc->xxxx -= ma->xxxx;
	mc->xyyy -= ma->xyyy;
	mc->xxxy -= ma->xxxy;
	mc->yyyy -= ma->yyyy;
	mc->xxxz -= ma->xxxz;
	mc->yyyz -= ma->yyyz;
	mc->xxyy -= ma->xxyy;
	mc->xxyz -= ma->xxyz;
	mc->xyyz -= ma->xyyz;
	mc->zz -= ma->zz;
	mc->xzz -= ma->xzz;
	mc->yzz -= ma->yzz;
	mc->zzz -= ma->zzz;
	mc->xxzz -= ma->xxzz;
	mc->xyzz -= ma->xyzz;
	mc->xzzz -= ma->xzzz;
	mc->yyzz -= ma->yyzz;
	mc->yzzz -= ma->yzzz;
	mc->zzzz -= ma->zzzz;
	}


/*
 ** This function subtracts the reduced moment ma from the reduced moment mr
 */
void momSubMomr(MOMR *mr,MOMR *ma)
{
	mr->m -= ma->m;
	mr->xx -= ma->xx;
	mr->yy -= ma->yy;
	mr->xy -= ma->xy;
	mr->xz -= ma->xz;
	mr->yz -= ma->yz;
	mr->xxx -= ma->xxx;
	mr->xyy -= ma->xyy;
	mr->xxy -= ma->xxy;
	mr->yyy -= ma->yyy;
	mr->xxz -= ma->xxz;
	mr->yyz -= ma->yyz;
	mr->xyz -= ma->xyz;
	mr->xxxx -= ma->xxxx;
	mr->xyyy -= ma->xyyy;
	mr->xxxy -= ma->xxxy;
	mr->yyyy -= ma->yyyy;
	mr->xxxz -= ma->xxxz;
	mr->yyyz -= ma->yyyz;
	mr->xxyy -= ma->xxyy;
	mr->xxyz -= ma->xxyz;
	mr->xyyz -= ma->xyyz;
	}


/*
 ** This function calculates a complete multipole from a single
 ** particle at position <x,y,z> from the center of mass.
 ** The strange order of evaluation reduces the number of 
 ** multiplications to a minimum.
 ** <x,y,z> := d := r(particle) - rcm.
 **
 ** OpCount (*,+) = (34,0)
 **
 */
void momMakeMomc(MOMC *mc,momFloat m,momFloat x,momFloat y,momFloat z)
{
	momFloat tx,ty,tz;

	mc->m = m;
	/*
	 ** Calculate the Quadrupole Moment.
	 */
	tx = m*x;
	ty = m*y;
	mc->xy = tx*y;
	mc->xz = tx*z;
	mc->yz = ty*z;
	tx *= x;
	ty *= y;
	tz = m*z*z;
	mc->xx = tx;
	mc->yy = ty;
	mc->zz = tz;
	/*
	 ** Calculate the Octopole Moment.
	 */
	mc->xxy = tx*y;
	mc->xxz = tx*z;
	mc->yyz = ty*z;
	mc->xyy = ty*x;
	mc->xzz = tz*x;
	mc->yzz = tz*y;
	mc->xyz = mc->xy*z;
	tx *= x;
	ty *= y;
	tz *= z;
 	mc->xxx = tx;
	mc->yyy = ty;
	mc->zzz = tz;
	/*
	 ** Calculate the Hexadecapole Moment.
	 */
	mc->xxxx = tx*x;
	mc->xxxy = tx*y;
	mc->xxxz = tx*z;
	mc->xyyy = ty*x;
	mc->yyyy = ty*y;
	mc->yyyz = ty*z;
	mc->xzzz = tz*x;
	mc->yzzz = tz*y;
	mc->zzzz = tz*z;
	mc->xxyy = mc->xxy*y;
	mc->xxyz = mc->xxy*z;
	mc->xyyz = mc->yyz*x;
	mc->yyzz = mc->yyz*z;
	mc->xxzz = mc->xzz*x;
	mc->xyzz = mc->xzz*y;
	}


/*
 ** This function calculates a reduced multipole from a single
 ** particle at position <x,y,z> from the center of mass.
 ** The strange order of evaluation reduces the number of 
 ** multiplications to a minimum.
 ** <x,y,z> := d := r(particle) - rcm.
 ** returns: d^2
 **
 ** OpCount (*,+) = (43,18) = ~60
 */
momFloat momMakeMomr(MOMR *mr,momFloat m,momFloat x,momFloat y,momFloat z)
{
	momFloat tx,ty,t,dx,dy;
	momFloat x2 = x*x;
	momFloat y2 = y*y;
	momFloat d2 = x2 + y2 + z*z;

	mr->m = m;
	/*
	 ** Calculate the Quadrupole Moment.
	 */
	tx = m*x;
	ty = m*y;
	mr->xy = tx*y;
	mr->xz = tx*z;
	mr->yz = ty*z;
	tx *= x;
	ty *= y;
	m *= d2;
	t = m/3;
	mr->xx = tx - t;
	mr->yy = ty - t;
	/*
	 ** Calculate the Octopole Moment.
	 */
	t = 0.2*m;
	dx = tx - t;
	dy = ty - t;
	mr->xxy = dx*y;
	mr->xxz = dx*z;
	mr->yyz = dy*z;
	mr->xyy = dy*x;
	mr->xyz = mr->xy*z;
	t *= 3;
 	mr->xxx = (tx - t)*x;
	mr->yyy = (ty - t)*y;
	/*
	 ** Calculate the Hexadecapole Moment.
	 */
	t = m/7;
	mr->xxyz = (tx - t)*y*z;
	mr->xyyz = (ty - t)*x*z;
	dx = (tx - 3*t)*x;
	dy = (ty - 3*t)*y;
	mr->xxxy = dx*y;
	mr->xxxz = dx*z;
	mr->xyyy = dy*x;
	mr->yyyz = dy*z;
	dx = t*(x2 - 0.1*d2);
	dy = t*(y2 - 0.1*d2);
	mr->xxxx = tx*x2 - 6*dx;
	mr->yyyy = ty*y2 - 6*dy;
	mr->xxyy = tx*y2 - dx - dy;

	return(d2);
	}


/*
 ** This function calculates a reduced multipole from a single
 ** particle at position <x,y,z> from the center of mass.
 ** This is the "straight forward" implementation which we 
 ** used in the original version of PKDGRAV. It remains a good
 ** test of more peculiar looking code.
 **
 ** <x,y,z> := d := r(particle) - rcm.
 ** 
 ** OpCount (*,+) = (115,20) = ~135
 */
void momOldMakeMomr(MOMR *mr,momFloat m,momFloat x,momFloat y,momFloat z)
{
	momFloat d2 = x*x + y*y + z*z;
	
	mr->xxxx = m*(x*x*x*x - 6.0/7.0*d2*(x*x - 0.1*d2));
	mr->xyyy = m*(x*y*y*y - 3.0/7.0*d2*x*y);
	mr->xxxy = m*(x*x*x*y - 3.0/7.0*d2*x*y);
	mr->yyyy = m*(y*y*y*y - 6.0/7.0*d2*(y*y - 0.1*d2));
	mr->xxxz = m*(x*x*x*z - 3.0/7.0*d2*x*z);
	mr->yyyz = m*(y*y*y*z - 3.0/7.0*d2*y*z);
	mr->xxyy = m*(x*x*y*y - 1.0/7.0*d2*(x*x + y*y - 0.2*d2));
	mr->xxyz = m*(x*x*y*z - 1.0/7.0*d2*y*z);
	mr->xyyz = m*(x*y*y*z - 1.0/7.0*d2*x*z);
	/*
	 ** Calculate reduced octopole moment...
	 */
	mr->xxx = m*(x*x*x - 0.6*d2*x);
	mr->xyy = m*(x*y*y - 0.2*d2*x);
	mr->xxy = m*(x*x*y - 0.2*d2*y);
	mr->yyy = m*(y*y*y - 0.6*d2*y);
	mr->xxz = m*(x*x*z - 0.2*d2*z);
	mr->yyz = m*(y*y*z - 0.2*d2*z);
	mr->xyz = m*x*y*z;
	/*
	 ** Calculate quadrupole moment...
	 */
	mr->xx = m*(x*x - 1.0/3.0*d2);
	mr->yy = m*(y*y - 1.0/3.0*d2);
	mr->xy = m*x*y;
	mr->xz = m*x*z;
	mr->yz = m*y*z;
	mr->m = m;
	}


/*
 ** This function shifts a complete multipole (MOMC) to a new center of mass.
 ** <x,y,z> := d := rcm(old) - rcm(new).
 **
 ** OpCount ShiftMomc   (*,+) = (111,84)
 **         MakeMomc    (*,+) = (34,0)
 **         MulAddMomc  (*,+) = (32,32)
 **         Total       (*,+) = (177,116) = 293
 */
void momShiftMomc(MOMC *m,momFloat x,momFloat y,momFloat z)
{
	MOMC f;

	momMakeMomc(&f,1,x,y,z);
	/*
	 ** Shift the Hexadecapole.
	 */
	m->xxxx += 4*m->xxx*x + 6*m->xx*f.xx;
	m->yyyy += 4*m->yyy*y + 6*m->yy*f.yy;
	m->zzzz += 4*m->zzz*z + 6*m->zz*f.zz;
	m->xyyy += m->yyy*x + 3*(m->xyy*y + m->yy*f.xy + m->xy*f.yy);
	m->xxxy += m->xxx*y + 3*(m->xxy*x + m->xx*f.xy + m->xy*f.xx);
	m->xxxz += m->xxx*z + 3*(m->xxz*x + m->xx*f.xz + m->xz*f.xx);
	m->yyyz += m->yyy*z + 3*(m->yyz*y + m->yy*f.yz + m->yz*f.yy);
	m->xzzz += m->zzz*x + 3*(m->xzz*z + m->zz*f.xz + m->xz*f.zz);
	m->yzzz += m->zzz*y + 3*(m->yzz*z + m->zz*f.yz + m->yz*f.zz);
	m->xxyy += 2*(m->xxy*y + m->xyy*x) + m->xx*f.yy + m->yy*f.xx + 4*m->xy*f.xy;
	m->xxzz += 2*(m->xxz*z + m->xzz*x) + m->xx*f.zz + m->zz*f.xx + 4*m->xz*f.xz;
	m->yyzz += 2*(m->yyz*z + m->yzz*y) + m->yy*f.zz + m->zz*f.yy + 4*m->yz*f.yz;
	m->xxyz += m->xxy*z + m->xxz*y + m->xx*f.yz + m->yz*f.xx + 2*(m->xyz*x + m->xy*f.xz + m->xz*f.xy);
	m->xyyz += m->xyy*z + m->yyz*x + m->yy*f.xz + m->xz*f.yy + 2*(m->xyz*y + m->xy*f.yz + m->yz*f.xy);
	m->xyzz += m->yzz*x + m->xzz*y + m->zz*f.xy + m->xy*f.zz + 2*(m->xyz*z + m->yz*f.xz + m->xz*f.yz);
	/*
	 ** Now shift the Octopole.
	 */
	m->xxx += 3*m->xx*x;
	m->yyy += 3*m->yy*y;
	m->zzz += 3*m->zz*z;
	m->xyy += 2*m->xy*y + m->yy*x;
	m->xzz += 2*m->xz*z + m->zz*x;
	m->yzz += 2*m->yz*z + m->zz*y;
	m->xxy += 2*m->xy*x + m->xx*y; 
	m->xxz += 2*m->xz*x + m->xx*z;
	m->yyz += 2*m->yz*y + m->yy*z;
	m->xyz += m->xy*z + m->xz*y + m->yz*x;
	/*
	 ** Now deal with the monopole terms.
	 */
	f.m = 0;
	momMulAddMomc(m,m->m,&f);
	}


/*
 ** This function shifts a reduced multipole (MOMR) to a new center of mass.
 ** <x,y,z> := d := rcm(old) - rcm(new).
 **
 ** OpCount ShiftMomr  (*,+) = (128,111)
 **         MakeMomr   (*,+) = (43,18)
 **         MulAddMomr (*,+) = (22,22)
 **         Total      (*,+) = (193,151) = 344
 */
void momShiftMomr(MOMR *m,momFloat x,momFloat y,momFloat z)
{
	MOMR f;
	momFloat t,tx,ty,tz,txx,tyy,txy,tyz,txz;
	const momFloat twosevenths = 2.0/7.0;

	momMakeMomr(&f,1,x,y,z);
	/*
	 ** Calculate the correction terms.
	 */
	tx = 0.4*(m->xx*x + m->xy*y + m->xz*z);
	ty = 0.4*(m->xy*x + m->yy*y + m->yz*z);
	tz = 0.4*(m->xz*x + m->yz*y - (m->xx + m->yy)*z);
	t = tx*x + ty*y + tz*z;
	txx = twosevenths*(m->xxx*x + m->xxy*y + m->xxz*z + 2*(m->xx*f.xx + m->xy*f.xy + m->xz*f.xz) - 0.5*t);
	tyy = twosevenths*(m->xyy*x + m->yyy*y + m->yyz*z + 2*(m->xy*f.xy + m->yy*f.yy + m->yz*f.yz) - 0.5*t);
	txy = twosevenths*(m->xxy*x + m->xyy*y + m->xyz*z + m->xy*(f.xx + f.yy) + (m->xx + m->yy)*f.xy + m->yz*f.xz + m->xz*f.yz);
	tyz = twosevenths*(m->xyz*x + m->yyz*y - (m->xxy + m->yyy)*z - m->yz*f.xx - m->xx*f.yz + m->xz*f.xy + m->xy*f.xz);
	txz = twosevenths*(m->xxz*x + m->xyz*y - (m->xxx + m->xyy)*z - m->xz*f.yy - m->yy*f.xz + m->yz*f.xy + m->xy*f.yz);
	/*
	 ** Shift the Hexadecapole.
	 */
	m->xxxx += 4*m->xxx*x + 6*(m->xx*f.xx - txx);
	m->yyyy += 4*m->yyy*y + 6*(m->yy*f.yy - tyy);
	m->xyyy += m->yyy*x + 3*(m->xyy*y + m->yy*f.xy + m->xy*f.yy - txy);
	m->xxxy += m->xxx*y + 3*(m->xxy*x + m->xx*f.xy + m->xy*f.xx - txy);
	m->xxxz += m->xxx*z + 3*(m->xxz*x + m->xx*f.xz + m->xz*f.xx - txz);
	m->yyyz += m->yyy*z + 3*(m->yyz*y + m->yy*f.yz + m->yz*f.yy - tyz);
	m->xxyy += 2*(m->xxy*y + m->xyy*x) + m->xx*f.yy + m->yy*f.xx + 4*m->xy*f.xy - txx - tyy;
	m->xxyz += m->xxy*z + m->xxz*y + m->xx*f.yz + m->yz*f.xx + 2*(m->xyz*x + m->xy*f.xz + m->xz*f.xy) - tyz;
	m->xyyz += m->xyy*z + m->yyz*x + m->yy*f.xz + m->xz*f.yy + 2*(m->xyz*y + m->xy*f.yz + m->yz*f.xy) - txz;
	/*
	 ** Now shift the Octopole.
	 */
	m->xxx += 3*(m->xx*x - tx);
	m->xyy += 2*m->xy*y + m->yy*x - tx;
	m->yyy += 3*(m->yy*y - ty);
	m->xxy += 2*m->xy*x + m->xx*y - ty;
	m->xxz += 2*m->xz*x + m->xx*z - tz;
	m->yyz += 2*m->yz*y + m->yy*z - tz;
	m->xyz += m->xy*z + m->xz*y + m->yz*x;
	/*
	 ** Now deal with the monopole terms.
	 */
	f.m = 0;
	momMulAddMomr(m,m->m,&f);
	}


/*
 ** This function converts a complete multipole (MOMC) to a reduced one (MOMR).
 */
void momReduceMomc(MOMC *mc,MOMR *mr)
{
	momFloat  t,tx,ty,tz,txx,txy,txz,tyy,tyz,tzz;

	/*
	 ** First reduce Hexadecapole.
	 */
	txx = (mc->xxxx + mc->xxyy + mc->xxzz)/7;
	txy = (mc->xxxy + mc->xyyy + mc->xyzz)/7;
	txz = (mc->xxxz + mc->xyyz + mc->xzzz)/7;
	tyy = (mc->xxyy + mc->yyyy + mc->yyzz)/7;
	tyz = (mc->xxyz + mc->yyyz + mc->yzzz)/7;
	tzz = (mc->xxzz + mc->yyzz + mc->zzzz)/7;
	t = 0.1*(txx + tyy + tzz);
	mr->xxxx = mc->xxxx - 6*(txx - t);
	mr->xyyy = mc->xyyy - 3*txy;
	mr->xxxy = mc->xxxy - 3*txy;
	mr->yyyy = mc->yyyy - 6*(tyy - t);
	mr->xxxz = mc->xxxz - 3*txz;
	mr->yyyz = mc->yyyz - 3*tyz;
	mr->xxyy = mc->xxyy - (txx + tyy - 2*t);
	mr->xxyz = mc->xxyz - tyz;
	mr->xyyz = mc->xyyz - txz;	
	/*
	 ** Now reduce the Octopole.
	 */
	tx = (mc->xxx + mc->xyy + mc->xzz)/5;
	ty = (mc->xxy + mc->yyy + mc->yzz)/5;
	tz = (mc->xxz + mc->yyz + mc->zzz)/5;
	mr->xxx = mc->xxx - 3*tx;
	mr->xyy = mc->xyy - tx;
	mr->xxy = mc->xxy - ty;
	mr->yyy = mc->yyy - 3*ty;
	mr->xxz = mc->xxz - tz;
	mr->yyz = mc->yyz - tz;
	mr->xyz = mc->xyz;
	/*
	 ** Now reduce the Quadrupole.
	 */
	t = (mc->xx + mc->yy + mc->zz)/3;
	mr->xx = mc->xx - t;
	mr->yy = mc->yy - t;
	mr->xy = mc->xy;
	mr->xz = mc->xz;
	mr->yz = mc->yz;
	/*
	 ** Finally the mass remains the same.
	 */
	mr->m = mc->m;
	}


/*
 ** This is a new fast version of QEVAL which evaluates
 ** the interaction due to the reduced moment 'm'.
 ** This version is nearly two times as fast as a naive
 ** implementation.
 **
 ** OpCount = (*,+) = (103,72) = 175 - 8 = 167
 */
void momEvalMomr(MOMR *m,momFloat dir,momFloat x,momFloat y,momFloat z,
				 momFloat *fPot,momFloat *ax,momFloat *ay,momFloat *az)
{
	const momFloat onethird = 1.0/3.0;
	momFloat xx,xy,xz,yy,yz,zz;
	momFloat xxx,xxy,xxz,xyy,yyy,yyz,xyz;
	momFloat tx,ty,tz,dir2,g2,g3,g4;

	dir = -dir;
	dir2 = dir*dir;
	g2 = 3*dir*dir2*dir2;
	g3 = -5*g2*dir2;
	g4 = -7*g3*dir2;
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
	dir2 *= -(dir + 5*g2 + 7*g3 + 9*g4);
	*fPot += dir + g2 + g3 + g4;
	*ax += xx + xxx + tx + x*dir2;
	*ay += xy + xxy + ty + y*dir2;
	*az += xz + xxz + tz + z*dir2;
	}


void momMomr2Momc(MOMR *ma,MOMC *mc)
{
	mc->m = ma->m;
	mc->xx = ma->xx;
	mc->yy = ma->yy;
	mc->xy = ma->xy;
	mc->xz = ma->xz;
	mc->yz = ma->yz;
	mc->xxx = ma->xxx;
	mc->xyy = ma->xyy;
	mc->xxy = ma->xxy;
	mc->yyy = ma->yyy;
	mc->xxz = ma->xxz;
	mc->yyz = ma->yyz;
	mc->xyz = ma->xyz;
	mc->xxxx = ma->xxxx;
	mc->xyyy = ma->xyyy;
	mc->xxxy = ma->xxxy;
	mc->yyyy = ma->yyyy;
	mc->xxxz = ma->xxxz;
	mc->yyyz = ma->yyyz;
	mc->xxyy = ma->xxyy;
	mc->xxyz = ma->xxyz;
	mc->xyyz = ma->xyyz;
	mc->zz = -(ma->xx + ma->yy);
	mc->xzz = -(ma->xxx + ma->xyy);
	mc->yzz = -(ma->xxy + ma->yyy);
	mc->zzz = -(ma->xxz + ma->yyz);
	mc->xxzz = -(ma->xxxx + ma->xxyy);
	mc->xyzz = -(ma->xxxy + ma->xyyy);
	mc->xzzz = -(ma->xxxz + ma->xyyz);
	mc->yyzz = -(ma->xxyy + ma->yyyy);
	mc->yzzz = -(ma->xxyz + ma->yyyz);
	mc->zzzz = -(mc->xxzz + mc->yyzz);
	}


void momPrintMomc(MOMC *m) {
	printf("MOMC:%20.15g\n",(double)m->m);
	printf("  xx:%20.15g   yy:%20.15g   zz:%20.15g\n",(double)m->xx,(double)m->yy,(double)m->zz);
    printf("  xy:%20.15g   yz:%20.15g   xz:%20.15g\n",(double)m->xy,(double)m->yz,(double)m->xz);
	printf(" xxx:%20.15g  xyy:%20.15g  xzz:%20.15g\n",(double)m->xxx,(double)m->xyy,(double)m->xzz);
    printf(" xxy:%20.15g  yyy:%20.15g  yzz:%20.15g\n",(double)m->xxy,(double)m->yyy,(double)m->yzz);
	printf(" xxz:%20.15g  yyz:%20.15g  zzz:%20.15g\n",(double)m->xxz,(double)m->yyz,(double)m->zzz);
	printf(" xyz:%20.15g\n",(double)m->xyz);
	printf("xxxx:%20.15g xxxy:%20.15g xxxz:%20.15g\n",(double)m->xxxx,(double)m->xxxy,(double)m->xxxz);
	printf("xyyy:%20.15g yyyy:%20.15g yyyz:%20.15g\n",(double)m->xyyy,(double)m->yyyy,(double)m->yyyz);
	printf("xzzz:%20.15g yzzz:%20.15g zzzz:%20.15g\n",(double)m->xzzz,(double)m->yzzz,(double)m->zzzz);
	printf("xxyy:%20.15g xxyz:%20.15g xyyz:%20.15g\n",(double)m->xxyy,(double)m->xxyz,(double)m->xyyz);
	printf("yyzz:%20.15g xxzz:%20.15g xyzz:%20.15g\n",(double)m->yyzz,(double)m->xxzz,(double)m->xyzz);
	}


void momPrintMomr(MOMR *m) {
	printf("MOMR:%20.15g\n",(double)m->m);
	printf("  xx:%20.15g   yy:%20.15g\n",(double)m->xx,(double)m->yy);
    printf("  xy:%20.15g   yz:%20.15g   xz:%20.15g\n",(double)m->xy,(double)m->yz,(double)m->xz);
	printf(" xxx:%20.15g  xyy:%20.15g\n",(double)m->xxx,(double)m->xyy);
    printf(" xxy:%20.15g  yyy:%20.15g\n",(double)m->xxy,(double)m->yyy);
	printf(" xxz:%20.15g  yyz:%20.15g\n",(double)m->xxz,(double)m->yyz);
	printf(" xyz:%20.15g\n",(double)m->xyz);
	printf("xxxx:%20.15g xxxy:%20.15g xxxz:%20.15g\n",(double)m->xxxx,(double)m->xxxy,(double)m->xxxz);
	printf("xyyy:%20.15g yyyy:%20.15g yyyz:%20.15g\n",(double)m->xyyy,(double)m->yyyy,(double)m->yyyz);
	printf("xxyy:%20.15g xxyz:%20.15g xyyz:%20.15g\n",(double)m->xxyy,(double)m->xxyz,(double)m->xyyz);
	}
