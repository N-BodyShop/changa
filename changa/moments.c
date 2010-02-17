/*
 * Mulitpole moment routines.
 * Main author: Joachim Stadel as part of PKGRAV
 */
#include <stdio.h>
#include <math.h>
#include "moments.h"

void momClearLocr(LOCR *l) {
    l->m = 0;
    l->x = 0;
    l->y = 0;
    l->z = 0;
    l->xx = 0;
    l->xy = 0;
    l->yy = 0;
    l->xz = 0;
    l->yz = 0;
    l->xxx = 0;
    l->xxy = 0;
    l->xyy = 0;
    l->yyy = 0;
    l->xxz = 0;
    l->xyz = 0;
    l->yyz = 0;
    l->xxxx = 0;
    l->xxxy = 0;
    l->xxyy = 0;
    l->xyyy = 0;
    l->yyyy = 0;
    l->xxxz = 0;
    l->xxyz = 0;
    l->xyyz = 0;
    l->yyyz = 0;
    l->xxxxx = 0;
    l->xxxxy = 0;
    l->xxxyy = 0;
    l->xxyyy = 0;
    l->xyyyy = 0;
    l->yyyyy = 0;
    l->xxxxz = 0;
    l->xxxyz = 0;
    l->xxyyz = 0;
    l->xyyyz = 0;
    l->yyyyz = 0;
    }

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
 ** This function shifts a reduced local expansion (LOCR) to a new center of expansion.
 ** <x,y,z> := d := rexp(new) - rexp(old).
 **
 ** Op Count (*,+) = (159,173) = 332
 */
/* IBM brain damage */
#undef hz
double momShiftLocr(LOCR *l,momFloat x,momFloat y,momFloat z) {
    const momFloat onethird = 1.0/3.0;
    momFloat hx,hy,hz,tx,ty,tz;
    momFloat L,Lx,Ly,Lz,Lxx,Lxy,Lxz,Lyy,Lyz,Lxxx,Lxxy,Lxxz,Lxyy,Lxyz,Lyyy,Lyyz;
    momFloat Lxxxx,Lxxxy,Lxxyy,Lxyyy,Lyyyy,Lxxxz,Lxxyz,Lxyyz,Lyyyz;

    hx = 0.5*x;
    hy = 0.5*y;
    hz = 0.5*z;
    tx = onethird*x;
    ty = onethird*y;
    tz = onethird*z;

    L = l->x*x + l->y*y + l->z*z;
    l->m += L;

    Lx = l->xx*x + l->xy*y + l->xz*z;
    Ly = l->xy*x + l->yy*y + l->yz*z;
    Lz = l->xz*x + l->yz*y - (l->xx + l->yy)*z;
    L = Lx*hx + Ly*hy + Lz*hz;
    l->x += Lx;
    l->y += Ly;
    l->z += Lz;
    l->m += L;

    Lxx = l->xxx*x + l->xxy*y + l->xxz*z;
    Lxy = l->xxy*x + l->xyy*y + l->xyz*z;
    Lxz = l->xxz*x + l->xyz*y - (l->xxx + l->xyy)*z;
    Lyy = l->xyy*x + l->yyy*y + l->yyz*z;
    Lyz = l->xyz*x + l->yyz*y - (l->xxy + l->yyy)*z;
    Lx = Lxx*hx + Lxy*hy + Lxz*hz;
    Ly = Lxy*hx + Lyy*hy + Lyz*hz;
    Lz = Lxz*hx + Lyz*hy - (Lxx + Lyy)*hz;
    L = Lx*tx + Ly*ty + Lz*tz;
    l->xx += Lxx;
    l->xy += Lxy;
    l->xz += Lxz;
    l->yy += Lyy;
    l->yz += Lyz;
    l->x += Lx;
    l->y += Ly;
    l->z += Lz;
    l->m += L;

    Lxxx = l->xxxx*x + l->xxxy*y + l->xxxz*z;
    Lxxy = l->xxxy*x + l->xxyy*y + l->xxyz*z;
    Lxxz = l->xxxz*x + l->xxyz*y - (l->xxxx + l->xxyy)*z;
    Lxyy = l->xxyy*x + l->xyyy*y + l->xyyz*z;
    Lxyz = l->xxyz*x + l->xyyz*y - (l->xxxy + l->xyyy)*z;
    Lyyy = l->xyyy*x + l->yyyy*y + l->yyyz*z;
    Lyyz = l->xyyz*x + l->yyyz*y - (l->xxyy + l->yyyy)*z;
    Lxx = Lxxx*hx + Lxxy*hy + Lxxz*hz;
    Lxy = Lxxy*hx + Lxyy*hy + Lxyz*hz;
    Lxz = Lxxz*hx + Lxyz*hy - (Lxxx + Lxyy)*hz;
    Lyy = Lxyy*hx + Lyyy*hy + Lyyz*hz;
    Lyz = Lxyz*hx + Lyyz*hy - (Lxxy + Lyyy)*hz;
    Lx = Lxx*tx + Lxy*ty + Lxz*tz;
    Ly = Lxy*tx + Lyy*ty + Lyz*tz;
    Lz = Lxz*tx + Lyz*ty - (Lxx + Lyy)*tz;
    L = Lx*hx + Ly*hy + Lz*hz;
    l->xxx += Lxxx;
    l->xxy += Lxxy;
    l->xxz += Lxxz;
    l->xyy += Lxyy;
    l->xyz += Lxyz;
    l->yyy += Lyyy;
    l->yyz += Lyyz;
    l->xx += Lxx;
    l->xy += Lxy;
    l->xz += Lxz;
    l->yy += Lyy;
    l->yz += Lyz;
    l->x += Lx;
    l->y += Ly;
    l->z += Lz;
    l->m += 0.5*L;

    Lxxxx = l->xxxxx*x + l->xxxxy*y + l->xxxxz*z;
    Lxxxy = l->xxxxy*x + l->xxxyy*y + l->xxxyz*z;
    Lxxyy = l->xxxyy*x + l->xxyyy*y + l->xxyyz*z;
    Lxyyy = l->xxyyy*x + l->xyyyy*y + l->xyyyz*z;
    Lyyyy = l->xyyyy*x + l->yyyyy*y + l->yyyyz*z;
    Lxxxz = l->xxxxz*x + l->xxxyz*y - (l->xxxxx + l->xxxyy)*z;
    Lxxyz = l->xxxyz*x + l->xxyyz*y - (l->xxxxy + l->xxyyy)*z;
    Lxyyz = l->xxyyz*x + l->xyyyz*y - (l->xxxyy + l->xyyyy)*z;
    Lyyyz = l->xyyyz*x + l->yyyyz*y - (l->xxyyy + l->yyyyy)*z;
    Lxxx = Lxxxx*hx + Lxxxy*hy + Lxxxz*hz;
    Lxxy = Lxxxy*hx + Lxxyy*hy + Lxxyz*hz;
    Lxyy = Lxxyy*hx + Lxyyy*hy + Lxyyz*hz;
    Lyyy = Lxyyy*hx + Lyyyy*hy + Lyyyz*hz;
    Lxxz = Lxxxz*hx + Lxxyz*hy - (Lxxxx + Lxxyy)*hz;
    Lxyz = Lxxyz*hx + Lxyyz*hy - (Lxxxy + Lxyyy)*hz;
    Lyyz = Lxyyz*hx + Lyyyz*hy - (Lxxyy + Lyyyy)*hz;
    Lxx = Lxxx*tx + Lxxy*ty + Lxxz*tz;
    Lxy = Lxxy*tx + Lxyy*ty + Lxyz*tz;
    Lyy = Lxyy*tx + Lyyy*ty + Lyyz*tz;
    Lxz = Lxxz*tx + Lxyz*ty - (Lxxx + Lxyy)*tz;
    Lyz = Lxyz*tx + Lyyz*ty - (Lxxy + Lyyy)*tz;
    Lx = Lxx*hx + Lxy*hy + Lxz*hz;
    Ly = Lxy*hx + Lyy*hy + Lyz*hz;
    Lz = Lxz*hx + Lyz*hy - (Lxx + Lyy)*hz;
    L = Lx*hx + Ly*hy + Lz*hz;
    l->xxxx += Lxxxx;
    l->xxxy += Lxxxy;
    l->xxyy += Lxxyy;
    l->xyyy += Lxyyy;
    l->yyyy += Lyyyy;
    l->xxxz += Lxxxz;
    l->xxyz += Lxxyz;
    l->xyyz += Lxyyz;
    l->yyyz += Lyyyz;
    l->xxx += Lxxx;
    l->xxy += Lxxy;
    l->xxz += Lxxz;
    l->xyy += Lxyy;
    l->xyz += Lxyz;
    l->yyy += Lyyy;
    l->yyz += Lyyz;
    l->xx += Lxx;
    l->xy += Lxy;
    l->xz += Lxz;
    l->yy += Lyy;
    l->yz += Lyz;
    l->x += 0.5*Lx;
    l->y += 0.5*Ly;
    l->z += 0.5*Lz;
    l->m += 0.2*L;

    return 332.0;
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
void momEvalMomr(MOMR *m,momFloat dir0,momFloat x,momFloat y,momFloat z,
				 momFloat *fPot,momFloat *ax,momFloat *ay,momFloat *az)
{
	const momFloat onethird = 1.0/3.0;
	momFloat xx,xy,xz,yy,yz,zz;
	momFloat xxx,xxy,xxz,xyy,yyy,yyz,xyz;
	momFloat tx,ty,tz,dir2,g2,g3,g4;

	momFloat dir = -dir0;
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

/*
** Op Count = (*,+) = (302,207) = 509
*/
double momLocrAddMomr5(LOCR *l,MOMR *m,momFloat dir,momFloat x,momFloat y,momFloat z,double *tax,double *tay,double *taz) {
    const momFloat onethird = 1.0/3.0;
    momFloat xx,xy,xz,yy,yz,zz;
    momFloat xxx,xxy,xyy,yyy,xxz,xyz,yyz;
    momFloat Ax,Ay,Az,A,Bxx,Bxy,Byy,Bxz,Byz,Bx,By,Bz,B,Cx,Cy,Cz,C;
    momFloat R1,R2,R3,T2,T3;
    momFloat g0,g1,g2,g3,g4,g5;
    momFloat g4xx,g4yy,g5xx,g5yy,fxx,fyy;

    x *= dir;
    y *= dir;
    z *= dir;
    g0 = -dir;
    g1 = -g0*dir;
    g2 = -3*g1*dir;
    g3 = -5*g2*dir;
    g4 = -7*g3*dir;
    g5 = -9*g4*dir;
    /*
    ** Calculate the funky distance terms.
    */
    xx = 0.5*x*x;
    xy = x*y;
    yy = 0.5*y*y;
    xz = x*z;
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

    Bxx = x*m->xxx + y*m->xxy + z*m->xxz;
    Bxy = x*m->xxy + y*m->xyy + z*m->xyz;
    Byy = x*m->xyy + y*m->yyy + z*m->yyz;
    Bxz = x*m->xxz + y*m->xyz - z*(m->xxx + m->xyy);
    Byz = x*m->xyz + y*m->yyz - z*(m->xxy + m->yyy);

    Cx = m->xxxx*xxx + m->xyyy*yyy + m->xxxy*xxy + m->xxxz*xxz + m->xxyy*xyy + m->xxyz*xyz + m->xyyz*yyz;
    Cy = m->xyyy*xyy + m->xxxy*xxx + m->yyyy*yyy + m->yyyz*yyz + m->xxyy*xxy + m->xxyz*xxz + m->xyyz*xyz;
    Cz = -m->xxxx*xxz - (m->xyyy + m->xxxy)*xyz - m->yyyy*yyz + m->xxxz*xxx + m->yyyz*yyy - m->xxyy*(xxz + yyz) + m->xxyz*xxy + m->xyyz*xyy;

    Ax = x*m->xx + y*m->xy + z*m->xz;
    Ay = x*m->xy + y*m->yy + z*m->yz;
    Az = x*m->xz + y*m->yz - z*(m->xx + m->yy);

    Bx = 0.5*(x*Bxx + y*Bxy + z*Bxz);
    By = 0.5*(x*Bxy + y*Byy + z*Byz);
    Bz = 0.5*(x*Bxz + y*Byz - z*(Bxx + Byy));

    C = 0.25*(x*Cx + y*Cy + z*Cz);

    A = 0.5*(x*Ax + y*Ay + z*Az);

    B = onethird*(x*Bx + y*By + z*Bz);

    xx = x*x;
    yy = y*y;

    l->m += g0*m->m + g2*A - g3*B + g4*C;
    R1 = g1*m->m + g3*A - g4*B + g5*C;
    R2 = g2*m->m + g4*A - g5*B;
    R3 = g3*m->m + g5*A;

    g1 *= dir;
    g2 *= dir;
    g3 *= dir;

    T2 = g1*m->m + g3*A;
    T3 = g2*m->m;

    *tax = -(g2*Ax - g3*Bx + x*R1);
    *tay = -(g2*Ay - g3*By + y*R1);
    *taz = -(g2*Az - g3*Bz + z*R1);

    g2 *= dir;

    g4xx = g4*xx;
    g4yy = g4*yy;

    l->xxxx += (3*g2 + (6*g3 + g4xx)*xx)*m->m;
    l->yyyy += (3*g2 + (6*g3 + g4yy)*yy)*m->m;
    fxx = (3*g3 + g4xx)*m->m;
    l->xxxy += fxx*xy;
    l->xxxz += fxx*xz;
    fyy = (3*g3 + g4yy)*m->m;
    l->xyyy += fyy*xy;
    l->yyyz += fyy*yz;
    fxx = (g3 + g4xx);
    l->xxyz += fxx*yz*m->m;
    l->xxyy += (g2 + g3*xx + fxx*yy)*m->m;
    l->xyyz += (g3 + g4yy)*xz*m->m;

    g4 *= dir;

    T2 -= g4*B;
    T3 += g4*A;

    *tax -= g4*Cx;
    *tay -= g4*Cy;
    *taz -= g4*Cz;
    l->x -= *tax;
    l->y -= *tay;
    l->z -= *taz;

    l->xx += g2*m->xx + T2 + R2*xx + 2*g3*Ax*x - 2*g4*Bx*x;
    l->yy += g2*m->yy + T2 + R2*yy + 2*g3*Ay*y - 2*g4*By*y;
    l->xy += g2*m->xy + R2*xy + g3*(Ax*y + Ay*x) - g4*(Bx*y + By*x);
    l->xz += g2*m->xz + R2*xz + g3*(Ax*z + Az*x) - g4*(Bx*z + Bz*x);
    l->yz += g2*m->yz + R2*yz + g3*(Ay*z + Az*y) - g4*(By*z + Bz*y);

    g3 *= dir;

    g5xx = g5*xx;
    g5yy = g5*yy;

    l->xx -= g3*Bxx;
    l->xy -= g3*Bxy;
    l->yy -= g3*Byy;
    l->xz -= g3*Bxz;
    l->yz -= g3*Byz;

    fxx = T3 + R3*xx;
    fyy = T3 + R3*yy;

    l->xxy += fxx*y + g4*(2*xy*Ax + xx*Ay) + g3*(Ay + 2*m->xy*x + m->xx*y);
    l->xxz += fxx*z + g4*(2*xz*Ax + xx*Az) + g3*(Az + 2*m->xz*x + m->xx*z);
    l->xyy += fyy*x + g4*(2*xy*Ay + yy*Ax) + g3*(Ax + 2*m->xy*y + m->yy*x);
    l->yyz += fyy*z + g4*(2*yz*Ay + yy*Az) + g3*(Az + 2*m->yz*y + m->yy*z);
    l->xyz += R3*xy*z + g4*(yz*Ax + xz*Ay + xy*Az) + g3*(m->xy*z + m->xz*y + m->yz*x);

    fxx += 2*T3;
    fyy += 2*T3;

    l->xxx += fxx*x + 3*(g4*xx*Ax + g3*(Ax + m->xx*x));
    l->yyy += fyy*y + 3*(g4*yy*Ay + g3*(Ay + m->yy*y));

    fxx = 3*g3 + (6*g4 + g5xx)*xx;
    fyy = 3*g3 + (6*g4 + g5yy)*yy;

    x *= m->m;
    y *= m->m;
    z *= m->m;

    l->xxxxx += (15*g3 + (10*g4 + g5xx)*xx)*x;
    l->yyyyy += (15*g3 + (10*g4 + g5yy)*yy)*y;
    l->xxxyy += (3*g3 + 3*g4*yy + (g4 + g5yy)*xx)*x;
    l->xxyyy += (3*g3 + 3*g4*xx + (g4 + g5xx)*yy)*y;
    l->xxyyz += (g3 + g4*xx + (g4 + g5xx)*yy)*z;
    l->xxxxy += fxx*y;
    l->xxxxz += fxx*z;
    l->xyyyy += fyy*x;
    l->yyyyz += fyy*z;
    l->xxxyz += (3*g4 + g5xx)*xy*z;
    l->xyyyz += (3*g4 + g5yy)*xy*z;

    return 509.0; /* a guess right now */
    }


void momEvalLocr(LOCR *l,momFloat x,momFloat y,momFloat z,
		 momFloat *fPot,momFloat *ax,momFloat *ay,momFloat *az) {
    const momFloat onethird = 1.0/3.0;
    momFloat xx,xy,xz,yy,yz,zz,xxx,xxz,yyy,yyz,xxy,xyy,xyz,xxxx,xxxy,xxyy,xyyy,yyyy,xxxz,xxyz,xyyz,yyyz;
    momFloat g1,A,Ax,Ay,Az,B,Bx,By,Bz,C,Cx,Cy,Cz,D,Dx,Dy,Dz;

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
    xxxx = 0.25*(x*xxx - z*xxz);
    xxxy = y*xxx;
    xxyy = xx*yy - 2*onethird*zz*zz;
    xyyy = x*yyy;
    yyyy = 0.25*(y*yyy - z*yyz);
    xxxz = onethird*x*z*xx;
    xxyz = y*xxz;
    xyyz = x*yyz;
    yyyz = y*yyz;
    /*
     ** Now calculate the interaction.
     */
    Dx = l->xxxxx*xxxx + l->xxxxy*xxxy + l->xxxyy*xxyy + l->xxyyy*xyyy + l->xyyyy*yyyy + l->xxxxz*xxxz + l->xxxyz*xxyz + l->xxyyz*xyyz + l->xyyyz*yyyz;
    Dy = l->xxxxy*xxxx + l->xxxyy*xxxy + l->xxyyy*xxyy + l->xyyyy*xyyy + l->yyyyy*yyyy + l->xxxyz*xxxz + l->xxyyz*xxyz + l->xyyyz*xyyz + l->yyyyz*yyyz;
    Dz = l->xxxxz*xxxx + l->xxxyz*xxxy + l->xxyyz*xxyy + l->xyyyz*xyyy + l->yyyyz*yyyy
	 - l->xxxxx*xxxz - l->xxxxy*xxyz - l->xxxyy*(xxxz + xyyz) - l->xxyyy*(xxyz + yyyz) + l->xyyyy*xyyz + l->yyyyy*yyyz;
    Cx = l->xxxx*xxx + l->xyyy*yyy + l->xxxy*xxy + l->xxxz*xxz + l->xxyy*xyy + l->xxyz*xyz + l->xyyz*yyz;
    Cy = l->xyyy*xyy + l->xxxy*xxx + l->yyyy*yyy + l->yyyz*yyz + l->xxyy*xxy + l->xxyz*xxz + l->xyyz*xyz;
    Cz = -l->xxxx*xxz - (l->xyyy + l->xxxy)*xyz - l->yyyy*yyz + l->xxxz*xxx + l->yyyz*yyy - l->xxyy*(xxz + yyz) + l->xxyz*xxy + l->xyyz*xyy;
    Bx = l->xxx*xx + l->xyy*yy + l->xxy*xy + l->xxz*xz + l->xyz*yz;
    By = l->xyy*xy + l->xxy*xx + l->yyy*yy + l->yyz*yz + l->xyz*xz;
    Bz = -(l->xxx + l->xyy)*xz - (l->xxy + l->yyy)*yz + l->xxz*xx + l->yyz*yy + l->xyz*xy;
    Ax = l->xx*x + l->xy*y + l->xz*z;
    Ay = l->yy*y + l->xy*x + l->yz*z;
    Az = -(l->xx + l->yy)*z + l->xz*x + l->yz*y;
    D = 0.2*(Dx*x + Dy*y + Dz*z);
    C = 0.25*(Cx*x + Cy*y + Cz*z);
    B = onethird*(Bx*x + By*y + Bz*z);
    A = 0.5*(Ax*x + Ay*y + Az*z);
    g1 = x*l->x + y*l->y + z*l->z;
    *ax -= l->x + Ax + Bx + Cx + Dx;
    *ay -= l->y + Ay + By + Cy + Dy;
    *az -= l->z + Az + Bz + Cz + Dz;
    *fPot += l->m + g1 + A + B + C + D;
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
