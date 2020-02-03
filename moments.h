#ifndef MOMENTS_INCLUDED
#define MOMENTS_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

#include "cosmoType.h"

/*************************************************
 *
 *  The precision of the non-scaled moments (MOMR,
 *  MOMC, and LOCR) are to always be double (or better) to
 *  prevent known over/underflow issues.
 *  This is the meaning of momFloat.
 *
 *              Do not change this
 *
 *************************************************/
#ifdef QUAD
typedef long double momFloat;
#define sqrt(x)	sqrtl(x)
#else
typedef double momFloat;
#endif

/**
 ** @brief moment tensor components for reduced multipoles.
 */
typedef struct momReduced {
	momFloat m;
	momFloat xx,yy,xy,xz,yz;
	momFloat xxx,xyy,xxy,yyy,xxz,yyz,xyz;
	momFloat xxxx,xyyy,xxxy,yyyy,xxxz,yyyz,xxyy,xxyz,xyyz;
	} MOMR;

/**
 ** @brief moment tensor components for complete multipoles.
 */
typedef struct momComplete {
	momFloat m;
	momFloat xx,yy,xy,xz,yz;
	momFloat xxx,xyy,xxy,yyy,xxz,yyz,xyz;
	momFloat xxxx,xyyy,xxxy,yyyy,xxxz,yyyz,xxyy,xxyz,xyyz;
	momFloat zz;
	momFloat xzz,yzz,zzz;
	momFloat xxzz,xyzz,xzzz,yyzz,yzzz,zzzz;
	} MOMC;

/**
 ** @brief moment tensor components for reduced local expansion.
 ** note that we have the 5th-order terms here now!
 */
typedef struct locReduced {
    momFloat m;
    momFloat x,y,z;
    momFloat xx,xy,yy,xz,yz;
    momFloat xxx,xxy,xyy,yyy,xxz,xyz,yyz;
    momFloat xxxx,xxxy,xxyy,xyyy,yyyy,xxxz,xxyz,xyyz,yyyz;
    momFloat xxxxx,xxxxy,xxxyy,xxyyy,xyyyy,yyyyy,xxxxz,xxxyz,xxyyz,xyyyz,yyyyz;
    } LOCR;

/*
** The next set of data structures are intended specifically for use with float
** precision. These moments are usually scaled to a characteristic size of the 
** cell or volume. The convention is to use the scaling factor u for the multipole
** moments and scaling factor v for the local expansion.
*/
typedef struct fmomReduced {
    cosmoType m;
    cosmoType xx, yy, xy, xz, yz;
    cosmoType xxx, xyy, xxy, yyy, xxz, yyz, xyz;
    cosmoType xxxx, xyyy, xxxy, yyyy, xxxz, yyyz, xxyy, xxyz, xyyz;
    } FMOMR;

typedef struct flocReduced {
    cosmoType m;
    cosmoType x, y, z;
    cosmoType xx, yy, xy, xz, yz;
    cosmoType xxx, xyy, xxy, yyy, xxz, yyz, xyz;
    cosmoType xxxx, xyyy, xxxy, yyyy, xxxz, yyyz, xxyy, xxyz, xyyz;
    cosmoType xxxxx, xyyyy, xxxxy, yyyyy, xxxxz, yyyyz, xxxyy, xxyyy, xxxyz,
        xyyyz, xxyyz;
    } FLOCR;

void momClearMomr(MOMR *mr);
void momClearFmomr(FMOMR *l);
void momAddMomc(MOMC *,MOMC *);
void momAddMomr(MOMR *,MOMR *);
void momAddFmomr(FMOMR *mr,FMOMR *ma);
void momScaledAddFmomr(FMOMR *mr, cosmoType ur, FMOMR *ma, cosmoType ua);
void momRescaleFmomr(FMOMR *mr, cosmoType unew, cosmoType uold);
void momMulAddMomc(MOMC *,momFloat,MOMC *);
void momMulAddMomr(MOMR *,momFloat,MOMR *);
void momMulAddFmomr(FMOMR *mr, cosmoType ur, cosmoType m, FMOMR *ma,
                    cosmoType ua);
void momSubMomc(MOMC *,MOMC *);
void momSubMomr(MOMR *,MOMR *);
void momScaledSubFmomr(FMOMR *mr, cosmoType ur, FMOMR *ma, cosmoType ua);
void momMakeMomc(MOMC *,momFloat,momFloat,momFloat,momFloat);
cosmoType momMakeFmomr(FMOMR *mr, cosmoType m, cosmoType u, cosmoType x,
                       cosmoType y, cosmoType z);
momFloat momMakeMomr(MOMR *,momFloat,momFloat,momFloat,momFloat);
void momOldMakeMomr(MOMR *,momFloat,momFloat,momFloat,momFloat);
void momShiftMomc(MOMC *,momFloat,momFloat,momFloat);
void momShiftMomr(MOMR *,momFloat,momFloat,momFloat);
void momShiftFmomr(FMOMR *m, cosmoType u, cosmoType x, cosmoType y,
                   cosmoType z);
double momShiftLocr(LOCR *,momFloat,momFloat,momFloat);
void momReduceMomc(MOMC *,MOMR *);
void momEvalMomr(MOMR *,momFloat,momFloat,momFloat,momFloat,
				 momFloat *,momFloat *,momFloat *,momFloat *);
void momEvalFmomrcm(FMOMR *m, cosmoType u, cosmoType dir, cosmoType x,
                    cosmoType y, cosmoType z, cosmoType *fPot, cosmoType *ax,
                    cosmoType *ay, cosmoType *az, cosmoType *magai);
void momMomr2Momc(MOMR *,MOMC *);
void momFmomr2Momc(FMOMR *ma,MOMC *mc);
void momPrintMomc(MOMC *);
void momPrintMomr(MOMR *);

void momClearLocr(LOCR *);
double momLocrAddMomr5(LOCR *,MOMR *,momFloat,momFloat,momFloat,momFloat,double *,double *,double *);
double momFlocrAddFmomr5cm(FLOCR *l, cosmoType v, FMOMR *m, cosmoType u,
                           cosmoType dir, cosmoType x, cosmoType y, cosmoType z,
                           cosmoType *tax, cosmoType *tay, cosmoType *taz);
void momEvalLocr(LOCR *,momFloat,momFloat,momFloat,
		 momFloat *,momFloat *,momFloat *,momFloat *);
double momLocrAddMomr(LOCR *,MOMR *,momFloat,momFloat,momFloat,momFloat);
#if defined(__cplusplus)
}
#endif
#endif
