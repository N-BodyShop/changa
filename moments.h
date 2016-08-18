#ifndef MOMENTS_INCLUDED
#define MOMENTS_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

#include "cosmoType.h"

#ifdef QUAD
typedef long double momFloat;
#define sqrt(x) sqrtl(x)
#else
typedef double momFloat;
#endif

/*************************************************
 *
 *  The precision of the non-scaled moments (MOMR,
 *  MOMC, and LOCR) are to always be double to
 *  prevent known over/underflow issues.
 *
 *              Do not change this
 *
 *************************************************/

/**
 ** @brief moment tensor components for reduced multipoles.
 */
typedef struct momReduced {
    double m;
    double xx, yy, xy, xz, yz;
    double xxx, xyy, xxy, yyy, xxz, yyz, xyz;
    double xxxx, xyyy, xxxy, yyyy, xxxz, yyyz, xxyy, xxyz, xyyz;
} MOMR;

/**
 ** @brief moment tensor components for complete multipoles.
 */
typedef struct momComplete {
    double m;
    double xx, yy, xy, xz, yz;
    double xxx, xyy, xxy, yyy, xxz, yyz, xyz;
    double xxxx, xyyy, xxxy, yyyy, xxxz, yyyz, xxyy, xxyz, xyyz;
    double zz;
    double xzz, yzz, zzz;
    double xxzz, xyzz, xzzz, yyzz, yzzz, zzzz;
} MOMC;

/**
 ** @brief moment tensor components for reduced local expansion.
 ** note that we have the 5th-order terms here now!
 */
typedef struct locReduced {
    double m;
    double x, y, z;
    double xx, xy, yy, xz, yz;
    double xxx, xxy, xyy, yyy, xxz, xyz, yyz;
    double xxxx, xxxy, xxyy, xyyy, yyyy, xxxz, xxyz, xyyz, yyyz;
    double xxxxx, xxxxy, xxxyy, xxyyy, xyyyy, yyyyy, xxxxz, xxxyz, xxyyz, xyyyz,
        yyyyz;
} LOCR;

/*
** These moments are usually scaled to a characteristic size of the
** cell or volume. The convention is to use the scaling factor u for the
** multipole moments and scaling factor v for the local expansion.
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
void momAddMomc(MOMC *, MOMC *);
void momAddMomr(MOMR *, MOMR *);
void momAddFmomr(FMOMR *mr, FMOMR *ma);
void momMulAddMomc(MOMC *, double, MOMC *);
void momMulAddMomr(MOMR *, double, MOMR *);
void momSubMomc(MOMC *, MOMC *);
void momSubMomr(MOMR *, MOMR *);
double momLocrAddMomr5(LOCR *, MOMR *, double, double, double, double, double *,
                       double *, double *);
void momEvalLocr(LOCR *, double, double, double, double *, double *, double *,
                 double *);
double momLocrAddMomr(LOCR *, MOMR *, double, double, double, double);
void momMakeMomc(MOMC *, double, double, double, double);
double momMakeMomr(MOMR *, double, double, double, double);
void momOldMakeMomr(MOMR *, double, double, double, double);
void momShiftMomc(MOMC *, double, double, double);
void momShiftMomr(MOMR *, double, double, double);
double momShiftLocr(LOCR *, double, double, double);
void momReduceMomc(MOMC *, MOMR *);
void momEvalMomr(MOMR *, double, double, double, double, double *, double *,
                 double *, double *);
void momMomr2Momc(MOMR *, MOMC *);
void momFmomr2Momc(FMOMR *ma, MOMC *mc);
void momPrintMomc(MOMC *);
void momPrintMomr(MOMR *);
void momClearLocr(LOCR *);

cosmoType momMakeFmomr(FMOMR *mr, cosmoType m, cosmoType u, cosmoType x,
                       cosmoType y, cosmoType z);
void momShiftFmomr(FMOMR *m, cosmoType u, cosmoType x, cosmoType y,
                   cosmoType z);
void momEvalFmomrcm(FMOMR *m, cosmoType u, cosmoType dir, cosmoType x,
                    cosmoType y, cosmoType z, cosmoType *fPot, cosmoType *ax,
                    cosmoType *ay, cosmoType *az, cosmoType *magai);
void momScaledAddFmomr(FMOMR *mr, cosmoType ur, FMOMR *ma, cosmoType ua);
void momRescaleFmomr(FMOMR *mr, cosmoType unew, cosmoType uold);
void momMulAddFmomr(FMOMR *mr, cosmoType ur, cosmoType m, FMOMR *ma,
                    cosmoType ua);
void momScaledSubFmomr(FMOMR *mr, cosmoType ur, FMOMR *ma, cosmoType ua);
double momFlocrAddFmomr5cm(FLOCR *l, cosmoType v, FMOMR *m, cosmoType u,
                           cosmoType dir, cosmoType x, cosmoType y, cosmoType z,
                           cosmoType *tax, cosmoType *tay, cosmoType *taz);

#if defined(__cplusplus)
}
#endif

#endif
