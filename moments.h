#ifndef MOMENTS_INCLUDED
#define MOMENTS_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

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
    double xxxxx, xxxxy, xxxyy, xxyyy, xyyyy, yyyyy, xxxxz, xxxyz, xxyyz, xyyyz, yyyyz;
} LOCR;

/*
** The next set of data structures are intended specifically for use with float
** precision. These moments are usually scaled to a characteristic size of the
** cell or volume. The convention is to use the scaling factor u for the multipole
** moments and scaling factor v for the local expansion.
*/
typedef struct fmomReduced {
    float m;
    float xx, yy, xy, xz, yz;
    float xxx, xyy, xxy, yyy, xxz, yyz, xyz;
    float xxxx, xyyy, xxxy, yyyy, xxxz, yyyz, xxyy, xxyz, xyyz;
} FMOMR;

typedef struct flocReduced {
    float m;
    float x, y, z;
    float xx, yy, xy, xz, yz;
    float xxx, xyy, xxy, yyy, xxz, yyz, xyz;
    float xxxx, xyyy, xxxy, yyyy, xxxz, yyyz, xxyy, xxyz, xyyz;
    float xxxxx, xyyyy, xxxxy, yyyyy, xxxxz, yyyyz, xxxyy, xxyyy, xxxyz, xyyyz, xxyyz;
} FLOCR;

void momClearMomr(MOMR *mr);
void momClearFmomr(FMOMR *l);
void momAddMomc(MOMC *, MOMC *);
void momAddMomr(MOMR *, MOMR *);
void momAddFmomr(FMOMR *mr, FMOMR *ma);
void momScaledAddFmomr(FMOMR *mr, float ur, FMOMR *ma, float ua);
void momRescaleFmomr(FMOMR *mr, float unew, float uold);
void momMulAddMomc(MOMC *, double, MOMC *);
void momMulAddMomr(MOMR *, double, MOMR *);
void momMulAddFmomr(FMOMR *mr, float ur, float m, FMOMR *ma, float ua);
void momSubMomc(MOMC *, MOMC *);
void momSubMomr(MOMR *, MOMR *);
void momScaledSubFmomr(FMOMR *mr, float ur, FMOMR *ma, float ua);
void momMakeMomc(MOMC *, double, double, double, double);
float momMakeFmomr(FMOMR *mr, float m, float u, float x, float y, float z);
double momMakeMomr(MOMR *, double, double, double, double);
void momOldMakeMomr(MOMR *, double, double, double, double);
void momShiftMomc(MOMC *, double, double, double);
void momShiftMomr(MOMR *, double, double, double);
void momShiftFmomr(FMOMR *m, float u, float x, float y, float z);
double momShiftLocr(LOCR *, double, double, double);
void momReduceMomc(MOMC *, MOMR *);
void momEvalMomr(MOMR *, double, double, double, double, double *, double *, double *, double *);
void momEvalFmomrcm(FMOMR *m, float u, float dir, float x, float y, float z, float *fPot, float *ax, float *ay,
                    float *az, float *magai);
void momMomr2Momc(MOMR *, MOMC *);
void momFmomr2Momc(FMOMR *ma, MOMC *mc);
void momPrintMomc(MOMC *);
void momPrintMomr(MOMR *);

void momClearLocr(LOCR *);
double momLocrAddMomr5(LOCR *, MOMR *, double, double, double, double, double *, double *, double *);
double momFlocrAddFmomr5cm(FLOCR *l, float v, FMOMR *m, float u, float dir, float x, float y, float z, float *tax,
                           float *tay, float *taz);
void momEvalLocr(LOCR *, double, double, double, double *, double *, double *, double *);
double momLocrAddMomr(LOCR *, MOMR *, double, double, double, double);

#if defined(__cplusplus)
}
#endif

#endif
