#ifndef MOMENTS_INCLUDED
#define MOMENTS_INCLUDED

#if defined(__cplusplus)
extern "C" {
#endif

#ifdef QUAD
typedef long double momFloat;
#define sqrt(x)	sqrtl(x)
#else
typedef double momFloat;
#endif

/*
 ** moment tensor components for reduced multipoles.
 */
typedef struct momReduced {
	momFloat m;
	momFloat xx,yy,xy,xz,yz;
	momFloat xxx,xyy,xxy,yyy,xxz,yyz,xyz;
	momFloat xxxx,xyyy,xxxy,yyyy,xxxz,yyyz,xxyy,xxyz,xyyz;
	} MOMR;

/*
 ** moment tensor components for complete multipoles.
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

/*
 ** moment tensor components for reduced local expansion.
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
    float m;
    float xx,yy,xy,xz,yz;
    float xxx,xyy,xxy,yyy,xxz,yyz,xyz;
    float xxxx,xyyy,xxxy,yyyy,xxxz,yyyz,xxyy,xxyz,xyyz;
    } FMOMR;

typedef struct flocReduced {
    float m;
    float x,y,z;
    float xx,yy,xy,xz,yz;
    float xxx,xyy,xxy,yyy,xxz,yyz,xyz;
    float xxxx,xyyy,xxxy,yyyy,xxxz,yyyz,xxyy,xxyz,xyyz;
    float xxxxx,xyyyy,xxxxy,yyyyy,xxxxz,yyyyz,xxxyy,xxyyy,xxxyz,xyyyz,xxyyz;
    } FLOCR;

void momClearMomr(MOMR *mr);
void momClearFmomr(FMOMR *l);
void momAddMomc(MOMC *,MOMC *);
void momAddMomr(MOMR *,MOMR *);
void momAddFmomr(FMOMR *mr,FMOMR *ma);
void momScaledAddFmomr(FMOMR *mr,float ur,FMOMR *ma,float ua);
void momRescaleFmomr(FMOMR *mr,float unew,float uold);
void momMulAddMomc(MOMC *,momFloat,MOMC *);
void momMulAddMomr(MOMR *,momFloat,MOMR *);
void momMulAddFmomr(FMOMR *mr,float ur,float m,FMOMR *ma,float ua);
void momSubMomc(MOMC *,MOMC *);
void momSubMomr(MOMR *,MOMR *);
void momScaledSubFmomr(FMOMR *mr,float ur,FMOMR *ma,float ua);
void momMakeMomc(MOMC *,momFloat,momFloat,momFloat,momFloat);
float momMakeFmomr(FMOMR *mr,float m,float u,float x,float y,float z);
momFloat momMakeMomr(MOMR *,momFloat,momFloat,momFloat,momFloat);
void momOldMakeMomr(MOMR *,momFloat,momFloat,momFloat,momFloat);
void momShiftMomc(MOMC *,momFloat,momFloat,momFloat);
void momShiftMomr(MOMR *,momFloat,momFloat,momFloat);
void momShiftFmomr(FMOMR *m,float u,float x,float y,float z);
double momShiftLocr(LOCR *,momFloat,momFloat,momFloat);
void momReduceMomc(MOMC *,MOMR *);
void momEvalMomr(MOMR *,momFloat,momFloat,momFloat,momFloat,
				 momFloat *,momFloat *,momFloat *,momFloat *);
void momEvalFmomrcm(FMOMR *m,float u,float dir,float x,float y,float z,
		    float *fPot,float *ax,float *ay,float *az,float *magai);
void momMomr2Momc(MOMR *,MOMC *);
void momFmomr2Momc(FMOMR *ma,MOMC *mc);
void momPrintMomc(MOMC *);
void momPrintMomr(MOMR *);

void momClearLocr(LOCR *);
double momLocrAddMomr5(LOCR *,MOMR *,momFloat,momFloat,momFloat,momFloat,double *,double *,double *);
double momFlocrAddFmomr5cm(FLOCR *l,float v,FMOMR *m,float u,float dir,float x,float y,float z,float *tax,float *tay,float *taz);
void momEvalLocr(LOCR *,momFloat,momFloat,momFloat,
		 momFloat *,momFloat *,momFloat *,momFloat *);
double momLocrAddMomr(LOCR *,MOMR *,momFloat,momFloat,momFloat,momFloat);
#if defined(__cplusplus)
}
#endif
#endif
