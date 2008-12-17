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

void momClearMomr(MOMR *mr);
void momAddMomc(MOMC *,MOMC *);
void momAddMomr(MOMR *,MOMR *);
void momMulAddMomc(MOMC *,momFloat,MOMC *);
void momMulAddMomr(MOMR *,momFloat,MOMR *);
void momSubMomc(MOMC *,MOMC *);
void momSubMomr(MOMR *,MOMR *);
void momMakeMomc(MOMC *,momFloat,momFloat,momFloat,momFloat);
momFloat momMakeMomr(MOMR *,momFloat,momFloat,momFloat,momFloat);
void momOldMakeMomr(MOMR *,momFloat,momFloat,momFloat,momFloat);
void momShiftMomc(MOMC *,momFloat,momFloat,momFloat);
void momShiftMomr(MOMR *,momFloat,momFloat,momFloat);
double momShiftLocr(LOCR *,momFloat,momFloat,momFloat);
void momReduceMomc(MOMC *,MOMR *);
void momEvalMomr(MOMR *,momFloat,momFloat,momFloat,momFloat,
				 momFloat *,momFloat *,momFloat *,momFloat *);
void momMomr2Momc(MOMR *,MOMC *);
void momPrintMomc(MOMC *);
void momPrintMomr(MOMR *);

void momClearLocr(LOCR *);
double momLocrAddMomr5(LOCR *,MOMR *,momFloat,momFloat,momFloat,momFloat,double *,double *,double *);
void momEvalLocr(LOCR *,momFloat,momFloat,momFloat,
		 momFloat *,momFloat *,momFloat *,momFloat *);
double momLocrAddMomr(LOCR *,MOMR *,momFloat,momFloat,momFloat,momFloat);
#if defined(__cplusplus)
}
#endif
#endif
