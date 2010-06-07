#ifndef PARAM_HINCLUDED
#define PARAM_HINCLUDED

/* Header for parameter parsing module.
 * First implemented by Joachim Stadel in PKDGRAV.
 */
#ifdef __cplusplus
extern "C" {
#endif
enum TypeParam {
    paramBool = 0,
    paramInt = 1,
    paramDouble = 2,
    paramString = 3
};

typedef struct prmNode {
	struct prmNode *pnNext;
	char *pszName;
	int iType;
	int bArg;
	int bFile;
	int iSize;
	void *pValue;
	char *pszArg;
	char *pszArgUsage;
	} PRM_NODE;

typedef struct prmContext {
	PRM_NODE *pnHead;
	void (*fcnLeader)(void);
	void (*fcnTrailer)(void);
	} * PRM;

#define PRM_LINE_SIZE	128

void prmInitialize(PRM *,void (*)(void),void (*)(void));
void prmFinish(PRM);
void prmAddParam(PRM,const char *,int,void *,int,const char *,const char *);
void prmArgUsage(PRM prm);
void prmLogParam(PRM prm, char *pszFile);
int prmParseParam(PRM,char *);
int prmArgProc(PRM,int,char **);
int prmSpecified(PRM,const char *);
int prmArgSpecified(PRM,const char *);
int prmFileSpecified(PRM,const char *);
#ifdef __cplusplus
}
#endif

#endif







