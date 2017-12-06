/* 
 * Module to parse parameter files and command line arguments.
 * First implemented by Joachim Stadel in PKDGRAV
 */
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "param.h"

void prmInitialize(PRM *pprm,void (*fcnLeader)(void),void (*fcnTrailer)(void))
{
	PRM prm;

	prm = (PRM)malloc(sizeof(struct prmContext));
	assert(prm != NULL);
	*pprm = prm;
	prm->pnHead = NULL;
	prm->pnTail = NULL;
	prm->fcnLeader = fcnLeader;
	prm->fcnTrailer = fcnTrailer;
	}


void prmFinish(PRM prm)
{
	PRM_NODE *pn,*pnKill;

	pn = prm->pnHead;
	while (pn) {
		pnKill = pn;
		pn = pn->pnNext;
		free(pnKill->pszName);
		if (pnKill->pszArg) free(pnKill->pszArg);
		if (pnKill->pszArgUsage) free(pnKill->pszArgUsage);
		free(pnKill);
		}
	free(prm);
	}


void prmAddParam(PRM prm,const char *pszName,int iType,void *pValue,
				 int iSize,const char *pszArg,const char *pszArgUsage)
{
	PRM_NODE *pn;

	pn = (PRM_NODE *)malloc(sizeof(PRM_NODE));
	assert(pn != NULL);
	pn->pszName = (char *)malloc(strlen(pszName)+1);
	assert(pn->pszName != NULL);
	strcpy(pn->pszName,pszName);
	pn->iType = iType;
	pn->iSize = iSize;
	pn->bArg = 0;
	pn->bFile = 0;
	pn->pValue = pValue;
	if (pszArg) {
		pn->pszArg = (char *)malloc(strlen(pszArg)+1);
		assert(pn->pszArg != NULL);
		strcpy(pn->pszArg,pszArg);
		}
	else pn->pszArg = NULL;
	if (pszArgUsage) {
		pn->pszArgUsage = (char *)malloc(strlen(pszArgUsage)+1);
		assert(pn->pszArgUsage != NULL);
		strcpy(pn->pszArgUsage,pszArgUsage);
		}
	else pn->pszArgUsage = NULL;
	pn->pnNext = NULL;
	if (!prm->pnHead) {
	        prm->pnHead = pn;
		prm->pnTail = pn; 
	        }
	else {	
		prm->pnTail->pnNext = pn;
		prm->pnTail = pn; 
		}
	}


void prmArgUsage(PRM prm)
{
	PRM_NODE *pn;

	if (prm->fcnLeader) (*prm->fcnLeader)();
	pn = prm->pnHead;
	while (pn) {
		if (pn->pszArg && pn->pszArgUsage) {
			if (pn->iType == paramBool) {
				printf("[+%s][-%s] %s\n",pn->pszArg,pn->pszArg,
					   pn->pszArgUsage);
				}
			else {
				printf("[-%s %s]\n",pn->pszArg,pn->pszArgUsage);
				}
			}
		pn = pn->pnNext;
		}
	if (prm->fcnTrailer) (*prm->fcnTrailer)();
	}

void prmLogParam(PRM prm, const char *pszFile)
{
    FILE *fpLog;
    PRM_NODE *pn;
    
    fpLog = fopen(pszFile,"a");
    if(fpLog == NULL) {
	fprintf(stderr, "prmLogParam: Can't open file %s\n", pszFile);
	return;
	}
    
    fprintf(fpLog, "# Parameters:\n");
    pn = prm->pnHead;
    while (pn) {
	switch (pn->iType) {
	case paramBool:
	case paramInt:
	    fprintf(fpLog, "# %s: %d\n", pn->pszName, *((int *)pn->pValue));
	    break;
	case paramDouble:
	    fprintf(fpLog, "# %s: %g\n", pn->pszName, *((double *)pn->pValue));
	    break;
	case paramString:
	    fprintf(fpLog, "# %s: %s\n", pn->pszName, (char *)pn->pValue);
	    }
	pn = pn->pnNext;
	}
    fclose(fpLog);
    }

int prmParseParam(PRM prm,char *pszFile)
{
	FILE *fpParam;
	PRM_NODE *pn;
	char achBuf[PRM_LINE_SIZE];
	char *p,*q,*pszCmd,t;
	int iLine,ret;

	// mark the end of the buffer with a special character
	//  if it gets overwritten, we will know the buffer was not large enough
	achBuf[PRM_LINE_SIZE - 1] = '\a'; 

	fpParam = fopen(pszFile,"r");
	if (!fpParam) {
		printf("Could not open file:%s\n",pszFile);
		return(0);
		}
	p = fgets(achBuf,PRM_LINE_SIZE,fpParam);
	assert(achBuf[PRM_LINE_SIZE - 1] == '\a'); 
	iLine = 1;
	while (p) {
		if (*p == 0) goto new_line;
		if (*p == '#') goto new_line;
		while (isspace((int) *p)) {
			++p;
			if (*p == 0) goto new_line;
			}
		if (isalpha((int) *p)) {
			pszCmd = p;
			++p;
			if (*p == 0) goto lookup_cmd;
			}
		else goto syntax_error;
		while (isalnum((int) *p)||strchr("_$",*p)) {
			++p;
			if (*p == 0) goto lookup_cmd;
			}
	lookup_cmd:
		t = *p;
		*p = 0;
		pn = prm->pnHead;
		while (pn) {
			if (!strcmp(pszCmd,pn->pszName)) break;
			pn = pn->pnNext;
			}
		if (!pn) goto cmd_error;
		*p = t;
		if (*p == 0) goto syntax_error;
		while (isspace((int) *p)) {
			++p;
			if (*p == 0) goto syntax_error;
			}
		if (*p != '=') goto syntax_error;
		++p;
		if (*p == 0) goto syntax_error;
		while (isspace((int) *p)) {
			++p;
			if (*p == 0) goto syntax_error;
			}
		switch (pn->iType) {
		case paramBool:
			assert(pn->iSize == sizeof(int));
			ret = sscanf(p,"%d",(int *)pn->pValue);
			if (ret != 1) goto syntax_error;
			break;
		case paramInt:
			assert(pn->iSize == sizeof(int));
			ret = sscanf(p,"%d",(int *)pn->pValue);
			if (ret != 1) goto syntax_error;
			break;
		case paramDouble:
			assert(pn->iSize == sizeof(double));
			ret = sscanf(p,"%lf",(double *)pn->pValue);
			if (ret != 1) goto syntax_error;
			break;
 		case paramString:
			/*
			 ** Make sure there is enough space to handle the string.
			 ** This is a CONSERVATIVE test.
			 */
			if(pn->iSize <= strlen(p))
			    printf("Line size problem in %s line %d\n",
				   pszFile, iLine);
			assert(pn->iSize > strlen(p));
			ret = sscanf(p,"%[^\n#]",(char *)pn->pValue);
			if (ret != 1) goto syntax_error;
			/*
			 ** Strip trailing whitespace. OKAY!
			 */
			p = pn->pValue;
			q = &p[strlen(p)];
			while (--q >= p) if (!isspace((int) *q)) break;
			++q;
			*q = 0;
			break;
		default:
			goto cmd_error;
			}
		pn->bFile = 1;
	new_line:
		p = fgets(achBuf,PRM_LINE_SIZE,fpParam);
		assert(achBuf[PRM_LINE_SIZE - 1] == '\a'); 
		++iLine;
		}
	fclose(fpParam);
	return(1);
 syntax_error:
	q = achBuf;
	while (*q) {
		if (*q == '\n') *q = 0;
		else ++q;
		}
	printf("Syntax error in %s(%d):\n%s",pszFile,iLine,achBuf);
	fclose(fpParam);
	return(0);
 cmd_error:
	q = achBuf;
	while (*q) {
		if (*q == '\n') *q = 0;
		else ++q;
		}
	printf("Unrecognized command in %s(%d):%s\n",pszFile,iLine,achBuf);
	fclose(fpParam);
	return(0);
	}


int prmArgProc(PRM prm,int argc,char **argv,int processSimfile)
{
	int i,ret;
	PRM_NODE *pn;

    if (argc < 2) return(1);
    /*
     ** Look for the sim_file name and process it first if requested. If no
     ** sim_file name then ignore and look at command line settings
     ** and options. The cases are a bit tricky.
     */
	if (*argv[argc-1] != '-' && *argv[argc-1] != '+') {
		if (*argv[argc-2] != '-') {
		        if (processSimfile) {
			       if (!prmParseParam(prm,argv[argc-1])) return(0);
			       }
			--argc;
			}
		else {
			pn = prm->pnHead;
			while (pn) {
				if (pn->pszArg) 
					if (!strcmp(&argv[argc-2][1],pn->pszArg)) break;
				pn = pn->pnNext;
				}
			if (pn) {
				if (pn->iType == paramBool) {
					/*
					 ** It's a boolean flag.
					 */
				        if (processSimfile) {
			                        if (!prmParseParam(prm,argv[argc-1])) return(0);
			                        }
					--argc;
					}
				}
			else {
				if (processSimfile) {
			                if (!prmParseParam(prm,argv[argc-1])) return(0);
			                }
				--argc;
				}
			}
		}
	for (i=1;i<argc;++i) {
		if (*argv[i] == '-' || *argv[i] == '+') {
			pn = prm->pnHead;
			while (pn) {
				if (pn->pszArg)
					if (!strcmp(&argv[i][1],pn->pszArg)) break;
				pn = pn->pnNext;
				}
			if (!pn) {
				printf("Unrecognized command line argument:%s\n",argv[i]);
				prmArgUsage(prm);
				return(0);
				}
			}
		else if (i == argc-1) return(1);
		else {
			printf("Unrecognized command line argument:%s\n",argv[i]);
			prmArgUsage(prm);
			return(0);
			}
		switch (pn->iType) {
		case paramBool:
			/*
			 ** It's a boolean.
			 */
			if (argv[i][0] == '-') *((int *)pn->pValue) = 0;
			else *((int *)pn->pValue) = 1;
			break;
		case paramInt:
			/*
			 ** It's an int
			 */
			++i;
			if (i == argc) {
				printf("Missing integer value after command line ");
			    printf("argument:%s\n",argv[i-1]);
				prmArgUsage(prm);
				return(0);
				}
			ret = sscanf(argv[i],"%d",(int *) pn->pValue);
			if (ret != 1) {
				printf("Expected integer after command line ");
			    printf("argument:%s\n",argv[i-1]);
				prmArgUsage(prm);
				return(0);
				}
			break;
		case paramDouble:
			/*
			 ** It's a DOUBLE
			 */
			++i;
			if (i == argc) {
				printf("Missing double value after command line ");
			    printf("argument:%s\n",argv[i-1]);
				prmArgUsage(prm);
				return(0);
				}
			ret = sscanf(argv[i],"%lf",(double *)pn->pValue);
			if (ret != 1) {
				printf("Expected double after command line ");
			    printf("argument:%s\n",argv[i-1]);
				prmArgUsage(prm);
				return(0);
				}
			break;
		case paramString:
			/*
			 ** It's a string
			 */
			++i;
			if (i == argc) {
				printf("Missing string after command line ");
			    printf("argument:%s\n",argv[i-1]);
				prmArgUsage(prm);
				return(0);
				}
			assert(pn->iSize > strlen(argv[i]));
			strcpy((char *)pn->pValue,argv[i]);
			break;
		default:
			assert(0);
			}
		pn->bArg = 1;
		}
	return(1);
	}

int prmArgSpecified(PRM prm,const char *pszName)
{
	PRM_NODE *pn;
	
	pn = prm->pnHead;	
	while (pn) {
		if (pn->pszArg)
			if (!strcmp(pn->pszName,pszName)) break;
		pn = pn->pnNext;
		}
	if (!pn) return(0);
	return(pn->bArg);
	}


int prmFileSpecified(PRM prm,const char *pszName)
{
	PRM_NODE *pn;
	
	pn = prm->pnHead;	
	while (pn) {
		if (!strcmp(pn->pszName,pszName)) break;
		pn = pn->pnNext;
		}
	if (!pn) return(0);
	return(pn->bFile);
	}


int prmSpecified(PRM prm,const char *pszName)
{
	return(prmArgSpecified(prm,pszName) || prmFileSpecified(prm,pszName));
	}
