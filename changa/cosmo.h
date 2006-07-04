#ifndef COSMO_HINCLUDED
#define COSMO_HINCLUDED

#if defined(__cplusplus)
extern "C" {
#endif
/* Cosmo routines originally written for PKDGRAV by Thomas Quinn */

typedef struct csmContext {
    int bComove;	   
    double dHubble0;
    double dOmega0;
    double dLambda;
    double dOmegaRad;
    double dQuintess; /* w = -1/2 equation of  state */
    double dOmegab;
    } * CSM;

void csmInitialize(CSM *pcsm);
double csmExp2Hub(CSM csm, double dExp);
double csmTime2Hub(CSM csm,double dTime);
double csmExp2Time(CSM csm,double dExp);
double csmTime2Exp(CSM csm,double dTime);
double csmComoveDriftInt(CSM csm, double dIExp);
double csmComoveKickInt(CSM csm, double dIExp);
double csmComoveDriftFac(CSM csm,double dTime,double dDelta);
double csmComoveKickFac(CSM csm,double dTime,double dDelta);
double csmComoveLookbackTime2Exp(CSM csm,double dComoveTime);
double csmGrowthFac(CSM csm, double dExp);
double csmGrowthFacDot(CSM csm, double dExp);
double csmExp2Om(CSM csm, double dExp);
#if defined(__cplusplus)
}
#endif
     
#endif
