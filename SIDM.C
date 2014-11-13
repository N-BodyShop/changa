//SIDM- self interacting dark matter by Alexander Fry
#include <math.h>
#ifdef HAVE_VALUES_H
#include <values.h>
#else
#include <float.h>
#endif
#include "ParallelGravity.h"
#include "smooth.h"
#include "SIDM.h"


void Main::doSIDM(double dTime,double dDelta, int activeRung) { 
    if (param.bSIDM) {
        SIDMSmoothParams pSIDM(TYPE_DARK, activeRung, param.csm, dTime,param.dSIDMSigma, param.dDelta );  
        treeProxy.startSmooth(&pSIDM, 1, param.nSmooth, 0, CkCallbackResumeThread());         
        }
    }

int SIDMSmoothParams::isSmoothActive(GravityParticle *p) {
    return (TYPETest(p, iType));
    }


void SIDMSmoothParams::initSmoothParticle(GravityParticle *p) {
    }

void SIDMSmoothParams::initTreeParticle(GravityParticle *p) {
   }

void SIDMSmoothParams::initSmoothCache(GravityParticle *p) {
    p->treeAcceleration = p->velocity;
    }

void SIDMSmoothParams::combSmoothCache(GravityParticle *p1, ExternalSmoothParticle *p2) {
    Vector3D<double> deltav; //vector template notation
    deltav=p2->velocity - p2->treeAcceleration;
    p1->velocity += deltav;
    }

void SIDMSmoothParams::fcnSmooth(GravityParticle *p, int nSmooth, pqSmoothNode *nnList) {
    GravityParticle *q;
    double fNorm,ih2,r2;
    double ran, probability,aDot, dvx,dvy,dvz,dvcosmo;
    int i;
    double mu, pxcm,pycm,pzcm,pcm0,ux,uy,uz;
    double norm,vxcm,vycm,vzcm, m1,m2;
    double uvar,vvar;

    aDot=a*H;
    ih2 = invH2(p); //smoothing length, 1/h^2
    fNorm = M_1_PI*ih2*sqrt(ih2); 

    for (i=0;i<nSmooth;++i)
        {
        q = nnList[i].p;
        ran=rand()/((double)RAND_MAX); //random number on [0,1]
        dvx = (-p->velocity.x + q->velocity.x)/a - aDot*nnList[i].dx.x; //v or vpred? 
        dvy = (-p->velocity.y + q->velocity.y)/a - aDot*nnList[i].dx.y; 
        dvz = (-p->velocity.z + q->velocity.z)/a - aDot*nnList[i].dx.z; 
        dvcosmo = sqrt(dvx*dvx + dvy*dvy + dvz*dvz); //relative velocity between particles
	r2 = nnList[i].fKey*ih2;
        //density in physical units is fNorm*KERNEL( r2 )*(p->mass)/(a*a*a)
        probability=dSIDMSigma*dvcosmo*dDelta*fNorm*KERNEL( r2 )*(q->mass)/(a*a*a*2.0);

        if (probability > .02) {
           CkPrintf("SIDM Warning! The probability is rather large: %g \n",probability);
           CkPrintf("iOrder: %i, dDelta: %g,  dvcosmo: %g, dSIDMSigma %d, \n",p->iOrder,dDelta,dvcosmo, dSIDMSigma);
           }

        if (probability > ran) {
            p->iNSIDMInteractions += 1;
            m1=p->mass; 
            m2=q->mass;
            mu=m1*m2/(m1+m2);
            pxcm = -mu*dvx; //momentum COM frame where vx1 is 0
            pycm = -mu*dvy; 
            pzcm = -mu*dvz; 
            vxcm=-pxcm/m1; //velocity in COM frame where vx1 0
            vycm=-pycm/m1;
            vzcm=-pzcm/m1;
            pcm0=sqrt(pxcm*pxcm+pycm*pycm+pzcm*pzcm);
            norm=666;
            while (norm>1.0){ //unit sphere point picking (Marsaglia 1972)
                uvar=2.0*(rand()/(double)RAND_MAX)-1.0;  //#random number on [-1,1]
                vvar=2.0*(rand()/(double)RAND_MAX)-1.0;
                norm=(uvar*uvar+vvar*vvar);
                }
            ux=2.0*uvar*sqrt(1.0-uvar*uvar-vvar*vvar);
            uy=2.0*vvar*sqrt(1.0-uvar*uvar-vvar*vvar);
            uz=1.0-2.0*(uvar*uvar+vvar*vvar);
            pxcm=pcm0*ux;
            pycm=pcm0*uy;
            pzcm=pcm0*uz;
            p->velocity.x+=a*(pxcm/m1); 
            p->velocity.y+=a*(pycm/m1);
            p->velocity.z+=a*(pzcm/m1);
            q->velocity.x+=a*(-pxcm/m2);
            q->velocity.y+=a*(-pycm/m2);
            q->velocity.z+=a*(-pzcm/m2);
            //to transform back to the lab frame we must add back in the COM velocity
            p->velocity.x+=a*(vxcm); 
            p->velocity.y+=a*(vycm);
            p->velocity.z+=a*(vzcm);
            q->velocity.x+=a*(vxcm-dvx);
            q->velocity.y+=a*(vycm-dvy);
            q->velocity.z+=a*(vzcm-dvz);
            }
        }    
    }
