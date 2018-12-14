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

///
/// @brief Main method to perform Self Interacting Dark Matter interactions
/// @param dTime current simulation time
/// @param dDelta timestep over which to calculate interactions
/// @param activeRung timestep rung corresponding to dDelta
///
void Main::doSIDM(double dTime,double dDelta, int activeRung) { 
    if (param.iSIDMSelect!=0) {
        double stime = CkWallTimer();
        if(verbosity > 0) CkPrintf("SIDM interactions ... ");
        SIDMSmoothParams pSIDM(TYPE_DARK, activeRung, param.csm, dTime,param.dSIDMSigma, param.dSIDMVariable,param.iSIDMSelect,param.dDelta );  
        treeProxy.startSmooth(&pSIDM, 1, param.nSmooth, 0, CkCallbackResumeThread());         
        if(verbosity > 0) CkPrintf("took %g seconds\n", CkWallTimer() - stime);
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
    double sigma_classic, beta, term;  //velocity depedence terms
    double Sigma = 0.0;

    aDot=a*H;
    ih2 = invH2(p); //smoothing length, 1/h^2
    fNorm = M_1_PI*ih2*sqrt(ih2); 

    for (i=0;i<nSmooth;++i)
        {
        q = nnList[i].p;
        if (q->iOrder != p->iOrder)  //don't interact with self
        {
        ran=rand()/((double)RAND_MAX); //random number on [0,1]
        dvx = (-p->velocity.x + q->velocity.x)/a - aDot*nnList[i].dx.x; //v or vpred? 
        dvy = (-p->velocity.y + q->velocity.y)/a - aDot*nnList[i].dx.y; 
        dvz = (-p->velocity.z + q->velocity.z)/a - aDot*nnList[i].dx.z; 
        dvcosmo = sqrt(dvx*dvx + dvy*dvy + dvz*dvz); //relative velocity between particles
	r2 = nnList[i].fKey*ih2;

        if (iSIDMSelect==1) 
            {
            Sigma=dSIDMSigma;
            }

        else if (iSIDMSelect==2) //classical cross section dSIDMVariable is break in km/s
            {
            beta=M_PI*dSIDMVariable*dSIDMVariable/(dvcosmo*dvcosmo);
            if (beta < .1){
              sigma_classic=(4.0*M_PI/ 22.7 ) * beta *beta* log(1.0+(1.0/beta));
              }
            if ((beta >=.1 ) && (beta<1000)) {
              sigma_classic=(8.0*M_PI/ 22.7 ) * beta *beta /( 1.0+1.5* pow(beta,1.65) );
              }
            if (beta >= 1000.0 ){
              term= log(beta)+ 1.0 - .5*(1.0/log(beta));
              sigma_classic=(M_PI/ 22.7 ) *term *term;
              }
            Sigma=sigma_classic*dSIDMSigma; //convert to simulation units
            }

        else if (iSIDMSelect==3) //resonant cross section, dSIDMVariable is exponent value
	    {
            Sigma=dSIDMSigma*pow(dvcosmo,dSIDMVariable);
	    }
        else
            CkAbort("SIDM done with unknown cross section type (iSIDMSelect");

        //density in physical units? fNorm*KERNEL( r2 )*(p->mass)/(a*a*a)
        probability=Sigma*dvcosmo*dDelta*fNorm*KERNEL( r2, nSmooth )*(q->mass)/(a*a*a*2.0);
    
        //if ( ran>0.999999) {
        //   CkPrintf("SIDM Diagnostics: %g \n",probability);
        //   CkPrintf("iOrder: %i, dDelta: %g,  dvcosmo: %g, dSIDMSigma: %g \n",p->iOrder,dDelta,dvcosmo, dSIDMSigma);
        //   CkPrintf("fnorm: %g, kern: %g, r2: %g, q->mass: %g, aaa: %g  \n",fNorm,KERNEL( r2 ),r2,(q->mass),(a*a*a*2.0));
        //   CkPrintf("probability: %g, dvcosmo: %g, sigma: %g \n", probability,dvcosmo,Sigma);
        //   CkPrintf("iSIDMSelect: %i \n",(int)iSIDMSelect);
	//   }
 
        if (probability > .1 && ran>0.99999) {
           CkPrintf("SIDM Warning! The probability is rather large: %g \n",probability);
           CkPrintf("iOrder: %i, dDelta: %g,  dvcosmo: %g, Sigma: %g  \n",p->iOrder,dDelta,dvcosmo, Sigma);
           CkPrintf("fnorm: %g, kern: %g, r2: %g, q->mass: %g, aaa: %g  \n",fNorm,KERNEL( r2, nSmooth ),r2,(q->mass),(a*a*a*2.0));
           }

        if (probability > ran) {
#ifdef SIDMINTERACT
            p->iNSIDMInteractions += 1;
#endif
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
    }
