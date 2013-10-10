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

//see dendvdsmoothparams

void Main::doSIDM(double dTime,double dDelta, int activeRung)
{
//construct bSIDM, dSIDMsigma
if param.bSIDM {
	SIDMSmoothParams pSIDM(TYPE_DARK, activeRung, param.csm, dTime,dSIDMSigma, dDelta ); //things needs to be in constructor?? 
        treeProxy.startSmooth(&pSIDM, 1, param.Nsmooth, 0, CkCallbackResumeThread());         
	}
}

int SIDMSmoothParams::isSmoothActive(GravityParticle *p) //is a particle active is it the right type?
{
    if(bActiveOnly && p->rung < activeRung)
        return 0;               // not active
    return (TYPETest(p, iType));
    }


//what about  postTreeParticle(GravityParticle *p) {}
//and  initSmoothParticle(GravityParticle *p);

//void SIDMSmoothParams::initSmoothParticle(GravityParticle *p)
//{
//    TYPEReset(p, TYPE_NbrOfACTIVE);
//    }

//void SIDMSmoothParams::initTreeParticle(GravityParticle *p)
//{
//   TYPEReset(p, TYPE_NbrOfACTIVE);
//   }

void SIDMSmoothParams::initSmoothCache(GravityParticle *p)
{
//        p->nSIDMCollisionCount = 0;
//        ((PARTICLE *)p)->curlv[0] = 0;
//        ((PARTICLE *)p)->curlv[1] = 0;
//        ((PARTICLE *)p)->curlv[2] = 0;
        }

void SIDMSmoothParams::combSmoothCache(GravityParticle *p1, ExternalSmoothParticle *p2)
{
    //p1->nSIDMCollisionCount += p2->nSIDMCollisionCount;
    //((PARTICLE *)p1)->v[0] += ((PARTICLE *)p2)->curlv[0];//we are hijacking the curlv variable and accumulating q->v to send back to home node
    //((PARTICLE *)p1)->v[1] += ((PARTICLE *)p2)->curlv[1];
    //((PARTICLE *)p1)->v[2] += ((PARTICLE *)p2)->curlv[2];
    //p1->iType |= p2->iType;
    }


void SIDMSmoothParams::fcnSmooth(GravityParticle *p, int nSmooth, pqSmoothNode *nnList)
    {
    GravityParticle *q;
    double fDensity;
    double fNorm,ih2,r2,rs;
    double ran, probability,v,fPhyDensity, aDot, dvx,dvy,dvz,dvcosmo;
    int i;
    double mu, pxcm,pycm,pzcm,ecm0,ux,uy,uz;
    double norm,vxcm,vycm,vzcm, m1,m2;
    //double ftempe; //only used in print diagnostics
    double uvar,vvar;
    
    //    vFac = 1./(a*a); /* converts v to xdot */
    aDot=a*H;
    ih2 = invH2(p); //grabs smoothing length, 1/h^2
    fNorm = M_1_PI*ih2*sqrt(ih2); //
    fDensity = 0.0;

    for (i=0;i<nSmooth;++i)
        {
        r2 = nnList[i].fKey*ih2; 
        rs=KERNEL(r2);
        q = nnList[i].p;
        fDensity += rs*q->mass;
        }

        fPhyDensity=fNorm*fDensity/(a*a*a); /*in cosmolgical units */

    for (i=0;i<nSmooth;++i)
        {
        q = nnList[i].p;
        ran=rand()/((FLOAT)RAND_MAX); /* generate random number between 0 and 1 */

        dvx = (-p->velocity.x + q->velocity.x/a - aDot*nnList[i].dx.x; //could reduce these to one statment using vector arithmetic
        dvy = (-p->velocity.y + q->velocity.y/a - aDot*nnList[i].dx.y; //v or vpred?
        dvz = (-p->velocity.z + q->velocity.z/a - aDot*nnList[i].dx.z;  //units? H? or adot?
        dvcosmo = sqrt(dvx*dvx + dvy*dvy + dvz*dvz); /*relative velocity between particles */
        probability=dSIDMSigma*dvcosmo*dDelta*fPhyDensity/2.0;
                /* units: (cm^2 g^-1) (cm s^-1) (s) (g cm^-3) =  none, however in reality these are all sim units...
 *                 mass sim unit = real unit in g * conversion from g to Msol / dMsolUnit
 *                                 distance sim unit = real unit in cm * conversion from cm to kpc / dKpcUnit */              
        if (probability > ran) {
            //p->nSIDMCollisionCount+=1;
            //q->nSIDMCollisionCount+=1;
            m1=p->mass; //is this physical mass?
            m2=q->mass;
            mu=m1*m2/(m1+m2);
            pxcm = -mu*dvx; //=mu*(vx1-vx2) in frame where vx1 is zero
            pycm = -mu*dvy;
            pzcm = -mu*dvz;

            vxcm=-pxcm/m1; /*where vx1 is zero in this frame */
            vycm=-pycm/m1;
            vzcm=-pzcm/m1;

            ecm0=sqrt(pxcm*pxcm+pycm*pycm+pzcm*pzcm);//for printing out ecm0 is usefull this way
            norm=666;
                        while (norm>1.0){ //unit sphere point picking (Marsaglia 1972)
                            uvar=2.0*(rand()/(double)RAND_MAX)-1.0;  //#random number on [-1,1]
                            vvar=2.0*(rand()/(double)RAND_MAX)-1.0;
                            norm=(uvar*uvar+vvar*vvar);
                            }
            ux=2.0*uvar*sqrt(1-uvar*uvar-vvar*vvar);
            uy=2.0*vvar*sqrt(1-uvar*uvar-vvar*vvar);
            uz=1.0-2.0*(uvar*uvar+vvar*vvar);

            pxcm=ecm0*ux;
            pycm=ecm0*uy;
            pzcm=ecm0*uz;

            /* vphy=a*xdot + xdot*a  vsim=a^2 * xdot */
            p->velocity.x+=a*(pxcm/m1+vxcm);
            p->velocity.y+=a*(pycm/m1+vycm);
            p->velocity.z+=a*(pzcm/m1+vzcm);
            q->velocity.x+=a*(-pxcm/m2+vxcm-dvx);
            q->velocity.y+=a*(-pycm/m2+vycm-dvy);
            q->velocity.z+=a*(-pzcm/m2+vzcm-dvz);

            //q->curlv[0]+=aFac*(-pxcm/m2+vxcm-dvx); //this accumulation is for sending q vel update back to home node
            //q->curlv[1]+=aFac*(-pycm/m2+vycm-dvy);
            //q->curlv[2]+=aFac*(-pzcm/m2+vzcm-dvz);
            }
        }    
    }
