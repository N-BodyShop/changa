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
//#include "physconst.h"

//void SIDM::AddParams(PRM prm, Parameters &params) 
//{
//}

void Main::doSIDM(double dTime,double dDelta, int activeRung)
{
//construct bSIDM, dSIDMsigma
if (param.bSIDM) {
	SIDMSmoothParams pSIDM(TYPE_DARK, activeRung, param.csm, dTime,param.dSIDMSigma, param.dDelta ); //things needs to be in constructor?? 
        treeProxy.startSmooth(&pSIDM, 1, param.nSmooth, 0, CkCallbackResumeThread());         
	}
}

int SIDMSmoothParams::isSmoothActive(GravityParticle *p) //is a particle active is it the right type?
{
//    if(bActiveOnly && p->rung < activeRung)
//        return 0;               // not active
    return (TYPETest(p, iType));
    }


//what about  postTreeParticle(GravityParticle *p) {}
//and  initSmoothParticle(GravityParticle *p);

void SIDMSmoothParams::initSmoothParticle(GravityParticle *p)
{
//    TYPEReset(p, TYPE_NbrOfACTIVE);
    }

void SIDMSmoothParams::initTreeParticle(GravityParticle *p)
{
//   TYPEReset(p, TYPE_NbrOfACTIVE);
   }

void SIDMSmoothParams::initSmoothCache(GravityParticle *p)
{
    //p->iNSIDMInteractions = 0;
    p->treeAcceleration = p->velocity; // This only works because
    //((extraSPHData*)p->extraData)->vPred() = p->velocity; // This only works because
			  // we know it is a cached particle,
			  // and therefore has "extraSPHData" even though it is dark matter.
			  // If the code at CacheInterface.C:125 gets changed, then this
			  // needs to be revisited.
    }

void SIDMSmoothParams::combSmoothCache(GravityParticle *p1, ExternalSmoothParticle *p2)
{
    Vector3D<double> deltav; //confirm vector template notation
    deltav=p2->velocity - p2->treeAcceleration;
    p1->velocity += deltav;
    //p1->iNSIDMInteractions += p2-> iNSIDMInteractions;
    //p1->iType |= p2->iType;
    }


//void SIDMSmoothParams::doSIDM(GravityParticle *p, int nSmooth, pqSmoothNode *nnList)
void SIDMSmoothParams::fcnSmooth(GravityParticle *p, int nSmooth, pqSmoothNode *nnList)
    {
    GravityParticle *q;
    double fDensity;
    double fNorm,ih2,r2,rs;
    double ran, probability,v,fPhyDensity, aDot, dvx,dvy,dvz,dvcosmo;
    int i;
    double mu, pxcm,pycm,pzcm,pcm0,ux,uy,uz;
    double norm,vxcm,vycm,vzcm, m1,m2;
    double uvar,vvar;
    //double e1,e2,etot_init,etot_final;
    double ran2,pvxcm,pvycm,pvzcm; //delete these, testing
    double qvxcm,qvycm,qvzcm; //delete thes

    //dSIDMSigma=dSIDMSigma*param.dMsolUnit*MSOLG/(param.dKpcUnit*KPCCM*param.dKpcUnit*KPCCM);

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
        ran=rand()/((double)RAND_MAX); /* generate random number between 0 and 1 */
        ran2=rand()/((double)RAND_MAX); /* generate random number between 0 and 1, for random diagnostic print outs */

        dvx = (-p->velocity.x + q->velocity.x)/a - aDot*nnList[i].dx.x; //could reduce these to one statment using vector arithmetic
        dvy = (-p->velocity.y + q->velocity.y)/a - aDot*nnList[i].dx.y; //v or vpred?
        dvz = (-p->velocity.z + q->velocity.z)/a - aDot*nnList[i].dx.z;  //units? H? or adot?
        dvcosmo = sqrt(dvx*dvx + dvy*dvy + dvz*dvz); /*relative velocity between particles */
	r2 = nnList[i].fKey*ih2;

        probability=dSIDMSigma*dvcosmo*dDelta*fNorm*KERNEL( r2 )*(p->mass)/(a*a*a*2.0);


        /*print out some diagnostics */

        if (p->iOrder%2000==0 && i==0) {
                CkPrintf("<SIDM a diagnostic>, %i, %i, %g, %g, %g, %g, %g, %g, %g, %g \n",p->iOrder, q->iOrder, fPhyDensity, dvcosmo, probability, nnList[i].dx.x,nnList[i].dx.y, nnList[i].dx.z, a, aDot);
                CkPrintf("<SIDM b diagnostic>, %i, %i, %g, %g, %g \n",p->iOrder, q->iOrder, fPhyDensity, dDelta, fNorm*KERNEL( r2 )*(p->mass));
                CkPrintf("<SIDM c diagnostic>, %g, %g, %g \n",r2,(-p->velocity.x + q->velocity.x)/a , aDot*nnList[i].dx.x );
	}


        if (probability > ran) {
            p->iNSIDMInteractions += 1;

            if (q->iOrder== p->iOrder)//turn this into asssert
                {
                CkPrintf("<bad meta-self-interaction occured on equal iorder> i: %i,  probability: %g \n",i,probability);
                CkPrintf("<bad meta-self-interaction occured> i: %i,  dvcosmo: %g \n",i,dvcosmo);
                CkPrintf("<bad meta-self-interaction occured> p->iOrder: %i,  q->iOrder: %i \n", p->iOrder,  q->iOrder);
                CkPrintf("<bad meta-self-interaction occured> p->iNSIDMInteractions : %i, q->iNSIDMInteractions: %i \n", p->iNSIDMInteractions , q->iNSIDMInteractions );
                //continue;
                }
            m1=p->mass; 
            m2=q->mass;
            mu=m1*m2/(m1+m2);
            pxcm = -mu*dvx; //the momentum of the particles in the COM frame
            pycm = -mu*dvy; //p1cm=mu*(vx1-vx2) in frame where vx1 is zero
            pzcm = -mu*dvz; //and p1cm=-p2cm

            vxcm=-pxcm/m1; /* velocity in COM frame. where vx1 is zero in this frame */
            vycm=-pycm/m1;
            vzcm=-pzcm/m1;

            pcm0=sqrt(pxcm*pxcm+pycm*pycm+pzcm*pzcm);
  	    //in non-cosmo runs check e conservation directly
           //     e1=0.5*m1*((p->velocity.x)*(p->velocity.x) + (p->velocity.y)*(p->velocity.y) + (p->velocity.z)*(p->velocity.z));
           //     e2=0.5*m2*((q->velocity.x)*(q->velocity.x) + (q->velocity.y)*(q->velocity.y) + (q->velocity.z)*(q->velocity.z));
           //     etot_init=e1+e2;
            if (ran2 < 0.01) {
                CkPrintf("<SIDM> i: %i,  dvcosmo: %g \n",i,dvcosmo);
                CkPrintf("<SIDM> i: %i,  fPhyDensity: %g \n",i,fPhyDensity);
                CkPrintf("<SIDM> kernelr, r2: %g, %g \n",fNorm*KERNEL( r2 )*(p->mass), r2);
           //     e1=0.5*m1*((p->velocity.x)*(p->velocity.x) + (p->velocity.y)*(p->velocity.y) + (p->velocity.z)*(p->velocity.z));
           //     e2=0.5*m2*((q->velocity.x)*(q->velocity.x) + (q->velocity.y)*(q->velocity.y) + (q->velocity.z)*(q->velocity.z));
           //     etot_init=e1+e2;
           //     CkPrintf("<SIDM> (noncosmo) etot_init: %g,  e1: %g, e2: %g \n",etot_init,e1,e2);
                CkPrintf("<SIDM> dt: %g \n",dDelta);
                CkPrintf("<SIDM> vx: %g, vy: %g, vz: %g \n", p->velocity.x,p->velocity.y,p->velocity.z );
                CkPrintf("<SIDM> dvx: %g, dvy: %g, dvz: %g \n", dvx,dvy,dvz );
                CkPrintf("<SIDM> ran: %g,  probability: %g \n",ran, probability);
                CkPrintf("<SIDM> dSIDMSigma: %g,  dDelta: %g \n",dSIDMSigma,dDelta);
                CkPrintf("<SIDM> mass p: %g,  mass q: %g \n",p->mass,q->mass);
                CkPrintf("<SIDM> momentum initial: %g,  density: %g \n",pcm0, fPhyDensity);
                CkPrintf("<SIDM> px: %g,  py: %g, pz: %g \n",pxcm,pycm,pzcm);
                CkPrintf("<SIDM> ninteractions: %i \n",p->iNSIDMInteractions);
                CkPrintf("<SIDM> a: %g, adot: %g  \n",a, aDot);
                }
         
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

            p->velocity.x+=a*(pxcm/m1); //equivalent to vcm1finalx
            p->velocity.y+=a*(pycm/m1);
            p->velocity.z+=a*(pzcm/m1);
            q->velocity.x+=a*(-pxcm/m2);
            q->velocity.y+=a*(-pycm/m2);
            q->velocity.z+=a*(-pzcm/m2);
 
            p->velocity.x+=a*(vxcm); //to transform back to the lab frame we must add back in the COM velocity
            p->velocity.y+=a*(vycm);
            p->velocity.z+=a*(vzcm);
            q->velocity.x+=a*(vxcm-dvx);
            q->velocity.y+=a*(vycm-dvy);
            q->velocity.z+=a*(vzcm-dvz);


            //these are diagnostic print outs
            //e1=0.5*m1*((p->velocity.x)*(p->velocity.x) + (p->velocity.y)*(p->velocity.y) + (p->velocity.z)*(p->velocity.z));
            //e2=0.5*m2*((q->velocity.x)*(q->velocity.x) + (q->velocity.y)*(q->velocity.y) + (q->velocity.z)*(q->velocity.z));
            //etot_final=e1+e2;
            if (ran2 < 0.01) {
                dvx = (-p->velocity.x + q->velocity.x)/a - aDot*nnList[i].dx.x; 
                dvy = (-p->velocity.y + q->velocity.y)/a - aDot*nnList[i].dx.y;
                dvz = (-p->velocity.z + q->velocity.z)/a - aDot*nnList[i].dx.z;
                dvcosmo = sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
                pxcm = -mu*dvx;
                pycm = -mu*dvy;
                pzcm = -mu*dvz;
                pcm0=sqrt(pxcm*pxcm+pycm*pycm+pzcm*pzcm);
                CkPrintf("<SIDM> i: %i,  dvcosmo: %g \n",i,dvcosmo);
            //    e1=0.5*m1*((p->velocity.x)*(p->velocity.x) + (p->velocity.y)*(p->velocity.y) + (p->velocity.z)*(p->velocity.z));
            //    e2=0.5*m2*((q->velocity.x)*(q->velocity.x) + (q->velocity.y)*(q->velocity.y) + (q->velocity.z)*(q->velocity.z));
            //    etot_final=e1+e2;
             //   CkPrintf("<SIDM> (noncosmo) etot_final: %g,  e1: %g, e2: %g \n",etot_final,e1,e2);
                CkPrintf("<SIDM> ran: %g,  probability: %g \n",ran, probability);
                CkPrintf("<SIDM> uvector: %g \n", sqrt(ux*ux + uy*uy + uz*uz) );
                CkPrintf("<SIDM> momentum final: %g,  density: %g \n",pcm0, fPhyDensity);
                CkPrintf("<SIDM> px: %g,  py: %g, pz: %g \n",pxcm,pycm,pzcm);
                }
         
            }
        }    
    }
