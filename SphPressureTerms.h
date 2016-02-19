/* 
   Code originally by James Wadsley.

   Macros for function: 
   void PressureSmoothParams::fcnSmooth(GravityParticle *p, int nSmooth,
                                        pqSmoothNode *nnList)
   in Sph.C

   This code is used 3 times.  Constrained to use Macros (not functions) because 
   Macro definitions change between uses as follows:
   1) with p and q active
     hash define PACTIVE(xxx) xxx
     hash define QACTIVE(xxx) xxx
   2) with just p active
     hash define PACTIVE(xxx) xxx
     hash define QACTIVE(xxx) 
   3) with just q active 
     hash define PACTIVE(xxx) 
     hash define QACTIVE(xxx) xxx

   All Macros and Variables not defined here are defined in Sph.C
 */

#define PRES_PDV(a,b) (a)
#define PRES_ACC(a,b) (a+b)
#define SWITCHCOMBINE(a,b) (0.5*(a->BalsaraSwitch()+b->BalsaraSwitch()))

#ifdef DRHODT
#define DRHODTACTIVE(xxx) xxx
#ifdef RTFORCE
#define RHO_DIVV(a,b) (b)
#else
#define RHO_DIVV(a,b) (a)
#endif
#else
#define DRHODTACTIVE(xxx) 
#endif

#ifdef DIFFUSION 

#ifdef FEEDBACKDIFFLIMIT
#define DIFFUSIONLimitTest() (diffSum == 0 || dTime < p->fTimeCoolIsOffUntil() || dTime < q->fTimeCoolIsOffUntil())
#else
#define DIFFUSIONLimitTest() (diffSum == 0)
#endif


#ifdef DIFFUSIONHARMONIC
#define DIFFUSIONBase() double diffSum = (p->diff()+q->diff());           \
                        double diffBase = (DIFFUSIONLimitTest() ? 0 : 4*p->diff()*q->diff()/diffSum);
#else
#define DIFFUSIONBase() double diffSum = (p->diff()+q->diff()); \
                        double diffBase = (DIFFUSIONLimitTest() ? 0 : diffSum);
#endif
#ifdef MASSDIFF
#define MASSDIFFFAC(pother_) ((pother_)->fMass)
#define DIFFUSIONMetalsBase() double diffMetalsBase = 4*smf->dMetalDiffusionCoeff*diffBase   \
     /((p->fDensity+q->fDensity)*(p->fMass+q->fMass)); 
#define DIFFUSIONMass() \
    { double diff = diffMetalsBase*(p->fMass - q->fMass); \
      PACTIVE( p->fMassDot += diff*p->fMass*rq ); \
      QACTIVE( q->fMassDot -= diff*q->fMass*rp ); \
    }
#define DIFFUSIONVelocity() \
    { double diff0 = diffMetalsBase*(p->v[0] - q->v[0]); \
      double diff1 = diffMetalsBase*(p->v[1] - q->v[1]); \
      double diff2 = diffMetalsBase*(p->v[2] - q->v[2]); \
      PACTIVE( ACCEL(p,0) += diff0*rq*MASSDIFFFAC(q) ); \
      QACTIVE( ACCEL(q,0) -= diff0*rp*MASSDIFFFAC(p) ); \
      PACTIVE( ACCEL(p,1) += diff1*rq*MASSDIFFFAC(q) ); \
      QACTIVE( ACCEL(q,1) -= diff1*rp*MASSDIFFFAC(p) ); \
      PACTIVE( ACCEL(p,2) += diff2*rq*MASSDIFFFAC(q) ); \
      QACTIVE( ACCEL(q,2) -= diff2*rp*MASSDIFFFAC(p) ); \
    }
#else
#define MASSDIFFFAC(pother_) 1
#define DIFFUSIONMetalsBase() double diffMetalsBase = 2*dMetalDiffusionCoeff*diffBase \
     /(p->fDensity+q->fDensity); 
#define DIFFUSIONMass()
#define DIFFUSIONVelocity()
#endif
#else
#define DIFFUSIONBase()
#define DIFFUSIONMetalsBase() 
#define DIFFUSIONMass()
#define DIFFUSIONVelocity()
#endif


#ifdef DIFFUSION
#ifdef UNONCOOL
#define DIFFUSIONThermaluNoncool() \
        { double diffuNc = diffTh*(p->uNoncoolPred-q->uNoncoolPred); \
        PACTIVE( p->uNoncoolDotDiff += diffuNc*rq );        \
        QACTIVE( q->uNoncoolDotDiff -= diffuNc*rp );        \
        }
#else
#define DIFFUSIONThermaluNoncool()  
#endif
#ifdef DIFFUSIONPRICE
#define DIFFUSIONThermal(dt_) \
    { double irhobar = 2/(p->fDensity+q->fDensity);     \
     double vsig = sqrt(fabs(qPoverRho2*q->fDensity*q->fDensity - pPoverRho2*p->fDensity*p->fDensity)*irhobar); \
     double diffTh = smf->dThermalDiffusionCoeff*0.5*(ph+sqrt(0.25*BALL2(q)))*irhobar*vsig; \
     double diffu = diffTh*(p->uPred-q->uPred);             \
     PACTIVE( p->uDotDiff += diffu*rq );                     \
     QACTIVE( q->uDotDiff-= diffu*rp );                     \
     DIFFUSIONThermaluNoncool(); }
#else
#ifndef NODIFFUSIONTHERMAL
/* Default -- thermal diffusion */
#ifdef THERMALCOND 
#if (0)
/* Harmonic average coeff */
#define DIFFUSIONThermalCondBase() double dThermalCondSum = p->fThermalCond + q->fThermalCond; \
    double dThermalCond = ( dThermalCondSum <= 0 ? 0 : 4*p->fThermalCond*q->fThermalCond/(dThermalCondSum*p->fDensity*q->fDensity) );
#else
/* Arithmetic average coeff */
#define DIFFUSIONThermalCondBase(dt_) \
      double dThermalCond = (p->fThermalCond + q->fThermalCond)/(p->fDensity*q->fDensity); \
      if (dThermalCond > 0 && (dt_diff = dtFacDiffusion*ph*p->fThermalLength/(dThermalCond*p->fDensity)) < dt_) dt_ = dt_diff; 

#endif
#else
#define DIFFUSIONThermalCondBase(dt_) double dThermalCond=0;
#endif

#define DIFFUSIONThermal(dt_) \
    { double diffTh = (2*dThermalDiffusionCoeff*diffBase/(p->fDensity+q->fDensity)); \
      double dt_diff, diffu;                                                  \
      DIFFUSIONThermalCondBase(dt_);          \
      if (diffTh > 0 && (dt_diff= dtFacDiffusion*ph*ph/(diffTh*p->fDensity)) < dt_) dt_ = dt_diff; \
      diffu = (diffTh+dThermalCond)*(p->uPred()-q->uPred());            \
      PACTIVE( p->PdV() += diffu*rq*MASSDIFFFAC(q) );                   \
      QACTIVE( q->PdV() -= diffu*rp*MASSDIFFFAC(p) );                \
      DIFFUSIONThermaluNoncool(); }
#else
#define DIFFUSIONThermal(dt_)
#endif
#endif

#define DIFFUSIONMetals() \
    { double diff = diffMetalsBase*(p->fMetals() - q->fMetals()); \
      PACTIVE( p->fMetalsDot() += diff*rq*MASSDIFFFAC(q) ); \
      QACTIVE( q->fMetalsDot() -= diff*rp*MASSDIFFFAC(p) ); }
#define DIFFUSIONMetalsOxygen() \
    { double diff = diffMetalsBase*(p->fMFracOxygen() - q->fMFracOxygen()); \
      PACTIVE( p->fMFracOxygenDot() += diff*rq*MASSDIFFFAC(q) ); \
      QACTIVE( q->fMFracOxygenDot() -= diff*rp*MASSDIFFFAC(p) ); }
#define DIFFUSIONMetalsIron() \
    { double diff = diffMetalsBase*(p->fMFracIron() - q->fMFracIron()); \
      PACTIVE( p->fMFracIronDot() += diff*rq*MASSDIFFFAC(q) ); \
      QACTIVE( q->fMFracIronDot() -= diff*rp*MASSDIFFFAC(p) ); }
#else /* No diffusion */
#define DIFFUSIONThermal(dt_)
#define DIFFUSIONMetals() 
#define DIFFUSIONMetalsOxygen() 
#define DIFFUSIONMetalsIron() 
#endif

#ifdef VARALPHA
#define ALPHA (smf->alpha*0.5*(p->alpha+q->alpha))
#define BETA  (smf->beta*0.5*(p->alpha+q->alpha))
#else
#define ALPHA (alpha)
#define BETA  (beta)
#endif
#define SETDTNEW_PQ(dt_)  { if (dt_ < p->dtNew()) p->dtNew()=dt_; \
                            if (dt_ < q->dtNew()) q->dtNew()=dt_; \
                            if (4*q->dt < p->dtNew()) p->dtNew() = 4*q->dt; \
                            if (4*p->dt < q->dtNew()) q->dtNew() = 4*p->dt; }
// Question: what is VSIGVISC?          
#ifdef VSIGVISC
#define ARTIFICIALVISCOSITY(visc_,dt_) { absmu = -dvdotdr*smf->a           \
            /sqrt(nnList[i].fDist2); /* mu multiply by a to be consistent with physical c */ \
        if (absmu>p->mumax) p->mumax=absmu; /* mu terms for gas time step */ \
		if (absmu>q->mumax) q->mumax=absmu; \
		visc_ = (ALPHA*(pc + q->c) + BETA*1.5*absmu); \
		dt_ = dtFacCourant*ph/(0.625*(pc + q->c())+0.375*visc_); \
		visc_ = SWITCHCOMBINE(p,q)*visc_ \
		    *absmu/(pDensity + q->fDensity); }
#else
/* Iryna's Artifical Viscosity */
#define ARTIFICIALVISCOSITY(visc_,dt_) { double hav=0.5*(ph+0.5*q->fBall);  /* h mean */ \
		absmu = -hav*dvdotdr*a  \
		  /(fDist2+0.01*hav*hav); /* mu multiply by a to be consistent with physical c*/ \
		if (absmu>p->mumax()) p->mumax()=absmu; /* mu terms for gas time step  */  \
		if (absmu>q->mumax()) q->mumax()=absmu;                 \
		visc_ = (-(pc + q->c()) + (BETA/ALPHA)*2*absmu);	\
		dt_ = dtFacCourant*hav/(0.625*(pc + q->c())+0.375*visc_); \
		visc_ = SWITCHCOMBINE(p,q)*visc_ \
		    *absmu/(pDensity + q->fDensity); }
//*/
#endif

    /* Force Calculation between particles p and q */
        DRHODTACTIVE( PACTIVE( p->fDivv_PdV -= rq/p->fDivv_Corrector/RHO_DIVV(pDensity,q->fDensity)*dvdotdr; )); 
        DRHODTACTIVE( QACTIVE( q->fDivv_PdV -= rp/p->fDivv_Corrector/RHO_DIVV(q->fDensity,pDensity)*dvdotdr; )); 
        DRHODTACTIVE( PACTIVE( p->fDivv_PdVcorr -= rq/RHO_DIVV(pDensity,q->fDensity)*dvdotdr; ));
        DRHODTACTIVE( QACTIVE( q->fDivv_PdVcorr -= rp/RHO_DIVV(q->fDensity,pDensity)*dvdotdr; ));
        PACTIVE( p->PdV() += rq*PRES_PDV(pPoverRho2,qPoverRho2)*dvdotdr; );
        QACTIVE( q->PdV() += rp*PRES_PDV(qPoverRho2,pPoverRho2)*dvdotdr; );
        PACTIVE( Accp = (PRES_ACC(pPoverRho2f,qPoverRho2f)); );
        QACTIVE( Accq = (PRES_ACC(qPoverRho2f,pPoverRho2f)); );
        if (dvdotdr>=0.0) {
            dt = dtFacCourant*ph/(2*(pc > q->c() ? pc : q->c()));
            }
        else {  
            hav=0.5*(ph+0.5*q->fBall);
            ARTIFICIALVISCOSITY(visc,dt); /* Calculate Artificial viscosity terms */		
            PACTIVE( p->PdV() += rq*(0.5*visc)*dvdotdr *p->CullenAlpha() / (hav*hav); );
            QACTIVE( q->PdV() += rp*(0.5*visc)*dvdotdr * q->CullenAlpha()/ (hav*hav); );
            PACTIVE( Accp += visc*p->CullenAlpha() / (hav*hav); );
            QACTIVE( Accq += visc*p->CullenAlpha()/ (hav*hav); );
            }
        PACTIVE( Accp *= rq*aFac; );/* aFac - convert to comoving acceleration */
        QACTIVE( Accq *= rp*aFac; );
        PACTIVE( p->treeAcceleration.x -= Accp * dx; );
        PACTIVE( p->treeAcceleration.y -= Accp * dy; );
        PACTIVE( p->treeAcceleration.z -= Accp * dz; );
        QACTIVE( q->treeAcceleration.x += Accq * dx; );
        QACTIVE( q->treeAcceleration.y += Accq * dy; );
        QACTIVE( q->treeAcceleration.z += Accq * dz; );

        DIFFUSIONBase();
        DIFFUSIONThermal(dt);
        DIFFUSIONMetalsBase();
        DIFFUSIONMetals();
        DIFFUSIONMetalsOxygen();
        DIFFUSIONMetalsIron();
        DIFFUSIONMass();
        DIFFUSIONVelocity(); 
#ifdef DTADJUST
        SETDTNEW_PQ(dt);
#endif
