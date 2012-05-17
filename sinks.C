/*
 * Sink particle module originally written for GASOLINE by James
 * Wadsley and Jillian Bellovary (Black holes).
 * Modified for inclusion in ChaNGa
 */
#include <math.h>
#include "ParallelGravity.h"
#include "smooth.h"
#include "sinks.h"

/// Utility function
static inline double sqr(double x) 
{
    return x*x;
    }

///
/// @brief initialize parameters for sinks
///

void Sinks::AddParams(PRM prm, Parameters &params) 
{
	//
	// Sink parameters
	//
	bDoSinks = 0;
	prmAddParam(prm,"bDoSinks",paramBool,&bDoSinks,sizeof(int),
		    "sinks","enable/disable sinks = -sinks");
	bBHSink = 0;
	prmAddParam(prm,"bBHSink",paramBool,&bBHSink,sizeof(int),
		    "bhsink","Bondi-Hoyle type sink = -bhsink");
	bSmallBHSmooth = 0;
	prmAddParam(prm,"bSmallBHSmooth",paramBool,&bSmallBHSmooth,
		    sizeof(int), "smallbhsmooth", "smooth BH feedback over blastwave or smoothing radius = -smallbhsmooth");
	bBHTurnOffCooling = 0;
	prmAddParam(prm,"bBHTurnOffCooling",paramBool,&bBHTurnOffCooling,
		    sizeof(int), "bhturnoffcooling",
		    "turn off cooling for BHs = -bhturnoffcooling");
	bDoBHKick = 0;
	prmAddParam(prm,"bDoBHKick",paramBool,&bDoBHKick,sizeof(int),
		    "dobhkick","turn on grav recoil for mergers = -dobhkick");
	dBHSinkEddEff = 0.1;
	prmAddParam(prm,"dBHSinkEddEff",paramDouble,&dBHSinkEddEff,
		    sizeof(double),"bhsinkeddeff",
		    "<BHSink Eddington Efficiency>");
	dBHSinkFeedbackEff = 0.05;
	prmAddParam(prm,"dBHSinkFeedbackEff",paramDouble,
		    &dBHSinkFeedbackEff,sizeof(double),"bhsinkfbeff",
		    "<BHSink Feedback Efficiency>");
	dBHSinkAlpha = 1;
	prmAddParam(prm,"dBHSinkAlpha",paramDouble,&dBHSinkAlpha,
		    sizeof(double),"bhsinkalpha", "<BHSink Alpha>");
	bBHMindv = 0;
	prmAddParam(prm,"bBHMindv",paramBool,&bBHMindv,sizeof(int),
		    "bhmindv","use mindeltav for BH accretion = -bhmindv");
	bDoSinksAtStart = 0;
	prmAddParam(prm,"bDoSinksAtStart",paramBool,&bDoSinksAtStart,
		    sizeof(int),
		    "sinksas","enable/disable sinks at start = -sinksas");
	bSinkThermal = 0;
	prmAddParam(prm,"bSinkThermal",paramBool,&bSinkThermal,
		    sizeof(int), "tsinks",
		    "enable/disable thermal energy in sink calcs = -tsinks");
	dSinkRadius = 0.0;
	prmAddParam(prm,"dSinkRadius",paramDouble,&dSinkRadius,
		    sizeof(double),"sinkr", "<Sink Radius>");
	dSinkBoundOrbitRadius = 0.0;
	prmAddParam(prm,"dSinkBoundOrbitRadius",paramDouble,
		    &dSinkBoundOrbitRadius,sizeof(double),"sinkbor",
		    "<Sink Bound Orbit Radius>");
	dSinkMustAccreteRadius = 0.0;
	prmAddParam(prm,"dSinkMustAccreteRadius",paramDouble,
		    &dSinkMustAccreteRadius,sizeof(double),"sinkr",
		    "<Sink Radius>");

	dDeltaSink = params.dDelta;
	prmAddParam(prm,"dDeltaSink", paramDouble, &dDeltaSink,
		    sizeof(double), "dDeltaSink",
		    "<Maximum sink timestep in years> = dDelta");
	dSinkMassMin = 0;  /* Default reset to FLT_MAX for BH sink */
	prmAddParam(prm,"dSinkMassMin", paramDouble, &dSinkMassMin,
		    sizeof(double), "dSinkMassMin",
		    "<Minimum Mass to act as a sink> = 0" );
	iSinkRung = 0; 
	prmAddParam(prm,"iSinkRung", paramInt, &iSinkRung, sizeof(int),
		    "iSinkRung", "<Sink Rung> = 0");
	bSinkForm = 0;
	prmAddParam(prm,"bSinkForm",paramBool,&bSinkForm,sizeof(int),
				"sinkform","enable/disable sinks = -sinkform");
	nJeans = 50; 
	prmAddParam(prm,"nJeans", paramInt, &nJeans, sizeof(int),
		    "nJeans",
		    "<no. particle masses inside Jeans mass for good resolution> = 50");
	dJeansConstant = 4/3.*pow(M_PI,2.5);
	prmAddParam(prm,"dJeansConstant", paramDouble, &dJeansConstant,
		    sizeof(double), "dJeansConstant",
		    "<Constant to multiply c_s^3 G^-3/2 rho^-1/2 to get Jeans Mass > = 4/3 pi^5/2");
	bSinkFormJeans = 0;
	prmAddParam(prm,"bSinkFormJeans",paramBool,&bSinkFormJeans,
		    sizeof(int),
		    "sinkformjeans","enable/disable sinks = -sinkformjeans");
	bSinkFormDivV = 0;
	prmAddParam(prm,"bSinkFormDivV",paramBool,&bSinkFormDivV,
		    sizeof(int),
		    "sinkformdivv","enable/disable sinks = -sinkformdivv");
	dSinkFormDivVCoeff = 0;
	prmAddParam(prm,"dSinkFormDivVCoeff",paramDouble,
		    &dSinkFormDivVCoeff,sizeof(double),
		    "sinkformdivvcoeff","sinks form div v coeff = -sinkformdivvcoeff");
	bSinkFormDivAcc = 0;
	prmAddParam(prm,"bSinkFormDivAcc",paramBool,&bSinkFormDivAcc,
		    sizeof(int),
		    "sinkformdivacc","enable/disable sinks = -sinkformdivacc");
	dSinkFormDivAccCoeff = 0;
	prmAddParam(prm,"dSinkFormDivAccCoeff",paramDouble,
		    &dSinkFormDivAccCoeff,sizeof(double),
		    "sinkformdivvacccoeff","sinks form div acc coeff = -sinkformdivacccoeff");
	bSinkFormDV = 0;
	prmAddParam(prm,"bSinkFormDV",paramBool,&bSinkFormDV,sizeof(int),
		    "sinkformdv","enable/disable sinks = -sinkformdv");
	bSinkFormPotMin = 0;
	prmAddParam(prm,"bSinkFormPotMin",paramBool,&bSinkFormPotMin,
		    sizeof(int),
		    "sinkformpotmin","enable/disable sinks = -sinkformpotmin");
	dSinkFormDensity = -1.0;
	prmAddParam(prm,"dSinkFormDensity",paramDouble,
		    &dSinkFormDensity,sizeof(double),
		    "sinkformdensity","sinks density = -sinkformdensity");
	dSinkTimeEligible = 0;
	prmAddParam(prm,"dSinkTimeEligible",paramDouble,
		    &dSinkTimeEligible,sizeof(double),
		    "sinktimeligible","sink time eligible = -sinktimeeligible");
	bSinkFormSimple = 0;
	prmAddParam(prm,"bSinkFormSimple",paramBool,&bSinkFormSimple,
		    sizeof(int),
		    "sinkformsimple","enable/disable sinks = -sinkformjeans");
	nSinkFormMin = 0;
	prmAddParam(prm,"nSinkFormMin",paramInt,&nSinkFormMin,
		    sizeof(int), "sinkformmin","nmin to form a sink");
    }

#include "physconst.h"

///
/// @brief check sink parameters
///

void Sinks::CheckParams(PRM prm, struct parameters &param) 
{
    if (!prmSpecified(prm, "nSinkFormMin")) nSinkFormMin = param.nSmooth;

    if (bBHSink) {
	assert (prmSpecified(prm, "dMsolUnit") &&
		prmSpecified(prm, "dKpcUnit"));

	/* For BH sinks -- default to negative tform as sink indicator */
	if(!prmSpecified(prm, "dSinkMassMin"))
	    dSinkMassMin = FLT_MAX;

	if(!param.bDoGas) {
	    ckerr << "Warning: BH sinks set without enabling SPH" << endl;
	    ckerr << "Enabling SPH" << endl;
	    param.bDoGas = 1;
	    }
	bDoSinks = 1;
	/* Units of inverse time -- code units */
	dBHSinkEddFactor = GCGS*4*M_PI*MHYDR/
	    (SIGMAT*LIGHTSPEED*dBHSinkEddEff)*param.dSecUnit;
	/* c^2 times efficiency factor (ergs per g) -- code units */
	dBHSinkFeedbackFactor = dBHSinkFeedbackEff
	    *dBHSinkEddEff*LIGHTSPEED*LIGHTSPEED/param.dErgPerGmUnit;
	}
    if (bDoSinks) {
	double testDelta;
	if (prmSpecified(prm,"iSinkRung")) {
	    /* Find associated timestep for iSinkRung */
	    int iRung;
	    
	    if (iSinkRung > param.iMaxRung) 
		iSinkRung = param.iMaxRung;

	    testDelta = param.dDelta;
	    for ( iRung = 0; iRung < iSinkRung ; iRung++ ) 
		testDelta *= 0.5;
	    }
	else if (prmSpecified(prm,"dDeltaSink")) {
	    /* Find associate Rung for dDeltaSink */
	    /* NB: dDeltaSink is always in code units */
	    iSinkRung = 0;
	    for ( testDelta = param.dDelta; testDelta > dDeltaSink ;
		  testDelta *= 0.5 ) {
		iSinkRung++;
		if (iSinkRung >= param.iMaxRung)
		    CkAbort("Sink timestep too small");
		}
	    }		  
	else {
	    testDelta = param.dDelta;
	    iSinkRung = 0;
	    }
	CkPrintf("dDeltaSink (set): %g, effectively: %g = %g yrs, iSinkRung: %i",
		 dDeltaSink, testDelta,
		 testDelta*param.dSecUnit/SECONDSPERYEAR, iSinkRung );
	dDeltaSink = testDelta;
	}
    }

///
/// @brief Initial identify sinks
///
void Main::SetSink()
{
    if (param.sinks.bDoSinks) {
	CkReductionMsg *msg;
	treeProxy.SetSink(param.sinks.dSinkMassMin,
			  CkCallbackResumeThread((void*&)msg));
	nSink = *(int *)msg->getData();
	if (verbosity) CkPrintf("Identified %d sink particles\n", nSink);
	delete msg;
	}
    }

///
/// @brief Initial identify sinks; processor specific method.
///
void TreePiece::SetSink(double dSinkMassMin, const CkCallback &cb)
{
    int i,nSink = 0;

    for(unsigned int i = 1; i <= myNumParticles; ++i) {
	GravityParticle *p = &myParticles[i];
#ifdef STARSINK
	if ((TYPETest(p,TYPE_STAR)))
#else
	if ((TYPETest(p,TYPE_STAR) && p->fTimeForm() < 0))
#endif
	    {
		TYPESet(p,TYPE_SINK);
		nSink++;
		}
	}

    contribute(sizeof(int), &nSink, CkReduction::sum_int, cb);
    }
    
///
/// @brief Form sink particles; main routine
///
void Main::FormSinks(double dTime, double dDelta, int iKickRung)
/* For momentum conservation, only particles on iKickRung or higher may
   contribute mass to a sink.
   This call is timed for just after those particles have completed a full KDK
*/
    {
    int64_t nStarOld,nSinkAdd;
    double sec,dsec;

    if (!param.sinks.bSinkForm) return;

    sec = CkWallTimer();
 
    nStarOld = nTotalStar;
    param.sinks.dSinkCurrentDelta = dDelta;
    param.sinks.iSinkCurrentRung = iKickRung;

    int bJeans = param.sinks.bSinkFormJeans;
    int bSimple = param.sinks.bSinkFormSimple;
    double dJConst2 = param.sinks.dJeansConstant/param.sinks.nJeans;
    dJConst2 *= dJConst2;
    double dDensityCut = param.sinks.dSinkFormDensity/param.dGmPerCcUnit;
    int bDensity = (dDensityCut > 0 ? 1 : 0);
    CkReductionMsg *msg;
    treeProxy.formSinks(bJeans, dJConst2, bDensity, dDensityCut, dTime,
			iKickRung, bSimple,
			CkCallbackResumeThread((void*&)msg));
    int nCandidates = *(int *)msg->getData();
    delete msg;
    
    nSinkAdd = 0;
    if (nCandidates) {
	if (param.sinks.bSinkFormSimple) {
	    nSinkAdd = nCandidates;
	    nSink += nSinkAdd;
	    } 
	else {
	    int iPhase = 0; // Keeps track of node cache use
	/* Note: Only gas particles are combined into sinks */
	    SinkFormTestSmoothParams pSFT(TYPE_GAS, iKickRung, param.sinks);
	    double dfBall2OverSoft2 = 4.0*param.dhMinOverSoft*param.dhMinOverSoft;
	    treeProxy.startIterationSmooth(&pSFT, 1, dfBall2OverSoft2,
					   CkCallbackResumeThread());
	    iPhase++;

	    SinkFormSmoothParams pSF(TYPE_GAS, iKickRung, param.csm, dTime,
				     param.sinks);
	    treeProxy.startIterationSmooth(&pSF, 1, dfBall2OverSoft2,
				   CkCallbackResumeThread());
	    iPhase++;
	    
	    addDelParticles();
	    nSink += nSinkAdd = nTotalStar-nStarOld;
	    }
	}


    if (verbosity)
	CkPrintf("Sinks Formation: %i candidates +%d sinks\n",nCandidates,
		 nSinkAdd);

    dsec = CkWallTimer() - sec;
    CkPrintf("Sink Form Done (+%d total) Calculated, Wallclock: %f secs\n\n",
	     nSinkAdd, dsec);

    }

///
/// @brief Form sink particles: per processor routine
///
void TreePiece::formSinks(int bJeans, double dJConst2, int bDensity,
			  double dDensityCut, double dTime, int iKickRung,
			  int bSimple, const CkCallback &cb)
{
    int i;
    GravityParticle *p;
    int n = myNumParticles;
    double Jval, Jvalmin = FLT_MAX;
    int nCandidates = 0;
    
    for(i = 1; i <= n; ++i) {
	GravityParticle *p = &myParticles[i];
    
        if(p->isGas() && p->rung >= iKickRung) {
#ifdef SINKING
	    if (TYPETest( p, TYPE_SINKING )) continue;
#endif
/* Jeans Mass compared to nJeans particle masses */
	    Jval =  dJConst2*(p->c()*p->c()*p->c()*p->c()*p->c()*p->c())
		/(p->mass*p->mass*p->fDensity);
	    if (Jval < Jvalmin) Jvalmin = Jval;
	    if ((bJeans && Jval < 1) ||
		(bDensity && p->fDensity >= dDensityCut)) {
		nCandidates++;
		if (bSimple) {
		    TYPESet(p, TYPE_SINK); /* Is now a sink! */
		    p->fMetals()= -dTime;
		    }
		else {
		    TYPESet(p, TYPE_SMOOTHACTIVE); /* Is now a candidate */
		    }
		}
	    }
	}
    contribute(cb);
}

#define fSinkRating(_a)    ((_a)->curlv()[0] )
#define fSinkRating_Ex(_a)    ((_a)->curlv[0] )
#define iOrderSink(_a)      (*((int *) (&(_a)->curlv()[1])))
#define iOrderSink_Ex(_a)      (*((int *) (&(_a)->curlv[1])))

void SinkFormTestSmoothParams::initSmoothCache(GravityParticle *p)
{
    fSinkRating(p) = -FLT_MAX;
    iOrderSink(p) = -1;
#ifdef SINKING
    TYPEReset(p, TYPE_NEWSINKING);
#endif
    }

void SinkFormTestSmoothParams::combSmoothCache(GravityParticle *p1,
					       ExternalSmoothParticle *p2)
{
/* Particle p1 belongs to candidate stored in iOrderSink of p1 initially but
   switch to p2's if that candidate is denser */
#ifdef SINKING
    if (TYPETest(p1, TYPE_SINKING)) return;
#endif
    if (fSinkRating_Ex(p2) > fSinkRating(p1)) {
	fSinkRating(p1) = fSinkRating_Ex(p2);
	iOrderSink(p1) = iOrderSink_Ex(p2);
	}
}

void SinkFormTestSmoothParams::fcnSmooth(GravityParticle *p,int nSmooth,
					 pqSmoothNode *nnList)
{
	int i;
	double dSinkRadius2 = s.dSinkRadius*s.dSinkRadius,r2;
	GravityParticle *q;

	/* Apply Bate tests in next phase
	   For now just decide which sink the particle belongs to:
   	          prefer joining a denser (or deeper potential) sink candidate

		  Must be local extremum to be a candidate or abort
	   Need to prevent double counting particles into two sinks
	*/
	for (i=0;i<nSmooth;++i) {
	    r2 = nnList[i].fKey;
	    if (r2 > 0 && r2 <= dSinkRadius2) {
		q = nnList[i].p;
		if (TYPETest( q, TYPE_GAS ) && ((s.bSinkFormPotMin ? q->potential < p->potential : q->fDensity > p->fDensity))) {
		    /* Abort without grabbing any particles -- this isn't an extremum particle */
/*		    CkPrintf("Sink aborted %d %g: Denser Neighbour %d %g\n",p->iOrder,p->fDensity,q->iOrder,q->fDensity);*/
		    return;
			
		    }
		}
	    }

	fSinkRating(p) =  (s.bSinkFormPotMin ? -p->potential : p->fDensity); /* rate yourself */
	
/*	printf("Sink %d %g: looking...\n",p->iOrder,p->fDensity);*/

	/* Identify all nbrs as to be eaten unless those nbrs already
	   belong to a "higher rated" sink 
	   Thus sink quality is rated by density or potential depth
	*/
	for (i=0;i<nSmooth;++i) {
	    r2 = nnList[i].fKey;
	    if (r2 > 0 && r2 <= dSinkRadius2) {
		q = nnList[i].p;
#ifdef SINKING
		if (TYPETest(q, TYPE_SINKING)) continue;
#endif
		if (TYPETest( q, TYPE_GAS ) && q->rung >= s.iSinkCurrentRung) {
		    if (s.bSinkFormPotMin) {
			if (-p->potential > fSinkRating(q)) {
			    fSinkRating(q) = -p->potential;
			    iOrderSink(q) = p->iOrder; /* Particle q belongs to sink p */
			    }

			}
		    else {
			if (p->fDensity > fSinkRating(q)) {
			    fSinkRating(q) = p->fDensity;
			    iOrderSink(q) = p->iOrder; /* Particle q belongs to sink p */
			    }
			}
		    }
		}
	    }
}

#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))

void tred3(double a[][3], double d[], double e[])
{
  int l,k,j,i,n=3;
  double scale,hh,h,g,f;

  for (i=n-1;i>=1;i--) {
    l=i-1;
    h=scale=0.0;
    if (l > 0) {
      for (k=0;k<=l;k++)
        scale += fabs(a[i][k]);
      if (scale == 0.0)
        e[i]=a[i][l];
      else {
        for (k=0;k<=l;k++) {
          a[i][k] /= scale;
          h += a[i][k]*a[i][k];
        }
        f=a[i][l];
        g=(f > 0.0 ? -sqrt(h) : sqrt(h));
        e[i]=scale*g;
        h -= f*g;
        a[i][l]=f-g;
        f=0.0;
        for (j=0;j<=l;j++) {
          a[j][i]=a[i][j]/h;
          g=0.0;
          for (k=0;k<=j;k++)
            g += a[j][k]*a[i][k];
          for (k=j+1;k<=l;k++)
            g += a[k][j]*a[i][k];
          e[j]=g/h;
          f += e[j]*a[i][j];
        }
        hh=f/(h+h);
        for (j=0;j<=l;j++) {
          f=a[i][j];
          e[j]=g=e[j]-hh*f;
          for (k=0;k<=j;k++)
            a[j][k] -= (f*e[k]+g*a[i][k]);
        }
      }
    } else
      e[i]=a[i][l];
    d[i]=h;
  }
  d[0]=0.0;
  e[0]=0.0;
  /* Contents of this loop can be omitted if eigenvectors not
      wanted except for statement d[i]=a[i][i]; */
  for (i=0;i<n;i++) {
    l=i-1;
    if (d[i]) {
      for (j=0;j<=l;j++) {
        g=0.0;
        for (k=0;k<=l;k++)
          g += a[i][k]*a[k][j];
        for (k=0;k<=l;k++)
          a[k][j] -= g*a[k][i];
      }
    }
    d[i]=a[i][i];
    a[i][i]=1.0;
    for (j=0;j<=l;j++) a[j][i]=a[i][j]=0.0;
  }
}

double pythag(double a, double b)
{
  double absa,absb;
  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) {
      double xxx = absb/absa;
      return absa*sqrt(1.0+xxx*xxx);
//      return absa*sqrt(1.0+SQR(absb/absa));
      }
  else if (absb == 0) return 0.;
  else {
      double xxx = absa/absb;
      return absb*sqrt(1.0+xxx*xxx);
//    return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
      }
}

void tqli3(double d[], double e[], double z[][3])
{
  double pythag(double a, double b);
  int m,l,iter,i,k,n=3;
  double s,r,p,g,f,dd,c,b;

  for (i=1;i<n;i++) e[i-1]=e[i];
  e[n-1]=0.0;
  for (l=0;l<n;l++) {
    iter=0;
    do {
      for (m=l;m<n-1;m++) {
        dd=fabs(d[m])+fabs(d[m+1]);
        if ((double)(fabs(e[m])+dd) == dd) break;
      }
      if (m != l) {
	  if (iter++ == 30) assert(iter<30); /* "Too many iterations in tqli" */
        g=(d[l+1]-d[l])/(2.0*e[l]);
        r=pythag(g,1.0);
        g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
        s=c=1.0;
        p=0.0;
        for (i=m-1;i>=l;i--) {
          f=s*e[i];
          b=c*e[i];
          e[i+1]=(r=pythag(f,g));
          if (r == 0.0) {
            d[i+1] -= p;
            e[m]=0.0;
            break;
          }
          s=f/r;
          c=g/r;
          g=d[i+1]-p;
          r=(d[i]-g)*s+2.0*c*b;
          d[i+1]=g+(p=s*r);
          g=c*r-b;
          for (k=0;k<n;k++) {
            f=z[k][i+1];
            z[k][i+1]=s*z[k][i]+c*f;
            z[k][i]=c*z[k][i]-s*f;
          }
        }
        if (r == 0.0 && i >= l) continue;
        d[l] -= p;
        e[l]=g;
        e[m]=0.0;
      }
    } while (m != l);
  }
}

void SinkFormSmoothParams::combSmoothCache(GravityParticle *p1,
					   ExternalSmoothParticle *p2)
{
    if (!(TYPETest( p1, TYPE_DELETED )) && TYPETest( p2, TYPE_DELETED ) ) {
	p1->mass = p2->mass;
	deleteParticle(p1);
	}
    }

void SinkFormSmoothParams::fcnSmooth(GravityParticle *p,int nSmooth,
				     pqSmoothNode *nnList)
{
	int i,j,nEaten,nEatNow;
	double mtot,im,Ek,Eth,Eg,dvx,dv2,vsink[3];
	double diva,ih2,r2,rs1;
	GravityParticle *q,*q1,*q2;
	GravityParticle sinkp;

	/* You are accreted */
	if ( iOrderSink(p) != -1 ) {
	    /* What if it was temporarily a candidate -- are it's nbrs confused? 
	       Shouldn't be -- they should go to a new better candidate if they are near one */

	    return;
	    }
#ifdef SINKING
	assert(!TYPETest(p,TYPE_SINKING));
#endif
	iOrderSink(p) = p->iOrder;

	nEaten=0;
	mtot = 0;
	for (i=0;i<nSmooth;++i) {
	    q1 = nnList[i].p;
	    if (iOrderSink(q1) == p->iOrder) {
		nEaten++;
		mtot += q1->mass;
		}
	    }

	if (nEaten < s.nSinkFormMin) {
	    return;
	    }

	if (s.bSinkFormDivAcc) {
	    diva = 0;
	    ih2 = 1./(s.dSinkRadius*s.dSinkRadius);
	    for (i=0;i<nSmooth;++i) {
		q1 = nnList[i].p;
		if (iOrderSink(q1) == p->iOrder) {
		    double dax,day,daz,dx,dy,dz;
		    
		    r2 = nnList[i].fKey*ih2;  /* Calculate div a (un-normalized) */
		    if (r2 < 4) {
			rs1 = DKERNEL(r2);
			rs1 *= q1->mass;
			dx = nnList[i].dx.x;
			dy = nnList[i].dx.y;
			dz = nnList[i].dx.z;
			dax = (-p->treeAcceleration[0]
			       + q1->treeAcceleration[0]);
			day = (-p->treeAcceleration[1]
			       + q1->treeAcceleration[1]);
			daz = (-p->treeAcceleration[2]
			       + q1->treeAcceleration[2]);
			diva += (dax*dx+day*dy+daz*dz)*rs1;
			}
		    }
		}
	    if (diva > 0) return;
	    }


	if (s.bSinkFormDV || s.bSinkFormDivV) {
	    double dvdx[3][3];
	    double d[3],e[3];
	    double vFac = 1./(a*a); /* converts v to xdot */

	    dvdx[0][0]=0; dvdx[0][1]=0; dvdx[0][2]=0;
	    dvdx[1][0]=0; dvdx[1][1]=0; dvdx[1][2]=0;
	    dvdx[2][0]=0; dvdx[2][1]=0; dvdx[2][2]=0;

	    ih2 = 1./(s.dSinkRadius*s.dSinkRadius);
	    for (i=0;i<nSmooth;++i) {
		q1 = nnList[i].p;
		if (iOrderSink(q1) == p->iOrder) {
		    double dvx,dvy,dvz,dx,dy,dz;
		    
		    r2 = nnList[i].fKey*ih2;  /* Calculate dv_i/dx_j (un-normalized) */
		    if (r2 < 4) {
			rs1 = DKERNEL(r2);
			rs1 *= q1->mass;
			dx = nnList[i].dx.x;
			dy = nnList[i].dx.y;
			dz = nnList[i].dx.z;
			dvx = ((-p->vPred()[0] + q1->vPred()[0])*vFac - dx*H)*rs1; /* NB: dx = px - qx */
			dvy = ((-p->vPred()[1] + q1->vPred()[1])*vFac - dy*H)*rs1;
			dvz = ((-p->vPred()[2] + q1->vPred()[2])*vFac - dz*H)*rs1;
			dvdx[0][0] += dvx*dx;
			dvdx[0][1] += dvx*dy;
			dvdx[0][2] += dvx*dz;
			dvdx[1][0] += dvy*dx;
			dvdx[1][1] += dvy*dy;
			dvdx[1][2] += dvy*dz;
			dvdx[2][0] += dvz*dx;
			dvdx[2][1] += dvz*dy;
			dvdx[2][2] += dvz*dz;
			}
		    }
		}
		
	    dvdx[0][1] = 0.5*(dvdx[0][1]+dvdx[1][0]); dvdx[1][0]=dvdx[0][1];
	    dvdx[2][1] = 0.5*(dvdx[2][1]+dvdx[1][2]); dvdx[1][2]=dvdx[2][1];
	    dvdx[2][0] = 0.5*(dvdx[2][0]+dvdx[0][2]); dvdx[0][2]=dvdx[2][0];

	    if (s.bSinkFormDivV) {
		if (dvdx[0][0]+dvdx[1][1]+dvdx[2][2] > 0) return; /* not converging */
		}
	    if (s.bSinkFormDV) {
		tred3(dvdx,d,e);
		tqli3(d,e,dvdx);
	    
		if (d[0] > 0 || d[1] > 0 || d[2] > 0) return; /* not converging in all directions */
		}
	    }

	Ek = 0;
	Eth = 0;
	Eg = 0;
	vsink[0] = 0; 	vsink[1] = 0;	vsink[2] = 0;

	for (i=0;i<nSmooth;++i) {
	    q1 = nnList[i].p;
	    if (iOrderSink(q1) == p->iOrder) {
		dvx = q1->velocity[0]-p->velocity[0];  
		vsink[0] += q1->mass*dvx;
		dv2 = dvx*dvx;
		dvx = q1->velocity[1]-p->velocity[1];
		vsink[1] += q1->mass*dvx;
		dv2 += dvx*dvx;
		dvx = q1->velocity[2]-p->velocity[2];
		vsink[2] += q1->mass*dvx;
		dv2 += dvx*dvx;
		Ek += 0.5*q1->mass*dv2;
		Eth += q1->mass*q1->u();
		for (j=i+1;j<nSmooth;j++) {
		    q2 = nnList[j].p;
		    if (iOrderSink(q2) == p->iOrder) {
			dvx = q1->position[0]-q2->position[0];
			dv2 = dvx*dvx;
			dvx = q1->position[1]-q2->position[1];
			dv2 += dvx*dvx;
			dvx = q1->position[2]-q2->position[2];
			dv2 += dvx*dvx;
			Eg -= q1->mass*q2->mass/sqrt(dv2);
			}
		    }
		}
	    }

	Ek -= 0.5*(vsink[0]*vsink[0]+vsink[1]*vsink[1]+vsink[2]*vsink[2])/mtot;

	/* Apply Bate tests here -- 
	   1. thermal energy < 1/2 |Grav|, 
	   2. thermal + rot E < |Grav|, 
	   3. total E < 0 (implies 2.)
	   4. div.acc < 0 (related to rate of change of total E I guess)  (optional see above)
           5. NEW eigenvalues of dvi/dxj all negative (Federrath) (optional see above)
           6. NEW div.v < 0 (Whitworth group)
	*/

	if (Eth < 0.5*fabs(Eg) && Ek + Eth + Eg < 0) {
	    /* Sink approved */	
	    if (s.dSinkTimeEligible > 0) {
		/* Check that formation conditions persisted for minimum time */
		/* This cannot work in its simplest form because the sink candidate
                   for a given density peak changes randomly from step to step as
                   densities etc.. fluctuate.  Need to code in some way to continuously
                   refer to the same pre-sink blob (inherit properties of particles inside 1/2 SinkRadius?) */
		if (p->fTimeForm() > 1e36) {
		    p->fTimeForm() = dTime;
		    return;
		    }
		if (dTime-p->fTimeForm() < s.dSinkTimeEligible) return;
		}

	    GravityParticle sinkp = *p;
	    sinkp.position[0] = 0;
	    sinkp.position[1] = 0;
	    sinkp.position[2] = 0;
	    sinkp.velocity[0] = 0;
	    sinkp.velocity[1] = 0;
	    sinkp.velocity[2] = 0;
	    sinkp.treeAcceleration[0] = 0;
	    sinkp.treeAcceleration[1] = 0;
	    sinkp.treeAcceleration[2] = 0;
	    sinkp.u() = 0;
#ifdef STARSINK
	    SINK_Lx(&sinkp) = 0;
	    SINK_Ly(&sinkp) = 0;
	    SINK_Lz(&sinkp) = 0;
#endif
	    sinkp.mass = 0;
#ifdef SINKING
	    sinkp.fTrueMass = 0;
	    sinkp.iSinkingOnto = 0; /* Must start at zero -- it records the number of assoc. sinking particles */
#endif
	    nEatNow = 0;
	    for (i=0;i<nSmooth;++i) {
		q = nnList[i].p;
		r2 = nnList[i].fKey;
		if (iOrderSink(q) == p->iOrder) {
#ifdef SINKING
		    if (r2 > smf->dSinkMustAccreteRadius*smf->dSinkMustAccreteRadius) {
			/* Make sinking -- can't force sink choice b/c sink iOrder isn't known yet */
			TYPESet(q, TYPE_NEWSINKING);
			}
		    else 
#endif /*SINKING*/
			{
			double mq,rx,ry,rz,vx,vy,vz;
			mq = q->mass;
			rx = q->position[0]; ry = q->position[1];
			rz = q->position[2];
			vx = q->velocity[0]; vy = q->velocity[1];
			vz = q->velocity[2];
			nEatNow++;
			sinkp.position[0] += mq*rx;
			sinkp.position[1] += mq*ry;
			sinkp.position[2] += mq*rz;
			sinkp.velocity[0] += mq*vx;
			sinkp.velocity[1] += mq*vy;
			sinkp.velocity[2] += mq*vz;
			sinkp.treeAcceleration[0] += mq*q->treeAcceleration[0];
			sinkp.treeAcceleration[1] += mq*q->treeAcceleration[1];
			sinkp.treeAcceleration[2] += mq*q->treeAcceleration[2];
			sinkp.u() += q->mass*q->u();
#ifdef STARSINK
			SINK_Lx(&sinkp) += mq*(ry*vz - rz*vy); /* Add to total L */
			SINK_Ly(&sinkp) += mq*(rz*vx - rx*vz);
			SINK_Lz(&sinkp) += mq*(rx*vy - ry*vx);
#endif
			sinkp.mass += mq; 
#ifdef SINKING
			sinkp.fTrueMass += mq;
#endif
			if (p!=q) {
			    q->mass = 0;
			    deleteParticle(q);
			    }
			}
		    }
		}   

	    im = 1/sinkp.mass;
	    sinkp.position[0] *= im;
	    sinkp.position[1] *= im;
	    sinkp.position[2] *= im;
	    sinkp.velocity[0] *= im;
	    sinkp.velocity[1] *= im;
	    sinkp.velocity[2] *= im;
	    sinkp.treeAcceleration[0] *= im;
	    sinkp.treeAcceleration[1] *= im;
	    sinkp.treeAcceleration[2] *= im;
	    sinkp.u() *= im;
#ifdef STARSINK
            /* Store Internal angular momentum only */
	    SINK_Lx(&sinkp) -= sinkp.fMass*(sinkp.r[1]*sinkp.v[2] - sinkp.r[2]*sinkp.v[1]); 
	    SINK_Ly(&sinkp) -= sinkp.fMass*(sinkp.r[2]*sinkp.v[0] - sinkp.r[0]*sinkp.v[2]);
	    SINK_Lz(&sinkp) -= sinkp.fMass*(sinkp.r[0]*sinkp.v[1] - sinkp.r[1]*sinkp.v[0]);
#endif
	    TYPEClear(&sinkp);
	    TYPESet(&sinkp,TYPE_SINK|TYPE_STAR);
	    sinkp.fTimeForm() = -dTime; /* -ve time is sink indicator */
#ifdef SINKEXTRADATA
	    for(j = 0; j < 3; j++) {
		sinkp.rForm[j] = sinkp.r[j];
		sinkp.vForm[j] = sinkp.v[j];
		}
#endif
	    printf("Sink Formed %d %g: np (%d) %d Mass (%g) %g Ek %g Eth %g Eg %g, %g %g\n",p->iOrder,p->fDensity,nEaten,nEatNow,mtot,sinkp.mass,Ek,Eth,Eg,4/3.*pow(M_PI,2.5)/50.*sqrt((p->c()*p->c()*p->c()*p->c()*p->c()*p->c())/(p->mass*p->mass*p->fDensity)),pow(Eth/fabs(Eg),1.5) );
#ifndef SINKING
	    assert(fabs(sinkp.mass/mtot-1) < 1e-4);
#else
/*	    printf("Sink Made %g %g\n",sinkp.fMass,sinkp.fTrueMass);*/
	    assert(fabs(sinkp.fMass/sinkp.fTrueMass-1) < 1e-7);
	    assert(!TYPETest(p, TYPE_SINKING));
#endif
	    p->mass = 0;
	    deleteParticle(p);
	    newParticle(&sinkp);    
	    }
	else {
	    printf("Sink Failed Tests %d %g: np (%d) Mass (%g) Ek %g Eth %g Eg %g, %g %g\n",p->iOrder,p->fDensity,nEaten,mtot,Ek,Eth,Eg, 4/3.*pow(M_PI,2.5)/50.*sqrt((p->c()*p->c()*p->c()*p->c()*p->c()*p->c())/(p->mass*p->mass*p->fDensity)),pow(Eth/fabs(Eg),1.5) );
	    }
}

/// Process sink particles
void Main::doSinks(double dTime, double dDelta, int iKickRung)
/* For momentum conservation, only particles on iKickRung or higher may
   contribute mass to a sink 
   This call is timed for just after those particles have completed a full KDK
*/
{
   double sec,sec1,dsec,dMass;
   int nAccreted,nSmoothTemp;

    /* I assume sink creation is rarer so the tree will be ok after this call most of the time */
    FormSinks(dTime, dDelta, iKickRung ); 
   
    if (param.sinks.bDoSinks == 0) return;
    if (nSink == 0) return;
    if (param.sinks.bBHSink && dDelta <= 0.0) return;

    sec = CkWallTimer();
    // dMass = msrMassCheck(msr, -1, "Accrete onto Sinks: Initial Value");

    nAccreted = nTotalSPH;
    if(nActiveGrav > 0) {
	param.sinks.dSinkCurrentDelta = dDelta;
	param.sinks.iSinkCurrentRung = iKickRung;
	double dfBall2OverSoft2 = 4.0*param.dhMinOverSoft*param.dhMinOverSoft;

	if (param.sinks.bBHSink) {
	    /* Smooth Bondi-Hoyle Accretion: radius set by nSmooth */
	    BHDensitySmoothParams pBHDen(TYPE_GAS, iKickRung, param.csm,
					 dTime, param.sinks);
	    treeProxy.startIterationSmooth(&pBHDen, 1, dfBall2OverSoft2,
					     CkCallbackResumeThread());
	    iPhase++;
	    BHAccreteSmoothParams pBHAcc(TYPE_GAS, iKickRung, param.csm,
					 dTime, param.sinks,
					 param.dConstGamma);
	    treeProxy.startIterationReSmooth(&pBHAcc,
					     CkCallbackResumeThread());
	    iPhase++;
	    /* build new tree of BHs for merging JMB 11/14/08  */
	    /* Search for BHs that need merging */
	    CkReductionMsg *msgCnt;
	    treeProxy.countType(TYPE_SINK,
				CkCallbackResumeThread((void *&)msgCnt));
	    nSink = *(int *) msgCnt->getData();
	    delete msgCnt;
	    if (nSink > 1) { /* no need to merge if there is only one! */
		BHIdentifySmoothParams pBHID(TYPE_SINK, iKickRung, param.csm,
					     dTime, param.sinks);
		treeProxy.startIterationSmooth(&pBHDen, 1, dfBall2OverSoft2,
					       nSink,
					       CkCallbackResumeThread());
		iPhase++;
		BHSinkMergeSmoothParams pBHM(TYPE_SINK, iKickRung, param.csm,
					     dTime, param.sinks,
					     param.dConstGamma);
		treeProxy.startIterationReSmooth(&pBHM,
						 CkCallbackResumeThread());
		iPhase++;
		}
	    }
	else {
	    /* Fixed Radius Accretion: particle by particle (cf. Bate) */
	    SinkAccreteTestSmoothParams pSinkTest(TYPE_GAS, iKickRung, dTime,
						  param.sinks);
	    treeProxy.startIterationSmooth(&pSinkTest, 1, dfBall2OverSoft2,
					   CkCallbackResumeThread());
	    iPhase++;
	    
#ifdef SINKINGAVERAGE
	    msrActiveType(msr,TYPE_NEWSINKING, TYPE_SMOOTHACTIVE);
	    msrSmooth(msr, dTime, SMX_SINKINGAVERAGE,0);
	    SinkingAverageSmoothParams pSinkAvg(TYPE_GAS, iKickRung, sinks);
	    treeProxy.startIterationSmooth(&pSinkAvg, 1, dfBall2OverSoft2,
					   CkCallbackResumeThread());
	    iPhase++;
#endif
	    SinkAccreteSmoothParams pSinkAcc(TYPE_GAS, iKickRung, dTime,
					     param.sinks);
	    treeProxy.startIterationSmooth(&pSinkAcc, 1, dfBall2OverSoft2,
					   CkCallbackResumeThread());
	    iPhase++;
	    }
	
	addDelParticles();
	if (param.sinks.bBHSink) { /* reset nSinks in case of merger */
	    CkReductionMsg *msgCnt;
	    treeProxy.countType(TYPE_SINK,
				CkCallbackResumeThread((void *&)msgCnt));
	    nSink = *(int *) msgCnt->getData();
	    delete msgCnt;
	    }
	}
    nAccreted -= nTotalSPH;
    
    sec1 = CkWallTimer();
    dsec = sec1 - sec;
    CkPrintf("Sinks Done (%d accreted) Calculated, Wallclock: %f secs\n\n",
	     nAccreted,dsec);
}

int BHDensitySmoothParams::isSmoothActive(GravityParticle *p) 
{
    if(p->rung < activeRung)
	return 0;
    return (TYPETest(p, TYPE_SINK));
    }

void BHDensitySmoothParams::initTreeParticle(GravityParticle *p)
{
    p->curlv()[1] = 0.0; /* total mass change */
    p->curlv()[0] = p->mass; /* initial mass */
    p->curlv()[2] = 0.0; /* flag if accreted */
    }

/* Cached Tree Active particles */
void BHDensitySmoothParams::initSmoothCache(GravityParticle *p)
{
    /*
     * Original particle curlv's is trashed (JW)
     */
    p->curlv()[1] = 0.0; /* total mass change */
    p->curlv()[0] = p->mass; /* initial mass */
    p->curlv()[2] = 0.0; /* flag if accreted */
    }

void BHDensitySmoothParams::combSmoothCache(GravityParticle *p1,
					    ExternalSmoothParticle *p2)
{
    if (p2->curlv[2] > p1->curlv()[2]) { /* only important if particle has been selected JMB 6/24/09 
					  * but necessary - if particle has 2 BH neighbors it can be
					  * overwritten */
	p1->curlv()[0] = p2->curlv[0];
	p1->curlv()[1] = p2->curlv[1]; /* total mass change */
	p1->curlv()[2] = p2->curlv[2]; 
	}

    /* this way the most recent one wins JMB 5/12/09 */
    /* don't sum up curlv[1], only let one BH eating be recorded. 
       curlv[1] will never be greater than curlv[0] this way */
}

/*
 * Calculate parameters for BH accretion for use in the BHSinkAccrete
 * function.  This is done separately to even competition between
 * neighboring Black Holes.
 * Gas particles are ID'd here to be eaten in BHSinkAccrete.
 */
void BHDensitySmoothParams::fcnSmooth(GravityParticle *p, int nSmooth,
				      pqSmoothNode *nnList)
{
	GravityParticle *q = NULL;

	double ih2,r2,rs,fDensity;
	double v[3],cs,fW,dv2,dv;
	double mdot, mdotEdd, mdotCurr, dm, dmq, dE, ifMass, dtEff;
	int i,iRung,naccreted, ieat, ivmin;
	double mdotsum, weat;
	double weight,wrs;
	double aFac, dCosmoDenFac,dCosmoVel2Fac;
	double dvmin,dvx,dvy,dvz;
	double ddvx,ddvy,ddvz,dvv; /* measure mindeltav */

	assert(p->rung >= s.iSinkCurrentRung);
	p->curlv()[1] = 0.0; /// XXX BH is a star
	naccreted = 0;
	aFac = a;
	dCosmoDenFac = aFac*aFac*aFac;
	dCosmoVel2Fac = aFac*aFac*aFac*aFac;

        mdotsum = 0.;
        weat = -1e37;

	ih2 = invH2(p);
	fDensity = 0.0; cs = 0;
	v[0] = 0; v[1] = 0; v[2] = 0;
	
	/* for mindv accretion method JMB 8/5/10 */
	if(s.bBHMindv == 1){
	  dvmin = FLT_MAX;
	  for (i=0;i<nSmooth;++i) {
	    ddvx = p->velocity[0] - nnList[i].p->velocity[0];
	    ddvy = p->velocity[1] - nnList[i].p->velocity[1];
	    ddvz = p->velocity[2] - nnList[i].p->velocity[2];
	    dvv = sqrt(ddvx*ddvx + ddvy*ddvy + ddvz*ddvz);
	    if (dvv < dvmin) {
	      dvmin=dvv;
	      ivmin = i;
	    }
	  }
	}

	for (i=0;i<nSmooth;++i) {
	    double fDist2 = nnList[i].fKey;
	    r2 = fDist2*ih2;
	    rs = KERNEL(r2);
	    q = nnList[i].p;
	    q->curlv()[0] = q->mass;
	    q->curlv()[1] = 0;
	    q->curlv()[2] = 0;

	    /* Sink should NOT be part of this list */
	    assert(TYPETest(q,TYPE_GAS));
	    fW = rs*q->mass;

	    if(s.bBHMindv == 1) weight = rs*pow(q->c()*q->c()+(dvmin*dvmin/dCosmoVel2Fac),-1.5)/dCosmoDenFac;
	    else {
	      dvx = p->velocity[0]-q->velocity[0];
	      dvy = p->velocity[1]-q->velocity[1];
	      dvz = p->velocity[2]-q->velocity[2];
	      weight = rs*pow(q->c()*q->c()+(dvx*dvx+dvy*dvy+dvz*dvz)/dCosmoVel2Fac,-1.5)/dCosmoDenFac; /* weight particles by mdot quantities */
	    /* cosmo factors put in 7/7/09  JMB */
	    }
	      if (weight > weat)  {   
		weat = weight;
		ieat = i; 
		wrs = rs;
	      }
	      mdotsum += weight*q->mass;


	    fDensity += fW;
	    v[0] += fW*q->velocity[0];
	    v[1] += fW*q->velocity[1];
	    v[2] += fW*q->velocity[2];
	    cs += fW*q->c();
	    }

	dv2 = 0;
	for (i=0;i<3;i++) {
	    dv = v[i]/fDensity-p->velocity[i];
	    dv2 += dv*dv;
	    }
	dv2 /= dCosmoVel2Fac;
	/*
	 * Store results in particle.
	 * XXX NB overloading "curlv" field of the BH particle.  I am
	 * assuming it is not used.
	 */
	p->c() = cs = cs/fDensity;
	p->fDensity = fDensity = M_1_PI*sqrt(ih2)*ih2*fDensity; 
	fDensity /= dCosmoDenFac;

	if(s.bBHMindv == 1)
	    CkPrintf("BHSink %d:  Time: %g Selected Particle Density: %g SPH: %g C_s: %g dv: %g dist: %g\n",
		     p->iOrder,dTime,nnList[ieat].p->fDensity, wrs,
		     nnList[ieat].p->c(),dvmin, sqrt(nnList[ieat].fKey));
	else CkPrintf("BHSink %d:  Time: %g Selected Particle Density: %g SPH: %g C_s: %g dv: %g dist: %g\n",
		      p->iOrder,dTime,nnList[ieat].p->fDensity,wrs,
		      nnList[ieat].p->c(),sqrt(dv2),sqrt(nnList[ieat].fKey));

	CkPrintf("BHSink %d:  Time: %g Mean Particle Density: %g C_s: %g dv: %g\n",p->iOrder,dTime,fDensity,cs,sqrt(dv2));


	mdot = mdotsum*s.dBHSinkAlpha*4*M_PI*p->mass*p->mass*M_1_PI*sqrt(ih2)*ih2; /* new mdot! */


	/* Eddington Limit Rate */
	mdotEdd = s.dBHSinkEddFactor*p->mass;
	CkPrintf("BHSink %d:  Time: %g mdot (BH): %g mdot (Edd): %g a: %g\n",
	       p->iOrder,dTime,mdot,mdotEdd, a);

	if (mdot > mdotEdd) mdot = mdotEdd;

	mdotCurr = p->divv() = mdot; /* store mdot in divv of sink XXX
					*/

        weight = -1e37;
	weat = -1e37;

	for (;;) {
	    double r2min = FLT_MAX;
	    q = NULL;
	    for (i=0;i<nSmooth;++i) {
		double fDist2 = nnList[i].fKey;
		r2 = fDist2;
		if(TYPETest(nnList[i].p,TYPE_DELETED)) continue;
		    /* don't accrete a deleted particle!  JMB 10/1/08 */
		rs = KERNEL(r2);
		fW = rs*nnList[i].p->mass;
		double fBall = nnList[i].p->fBall;
		if(r2 > 0.25*fBall*fBall) continue; 
		/* don't accrete gas that doesn't have the BH
		 * in its smoothing length  JMB 10/22/08 */
		if (nnList[i].p->rung < s.iSinkCurrentRung) continue; /* JMB 7/9/09 */

		if(s.bBHMindv == 1) weight = rs*pow(nnList[i].p->c()*nnList[i].p->c()+(dvmin*dvmin/dCosmoVel2Fac),-1.5)/dCosmoDenFac;
		else {
		  dvx = p->velocity[0]-nnList[i].p->velocity[0];
		  dvy = p->velocity[1]-nnList[i].p->velocity[1];
		  dvz = p->velocity[2]-nnList[i].p->velocity[2];
		  weight = rs*pow(nnList[i].p->c()*nnList[i].p->c()+(dvx*dvx+dvy*dvy+dvz*dvz)/dCosmoVel2Fac,-1.5)/dCosmoDenFac; /* weight particles by mdot quantities */
	    /* cosmo factors put in 7/7/09  JMB */
		}
	
		if (weight > weat && nnList[i].p->curlv()[2] != p->iOrder) {
		  weat = weight;
		  q = nnList[i].p;
		}
	    }
	    assert(q != p); /* shouldn't happen because p isn't in
			       tree */

	    if(q != NULL) {	    
	      /* Timestep for accretion is larger of sink and victim timestep */
	      iRung = q->rung; 
	      if (iRung > p->rung) iRung = p->rung;
	      dtEff = s.dSinkCurrentDelta*pow(0.5,iRung-s.iSinkCurrentRung);
	      /* If victim has unclosed kick -- don't actually take the mass
	       If sink has unclosed kick we shouldn't even be here!
	       When victim is active use his timestep if longer 
	       Statistically expect to get right effective mdot on average */
	      dmq = mdotCurr*dtEff;
	      }
	    else {
	      if(naccreted != 0) CkPrintf("BHSink %d:  only %i accreted, %g uneaten mass\n",p->iOrder,naccreted,mdotCurr*dtEff);
	      /* JMB 6/29/09 */
	      if(naccreted == 0) p->curlv()[1] = 0.0; /* should be anyway XXX*/
	      return;  	     	      
	    }
	    assert(q != NULL);
	    /* We have our victim */

	    if (q->rung >= s.iSinkCurrentRung) { /* moved up from below JMB 5/4/09 */
	                                             /* no longer needed, moved above */

	      if (dmq < q->mass) {
		mdotCurr = 0.0;
	      }
	      else {
		mdotCurr -= mdotCurr*(q->mass/dmq); /* need an additional victim */
		assert(mdotCurr >= 0.0);
		dmq = q->mass;
	      }
	      q->curlv()[2] = p->iOrder; /* flag for pre-used victim particles */
	      /* temporarily store mass lost -- to later check for double dipping */
	      q->curlv()[1] += dmq;
	      p->curlv()[1] += dmq; /* XXX */
	      naccreted += 1;
	      weat = -1e37; /* reset so other particles can be selected JMB 7/9/09 */

	      CkPrintf("BHSink %d:  Time %g %d dmq %g %g %g\n",p->iOrder,dTime,q->iOrder,dmq,q->curlv()[1],p->curlv()[1]);

		}
	    
	    if (mdotCurr == 0.0) break;
	    }
    }

void BHAccreteSmoothParams::initTreeParticle(GravityParticle *p1)
{
    /* Convert energy and metals to non-specific quantities (not per mass)
     * to make it easier to divvy up BH energy.  
     */
    
    if(TYPETest(p1, TYPE_GAS))  p1->u() *= p1->mass;    
    if(TYPETest(p1, TYPE_GAS))  p1->fNSN() = 0.0;
    /* fNSN will hold the place of FB energy because uPred is needed for
     *  calculations. it should be zero anyways.  JMB 10/5/09  */
    }

/* Cached Tree Active particles */
void BHAccreteSmoothParams::initSmoothCache(GravityParticle *p)
{
    if (TYPETest(p, iType)) p->u() = 0.0; /*added 6/10/08*/    
    if (TYPETest(p, iType)) p->fNSN() = 0.0; 
    }

static inline double max(double x, double y) 
{
    if (x > y) return x;
    else return y;
    }

void BHAccreteSmoothParams::combSmoothCache(GravityParticle *p1,
					     ExternalSmoothParticle *p2)
{
    if (!(TYPETest( p1, TYPE_DELETED )) &&
        TYPETest( p2, TYPE_DELETED ) ) {
	p1->mass = p2->mass;
	deleteParticle(p1);
	}
    else if (TYPETest(p1, iType)) {  // TreeActive
	/*
	 * See kludgery notice above: record eaten mass in original
	 * particle.
	 */
	double fEatenMass = p2->curlv[0] - p2->mass;
	p1->mass -= fEatenMass;
	if(p1->mass < p1->curlv()[0]*1e-3) {
	    /* This could happen if BHs on two
	       different processors are eating
	       gas from a third processor */
	    CkError("ERROR: Overeaten gas particle %d: %g %g\n",
		    p1->iOrder, p1->mass, fEatenMass);
	    if (!(TYPETest( p1, TYPE_DELETED ))) {
		deleteParticle( p1 );
		}
	    return;
	    }
	
	if(!(TYPETest(p1,TYPE_DELETED)))  { /*added 8/21/08 to check -0 mass particles  */
	    assert(p1->mass > 0.0); 
	    }
	p1->u() += p2->u; /*added 6/10/08 for BH blastwave FB*/
	p1->fNSN() += p2->fNSN;  /* (this is uPred) JMB 10/5/09  */
        p1->fTimeCoolIsOffUntil() = max( p1->fTimeCoolIsOffUntil(),
					 p2->fTimeCoolIsOffUntil );
	}
}

void BHAccreteSmoothParams::fcnSmooth(GravityParticle *p,int nSmooth,
				      pqSmoothNode *nnList)
{
	GravityParticle *q = NULL;

	double ih2,r2,rs,fDensity;
	double fW;
	double mdot, mdotCurr, dmAvg, dm, dmq, dE, ifMass, dtEff, dEwave;
	double fNorm,rstot,fNorm_u,fNorm_Pres,fAveDens;
        double fBHBlastRadius,fBHShutoffTime,fmind;
	int i,iRung,counter,imind,naccreted,ieat;
	double weat;
	double weight,fbweight; /* weight is for accretion, fbweight is for FB  */
	double aFac, dCosmoDenFac,dCosmoVel2Fac;
	double dvmin, dvx,dvy,dvz,dvv;
	int ivmin;

        weat = -1e37;
	aFac = a;
	dCosmoDenFac = aFac*aFac*aFac;
        dCosmoVel2Fac = aFac*aFac*aFac*aFac;

	mdot = p->divv();	
	if (p->curlv()[1] == 0.0) {
	    dtEff = s.dSinkCurrentDelta*pow(0.5,p->rung-s.iSinkCurrentRung);
	    dmAvg = mdot*dtEff;
	    CkPrintf("BHSink %d:  Delta: %g dm: 0 ( %g ) (victims on wrong step)\n",p->iOrder,dtEff,dmAvg);
	    return;
	    }

	mdotCurr = mdot;
	dm = 0;
	naccreted = 0;
	ih2 = invH2(p);

	for (i=0;i<nSmooth;++i){
	    double r2min = FLT_MAX;
	    q = NULL;
	    if (nnList[i].p->curlv()[2] == 0) continue;
	    if(TYPETest(nnList[i].p,TYPE_DELETED)) continue;
	    if (nnList[i].p->curlv()[2] == p->iOrder) q = nnList[i].p;
	    else continue; /* won't work if selected by another BH */
	    assert(q != NULL);
	    r2min = nnList[i].fKey;

	    /* Timestep for accretion is larger of sink and victim timestep */
	    iRung = q->rung;
	    if (iRung > p->rung) iRung = p->rung;
	    dtEff = s.dSinkCurrentDelta*pow(0.5,iRung-s.iSinkCurrentRung);
	      /* If victim has unclosed kick -- don't actually take the mass
	       If sink has unclosed kick we shouldn't even be here!
	       When victim is active use his timestep if longer 
	       Statistically expect to get right effective mdot on average */
	    dmq = mdotCurr*dtEff;

	    /* update mdotCurr */
	    if( fabs(dmq - q->curlv()[1])/dmq < 0.0001) mdotCurr = 0.0;
	    /* some floating point issues here JMB 7/9/09 */
	    else {
		mdotCurr -= q->curlv()[1]/dtEff;
		dmq = q->curlv()[1];
		}
	    CkAssert(mdotCurr >= 0.0);

	    CkPrintf("BHSink %d:  Time %g %d dmq %g %g %g\n",p->iOrder, dTime,
		     q->iOrder,dmq,q->curlv()[1],p->curlv()[1]);
	    if (q->rung >= s.iSinkCurrentRung) {
		ifMass = 1./(p->mass + dmq);
		/* to record angular momentum JMB 11/9/10 */
		CkPrintf("BHSink %d:  Gas: %d  dx: %g dy: %g dz: %g \n",
			 p->iOrder,q->iOrder,p->position[0]-q->position[0],
			 p->position[1]-q->position[1],
			 p->position[2]-q->position[2]);
		CkPrintf("BHSink %d:  Gas: %d  dvx: %g dvy: %g dvz: %g \n",
			 p->iOrder,q->iOrder,p->velocity[0]-q->velocity[0],
			 p->velocity[1]-q->velocity[1],
			 p->velocity[2]-q->velocity[2]);
		/* Adjust sink properties (conserving momentum etc...) */
		p->position[0] = ifMass*(p->mass*p->position[0]
					 +dmq*q->position[0]);
		p->position[1] = ifMass*(p->mass*p->position[1]
					 +dmq*q->position[1]);
		p->position[2] = ifMass*(p->mass*p->position[2]
					 +dmq*q->position[2]);
		p->velocity[0] = ifMass*(p->mass*p->velocity[0]
					 +dmq*q->velocity[0]);
		p->velocity[1] = ifMass*(p->mass*p->velocity[1]
					 +dmq*q->velocity[1]);
		p->velocity[2] = ifMass*(p->mass*p->velocity[2]
					 +dmq*q->velocity[2]);
		p->treeAcceleration[0] = ifMass*(p->mass*p->treeAcceleration[0]
						 +dmq*q->treeAcceleration[0]);
		p->treeAcceleration[1] = ifMass*(p->mass*p->treeAcceleration[1]
						 +dmq*q->treeAcceleration[1]);
		p->treeAcceleration[2] = ifMass*(p->mass*p->treeAcceleration[2]
						 +dmq*q->treeAcceleration[2]);
		p->fMetals() = ifMass*(p->mass*p->fMetals()+dmq*q->fMetals());
		p->mass += dmq;
		dm += dmq;
		q->mass -= dmq;
		/*assert(q->mass >= 0.0);*/
		naccreted += 1;  /* running tally of how many are accreted JMB 10/23/08 */
		CkPrintf("BHSink %d:  Time %g dist2: %d %g gas smooth: %g eatenmass %g \n",
			 p->iOrder,dTime,q->iOrder,r2min,q->fBall,dmq);
		if (q->mass <= 1e-3*dmq) { /* = added 8/21/08 */
		    q->mass = 0;
		    if(!(TYPETest(q,TYPE_DELETED))) deleteParticle(q);
		    /* Particles are getting deleted twice, which is messing
		      up the bookkeeping.  I think this will make them only be 
		      deleted once.  JMB 9/23/08*/
		    }

		}

	}

	if (mdotCurr == 0.0) goto dofeedback;
	else{	  
	  /* go looking for more particles, sharing must happen */
	    /****** SHARING CODE HERE ***********/ 

	  if(s.bBHMindv == 1){
	    dvmin = FLT_MAX;
	    for (i=0;i<nSmooth;++i) {
	      dvx = p->velocity[0] - nnList[i].p->velocity[0];
	      dvy = p->velocity[1] - nnList[i].p->velocity[1];
	      dvz = p->velocity[2] - nnList[i].p->velocity[2];
	      dvv = sqrt(dvx*dvx + dvy*dvy + dvz*dvz);
	      if (dvv < dvmin) {
		dvmin=dvv;
		ivmin = i;
	      }
	    }
	  }

	  for (;;) {
	    double r2min = FLT_MAX;
	    q = NULL;
	    for (i=0;i<nSmooth;++i) {
	      if(TYPETest(nnList[i].p,TYPE_DELETED)) continue;
	      if(nnList[i].p->curlv()[2] == p->iOrder) continue;
	      /* don't choose a pre-accreted particle 
	         but it can have been accreted by another BH */
	      if(nnList[i].p->rung < s.iSinkCurrentRung) continue;
	      /* has to be on the right timestep */
	      r2 = nnList[i].fKey;
	      if(r2 > 0.25*(nnList[i].p->fBall)*(nnList[i].p->fBall))
		  continue;
	      /* has to be nearby! */
	      rs = KERNEL(r2);

	      if(s.bBHMindv == 1)
		  weight = rs*pow(nnList[i].p->c()*nnList[i].p->c()
				  +(dvmin*dvmin/dCosmoVel2Fac),-1.5)/dCosmoDenFac;
	      else {
		dvx = p->velocity[0]-nnList[i].p->velocity[0];
		dvy = p->velocity[1]-nnList[i].p->velocity[1];
		dvz = p->velocity[2]-nnList[i].p->velocity[2];
		weight = rs*pow(nnList[i].p->c()*nnList[i].p->c()
				+(dvx*dvx+dvy*dvy+dvz*dvz)/dCosmoVel2Fac,-1.5)
		    /dCosmoDenFac; /* weight particles by mdot quantities */
		/* cosmo factors put in 7/7/09  JMB */
	      }	    	      
	      if (weight > weat) {
		r2min = r2; /* note r2min is not really the min r2 anymore */
		weat = weight;
		ieat = i;
		q = nnList[i].p;
	      } 
	    }

	    weat = -1e37; /* reset so other particles can be selected JMB 7/9/09 */

	    if(q == NULL) {
	      dtEff = s.dSinkCurrentDelta; /**pow(0.5,iRung-smf->iSinkCurrentRung);*/
	      CkPrintf("BHSink %d:  WARNING!! Not enough edible particles.  Time: %g Mass not eaten: %g \n",p->iOrder,dTime,mdotCurr*dtEff);
	      if(naccreted == 0) return;
	      else goto dofeedback;
	    }
	    else {
	      assert(q != NULL);
	      /* Timestep for accretion is larger of sink and victim timestep */
	      iRung = q->rung;
	      if (iRung > p->rung) iRung = p->rung;
	      dtEff = s.dSinkCurrentDelta
		  *pow(0.5,iRung-s.iSinkCurrentRung);
	      dmq = mdotCurr*dtEff;


	      if (dmq < q->curlv()[0]) mdotCurr = 0.0; /* Original mass in q->curlv[0] */	       
	      else {	       
		  mdotCurr -= mdotCurr*(q->curlv()[0]/dmq); /* need an additional victim */
		  dmq = q->curlv()[0];
	      }

	    
	    CkPrintf("BHSink %d:  Time: %g %d dmq %g %g %g sharing \n",
		     p->iOrder,dTime,q->iOrder,dmq,q->curlv()[1],
		     p->curlv()[1]);
	      if (q->rung >= s.iSinkCurrentRung) {
		ifMass = 1./(p->mass + dmq);
		/* to record angular momentum JMB 11/9/10 */
		CkPrintf("BHSink %d:  Gas: %d  dx: %g dy: %g dz: %g \n",
			 p->iOrder,q->iOrder,p->position[0]-q->position[0],
			 p->position[1]-q->position[1],
			 p->position[2]-q->position[2]);
		CkPrintf("BHSink %d:  Gas: %d  dvx: %g dvy: %g dvz: %g \n",
			 p->iOrder,q->iOrder,p->velocity[0]-q->velocity[0],
			 p->velocity[1]-q->velocity[1],
			 p->velocity[2]-q->velocity[2]);
		/* Adjust sink properties (conserving momentum etc...) */
		p->position[0] = ifMass*(p->mass*p->position[0]
					 +dmq*q->position[0]);
		p->position[1] = ifMass*(p->mass*p->position[1]
					 +dmq*q->position[1]);
		p->position[2] = ifMass*(p->mass*p->position[2]
					 +dmq*q->position[2]);
		p->velocity[0] = ifMass*(p->mass*p->velocity[0]
					 +dmq*q->velocity[0]);
		p->velocity[1] = ifMass*(p->mass*p->velocity[1]
					 +dmq*q->velocity[1]);
		p->velocity[2] = ifMass*(p->mass*p->velocity[2]
					 +dmq*q->velocity[2]);
		p->treeAcceleration[0] = ifMass*(p->mass*p->treeAcceleration[0]
						 +dmq*q->treeAcceleration[0]);
		p->treeAcceleration[1] = ifMass*(p->mass*p->treeAcceleration[1]
						 +dmq*q->treeAcceleration[1]);
		p->treeAcceleration[2] = ifMass*(p->mass*p->treeAcceleration[2]
						 +dmq*q->treeAcceleration[2]);
		p->fMetals() = ifMass*(p->mass*p->fMetals()+dmq*q->fMetals());
		p->mass += dmq;
		dm += dmq;
		q->mass -= dmq;		   
		naccreted += 1;  /* running tally of how many are accreted JMB 10/23/08 */
		CkPrintf("BHSink %d:  Time %g dist2 %d %g gas smooth: %g eatenmass %g\n",
			 p->iOrder,dTime,q->iOrder,r2min,q->fBall,dmq);
		if (q->mass <= 1e-3*dmq) { /* = added 8/21/08 */
		    q->mass = 0;
		    if(!(TYPETest(q,TYPE_DELETED))) deleteParticle(q);
		    /* Particles are getting deleted twice, which is messing
		      up the bookkeeping.  I think this will make them only be 
		      deleted once.  JMB 9/23/08*/
		    }
		}
	    if (mdotCurr == 0.0) break;

	    }
	  }
	}

	/********************************************************************/
        dofeedback:  

	if(s.dBHSinkFeedbackFactor != 0.0) {  /* allows FB to be turned off JMB 7/20/09 */
	  dE = s.dBHSinkFeedbackFactor*dm; /* dE based on actual mass eaten */
	  dEwave =  dE / s.dBHSinkFeedbackEff; /* dE for blastwave equations.  assume excess is radiated away. JMB 9/10/09 */

	  dtEff = s.dSinkCurrentDelta*pow(0.5,p->rung-s.iSinkCurrentRung);
	  dmAvg = mdot*dtEff;    
	  CkPrintf("BHSink %d:  Delta: %g Time: %g dm: %g dE %g\n",
		   p->iOrder,dtEff,dTime,dm,dE);

	  /* Recalculate Normalization */
	  ih2 = invH2(p);
	  fDensity = 0.0; 
	  rstot = 0.0;  
	  fNorm_u = 0.0;
	  fNorm_Pres = 0.0;
	  fAveDens = 0.0;
	  fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;

	  for (i=0;i<nSmooth;++i) {
	    r2 = nnList[i].fKey*ih2;
	    rs = KERNEL(r2);
	    fDensity += rs*nnList[i].p->mass;
	    q = nnList[i].p;
            fNorm_u += q->mass*rs;
            rs *= fNorm;
            fAveDens += q->mass*rs;
            fNorm_Pres += q->mass*(q->uPred()/q->curlv()[0])*rs;
	    assert(TYPETest(q, TYPE_GAS));
	  }

	  fNorm_Pres *= (gamma-1.0)/dCosmoDenFac;
	  fAveDens /= dCosmoDenFac;


	  /* eliminating the blast wave feedback here.  just give dE
	     to nearest particle and shut off cooling for deltat.
	     JMB 8/25/10 */
	
	  fmind = p->fBall*p->fBall;
	  imind = 0;
	  rstot = 0.0;  
	  fNorm_u = 0.0;

	  /* find the closest gas particle */
	  for (i=0;i<nSmooth;++i) {
	      q = nnList[i].p;
	      if(TYPETest(q,TYPE_DELETED)) continue;
	      if(q->mass == 0.0) continue;
	      if ( nnList[i].fKey < fmind ){imind = i; fmind = nnList[i].fKey;}
	      }

          /* give mass and energy to nearest gas particle. */

          r2 = nnList[imind].fKey*ih2;            
	  rs = KERNEL(r2);
	    /*
	     * N.B. This will be NEGATIVE, but that's OK since it will
	     * cancel out down below.
	     */
#ifdef VOLUMEFEEDBACK
	  fNorm_u = nnList[imind].p->mass/nnList[imind].p->fDensity*rs;
#else
	  fNorm_u = nnList[imind].p->mass*rs;
#endif
	 
	  assert(fNorm_u != 0.0);
	  fNorm_u = 1./fNorm_u;
	  counter=0;
	  q = nnList[imind].p;
	  fBHShutoffTime = q->dt;
	  if(s.bBHTurnOffCooling)
	      q->fTimeCoolIsOffUntil() = max(q->fTimeCoolIsOffUntil(),
					     dTime + fBHShutoffTime);
	  CkPrintf("BHSink %d: Time %g tCoolOffUntil %g \n",
		   p->iOrder,dTime,q->fTimeCoolIsOffUntil());
	  
	  q->fMassForm() = dTime;
	  /* track BHFB time in MassForm
		     JMB 11/19/08 */
	  counter++;  

	  /* Remember: We are dealing with total energy rate and total metal
	   * mass, not energy/gram or metals per gram.  
	   * q->fMass is in product to make units work for fNorm_u.
	   */
#ifdef VOLUMEFEEDBACK
	  fbweight = rs*fNorm_u*q->mass/q->fDensity;
#else
	  fbweight = rs*fNorm_u*q->mass;
#endif
	  q->u() += fbweight*dE;
	  q->fNSN() += fbweight*dE;
	  CkPrintf("BHSink %i:  FB Energy to %i dE %g\n",
		   p->iOrder,q->iOrder,fbweight*dE);	
	}
}

void BHAccreteSmoothParams::postTreeParticle(GravityParticle *p1)
{
    /* Convert energy  back to specific quantities (per mass)
       because we are done with our conservative calculations */
    if(TYPETest(p1, TYPE_GAS) && p1->mass != 0
       && !(TYPETest(p1,TYPE_DELETED))) {
	p1->u() /= p1->mass;  
	p1->fNSN() /= p1->mass;
	p1->uPred() += p1->fNSN(); /* now combine fNSN to uPred.  JMB
				      10/5/09  */
	}  
}

void BHIdentifySmoothParams::initSmoothCache(GravityParticle *p)
{
    p->fNSNtot() = 0.0;
    /*  fNSNtot is a placeholder for merging info*/
}

void BHIdentifySmoothParams::combSmoothCache(GravityParticle *p1,
					     ExternalSmoothParticle *p2)
{
    if(p2->fNSNtot > p1->fNSNtot()) p1->fNSNtot = p2->fNSNtot;
}

void BHIdentifySmoothParams::fcnSmooth(GravityParticle *p, int nSmooth,
				       pqSmoothNode *nnList)
{
  /* Identify BHs for merging JMB */

    GravityParticle *q = NULL;
    int i;
    double ifMass, deltaa, deltar, deltav;

    for (i=0;i<nSmooth;++i) {
	q = nnList[i].p;
	         
	if(nnList[i].fKey > 4.0*p->soft*p->soft || q->iOrder == p->iOrder)
	    continue; 
	/* are they close together? (within TWO softenings JMB 10/21/10*/
	/* don't include yourself!  JMB 12/11/08 */
	deltaa=sqrt(sqr(p->treeAcceleration[0] - q->treeAcceleration[0])
		    + sqr(p->treeAcceleration[1] - q->treeAcceleration[1])
		    + sqr(p->treeAcceleration[2] - q->treeAcceleration[2]));
	deltar=sqrt(sqr(p->position[0] - q->position[0])
		    + sqr(p->position[1] - q->position[1])
		    + sqr(p->position[2] - q->position[2]));
	deltav=sqrt(sqr(p->velocity[0] - q->velocity[0])
		    + sqr(p->velocity[1] - q->velocity[1])
		    + sqr(p->velocity[2] - q->velocity[2]));
	  
	if ( deltaa*deltar < 0.5*deltav*deltav ) continue;
	/* Selects other BH particles that are 
	 * within  the criteria 
	 * delta_a*delta_r < .5*delta_v^2  
	 AND 
	 * within the softening  */

	if(p->iOrder > q->iOrder) {
	    if(p->iOrder > q->fNSNtot()) {
		q->fNSNtot() = p->iOrder;
		CkPrintf("BHSink MergeID %d will be eaten by %d \n",
			 q->iOrder,p->iOrder);

		/* fNSNtot is the place holder for the iord of the 
		   eating black hole.  the victim remembers who eats
		   it for the merging step which is next.  JMB 4/30/09 */
		}
	    }
	}
}

void BHSinkMergeSmoothParams::fcnSmooth(GravityParticle *p,int nSmooth,
					pqSmoothNode *nnList)
{
  /* Makes BHs merge together, based on the criteria
     set in BHSinkIdentify above  */

    GravityParticle *q = NULL;
    double ifMass, deltaa, deltar, deltav;
    int i;
    /*  for kicks */
    const double A = 12000.0;
    const double B = -0.93;
    const double H = 7300.0;
    const double K = 60000.0;  /* BH kick parameters in km/s  */
    double vkick, vm, vpar, vkx, vky, vkz;
    double vperp, mratio, mfactor, cosine, spin1, spin2, angle1, angle2, xi;
    double xrand, yrand, zrand;

    for (i=0;i<nSmooth;++i) {
	q = nnList[i].p;
		
	/* Give mass to higher iOrd particle, delete lower. */
	if(q->fNSNtot() == p->iOrder) {
	    /* q has already been identified as a victim of p */	  
	    ifMass = 1./(p->mass + q->mass);
	    mratio = p->mass/q->mass;
	    if(mratio > 1.0) mratio = q->mass/p->mass;
	    /* Adjust sink properties (conserving momentum etc...) */
	    p->position[0] = ifMass*(p->mass*p->position[0]
				     +q->mass*q->position[0]);
	    p->position[1] = ifMass*(p->mass*p->position[1]
				     +q->mass*q->position[1]);
	    p->position[2] = ifMass*(p->mass*p->position[2]
				     +q->mass*q->position[2]);
	    p->velocity[0] = ifMass*(p->mass*p->velocity[0]
				     +q->mass*q->velocity[0]);
	    p->velocity[1] = ifMass*(p->mass*p->velocity[1]
				     +q->mass*q->velocity[1]);
	    p->velocity[2] = ifMass*(p->mass*p->velocity[2]
				     +q->mass*q->velocity[2]);
	    p->treeAcceleration[0] = ifMass*(p->mass*p->treeAcceleration[0]
					     +q->mass*q->treeAcceleration[0]);
	    p->treeAcceleration[1] = ifMass*(p->mass*p->treeAcceleration[1]
					     +q->mass*q->treeAcceleration[1]);
	    p->treeAcceleration[2] = ifMass*(p->mass*p->treeAcceleration[2]
					     +q->mass*q->treeAcceleration[2]);
	    p->mass += q->mass;


	      /**** Gravitational Recoil  ****/
	    if(s.bDoBHKick == 1
	       && (fabs(-1.0*p->fTimeForm() - dTime)
		   > s.dDeltaStarForm*SECONDSPERYEAR/smf->dSecUnit)) {
	      /* Turn off recoil if BH has just formed. Ejecting them
		 immediately is not helpful. JMB 8/5/09  */
	      mfactor = pow(mratio,2)/pow((1+mratio),5);
	      vm = A * (1.0-mratio) * mfactor * (1.0 + B*mratio/pow((1.0+mratio),2));
	      cosine = (rand()/((double) RAND_MAX )) * 2.*M_PI; /* angle */
	      spin1 = (rand()/((double) RAND_MAX )); /* spins of BHs 1 and 2 */
	      spin2 = (rand()/((double) RAND_MAX ));
	      angle1 = (rand()/((double) RAND_MAX )) * 2.*M_PI;  /* orientations of BH spins */
	      angle2 = (rand()/((double) RAND_MAX )) * 2.*M_PI;
	      xi = 1.5358897; /* 88.0 * 2. * 3.14159 / 360.0  */
	      /* in Campanelli et al xi = 88 degrees */

	      vperp = H * mfactor * (spin2*cos(angle2) - mratio*spin1*cos(angle1));  /* check sin cos! */
	      vpar = K * cos(cosine) * mfactor * (spin2*sin(angle2) - mratio*spin1*sin(angle1));

	      vkick = sqrt ( (vm + vperp*cos(xi))*(vm + vperp*cos(xi)) + (vperp*sin(xi))*(vperp*sin(xi)) + vpar*vpar);
	      vkick = vkick / smf->dKmPerSecUnit;
	      /* comoving coords are important JMB 5/15/09 */
	      vkick = vkick*a*a;
	    
	      /* random direction */

	      for(;;){
		xrand = (rand()/((double) RAND_MAX )) * 2.0 - 1.0;
		yrand = (rand()/((double) RAND_MAX )) * 2.0 - 1.0;
		zrand = (rand()/((double) RAND_MAX )) * 2.0 - 1.0;
		if(sqrt(xrand*xrand+yrand*yrand+zrand*zrand) < 1.0) break;
	      }
	      vkx = vkick * xrand / sqrt(xrand*xrand+yrand*yrand+zrand*zrand);
	      vky = vkick * yrand / sqrt(xrand*xrand+yrand*yrand+zrand*zrand);
	      vkz = vkick * zrand / sqrt(xrand*xrand+yrand*yrand+zrand*zrand);

	      p->velocity[0] += vkx;
	      p->velocity[1] += vky;
	      p->velocity[2] += vkz;
	    }
	    else vkick = 0.0;

	    CkPrintf("BHSink Merge %d eating %d  Time %g kick velocity %g mass ratio %g \n",
		     p->iOrder,q->iOrder,dTime,vkick,mratio);

	    if(q->fNSNtot() > 0 && !(TYPETest(q,TYPE_DELETED))) {
	      CkPrintf("BHSink Merge %d eaten Time %g \n", q->iOrder,dTime);
	      deleteParticle(q);
	    }
	  }
	}
}

int SinkAccreteTestSmoothParams::isSmoothActive(GravityParticle *p)
{
    if(p->rung < activeRung)
	return 0;
    return (TYPETest(p, TYPE_SINK));
    }

static inline double & fBindingEnergy(GravityParticle *a) {return a->curlv()[0];}
static inline double & fBindingEnergy(ExternalSmoothParticle *a)
{return a->curlv[0];}

/* Indicator for r,v,a update */
#define bRVAUpdate(_a)         ((_a)->curlv()[2])

void SinkAccreteTestSmoothParams::initSmoothCache(GravityParticle *p) 
{
    fBindingEnergy(p) = FLT_MAX;
    iOrderSink(p) = -1;
#ifdef SINKING
    if (TYPETest(p, TYPE_SINKING)) {
	iOrderSink(p) = p->iSinkingOnto;
	}
#endif
}

void SinkAccreteTestSmoothParams::combSmoothCache(GravityParticle *p1,
						  ExternalSmoothParticle *p2)
{
/* Particle p1 belongs to sink iOrderSink(p1) initially but will
   switch to iOrderSink(p2) if more bound to that sink */
#ifdef SINKING
    if (TYPETest(p1, TYPE_SINKING)) return;
#endif
    if (fBindingEnergy(p2) < fBindingEnergy(p1)) {
	fBindingEnergy(p1) = fBindingEnergy(p2);
	iOrderSink(p1) = iOrderSink_Ex(p2);
	}

    if (TYPETest( p2, TYPE_NEWSINKING)) TYPESet( p1, TYPE_NEWSINKING);
}

#ifdef SINKING
#define TRUEMASS(p__) (p__)->fTrueMass
#else
#define TRUEMASS(p__) (p__)->mass
#endif

void SinkAccreteTestSmoothParams::fcnSmooth(GravityParticle *p,int nSmooth,
					    pqSmoothNode *nnList)
{
	int i;
	double dSinkRadius2 = s.dSinkRadius*s.dSinkRadius, 
	       EBO,Eq,r2,dvx,dv2,ifMass;
	GravityParticle *q;

#ifdef SINKING
	assert(TYPETest(p, TYPE_SINK));
	assert(!TYPETest(p, TYPE_SINKING));
#endif
	/* G = 1 
	 p is sink particle
	 q is gas particle */
        if (s.dSinkBoundOrbitRadius > 0)
            EBO = -0.5*TRUEMASS(p)/s.dSinkBoundOrbitRadius;
        else
            EBO = FLT_MAX;

	for (i=0;i<nSmooth;++i) {
	    r2 = nnList[i].fKey;
	    if (r2 > 0 && r2 <= dSinkRadius2) {
		q = nnList[i].p;
#ifdef SINKING
		if (TYPETest(q, TYPE_SINKING)) {
		    FLOAT r0 = q->rSinking0Mag;
		    FLOAT r2 = r0 + q->vSinkingr0*(smf->dTime-q->fSinkingTime);
		    continue;
		    }
#endif
		if (TYPETest( q, TYPE_GAS )
		    && (q->rung >= s.iSinkCurrentRung
			|| TYPETest(q, TYPE_NEWSINKING))) {
		    dvx = p->velocity[0]-q->velocity[0];
		    dv2 = dvx*dvx;
		    dvx = p->velocity[1]-q->velocity[1];
		    dv2 += dvx*dvx;
		    dvx = p->velocity[2]-q->velocity[2];
		    dv2 += dvx*dvx;
		    Eq = -TRUEMASS(p)/sqrt(r2) + 0.5*dv2;
		    if (s.bSinkThermal) Eq+= q->u();
		    /* Being labelled NEWSINKING forces accretion to somebody even if unbound */
		    if (Eq < EBO || r2 < sqr(s.dSinkMustAccreteRadius)
			|| TYPETest(q, TYPE_NEWSINKING)) {
			if (Eq < fBindingEnergy(q)) {
			    fBindingEnergy(q) = Eq;
			    iOrderSink(q) = p->iOrder; /* Particle q belongs to sink p */
			    TYPESet(q, TYPE_NEWSINKING);
			    }
			}
		    }
		}   
	    }
}

int SinkingAverageSmoothParams::isSmoothActive(GravityParticle *p)
{
    return (TYPETest(p, TYPE_NEWSINKING));
    }

void SinkingAverageSmoothParams::fcnSmooth(GravityParticle *p,int nSmooth,
					   pqSmoothNode *nnList)
{
#ifdef SINKING
	int i,j;
	PARTICLE *q;
	double dv[3],wt,ih2,rs,r2,norm;

	assert(TYPETest(p,TYPE_NEWSINKING));
	wt = 0;
	for (j=0;j<3;j++) dv[j]=0;

	for (i=0;i<nSmooth;++i) {
	    q = nnList[i].p;
	    ih2 = 4.0/BALL2(p);
	    r2 = nnList[i].fDist2*ih2;
	    KERNEL(rs,r2);
	    rs*= q->fMass;
	    wt+=rs;
	    for (j=0;j<3;j++) dv[j]+=rs*(q->v[j]-p->v[j]);
	    }

	norm=1./wt;
        /* smoothed velocity */
	for (j=0;j<3;j++) {
	    p->vSinkingTang0Unit[j] = dv[j]*norm;
	    }
#endif
}

int SinkAccreteSmoothParams::isSmoothActive(GravityParticle *p)
{
    if(p->rung < activeRung)
	return 0;
    return (TYPETest(p, TYPE_SINK));
    }

void SinkAccreteSmoothParams::initSmoothCache(GravityParticle *p)
    { /* cached copies only */
#ifdef SINKING
    if (TYPETest( ((PARTICLE *) p), TYPE_SINKING )) {
	bRVAUpdate(p) = 0; /* Indicator for r,v,a update */
	}
#endif
    }

void SinkAccreteSmoothParams::combSmoothCache(GravityParticle *p1,
					      ExternalSmoothParticle *p2)
{
    if (!(TYPETest( p1, TYPE_DELETED )) &&
        TYPETest( p2, TYPE_DELETED ) ) {
	p1->mass = p2->mass;
	TYPEReset(p1,TYPE_SINKING);
	deleteParticle(p1);
	}
#ifdef SINKING
    else if (TYPETest(p2, TYPE_SINKING ) && !TYPETest( p1, TYPE_SINKING ))  {
	GravityParticle *p=p1,*q=p2;
	int j;
#ifdef SINKDBG
	printf("SINKING Sinking particle on other processor -- sinking properties combined %d %g\n",p->iOrder,p->fMetals);
#endif
	bRVAUpdate(p) = 1;
	assert(!TYPETest(p,TYPE_SINK));
	assert(!TYPETest(q,TYPE_SINK));
	p->fTrueMass = q->fTrueMass;
	TYPESet(p, TYPE_SINKING);
	p->fMetals = -1; /* HACK -- flag in output for sinking state */
	p->iSinkingOnto = q->iSinkingOnto;
	p->dt = q->dt;
	p->rung = q->iRung;
	p->fSinkingTime = q->fSinkingTime;

	p->vSinkingr0 = q->vSinkingr0;
        p->rSinking0Mag = q->rSinking0Mag;
	p->vSinkingTang0Mag = q->vSinkingTang0Mag;
	assert(q->vSinkingr0 < 0);
	for (j=0;j<3;j++) {
	    p->rSinking0Unit[j] = q->rSinking0Unit[j];
	    p->vSinkingTang0Unit[j] = q->vSinkingTang0Unit[j];
	    }
	}
    if (TYPETest( p1, TYPE_SINKING )) {
#ifdef SINKINGAVERAGE
	TYPEReset(p1,TYPE_NEWSINKING);
#endif
	if (bRVAUpdate(p2)) {
	    GravityParticle *p=p1,*q=p2;
	    p->r[0] = q->r[0];
	    p->r[1] = q->r[1];
	    p->r[2] = q->r[2];
	    p->v[0] = q->v[0];
	    p->v[1] = q->v[1];
	    p->v[2] = q->v[2];
	    p->vPred[0] = q->vPred[0];
	    p->vPred[1] = q->vPred[1];
	    p->vPred[2] = q->vPred[2];
	    p->a[0] = q->a[0];
	    p->a[1] = q->a[1];
	    p->a[2] = q->a[2];
#ifdef SINKDBG
	    if (((PARTICLE *)p1)->iOrder == 80) printf("FORCESINKACCRETECOMB %d with %d \n",-1,((PARTICLE *)p1)->iOrder);
#endif
	    }
	}
#endif
    }

/*#define DBGIORDER 2000000 */
#define DBGIORDER -1

void SinkAccreteSmoothParams::fcnSmooth(GravityParticle *p,int nSmooth,
					pqSmoothNode *nnList)
{
	int i,j,bEat=0;
	double ifMass;
	GravityParticle *q;
#ifdef SINKING
	double drSink[3];
	double r2;

	for (j=0;j<3;j++) drSink[j] = -p->position[j];
#endif

#ifdef STARSINK
	/* Convert to total angular momentum for additions */
	SINK_Lx(p) += TRUEMASS(p)*(p->position[1]*p->velocity[2]
				   - p->position[2]*p->velocity[1]); 
	SINK_Ly(p) += TRUEMASS(p)*(p->position[2]*p->velocity[0]
				   - p->position[0]*p->velocity[2]);
	SINK_Lz(p) += TRUEMASS(p)*(p->position[0]*p->velocity[1]
				   - p->position[1]*p->velocity[0]);
#endif
	for (i=0;i<nSmooth;++i) {
	    q = nnList[i].p;
	    if ( iOrderSink(q) == p->iOrder) {
#ifndef SINKING
		double mp,mq,rx,ry,rz,vx,vy,vz;
		mp = p->mass; mq = q->mass;
		rx = q->position[0]; ry = q->position[1]; rz = q->position[2];
		vx = q->velocity[0]; vy = q->velocity[1]; vz = q->velocity[2];
		ifMass = 1./(mp+mq);
		for (j=0;j<3;j++) {
		    p->position[j] = ifMass*(mp*p->position[j]
					     +mq*q->position[j]);
		    p->velocity[j] = ifMass*(mp*p->velocity[j]
					     +mq*q->velocity[j]);
		    p->treeAcceleration[j] = ifMass*(mp*p->treeAcceleration[j]
						     +mq*q->treeAcceleration[j]);
		    }
#ifdef STARSINK
		SINK_Lx(p) += mq*(ry*vz - rz*vy); /* add to total L */
  	        SINK_Ly(p) += mq*(rz*vx - rx*vz);
		SINK_Lz(p) += mq*(rx*vy - ry*vx);
#endif
		p->u() += ifMass*(mp*p->u()+mq*q->u());
		p->mass += mq;
		bEat = 1;
		assert(q->mass != 0);
		q->mass = 0;
		deleteParticle(q);
#endif
#ifdef SINKING
		assert(TYPETest(q, TYPE_GAS));
		assert(!TYPETest(q, TYPE_SINK));
		r2 = nnList[i].fKey;
		if (r2 < sinks.dSinkMustAccreteRadius*sinks.dSinkMustAccreteRadius || 
		    (TYPETest(q, TYPE_SINKING) && 
		     q->rSinking0Mag + q->vSinkingr0*(dTime-q->fSinkingTime) < smf->dSinkMustAccreteRadius)) {
#ifdef SINKDBG
		    if (q->iOrder == 55) printf("FORCESINKACCRETE0 %d with %d \n",p->iOrder,q->iOrder);
#endif
		    if (!TYPETest(q, TYPE_SINKING)) {
			double mp,mq,rx,ry,rz,vx,vy,vz;
			mp = p->fTrueMass; mq = q->fMass;
			rx = q->r[0]; ry = q->r[1]; rz = q->r[2];
			vx = q->v[0]; vy = q->v[1]; vz = q->v[2];
			/* Do a standard accrete -- no sinking stage */
			ifMass = 1./(mp+mq);
			for (j=0;j<3;j++) {
			    p->r[j] = ifMass*(mp*p->r[j]+mq*q->r[j]);
			    p->v[j] = ifMass*(mp*p->v[j]+mq*q->v[j]);
			    p->a[j] = ifMass*(mp*p->a[j]+mq*q->a[j]);
			    }
#ifdef STARSINK
			SINK_Lx(p) += mq*(ry*vz - rz*vy); /* add to total L */
			SINK_Ly(p) += mq*(rz*vx - rx*vz);
			SINK_Lz(p) += mq*(rx*vy - ry*vx);
#endif
			p->u += ifMass*(mp*p->u+mq*q->u);
			p->fMass += mq;
			p->fTrueMass += q->fMass;
			if (p->iOrder==DBGIORDER) printf("SINKING %d Direct Delete of %d +%g %g %g\n",p->iOrder,q->iOrder,q->fMass,p->fMass,p->fTrueMass);
			bEat = 1;
			assert(q->fMass != 0);
			q->fMass = 0;
			pkdDeleteParticle(smf->pkd, q);
			}
		    else if (q->iSinkingOnto == p->iOrder) {
			/* Already sinking -- now do final accrete */
			p->fMass += q->fMass; /* Note: p->fTrueMass NOT modified */
			p->iSinkingOnto--; /* One less sinking onto this sink */
			if (p->iOrder==DBGIORDER) printf("SINKING %d Final Delete of %d +%g %g %g %d\n",p->iOrder,q->iOrder,q->fMass,p->fMass,p->fTrueMass,p->iSinkingOnto);
			bEat = 1;
			assert(q->fMass != 0);
			q->fMass = 0;
			TYPEReset(q,TYPE_SINKING); /* Don't want it to participate later */
			pkdDeleteParticle(smf->pkd, q);
			}
		    }
		else if (!TYPETest(q, TYPE_SINKING)) {
		    FLOAT mp,mq,rx,ry,rz,vx,vy,vz;
		    mp = p->fTrueMass; mq = q->fMass;
		    rx = q->r[0]; ry = q->r[1]; rz = q->r[2];
		    vx = q->v[0]; vy = q->v[1]; vz = q->v[2];
		    /* Enter sinking stage */
		    assert(r2 < smf->dSinkRadius*smf->dSinkRadius);
		    /* If this code is active, ignore non-infalling particles */
/*		    vr = 0;
		    for (j=0;j<3;j++) {
			vr +=  (q->r[j]-p->r[j])*(q->v[j]-p->v[j]);
			}
			if (vr >= 0) continue;*/

		    /* All force, velocity, position info associated
		       with particle now belong to sink instead.
		       Note: If we want accurate estimate of L for the sink
		       we should adjust the sink L now */
		    /* Everything short of eating the particle 
		       If many additional particles are added sink position could vary wildly 
		       -- defer trajectory calculation until after */
		    ifMass = 1./(p->fTrueMass+q->fMass);
		    for (j=0;j<3;j++) {
			p->r[j] = ifMass*(mp*p->r[j]+q->fMass*q->r[j]);
			p->v[j] = ifMass*(mp*p->v[j]+q->fMass*q->v[j]);
			p->a[j] = ifMass*(mp*p->a[j]+q->fMass*q->a[j]);
			}
#ifdef STARSINK
		    SINK_Lx(p) += mq*(ry*vz - rz*vy); /* Add to total L */
		    SINK_Ly(p) += mq*(rz*vx - rx*vz);
		    SINK_Lz(p) += mq*(rx*vy - ry*vx);
#endif
		    p->u += ifMass*(mp*p->u+mq*q->u);
		    p->fTrueMass += q->fMass;
		    p->iSinkingOnto++; /* One more sinking onto this sink */
		    q->fTrueMass = 0; /* Sinking particles retain mass to act as ghost particles for forces g+SPH */
		    assert(!TYPETest(q,TYPE_SINK));
		    assert(TYPETest(q,TYPE_NEWSINKING));
		    TYPESet(q, TYPE_SINKING);
		    q->fMetals = -1; /* HACK -- flag for sinking state */
		    q->iSinkingOnto = p->iOrder;
		    q->dt = p->dt;
		    q->iRung = p->iRung;
		    q->fSinkingTime = smf->dTime;
		    q->vSinkingr0 = 1e30;
		    if (p->iOrder==DBGIORDER) printf("SINKING %d Add Sinking %d +%g %g %g %d\n",p->iOrder,q->iOrder,q->fMass,p->fMass,p->fTrueMass,p->iSinkingOnto);
		    bEat = 1;
		    }
#endif /*SINKING*/
		}
	    }
#ifdef STARSINK
	/* Store Internal angular momentum only as it is invariant under motions  */
	SINK_Lx(p) -= TRUEMASS(p)*(p->position[1]*p->velocity[2]
				   - p->position[2]*p->velocity[1]); 
	SINK_Ly(p) -= TRUEMASS(p)*(p->position[2]*p->velocity[0]
				   - p->position[0]*p->velocity[2]);
	SINK_Lz(p) -= TRUEMASS(p)*(p->position[0]*p->velocity[1]
				   - p->position[1]*p->velocity[0]);
#endif

#ifdef SINKING
	if (bEat) {
	    for (j=0;j<3;j++) drSink[j] += p->r[j]; /* If sink ate it moved: old sinking particles must follow */

	    for (i=0;i<nSmooth;++i) {
		q = nnList[i].p;
		if (TYPETest( q, TYPE_SINKING ) && q->iSinkingOnto == p->iOrder) {
		    if (TYPETest(q,TYPE_NEWSINKING)) {
			double r0, dx[3], dv[3], vr, vrraw, vrmin, vt, norm, dot;
			/* Initialize sinking trajectory */
			r0 = 0;
			for (j=0;j<3;j++) {
			    dx[j] = (q->r[j]-p->r[j]);
			    r0 += dx[j]*dx[j];
			    }
			r0 = sqrt(r0);
			norm = 1/r0;
			for (j=0;j<3;j++) dx[j] *= norm;

			vr = 0;
			for (j=0;j<3;j++) {
			    dv[j] = q->v[j]-p->v[j];  /* dv = raw relative velocity */
			    vr += dv[j]*dx[j];
			    } /* raw infall v */
#ifdef SINKINGAVERAGE
			vrraw = vr;
			    /* See if smoothed velocity is still infalling 
			       q->vSinkingTang0Unit temporarily contains mean of (q'-q) */
			vr = 0;
			for (j=0;j<3;j++) {
			    dv[j] = q->vSinkingTang0Unit[j]+q->v[j]-p->v[j]; /* dv now smoothed relative vel */
			    vr += dv[j]*dx[j];
			    }
			if (p->iOrder==DBGIORDER) printf("SINKING %d Add Sinking Trajectory vr %d %g %g %g\n",p->iOrder,q->iOrder,vr,vrraw,-q->c);
/*                    printf("Infall vr avg %g ",vr);*/
#endif
			vt = 0;
			for (j=0;j<3;j++) {
			    dv[j] -= vr*dx[j]; /*subtract off radial part*/
			    vt += dv[j]*dv[j];
			    }
			vt = sqrt(vt); /* Averaged tangential motion */
#ifdef SINKINGAVERAGE
			if (vr >= 0) vr=vrraw;
#endif
			vrmin = -q->c; /* Bondi-Hoyle minimum inflow speed -- perhaps not ideal for rotating case */
			if (vr > vrmin) vr = vrmin;
/*                    printf("min %g : %g   vt %g \n",-vrmin,vr,vt);*/
			
			q->rSinking0Mag = r0;
			if (r0 > smf->dSinkRadius) printf("WARNING: %i outside sink r=%g %g %g\n",q->iOrder,r0,sqrt(nnList[i].fDist2),smf->dSinkRadius);
			q->vSinkingr0 = vr;
			q->vSinkingTang0Mag = vt;
			norm = 1/(vt+1e-37);
			dot = 0;
			for (j=0;j<3;j++) {
			    q->rSinking0Unit[j] = dx[j];
			    q->vSinkingTang0Unit[j] = dv[j]*norm;
			    dot += q->rSinking0Unit[j]*q->vSinkingTang0Unit[j]; /* sanity check */
			    q->a[j] = p->a[j];
			    q->v[j] = p->v[j]; /* move with sink now */
			    }
			assert(fabs(dot) < 1e-4);
			}
		    else {
			assert(q->vSinkingr0 < 0);
			for (j=0;j<3;j++) {
			    q->r[j] += drSink[j];
			    q->a[j] = p->a[j];
			    q->v[j] = p->v[j];
			    }
			}
		    }
		}
	    }    
	
	/* Make sure sinking particles are still on their trajectories */
	for (i=0;i<nSmooth;++i) {
	    q = nnList[i].p;
	    if (TYPETest(q,TYPE_SINKING) && q->iSinkingOnto == p->iOrder) {
		FLOAT r0 = q->rSinking0Mag;
		FLOAT r2 = r0 + q->vSinkingr0*(smf->dTime-q->fSinkingTime);
		FLOAT thfac, th2, costh2, sinth2, dr2, sqr02;
		FLOAT dr[3];

		assert(q->vSinkingr0 < 0);
		if (r2 < 0.1*r0) r2 = 0.1*r0; /* HACK */
		thfac = q->vSinkingTang0Mag*2/(q->vSinkingr0);
		sqr02 = sqrt(r0/r2);
		th2 = thfac*(1-sqr02);
		costh2 = cos(th2);
		sinth2 = sin(th2);
		for (j=0;j<3;j++) {
		    q->r[j] = p->r[j]+r2*costh2*q->rSinking0Unit[j]+r2*sinth2*q->vSinkingTang0Unit[j];
		    q->vPred[j] = p->v[j]+q->vSinkingTang0Mag*sqr02*
			(-sinth2*q->rSinking0Unit[j]+costh2*q->vSinkingTang0Unit[j])
			+q->vSinkingr0*(costh2*q->rSinking0Unit[j]+sinth2*q->vSinkingTang0Unit[j]);
		    }
		bRVAUpdate(q) = 1; /* Indicator for r,v,a update */
#ifdef SINKDBG
		if (q->iOrder == 160460) printf("SINKINGFORCESHARE %d  %g %g %g  %g %g %g  %g %g %g (%g)\n",q->iOrder,smf->dTime,r2,sqrt(pow(q->r[0]-p->r[0],2.)+pow(q->r[1]-p->r[1],2.)+pow(q->r[2]-p->r[2],2.)), q->r[0]-p->r[0],q->r[1]-p->r[1],q->r[2]-p->r[2],q->vPred[0]-p->v[0],q->vPred[1]-p->v[1],q->vPred[2]-p->v[2],(q->vPred[0]-p->v[0])*(q->r[0]-p->r[0])+(q->vPred[1]-p->v[1])*(q->r[1]-p->r[1])+(q->vPred[2]-p->v[2])*(q->r[2]-p->r[2]));
#endif
		TYPEReset(q,TYPE_NEWSINKING);
		}

	    }
#endif
}
