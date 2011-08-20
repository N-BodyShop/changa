///
/// @brief initialize parameters for sinks
///

void Sink::AddParams(PRM prm) 
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

	dDeltaSink = param.dDelta;
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

///
/// Form sink particles; main routine
///
void Main::FormSinks(double dTime, double dDelta, int iKickRung)
/* For momentum conservation, only particles on iKickRung or higher may
   contribute mass to a sink.
   This call is timed for just after those particles have completed a full KDK
*/
    {
    int nStarOld,nSinkAdd;
    double sec,dsec;

    if (!param.bSinkForm) return;

    sec = CkWallTimer();
 
    nStarOld = msr->nStar;
    msr->param.dSinkCurrentDelta = dDelta;
    msr->param.iSinkCurrentRung = iKickRung;

    int bJeans = param.bSinkFormJeans;
    int bSimple = param.bSinkFormSimple;
    double dJConst2 = msr->param.dJeansConstant/msr->param.nJeans;
    dJConst2 *= dJConst2;
    double dDensityCut = param.dSinkFormDensity/param.dGmPerCcUnit;
    int bDensity = (dDensityCut > 0 ? 1 : 0);
    // msrResetType(msr, TYPE_ALL, TYPE_SMOOTHACTIVE);
    /* set sink candidates as SMOOTHACTIVE */
    CkReductionMsg *msg;
    treeProxy.formSinks(bJeans, dJConst2, bDensity, dDensityCut, dTime,
			iKickRung, bSimple,
			CkCallbackResumeThread((void*&)msg));
    int nCandidates = *(int *)msg->getData();
    delete msg;
    
    nSinkAdd = 0;
    if (nCandidates) {
	if (param.bSinkFormSimple) {
	    nSinkAdd = nCandidates;
	    nSink += nSinkAdd;
	    } 
	else {
	    int iPhase = 0; // Keeps track of node cache use
	/* Note: Only gas particles are combined into sinks */
	    SinkFormTestSmoothParams pSFT(TYPE_GAS);
	    // msrSmooth(msr, dTime, SMX_SINKFORMTEST, 1); /* should exclude some who fail this test? */
	    double dfBall2OverSoft2 = 4.0*param.dhMinOverSoft*param.dhMinOverSoft;
	    treeProxy.startIterationSmooth(&pSFT, 1, dfBall2OverSoft2,
					   CkCallbackResumeThread());
	    iPhase++;

	    SinkFormSmoothParams pSF(TYPE_GAS);
	    // msrSmooth(msr, dTime, SMX_SINKFORM, 1);
	    double dfBall2OverSoft2 = 4.0*param.dhMinOverSoft*param.dhMinOverSoft;
	    treeProxy.startIterationSmooth(&pDSFB, 1, dfBall2OverSoft2,
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
			  int bSimple)
{
    int i;
    PARTICLE *p;
    int n = myNumParticles;
    double Jval, Jvalmin = FLT_MAX;
    int nCandidates = 0;
    
    for(i = 1; i <= n; ++i) {
	GravityParticle *p = &myParticles[i];
    
        if(p->isGas() && p->iRung >= iKickRung) {
#ifdef SINKING
	    if (TYPETest( p, TYPE_SINKING )) continue;
#endif
/* Jeans Mass compared to nJeans particle masses */
	    Jval =  dJConst2*(p->c*p->c*p->c*p->c*p->c*p->c)
		/(p->fMass*p->fMass*p->fDensity);
	    if (Jval < Jvalmin) Jvalmin = Jval;
	    if ((bJeans && Jval < 1) ||
		(bDensity && p->fDensity >= dDensityCut)) {
		nCandidates++;
		if (bSimple) {
		    TYPESet(p, TYPE_SINK); /* Is now a sink! */
		    p->fMetals = -dTime;
		    }
		else {
		    TYPESet(p, TYPE_SMOOTHACTIVE); /* Is now a candidate */
		    }
		}
	    }
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
    formSinks(dTime, dDelta, iKickRung ); 
   
    if (msr->param.bDoSinks == 0) return;
    if (msr->nSink == 0) return;
    if (msr->param.bBHSink && dDelta <= 0.0) return;

    sec = msrTime();
    dMass = msrMassCheck(msr, -1, "Accrete onto Sinks: Initial Value");

    /* Note: Only gas particles are accreted by sinks */
/*    printf("Tree: %d %d\n",msr->iTreeType,MSR_TREE_DENSITY);*/
    if (msr->iTreeType != MSR_TREE_DENSITY) {
	    msrActiveType(msr,TYPE_GAS,TYPE_TREEACTIVE);
	    msrBuildTree(msr,1,-1.0,1);  /* bTreeActive */
	    }

    msrResetType(msr,TYPE_SINK,TYPE_SMOOTHDONE);
    msrActiveTypeRung(msr,TYPE_SINK,TYPE_ACTIVE|TYPE_SMOOTHACTIVE,iKickRung,1);

    nAccreted = msr->nGas;
    if (msr->nActive > 0) {
	msr->param.dSinkCurrentDelta = dDelta;
	msr->param.iSinkCurrentRung = iKickRung;

	if (msr->param.bBHSink) {
	    /* Smooth Bondi-Hoyle Accretion: radius set by nSmooth */
	    msrSmooth(msr, dTime, SMX_BHDENSITY, 1);
	    msrResetType(msr,TYPE_SINK,TYPE_SMOOTHDONE);
	    msrSmooth(msr, dTime, SMX_BHSINKACCRETE,1);
	    /* build new tree of BHs for merging JMB 11/14/08  */
	    msrActiveType(msr,TYPE_SINK,TYPE_TREEACTIVE);
	    if (msr->nTreeActive > 1) { /* no need to merge if there is only one! */
	      msrBuildTree(msr,1,-1.0,1);  /* bTreeActive */
	      msrResetType(msr,TYPE_SINK,TYPE_SMOOTHDONE);
	      msrActiveTypeRung(msr,TYPE_SINK,TYPE_ACTIVE|TYPE_SMOOTHACTIVE,iKickRung,1);
	      /* need to change nSmooth to number of BHs.  */
	      nSmoothTemp = msr->param.nSmooth;
	      msr->param.nSmooth = msr->nTreeActive;
	      msrSmooth(msr,dTime, SMX_BHSINKIDENTIFY,1);
      	      msrResetType(msr,TYPE_SINK,TYPE_SMOOTHDONE);
	      msrSmooth(msr,dTime, SMX_BHSINKMERGE,0);
	      /* now change it back to what it was before JMB 12/10/08  */
	      msr->param.nSmooth = nSmoothTemp;
	    }
	    }
	else {
	    /* Fixed Radius Accretion: particle by particle (cf. Bate) */
	    msrSmooth(msr, dTime, SMX_SINKACCRETETEST,1);
	    msrResetType(msr,TYPE_SINK,TYPE_SMOOTHDONE);
#ifdef SINKINGAVERAGE
	    msrActiveType(msr,TYPE_NEWSINKING, TYPE_SMOOTHACTIVE);
	    msrSmooth(msr, dTime, SMX_SINKINGAVERAGE,0);
#endif
	    msrActiveTypeRung(msr,TYPE_SINK,TYPE_ACTIVE|TYPE_SMOOTHACTIVE,iKickRung,1);
	    msrSmooth(msr, dTime, SMX_SINKACCRETE,1);
	    }
	
	msrMassCheck(msr, dMass, "Accrete onto Sinks: before particle adjustment");
	
	msrAddDelParticles(msr);
	if (msr->param.bBHSink) { /* reset nSinks in case of merger */
	    msrActiveType(msr,TYPE_SINK,TYPE_TREEACTIVE);
	    msr->nSink = msr->nTreeActive;
	    }
	msrMassCheck(msr, dMass, "Accrete onto Sinks: after particle adjustment");
	}
    nAccreted -= msr->nGas;
    
    sec1 = msrTime();
    dsec = sec1 - sec;
    printf("Sinks Done (%d accreted) Calculated, Wallclock: %f secs\n\n",nAccreted,dsec);
    LOGTIMINGUPDATE( dsec, TIMING_Sink );
}
