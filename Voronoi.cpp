#ifdef VORONOI
#include "ParallelGravity.h"
#include "DataManager.h"
#include "smooth.h"
#include "physconst.h"
#include "Sph.h"
#include "Neighbor.h"
#include "HydroParticle.h"
#include "VoronoiSmoothParams.h"
#include "parameters.h"
#include "Reductions.h"
#include "UserGravity.h"
#include "UserSources.h"
#include "FixedBC.h"
#include "UserProblemGenerator.h"
#ifdef RADIATION
#include "UserRadiationSources.h"
#include "RadiationSource.h"
#include "Opacity.h"
#include "SimpleOpacity.h"
#include "DustOpacity.h"
#endif

void Main::setupVoroICs() {
  bool arrayFileExists(const std::string filename, const int64_t count) ;
  if( param.bDoGas) {
    double dvFac;
    if(param.csm->bComove) {
      double dOutTime = csmTime2Exp(param.csm, dTime);
      dvFac = 1.0/(dOutTime*dOutTime);
    }
    else {
      double dOutTime = dTime;
      dvFac = 1.0;
    }
    if(arrayFileExists(basefilename + ".rho", nTotalParticles) &&
       arrayFileExists(basefilename + ".ie", nTotalParticles) &&
       arrayFileExists(basefilename + ".vx", nTotalParticles) &&
       arrayFileExists(basefilename + ".vy", nTotalParticles) &&
       arrayFileExists(basefilename + ".vz", nTotalParticles)) {
      ckout << "Reading Primitives" << endl;
      RhoOutputParams rhoOut(basefilename, param.iBinaryOut, 0.0);
      IeOutputParams iEOut(basefilename, param.iBinaryOut, 0.0);
      VxOutputParams vxOut(basefilename, param.iBinaryOut, 0.0, dvFac);
      VyOutputParams vyOut(basefilename, param.iBinaryOut, 0.0, dvFac);
      VzOutputParams vzOut(basefilename, param.iBinaryOut, 0.0, dvFac);

      useVoronoiICs = true;
#ifdef MHD
      BxOutputParams BxOut(basefilename, param.iBinaryOut, 0.0);
      ByOutputParams ByOut(basefilename, param.iBinaryOut, 0.0);
      BzOutputParams BzOut(basefilename, param.iBinaryOut, 0.0);
      AxOutputParams AxOut(basefilename, param.iBinaryOut, 0.0);
      AyOutputParams AyOut(basefilename, param.iBinaryOut, 0.0);
      AzOutputParams AzOut(basefilename, param.iBinaryOut, 0.0);
      treeProxy.readTipsyArray(AxOut, CkCallbackResumeThread());
      treeProxy.readTipsyArray(AyOut, CkCallbackResumeThread());
      treeProxy.readTipsyArray(AzOut, CkCallbackResumeThread());
      treeProxy.readTipsyArray(BxOut, CkCallbackResumeThread());
      treeProxy.readTipsyArray(ByOut, CkCallbackResumeThread());
      treeProxy.readTipsyArray(BzOut, CkCallbackResumeThread());
#endif

      treeProxy.readTipsyArray(rhoOut, CkCallbackResumeThread());
      treeProxy.readTipsyArray(iEOut, CkCallbackResumeThread());
      treeProxy.readTipsyArray(vxOut, CkCallbackResumeThread());
      treeProxy.readTipsyArray(vyOut, CkCallbackResumeThread());
      treeProxy.readTipsyArray(vzOut, CkCallbackResumeThread());
      ckout << "Finished Reading Primitives " << endl;
      // TODO: Modify OutputParams to allow for setting Vector Values for reading in vectors
    }
  }
}

void Main::initVoro() {
  if(param.bDoGas) {

    dMProxy.setupMovingMesh( param, CkCallbackResumeThread());

    //if( !bVoroICs) {
      // The following smooths all GAS, and also marks neighbors of
      // actives, and those who have actives as neighbors.
      ckout << "Calculating densities/divv ...";
      DenDvDxSmoothParams pDen(TYPE_GAS, 0, param.csm, dTime, 0,
			       param.bConstantDiffusion, 1, bHaveAlpha,
                               param.dConstAlphaMax);
      double startTime = CkWallTimer();
      double dfBall2OverSoft2 = 4.0*param.dhMinOverSoft*param.dhMinOverSoft;
      treeProxy.startSmooth(&pDen, 1, param.nSmooth, dfBall2OverSoft2,
                            CkCallbackResumeThread());
      ckout << " took " << (CkWallTimer() - startTime) << " seconds."
            << endl;
      if(verbosity > 1 && !param.bConcurrentSph)
          memoryStatsCache();

      double dTuFac = param.dGasConst/(param.dConstGamma-1)
    	    /param.dMeanMolWeight;
      double z = 1.0/csmTime2Exp(param.csm, dTime) - 1.0;

      if(param.bGasCooling) {
        // Update cooling on the datamanager
        dMProxy.CoolingSetTime(z, dTime, CkCallbackResumeThread());
        if(!bIsRestarting)  // Energy is already OK from checkpoint.
          treeProxy.InitEnergy(dTuFac, z, dTime, (param.dConstGamma-1), CkCallbackResumeThread());
      }
    //}
    ckout << "Starting doVoroStart" << endl;
    Main::doVoroStart(0, 0., true);
    int nMaxRungs = adjust(0);
    treeProxy.setupHydroParticleRungs( param.bFastGas ? -1 : nMaxRungs, CkCallbackResumeThread());
  }
}

void
Main::doVoroStart(int activeRung, double dt, bool initialize)
{
  bool activeOnly = param.bFastGas;// && nActiveSPH < nTotalSPH*param.dFracFastGas;
  if( initialize) activeOnly = false;
  if(param.bDoGas) {
    double startTime = CkWallTimer();

    if (initialize) ckout << "Initializing and ";

    // reset the hydroparticles
    treeProxy.resetHydroParticle( CkCallbackResumeThread());

    // build the voronoi cells.
#ifdef MHD
    int mode = initialize ? VoronoiSmoothParams::INITIALIZE : VoronoiSmoothParams::MAP_A_TO_B;
    if( initialize){
      double convFactor = 1./sqrt(param.dErgPerGmUnit*param.dGmPerCcUnit);
      treeProxy.updateHydroParticleBfield( param.dBx*convFactor, param.dBy*convFactor, param.dBz*convFactor,
                                           CkCallbackResumeThread());
    }
#else
    int mode = initialize ? VoronoiSmoothParams::INITIALIZE : VoronoiSmoothParams::BUILD_MESH;
#endif //MHD
    int numIterations = activeOnly && !initialize ? 2 : 1;
    int iterations[2] = {numIterations == 2 ? VSP_ACTIVEONLY : VSP_EVERYTHING, VSP_NEIGHBORONLY};
    for( int iter = 0; iter < numIterations; iter++) {

#ifdef MHD
      ckout << "Mapping A to B ";
#else
      ckout << "Computing gradient and half time predictions on Rung:" << activeRung << "on iteration: " << iter;
#endif //MHD

      VoronoiSmoothParams voro(TYPE_GAS,
                               activeRung,
                               param.csm,
                               dTime,
                               iterations[iter],
                               dt,
                               mode,
                               false);
    
      if( useVoronoiICs ) {
        voro.useVoronoiICs();
      }
      if( UserProblemGenerator::instance != 0) {
        voro.useProblemGenerator();
      }

      double dfBall2OverSoft2 = 4.0 * param.dhMinOverSoft * param.dhMinOverSoft;
      //double dfBall2OverSoft2 = 100.0 * param.dhMinOverSoft * param.dhMinOverSoft;
      if( initialize) {
        treeProxy.startSmooth(&voro, 1, param.nSmooth, dfBall2OverSoft2,
                              CkCallbackResumeThread());
      }
      else {
        treeProxy.startReSmooth(&voro, CkCallbackResumeThread());
      }
      // treeProxy.startReSmooth(&voro, CkCallbackResumeThread());
      ckout << " took " << (CkWallTimer() - startTime) << " seconds."
            << endl;

      startTime = CkWallTimer();
      ckout << "Clean up Voronoi on " << iter;
      treeProxy.setFBallHydroParticle( activeRung, CkCallbackResumeThread());
      VoronoiSmoothParams voroClean(TYPE_GAS,
                               activeRung,
                               param.csm,
                               dTime,
                               iterations[iter],
                               dt,
                               mode,
                               true);
      if( useVoronoiICs ) {
        voroClean.useVoronoiICs();
      }
      if( UserProblemGenerator::instance != 0) {
        voroClean.useProblemGenerator();
      }

      treeProxy.startReSmooth(&voroClean, CkCallbackResumeThread());

      ckout << " took " << (CkWallTimer() - startTime) << " seconds."
            << endl;

      getGlobalInfo(activeRung);

    }

    ckout << "call checkHydroParticle" << endl;
    treeProxy.checkHydroParticle( CkCallbackResumeThread());

    //ckout << "call infoFromHydroParticle" << endl;
    //treeProxy.infoFromHydroParticle( dt, param.dVoroKick, CkCallbackResumeThread());
    
    treeProxy.setFBallHydroParticle( activeRung, CkCallbackResumeThread());
    getGlobalInfo( activeRung);
    
    treeProxy.ballMax(activeRung, 1.0 + param.ddHonHLimit,
                      CkCallbackResumeThread());
    ckout << "finished ballMax" << endl;

#ifdef MHD
    if( !initialize) {
      // We need to build the mesh
      ckout << "Computing gradient and half time predictions";
      startTime = CkWallTimer();
      VoronoiSmoothParams voroBuild(TYPE_GAS, activeRung, param.csm, dTime, 
                              activeOnly ? VSP_ACTIVE_AND_NEIGHBOR : VSP_EVERYTHING,
                              dt, VoronoiSmoothParams::BUILD_MESH, false);
      treeProxy.startReSmooth(&voroBuild, CkCallbackResumeThread());
      ckout << " took " << (CkWallTimer() - startTime) << " seconds."
            << endl;
    }
#endif//MHD
  }
}

void
Main::doVoroSolve(int activeRung, double dt)
{
  bool activeOnly = param.bFastGas;// && nActiveSPH < nTotalSPH*param.dFracFastGas;
  if(param.bDoGas) {
  
    // reset the hydroparticles
    treeProxy.resetHydroParticle( CkCallbackResumeThread());

    int mode = VoronoiSmoothParams::HALF_STEP;
    //int mode = VoronoiSmoothParams::BUILD_MESH;
    int numIterations = activeOnly ? 2 : 1;
    int iterations[2] = {numIterations == 2 ? VSP_ACTIVEONLY : VSP_EVERYTHING, VSP_NEIGHBORONLY};
    for( int iter = 0; iter < numIterations; iter++) {
      double startTime = CkWallTimer();

      ckout << "Performing half step gradients on " << iter;
      VoronoiSmoothParams voro(TYPE_GAS,
                               activeRung,
                               param.csm,
                               dTime,
                               iterations[iter],
                               dt,
                               mode,
                               false);
      treeProxy.startReSmooth(&voro, CkCallbackResumeThread());
      // treeProxy.startReSmooth(&voro, CkCallbackResumeThread());
      ckout << " took " << (CkWallTimer() - startTime) << " seconds."
            << endl;
      //getGlobalInfo( activeRung);
      startTime = CkWallTimer();
      ckout << "Clean up Voronoi Half Step on " << iter;
      //treeProxy.setFBallHydroParticle( activeRung, CkCallbackResumeThread());

      VoronoiSmoothParams voroClean(TYPE_GAS,
                               activeRung,
                               param.csm,
                               dTime,
                               iterations[iter],
                               dt,
                               mode,
                               true);
      treeProxy.startReSmooth(&voroClean, CkCallbackResumeThread());

      ckout << " took " << (CkWallTimer() - startTime) << " seconds."
            << endl;

      //ckout << "call checkHydroParticle" << endl;
      //treeProxy.checkHydroParticle( CkCallbackResumeThread());
    
      treeProxy.setFBallHydroParticle( activeRung, CkCallbackResumeThread());
      getGlobalInfo( activeRung);
    
      treeProxy.ballMax(activeRung, 1.0 + param.ddHonHLimit,
                        CkCallbackResumeThread());
    }
    ckout << "Solving Riemann ";

    // proceed to solve Riemann problem
    double startTime = CkWallTimer();
    VoronoiSmoothParams voroSolve(TYPE_GAS, activeRung, param.csm, dTime, 
                            activeOnly ? VSP_ACTIVE_AND_NEIGHBOR : VSP_EVERYTHING,
                            dt, VoronoiSmoothParams::FULL_STEP, false);
    treeProxy.startReSmooth(&voroSolve, CkCallbackResumeThread());
    ckout << " took " << (CkWallTimer() - startTime) << " seconds."
          << endl;

    // put back the corrected acceleration
    //getGlobalInfo( activeRung);
#ifdef MHD
    double bx, by, bz;
    ckout << "Computing global B-field ";
    startTime = CkWallTimer();
    getGlobalBfield( dt, bx, by, bz);
    //treeProxy.updateHydroParticleBfield( bx, by, bz, CkCallbackResumeThread());
 
    ckout << " took " << (CkWallTimer() - startTime) << " seconds."
          << endl;

#endif
      startTime = CkWallTimer();
      if( false) {
        ckout << "Checking Solution ";
        VoronoiSmoothParams voroCheck(TYPE_GAS, activeRung, param.csm, dTime, VSP_EVERYTHING,
                                      dt, VoronoiSmoothParams::CHECK_SOLUTION, true);
        treeProxy.startReSmooth(&voroCheck, CkCallbackResumeThread());
        ckout << " took " << (CkWallTimer() - startTime) << " seconds."
              << endl;
      }
    //treeProxy.infoFromHydroParticle( dt, param.dVoroKick, CkCallbackResumeThread());



  }
}


void
Main::doVoroEnd(int activeRung, double dt)
{
  if(param.bDoGas) {
    bool activeOnly = param.bFastGas;//&& nActiveSPH < nTotalSPH*param.dFracFastGas;
    treeProxy.infoFromHydroParticle( activeRung, param.dVoroKick, CkCallbackResumeThread(), activeOnly);
    ckout << "call checkHydroParticle after update" << endl;
    treeProxy.checkHydroParticle( CkCallbackResumeThread());
  }
}

void Main::getGlobalInfo( int activeRung) {
   CkReductionMsg *msg;
   treeProxy.getGlobalInfo(CkCallbackResumeThread((void*&)msg));
   double *dGlobalInfo = (double *) msg->getData();

   printf("total mass = %7.5e %7.5e\n", dGlobalInfo[0], dGlobalInfo[1]);
   printf("total KE = %7.5e\n", dGlobalInfo[5]);
#ifdef RADIATION
   printf("total radiation = %7.5e\n", dGlobalInfo[6]);
#endif
   CkReductionMsg *msg2;
   treeProxy.getGlobalStatistics( activeRung, CkCallbackResumeThread((void*&)msg2));
   int *iStats = (int *) msg2->getData();

   printf("validCells half= %d %d\n", iStats[0], iStats[1]);
   printf("validCells full= %d %d\n", iStats[2], iStats[3]);
   printf("neighbors %d full %d half %d\n", iStats[4], iStats[5], iStats[6]);

   CkReductionMsg *msg3;
   treeProxy.getGlobalDistributions(CkCallbackResumeThread((void*&)msg3));
   int *iDist = (int *) msg3->getData();

   printf("max #search = %d \n", iDist[0]);
}

#ifdef MHD
void Main::getGlobalBfield( double dt, double &bx, double &by, double &bz) {
   CkReductionMsg *msg;
   treeProxy.getGlobalBfield( dt, CkCallbackResumeThread((void*&)msg));
   double *dGlobalBfield = (double *) msg->getData();
   double vol = dGlobalBfield[6];
    
   double convFactor = sqrt(param.dErgPerGmUnit*param.dGmPerCcUnit);
   bx = dGlobalBfield[0]/vol*convFactor;
   by = dGlobalBfield[1]/vol*convFactor;
   bz = dGlobalBfield[2]/vol*convFactor;
   
   double dbx = dGlobalBfield[3]/vol*convFactor;
   double dby = dGlobalBfield[4]/vol*convFactor;
   double dbz = dGlobalBfield[5]/vol*convFactor;

   double totalKE = dGlobalBfield[7];
   double totalMass = dGlobalBfield[8]; 
   double v = sqrt(2.*totalKE/totalMass);
   double vA = sqrt(dGlobalBfield[9]/totalMass);
   double cs = sqrt(dGlobalBfield[10]/totalMass)*(3.0857e21*param.dKpcUnit)/param.dSecUnit;
   printf( "Global Bfield: %5.3e %5.3e %5.3e %5.3e\n", bx, by, bz, vol);
   printf( "Global Bfield Change: %5.3e %5.3e %5.3e \n", dbx, dby, dbz);
   printf( "<v>: %5.3e %5.3e %5.3e \n", v*(3.0857e21*param.dKpcUnit)/param.dSecUnit, totalKE, totalMass);
   printf( "<vA>, <cs>: %5.3e %5.3e \n", vA*(3.0857e21*param.dKpcUnit)/param.dSecUnit, cs);

}
#endif //MHD

void TreePiece::resetHydroParticle( const CkCallback &cb) {
  GravityParticle *p;
#ifndef COOLING_NONE
  HydroUtils::setupEOS( dm->Cool);
#endif
  for(int i=1; i<= myNumParticles; ++i) {
    p = &myParticles[i];
    TYPEReset( p, TYPE_NbrOfACTIVE);

    if (TYPETest(p, TYPE_GAS)) {
      p->hydroParticle().reset();
      p->hydroParticle().iOrder = p->iOrder;
      p->hydroParticle().updateParticlePosition(p->position.x,p->position.y,p->position.z);
      p->hydroParticle().updateParticleVelocity(p->velocity.x,p->velocity.y,p->velocity.z);
      p->hydroParticle().updateParticleAcceleration( p->treeAcceleration.x,
                                                     p->treeAcceleration.y,
                                                     p->treeAcceleration.z);
    }
  }
  // Use shadow array to avoid reduction conflict
  smoothProxy[thisIndex].ckLocal()->contribute(cb);
}
#ifdef MHD
void TreePiece::updateHydroParticleBfield( const double bx, const double by, const double bz,
                                           const CkCallback &cb) {
  GravityParticle *p;
  for(int i=1; i<= myNumParticles; ++i) {
    p = &myParticles[i];
    if (TYPETest(p, TYPE_GAS)) {
      p->hydroParticle().updateGlobalBfield( bx, by, bz);
    }
  }
  // Use shadow array to avoid reduction conflict
  smoothProxy[thisIndex].ckLocal()->contribute(cb);
}
#endif

void TreePiece::checkHydroParticle( const CkCallback &cb) {
  GravityParticle *p;
  int numBadParticles = 0;
  for(int i=1; i<= myNumParticles; ++i) {
    p = &myParticles[i];
    if (TYPETest(p, TYPE_GAS))
      p->hydroParticle().check() ? 0 : numBadParticles++;
  }

  //ckout << "Num Bad Particles = " << numBadParticles << endl;
  // Use shadow array to avoid reduction conflict
  smoothProxy[thisIndex].ckLocal()->contribute(cb);
}

void TreePiece::setFBallHydroParticle( int activeRung, const CkCallback &cb) {
  GravityParticle *p;
  for(int i=1; i<= myNumParticles; ++i) {
    p = &myParticles[i];
    if (TYPETest(p, TYPE_GAS)) {
      if(p->hydroParticle().isValid()) p->fBall = p->hydroParticle().maxSearchRadius*1.05;
      else p->fBall = p->hydroParticle().maxSearchRadius*1.5;
    }
  }
  // Use shadow array to avoid reduction conflict
  smoothProxy[thisIndex].ckLocal()->contribute(cb);
}

void TreePiece::getGlobalInfo( const CkCallback &cb) {
  double dGlobalValues[8] = {0};
  varVector Q, deltaQ;
  // Loop over all particles to sum mass and mass*position for only star
  // particles
  for (unsigned int i = 0; i < myNumParticles; ++i) {
    GravityParticle *p = &myParticles[i+1];
    if (TYPETest(p, TYPE_GAS)) {
      p->hydroParticle().getQ( Q);
      p->hydroParticle().getDeltaQ( deltaQ);
      dGlobalValues[0] += Q[iRho];
      dGlobalValues[1] += deltaQ[iRho];
      dGlobalValues[2] += Q[ipx];
      dGlobalValues[3] += Q[ipy];
      dGlobalValues[4] += Q[ipz];
      dGlobalValues[5] += 0.5*(Q[ipx]*Q[ipx] + Q[ipy]*Q[ipy] + Q[ipz]*Q[ipz])/Q[iRho];
#ifdef RADIATION
      dGlobalValues[6] += p->hydroParticle().getRadEnergy()*p->hydroParticle().getVolume();
#endif
    }
  }
  contribute(8*sizeof(double), dGlobalValues, CkReduction::sum_double, cb);

}

void TreePiece::getGlobalStatistics( int activeRung, const CkCallback &cb) {
  int iStats[7] = {0};
  // Loop over all particles to sum mass and mass*position for only star
  // particles
  for (unsigned int i = 0; i < myNumParticles; ++i) {
    GravityParticle *p = &myParticles[i+1];
    if (TYPETest(p, TYPE_GAS)) {
      if(TYPETest( p, TYPE_NbrOfACTIVE) || p->hydroParticle().isActiveHalf( activeRung)) 
        p->hydroParticle().isValid() ? iStats[0]++ : iStats[1]++;
      if (TYPETest( p, TYPE_NbrOfACTIVE) || p->hydroParticle().isActive( activeRung))
        p->hydroParticle().isValid() ? iStats[2]++ : iStats[3]++;
      if( TYPETest( p, TYPE_NbrOfACTIVE)) iStats[4]++;
      if( p->hydroParticle().isActive( activeRung)) iStats[5]++;
      if( p->hydroParticle().isActiveHalf( activeRung)) iStats[6]++;
    }
    
  }
  contribute(7*sizeof(int), iStats, CkReduction::sum_int, cb);

}

void TreePiece::getGlobalDistributions( const CkCallback &cb) {
  int iStats[1] = {0};
  // Loop over all particles to sum mass and mass*position for only star
  // particles
  int nmin = 1000, nmax = 0, nsearch = 0;
  for (unsigned int i = 0; i < myNumParticles; ++i) {
    GravityParticle *p = &myParticles[i+1];
    if (TYPETest(p, TYPE_GAS)) {
      nsearch = p->hydroParticle().getNumSearch();
      nmax = max(nsearch, nmax);
      nmin = min(nsearch, nmin);
    }

  }
  iStats[0] = nmax;
  contribute(1*sizeof(int), iStats, CkReduction::max_int, cb);
}

#ifdef MHD
void TreePiece::getGlobalBfield( double dt, const CkCallback &cb) {
  double dGlobalBfield[11] = {0};
  varVector Q, deltaQ, conserved;
  dVector cGradients[NVARS];
  // Loop over all particles to sum mass and mass*position for only star
  // particles
  for (unsigned int i = 0; i < myNumParticles; ++i) {
    GravityParticle *p = &myParticles[i+1];
    if (TYPETest(p, TYPE_GAS)) {
      p->hydroParticle().getQ( Q);
      p->hydroParticle().getGradients( cGradients);
      p->hydroParticle().getConserved( conserved);
      double volume = p->hydroParticle().getVolume();
      dVector v(conserved[ipx], conserved[ipy], conserved[ipz]);
      dVector Bfield(conserved[iBx], conserved[iBy], conserved[iBz]);
      dVector rhov = v;
      v /= conserved[iRho];
      double totalEnergy = conserved[ipx]*v[0] + conserved[ipy]*v[1] + conserved[ipz]*v[2];
      totalEnergy *= 0.5*volume;
      double totalMass =conserved[iRho] * volume;
      double v2A = dotProduct(Bfield, Bfield) * volume;
      double cs2 = p->hydroParticle().getSoundSpeed(); cs2 = cs2*cs2*volume*conserved[iRho];
      dVector vdotGradB(0.,0.,0.);
      vdotGradB[0] = dotProduct(v,cGradients[iBx]);
      vdotGradB[1] = dotProduct(v,cGradients[iBy]);
      vdotGradB[2] = dotProduct(v,cGradients[iBz]);

      dVector BdotGradv(0.,0.,0.);
      dVector Gradv[3];
      for( int i = 0; i < 3; i++)
        Gradv[i] = (cGradients[ipx+i] - v[i]*cGradients[iRho])/conserved[iRho];

      BdotGradv[0] = dotProduct(Bfield,Gradv[0]);
      BdotGradv[1] = dotProduct(Bfield,Gradv[1]);
      BdotGradv[2] = dotProduct(Bfield,Gradv[2]);
      double divRhoV =  cGradients[ipx][0] + cGradients[ipy][1] + cGradients[ipz][2];
      double div_vB = (divRhoV - dotProduct(v,cGradients[iRho]))/conserved[iRho];
      
      dGlobalBfield[3] += -dt*volume*(vdotGradB[0] - BdotGradv[0] + Bfield[0]*div_vB);
      dGlobalBfield[4] += -dt*volume*(vdotGradB[1] - BdotGradv[1] + Bfield[1]*div_vB);
      dGlobalBfield[5] += -dt*volume*(vdotGradB[2] - BdotGradv[2] + Bfield[2]*div_vB);
      dGlobalBfield[0] += Q[iBx];
      dGlobalBfield[1] += Q[iBy];
      dGlobalBfield[2] += Q[iBz];
      dGlobalBfield[6] += volume;
      dGlobalBfield[7] += totalEnergy;
      dGlobalBfield[8] += totalMass;
      dGlobalBfield[9] += v2A;
      dGlobalBfield[10] += cs2;

    }
  }
  contribute(11*sizeof(double), dGlobalBfield, CkReduction::sum_double, cb);

}
#endif

void TreePiece::infoFromHydroParticle( int activeRung, double voroKick, const CkCallback &cb, const bool activeOnly) {
  GravityParticle *p;
  voroKick = min(max(voroKick,0.), 1.);
#ifndef COOLING_NONE
  HydroUtils::setupEOS( dm->Cool);
#endif
  for(int i=1; i<= myNumParticles; ++i) {
    p = &myParticles[i];
    if (TYPETest(p, TYPE_GAS) && (!activeOnly || p->hydroParticle().isActive( activeRung))) {
#ifdef DAMPING
      p->hydroParticle().setRelaxVelocity();
#endif
      p->hydroParticle().update();
      dVector vel2 = p->hydroParticle().getVelocity();
      std::isnan( vel2[0]) && printf( "vel2[0] nan: %5.3e %5.3e\n", vel2[0], voroKick);

      p->hydroParticle().kick( voroKick);
      dVector vel = p->hydroParticle().getVelocity();
      std::isnan( vel[0]) && printf( "vel[0] nan: %5.3e %5.3e\n", vel[0], voroKick);

      p->mass = p->hydroParticle().getMass();
      varVector primitives;
      p->hydroParticle().getPrimitives( primitives);
      p->u() = primitives[iE];
      p->fDensity = primitives[iRho];
#ifdef COOLING_SROEOS
      p->fMetals() = primitives[iye];
#endif
#ifdef COOLING_MESA
      p->fMetals() = primitives[iX];
#endif

      for( int j=0; j<NDIM; j++) {
        p->velocity[j] = vel[j];
      }

      p->hydroParticle().updateParticleVelocity(p->velocity.x,p->velocity.y,p->velocity.z);

    }
  }
  // Use shadow array to avoid reduction conflict
  smoothProxy[thisIndex].ckLocal()->contribute(cb);
}

// Setup the HydroParticle Rungs 
void TreePiece::setupHydroParticleRungs( int _iRung, const CkCallback &cb) {
  GravityParticle *p;
  for(int i=1; i<= myNumParticles; ++i) {
    p = &myParticles[i];
    if (TYPETest(p, TYPE_GAS)) {
      if( _iRung < 0) { 
        p->hydroParticle().setRung( p->rung);
      }
      else { 
        p->hydroParticle().setRung( _iRung);
        p->rung = _iRung;
      }
      //p->rung = p->hydroParticle().getRung();
      //p->rung = p->hydroParticle().getHalfRung();
    }
  }
  // Use shadow array to avoid reduction conflict
  smoothProxy[thisIndex].ckLocal()->contribute(cb);
}

void DataManager::setupMovingMesh( Parameters param, const CkCallback &cb) {
  HydroUtils::setupBC( param.vPeriod.x, param.vPeriod.y, param.vPeriod.z);
  HydroUtils::instance->particleTrack = param.iParticleTrack;
  HydroUtils::instance->firstOrder = param.bFirstOrder;
  HydroUtils::instance->smallRho = param.dSmallRho;
  HydroUtils::instance->useEntropy = param.bUseEntropy;
  HydroUtils::instance->useOldGravMethod = param.bUseOldGravMethod;
  HydroUtils::instance->useVanLeer = param.bUseVanLeer != 0;
  HydroUtils::instance->machThreshold = param.dMachThreshold;
  HydroUtils::instance->theta = param.dHydroTheta;
  HydroUtils::instance->invRelaxTau = param.dGlassDamper;
  HydroUtils::instance->vPsi = param.dVPsi*param.dSecUnit/(3.0857e21*param.dKpcUnit);
  HydroUtils::instance->voroMesh = fmin(1.,fmax(0., param.dVoroMesh));
  HydroUtils::instance->minTimeStep = fmin(fmax( 0., param.dDeltaMin), param.dDelta);
  HydroUtils::instance->useEOSSpeedup = param.bUseEOSSpeedup;
  HydroUtils::instance->checkNaN = param.bDoCheckNaN;
  HydroUtils::defaultGamma = param.dConstGamma;
  HydroUtils::defaultMu = param.dMeanMolWeight;
  HydroUtils::riemannSolver = param.iRiemannSolver;
  HydroUtils::instance->cInCodeUnits = SPEED_OF_LIGHT*param.dSecUnit/(3.0857e21*param.dKpcUnit);
#ifdef RADIATION
  if( param.dRedSpeedOfLight > 0.) 
    HydroUtils::instance->redCInCodeUnits = param.dRedSpeedOfLight*param.dSecUnit/(3.0857e21*param.dKpcUnit);
  else 
    HydroUtils::instance->redCInCodeUnits = HydroUtils::instance->cInCodeUnits;
  printf("iAngleTrack: %d\n", param.iAngleTrack);
  HydroUtils::instance->angleTrack = param.iAngleTrack;
  HydroUtils::instance->initRadNormals();
  UserSources::registerSource( new RadiationSource());
  UserSources::registerSource( new UserRadiationSources());
  std::string opacityName( param.opacityName);
  Opacity::registerOpacity( opacityName, new SimpleOpacity());
  Opacity::registerOpacity( opacityName, new DustOpacity());
  Opacity::initOpacityParameters( param);
  const double Ti = param.dRadTinit;
  const double Ti4 = Ti*Ti*Ti*Ti;
  const double Iinit = 7.5657e-15*Ti4/(param.dErgPerGmUnit*param.dGmPerCcUnit)/(4.*3.1415);
  HydroUtils::instance->initRadIntensity = Iinit;
  HydroUtils::instance->initRadMu = max(param.dRadMuinit, -1.);
#endif  
  UserSources::registerSource( new FixedBC());
  UserSources::initSourceParameters( param);
  std::string NoneProblemGeneratorName = "None";
  if( NoneProblemGeneratorName.compare( param.problemGeneratorName) != 0) { 
    std::string problemName( param.problemGeneratorName);
    UserProblemGenerator::registerUserProblemGenerator( problemName, new WaveProblemGenerator());
    UserProblemGenerator::registerUserProblemGenerator( problemName, new DynDiffusionProblemGenerator());
    UserProblemGenerator::initUserProblemGeneratorParameters( param);
  }
  contribute(cb);
}
#endif //VORONOI
