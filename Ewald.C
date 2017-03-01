#include <assert.h>
#include "ParallelGravity.h"
// Ewald summation code.
// First implemented by Thomas Quinn and Joachim Stadel in PKDGRAV.

#ifdef HEXADECAPOLE
inline
static
void QEVAL(MOMC mom, double gam[], double dx, double dy,
		  double dz, double &ax, double &ay, double &az, double &fPot)
{
    double Qmirx,Qmiry,Qmirz,Qmir,Qta;
    Qta = 0.0;
    Qmirx = (1.0/6.0)*(mom.xzzz*dz*dz*dz + 3*mom.xyzz*dy*dz*dz + 3*mom.xyyz*dy*dy*dz + mom.xyyy*dy*dy*dy + 3*mom.xxzz*dx*dz*dz + 6*mom.xxyz*dx*dy*dz + 3*mom.xxyy*dx*dy*dy + 3*mom.xxxz*dx*dx*dz + 3*mom.xxxy*dx*dx*dy + mom.xxxx*dx*dx*dx);
    Qmiry = (1.0/6.0)*(mom.yzzz*dz*dz*dz + 3*mom.xyzz*dx*dz*dz + 3*mom.xxyz*dx*dx*dz + mom.xxxy*dx*dx*dx + 3*mom.yyzz*dy*dz*dz + 6*mom.xyyz*dx*dy*dz + 3*mom.xxyy*dx*dx*dy + 3*mom.yyyz*dy*dy*dz + 3*mom.xyyy*dx*dy*dy + mom.yyyy*dy*dy*dy);
    Qmirz = (1.0/6.0)*(mom.yyyz*dy*dy*dy + 3*mom.xyyz*dx*dy*dy + 3*mom.xxyz*dx*dx*dy + mom.xxxz*dx*dx*dx + 3*mom.yyzz*dy*dy*dz + 6*mom.xyzz*dx*dy*dz + 3*mom.xxzz*dx*dx*dz + 3*mom.yzzz*dy*dz*dz + 3*mom.xzzz*dx*dz*dz + mom.zzzz*dz*dz*dz);
    Qmir = (1.0/4.0)*(Qmirx*dx + Qmiry*dy + Qmirz*dz);
    fPot -= gam[4]*Qmir;
    Qta += gam[5]*Qmir;
    ax += gam[4]*Qmirx;
    ay += gam[4]*Qmiry;
    az += gam[4]*Qmirz;
    Qmirx = (1.0/2.0)*(mom.xzz*dz*dz + 2*mom.xyz*dy*dz + mom.xyy*dy*dy + 2*mom.xxz*dx*dz + 2*mom.xxy*dx*dy + mom.xxx*dx*dx);
    Qmiry = (1.0/2.0)*(mom.yzz*dz*dz + 2*mom.xyz*dx*dz + mom.xxy*dx*dx + 2*mom.yyz*dy*dz + 2*mom.xyy*dx*dy + mom.yyy*dy*dy);
    Qmirz = (1.0/2.0)*(mom.yyz*dy*dy + 2*mom.xyz*dx*dy + mom.xxz*dx*dx + 2*mom.yzz*dy*dz + 2*mom.xzz*dx*dz + mom.zzz*dz*dz);
    Qmir = (1.0/3.0)*(Qmirx*dx + Qmiry*dy + Qmirz*dz);
    fPot -= gam[3]*Qmir;
    Qta += gam[4]*Qmir;
    ax += gam[3]*Qmirx;
    ay += gam[3]*Qmiry;
    az += gam[3]*Qmirz;
    Qmirx = (1.0/1.0)*(mom.xz*dz + mom.xy*dy + mom.xx*dx);
    Qmiry = (1.0/1.0)*(mom.yz*dz + mom.xy*dx + mom.yy*dy);
    Qmirz = (1.0/1.0)*(mom.yz*dy + mom.xz*dx + mom.zz*dz);
    Qmir = (1.0/2.0)*(Qmirx*dx + Qmiry*dy + Qmirz*dz);
    fPot -= gam[2]*Qmir;
    Qta += gam[3]*Qmir;
    ax += gam[2]*Qmirx;
    ay += gam[2]*Qmiry;
    az += gam[2]*Qmirz;
    fPot -= gam[0]*mom.m;
    Qta += gam[1]*mom.m;
    ax -= dx*Qta;
    ay -= dy*Qta;
    az -= dz*Qta;
}
#else
inline
static
void QEVAL(MultipoleMoments mom, double gam[], double dx, double dy,
		  double dz, double &ax, double &ay, double &az, double &fPot)
{
    double Qmirx,Qmiry,Qmirz,Qmir,Qta;
    Qta = 0.0;
    Qmirx = (1.0/1.0)*(mom.xz*dz + mom.xy*dy + mom.xx*dx);
    Qmiry = (1.0/1.0)*(mom.yz*dz + mom.xy*dx + mom.yy*dy);
    Qmirz = (1.0/1.0)*(mom.yz*dy + mom.xz*dx + mom.zz*dz);
    Qmir = (1.0/2.0)*(Qmirx*dx + Qmiry*dy + Qmirz*dz);
    fPot -= gam[2]*Qmir;
    Qta += gam[3]*Qmir;
    ax += gam[2]*Qmirx;
    ay += gam[2]*Qmiry;
    az += gam[2]*Qmirz;
    fPot -= gam[0]*mom.totalMass;
    Qta += gam[1]*mom.totalMass;
    ax -= dx*Qta;
    ay -= dy*Qta;
    az -= dz*Qta;
}
#endif

void TreePiece::BucketEwald(GenericTreeNode *req, int nReps,double fEwCut)
{
#ifndef BENCHMARK_NO_WORK
	GravityParticle *p;
#ifdef HEXADECAPOLE
	MOMC mom = momcRoot;
	MultipoleMoments momQuad = root->moments;
	double xx,xxx,xxy,xxz,yy,yyy,yyz,xyy,zz,zzz,xzz,yzz,xy,xyz,xz,yz;
	double Q4mirx,Q4miry,Q4mirz,Q4mir,Q4x,Q4y,Q4z;
	double Q4xx,Q4xy,Q4xz,Q4yy,Q4yz,Q4zz,Q4,Q3x,Q3y,Q3z;
	double Q3mirx,Q3miry,Q3mirz,Q3mir;
	const double onethird = 1.0/3.0;
#else
	MultipoleMoments mom = root->moments;
#endif
	double Q2;
	double L,fEwCut2,fInner2,alpha,alpha2,alphan,k1,ka;
	double fPot,ax,ay,az;
	double dx,dy,dz,x,y,z,r2,dir,dir2,a;
	double Q2mirx,Q2miry,Q2mirz,Q2mir,Qta;
	double g0,g1,g2,g3,g4,g5;
	double hdotx,s,c;
	int i,j,n,ix,iy,iz,nEwReps,bInHole,bInHolex,bInHolexy;
	int nLoop = 0;

	/*
	 ** Set up traces of the complete multipole moments.
	 */
#ifdef HEXADECAPOLE
        Q4xx = 0.5*(mom.xxxx + mom.xxyy + mom.xxzz);
        Q4xy = 0.5*(mom.xxxy + mom.xyyy + mom.xyzz);
        Q4xz = 0.5*(mom.xxxz + mom.xyyz + mom.xzzz);
        Q4yy = 0.5*(mom.xxyy + mom.yyyy + mom.yyzz);
        Q4yz = 0.5*(mom.xxyz + mom.yyyz + mom.yzzz);
        Q4zz = 0.5*(mom.xxzz + mom.yyzz + mom.zzzz);
        Q4 = 0.25*(Q4xx + Q4yy + Q4zz);
        Q3x = 0.5*(mom.xxx + mom.xyy + mom.xzz);
        Q3y = 0.5*(mom.xxy + mom.yyy + mom.yzz);
        Q3z = 0.5*(mom.xxz + mom.yyz + mom.zzz);
#endif
	Q2 = 0.5*(mom.xx + mom.yy + mom.zz);	

	n = req->lastParticle - req->firstParticle + 1;
	p = &myParticles[req->firstParticle];
	nEwReps = (int) ceil(fEwCut);
	L = fPeriod.x;
	fEwCut2 = fEwCut*fEwCut*L*L;
	fInner2 = 1.2e-3*L*L;
	nEwReps = nEwReps > nReps ? nEwReps : nReps;
	alpha = 2.0/L;
	alpha2 = alpha*alpha;
	k1 = M_PI/(alpha2*L*L*L);
	ka = 2.0*alpha/sqrt(M_PI);
	for(j=0;j<n;++j) {
                if (p[j].rung < activeRung) continue;
#ifdef HEXADECAPOLE
		fPot = momQuad.totalMass*k1;
#else
		fPot = mom.totalMass*k1;
#endif
		ax = 0.0;
		ay = 0.0;
		az = 0.0;
#ifdef HEXADECAPOLE
		dx = p[j].position.x - momQuad.cm.x;
		dy = p[j].position.y - momQuad.cm.y;
		dz = p[j].position.z - momQuad.cm.z;
#else
		dx = p[j].position.x - mom.cm.x;
		dy = p[j].position.y - mom.cm.y;
		dz = p[j].position.z - mom.cm.z;
#endif
		for (ix=-nEwReps;ix<=nEwReps;++ix) {
			bInHolex = (ix >= -nReps && ix <= nReps);
			x = dx + ix*L;
			for(iy=-nEwReps;iy<=nEwReps;++iy) {
				bInHolexy = (bInHolex && iy >= -nReps && iy <= nReps);
				y = dy + iy*L;
				for(iz=-nEwReps;iz<=nEwReps;++iz) {
					bInHole = (bInHolexy && iz >= -nReps && iz <= nReps);
					/*
					 ** Scoring for Ewald inner stuff = (+,*)
					 **		Visible ops 		= (104,161)
					 **		sqrt, 1/sqrt est. 	= (6,11)
					 **     division            = (6,11)  same as sqrt.
					 **		exp est.			= (6,11)  same as sqrt.
					 **		erf/erfc est.		= (12,22) twice a sqrt.	
					 **		Total			= (128,205) = 333
					 **     Old scoring				    = 447
					 */
					z = dz + iz*L;
					r2 = x*x + y*y + z*z;
					if (r2 > fEwCut2 && !bInHole) continue;
					if (r2 < fInner2) {
						/*
						 * For small r, series expand about
						 * the origin to avoid errors caused
						 * by cancellation of large terms.
						 */
						alphan = ka;
						r2 *= alpha2;
						g0 = alphan*((1.0/3.0)*r2 - 1.0);
						alphan *= 2*alpha2;
						g1 = alphan*((1.0/5.0)*r2 - (1.0/3.0));
						alphan *= 2*alpha2;
						g2 = alphan*((1.0/7.0)*r2 - (1.0/5.0));
						alphan *= 2*alpha2;
						g3 = alphan*((1.0/9.0)*r2 - (1.0/7.0));
						alphan *= 2*alpha2;
						g4 = alphan*((1.0/11.0)*r2 - (1.0/9.0));
						alphan *= 2*alpha2;
						g5 = alphan*((1.0/13.0)*r2 - (1.0/11.0));
						}
					else {
					    dir = 1/sqrt(r2);
					    dir2 = dir*dir;
					    a = exp(-r2*alpha2);
					    a *= ka*dir2;
					    if (bInHole) g0 = -erf(alpha/dir);
					    else g0 = erfc(alpha/dir);
					    g0 *= dir;
					    g1 = g0*dir2 + a;
					    alphan = 2*alpha2;
					    g2 = 3*g1*dir2 + alphan*a;
					    alphan *= 2*alpha2;
					    g3 = 5*g2*dir2 + alphan*a;
					    alphan *= 2*alpha2;
						g4 = 7*g3*dir2 + alphan*a;
					    alphan *= 2*alpha2;
					    g5 = 9*g4*dir2 + alphan*a;
					    }
#ifdef HEXADECAPOLE
					xx = 0.5*x*x;
					xxx = onethird*xx*x;
					xxy = xx*y;
					xxz = xx*z;
					yy = 0.5*y*y;
					yyy = onethird*yy*y;
					xyy = yy*x;
					yyz = yy*z;
					zz = 0.5*z*z;
					zzz = onethird*zz*z;
					xzz = zz*x;
					yzz = zz*y;
					xy = x*y;
					xyz = xy*z;
					xz = x*z;
					yz = y*z;
					Q2mirx = mom.xx*x + mom.xy*y + mom.xz*z;
					Q2miry = mom.xy*x + mom.yy*y + mom.yz*z;
					Q2mirz = mom.xz*x + mom.yz*y + mom.zz*z;
					Q3mirx = mom.xxx*xx + mom.xxy*xy + mom.xxz*xz + mom.xyy*yy + mom.xyz*yz + mom.xzz*zz;
					Q3miry = mom.xxy*xx + mom.xyy*xy + mom.xyz*xz + mom.yyy*yy + mom.yyz*yz + mom.yzz*zz;
					Q3mirz = mom.xxz*xx + mom.xyz*xy + mom.xzz*xz + mom.yyz*yy + mom.yzz*yz + mom.zzz*zz;
					Q4mirx = mom.xxxx*xxx + mom.xxxy*xxy + mom.xxxz*xxz + mom.xxyy*xyy + mom.xxyz*xyz +
						mom.xxzz*xzz + mom.xyyy*yyy + mom.xyyz*yyz + mom.xyzz*yzz + mom.xzzz*zzz;
					Q4miry = mom.xxxy*xxx + mom.xxyy*xxy + mom.xxyz*xxz + mom.xyyy*xyy + mom.xyyz*xyz +
						mom.xyzz*xzz + mom.yyyy*yyy + mom.yyyz*yyz + mom.yyzz*yzz + mom.yzzz*zzz;
					Q4mirz = mom.xxxz*xxx + mom.xxyz*xxy + mom.xxzz*xxz + mom.xyyz*xyy + mom.xyzz*xyz +
						mom.xzzz*xzz + mom.yyyz*yyy + mom.yyzz*yyz + mom.yzzz*yzz + mom.zzzz*zzz;
					Q4x = Q4xx*x + Q4xy*y + Q4xz*z;
					Q4y = Q4xy*x + Q4yy*y + Q4yz*z;
					Q4z = Q4xz*x + Q4yz*y + Q4zz*z;
					Q2mir = 0.5*(Q2mirx*x + Q2miry*y + Q2mirz*z) - (Q3x*x + Q3y*y + Q3z*z) + Q4;
					Q3mir = onethird*(Q3mirx*x + Q3miry*y + Q3mirz*z) - 0.5*(Q4x*x + Q4y*y + Q4z*z);
					Q4mir = 0.25*(Q4mirx*x + Q4miry*y + Q4mirz*z);
					Qta = g1*mom.m - g2*Q2 + g3*Q2mir + g4*Q3mir + g5*Q4mir;
					fPot -= g0*mom.m - g1*Q2 + g2*Q2mir + g3*Q3mir + g4*Q4mir;
					ax += g2*(Q2mirx - Q3x) + g3*(Q3mirx - Q4x) + g4*Q4mirx - x*Qta;
					ay += g2*(Q2miry - Q3y) + g3*(Q3miry - Q4y) + g4*Q4miry - y*Qta;
					az += g2*(Q2mirz - Q3z) + g3*(Q3mirz - Q4z) + g4*Q4mirz - z*Qta;
#else
					Q2mirx = mom.xx*x + mom.xy*y + mom.xz*z;
					Q2miry = mom.xy*x + mom.yy*y + mom.yz*z;
					Q2mirz = mom.xz*x + mom.yz*y + mom.zz*z;
					Q2mir = 0.5*(Q2mirx*x + Q2miry*y + Q2mirz*z);
					Qta = g1*mom.totalMass - g2*Q2 + g3*Q2mir;
					fPot -= g0*mom.totalMass - g1*Q2 + g2*Q2mir;
					ax += g2*(Q2mirx) - x*Qta;
					ay += g2*(Q2miry) - y*Qta;
					az += g2*(Q2mirz) - z*Qta;
#endif
					++nLoop;
					}
				}
			}
		/*
		 ** Scoring for the h-loop (+,*)
		 ** 	Without trig = (10,14)
		 **	    Trig est.	 = 2*(6,11)  same as 1/sqrt scoring.
		 **		Total        = (22,36)
		 **					 = 58
		 */
		for (i=0;i<nEwhLoop;++i) {
			hdotx = ewt[i].hx*dx + ewt[i].hy*dy + ewt[i].hz*dz;
			c = cos(hdotx);
			s = sin(hdotx);
			fPot += ewt[i].hCfac*c + ewt[i].hSfac*s;
			ax += ewt[i].hx*(ewt[i].hCfac*s - ewt[i].hSfac*c);
			ay += ewt[i].hy*(ewt[i].hCfac*s - ewt[i].hSfac*c);
			az += ewt[i].hz*(ewt[i].hCfac*s - ewt[i].hSfac*c);
			}
		p[j].potential += fPot;
		p[j].treeAcceleration.x += ax;
		p[j].treeAcceleration.y += ay;
		p[j].treeAcceleration.z += az;
	    }
	return;
#endif
}

// Set up table for Ewald h (Fourier space) loop

void TreePiece::EwaldInit()
{
	int i,hReps,hx,hy,hz,h2;
	double alpha,k4,L;
	double gam[6],mfacc,mfacs;
	double ax,ay,az;

#ifdef HEXADECAPOLE
	/* convert to complete moments */
	momRescaleFmomr(&(root->moments.mom),1.0f,root->moments.getRadius());
	momFmomr2Momc(&(root->moments.mom), &momcRoot);
	/* XXX note that we could leave the scaling as is and change
	   the radius of the root. */
	momRescaleFmomr(&(root->moments.mom),root->moments.getRadius(),1.0f);
#endif
	/*
	 ** Now setup stuff for the h-loop.
	 */
	hReps = (int) ceil(dEwhCut);
	L = fPeriod.x;
	alpha = 2.0/L;
	k4 = M_PI*M_PI/(alpha*alpha*L*L);
	i = 0;
	for (hx=-hReps;hx<=hReps;++hx) {
		for (hy=-hReps;hy<=hReps;++hy) {
			for (hz=-hReps;hz<=hReps;++hz) {
				h2 = hx*hx + hy*hy + hz*hz;
				if (h2 == 0) continue;
				if (h2 > dEwhCut*dEwhCut) continue;
				if (i == nMaxEwhLoop) {
				    nMaxEwhLoop *= 2;
				    /* avoid realloc() */
				    EWT *ewtTmp = new EWT[nMaxEwhLoop];
				    assert(ewtTmp != NULL);
				    for(int j = 0; j < i; j++) {
					ewtTmp[j] = ewt[j];
					}
				    delete[] ewt;
				    ewt = ewtTmp;
				    }
				gam[0] = exp(-k4*h2)/(M_PI*h2*L);
				gam[1] = 2*M_PI/L*gam[0];
				gam[2] = -2*M_PI/L*gam[1];
				gam[3] = 2*M_PI/L*gam[2];
				gam[4] = -2*M_PI/L*gam[3];
				gam[5] = 2*M_PI/L*gam[4];
				gam[1] = 0.0;
				gam[3] = 0.0;
				gam[5] = 0.0;
				ax = 0.0;
				ay = 0.0;
				az = 0.0;
				mfacc = 0.0;
#ifdef HEXADECAPOLE
				QEVAL(momcRoot, gam, hx, hy, hz,
				      ax, ay, az, mfacc);
#else
				QEVAL(root->moments, gam, hx, hy, hz,
				      ax, ay, az, mfacc);
#endif
				gam[0] = exp(-k4*h2)/(M_PI*h2*L);
				gam[1] = 2*M_PI/L*gam[0];
				gam[2] = -2*M_PI/L*gam[1];
				gam[3] = 2*M_PI/L*gam[2];
				gam[4] = -2*M_PI/L*gam[3];
				gam[5] = 2*M_PI/L*gam[4];
				gam[0] = 0.0;
				gam[2] = 0.0;
				gam[4] = 0.0;
				ax = 0.0;
				ay = 0.0;
				az = 0.0;
				mfacs = 0.0;
#ifdef HEXADECAPOLE
				QEVAL(momcRoot, gam,hx,hy,hz,
				      ax,ay,az,mfacs);
#else
				QEVAL(root->moments, gam,hx,hy,hz,
				      ax,ay,az,mfacs);
#endif
				ewt[i].hx = 2*M_PI/L*hx;
				ewt[i].hy = 2*M_PI/L*hy;
				ewt[i].hz = 2*M_PI/L*hz;
				ewt[i].hCfac = mfacc;
				ewt[i].hSfac = mfacs;
				++i;
				}
			}
		}
	nEwhLoop = i;

	//contribute(cb);
	dummyMsg *msg = new (8*sizeof(int)) dummyMsg;
	*((int *)CkPriorityPtr(msg)) = numTreePieces * numChunks + numTreePieces + thisIndex + 1;
	CkSetQueueing(msg,CK_QUEUEING_IFIFO);
	thisProxy[thisIndex].calculateEwald(msg);
}


void TreePiece::EwaldGPU() {
  /* when not using CUDA, definition is required because
     EwaldGPU is an entry method
  */
  
#ifdef SPCUDA
  int NumberOfGPUParticles = 0;
  if(FirstGPUParticleIndex != -1){
  	NumberOfGPUParticles = LastGPUParticleIndex - FirstGPUParticleIndex - 1; 
  } 
  else{
  	for (int i=0; i<numBuckets; i++){bucketReqs[i].finished = 1;}
  	return;
  }
  h_idata = (EwaldData*) malloc(sizeof(EwaldData));
  EwaldHostMemorySetup(h_idata, NumberOfGPUParticles, nEwhLoop);
  EwtData *ewtTable; 
  EwaldReadOnlyData *roData; 
  MultipoleMoments mm = root->moments;

  cudatype L = fPeriod.x;
  cudatype alpha = 2.0f/L;
  
  ewtTable = (EwtData*) h_idata->ewt;
  roData = (EwaldReadOnlyData*) h_idata->cachedData; 

  int nActive = NumberOfGPUParticles;
  h_idata->EwaldRange[0] = FirstGPUParticleIndex;
  h_idata->EwaldRange[1] = LastGPUParticleIndex; 

  for (int i=0; i<nEwhLoop; i++) {
    ewtTable[i].hx = (cudatype) ewt[i].hx; 
    ewtTable[i].hy = (cudatype) ewt[i].hy; 
    ewtTable[i].hz = (cudatype) ewt[i].hz; 
    ewtTable[i].hCfac = (cudatype) ewt[i].hCfac; 
    ewtTable[i].hSfac = (cudatype) ewt[i].hSfac; 
  }

#ifdef HEXADECAPOLE
  roData->momcRoot.m    = (cudatype)     momcRoot.m;
  roData->momcRoot.xx   = (cudatype)    momcRoot.xx;
  roData->momcRoot.yy   = (cudatype)    momcRoot.yy;
  roData->momcRoot.xy   = (cudatype)    momcRoot.xy;
  roData->momcRoot.xz   = (cudatype)    momcRoot.xz;
  roData->momcRoot.yz   = (cudatype)    momcRoot.yz;
  roData->momcRoot.xxx  = (cudatype)   momcRoot.xxx;
  roData->momcRoot.xyy  = (cudatype)   momcRoot.xyy;
  roData->momcRoot.xxy  = (cudatype)   momcRoot.xxy;
  roData->momcRoot.yyy  = (cudatype)   momcRoot.yyy;
  roData->momcRoot.xxz  = (cudatype)   momcRoot.xxz;
  roData->momcRoot.yyz  = (cudatype)   momcRoot.yyz;
  roData->momcRoot.xyz  = (cudatype)   momcRoot.xyz;
  roData->momcRoot.xxxx = (cudatype)  momcRoot.xxxx;
  roData->momcRoot.xyyy = (cudatype)  momcRoot.xyyy;
  roData->momcRoot.xxxy = (cudatype)  momcRoot.xxxy;
  roData->momcRoot.yyyy = (cudatype)  momcRoot.yyyy;
  roData->momcRoot.xxxz = (cudatype)  momcRoot.xxxz;
  roData->momcRoot.yyyz = (cudatype)  momcRoot.yyyz;
  roData->momcRoot.xxyy = (cudatype)  momcRoot.xxyy;
  roData->momcRoot.xxyz = (cudatype)  momcRoot.xxyz;
  roData->momcRoot.xyyz = (cudatype)  momcRoot.xyyz;
  roData->momcRoot.zz   = (cudatype)    momcRoot.zz;
  roData->momcRoot.xzz  = (cudatype)   momcRoot.xzz;
  roData->momcRoot.yzz  = (cudatype)   momcRoot.yzz;
  roData->momcRoot.zzz  = (cudatype)   momcRoot.zzz;
  roData->momcRoot.xxzz = (cudatype)  momcRoot.xxzz;
  roData->momcRoot.xyzz = (cudatype)  momcRoot.xyzz;
  roData->momcRoot.xzzz = (cudatype)  momcRoot.xzzz;
  roData->momcRoot.yyzz = (cudatype)  momcRoot.yyzz;
  roData->momcRoot.yzzz = (cudatype)  momcRoot.yzzz;
  roData->momcRoot.zzzz = (cudatype)  momcRoot.zzzz;
#else
  roData->mm.xx = (cudatype) mm.xx;
  roData->mm.xy = (cudatype) mm.xy;
  roData->mm.xz = (cudatype) mm.xz;
  roData->mm.yy = (cudatype) mm.yy;
  roData->mm.yz = (cudatype) mm.yz;
  roData->mm.zz = (cudatype) mm.zz;
#endif
  roData->mm.totalMass = (cudatype) mm.totalMass;
  roData->mm.cmx = (cudatype) mm.cm.x;
  roData->mm.cmy = (cudatype) mm.cm.y;
  roData->mm.cmz = (cudatype) mm.cm.z;
  roData->n = nActive;
  roData->fEwCut = (cudatype) fEwCut;
  roData->nReps = nReplicas ;
  roData->nEwReps = (int) ceil(fEwCut);
  roData->nEwhLoop = nEwhLoop;
  roData->L = (cudatype) L;
  roData->alpha = alpha;
  roData->alpha2 = (cudatype) alpha*alpha;
  roData->k1 = (cudatype) M_PI/(alpha*alpha*L*L*L);
  roData->ka = (cudatype) 2.0*alpha/sqrt(M_PI);
  roData->fEwCut2 = (cudatype) fEwCut*fEwCut*L*L;
/*
  Break between Taylor expansion for small r and multipole expansion.
  This value is for double precision.
  roData->fInner2 = (cudatype) 1.2e-3*L*L;
  The following is for single precision.  The CUDA version currently uses
  erff() and erfcf().  If these ever get changed then the following line
  needs to be changed accordingly.
 */
  roData->fInner2 = (cudatype) 1.1e-2*L*L;

  CkCallback *cb; 
  CkArrayIndex1D myIndex = CkArrayIndex1D(thisIndex); 
  cb = new CkCallback(CkIndex_TreePiece::EwaldGPUComplete(), myIndex, 
		      thisArrayID); 

  
  //CkPrintf("[%d] in EwaldGPU, calling EwaldHost\n", thisIndex);
#ifdef CUDA_INSTRUMENT_WRS
  EwaldHost(h_idata, (void *) cb, instrumentId, activeRung); 
#else
  int myLocalIndex;
  for(myLocalIndex = 0; this != dm->registeredTreePieces[myLocalIndex].treePiece;
      myLocalIndex++);
  CkAssert(myLocalIndex < dm->registeredTreePieces.length());
  
  EwaldHost(h_idata, (void *) cb, myLocalIndex); 
#endif

#endif
}


void TreePiece::EwaldGPUComplete() {
  /* when not using CUDA, definition is required because
     EwaldGPUComplete is an entry method
  */
#ifdef SPCUDA
  EwaldHostMemoryFree(h_idata); 
  free(h_idata); 

  /* indicate completion of ewald */
  
  for (int i=0; i<numBuckets; i++) {  
    bucketReqs[i].finished = 1; 
    finishBucket(i); 
  }
  //CkPrintf("[%d] in EwaldGPUComplete, completed book-keeping\n", thisIndex);
#endif 
}


