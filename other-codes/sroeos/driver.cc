#include<stdio.h>
#include<stdlib.h>
#include "nuc_eos.hh"
#include<cmath>

int main(void) {

  using namespace nuc_eos;

  //nuc_eos_C_ReadTable("LS220_234r_136t_50y_analmu_20091212_SVNr26.h5");
  nuc_eos_readtable_("LS220_234r_136t_50y_analmu_20091212_SVNr26.h5");
  printf("nrho: %d\n",nrho);
  printf("ntemp: %d\n",ntemp);
  printf("nye: %d\n",nye);

  double xrho = 2.74405e14 * RHOGF;
  double xtemp = 8.12426;
  double xye =  0.267428;
  double xeps = 0.0;
  double xprs = 0.0;
  int keyerr = 0;
  int anyerr = 0;
  int n = 1;
  
  nuc_eos_m_kt1_press_eps(&n,&xrho,&xtemp,&xye,&xeps,&xprs,
			  &keyerr,&anyerr);

  fprintf(stderr,"energy_shift: %15.6E\n",energy_shift);
  fprintf(stderr,"%15.6E %15.6E %15.6E\n",xrho,xtemp,xye);
  fprintf(stderr,"%15.6E %15.6E %d %d\n",xprs,xeps,keyerr,anyerr);

  // now see what we get
  const double prec = 1.0e-12;
  //  xeps = 4.0e-4;
  nuc_eos_m_kt0_press(&n,&xrho,&xtemp,&xye,&xeps,&xprs,
		      &prec,&keyerr,&anyerr);

  fprintf(stderr,"*********************\n");
  fprintf(stderr,"%15.6E %15.6E %15.6E\n",xrho,xtemp,xye);
  fprintf(stderr,"%15.6E %15.6E %d %d\n",xprs,xeps,keyerr,anyerr);

  fprintf(stderr,"******************************************\n");
  fprintf(stderr,"Short EOS Call:\n");
  xrho = 2.0e14 * RHOGF;
  xtemp = 5.0;
  xye = 0.25;

  // we assume we know our temperature
  // declare additional vars
  double xent,xcs2,xdedt,xdpderho,xdpdrhoe,xmunu;
  nuc_eos_m_kt1_short(&n,&xrho,&xtemp,&xye,
		      &xeps,&xprs,&xent,&xcs2,&xdedt,
		      &xdpderho,&xdpdrhoe,&xmunu,
		      &keyerr,&anyerr);
  fprintf(stderr,"     rho: %15.6E     temp: %15.6E    ye: %15.6E\n",xrho*INVRHOGF,xtemp,xye);
  fprintf(stderr,"    xprs: %15.6E     xeps: %15.6E\n",xprs*INVPRESSGF,xeps);
  fprintf(stderr,"  keyerr: %15d   anyerr: %15d\n",keyerr,anyerr);
  fprintf(stderr,"    xent: %15.6E      xcs: %15.6E xdedt: %15.6E\n",xent,sqrt(xcs2),xdedt);
  fprintf(stderr,"xdpderho: %15.6E xdpdrhoe: %15.6E \n",xdpderho,xdpdrhoe);
  fprintf(stderr,"   xmunu: %15.6E\n",xmunu);
  

# if 0
  double xxa,xxh,xxn,xxp;
  double xabar,xzbar,xmue,xmun,xmup,xmuhat;
  nuc_eos_m_kt1_full(&n,&xrho,&xtemp,&xye,&xeps,
		     &xprs,&xent,&xcs2,&xdedt,&xdpderho,&xdpdrhoe,
		     &xxa,&xxh,&xxn,&xxp,&xabar,&xzbar,
		     &xmue,&xmun,&xmup,&xmuhat,&keyerr,&anyerr);

  fprintf(stderr,"*********************\n");
  fprintf(stderr,"%15.6E %15.6E %15.6E\n",xrho,xtemp,xye);
  fprintf(stderr,"%15.6E %15.6E %d %d\n",xprs,xeps,keyerr,anyerr);
  fprintf(stderr,"xent: %15.6E xcs2: %15.6E xdedt: %15.6E\n",xent,xcs2,xdedt);
  fprintf(stderr,"xdpderho: %15.6E xdpdrhoe: %15.6E \n",xdpderho,xdpdrhoe);
  fprintf(stderr,"xxa: %15.6E xxh: %15.6E xxn: %15.6E xxp: %15.6e\n",
	  xxa,xxh,xxn,xxp);
  fprintf(stderr,"xabar: %15.6E xzbar: %15.6E \n",xabar,xzbar);
  fprintf(stderr,"xmue: %15.6E xmun: %15.6E  xmup: %15.6E  xmuhat: %15.6E\n",
	  xmue,xmun,xmup,xmuhat);


#endif

  //  testeos();

  return 0;
}
