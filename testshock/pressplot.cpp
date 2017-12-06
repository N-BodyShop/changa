/* 
 * Analysis program to generate a run of pressure vs. z for the test shock
 * tube.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "tipsydefs.h"
#include <rpc/types.h>
#include <rpc/xdr.h>

int
main(argc, argv)
     int argc;
     char **argv;
{
  double zmin;
  double zmax;
  struct dump h;
  double deltaz;
  int i;
  int nz;
  double *pressa;
  int *na;
  XDR xdrs;

  if(argc != 4) {
      fprintf(stderr, "Usage: pressplot zmin zmax nz\n");
      return 1;
      }
  
  zmin = atof(argv[1]);
  zmax = atof(argv[2]);
  nz = atoi(argv[3]);
  pressa = malloc(nz*sizeof(double)) ;
  na = malloc(nz*sizeof(int));
  for(i = 0; i < nz; i++) {
	pressa[i] = 0.0;
	na[i] = 0;
	}

  xdrstdio_create(&xdrs, stdin, XDR_DECODE);

  xdr_header(&xdrs, &h);

  deltaz = (zmax - zmin)/nz;
  for(i = 0; i < h.nsph; i++)
    {
	  int iz;
	  struct gas_particle gp;
	  xdr_gas(&xdrs, &gp);

	  iz = (gp.pos[2] - zmin)/deltaz;
	  if(iz < 0) continue;
	  if(iz >= nz) continue;
	  na[iz]++;
          pressa[iz] +=  gp.rho*gp.temp;
	}
  xdr_destroy(&xdrs);
  for(i = 0; i < nz; i++) {
	printf("%d %g %d %g\n", i, zmin + deltaz*i, na[i], pressa[i]/na[i]);
	}
  return(0);
}
