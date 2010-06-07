/* program to create an initial Poisson distribution particles */
/* This is stripped down from the version in "initial" to make it
   easier to use for a non-cosmologist. */
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
  double mass;
  double xmin;
  double xmax;
  double eps;
  struct dump h;
  double deltax;
  int i;
  int seed;
  XDR xdrs;

  if(argc != 5) {
      fprintf(stderr, "Usage: ppartt ndark xmin xmax seed\n");
      return 1;
      }
  
  xmin = atof(argv[2]);
  xmax = atof(argv[3]);
  seed = atoi(argv[4]);

  srand(seed);
  
  h.ndark = atoi(argv[1]);
  h.nbodies = h.ndark;
  h.nsph = 0;
  h.nstar = 0;
  h.ndim = 3;
  h.time = 0.1;

  mass = 1.0/h.ndark; /* scale so the total mass is 1  (Omega = 1.0) */
  eps = pow(h.ndark, -1.0/3.0)*(xmax - xmin)/20.0;
  
  xdrstdio_create(&xdrs, stdout, XDR_ENCODE);

  xdr_header(&xdrs, &h);

  deltax = (xmax - xmin);
  for(i = 0; i < h.ndark; i++)
    {
	  struct dark_particle dp;

	  dp.mass = mass;
	  dp.pos[0] = xmin + rand()/((double) RAND_MAX)*deltax;
	  dp.pos[1] = xmin + rand()/((double) RAND_MAX)*deltax;
	  dp.pos[2] = xmin + rand()/((double) RAND_MAX)*deltax;
	  dp.vel[0] = 0.0;
	  dp.vel[1] = 0.0;
	  dp.vel[2] = 0.0;
	  dp.eps = eps;
	  dp.phi = 0.0;

	  xdr_dark(&xdrs, &dp);
	}
  xdr_destroy(&xdrs);
  return(0);
}
