/* program to create an initial Poisson distribution particles */
/* This is stripped down from the version in "initial" to make it
   easier to use for a non-cosmologist. */
/* This version produces gas to test SPH */
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
      fprintf(stderr, "Usage: pgaspartt ngas xmin xmax seed\n");
      return 1;
      }
  
  xmin = atof(argv[2]);
  xmax = atof(argv[3]);
  seed = atoi(argv[4]);

  srand(seed);
  
  h.nsph = atoi(argv[1]);
  h.ndark = 0;
  h.nbodies = h.nsph;
  h.nstar = 0;
  h.ndim = 3;
  h.time = 0.1;

  mass = 1.0/h.nsph; /* scale so the total mass is 1  (Omega = 1.0) */
  eps = pow(h.nsph, -1.0/3.0)*(xmax - xmin)/20.0;
  
  xdrstdio_create(&xdrs, stdout, XDR_ENCODE);

  xdr_header(&xdrs, &h);

  deltax = (xmax - xmin);
  for(i = 0; i < h.nsph; i++)
    {
	  struct gas_particle gp;

	  gp.mass = mass;
	  gp.pos[0] = xmin + rand()/((double) RAND_MAX)*deltax;
	  gp.pos[1] = xmin + rand()/((double) RAND_MAX)*deltax;
	  gp.pos[2] = xmin + rand()/((double) RAND_MAX)*deltax;
	  gp.vel[0] = 0.0;
	  gp.vel[1] = 0.0;
	  gp.vel[2] = 0.0;
	  gp.hsmooth = eps;
	  gp.rho = 0.0;
	  gp.temp = 100.0;
	  gp.phi = 0.0;

	  xdr_gas(&xdrs, &gp);
	}
  xdr_destroy(&xdrs);
  return(0);
}
