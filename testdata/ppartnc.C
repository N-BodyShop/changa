/* program to create an initial Poisson distribution particles */
/* This is stripped down from the version in "initial" to make it
   easier to use for a non-cosmologist. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include "xdr_template.h"
#include "tree_xdr.h"

int
main(int argc, char **argv)
{
  float mass;
  double xmin;
  double xmax;
  float eps;
  double deltax;
  int64_t i;
  int seed;
  XDR xdrs;
  FieldHeader fh;

  if(argc != 6) {
      fprintf(stderr, "Usage: ppartt ndark xmin xmax seed directory\n");
      return 1;
      }
  
  xmin = atof(argv[2]);
  xmax = atof(argv[3]);
  seed = atoi(argv[4]);

  srand(seed);
  
  fh.numParticles = atol(argv[1]);
  fh.dimensions = 1;
  fh.time = 0.1;
  fh.code = TypeHandling::float32;

  mass = 1.0/fh.numParticles; /* scale so the total mass is 1  (Omega = 1.0) */
  eps = pow(fh.numParticles, -1.0/3.0)*(xmax - xmin)/20.0;

  //make directory for files
  if(mkdir(argv[5], 0775) || chdir(argv[5])) {
      fprintf(stderr, "Could not create and move into a directory for the converted files, maybe you don't have permission?");
      return -1;
  }
  if(mkdir("dark", 0755))
      return -1;
  
  FILE* outfile = fopen("dark/mass", "wb");
  xdrstdio_create(&xdrs, outfile, XDR_ENCODE);

  xdr_template(&xdrs, &fh);
  for(i = 0; i < fh.numParticles+2; i++)
      xdr_template(&xdrs, &mass);
  xdr_destroy(&xdrs);
  fclose(outfile);

  outfile = fopen("dark/soft", "wb");
  xdrstdio_create(&xdrs, outfile, XDR_ENCODE);

  xdr_template(&xdrs, &fh);
  for(i = 0; i < fh.numParticles+2; i++)
      xdr_template(&xdrs, &eps);
  xdr_destroy(&xdrs);
  fclose(outfile);

  fh.dimensions = 3;
  
  outfile = fopen("dark/vel", "wb");
  xdrstdio_create(&xdrs, outfile, XDR_ENCODE);

  xdr_template(&xdrs, &fh);
  Vector3D<float> vel(0.0);
  for(i = 0; i < fh.numParticles+2; i++)
      xdr_template(&xdrs, &vel);
  xdr_destroy(&xdrs);
  fclose(outfile);

  outfile = fopen("dark/pos", "wb");
  xdrstdio_create(&xdrs, outfile, XDR_ENCODE);

  xdr_template(&xdrs, &fh);
  Vector3D<float> pos(xmin);
  xdr_template(&xdrs, &pos);
  pos = xmax;
  xdr_template(&xdrs, &pos);

  deltax = (xmax - xmin);
  for(i = 0; i < fh.numParticles; i++) {
      pos.x = xmin + rand()/((double) RAND_MAX)*deltax;
      pos.y = xmin + rand()/((double) RAND_MAX)*deltax;
      pos.z = xmin + rand()/((double) RAND_MAX)*deltax;
      xdr_template(&xdrs, &pos);
      }
  
  xdr_destroy(&xdrs);
  fclose(outfile);
  
  return(0);
}
