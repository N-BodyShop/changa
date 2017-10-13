/* 
 * Routine to take one "block" of a gadget file
 * Based on tipsy converter Obtained by trq from Rubert Croft via
 * Tiziana De Mateo.
 * Modified significantly by trq.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <endian.h>

/*
 * In place swap of data
 */
void swapEndian(void *data, int size, int count)
{
    char *cdata = (char *) data;
    int iCount;
    
    for(iCount = 0; iCount < count; iCount++) {
	int i;
	char temp;
	
	for(i = 0; i < size/2; i++) {
	    temp = cdata[size-i-1];
	    cdata[size-i-1] = cdata[i];
	    cdata[i] = temp;
	    }
	cdata += size;
	}
    }
	
struct io_header_1
{
  int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
} header1;


#define TYPE_GAS 0
#define TYPE_HALO 1
#define TYPE_DISK 2
#define TYPE_BULGE 3
#define TYPE_STAR 4
#define TYPE_BNDRY 5

int     NumPart, NumPartFiltered;


struct particle_data 
{
  float  Pos[3];		/* This could be any of the variables */
} *P;

void free_memory(void);
void allocate_memory(void);
/* returns dimensions (1 or 3) of data */
int load_snapshot(char *, int block_nr, int files);

/* Here we load one block of a snapshot file. It can be distributed
 * onto several files (for files>1).  It then gets output as an array file.
 */
int main(int argc, char **argv)
{
  char basename[200];
  int  type, files;
  int block_nr;
  int iDim;
  int i;
  int j;
  
  if(argc != 3 && argc != 4) 
    {
      fprintf(stderr,
	      "usage: snap2arr block snapbase [N_files]\n");
      exit(-1);
    }

  block_nr = atoi(argv[1]);
  strcpy(basename, argv[2]);
  
  if(argc == 4 ) 
    {
      files = atoi(argv[3]);	/* # of files per snapshot */
    }
  else 
    files = 1;

  iDim = load_snapshot(basename, block_nr, files);

  for(i=1; i<=NumPart; i++) {
      for(j = 0; j < iDim; j++) {
	  printf("%g\n", P[i].Pos[j]);
	  }
      }
  free_memory(); 

  return 0;
}


/* this routine loads particle data from one block of Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.
 */
int load_snapshot(char *fname, int block, int files)
{
  char   buf[200];
  int    i,k,dummy,ntot_withmasses;
  int    n,pc,pc_new,pc_sph;
  int nread;
  int swap = 0;
  int type;
  int iDim;
  FILE *fd;
  int iBlock;
  int nskip;
  

#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

  for(i=0, pc=1; i<files; i++, pc=pc_new)
    {
      if(files>1)
	sprintf(buf,"%s.%d",fname,i);
      else
	sprintf(buf,"%s",fname);

      if(!(fd=fopen(buf,"r")))
	{
	  fprintf(stderr,"can't open file `%s`\n",buf);
	  exit(0);
	}

      fread(&dummy, sizeof(dummy), 1, fd);
      if(dummy!=sizeof(header1)) {
	  swap = 1;
	  fprintf(stderr, "Trying endian swap\n");
	  swapEndian(&dummy, sizeof(dummy), 1);
	  assert(dummy == sizeof(header1));
	  }
      
      nread = fread(&header1, sizeof(header1), 1, fd);
      if(nread != 1) {
	fprintf(stderr, "Bad header read of %s\n", buf);
	exit(-1);
	}
      if(swap) {
	  swapEndian(&header1, 4, 6);           // 6 integers
	  swapEndian(&header1.mass, 8, 8);     // 8 doubles
	  swapEndian(&header1.flag_sfr, 4, 10);  // 10 more integers
	  swapEndian(&header1.BoxSize,8, 4);   // 4 more doubles
	  }
      
      fprintf(stderr, "BoxSize: %g\n", header1.BoxSize);
      fprintf(stderr, "Hubble parameter: %g\n", header1.HubbleParam);
      fread(&dummy, sizeof(dummy), 1, fd);

      NumPart = 0;
      if(files==1)
	{
	    for(type = 0; type < 6; type++)
		NumPart += header1.npart[type];
	}
      else
	{
	    for(type = 0; type < 6; type++)
		NumPart += header1.npartTotal[type];
	}


      if(i==0)
	allocate_memory();

      for(iBlock = 0; iBlock < block; iBlock++) {
	  SKIP;
	  nskip = dummy;
	  fseek(fd, nskip, SEEK_CUR);
	  SKIP;
	  assert(dummy == nskip);
	  }
      
      SKIP;

      if(dummy == NumPart*sizeof(float))
	  iDim = 1;
      else if(dummy == 3*NumPart*sizeof(float))
	  iDim = 3;
      else
	  assert(0);
      
      for(n=0;n<NumPart;n++)
	  {
	      if(iDim == 3) {
		  nread = fread(&P[pc_new].Pos[0], sizeof(float), 3, fd);
		  assert(nread == 3);
		  if(swap)
		      swapEndian(&P[pc_new].Pos[0], sizeof(float), 3);
		  }
	      if(iDim == 1) {
		  nread = fread(&P[pc_new].Pos[0], sizeof(float), 1, fd);
		  assert(nread == 1);
		  if(swap)
		      swapEndian(&P[pc_new].Pos[0], sizeof(float), 1);
		  }
	      pc_new++;
	      }
	}

      fclose(fd);
      return iDim;
    }

/* this routine allocates the memory for the 
 * particle data.
 */
void allocate_memory(void)
{
/*  fprintf(stderr,"allocating memory...\n");*/

  fprintf(stderr,"NumPart=%d \n",NumPart);

  if(!(P=malloc(NumPart*sizeof(struct particle_data))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
  
  P--;   /* start with offset 1 */

  
}


void free_memory(void)
{
  P++;
  free(P);
}
