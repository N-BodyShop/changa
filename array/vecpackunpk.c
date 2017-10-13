#include <stdio.h>
#include <math.h>
#include <assert.h>

struct PART 
{
    double vec[3];
    };

int main(int argc,char **argv)
{
	FILE *fpv,*fpi;
	int nVec;
	struct PART *part;
	int i, j;

	if(argc != 1) {
	    fprintf(stderr, "Usage: vecpackup < packed > unpacked\n");
	    return 1;
	    }
	scanf("%d",&nVec);
	part = malloc(nVec*sizeof(struct PART));
	for (i = 0; i < nVec; i++) {	/* read in packed order */
	    for(j = 0; j < 3; j++) {
		scanf("%lg", &(part[i].vec[j]));
		}
	    }
	printf("%d\n", nVec);
	for(j = 0; j < 3; j++) { 	/* write out in unpacked order */
	    for (i = 0; i < nVec; i++) {
		printf("%.15g\n", part[i].vec[j]);
		}
	    }
	return 0;
    }
