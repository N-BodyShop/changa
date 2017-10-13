#include <stdio.h>
#include <math.h>
#include <assert.h>

struct PART 
{
    double vec[3];
    int iOrd;
    };

int compare(const struct PART *p1, const struct PART *p2)
{
    if(p1->iOrd < p2->iOrd)
	return -1;
    if(p1->iOrd > p2->iOrd)
	return 1;
    assert(0);			/* Two particles with the same iOrder? */
    return 0;
    }

int main(int argc,char **argv)
{
	FILE *fpv,*fpi;
	int nVec;
	int nOrd;
	struct PART *part;
	int i, j;

	if(argc != 3) {
	    fprintf(stderr, "Usage: sortarr vecfile indexfile\n");
	    return 1;
	    }
	fpv = fopen(argv[1],"r");
	fpi = fopen(argv[2],"r");
	fscanf(fpv,"%d",&nVec);
	fscanf(fpi,"%d",&nOrd);
	assert(nOrd == nVec);
	part = malloc(nOrd*sizeof(struct PART));
	for(j = 0; j < 3; j++) {
	    for (i = 0; i < nVec; i++) {
		fscanf(fpv, "%lg", &(part[i].vec[j]));
		}
	    }
	for (i = 0; i < nVec; i++) {
	    fscanf(fpi, "%d", &part[i].iOrd);
	    }
	
	qsort(part, nVec, sizeof(struct PART), compare);
	printf("%d\n", nVec);
	for(j = 0; j < 3; j++) {
	    for (i = 0; i < nVec; i++) {
		printf("%g\n", part[i].vec[j]);
		}
	    }
	return 0;
    }
