#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define LINE_SIZE 256

int main(void)
{
	char ach[LINE_SIZE];
	int a,n,i;
	double da;
	double *t;

	fgets(ach, LINE_SIZE, stdin);
	fputs(ach, stdout);
	sscanf(ach,"%d",&n);
	t = (double *)malloc(n*sizeof(double));
	if (!t) {
		fprintf(stderr, "Can't allocate memory for %d particles\n", n);
		return;
		}
	for (i=0;i<n;++i) {
		a = scanf("%lg",&da);
		if (a != 1) break;
		t[i] = da*da;
		}
	for (i=0;i<n;++i) {
		a = scanf("%lg",&da);
		if (a != 1) break;
		t[i] += da*da;
		}
	for (i=0;i<n;++i) {
		a = scanf("%lg",&da);
		if (a != 1) break;
		t[i] += da*da;
		printf("%.17g\n",sqrt(t[i]));
		}
	return 0;
	}
