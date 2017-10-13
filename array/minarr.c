#include <stdio.h>
#include <math.h>

#define LINE_SIZE 256

int main(void)
{
	char ach[LINE_SIZE];
	int a;
	double da,dmin;

	fgets(ach, LINE_SIZE, stdin);
	a = scanf("%lg",&da);
	if (a != 1) return;
	dmin = da;
	while (1) {
		a = scanf("%lg",&da);
		if (a != 1) break;
		if (da < dmin) dmin = da;
		}
	printf("%.17g\n",dmin);
	return 0;
	}







