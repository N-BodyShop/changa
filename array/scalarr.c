/*
 * routine to scale a tipsy array by a constant factor.
 */
#include <stdio.h>
#include <stdlib.h>

int main(int argc,char **argv)
{
	char ach[80];
	int a;
	double da,dScale;
	
	dScale = atof(argv[1]);
	gets(ach);
	puts(ach);
	while (1) {
		a = scanf("%lg",&da);
		if (a != 1) break;  /* End of File */
		printf("%.17g\n",da*dScale);
		}
	return 0;
	}
