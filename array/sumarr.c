#include <stdio.h>
#include <math.h>

#define LINE_SIZE 256

int main(void)
{
	char ach[LINE_SIZE];
	int a;
	double da,dsum;

	fgets(ach, LINE_SIZE, stdin);
	dsum = 0.0;
	while (1) {
		a = scanf("%lg",&da);
		if (a != 1) break;
		dsum += da;
		}
	printf("%.17g\n",dsum);
	return 0;
	}







