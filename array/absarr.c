#include <stdio.h>
#include <math.h>

#define LINE_SIZE 256

int main(int argc,char **argv)
{
	char ach[LINE_SIZE];
	int a;
	double da;

	fgets(ach, LINE_SIZE, stdin);
	puts(ach);
	while (1) {
		a = scanf("%lg",&da);
		if (a != 1) break;
		printf("%.17g\n",fabs(da));
		}
	return 0;
	}







