#include <stdio.h>
#include <math.h>

#define LINE_SIZE 256

int main(void)
{
	char ach[LINE_SIZE];
	int a,cnt;
	double da,dsum;

	fgets(ach, LINE_SIZE, stdin);
	dsum = 0.0;
	cnt = 0;
	while (1) {
		a = scanf("%lg",&da);
		if (a != 1) break;
		dsum += da*da;
		++cnt;
		}
	if (cnt) printf("%.17g\n",sqrt(dsum/cnt));
	return 0;
	}







