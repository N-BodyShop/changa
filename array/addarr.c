#include <stdio.h>
#include <math.h>

int main(int argc,char **argv)
{
	FILE *fpa,*fpb;
	char ach[80];
	int a,b;
	double da,db;

	fpa = fopen(argv[1],"r");
	fpb = fopen(argv[2],"r");
	fgets(ach,80,fpa);
	fgets(ach,80,fpb);
	fputs(ach,stdout);
	while (1) {
		a = fscanf(fpa,"%lg",&da);
		b = fscanf(fpb,"%lg",&db);
		if (a != 1 || b != 1) break;
		printf("%.17g\n",da+db);
		}
	fclose(fpa);
	fclose(fpb);
	return 0;
	}







