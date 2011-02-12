
#if defined(__cplusplus)
extern "C" {
#endif

double
dRombergO(void *CTX,double (*func)(void *, double),double a,double b,
	  double eps);
double
dRombergC(void *CTX,double (*func)(void *, double),double a,double b,
	  double eps);

#if defined(__cplusplus)
}

#endif
