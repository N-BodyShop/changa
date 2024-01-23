
#if defined(__cplusplus)
extern "C" {
#endif

double
dRombergO(const void *CTX,double (*func)(const void *, double),double a,double b,
	  double eps);
double
dRombergC(const void *CTX,double (*func)(const void *, double),double a,double b,
	  double eps);

#if defined(__cplusplus)
}

#endif
