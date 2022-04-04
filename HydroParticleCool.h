#ifndef HYDROPARTICLE_COOL_HH
#define HYDROPARTICLE_COOL_HH

#if defined(COOLING_MESA)
#ifdef USE_APPROX_21
#else
#define NSCALARS 1
#define iX NHYDROVARS
#endif
#elif defined(COOLING_SROEOS)
#define NSCALARS 1
#define iye NHYDROVARS
#else
#define NSCALARS 0
#endif
#endif
