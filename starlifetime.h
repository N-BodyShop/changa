#ifndef STARTIME_H_INCLUDED
#define STARTIME_H_INCLUDED
/*
  Uses Raiteri, Villata, and Navarro (A&A, 315, 105, 1996) fit to 
  Padova group stellar models to find stellar lifetimes as a 
  function of mass, M,and metallicity, Z.

  log t = a0(Z) + a1(Z)*log(M) + a2(Z)*(log(M))**2

  where t is stellar lifetime in years, M is stellar mass 
  in solar masses, and Z, defined as the total mass fracion in
  elements heavier than He, is in the range 0.0004-0.05.  
  Coefficients are given as follows:

  a0(Z) = 10.13 + 0.07547*log(Z) - 0.008084*(log(Z))**2

  a1(Z) = -4.424 + 0.7939*log(Z) - 0.1187*(log(Z))**2

  a2(Z) = 1.262 + 0.3385*log(Z) - 0.1187*(log(Z))**2

  zmin, zmax = minimum and maximum metallicity in stellar lifetimes equation 
  zsol = solar metal abundance 
  xoxsol = solar oxygen abundance */

class Padova /*: public PUP::able*/ {
    double a00, a01, a02;
    double a10, a11, a12;
    double a20, a21, a22;
    double a0, a1, a2;    
    double zmin, zmax, zsol, xoxsol;    
    void CoefInit (double dMetals);
 public:
 Padova() : a00(10.13),	a01(0.07547), a02(-0.008084), 
      a10(-4.424),a11(-0.7939), a12(-0.1187),
      a20(1.262),a21(0.3385),a22(0.05417),
      a0(0.0), a1(0.0), a2(0.0),
      zmin(7e-5), zmax(3e-2), zsol(0.02), xoxsol(9e-3) {}

    double Lifetime(double dStarMass, double dMetals);
    double StarMass(double dStarLtime, double dMetals);
    /*PUPable_decl(Padova);
 Padova(CkMigrateMessage *m) : PUP::able(m) {}
    virtual*/ void pup(PUP::er& p) {
	/*      PUP::able::pup(p);*/
	p | a00;
	p | a01;
	p | a02;
	p | a10;
	p | a11;
	p | a12;
	p | a20;
	p | a21;
	p | a22;
	p | a0;
	p | a1;
	p | a2;
	p | zmin;
	p | zmax;
	p | zsol;
	p | xoxsol;
	}
    };


#endif
