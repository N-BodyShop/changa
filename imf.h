#ifndef IMF_HINCLUDED
#define IMF_HINCLUDED

// Please see 
// http://charm.cs.uiuc.edu/manuals/html/charm++/3_17.html#SECTION000317400000000000000
// for details about PUPing child and parent classes.

// http://www.parashift.com/c++-faq-lite/abcs.html#faq-22.5
// for copy constructors in inherited abstract classes and a
// description of pure virtual functions

/**
 * @brief Interface class for initial mass function.
 */

/**
* Modified by Elaad Applebaum to implement a stochastic IMF.
* Added DrawStar methods
* NB: So far, only added for Kroupa01
*/
class IMF : public PUP::able {

 public:
    IMF() {};
    /** @brief return stars per unit logarithmic mass
	@param mass in solar masses.
     */
    virtual double returnimf(double mass) = 0;
    /** @brief Charm++ requirement for passing polymorphic objects. */
    PUPable_abstract(IMF);
    /** @brief Charm++ migrate constructor */
    IMF(CkMigrateMessage *m) : PUP::able(m) {}
    /** @brief Charm++ Pack-UnPack method */
    virtual void pup(PUP::er &p) = 0;
    /** @brief Cumulative number of stars with mass greater than mass.
	@param mass in solar masses */
    virtual double CumNumber(double mass) = 0;
    /** @brief CumNumber for use with stochastic IMF.
    * NOTA BENE - CumNumber for stochastic use returns the actual number for
    * the star particle - no renormalization necessary
    */
    virtual double CumNumberStoch(double mass, double lownorm, double *hmstars, double cutmass) = 0;
    /** @brief Cumulative mass of stars with mass greater than mass.
	@param mass in solar masses */
    virtual double CumMass(double mass) = 0;
    /** @brief CumMass for use with stochastic IMF.
    * NOTA BENE - CumMass for stochastic use returns the actual mass for
    * the star particle - no renormalization necessary
    */
    virtual double CumMassStoch(double mass, double lownorm, double *hmstars, double cutmass) = 0;
    /** @brief inverse CDF of IMF to draw stars stochastically from the IMF */
    virtual double DrawStar(double num) = 0;
    /** @brief copy IMF object */
    virtual IMF* clone() const = 0;
    ~IMF() {};
};


/** @brief Implement Miller-Scalo IMF
    Uses the 3 segment power law fit for the Miller-Scalo IMF
    (Ap.J. Supp., 41,1979).

                                a1*(M**(b1))          0.1<M<1.
              IMF(Log10(M))=    a2*(M**(b2))         1.<M<10.
                                a3*(M**(b3))         10.<M
   Miller-Scalo IMF (Miller & Scalo, Ap.J. Supp., 41, 513, 1979) in
   stars per unit logarithmic mass.  Divide by M (mass) for IMF in
   stars per unit mass.  Also IMF is defined per yer per pc^2,
   integrated over a cylinder that extends "several hundred parsecs on
   either side of the plane of the galaxy"
*/

class MillerScalo : public IMF {

    double a1, b1, m1;
    double a2, b2, m2;
    double a3, b3, m3;
    double mmax;

 public:
    /* normalization, index, minimum mass */
    MillerScalo() {
	a1=42.0; b1=-0.4;m1=0.1; /* parameters from Ap.J. Supp., 41,1979 */
	a2=42.0;b2=-1.5;m2=1.0; /* This is discontinuous, but is what */
	a3=240.0;b3=-2.3;m3=10.0;/* they report in paper, so we leave it.*/
	mmax=100.0;
	} 
    /** @brief Charm++ method for migrating derived classes */
    PUPable_decl(MillerScalo);
    /** @brief Charm++ migration constructor */
    MillerScalo(CkMigrateMessage *m) : IMF(m) {}
    virtual double returnimf(double mass);
    virtual double CumNumber(double mass);
    virtual double CumMass(double mass);
    virtual double CumNumberStoch(double mass, double lownorm, double *hmstars, double cutmass);
    virtual double CumMassStoch(double mass, double lownorm, double *hmstars, double cutmass);
    virtual double DrawStar(double num);
    virtual MillerScalo* clone() const;
    virtual void pup(PUP::er &p) {
	PUP::able::pup(p);
	p|a1; p|b1; p|m1;
	p|a2; p|b2; p|m2;
	p|a3; p|b3; p|m3;
	p|mmax;
	}
    };

/**
 * @brief Implement IMF from Kroupa, Tout & Gilmore, 1993
 */
class Kroupa93 : public IMF {
    double a1, b1, m1;
    double a2, b2, m2;
    double a3, b3, m3;
    double mmax;
 public:
/* parameters from Raiteri et. al. A&A, 315,1996, eq. 2;  See also the
   conclusions of Kroupa, Tout & Gilmore, 1993. */
/* To convert to the IMF(log10(M)) convention of Miller-Scale, we
    increase the power law by 1 and multiply the coefficient by
    ln(10.0). See, eg., Chabrier 2003, eq. 2 */
    Kroupa93() {a1=0.3029*1.86606*log(10.0);b1=-0.3;m1=.08; 
	a2=0.3029*log(10.0);b2=-1.2;m2=0.5; 
	a3=0.3029*log(10.0); b3=-1.7; m3=1.0; 
	mmax=100.0; }
    /** @brief Charm++ method for migrating derived classes */
    PUPable_decl(Kroupa93);
    /** @brief Charm++ migration constructor */
    Kroupa93(CkMigrateMessage *m) : IMF(m) {}
    virtual double returnimf(double mass);
    virtual double CumNumber(double mass);
    virtual double CumMass(double mass);
    virtual double CumNumberStoch(double mass, double lownorm, double *hmstars, double cutmass);
    virtual double CumMassStoch(double mass, double lownorm, double *hmstars, double cutmass);
    virtual double DrawStar(double num);
    virtual Kroupa93* clone() const;
    virtual void pup(PUP::er &p) {
	PUP::able::pup(p);
	p|a1; p|b1; p|m1;
	p|a2; p|b2; p|m2;
	p|a3; p|b3; p|m3;
	p|mmax;
	}
};

/**
 * @brief Implement IMF from Kroupa 2001.
 */
class Kroupa01 : public IMF {
    double a1, b1, m1;
    double a2, b2, m2;
    double mmax;
 public:
/* parameters from Kroupa 2001, equation 2, and ignoring brown dwarfs,
   Also normalized so that the mass integral is 1. */
/* NOTE BENE: Kroupa 2001 has a revised IMF in section 6.2 which is
   different than this; however, below is what is used as the default in
   Starburst99
   (http://www.stsci.edu/science/starburst99/mappings/docs/run.html)
   with the exception that the low mass cutoff is .1 instead of the .08
   below and in the Kroupa paper.
 */
/* To convert to the IMF(log10(M)) convention of Miller-Scalo, we
    increase the power law by 1 and multiply the coefficient by
    ln(10.0). See, eg., Chabrier 2003, eq. 2 */
    Kroupa01() {
	a1=0.22038*2.0*log(10.0);b1=-0.3;m1=.08; 
	a2=0.22038*log(10.0);b2=-1.3;m2=0.5; 
	mmax=100.0; }
    /** @brief Charm++ method for migrating derived classes */
    PUPable_decl(Kroupa01);
    /** @brief Charm++ migration constructor */
    Kroupa01(CkMigrateMessage *m) : IMF(m) {}
    virtual double returnimf(double mass);
    virtual double CumNumber(double mass);
    virtual double CumMass(double mass);
    virtual double CumNumberStoch(double mass, double lownorm, double *hmstars, double cutmass);
    virtual double CumMassStoch(double mass, double lownorm, double *hmstars, double cutmass);
    virtual double DrawStar(double num);
    virtual Kroupa01* clone() const;
    virtual void pup(PUP::er &p) {
	PUP::able::pup(p);
	p|a1; p|b1; p|m1;
	p|a2; p|b2; p|m2;
	p|mmax;
	}
};

/** @brief Implement Chabrier IMF
  Use the log normal + power law fit of Chabrier, 2003, Galactic
  Stellar and Substellar Initial Mass Function", PASP 115, 763.
*/

class Chabrier : public IMF {
    double a1, b1, m1;
    double a2, b2, m2;
    double mmax;
 public:
    /*
      Chabrier low mass formula:
      \xi(log m) = A exp [ - (log m - log m_c)^2/2 \sigma^2]
      double a1, sigma, mc;  (b1 is sigma, m1 is mc)
    */
    Chabrier() {
	a1=0.158;b1=0.69;m1=.079;
	/* For high mass: normalization, index, minimum mass */
	/* Parameters from Table 1 of Chabrier, 2003. */
	a2=4.43e-2;b2=-1.3; m2=1.0;
	mmax=100.0;
	}
    virtual double returnimf(double mass);
    virtual double CumNumber(double mass);
    virtual double CumMass(double mass);
    virtual double CumNumberStoch(double mass, double lownorm, double *hmstars, double cutmass);
    virtual double CumMassStoch(double mass, double lownorm, double *hmstars, double cutmass);
    virtual double DrawStar(double num);
    virtual Chabrier* clone() const;
    /** @brief Charm++ method for migrating derived classes */
    PUPable_decl(Chabrier);
    /** @brief Charm++ migration constructor */
    Chabrier(CkMigrateMessage *m) : IMF(m) {}
    virtual void pup(PUP::er &p) {
	PUP::able::pup(p);
	p|a1; p|b1; p|m1;
	p|a2; p|b2; p|m2;
	p|mmax;
	}
};

#endif
