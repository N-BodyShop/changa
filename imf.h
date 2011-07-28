#ifndef IMF_HINCLUDED
#define IMF_HINCLUDED

// Please see 
// http://charm.cs.uiuc.edu/manuals/html/charm++/3_17.html#SECTION000317400000000000000
// for details about PUPing child and parent classes.

// http://www.parashift.com/c++-faq-lite/abcs.html#faq-22.5
// for copy constructors in inherited abstract classes and a
// description of pure virtual functions

class IMF : public PUP::able {

 public:
    IMF() {};
    virtual double returnimf(double mass) = 0;
    PUPable_abstract(IMF);
    IMF(CkMigrateMessage *m) : PUP::able(m) {}
    virtual double CumNumber(double mass) = 0;
    virtual double CumMass(double mass) = 0;
    virtual IMF* clone() const = 0;
    ~IMF() {};
};


/*
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
    PUPable_decl(MillerScalo);
    MillerScalo(CkMigrateMessage *m) : IMF(m) {}
    virtual double returnimf(double mass);
    virtual double CumNumber(double mass);
    virtual double CumMass(double mass);
    virtual MillerScalo* clone() const;
    virtual void pup(PUP::er &p) {
	PUP::able::pup(p);
	p|a1; p|b1; p|m1;
	p|a2; p|b2; p|m2;
	p|a3; p|b3; p|m3;
	p|mmax;
	}
    };

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
    PUPable_decl(Kroupa93);
    Kroupa93(CkMigrateMessage *m) : IMF(m) {}
    virtual double returnimf(double mass);
    virtual double CumNumber(double mass);
    virtual double CumMass(double mass);
    virtual Kroupa93* clone() const;
    virtual void pup(PUP::er &p) {
	PUP::able::pup(p);
	p|a1; p|b1; p|m1;
	p|a2; p|b2; p|m2;
	p|a3; p|b3; p|m3;
	p|mmax;
	}
};

/*
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
    virtual Chabrier* clone() const;
    PUPable_decl(Chabrier);
    Chabrier(CkMigrateMessage *m) : IMF(m) {}
    virtual void pup(PUP::er &p) {
	PUP::able::pup(p);
	p|a1; p|b1; p|m1;
	p|a2; p|b2; p|m2;
	p|mmax;
	}
};

#endif
