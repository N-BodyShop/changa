/*
 ** see (A1) and (A2) of TREESPH: A UNIFICATION OF SPH WITH THE 
 ** HIERARCHICAL TREE METHOD by Lars Hernquist and Neal Katz.
 ** APJ Supplemant Series 70:416-446, 1989
 ** 
 ** Higher derivative terms c and d for use with quadrupole spline
 ** softening (Joachim Stadel, Dec. 94).
 */
inline
void SPLINEQ(double invr,double r2,double twoh,double& a,double& b,
	     double& c,double& d)
{
	double u,dih,dir=(invr);
	if ((r2) < (twoh)*(twoh)) {
		dih = 2.0/(twoh);
		u = dih/dir;
		if (u < 1.0) {
			a = dih*(7.0/5.0 - 2.0/3.0*u*u + 3.0/10.0*u*u*u*u
					 - 1.0/10.0*u*u*u*u*u);
			b = dih*dih*dih*(4.0/3.0 - 6.0/5.0*u*u + 1.0/2.0*u*u*u);
		    c = dih*dih*dih*dih*dih*(12.0/5.0 - 3.0/2.0*u);
			d = 3.0/2.0*dih*dih*dih*dih*dih*dih*dir;
			}
		else {
			a = -1.0/15.0*dir + dih*(8.0/5.0 - 4.0/3.0*u*u + u*u*u
			              - 3.0/10.0*u*u*u*u + 1.0/30.0*u*u*u*u*u);
			b = -1.0/15.0*dir*dir*dir + dih*dih*dih*(8.0/3.0 - 3.0*u + 6.0/5.0*u*u - 1.0/6.0*u*u*u);
			c = -1.0/5.0*dir*dir*dir*dir*dir + 3.0*dih*dih*dih*dih*dir
				+ dih*dih*dih*dih*dih*(-12.0/5.0 + 1.0/2.0*u);
			d = -dir*dir*dir*dir*dir*dir*dir
				+ 3.0*dih*dih*dih*dih*dir*dir*dir
					- 1.0/2.0*dih*dih*dih*dih*dih*dih*dir;
			}
		}
	else {
		a = dir;
		b = a*a*a;
		c = 3.0*b*a*a;
		d = 5.0*c*a*a;
		}
	}


inline void SPLINEM(double invr,double r2,double twoh,double& a,double& b)
{
	double u,dih,dir=(invr);
	if ((r2) < (twoh)*(twoh)) {
		dih = 2.0/(twoh);
		u = dih/dir;
		if (u < 1.0) {
			a = dih*(7.0/5.0 - 2.0/3.0*u*u + 3.0/10.0*u*u*u*u
					 - 1.0/10.0*u*u*u*u*u);
			b = dih*dih*dih*(4.0/3.0 - 6.0/5.0*u*u + 1.0/2.0*u*u*u);
			}
		else {
			a = -1.0/15.0*dir + dih*(8.0/5.0 - 4.0/3.0*u*u + u*u*u
			              - 3.0/10.0*u*u*u*u + 1.0/30.0*u*u*u*u*u);
			b = -1.0/15.0*dir*dir*dir + dih*dih*dih*(8.0/3.0 - 3.0*u + 6.0/5.0*u*u - 1.0/6.0*u*u*u);
			}
		}
	else {
		a = dir;
		b = a*a*a;
		}
	}


inline
void SPLINE(double r2, double twoh, double &a, double &b)
{
	double r, u,dih,dir;
	r = sqrt(r2);
	if (r < (twoh)) {
		dih = 2.0/(twoh);
		u = r*dih;
		if (u < 1.0) {
			a = dih*(7.0/5.0 - 2.0/3.0*u*u + 3.0/10.0*u*u*u*u
					 - 1.0/10.0*u*u*u*u*u);
			b = dih*dih*dih*(4.0/3.0 - 6.0/5.0*u*u + 1.0/2.0*u*u*u);
			}
		else {
			dir = 1.0/r;
			a = -1.0/15.0*dir + dih*(8.0/5.0 - 4.0/3.0*u*u + u*u*u
			              - 3.0/10.0*u*u*u*u + 1.0/30.0*u*u*u*u*u);
			b = -1.0/15.0*dir*dir*dir + dih*dih*dih*(8.0/3.0 - 3.0*u + 6.0/5.0*u*u - 1.0/6.0*u*u*u);
			}
		}
	else {
		a = 1.0/r;
		b = a*a*a;
		}
	}
