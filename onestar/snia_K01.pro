function snia_K01, m 
M_s = 459040.5 ;mass of my ssp (mass of star particle in M_sun)
A = 0.22038
frac_b = 0.16
 if 2.*m lt 3 then x = 3. else x=2.*m
 return, M_s*A*frac_b*(-1./3.3)*m^2.*((m+8.)^(-3.3)-x^(-3.3))
end

