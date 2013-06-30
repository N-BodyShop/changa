function snia_MS, m 
M_s = 459040.5 ;mass of my ssp (mass of star particle in M_sun)
A = 8.429e-3
frac_b = 0.16
 if 2.*m lt 3 then x = 3. else x=2.*m
; if m ge 10 then return, M_s*A*frac_b*(-240./4.3)*m^2.*((m+8.)^(-4.3)-x^(-4.3)) 
 return, M_s*A*frac_b*(-42./3.5)*m^2.*((m+8.)^(-3.5)-x^(-3.5))
end



