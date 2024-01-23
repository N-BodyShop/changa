function snia_K, m 
M_s = 459040.5 ;mass of my ssp (mass of star particle in M_sun)
A = 0.3029
frac_b = 0.16
 if 2.*m lt 3 then x = 3. else x=2.*m
 return, M_s*A*frac_b*(-1./3.7)*m^2.*((m+8.)^(-3.7)-x^(-3.7))
end

function snii_fe_K, mass

M_s = 459040.5 ;mass of my ssp (mass of star particle in M_sun)
A = 0.3029
;phi_1 = M_s*A*(2^0.9)*mass^(-1.3)       ; 0.08 leq M leq 0.5
;phi_2 = M_s*A*mass^(-2.2)               ; 0.5 leq M leq 1
phi_3 = M_s*A*mass^(-2.7)               ; M geq 1

; For this code, never use values < 1 
return, 2.802e-4*mass^(1.864)*phi_3

end


function snii_ox_K, mass

M_s = 459040.5 ;mass of my ssp (mass of star particle in M_sun)
A = 0.3029
;phi_1 = M_s*A*(2^0.9)*mass^(-1.3)       ; 0.08 leq M leq 0.5
;phi_2 = M_s*A*mass^(-2.2)               ; 0.5 leq M leq 1
phi_3 = M_s*A*mass^(-2.7)               ; M geq 1

; For this code, never use values < 1 
return, 4.586e-4*mass^(2.721)*phi_3

end

function snii_linfe_K, mass

M_s = 459040.5 ;mass of my ssp (mass of star particle in M_sun)
A = 0.3029
;phi_1 = M_s*A*(2^0.9)*mass^(-1.3)       ; 0.08 leq M leq 0.5
;phi_2 = M_s*A*mass^(-2.2)               ; 0.5 leq M leq 1
phi_3 = M_s*A*mass^(-2.7)               ; M geq 1

; For this code, never use values < 1 
return, (mass/300 + 0.017)*phi_3

end


function snii_linox_K, mass

M_s = 459040.5 ;mass of my ssp (mass of star particle in M_sun)
A = 0.3029
;phi_1 = M_s*A*(2^0.9)*mass^(-1.3)       ; 0.08 leq M leq 0.5
;phi_2 = M_s*A*mass^(-2.2)               ; 0.5 leq M leq 1
phi_3 = M_s*A*mass^(-2.7)               ; M geq 1

; For this code, never use values < 1 
return, (0.210714*mass - 2)*phi_3

end


