function ejecta_K01, mass

M_s = 459040.5 ;mass of my ssp (mass of star particle in M_sun)
A = 0.2204
;phi_1 = M_s*A*(2^1.0)*mass^(-1.3)       ; 0.08 leq M leq 0.5
;phi_2 = M_s*A*mass^(-2.3)               ; 0.5 leq M leq 1
phi_3 = M_s*A*mass^(-2.3)               ; M geq 1

; For this code, never use values < 1 
return, 0.7682*mass^(1.056)*phi_3

end


