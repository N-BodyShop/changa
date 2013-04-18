function ejecta_MS, mass

M_s = 459040.5 ;mass of my ssp (mass of star particle in M_sun)
A = 8.429e-3
;phi_1 = M_s*42*mass^(-1.4)                  ; 0.1 leq M leq 1
phi_2 = M_s*A*42*mass^(-2.5)                  ; 1 leq M leq 10
phi_3 = M_s*A*240*mass^(-3.3)                 ; M geq 10
; Note this is discontinuous at 10 M_sun, but this is what's in the paper.

; For this code, never use values < 1 
if mass ge 10 then return, 0.7682*mass^(1.056)*phi_3
if mass lt 10 then return, 0.7682*mass^(1.056)*phi_2

end


