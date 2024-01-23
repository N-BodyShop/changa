function snii_oxlin_K, mass

M_s = 459040.5 ;mass of my ssp (mass of star particle in M_sun)
A = 0.3029
;phi_1 = M_s*A*(2^0.9)*mass^(-1.3)       ; 0.08 leq M leq 0.5
;phi_2 = M_s*A*mass^(-2.2)               ; 0.5 leq M leq 1
phi_3 = M_s*A*mass^(-2.7)               ; M geq 1

; For this code, never use values < 1 
return, (0.210714*mass -2.0)*phi_3

end

function lin_K, mass, m,b

M_s = 459040.5 ;mass of my ssp (mass of star particle in M_sun)
A = 0.3029
;phi_1 = M_s*A*(2^0.9)*mass^(-1.3)       ; 0.08 leq M leq 0.5
;phi_2 = M_s*A*mass^(-2.2)               ; 0.5 leq M leq 1
phi_3 = M_s*A*mass^(-2.7)               ; M geq 1

; For this code, never use values < 1 
if (m*mass +b LT 0) then return, 0
return, (m*mass +b)*phi_3

end

function snii_mglin_K, m
  return,lin_K(m,0.01,-0.1)
end
function snii_silin_K, m
  return,lin_K(m,0.01,-0.1)
end
function snii_nelin_K, m
  return,lin_K(m,0.08,-2.0)
end
function snii_clin_K, m
  return,lin_K(m,0.01,0.0)
end
function snii_nlin_K, m
  return,lin_K(m,3.67e-3,-6.7e-3)
end
function agb_oxlin_K, m
  return,lin_K(m,0.0,0.001)
end
