
@dummy.pro

pro findmaxmass

maxstars = [15.,20.,25.,30.,40.,80.,100.,200.]


for i=0,n_elements(maxstars)-1 do begin
    maxstar = maxstars[i]
    N_SNIa = qromb('snia_K', 1.5, 8.0)
    M_SNIa_ej = N_SNIa*1.40
    Fe_SNII = qromb('snii_fe_K', 8.0, maxstar)
    Ox_SNII = qromb('snii_ox_K', 8.0, maxstar)
    LinFe_SNII = qromb('snii_linfe_K', 8.0, maxstar)
    LinOx_SNII = qromb('snii_linox_K', 8.0, maxstar)
    Fe_SNIa = N_SNIa*0.63
    Ox_SNIa = N_SNIa*0.13
    
    iron = Fe_SNII + Fe_SNIa
    ox = Ox_SNII + Ox_SNIa
    liniron = LinFe_SNII + Fe_SNIa
    linox = LinOx_SNII + Ox_SNIa

    print,maxstar,'  [O/Fe]=',alog10(ox/iron)-alog10(9.6/1.17), $
      '  [O/Fe]_II=',alog10(Ox_SNII/Fe_SNII)-alog10(9.6/1.17)
    print,'Raiteri Fe:',iron,'  Ox:',ox
    print,'Lin fit Fe:',liniron,'  Ox:',linox
endfor

end
