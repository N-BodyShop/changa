; Authors: Alyson Brooks + Greg Stinson
; Last upate: Mar 2011
; Purpose: to check various aspects of the metal output from gasoline 

; Inputs:
; IMF: either Miller-Scale or Kroupa is used in gasoline currently. 
pro metalcheck

IMF = ' '
get_lun, lun
openr, lun, 'imf.dat'
readf, lun, IMF
close, lun

; Get metallicity of the one star
rtipsy, 'onestar.tbin', h,g,d,s
z = s[0].metals
maxstar = 40.0
minsecondary = 0.9

; find units from param file via eff.out file
effstring = ' '
openr, lun, 'eff.out',/get_lun
readf, lun, effstring
line = strsplit(effstring,/extract)
dDSF = float(line[n_elements(line)-1])
readf, lun, effstring
line = strsplit(effstring,/extract)
dKpcUnit = float(line[n_elements(line)-1])
readf, lun, effstring
line = strsplit(effstring,/extract)
dMsolUnit = float(line[n_elements(line)-1])
close, lun
;You have: 1 / (G sunmass /kpc^3)^0.5
;You want: yr
;	* 4.7141126e+11
eff = 4.7141126e11*dDSF/(dMsolUnit/dKpcUnit^3.0)^0.5

if z lt 7e-5 then z = 7e-5 
if z gt 3e-2 then z = 3e-2 
; Padova lifetimes
a0 = 10.13+0.07547*alog10(z)-0.008084*(alog10(z))^2.
a1 = -4.424-0.7939*alog10(z)-0.1187*(alog10(z))^2.
a2 = 1.262+0.3385*alog10(z)+0.05417*(alog10(z))^2.
; t = 10.^(a0+a1*alog10(mass)+a2*(alog10(mass))^2.)

;------------------------------------------------------------------
; Do the integrals relating to each IMF here.

if IMF eq 'MS' then goto, MS
if IMF eq 'K' then goto, Kroupa
if IMF eq 'K01' then goto, Kroupa01

;Kroupa IMF
; fit to M, not log M!
; Analytic solutions from Raiteri et al. (1996)
Kroupa: M_SNII_ej = qromb('ejecta_K', 8.0, maxstar)
N_SNIa = qromb('snia_K', minsecondary, 8.0)
M_SNIa_ej = N_SNIa*1.40
M_winds_ej = qromb('wind_ejecta_K',1.0,8.0)
Fe_SNII = qromb('snii_fe_K', 8.0, maxstar)
Ox_SNII = qromb('snii_ox_K', 8.0, maxstar)
Mg_SNII = qromb('snii_mglin_K', 8.0, maxstar)
Si_SNII = qromb('snii_silin_K', 8.0, maxstar)
Neon_SNII = qromb('snii_Nelin_K', 8.0, maxstar)
C_SNII = qromb('snii_clin_K', 8.0, maxstar)
N_SNII = qromb('snii_nlin_K', 8.0, maxstar)
Ox_agb = qromb('agb_oxlin_K', 1.0, 8.0)
Fe_SNIa = N_SNIa*0.63
Ox_SNIa = N_SNIa*0.13

goto, jump1

;Kroupa01 IMF
; fit to M, not log M!
; Analytic solutions from Raiteri et al. (1996)
Kroupa01: M_SNII_ej = qromb('ejecta_K01', 8.0, maxstar)
N_SNIa = qromb('snia_K01', minsecondary, 8.0)
M_SNIa_ej = N_SNIa*1.40
M_winds_ej = qromb('wind_ejecta_K01',1.0,8.0)
Fe_SNII = qromb('snii_felin_K01', 8.0, maxstar)
Ox_SNII = qromb('snii_oxlin_K01', 8.0, maxstar)
Mg_SNII = qromb('snii_mglin_K01', 8.0, maxstar)
Si_SNII = qromb('snii_silin_K01', 8.0, maxstar)
Neon_SNII = qromb('snii_Nelin_K01', 8.0, maxstar)
C_SNII = qromb('snii_clin_K01', 8.0, maxstar)
N_SNII = qromb('snii_nlin_K01', 8.0, maxstar)
Ox_agb = qromb('agb_oxlin_K01', 1.0, 8.0)
Fe_SNIa = N_SNIa*0.63
Ox_SNIa = N_SNIa*0.13

goto, jump1

;Miller-Scalo IMF
; this is for fit to M, not log M! 
MS: M_SNII_ej = qromb('ejecta_MS', 8.0, maxstar)
N_SNIa = qromb('snia_MS', minsecondary, 8.0)
M_SNIa_ej = N_SNIa*1.40
Fe_SNII = qromb('snii_fe_ms', 8.0, maxstar)
Ox_SNII = qromb('snii_ox_ms', 8.0, maxstar)
Fe_SNIa = N_SNIa*0.63
Ox_SNIa = N_SNIa*0.13

jump1: !quiet = 1
; Compare to output from simulations
rdfloat, 'SNII.out', massII, EII, metalsII
rdfloat, 'SNIa.out', massIa, EIa, metalsIa
rdfloat, 'winds.out', massw, Ew, metalsw
m_snII_sim = total(massII)*dMsolUnit
m_snIa_sim = total(massIa)*dMsolUnit
m_winds_sim = total(massw)*dMsolUnit
Fe_winds = m_winds_sim * z * 0.13
Ox_winds = m_winds_sim * z * 0.58
; Compare mass of ejecta
print, 'Checking ejected mass:'
print, '     expected (Msun)','   simulation   ',' Difference (%)'
ejecta_II_diff = abs(((m_snII_sim - M_SNII_ej)/M_SNII_ej)*100)
print, 'SNII:',M_SNII_ej, m_snII_sim,  ejecta_II_diff
ejecta_Ia_diff = abs(((m_snIa_sim - M_SNIa_ej)/M_SNIa_ej)*100) 
print, 'SNIa:',M_SNIa_ej, m_snIa_sim, ejecta_Ia_diff
print, 'winds:',M_winds_ej, m_winds_sim, abs((m_winds_sim - M_winds_ej)/M_winds_ej)*100
; Compare yields
print, ""
print, 'Checking yields:'
file = 'onestar.000130'
fe = read_ascii_array(file+'.FeMassFrac')
ox = read_ascii_array(file+'.OxMassFrac')
rtipsy, file, h,g,d,s
mass_fe = total(g.mass*fe)*1.665e9
mass_ox = total(g.mass*ox)*1.665e9
if (file_test(file+'.CMassFrac')) then begin
   C = read_ascii_array(file+'.CMassFrac')
   N = read_ascii_array(file+'.NMassFrac')
   Neon = read_ascii_array(file+'.NeMassFrac')
   Mg = read_ascii_array(file+'.MgMassFrac')
   Si = read_ascii_array(file+'.SiMassFrac')
   mass_C = total(g.mass*C)*1.665e9
   mass_N = total(g.mass*N)*1.665e9
   mass_Neon = total(g.mass*Neon)*1.665e9
   mass_Mg = total(g.mass*Mg)*1.665e9
   mass_Si = total(g.mass*Si)*1.665e9
   moremetals=1
endif else moremetals=0

print, '  SNII pred (M_sun)',"     SNIa pred  ",'winds','   Simulation   ','Difference (%)'
ejecta_Fe_diff = abs(((Fe_SNII+Fe_SNIa+Fe_winds-mass_fe)/ $
                      (Fe_SNII+Fe_SNIa+Fe_winds))*100)
print, format='("Fe:",4(G12))', Fe_SNII, Fe_SNIa, Fe_winds, $
       mass_fe, ejecta_Fe_diff
ejecta_Ox_diff = abs(((Ox_SNII+Ox_SNIa+Ox_winds-mass_ox)/ $
                      (Ox_SNII+Ox_SNIa+Ox_winds))*100)
print, format='("Ox:",4(G12))', Ox_SNII, Ox_SNIa, Ox_winds, $
       mass_ox, ejecta_Ox_diff
if (moremetals) then begin
   print, 'C:', C_SNII, 'C_agb', mass_C, abs(((C_SNII-mass_C)/(C_SNII))*100)
   print, 'N:', N_SNII, 'N_agb', mass_N, abs(((N_SNII-mass_N)/(N_SNII))*100)
   print, 'Neon:',Neon_SNII, 'Neon_agb', mass_Neon, abs(((Neon_SNII-mass_Neon)/(Neon_SNII))*100)
   print, 'Mg:', Mg_SNII,'Mg_agb',mass_Mg,abs(((Mg_SNII-mass_Mg)/(Mg_SNII))*100)
   print, 'Si:', Si_SNII,'Si_agb', mass_Si, abs(((Si_SNII-mass_Si)/(Si_SNII))*100)
endif

;------------------------------------------------------------------
; Check that the timescales of ejecta match the analytic timescales

; Find first and last timesteps that emit SNII ejecta
nonzero = where(massII ne 0)
firstII = [min(nonzero)+1.]*eff
lastII = [max(nonzero)+1.]*eff
; Find first and last timesteps that emit SNIa ejecta
nonzero = where(massIa ne 0)
firstIa = [min(nonzero)+1.]*eff
lastIa = [max(nonzero)+1.]*eff
; Find first and last timesteps that emit wind ejecta
nonzero = where(massw ne 0)
firstw = [min(nonzero)+1.]*eff
lastw = [max(nonzero)+1.]*eff
; Compare to expected theoretical lifetimes
mass = maxstar
t = 10.^(a0+a1*alog10(mass)+a2*(alog10(mass))^2.)  ;Padua ischrones
mass = 8.
t_eight = 10.^(a0+a1*alog10(mass)+a2*(alog10(mass))^2.)
print, ""
print, 'Checking timescales:'
print, 'First exp [Myrs]    1st in sim     diff(%)    last exp   last in sim    diff(%)'
time_II_start = abs((((t/1e6)-(firstII/1e6))/(t/1e6))*100)
time_II_last = abs((((t_eight/1e6)-(lastII/1e6))/(t_eight/1e6))*100)
print, format='("SNII:", 6(G12.4))',t/1e6, firstII/1e6, time_II_start, t_eight/1e6, lastII/1e6, time_II_last
time_Ia_start = abs((((t_eight/1e6)-(firstIa/1e6))/(t_eight/1e6))*100)
mass = minsecondary
t_last_snia = 10.^(a0+a1*alog10(mass)+a2*(alog10(mass))^2.)
time_Ia_last = abs((((t_last_snia/1e9)-(lastIa/1e9))/(t_last_snia/1e9))*100)
print, format='("SNIa:", 4(G12.4),2(" Gyr",(G8.4)))', t_eight/1e6, firstIa/1e6, time_Ia_start, t_last_snia/1e9, lastIa/1e9, time_Ia_last
time_w_start = abs((((t_eight/1e6)-(firstw/1e6))/(t_eight/1e6))*100)
mass = 1.0
t = 10.^(a0+a1*alog10(mass)+a2*(alog10(mass))^2.)
time_w_last = abs((((t/1e9)-(lastw/1e9))/(t/1e9))*100)
print, format='("wind:", 4(G12.4),2(" Gyr",(G8.4)))', t_eight/1e6, firstw/1e6, time_w_start, t/1e9, lastw/1e9, time_w_last 
; Get percentages for final verdict.
; Yields first.
print, 'SUMMARY: YIELDS'
if (ejecta_II_diff lt 2 and ejecta_Ia_diff lt 2 and ejecta_Fe_diff lt 2 and ejecta_Ox_diff lt 2) then print, 'Congrats!  All the ejecta and metal yields in the simulation agree with theoretical predictions to less than 2%.  The code seems to be working.' 
if ejecta_II_diff gt 2 and ejecta_II_diff lt 5 then print, 'The mass of SNII ejecta varies from analytic predictions by 2-5%.  BEWARE! Something may be wrong.'
if ejecta_Ia_diff gt 2 and ejecta_Ia_diff lt 5 then print, 'The mass of SNIa ejecta varies from analytic predictions by 2-5%.  BEWARE! Something may be wrong.'
if ejecta_Fe_diff gt 2 and ejecta_Fe_diff lt 5 then print, 'The total mass of Fe ejecta varies from analytic predictions by 2-5%.  BEWARE! Something may be wrong.  It is uncertain if SNIa or SNII are the problem.'
if ejecta_Ox_diff gt 2 and ejecta_Ox_diff lt 5 then print, 'The total mass of Ox ejecta varies from analytic predictions by 2-5%.  BEWARE! Something may be wrong.  It is uncertain if SNIa or SNII are the problem.'
if ejecta_II_diff gt 5 then print, 'The mass of SNII ejecta varies from analytic predictions by more than 5%.  SOMETHING APPEARS TO BE WRONG WITH THE CODE!'
if ejecta_Ia_diff gt 5 then print, 'The mass of SNIa ejecta varies from analytic predictions by more than 5%.  SOMETHING APPEARS TO BE WRONG WITH THE CODE!'
if ejecta_Fe_diff gt 5 then print, 'The total mass of Fe ejecta varies from analytic predictions by more than 5%.  It is uncertain if SNIa or SNII are the problem.  SOMETHING APPEARS TO BE WRONG WITH THE CODE!'
if ejecta_Ox_diff gt 5 then print, 'The total mass of Ox ejecta varies from analytic predictions by more than 5%.  It is uncertain if SNIa or SNII are the problem.  SOMETHING APPEARS TO BE WRONG WITH THE CODE!'
; Now timescales.
print, 'SUMMARY: TIMESCALES'
if (time_II_start lt 5 and time_II_last lt 2 and time_Ia_start lt 2 and time_Ia_last lt 2 and time_w_start lt 2 and time_w_last lt 2) then print, 'Congrats!  All the ejecta timescales in the simulation agree with theoretical predictions to less than 5%.  The code seems to be working.' 
if time_II_start gt 5 and time_II_start lt 10 then print, 'The time that SNII ejecta is first produced in the simulations is different from analytic predictions by 5-10%.  BEWARE! Check the numbers above to see the discrepancy.'
if time_II_last gt 2 and time_II_last lt 5 then print, 'The time that SNII ejecta stops being produced in the simulations is different from analytic predictions by 2-5%.  BEWARE! Check the numbers above to see the discrepancy.'
if time_Ia_start gt 2 and time_Ia_start lt 5 then print, 'The time that SNIa ejecta is first produced in the simulations is different from analytic predictions by 2-5%.  BEWARE! Check the numbers above to see the discrepancy.'
if time_Ia_last gt 2 and time_Ia_last lt 5 then print, 'The time that SNIa ejecta stops being produced in the simulations is different from analytic predictions by 2-5%.  BEWARE! Check the numbers above to see the discrepancy.'
if time_w_start gt 2 and time_w_start lt 5 then print, 'The time that wind ejecta is first produced in the simulations is different from analytic predictions by 2-5%.  BEWARE! Check the numbers above to see the discrepancy.'
if time_w_last gt 2 and time_w_last lt 5 then print, 'The time that wind ejecta stops being produced in the simulations is different from analytic predictions by 2-5%.  BEWARE! Check the numbers above to see the discrepancy.'
if time_II_start gt 10 then print, 'The time that SNII ejecta is first produced in the simulations is different from analytic predictions more than 10%.  SOMETHING APPEARS TO BE WRONG WITH THE CODE!'
if time_Ia_start gt 5 then print, 'The time that SNIa ejecta is first produced in the simulations is different from analytic predictions more than 5%.  SOMETHING APPEARS TO BE WRONG WITH THE CODE!'
if time_w_start gt 5 then print, 'The time that wind ejecta is first produced in the simulations is different from analytic predictions more than 5%.  SOMETHING APPEARS TO BE WRONG WITH THE CODE!'
if time_II_last gt 5 then print, 'The time that SNII ejecta stops being produced in the simulations is different from analytic predictions more than 5%.  SOMETHING APPEARS TO BE WRONG WITH THE CODE!'
if time_Ia_last gt 5 then print, 'The time that SNIa ejecta stops being produced in the simulations is different from analytic predictions more than 5%.  SOMETHING APPEARS TO BE WRONG WITH THE CODE!'
if time_w_last gt 5 then print, 'The time that wind ejecta stops being produced in the simulations is different from analytic predictions more than 5%.  SOMETHING APPEARS TO BE WRONG WITH THE CODE!'

!quiet = 0

end

