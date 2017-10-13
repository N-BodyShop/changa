; Author: Alyson Brooks
; Date: Jan 2006
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

if z lt 7d-5 then z = 7d-5 
if z gt 3d-2 then z = 3d-2 
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
N_SNIa = qromb('snia_K', 1.5, 8.0)
M_SNIa_ej = N_SNIa*1.40
Fe_SNII = qromb('snii_fe_K', 8.0, maxstar)
Ox_SNII = qromb('snii_ox_K', 8.0, maxstar)
Fe_SNIa = N_SNIa*0.63
Ox_SNIa = N_SNIa*0.13

goto, jump1

;Kroupa01 IMF
; fit to M, not log M!
Kroupa01: M_SNII_ej = qromb('ejecta_K01', 8.0, maxstar)
N_SNIa = qromb('snia_K01', 1.5, 8.0)
M_SNIa_ej = N_SNIa*1.40
Fe_SNII = qromb('snii_fe_K01', 8.0, maxstar)
Ox_SNII = qromb('snii_ox_K01', 8.0, maxstar)
Fe_SNIa = N_SNIa*0.63
Ox_SNIa = N_SNIa*0.13

goto, jump1

;Miller-Scalo IMF
; this is for fit to M, not log M! 
MS: M_SNII_ej = qromb('ejecta_MS', 8.0, maxstar)
N_SNIa = qromb('snia_MS', 1.5, 8.0)
M_SNIa_ej = N_SNIa*1.40
Fe_SNII = qromb('snii_fe_ms', 8.0, maxstar)
Ox_SNII = qromb('snii_ox_ms', 8.0, maxstar)
Fe_SNIa = N_SNIa*0.63
Ox_SNIa = N_SNIa*0.13

jump1: !quiet = 1
; Compare to output from simulations
res = sntots()
m_snII_sim = res[0]
m_snIa_sim = res[1]
m_wind_sim = res[2]
; Compare mass of ejecta
print, 'Checking the yields:'
print, 'Mass of total SNII ejecta expected is ', M_SNII_ej, ' M_sun.'
print, 'Mass of total SNII ejecta in simulation is', m_snII_sim, ' M_sun.'
print, 'Percent difference: ', abs(((m_snII_sim - M_SNII_ej)/M_SNII_ej)*100), ' %'
print, 'Mass of total SNIa ejecta expected is ', M_SNIa_ej, ' M_sun.'
print, 'Mass of total SNIa ejecta in simulation is', m_snIa_sim, ' M_sun.'
print, 'Percent difference: ', abs(((m_snIa_sim - M_SNIa_ej)/M_SNIa_ej)*100), ' %'
; Compare yields
fe = read_ascii_array('onestar.00320.FeMassFrac')
ox = read_ascii_array('onestar.00320.OxMassFrac')
rtipsy, 'onestar.00320', h,g,d,s
mass_fe = total(g.mass*fe)*1.665e9
mass_ox = total(g.mass*ox)*1.665e9

print, 'Total mass of Fe expected in SNII is ', Fe_SNII, ' M_sun.'
print, 'Total mass of Fe expected in SNIa is ', Fe_SNIa, ' M_sun.'
print, 'Total mass of Fe produced in simulation is ', mass_fe, ' M_sun.'
print, 'Percent difference: ', abs(((Fe_SNII+Fe_SNIa-mass_fe)/(Fe_SNII+Fe_SNIa))*100), ' %'
print, 'Total mass of Ox expected in SNII is ', Ox_SNII, ' M_sun.'
print, 'Total mass of Ox expected in SNIa is ', Ox_SNIa, ' M_sun.'
print, 'Total mass of Ox produced in simulation is ', mass_ox, ' M_sun.'
print, 'Percent difference: ', abs(((Ox_SNII+Ox_SNIa-mass_ox)/(Ox_SNII+Ox_SNIa))*100), ' %'


;------------------------------------------------------------------
; Check that the timescales of ejecta match the analytic timescales

; Assuming that onestar is run as I have created it, the current effective 
; dDeltaStarForm is 807347 yrs.  Adopt this for now, as there is no easy way 
; I can find to extract this info from the log file and read it into IDL.  If 
; anyone changes my param file, this number may be affected!  
; Check that the numbers match:
effstring = ' '
get_lun, lun
openr, lun, 'eff.out'
readf, lun, effstring
close, lun
pos = stregex(effstring,'[0-9]+ yrs',length=len)
eff = strmid(effstring, pos, len-4)
if long(eff) ne 807347 then read, 'dDeltaStarForm appears to have been modified.  Please enter the effective dDeltaStarForm used in yrs (see .log file to obtain this value): ', eff

rdfloat, 'SNII.out', massII, EII, metalsII
rdfloat, 'SNIa.out', massIa, EIa, metalsIa
rdfloat, 'winds.out', massw, Ew, metalsw
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
print, 'Check that the timescales of ejecta match the analytic timescales:'
print, 'First SNII ejecta expected at ', t/1d6, ' Myrs.'
print, 'First SNII ejecta in simulation occurs at ', long(firstII)/1d6, ' Myrs.'
time_II_start = abs((((t/1d6)-(firstII/1d6))/(t/1d6))*100)
print, 'Percent difference: ', time_II_start
mass = 8.
t = 10.^(a0+a1*alog10(mass)+a2*(alog10(mass))^2.)
print, 'Last SNII ejecta expected at ', t/1d6, ' Myrs.'
print, 'Last SNII ejecta in simulation occurs at ', long(lastII)/1d6, ' Myrs.'
time_II_last = abs((((t/1d6)-(lastII/1d6))/(t/1d6))*100)
print, 'Percent difference: ', time_II_last
print, 'First SNIa ejecta expected at ', t/1d6, ' Myrs.'
print, 'First SNIa ejecta in simulation occurs at ', long(firstIa)/1d6, ' Myrs.'
time_Ia_start = abs((((t/1d6)-(firstIa/1d6))/(t/1d6))*100)
print, 'Percent difference: ', time_Ia_start
print, 'First wind ejecta expected at ', t/1d6, ' Myrs.'
print, 'First wind ejecta in simulation occurs at ', long(firstw)/1d6, ' Myrs.'
time_w_start = abs((((t/1d6)-(firstw/1d6))/(t/1d6))*100)
print, 'Percent difference: ', time_w_start
mass = 1.5
t = 10.^(a0+a1*alog10(mass)+a2*(alog10(mass))^2.)
print, 'Last SNIa ejecta expected at ', t/1d9, ' Gyrs.'
print, 'Last SNIa ejecta in simulation occurs at ', long(lastIa)/1d9, ' Gyrs.'
time_Ia_last = abs((((t/1d9)-(lastIa/1d9))/(t/1d9))*100)
print, 'Percent difference: ', time_Ia_last
mass = 1.0
t = 10.^(a0+a1*alog10(mass)+a2*(alog10(mass))^2.)
print, 'Last wind ejecta expected at ', t/1d9, ' Gyrs.'
print, 'Last wind ejecta in simulation occurs at ', lastw/1d9, ' Gyrs.'
time_w_last = abs((((t/1d9)-(lastw/1d9))/(t/1d9))*100)
print, 'Percent difference: ', time_w_last 
; Get percentages for final verdict.
; Yields first.
print, 'SUMMARY: YIELDS'
ejecta_II_diff = abs(((m_snII_sim - M_SNII_ej)/M_SNII_ej)*100)
ejecta_Ia_diff = abs(((m_snIa_sim - M_SNIa_ej)/M_SNIa_ej)*100) 
ejecta_Fe_diff = abs(((Fe_SNII+Fe_SNIa-mass_fe)/(Fe_SNII+Fe_SNIa))*100)
ejecta_Ox_diff = abs(((Ox_SNII+Ox_SNIa-mass_ox)/(Ox_SNII+Ox_SNIa))*100)
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
if (time_II_start lt 2 and time_II_last lt 2 and time_Ia_start lt 2 and time_Ia_last lt 2 and time_w_start lt 2 and time_w_last lt 2) then print, 'Congrats!  All the ejecta timescales in the simulation agree with theoretical predictions to less than 2%.  The code seems to be working.' 
if time_II_start gt 2 and time_II_start lt 5 then print, 'The time that SNII ejecta is first produced in the simulations is different from analytic predictions by 2-5%.  BEWARE! Check the numbers above to see the discrepancy.'
if time_II_last gt 2 and time_II_last lt 5 then print, 'The time that SNII ejecta stops being produced in the simulations is different from analytic predictions by 2-5%.  BEWARE! Check the numbers above to see the discrepancy.'
if time_Ia_start gt 2 and time_Ia_start lt 5 then print, 'The time that SNIa ejecta is first produced in the simulations is different from analytic predictions by 2-5%.  BEWARE! Check the numbers above to see the discrepancy.'
if time_Ia_last gt 2 and time_Ia_last lt 5 then print, 'The time that SNIa ejecta stops being produced in the simulations is different from analytic predictions by 2-5%.  BEWARE! Check the numbers above to see the discrepancy.'
if time_w_start gt 2 and time_w_start lt 5 then print, 'The time that wind ejecta is first produced in the simulations is different from analytic predictions by 2-5%.  BEWARE! Check the numbers above to see the discrepancy.'
if time_w_last gt 2 and time_w_last lt 5 then print, 'The time that wind ejecta stops being produced in the simulations is different from analytic predictions by 2-5%.  BEWARE! Check the numbers above to see the discrepancy.'
if time_II_start gt 5 then print, 'The time that SNII ejecta is first produced in the simulations is different from analytic predictions more than 5%.  SOMETHING APPEARS TO BE WRONG WITH THE CODE!'
if time_Ia_start gt 5 then print, 'The time that SNIa ejecta is first produced in the simulations is different from analytic predictions more than 5%.  SOMETHING APPEARS TO BE WRONG WITH THE CODE!'
if time_w_start gt 5 then print, 'The time that wind ejecta is first produced in the simulations is different from analytic predictions more than 5%.  SOMETHING APPEARS TO BE WRONG WITH THE CODE!'
if time_II_last gt 5 then print, 'The time that SNII ejecta stops being produced in the simulations is different from analytic predictions more than 5%.  SOMETHING APPEARS TO BE WRONG WITH THE CODE!'
if time_Ia_last gt 5 then print, 'The time that SNIa ejecta stops being produced in the simulations is different from analytic predictions more than 5%.  SOMETHING APPEARS TO BE WRONG WITH THE CODE!'
if time_w_last gt 5 then print, 'The time that wind ejecta stops being produced in the simulations is different from analytic predictions more than 5%.  SOMETHING APPEARS TO BE WRONG WITH THE CODE!'

!quiet = 0

end

