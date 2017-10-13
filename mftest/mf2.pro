pro mf
!p.thick=3
!quiet=1

; For test purposes, this is a DM only volume.  I believe there is a 
; uniform resolution for these volumes.  The file .sovcirc contains 
; M_vir for each of the groups.  

rdfloat, 'cube300.00128.so.sovcirc', grp, mvir, rvir, skipline = 31 
volume = 300^3. ; This test run is 300 Mpc comoving on a side.
mvir = mvir(where(mvir GT 0 and rvir GT 0))

; Here I check if the halo is in the high res region.
; I think there is a uniform resolution for these volumes, though.     

; Set the number of bins; the x-axis is logarithmic in mass
binsize = 0.3   ; provides a fairly smooth curve in log(mass) 
upper_mass = max(alog10(mvir))
lower_mass = min(alog10(mvir))
nbins = (ceil(upper_mass)-floor(lower_mass))/binsize
bin=fltarr(nbins)
cum=fltarr(nbins)
mtest = alog10(mvir)-lower_mass  

; Get the indices to each m_vir in a given bin.
for i=0, n_elements(bin)-1 do begin
	ind1 = where(mtest GE (binsize*i) and mtest LT (binsize*(i+1)) )
	if ind1(0) EQ -1 then bin(i) = 0 else bin(i) = n_elements(ind1)
endfor
for i=0, n_elements(bin)-1 do cum(i) = total(bin[i:n_elements(bin)-1])
x = (indgen(n_elements(cum))*binsize)+lower_mass
y = alog10(cum/volume)

; Import values for previously determined mass function.  This MF has been 
; compared initially to PS expectations, and is now the standard to compare 
; future runs to, given that PS theory is slightly wrong.  
logm = [14.1744, 14.4744, 14.7744, 15.0744, 15.3744]
logn = [-5.1217, -5.5176, -6.0164, -6.4314, -7.4314]

test = spline(logm, logn, x) 

diff = abs(((test-y[0:n_elements(logm)-1])/test)*100.)
bad = where(diff GE 5.)

if bad EQ -1 then print, 'Congrats! This simulation produces a mass function as predicted!' 

if bad NE -1 then begin
print, 'PROBLEM!  This mass function differs from predicted by more than 5%!'
badmass = x(bad)
badn = y(bad)
pred = logn(bad)
print, 'At a value of Log(mass) = ', badmass, ' this simulation has Log (n) = ', badn, ' but expect Log(n) = ', pred
set_plot, 'ps'
device, filename='test.ps', /color
plot, x, y, xtitle = 'Log Mass (M_sun)', ytitle = 'Log (N/Vol) (Mpc^-3)' 
oplot, logm, logn, linestyle=2 
legend, ['Expected', 'Simulated'], linestyle=[1,2], /right
device, /close
print, 'See plot test.ps for visual comparison'

set_plot, 'x'
endif

end



