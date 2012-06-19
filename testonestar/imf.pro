pro imf

IMF = ' '
read, 'Enter type of IMF used: MS for Miller-Scalo, K01 for Kroupa01, or K for Kroupa: ', IMF
 
get_lun, lun  
openw, lun, 'imf.dat'
printf, lun, IMF
close, lun

end
 

