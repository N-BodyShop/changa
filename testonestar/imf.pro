pro imf

IMF = ' '
read, 'Enter type of IMF used: MS for Miller-Scalo, or K for Kroupa: ', IMF
 
get_lun, lun  
openw, lun, 'imf.dat'
printf, lun, IMF
close, lun

end
 

