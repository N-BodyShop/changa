;******************************************************************************
;This functions reads in ascii arrays output by the std_array_convert
;program.  The default output is float, but /DOUBLE or /LONG will output
;arrays of the appropriate type.
;******************************************************************************
FUNCTION read_ascii_array, filename,DOUBLE=double,LONG=long

openr, lun, filename,/get_lun
readf, lun, numlines
numlines=LONG(numlines)

;Get the array into the proper format
IF KEYWORD_SET(double) THEN BEGIN
    array=DBLARR(numlines) 
ENDIF ELSE BEGIN
    IF KEYWORD_SET(long) THEN BEGIN
        array=LONARR(numlines)
    ENDIF ELSE BEGIN
        array=FLTARR(numlines)
    ENDELSE
ENDELSE

readf, lun, array
close, lun
free_lun, lun

return, array

END

