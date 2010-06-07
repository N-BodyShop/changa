pro rdfloat,name,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15,v16,v17, $
            v18,v19,SKIPLINE = skipline, NUMLINE = numline,DOUBLE=double, $
	    LONG=long
;+
; NAME:
;      RDFLOAT
; PURPOSE:
;      Quickly read a numeric ASCII data file into IDL floating pt. vectors.  
; EXPLANATION:
;      Columns of data may be separated by commas or spaces.      This 
;      program is fast but is restricted to data files where all columns can 
;      be read as floating point (or all double precision).   
;
;      Use READCOL if  greater flexibility is desired.   Use READFMT to read a 
;      fixed-format ASCII file.   Use FORPINT to print columns of data.
;
; CALLING SEQUENCE:
;      RDFLOAT, name, v1, [ v2, v3, v4, v5, ...  v19] 
;                         /DOUBLE, SKIPLINE = , NUMLINE = ]
;
; INPUTS:
;      NAME - Name of ASCII data file, scalar string.  In VMS, an extension of 
;              .DAT is assumed, if not supplied.
;
; OPTIONAL INPUT KEYWORDS:
;      SKIPLINE - Integer scalar specifying number of lines to skip at the top
;              of file before reading.   Default is to start at the first line.
;      NUMLINE - Integer scalar specifying number of lines in the file to read.  
;             Default is to read the entire file
;      /DOUBLE - If this keyword is set, then all variables are read in as
;              double precision.
;
; OUTPUTS:
;      V1,V2,V3,...V19 - IDL vectors to contain columns of data.
;               Up to 19 columns may be read.  All output vectors are of type
;               float, unless the /DOUBLE keyword is set, 
;
; EXAMPLES:
;      Each row in a file POSITION.DAT contains a star number and 6 columns
;      of data giving an RA and Dec in sexigesimal format.   Read into IDL 
;      variables.     
;
;       IDL> rdfloat,'POSITION',ID,hr,min,sec,deg,dmin,dsec  
;
;       All output vectors will be floating point
;
; RESTRICTIONS:
;      (1) All rows in the file must be formatted identically (except for 
;          those skipped by SKIPLINE).    RDFLOAT reads the first line of 
;          the data (after SKIPLINE) to determine the number of columns of 
;          data.
;      (2) Cannot be used to read strings
; PROCEDURES USED:
;      STR_SEP(), NUMLINES()
; REVISION HISTORY:
;      Written         W. Landsman                 September 1995
;      Call NUMLINES() function                    February 1996
;      Read up to 19 columns                       August 1997
;      Converted to IDL V5.0   W. Landsman         September 1997
;      Allow to skip more than 32767 lines  W. Landsman  June 2001
;-
  On_error,2                           ;Return to caller

  if N_params() lt 2 then begin
     print,'Syntax - RDFLOAT, name, v1, [ v2, v3,...v19 '
     print,'                    /DOUBLE, SKIPLINE =, NUMLINE = ]'
     return
  endif

; Get number of lines in file

   nlines = NUMLINES( name )
   if nlines LT 0 then return

   if not keyword_set( SKIPLINE ) then skipline = 0
   nlines = nlines - skipline
   if keyword_set( NUMLINE) then nlines = numline < nlines

;Read first line, and determine number of columns of data

   openr, lun, name, /GET_LUN
   temp = ''
   if skipline GT 0 then $
        for i=0L,skipline-1 do readf, lun, temp
   readf,lun,temp
   colval = strsplit( strtrim( strcompress(temp),2),' ')
   ncol = N_elements(colval)

;Create big output array and read entire file into the array

   if keyword_set(DOUBLE) then bigarr = dblarr(ncol, nlines, /NOZERO) $
   else if keyword_set(LONG) then bigarr = lonarr(ncol, nlines, /NOZERO) $ 
                          else bigarr = fltarr(ncol, nlines, /NOZERO) 

   close,lun
   openr, lun, name
   if skipline GT 0 then $
        for i=0L,skipline-1 do readf, lun, temp

   readf, lun, bigarr
   free_lun, lun

   message, strtrim(nlines,2) + ' lines of data read',/INF

   Nvector = (N_params()-1) < ncol
   v1 = reform( bigarr[0,*])

   if Nvector GT 1 then v2 = reform( bigarr[1,*]) else return
   if Nvector GT 2 then v3 = reform( bigarr[2,*]) else return
   if Nvector GT 3 then v4 = reform( bigarr[3,*]) else return
   if Nvector GT 4 then v5 = reform( bigarr[4,*]) else return
   if Nvector GT 5 then v6 = reform( bigarr[5,*]) else return
   if Nvector GT 6 then v7 = reform( bigarr[6,*]) else return
   if Nvector GT 7 then v8 = reform( bigarr[7,*]) else return
   if Nvector GT 8 then v9 = reform( bigarr[8,*]) else return
   if Nvector GT 9 then v10 = reform( bigarr[9,*]) else return
   if Nvector GT 10 then v11 = reform( bigarr[10,*]) else return
   if Nvector GT 11 then v12 = reform( bigarr[11,*]) else return
   if Nvector GT 12 then v13 = reform( bigarr[12,*]) else return
   if Nvector GT 13 then v14 = reform( bigarr[13,*]) else return
   if Nvector GT 14 then v15 = reform( bigarr[14,*]) else return
   if Nvector GT 15 then v16 = reform( bigarr[15,*]) else return
   if Nvector GT 16 then v17 = reform( bigarr[16,*]) else return
   if Nvector GT 17 then v18 = reform( bigarr[17,*]) else return
   if Nvector GT 18 then v19 = reform( bigarr[18,*]) 

  return
  end
