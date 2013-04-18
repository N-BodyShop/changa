pro rtipsy,file,header,catg,catd,cats,TIME = time,VERBOSE = verbose, JUSTHEAD=justhead
;;; RTIPSY:  Tipsy reader for IDL
;;; Author:  James Wadsley
;;; 
if (N_PARAMS() eq 0) then begin
  print, "rtipsy.pro  Reads tipsy files detecting the format: "
  print, "big endian, little endian, padded (standard) or non-padded header "
  print
  print, "Usage: "
  print, "        rtipsy, filename ,header [,g] [,d] [,s] [,TIME=time] [,/VERBOSE]"
  print
  print, "Input parameters: "
  print, "  filename  filename string"
  print, "  time      desired output time (optional)"
  print, "  /VERBOSE  print messages (optional)"
  print, "Return values:"
  print, "  header    tipsy header struct"
  print, "  g,d,s     gas, dark and star structures"
  print, "Please read rtipsy.pro for the structure definitions"
  print
  print, "Example: "
  print, "  rtipsy, '/home/wadsley/usr5/mihos/mihos.std',h,g,d"
  print, "  print, h.ndark"
  print, "  plot, d.x, d.y, psym=3"
  return
endif

;;; Note: IDL structures are never paddded 
header = { time:double(0.0), n:0L, ndim:0L, ngas:0L, ndark:0L, nstar:0L }

close,1
openr,1,file

Loop:  

readu,1,header
endianswap = 0
if (header.ndim lt 1 or header.ndim gt 3) then begin
  endianswap = 1
  header = swap_endian(header)
  if (keyword_set(verbose)) then print,"SWAP_ENDIAN"
endif

if (keyword_set(verbose)) then print,"Read time,n,ngas,nstar,ndark: ",header.time,header.n,header.ngas,header.ndark,header.nstar

fs = fstat(1)
;;; Explicitly pad header if required 
if (fs.size eq 32L+header.ngas*48+header.ndark*36+header.nstar*44) then begin
  dummy = 1L
  readu,1,dummy
endif else if (fs.size ne 28L+header.ngas*48+header.ndark*36+header.nstar*44) then begin  
  print, "RTIPSY ERROR: Header and file size inconsistent"
  print, "Estimates: Header bytes:  28 or 32 (either is OK)"
  print, "     ngas: ",header.ngas," bytes:",48*header.ngas
  print, "    ndark: ",header.ndark," bytes:",36*header.ndark
  print, "    nstar: ",header.nstar," bytes:",44*header.nstar
  print, "Actual File bytes:",fs.size,"  not one of:",32L+header.ngas*48+header.ndark*36+header.nstar*44,28L+header.ngas*48+header.ndark*36+header.nstar*44
  close,1
  return
endif

IF (KEYWORD_SET(justhead) EQ 0) THEN BEGIN 
if (header.ngas ne 0) then begin 
	catg = replicate({mass: 1.,x: 1.,y : 1., z:1.,vx:1.,vy:1.,vz:1.,dens:1.,tempg:1.,h : 1. , zmetal : 1., phi : 1.},header.ngas)
	readu,1,catg
    if (endianswap eq 1) then catg=swap_endian(catg)
endif
if (header.ndark ne 0) then begin
    catd = replicate({mass: 1.,x: 1.,y : 1., z:1.,vx:1.,vy:1.,vz:1.,eps: 1.,phi: 1.},header.ndark)
    readu,1,catd
    if (endianswap eq 1) then catd=swap_endian(catd)
endif
if (header.nstar ne 0) then begin
    cats = replicate({mass: 1.,x: 1.,y : 1., z:1.,vx:1.,vy:1.,vz:1.,metals:1.,tform:1.,eps: 1.,phi: 1.},header.nstar)
    readu,1,cats
    if (endianswap eq 1) then cats=swap_endian(cats)
endif
  
;;; Loop over output times if requested
if (keyword_set(time)) then begin
  if (abs(time-header.time) gt 1e-3) then begin
    on_ioerror, ReadError
    goto, Loop
  endif
endif

close,1
ENDIF ;end the justhead if statement
return

ReadError:
print,"RTIPSY ERROR: Output time not found ",time
on_ioerror,NULL

close,1
return

end
