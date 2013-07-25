;
; NAME:
;   getCommonMangaWave
;   
; PURPOSE:  
;   puts input arrays onto a common MANGA wavelength grid, provided by the calibration files 
;
; INPUTS:
;   wset - the traceset of the wavelength solution associated with input array
;   arr  - the input array
;
; OPTIONAL INPUTS:
;   arr2:arr6 - optional additional arrays of same size as arr
;   FILE      - a string indicating the fits filename to write the new 
;
; PROCEDURES CALLED:
;   mlsetwcalib
;
; REVISION HISTORY:
;   24-Jul-2013 written by B. Cherinka 
;

pro getCommonMangaWave, wset, arr, arr2, arr3, arr4, arr5, arr6, FILE=FILE

  on_error,0
  compile_opt idl2
  
  ;number of parameters passed
  npar = n_params()
  if npar lt 2 then begin
     print, 'At least two inputs are required!'
     return
  endif
  
  traceset2xy, wset, junk, loglam
  lam = 10.^loglam
  status = mlsetwcalib(wave)  ;Put everythong on a fiducial wavelength calibration from the MANGA calibration files
  maxfib = (size(arr))[2]     ;Grab total number of fibers
  wavearray = wave#(fltarr(maxfib)+1)
  
  ;create new arrays
  newarr = fltarr((size(wave))[1],maxfib, npar-1)
  
  ;interpolate arrays onto new wavelength grid
   for i=0,maxfib-1 do begin
        newarr[*,i,0] = interpol(arr[*,i], lam[*,i], wave)
        if npar eq 3 then newarr[*,i,1] = interpol(arr2[*,i], lam[*,i], wave)
        if npar eq 4 then newarr[*,i,2] = interpol(arr3[*,i], lam[*,i], wave)
        if npar eq 5 then newarr[*,i,3] = interpol(arr4[*,i], lam[*,i], wave)
        if npar eq 6 then newarr[*,i,4] = interpol(arr5[*,i], lam[*,i], wave)
        if npar eq 7 then newarr[*,i,5] = interpol(arr6[*,i], lam[*,i], wave)
   endfor    
  
  ;if filename specified, then write new arrays to the existing file
  if n_elements(file) ne 0 then begin
    mwrfits, newarr[*,*,0], file
    if npar gt 2 then for j=1, npar-2 do mwrfits, newarr[*,*,j], file
    mwrfits, wavewarray, file
  endif
  
  return

end

