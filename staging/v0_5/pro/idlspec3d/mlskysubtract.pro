;+
; function mlskysubtract
;
; This is the high-level sky subtraction routine, analagous to what is
; done in BOSS extract_object.pro.  It calls mlcalcsky to actually
; figure out the sky values.;
; Returns 0 if everything ok, returns an error code if there was a problem.
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 01/28/2013
;   Last modified: 02/27/2013
;
; REVISION HISTORY:
;   v1.0: 28-Jan-2013  D. Law
;       Imported code from BOSS extract_object.pro, start integrating
;       with MaNGA algorithms and calling
;   v1.1: 05-Feb-2013  D. Law
;       Revised call to require output filename
;   v1.2: 15-Feb-2013  D. Law
;       Revised to require keyword 'wave', and handle 2'', 3'', 5''
;       fibers seperately.  Outputs extra extensions to spSFrame
;       where everything is rectified to the same wave scale.
;   v1.3: 27-Feb-2013 D. Law
;       Added ximg and superflat to output extensions again.  Now
;       mimics spFrame except with slitmap in place of plugmap.
;-

function mlskysubtract, spframefile, spsframefile, obsparam, slitmap, camera, wave, visual=visual

;;;;;;;;;;;;;;;;;;;;;;;;;;

; Read in files
mlframe_read,spframefile,objflux=flux,objivar=fluxivar,wset=vacset,mask=pixelmask, $
   dispset=dispset,ximg=ximg,superflat=superflat,hdr=objhdr

nx = (size(flux,/dim))[0] 
ny = (size(flux,/dim))[1] 

; Set up the blue/red frame parameters
if (camera eq 'b1') then nbkpt = 3*nx/4
if (camera eq 'r1') then nbkpt = nx

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Define which fibers are sky fibers of 2, 3, 5'' flavor

ftype2=where((slitmap.fsize eq 2.))
ftype3=where((slitmap.fsize eq 3.))
ftype5=where((slitmap.fsize eq 5.))

; Split up structures accordingly
flux2=flux[*,ftype2]
flux3=flux[*,ftype3]
flux5=flux[*,ftype5]

fluxivar2=fluxivar[*,ftype2]
fluxivar3=fluxivar[*,ftype3]
fluxivar5=fluxivar[*,ftype5]

pixelmask2=pixelmask[*,ftype2]
pixelmask3=pixelmask[*,ftype3]
pixelmask5=pixelmask[*,ftype5]

traceset2xy,vacset,junk,loglam
loglam2=loglam[*,ftype2]
loglam3=loglam[*,ftype3]
loglam5=loglam[*,ftype5]

vacset2=traceset_trim(vacset,ftype2)
vacset3=traceset_trim(vacset,ftype3)
vacset5=traceset_trim(vacset,ftype5)

traceset2xy,dispset,junk,dispval
dispval2=dispval[*,ftype2]
dispval3=dispval[*,ftype3]
dispval5=dispval[*,ftype5]

slitmap2=slitmap[ftype2]
slitmap3=slitmap[ftype3]
slitmap5=slitmap[ftype5]

; Define where the skies are within each subgroup
iskies2=where((slitmap2.ifuname eq 'SKY2') AND (slitmap2.plugstatus eq 1), nskies2)
iskies3=where((slitmap3.ifuname eq 'SKY3') AND (slitmap3.plugstatus eq 1), nskies3)
iskies5=where((slitmap5.ifuname eq 'SKY5') AND (slitmap5.plugstatus eq 1), nskies5)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; First pass
tai_mid=double(obsparam.taistart)+obsparam.exptime/2.

; Do 2'' fibers
skystruct2 = mlcalcsky(flux2, fluxivar2, vacset2, slitmap2, $
  skysub2, skysubivar2, iskies=iskies2, pixelmask=pixelmask2, $
  fibermask=fibermask2, upper=3.0, lower=3.0, tai=tai_mid, nbkpt=nbkpt)

if (NOT keyword_set(skystruct2)) then begin
  splog,'Problem with skysubtract2- quit!'
  mlquitmanga3d,-200L 
endif

; Do 3'' fibers
skystruct3 = mlcalcsky(flux3, fluxivar3, vacset3, slitmap3, $
  skysub3, skysubivar3, iskies=iskies3, pixelmask=pixelmask3, $
  fibermask=fibermask3, upper=3.0, lower=3.0, tai=tai_mid, nbkpt=nbkpt)

if (NOT keyword_set(skystruct3)) then begin
  splog,'Problem with skysubtract3- quit!'
  mlquitmanga3d,-200L 
endif

; Do 5'' fibers
skystruct5 = mlcalcsky(flux5, fluxivar5, vacset5, slitmap5, $
  skysub5, skysubivar5, iskies=iskies5, pixelmask=pixelmask5, $
  fibermask=fibermask5, upper=3.0, lower=3.0, tai=tai_mid, nbkpt=nbkpt)

if (NOT keyword_set(skystruct5)) then begin
  splog,'Problem with skysubtract5- quit!'
  mlquitmanga3d,-200L 
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; If any of the sky-fibers are bad, then re-do sky-subtraction.

; Do 2'' fibers
ibadfib2 = where(djs_median(skysub2[*,iskies2]^2 * $
  skysubivar2[*,iskies2], 1) GT 2.0)               
if (ibadfib2[0] NE -1) then begin
  fibermask2[iskies2[ibadfib2]] = fibermask2[iskies2[ibadfib2]] OR fibermask_bits('BADSKYFIBER')

  splog, 'Calling skysubtract2 again; masked skyfibers',string(iskies2[ibadfib2])
  skystruct2 = mlcalcsky(flux2, fluxivar2, vacset2, slitmap2, $
    skysub2, skysubivar2, iskies=iskies2, pixelmask=pixelmask2, $
     fibermask=fibermask2, upper=10.0, lower=10.0, tai=tai_mid, nbkpt=nbkpt)

  if (NOT keyword_set(skystruct2)) then begin
   splog,'Problem with skysubtract2- quit!'
   mlquitmanga3d,-200L 
  endif
endif

; Do 3'' fibers
ibadfib3 = where(djs_median(skysub3[*,iskies3]^2 * $
  skysubivar3[*,iskies3], 1) GT 2.0)               
if (ibadfib3[0] NE -1) then begin
  fibermask3[iskies3[ibadfib3]] = fibermask3[iskies3[ibadfib3]] OR fibermask_bits('BADSKYFIBER')

  splog, 'Calling skysubtract3 again; masked skyfibers',string(iskies3[ibadfib3])
  skystruct3 = mlcalcsky(flux3, fluxivar3, vacset3, slitmap3, $
    skysub3, skysubivar3, iskies=iskies3, pixelmask=pixelmask3, $
     fibermask=fibermask3, upper=10.0, lower=10.0, tai=tai_mid, nbkpt=nbkpt)

  if (NOT keyword_set(skystruct3)) then begin
   splog,'Problem with skysubtract3- quit!'
   mlquitmanga3d,-200L 
  endif
endif

; Do 5'' fibers
ibadfib5 = where(djs_median(skysub5[*,iskies5]^2 * $
  skysubivar5[*,iskies5], 1) GT 2.0)               
if (ibadfib5[0] NE -1) then begin
  fibermask5[iskies5[ibadfib5]] = fibermask5[iskies5[ibadfib5]] OR fibermask_bits('BADSKYFIBER')

  splog, 'Calling skysubtract5 again; masked skyfibers',string(iskies5[ibadfib5])
  skystruct5 = mlcalcsky(flux5, fluxivar5, vacset5, slitmap5, $
    skysub5, skysubivar5, iskies=iskies5, pixelmask=pixelmask5, $
     fibermask=fibermask5, upper=10.0, lower=10.0, tai=tai_mid, nbkpt=nbkpt)

  if (NOT keyword_set(skystruct5)) then begin
   splog,'Problem with skysubtract5- quit!'
   mlquitmanga3d,-200L 
  endif
endif

; QA plots for chi^2 from 1D sky-subtraction.
;if keyword_set(VISUAL) then getwindow,/open
qaplot_skysub, flux2, fluxivar2, skysub2, skysubivar2, vacset2, iskies2, title=' 1D Sky-subtraction'
qaplot_skysub, flux3, fluxivar3, skysub3, skysubivar3, vacset3, iskies3, title=' 1D Sky-subtraction'
qaplot_skysub, flux5, fluxivar5, skysub5, skysubivar5, vacset5, iskies5, title=' 1D Sky-subtraction'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;----------
; Sky-subtract one final time, this time with dispset (PSF subtraction)
; (rejected sky fibers from above remain rejected).
; Modify pixelmask in this call.
; DRL - these results looks BAD.  Turn it off for now.
;nskypoly = 3L
;skystruct = mlcalcsky(flux, fluxivar, vacset, slitmap, $
;  skysub, skysubivar, iskies=iskies, pixelmask=pixelmask, $
;  fibermask=fibermask, upper=10.0, lower=10.0, tai=tai_mid, $
;  ; dispset=dispset, $ ; Why is this commented out??
;  npoly=nskypoly, nbkpt=nbkpt, $
;  relchi2set=relchi2set, newmask=newmask)
;pixelmask = newmask

;if (NOT keyword_set(skystruct)) then begin
;  splog,'Problem with skysubtract- quit!'
;  mlquitmanga3d,-200L 
;endif

; QA plots for chi^2 from 2D sky-subtraction.
;if keyword_set(VISUAL) then getwindow,/open
;qaplot_skysub, flux, fluxivar, skysub, skysubivar, $
;  vacset, iskies, title=' 2D Sky-subtraction'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Re-assemble the different fiber sizes back into the same structures
skysub=fltarr(nx,ny)
skysubivar=fltarr(nx,ny)
pixelmask=fltarr(nx,ny)
fibermask=fltarr(ny)

skysub[*,ftype2]=skysub2[*,*]
skysub[*,ftype3]=skysub3[*,*]
skysub[*,ftype5]=skysub5[*,*]

skysubivar[*,ftype2]=skysubivar2[*,*]
skysubivar[*,ftype3]=skysubivar3[*,*]
skysubivar[*,ftype5]=skysubivar5[*,*]

pixelmask[*,ftype2]=pixelmask2[*,*]
pixelmask[*,ftype3]=pixelmask3[*,*]
pixelmask[*,ftype5]=pixelmask5[*,*]

fibermask[ftype2]=fibermask2[*]
fibermask[ftype3]=fibermask3[*]
fibermask[ftype5]=fibermask5[*]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Save the sky-subtracted flux values as is, and now modify flambda.
flambda = skysub
flambdaivar = skysubivar
skyimg = flux - flambda

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; DRL- telluric correction would be here, but not done yet



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;----------
; Interpolate over masked pixels, just for aesthetic purposes

flambda = djs_maskinterp(flambda, flambdaivar EQ 0, /const, iaxis=0 )

;----------
; Combine FIBERMASK and PIXELMASK to FINALMASK
finalmask = pixelmask
ntrace=ny ;???
for itrace=0, ntrace-1 do $
  finalmask[*,itrace] = finalmask[*,itrace] OR fibermask[itrace]

;----------
; Disable some mask bits in regions where 'NODATA' is set
q_nodata = (finalmask AND sdss_flagval('SPPIXMASK','NODATA')) NE 0
discards = ['NEARBADPIXEL','LOWFLAT','SCATTEREDLIGHT','NOSKY']
for j=0, n_elements(discards)-1 do $
  finalmask = finalmask - q_nodata $
  * (finalmask AND sdss_flagval('SPPIXMASK',discards[j]))

;----------
; Get an estimate of the relative chi^2 at each pixel.
; Do this with a simple linear interpolation.
if (keyword_set(relchi2set)) then begin
  xx = 0
  traceset2xy, vacset, xx, loglam
;   rchi2img = interpol(relchi2struct.chi2, relchi2struct.wave, loglam)
  rchi2img = bspline_valu(loglam, relchi2set)
   ; Compute the mean relative chi2 of sky-subtraction, after masking
   ; bad regions of the CCD
  fval = sdss_flagval('SPPIXMASK','NOPLUG') $
    + sdss_flagval('SPPIXMASK','BADTRACE') $
    + sdss_flagval('SPPIXMASK','BADFLAT') $
    + sdss_flagval('SPPIXMASK','BADARC') $
    + sdss_flagval('SPPIXMASK','LOWFLAT') $
    + sdss_flagval('SPPIXMASK','NOSKY') $
    + sdss_flagval('SPPIXMASK','NODATA') $
    + sdss_flagval('SPPIXMASK','BADFLUXFACTOR')
  indx = where((finalmask AND fval) EQ 0 AND flambdaivar NE 0, ct)
  if (ct EQ 0) then skychi2 = 0. $
  else skychi2 = mean(rchi2img[indx])
endif else begin
  rchi2img = 0 * flambda + 1.
  skychi2 = 0.
endelse
sxaddpar, objhdr, 'SKYCHI2', skychi2, ' Mean chi^2 of sky-subtraction'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


sxaddpar, objhdr, 'MANGADRP_VER', getenv('MANGADRP_VER')

; Determine output details from provided spsframefile name
outname=spsframefile
len=strlen(outname)
suffix=strmid(outname,len-3,3)
; If filename ended in .gz suffix, remove it.
; Will be added back when we gzip later
if (suffix eq '.gz') then outname=strmid(outname,0,len-3)

mwrfits, flambda, outname, objhdr, /create;sky subtracted flux
mwrfits, flambdaivar, outname   ; sky subtracted inverse variance
mwrfits, finalmask, outname ; final pixel mask
mwrfits, vacset, outname    ;trace-set for wavelength sol; wset
mwrfits, dispset, outname   ;trace-set for dispersion sol
mwrfits, slitmap, outname  ;slitmap used
mwrfits, skyimg, outname    ;sky flux
; Not sure what these last two are, just spit out same as was read in
mwrfits, ximg, outname      ;x pos on CCD
mwrfits, superflat, outname  ;superflat vector from quartz lamps

; Put everything on the same wavelength solution just to be helpful
; for Bershady analysis (not to be used for science)
lambdaber=10.^loglam
flamber=fltarr((size(wave))[1],ny)
flamivarber=fltarr((size(wave))[1],ny)
maskber=fltarr((size(wave))[1],ny)
wavearray=fltarr((size(wave))[1],ny)

for bershady=0,ny-1 do begin
  flamber[*,bershady]=interpol(flambda[*,bershady],lambdaber[*,bershady],wave)
  flamivarber[*,bershady]=interpol(flambdaivar[*,bershady],lambdaber[*,bershady],wave)
  maskber[*,bershady]=interpol(finalmask[*,bershady],lambdaber[*,bershady],wave)
  wavearray[*,bershady]=wave[*]
endfor
mwrfits, flamber, outname
mwrfits, flamivarber, outname
mwrfits, maskber, outname
mwrfits, wavearray, outname

; gzip output file
spawn, ['gzip', '-f', outname], /noshell

return,0
end
