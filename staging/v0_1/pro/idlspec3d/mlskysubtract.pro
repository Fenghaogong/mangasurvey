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
;   Last modified: 01/28/2013
;
; REVISION HISTORY:
;   v1: 28-Jan-2013  D. Law
;       Imported code from BOSS extract_object.pro, start integrating
;       with MaNGA algorithms and calling
;-

function mlskysubtract, spframefile, obsparam, slitmap, camera, visual=visual

;;;;;;;;;;;;;;;;;;;;;;;;;;

; Read in files
mlframe_read,spframefile,objflux=flux,objivar=fluxivar,wset=vacset,mask=pixelmask, $
   dispset=dispset,hdr=objhdr,ximg=xnow,superflat=superfit

nx = (size(flux,/dim))[0] 
ny = (size(flux,/dim))[1] 

; Set up the blue/red frame parameters
if (camera eq 'b1') then nbkpt = 3*nx/4
if (camera eq 'r1') then nbkpt = nx

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Define which fibers are sky fibers.  If this isn't done here it's
; set to a default in mlcalcsky.pro.
; Anne-Marie: here is the easiest place to tweak this

iskies=where((slitmap.ifuname eq 'SKY2') AND (slitmap.plugstatus eq 1), nskies)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; First pass
tai_mid=double(obsparam.taistart)+obsparam.exptime/2.
skystruct = mlcalcsky(flux, fluxivar, vacset, slitmap, $
  skysub, skysubivar, iskies=iskies, pixelmask=pixelmask, $
  fibermask=fibermask, upper=3.0, lower=3.0, tai=tai_mid, nbkpt=nbkpt)

if (NOT keyword_set(skystruct)) then begin
  splog,'Problem with skysubtract- quit!'
  mlquitmanga3d,-200L 
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; If any of the sky-fibers are bad, then re-do sky-subtraction.
ibadfib = where(djs_median(skysub[*,iskies]^2 * $
  skysubivar[*,iskies], 1) GT 2.0)               
if (ibadfib[0] NE -1) then begin
  fibermask[iskies[ibadfib]] = fibermask[iskies[ibadfib]] OR fibermask_bits('BADSKYFIBER')

  splog, 'Calling skysubtract again; masked skyfibers',string(iskies[ibadfib])
  skystruct = mlcalcsky(flux, fluxivar, vacset, slitmap, $
    skysub, skysubivar, iskies=iskies, pixelmask=pixelmask, $
     fibermask=fibermask, upper=10.0, lower=10.0, tai=tai_mid, nbkpt=nbkpt)

  if (NOT keyword_set(skystruct)) then begin
   splog,'Problem with skysubtract- quit!'
   mlquitmanga3d,-200L 
  endif
endif

; QA plots for chi^2 from 1D sky-subtraction.
;if keyword_set(VISUAL) then getwindow,/open
qaplot_skysub, flux, fluxivar, skysub, skysubivar, vacset, iskies, title=' 1D Sky-subtraction'

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

if (NOT keyword_set(skystruct)) then begin
  splog,'Problem with skysubtract- quit!'
  mlquitmanga3d,-200L 
endif

; QA plots for chi^2 from 2D sky-subtraction.
;if keyword_set(VISUAL) then getwindow,/open
;qaplot_skysub, flux, fluxivar, skysub, skysubivar, $
;  vacset, iskies, title=' 2D Sky-subtraction'

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


sxaddpar, objhdr, 'MANGAVER', getenv('MANGAVER')

; Determine output filename and write files
spsframefile=str_replace(spframefile,'spFrame','spSFrame')
len=strlen(spsframefile)
suffix=strmid(spsframefile,len-3,3)
; If old file ended in .gz suffix, remove it
if (suffix eq '.gz') then spsframefile=strmid(spsframefile,0,len-3)

mwrfits, flambda, spsframefile, objhdr, /create;sky subtracted flux
mwrfits, flambdaivar, spsframefile   ; sky subtracted inverse variance
mwrfits, finalmask, spsframefile ; final pixel mask
mwrfits, vacset, spsframefile    ;trace-set for wavelength sol; wset
mwrfits, dispset, spsframefile   ;trace-set for dispersion sol
mwrfits, slitmap, spsframefile  ;slitmap used
mwrfits, skyimg, spsframefile    ;sky flux
; Not sure what these last two are, just spit out same as was read in
mwrfits, xnow, spsframefile      ;x pos on CCD
mwrfits, superfit, spsframefile  ;superflat vector from quartz lamps

; gzip output file
spawn, ['gzip', '-f', spsframefile], /noshell

return,0
end
