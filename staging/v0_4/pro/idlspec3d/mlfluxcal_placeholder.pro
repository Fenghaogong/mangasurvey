; function mlfluxcal_placeholder
;
; Returns 0 if everything ok, returns an error code if there was a problem.
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 01/23/2013
;   Last modified: 01/23/2013
;
; REVISION HISTORY:
;   v1: 23-Jan-2013  D. Law
;       Placeholder function written

function mlfluxcal_placeholder, fileb1, filer1, wave, spec,assignorder,sentry,ifusize

  ; Read all data from files
  b1sci=mrdfits(fileb1,0)
  b1wset=mrdfits(fileb1,3); Array of legendre coefficients
  r1sci=mrdfits(filer1,0)
  r1wset=mrdfits(filer1,3); Array of legendre coefficients

  hdr = headfits(fileb1)
nwave=(size(wave))[1]

  ; Convert from legendre coefficient to CCD wavelength solution
  traceset2xy,b1wset,dummy,b1loglam; Convert from coefficients to log wavelengths
  b1lam=10.^b1loglam; Array of wavelengths for each fiber
  traceset2xy,r1wset,dummy2,r1loglam; Convert from coefficients to log wavelengths
  r1lam=10.^r1loglam; Array of wavelengths for each fiber

; Note that BOSS routine rebin_spectrum changes units of spectrum,
; while interpol does not.  Therefore if flux units were erg/s/cm2/Ang
; originally, they still will be after running interpol.  They will
; not after running rebin_spectrum

  tempspecb1=fltarr(nwave)
  tempspecr1=fltarr(nwave)

  for j=0,ifusize-1 do begin
    tempspecb1=interpol(b1sci[*,sentry[assignorder[j]].fiberid-1],b1lam[*,sentry[assignorder[j]].fiberid-1],wave)
    tempspecr1=interpol(r1sci[*,sentry[assignorder[j]].fiberid-1],r1lam[*,sentry[assignorder[j]].fiberid-1],wave)

    tempspecr1[0:2283]=0.
    tempspecb1[2284:nwave-1]=0.

  ; Quick version just adds blue and red
  spec[*,j]=tempspecb1+tempspecr1
  endfor

; Determine output filename and write files
spcframefile=str_replace(fileb1,'spSFrame-b1-','spCFrame-')
len=strlen(spcframefile)
suffix=strmid(spcframefile,len-3,3)
; If old file ended in .gz suffix, remove it
if (suffix eq '.gz') then spcframefile=strmid(spcframefile,0,len-3)

mwrfits, spec, spcframefile, hdr, /create;sky subtracted flux

return,0
end
