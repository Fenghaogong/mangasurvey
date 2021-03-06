;+
;
; function mlfluxcal_nostds
;
; Try to roughly flux calibrate the data in the case where there are
; no standards available.  This won't be dreadfully good, as it won't
; account for any actual observing conditions, and just uses an
; ideal calibration vector based on the BOSS throughput model
; used for the MaNGA simulator.
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 2/18/2013
;   Last modified: 2/24/2013
;
; REVISION HISTORY:
;   v1: 18-Feb-2013  D. Law
;       First written.
;   v1.1: 24-Feb-2013 D. Law
;       Fixed inverse variance scaling
;-

function mlfluxcal_nostds,fileb1,filer1,Cfile1,obsparam,slitmap,waveset
  nwave=(size(waveset))[1]

  ; Read all data from files
  b1sci=mrdfits(fileb1,0)
  b1sciivar=mrdfits(fileb1,1)
  b1mask=mrdfits(fileb1,2)
  b1wset=mrdfits(fileb1,3); Array of legendre coefficients
  b1dispwset=mrdfits(fileb1,4)
  b1sky=mrdfits(fileb1,6)
  b1superflat=mrdfits(fileb1,8)

  r1sci=mrdfits(filer1,0)
  r1sciivar=mrdfits(filer1,1)
  r1mask=mrdfits(filer1,2)
  r1wset=mrdfits(filer1,3); Array of legendre coefficients
  r1dispwset=mrdfits(filer1,4)
  r1sky=mrdfits(filer1,6)
  r1superflat=mrdfits(filer1,8)

  ; Multiply by superflat
  ; (this takes out the DIVISION by superflat applied in extract_object)
  b1sci=b1sci*b1superflat
  r1sci=r1sci*r1superflat

  ; Convert from legendre coefficient to CCD wavelength solution
  traceset2xy,b1wset,dummy,b1loglam; Convert from coefficients to log wavelengths
  b1lam=10.^b1loglam; Array of wavelengths for each fiber
  traceset2xy,r1wset,dummy2,r1loglam; Convert from coefficients to log wavelengths
  r1lam=10.^r1loglam; Array of wavelengths for each fiber

  ; Convert dispersion from traceset
  traceset2xy,b1dispwset,dummy3,b1disp
  traceset2xy,r1dispwset,dummy4,r1disp

  ; Put blue and red on a common wavelength grid.
  ; Note that BOSS routine rebin_spectrum changes units of spectrum,
  ; while interpol does not.  Therefore if flux units were erg/s/cm2/Ang
  ; originally, they still will be after running interpol.  They will
  ; not after running rebin_spectrum

  tempspecb1=fltarr(nwave)
  tempspecr1=fltarr(nwave)
  tempspecb1ivar=fltarr(nwave)
  tempspecr1ivar=fltarr(nwave)
  tempspecb1mask=fltarr(nwave)
  tempspecr1mask=fltarr(nwave)
  tempspecb1disp=fltarr(nwave)
  tempspecr1disp=fltarr(nwave)
  tempspecb1sky=fltarr(nwave)
  tempspecr1sky=fltarr(nwave)

  nfiber=n_elements(slitmap)
  finalflux=fltarr(nwave,nfiber)
  finalivar=fltarr(nwave,nfiber)
  finalmask=fltarr(nwave,nfiber)
  finaldispersion=fltarr(nwave,nfiber)
  finalsky=fltarr(nwave,nfiber)
  finalwave=fltarr(nwave,nfiber)

  ; Define approximate calibration vector from the blue/red throughput
  ; values in the MaNGA/BOSS simulator calibration file
  calibfile=strcompress(getenv('MANGACORE_DIR')+'/'+getenv('MANGACORE_VER')+'/bosscal/56280/calibmatrix.fits',/remove_all
  ; Check whether calibration file exists, return error code -120 if not
  junk=findfile(calibfile,count=ct)
  if (ct EQ 0) then return,-120L
  calibmatrix=readfits(calibfile)
  bluecalvec=dblarr(nwave)
  redcalvec=dblarr(nwave)
  bluecalvec[*]=calibmatrix[4,*]
  redcalvec[*]=calibmatrix[5,*]

  ; Avoid divide by zero errors at edges
  bluecalvec[where(bluecalvec le 0.)]=1.
  redcalvec[where(redcalvec le 0.)]=1.

  for j=0,nfiber-1 do begin
    tempspecb1=interpol(b1sci[*,j],b1lam[*,j],waveset)/obsparam.exptime/bluecalvec
    tempspecr1=interpol(r1sci[*,j],r1lam[*,j],waveset)/obsparam.exptime/redcalvec

    tempspecb1ivar=interpol(b1sciivar[*,j],b1lam[*,j],waveset)*(obsparam.exptime*bluecalvec)^2
    tempspecr1ivar=interpol(r1sciivar[*,j],r1lam[*,j],waveset)*(obsparam.exptime*redcalvec)^2

    tempspecb1mask=interpol(b1mask[*,j],b1lam[*,j],waveset)
    tempspecr1mask=interpol(r1mask[*,j],r1lam[*,j],waveset)

    tempspecb1disp=interpol(b1disp[*,j],b1lam[*,j],waveset)
    tempspecr1disp=interpol(r1disp[*,j],r1lam[*,j],waveset)

    tempspecb1sky=interpol(b1sky[*,j],b1lam[*,j],waveset)/obsparam.exptime/bluecalvec
    tempspecr1sky=interpol(r1sky[*,j],r1lam[*,j],waveset)/obsparam.exptime/redcalvec

    ; Break blue and red cameras at channel 2326
    ktrans=2326;
    tempspecr1[0:ktrans]=0.
    tempspecb1[ktrans:nwave-1]=0.
    tempspecr1ivar[0:ktrans]=0.
    tempspecb1ivar[ktrans:nwave-1]=0.
    tempspecr1mask[0:ktrans]=0
    tempspecb1mask[ktrans:nwave-1]=0
    tempspecr1disp[0:ktrans]=0.
    tempspecb1disp[ktrans:nwave-1]=0.
    tempspecr1sky[0:ktrans]=0.
    tempspecb1sky[ktrans:nwave-1]=0.

    ; Quick version just adds blue and red
    finalflux[*,j]=tempspecb1+tempspecr1
    finalivar[*,j]=tempspecb1ivar+tempspecr1ivar
    finalmask[*,j]=tempspecb1mask+tempspecr1mask
    finaldispersion[*,j]=tempspecb1disp+tempspecr1disp
    finalsky[*,j]=tempspecb1sky+tempspecr1sky
  endfor
  finalwave=waveset

   outname = Cfile1
   len=strlen(outname)
   suffix=strmid(outname,len-3,3)
   ; If filename ended in .gz suffix, remove it.
   ; Will be added back when we gzip later
   if (suffix eq '.gz') then outname=strmid(outname,0,len-3)

   ; HDU #0 is flux
   sxaddpar, bighdr, 'BUNIT', '1E-17 erg/cm^2/s/Ang'
   mwrfits, finalflux, outname, bighdr, /create

   ; HDU #1 is inverse variance
   sxaddpar, hdrfloat, 'BUNIT', '1/(1E-17 erg/cm^2/s/Ang)^2'
   mwrfits, finalivar, outname, hdrfloat

   ; HDU #2 is AND-pixelmask
   mwrfits, finalmask, outname, hdrlong

   ; HDU #3 is OR-pixelmask
   mwrfits, finalmask, outname, hdrlong

   ; HDU #4 is dispersion map
   sxaddpar, hdrfloat, 'BUNIT', 'pixels'
   mwrfits, finaldispersion, outname, hdrfloat

   ; HDU #5 is slitmap
   mwrfits, slitmap, outname

   ; HDU #6 is the sky
   mwrfits, finalsky, outname

   ; HDU #7 is the wavelength solution
   mwrfits, finalwave, outname

; gzip output file
spawn, ['gzip', '-f', outname], /noshell

end
