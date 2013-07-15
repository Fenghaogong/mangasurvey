;+
; function mlfluxcal
;
; Imports RSS-format data from the Phase 0 BOSS pipeline spFrame
; outputs into a MaNGA friendly set of vectors, doing proper
; flux calibration and joining of blue and red channels.
;
; This would use BOSS routines spframe_read.pro and spflux_v5.pro,
; but we've got to reimplement all of this to handle the fact that 
; we've got quite different rules and are using a slitmap .slm
; file in addition to the standard BOSS plPlugMapM file.  Code will
; crib heavily from the BOSS versions though.
;
; As per spframe_read.pro, note that extensions of spFrame are
;   HDU #0:  Flux
;   HDU #1:  Invvar
;   HDU #2:  mask
;   HDU #3:  wset
;   HDU #4:  dispset
;   HDU #5:  plugmap
;   HDU #6:  sky
;   HDU #7:  ximg
;   HDU #8:  superflat
;   HDU #9:  skystruct
;
; Returns 0 if everything ok, returns an error code if there was a problem.
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 12/12/2012
;   Last modified: 01/22/2013
;
; REVISION HISTORY:
;   v1: 12-Dec-2012  D. Law
;       First written.
;   v1.1: 30-Jan-2013 D. Law
;       Heavy revision for functionality
;-

;------------------------------------------------------------------------------
; Create a mask of 1's and 0's, where wavelengths that should not be used
; for fluxing (like near stellar features) are masked.
; 0 = not near lines, 1 = near lines
; HWIDTH = half width in log-wavelength for masking stellar lines

function spflux_masklines, loglam, hwidth=hwidth, stellar=stellar, $
 telluric=telluric

   if (NOT keyword_set(hwidth)) then $
    hwidth = 5.7e-4 ; Default is to mask +/- 5.7 pix = 400 km/sec

   mask = bytarr(size(loglam,/dimens))

   if (keyword_set(stellar)) then begin
      starwave = [ $
       3830.0 , $ ; ? (H-7 is at 3835 Ang)
       3889.0 , $ ; H-6
       3933.7 , $ ; Ca_k
       3968.5 , $ ; Ca_H (and H-5 at 3970. Ang)
       4101.7 , $ ; H-delta
       4300.  , $ ; G-band
       4305.  , $ ; G-band
       4310.  , $ ; more G-band
       4340.5 , $ ; H-gamma
       4861.3 , $ ; H-beta
       5893.0 , $ ; Mg
       6562.8 , $ ; H-alpha
       8500.8 , $
       8544.6 , $
       8665.0 , $
       8753.3 , $
       8866.1 , $
       9017.5 , $
       9232.0 ]
      airtovac, starwave

      for i=0L, n_elements(starwave)-1 do begin
         mask = mask OR (loglam GT alog10(starwave[i])-hwidth $
          AND loglam LT alog10(starwave[i])+hwidth)
      endfor
   endif

   if (keyword_set(telluric)) then begin
      tellwave1 = [6850., 7150., 7560., 8105., 8930.]
      tellwave2 = [6960., 7350., 7720., 8240., 9030.]
      for i=0L, n_elements(tellwave1)-1 do begin
         mask = mask OR (loglam GT alog10(tellwave1[i]) $
          AND loglam LT alog10(tellwave2[i]))
      endfor
   endif

   return, mask
end
;------------------------------------------------------------------------------

;------------------------------------------------------------------------------
; Divide the spectrum by a median-filtered spectrum.
; The median-filtered version is computed ignoring stellar absorp. features.

function spflux_medianfilt, loglam, objflux, objivar, width=width, $
 newivar=newivar, _EXTRA=KeywordsForMedian

   ndim = size(objflux, /n_dimen)
   dims = size(objflux, /dimens)
   npix = dims[0]
   if (ndim EQ 1) then nspec = 1 $
    else nspec = dims[1]

   ;----------
   ; Loop over each spectrum

   medflux = 0 * objflux
   if (arg_present(objivar)) then newivar = 0 * objivar
   for ispec=0L, nspec-1 do begin

      ; For the median-filter, ignore points near stellar absorp. features,
      ; but keep points near telluric bands.
      qgood = 1 - spflux_masklines(loglam[*,ispec], /stellar)

      ; Median-filter, but skipping masked points
      igood = where(qgood, ngood)
      thisback = fltarr(dims[0])
      if (ngood GT 1) then begin
         thisback[igood] = djs_median(objflux[igood,ispec], width=width, $
          _EXTRA=KeywordsForMedian)
      endif
      thisback = djs_maskinterp(thisback, (qgood EQ 0), /const)

      ; Force the ends of the background to be the same as the spectrum,
      ; which will force the ratio of the two to be unity.
      hwidth = ceil((width-1)/2.)
      thisback[0:hwidth] = objflux[0:hwidth,ispec]
      thisback[npix-1-hwidth:npix-1] = objflux[npix-1-hwidth:npix-1,ispec]
      czero2 = where(thisback eq 0., count2)
      if count2 gt 0 then thisback[czero2] = 1.
      medflux[*,ispec] = objflux[*,ispec] / thisback
      if (arg_present(objivar)) then $
      newivar[*,ispec] = objivar[*,ispec] * thisback^2
   endfor

   return, medflux
end
;------------------------------------------------------------------------------

;------------------------------------------------------------------------------
function spflux_bestmodel, loglam, objflux, objivar, dispimg, kindx=kindx1, $
 plottitle=plottitle

   filtsz = 99 ; ???
   cspeed = 2.99792458e5

   ndim = size(objflux, /n_dimen)
   dims = size(objflux, /dimens)
   npix = dims[0]
   if (ndim EQ 1) then nspec = 1 $
    else nspec = dims[1]

   ;----------
   ; Median-filter the object fluxes

   medflux = spflux_medianfilt(loglam, objflux, objivar, $
    width=filtsz, /reflect, newivar=medivar)
   sqivar = sqrt(medivar)

   ;----------
   ; Mask out the telluric bands

   sqivar = sqivar * (1 - spflux_masklines(loglam, /telluric))

   ;----------
   ; Load the Kurucz models into memory

   junk = spflux_read_kurucz(kindx=kindx)
   nmodel = n_elements(kindx)

   ;----------
   ; Fit the redshift just by using a canonical model

   ifud = where(kindx.teff EQ 6000 AND kindx.g EQ 4 AND kindx.feh EQ -1.5)
   if (ifud[0] EQ -1) then $
    message, 'Could not find fiducial model!'
   nshift = 20
   logshift = (-nshift/2. + findgen(nshift)) * 1.d-4
   chivec = fltarr(nshift)
   for ishift=0L, nshift-1 do begin
      modflux = spflux_read_kurucz(loglam-logshift[ishift], $
       dispimg, iselect=ifud)
      ; Median-filter this model
      medmodel = spflux_medianfilt(loglam, modflux, $
       width=filtsz, /reflect)
      for ispec=0L, nspec-1 do begin
         chivec[ishift] = chivec[ishift] + computechi2(medflux[*,ispec], $
          sqivar[*,ispec], medmodel[*,ispec])
      endfor
   endfor
   zshift = (10.d^logshift - 1) ; Convert log-lambda shift to redshift
   zpeak = find_nminima(chivec, zshift, errcode=errcode)
   splog, 'Best-fit velocity for std star = ', zpeak * cspeed, ' km/s'
   if (errcode NE 0) then $
    splog, 'Warning: Error code ', errcode, ' fitting std star'

   ;----------
   ; Generate the Kurucz models at the specified wavelengths + dispersions,
   ; using the best-fit redshift

   modflux = spflux_read_kurucz(loglam-alog10(1.+zpeak), dispimg)

   ;----------
   ; Loop through each model, computing the best chi^2
   ; as the sum of the best-fit chi^2 to each of the several spectra
   ; for this same object.
   ; We do this after a median-filtering of both the spectra + the models.

   chiarr = fltarr(nmodel,nspec)
   chivec = fltarr(nmodel)
   for imodel=0L, nmodel-1 do begin
      ; Median-filter this model
      medmodel = spflux_medianfilt(loglam, modflux[*,*,imodel], $
       width=filtsz, /reflect)

      for ispec=0L, nspec-1 do begin
         chiarr[imodel,ispec] = computechi2(medflux[*,ispec], $
          sqivar[*,ispec], medmodel[*,ispec])
      endfor
      chivec[imodel] = total(chiarr[imodel,*])
   endfor

   ;----------
   ; Return the best-fit model

   minchi2 = min(chivec, ibest)
   dof = total(sqivar NE 0)
   splog, 'Best-fit total chi2/DOF = ', minchi2/(dof>1)
   bestflux = modflux[*,*,ibest]

   ;----------
   ; Compute the chi^2 just around the stellar absorp. lines
   ; for the best-fit model star

   mlines = spflux_masklines(loglam, hwidth=12e-4, /stellar)
   linesqivar = sqivar * mlines
   linechi2 = 0.
   for ispec=0L, nspec-1 do begin
      thismodel = spflux_medianfilt(loglam, modflux[*,ispec,ibest], $
       width=filtsz, /reflect)
      linechi2 = linechi2 + computechi2(medflux[*,ispec], $
       linesqivar[*,ispec], thismodel)
   endfor
   linedof = total(linesqivar NE 0)
   splog, 'Best-fit line chi2/DOF = ', linechi2/(linedof>1)

   ;----------
   ; Compute the median S/N for all the spectra of this object,
   ; and for those data just near the absorp. lines

   sn_median = median(objflux * sqrt(objivar))
   indx = where(mlines, ct)
   if (ct GT 1) then $
    linesn_median = median(objflux[indx] * sqrt(objivar[indx])) $
   else $
    linesn_median = 0.
   splog, 'Full median S/N = ', sn_median
   splog, 'Line median S/N = ', linesn_median

   kindx1 = create_struct(kindx[ibest], $
    'IMODEL', ibest, $
    'Z', zpeak, $
    'SN_MEDIAN', sn_median, $
    'CHI2', minchi2, $
    'DOF', dof, $
    'LINESN_MEDIAN', linesn_median, $
    'LINECHI2', linechi2, $
    'LINEDOF', linedof)

   ;----------
   ; Plot the filtered object spectrum, overplotting the best-fit Kurucz model

   ; Select the observation to plot that has the highest S/N,
   ; and one that goes blueward of 4000 Ang.
   snvec = total(objflux * sqrt(objivar), 1) $
    * (10.^loglam[0,*] LT 4000 OR 10.^loglam[npix-1,*] LT 4000)
   junk = max(snvec, iplot) ; Best blue exposure

   snvec = total(objflux * sqrt(objivar), 1) $
    * (10.^loglam[0,*] GT 8600 OR 10.^loglam[npix-1,*] GT 8600)
   junk = max(snvec, jplot) ; Best red exposure

   csize = 0.85
   djs_plot, [3840., 4120.], [0.0, 1.4], /xstyle, /ystyle, /nodata, $
    xtitle='Wavelength [Ang]', ytitle='Normalized Flux', $
    title=plottitle
   if (iplot[0] NE -1) then begin
      djs_oplot, 10^loglam[*,iplot], medflux[*,iplot]
      djs_oplot, 10^loglam[*,iplot], medmodel[*,iplot], color='red'
   endif
   xyouts, 3860, 1.25, kindx1.model, charsize=csize
   djs_xyouts, 4000, 0.3, charsize=csize, $
    string(minchi2/(dof>1), format='("Total \chi^2/DOF=",f5.2)')
   djs_xyouts, 4000, 0.2, charsize=csize, $
    string(linechi2/(linedof>1), format='("Lines \chi^2/DOF=",f5.2)')
   djs_xyouts, 3860, 0.1, string(kindx1.feh, kindx1.teff, kindx1.g, $
    zpeak*cspeed, $
    format='("Fe/H=", f4.1, "  T_{eff}=", f6.0, "  g=", f3.1, "  cz=",f5.0)'), $
    charsize=csize

   djs_plot, [8440., 9160.], [0.0, 1.4], /xstyle, /ystyle, /nodata, $
    xtitle='Wavelength [Ang]', ytitle='Normalized Flux'
   if (jplot[0] NE -1) then begin
      djs_oplot, 10^loglam[*,jplot], medflux[*,jplot]
      djs_oplot, 10^loglam[*,jplot], medmodel[*,jplot], color='red'
   endif

   return, bestflux
end
;------------------------------------------------------------------------------



;------------------------------------------------------------------------------
function spflux_goodfiber, pixmask
   qgood = ((pixmask AND pixelmask_bits('NOPLUG')) EQ 0) $
       AND ((pixmask AND pixelmask_bits('BADTRACE')) EQ 0) $
       AND ((pixmask AND pixelmask_bits('BADFLAT')) EQ 0) $
       AND ((pixmask AND pixelmask_bits('BADARC')) EQ 0) $
       AND ((pixmask AND pixelmask_bits('MANYBADCOLUMNS')) EQ 0) $
       AND ((pixmask AND pixelmask_bits('NEARWHOPPER')) EQ 0) $
       AND ((pixmask AND pixelmask_bits('MANYREJECTED')) EQ 0)
   return, qgood
end

;------------------------------------------------------------------------------

;------------------------------------------------------------------------------
function spflux_bspline, loglam, mratio, mrativar, outmask=outmask, $
 everyn=everyn, airmass=airmass

   isort = sort(loglam)
   nord = 3

   ; Choose the break points using the EVERYN option, but masking
   ; out more pixels near stellar features just when selecting them.
   mask1 = 1 - spflux_masklines(loglam, hwidth=12.e-4, /stellar)
   ii = where(mrativar[isort] GT 0 AND mask1[isort] EQ 1)
   bkpt = 0
   fullbkpt = bspline_bkpts(loglam[isort[ii]], everyn=everyn, $
    bkpt=bkpt, nord=nord)

   outmask1 = 0
   if (keyword_set(airmass)) then begin
      x2 = airmass[isort]
   endif
   sset = bspline_iterfit(loglam[isort], mratio[isort], $
    invvar=mrativar[isort], lower=3, upper=3, fullbkpt=fullbkpt, $
    maxrej=ceil(0.05*n_elements(indx)), outmask=outmask1, nord=nord, $
    x2=x2, npoly=2*keyword_set(airmass), requiren=(everyn-1)>1)
   if (max(sset.coeff) EQ 0) then $
    message, 'B-spline fit failed!!'
   if (arg_present(outmask)) then begin
      outmask = bytarr(size(loglam,/dimens))
      outmask[isort] = outmask1
   endif

   return, sset
end

;------------------------------------------------------------------------------
function spflux_mratio_flatten, loglam1, mratio1, mrativar1, pres=pres

   ;--------
   ; Re-form the input data arrays from multi-dimensional to N x M

   ndim = size(loglam1, /n_dimen)
   dims = size(loglam1, /dimens)
   npix = dims[0]
   nobj = n_elements(loglam1) / npix
   loglam = reform(loglam1, npix, nobj)
   mratio = reform(mratio1, npix, nobj)
   mrativar = reform(mrativar1, npix, nobj)

   ;--------
   ; Re-bin the spectra to the same spacing

   minlog1 = min(loglam, max=maxlog1)
   newloglam = wavevector(minlog1, maxlog1)
   nnewpix = n_elements(newloglam)

   newratio = fltarr(nnewpix, nobj)
   newivar = fltarr(nnewpix, nobj)

   for iobj=0L, nobj-1 do begin
      isort = sort(loglam[*,iobj])
      combine1fiber, loglam[isort,iobj], mratio[isort,iobj], $
       mrativar[isort,iobj], $
       newloglam=newloglam, newflux=newratio1, newivar=newivar1
      newratio[*,iobj] = newratio1
      newivar[*,iobj] = newivar1
   endfor

   ;--------
   ; Compute the straight weighted mean at each wavelength
   ; (Avoid divide-by- zeros.)

   if (ndim EQ 1) then begin
      meanratio = (newratio * newivar) / (newivar + (newivar EQ 0))
   endif else begin
      denom = total(newivar, 2)
      meanratio = total(newratio * newivar, 2) / (denom + (denom EQ 0))
   endelse

   qbadpix = meanratio LE 0
   ibadpix = where(qbadpix, nbadpix)
   if (nbadpix GT 0) then newivar[ibadpix,*] = 0

   ;--------
   ; Actually take this "mean" and turn it into something more like
   ; a median, to protect us against standard stars that have bad
   ; magnitudes from the imaging.

; Comment-out ???
;   igoodpix = where(qbadpix EQ 0)
;   if (ndim EQ 1) then medratio = newratio $
;    else medratio = djs_median(newratio, 2)
;   rescale = median( medratio[igoodpix] / meanratio[igoodpix] )
;   if (rescale LE 0) then begin
;      splog, 'Warning: RESCALE = ', rescale
;   endif else begin
;      meanratio = rescale * meanratio
;      splog, 'Rescale factor median/mean = ', rescale
;   endelse

   ;--------
   ; Now for each object, compute the polynomial fit of it relative to the mean

   npoly = 3 ; ???
   flatarr = fltarr(npix, nobj)
   pres = fltarr(npoly, nobj)
   for iobj=0L, nobj-1 do begin
      ii = where(newivar[*,iobj] GT 0, ct)
      if (ct GT npoly+1) then begin ; At least NPOLY+1 pixels for a fit...
         thisloglam = newloglam[ii]
         thisratio = newratio[ii,iobj] / meanratio[ii]
         thisivar = newivar[ii,iobj] * meanratio[ii]^2

         ; This fit requires no rejection, because this function falls
         ; within an iteration loop that rejects points.

         ; The following is a weighted fit...
;         pres1 = poly_fit(thisloglam-3.5d0, thisratio, npoly-1, $
;          measure_errors=1./sqrt(thisivar))

         ; The following would be an unweighted fit...
;         pres1 = poly_fit(thisloglam-3.5d0, thisratio, npoly-1)

         ; The following is an unweighted fit but with outlier-rejection...
         poly_iter, thisloglam-3.5d0, thisratio, npoly-1, 3., coeff=pres1

         flatarr[*,iobj] = poly(loglam[*,iobj]-3.5d0, pres1)
         pres[*,iobj] = reform(pres1, npoly)
       endif else begin
         flatarr[*,iobj] = 1
         pres[*,iobj] = 0
         pres[0,iobj] = 1
       endelse
   endfor

   if (ndim GT 1) then $
    pres = reform(pres, [npoly, dims[1:ndim-1]])
   return, reform(flatarr, dims)
end

;------------------------------------------------------------------------------
pro spflux_plotcalib, mratiologlam, mratioflux, mrativar, $
 fitloglam, fitflux, fitflux2, logrange=logrange, plottitle=plottitle

   xrange = 10.^logrange
   ii = where(fitloglam GE logrange[0] AND fitloglam LE logrange[1])
   yrange = [0.9 * min(fitflux[ii]), 1.1 * max(fitflux[ii])]
   if (size(mratioflux, /n_dimen) EQ 1) then nfinal = 1 $
    else nfinal = (size(mratioflux, /dimens))[2]

   djs_plot, xrange, yrange, /xstyle, /ystyle, /nodata, $
    xtitle='Wavelength [Ang]', ytitle='Counts/(10^{-17}erg/cm^2/s/Ang', $
    title=plottitle
   for k=0, nfinal-1 do begin
      jj = where(mratiologlam[*,0,k] GE logrange[0] $
       AND mratiologlam[*,0,k] LE logrange[1] $
       AND mrativar[*,0,k] GT 0, ct)
      if (ct GT 1) then $
       djs_oplot, 10.^mratiologlam[jj,0,k], mratioflux[jj,0,k], psym=3
   endfor
   djs_oplot, 10.^fitloglam[ii], fitflux[ii], color='green'
   if (total(fitflux2) GT 0) then $
    djs_oplot, 10.^fitloglam[ii], fitflux2[ii], color='red'

   return
end

;------------------------------------------------------------------------------

function makelabel, hdr

   camera = strtrim(sxpar(hdr, 'CAMERAS'),2)
   expos =  strtrim(string(sxpar(hdr, 'EXPOSURE')),2)
   flat  =  strmid(sxpar(hdr, 'FLATFILE'),7,8)
   arc   =  strmid(sxpar(hdr, 'ARCFILE'),7,8)

   label = string(camera, expos, flat, arc, $
    format='(a2,"-",i8.8,"-",a8,"-",a8)')

   return, label
end

;------------------------------------------------------------------------------
pro add_iraf_keywords, hdr, wavemin, binsz

   sxaddpar, hdr, 'WAT0_001', 'system=linear'
   sxaddpar, hdr, 'WAT1_001', $
    'wtype=linear label=Wavelength units=Angstroms'
   sxaddpar, hdr, 'CRVAL1', wavemin, $
    ' Central wavelength (log10) of first pixel'
   sxaddpar, hdr, 'CD1_1', binsz, ' Log10 dispersion per pixel'
   sxaddpar, hdr, 'CRPIX1', 1, ' Starting pixel (1-indexed)'
   sxaddpar, hdr, 'CTYPE1', 'LINEAR'
   sxaddpar, hdr, 'DC-FLAG', 1, ' Log-linear flag'

   return
end


;------------------------------------------------------------------------------

function mlfluxcal,fileb1,filer1,obsparam,slitmap,waveset

; Status starts off OK
status=0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Start defining parameters for the extraction and calibration
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

splog,'Full reduction set- assuming spSFrame inputs.'
splog,'Reading science + wavelength data from files:'
splog,fileb1
splog,filer1

adderr=0.03
minfracthresh=0.80
objname=strarr(2)
objname[0]=fileb1
objname[1]=filer1

; Populate some variables with header information
nfile=n_elements(objname)
plateid=lonarr(nfile)
mjd=lonarr(nfile)
camname=strarr(nfile); b1,r1,b2, or r2
expnum=lonarr(nfile)
spectroid=lonarr(nfile); Spectrograph id (1 or 2)
npixarr=lonarr(nfile); Number of spectral pixels

for ifile=0,nfile-1 do begin
  hdr=headfits(objname[ifile])
  plateid[ifile] = strtrim(sxpar(hdr, 'PLATEID'),2)
  mjd[ifile] = strtrim(sxpar(hdr, 'MJD'),2)
  camname[ifile] = strtrim(sxpar(hdr, 'CAMERAS'),2)
  spectroid[ifile] = strmid(camname[ifile],1,1)
  expnum[ifile] = sxpar(hdr, 'EXPOSURE')
  npixarr[ifile] = sxpar(hdr, 'NAXIS1')
endfor
maxmjd = max(mjd)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Figure out which objects are F stars
; (i.e., which are 2'' cal fibers on std stars)
; BOSS software would look through extension 5 of the spFrame file,
; which I think is a copy of the plPlugMapM file, for the
; keyword SPECTROPHOTO_STD.  I'll use slitmap instead
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

mlframe_read, objname[0], hdr=hdr

; Determine whether standards are STD2D1, STD2D2, or STD2D3
stdname=strcompress('STD2D'+string(obsparam.dposn),/remove_all)

iphoto=where(slitmap.ifuname eq stdname,nphoto)

; This has defined iphoto, which indexes where the std spectra
; are along the slit, and nphoto, which is the total number of
; such standards
if (nphoto EQ 0) then begin
  splog,strcompress('Error in flux calibration; no standards found!')
  status=-140L
  mlquitmanga3d,status
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Read the raw F-star spectra
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

npix = max(npixarr); Number of spectral pixels
nfiber = n_elements(slitmap) ; Number of fibers in one spectrograph
; WARNING- the above will be FALSE when we move to using two spectrographs,
; will have to come up with something more clever

; Define wavelength, flux, inverse variance, dispersion, and airmass
; grids for all standards, in all files, at all spectra elements
loglam = fltarr(npix, nfile, nphoto)
objflux = fltarr(npix, nfile, nphoto)
objivar = fltarr(npix, nfile, nphoto)
dispimg = fltarr(npix, nfile, nphoto)
; Airmass is for all fibers
airmass = fltarr(npix, nfile, nfiber)

for ifile=0L, nfile-1 do begin
  ; Read info for the standard fibers
  mlframe_read, objname[ifile], iphoto, wset=wset1, loglam=loglam1, $
    objflux=objflux1, objivar=objivar1, dispimg=dispimg1, $
    mask=mask1, hdr=hdr1, adderr=adderr

  ; Compute the airmass for every pixel of every object
  ; (every pixel is the same, of course)
  get_tai, hdr1, tai_beg, tai, tai_end
  ; Defines tai, which is time in seconds since Nov 17 1858
  ; Calls RA, Dec from the .slm file
  for j=0, nfiber-1 do $
    airmass[0:npixarr[ifile]-1,ifile,j] = tai2airmass(slitmap[j].ra,slitmap[j].dec, tai=tai)

  ; Make a map of the size of each pixel in delta-(log10-Angstroms).
  ; Re-normalize the flux to ADU/(dloglam).
  ; Re-normalize the dispersion from /(raw pixel) to /(new pixel).
  correct_dlam, objflux1, objivar1, wset1, dlam=dloglam
  correct_dlam, dispimg1, 0, wset1, dlam=dloglam, /inverse

  ; Mask pixels on bad fibers
  objivar1 = objivar1 * spflux_goodfiber(mask1)

  loglam[0:npixarr[ifile]-1,ifile,*] = loglam1
  ;it wont do having a tail of zeros in the wavelength so add some dummy values
  if (npix GT npixarr[ifile]) then begin
    dllam=loglam1[npixarr[ifile]-1,*]-loglam1[npixarr[ifile]-2,*]
    for j=0, nphoto-1 do $
      loglam[npixarr[ifile]:*,ifile,j] = loglam1[npixarr[ifile]-1,j]+dllam[0,j]*(1+findgen(npix-npixarr[ifile]))
  endif
  objflux[0:npixarr[ifile]-1,ifile,*] = objflux1
  ;hopefully the inverse variance of 0 of non-filled objects will indicate the uselessness
  ; of the extra
  ; DRL- not sure I'm happy with what skymask actually DOES...
  ; Ok- it just masks out where sky values were flagged as bad
  objivar[0:npixarr[ifile]-1,ifile,*] = skymask(objivar1, mask1, mask1)
  dispimg[0:npixarr[ifile]-1,ifile,*] = dispimg1
endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Keep track of which F stars are good
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
qfinal = bytarr(nphoto) + 1B

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; For each star, find the best-fit model.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; !p is system structure for plotting variables
!p.multi = [0,2,3]
modflux = 0 * objflux

; Define sfdebv, which is dust contribution
; Convert from ra/dec to galactic long/lat
euler, slitmap.ra, slitmap.dec, ll, bb, 1
; Interpolate dust values from dust maps
sfdebv=dust_getval(ll,bb,/interp)

; Define correction, which is the AB/SDSS magnitude correction
correction = [-0.042, 0.036, 0.015, 0.013, -0.002]

; Loop over all of the standards
for ip=0L, nphoto-1 do begin
  thisfiber = iphoto[ip] + 1 + nfiber * (spectroid[0] - 1)
  splog, prelog='Fiber '+string(thisfiber,format='(I4)')

  plottitle = 'PLATE=' + string(plateid[0], format='(i4.4)') $
   + ' MJD=' + string(maxmjd, format='(i5.5)') $
   + ' Spectro-Photo Star' $
   + ' Fiber ' + strtrim(thisfiber,2)

  ; Find the best-fit model -- evaluated for each exposure [NPIX,NEXP]
  thismodel = spflux_bestmodel(loglam[*,*,ip], objflux[*,*,ip], $
   objivar[*,*,ip], dispimg[*,*,ip], kindx=thisindx, plottitle=plottitle)

  ; Also evaluate this model over a big wavelength range [3006,10960] Ang.
  tmploglam = 3.4780d0 + lindgen(5620) * 1.d-4
  tmpdispimg = 0 * tmploglam + 1.0 ; arbitrarily select this resolution
  tmpflux = spflux_read_kurucz(tmploglam, tmpdispimg, $
   iselect=thisindx.imodel)

  ; The returned models are redshifted, but not fluxed or
  ; reddened.  Do that now...  we compare data vs. model reddened.
  extcurve1 = ext_odonnell(10.^loglam[*,*,ip], 3.1)
  thismodel = thismodel $
   * 10.^(-extcurve1 * 3.1 * sfdebv[iphoto[ip]] / 2.5)
  extcurve2 = ext_odonnell(10.^tmploglam, 3.1)
  tmpflux = tmpflux $
   * 10.^(-extcurve2 * 3.1 * sfdebv[iphoto[ip]] / 2.5)

  ; Now integrate the apparent magnitude for this spectrum,
  ; The units of FTHRU are such that m = -2.5*alog10(FTHRU) + (48.6-2.5*17)
  ; Note that these computed magnitudes, THISMAG, should be equivalent
  ; to THISINDX.MAG in the case of no reddening.
  wavevec = 10.d0^tmploglam
  flambda2fnu = wavevec^2 / 2.99792e18
  fthru = filter_thru(tmpflux * flambda2fnu, waveimg=wavevec, /toair)
  thismag = -2.5 * alog10(fthru) - (48.6-2.5*17)

  ; Compute SCALEFAC = (plugmap flux) / (uncalibrated flux)
  ; We aren't using the plugmap, so we reconstruct it instead from
  ; calibflux = 10.^((22.5 - slitmap.mag-correction) /2.5)
  ; where correction is an AB conversion term
  ; correction = [-0.042, 0.036, 0.015, 0.013, -0.002]
  calibflux = 10.^((22.5 - slitmap[iphoto[ip]].mag-correction) /2.5)
  scalefac=calibflux[2]/ 10.^((22.5-thismag[2])/2.5)
  ; Reject this star if we don't know its flux.
  if (calibflux[2] LE 0) then begin
    splog, 'Warning: Rejecting std star in fiber = ', $
      iphoto[ip] + 1 + nfiber * (spectroid[0] - 1), $
      ' with unknown calibObj flux'
    qfinal[ip] = 0
  endif
  thismodel = thismodel * scalefac

  modflux[*,*,ip] = thismodel
  if (ip EQ 0) then kindx = replicate( create_struct( $
   'PLATE', 0L, $
   'MJD', 0L, $
   'FIBERID', 0L, $
   'QGOOD', 0, $
   thisindx, $
   'MODELFLUX', fltarr(npix)), $
   nphoto)
  copy_struct_inx, thisindx, kindx, index_to=ip
  kindx[ip].plate = plateid[0]
  kindx[ip].mjd = maxmjd
  kindx[ip].fiberid = thisfiber
  splog, prelog=''
endfor
!p.multi = 0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Reject any stars where more than 20% of the pixels marked are bad
; in any observation.
; DRL- this will FAIL as is for all MaNGA observations, about 39%
; will always be bad just from contribution from the crap on the ends
; of the spectrum...  I'm not sure why it works for BOSS, but kludge
; by resetting minfracthresh to 0.5 instead of 0.8
; Work out details later...
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
minfracthresh=0.50; DRL kludge

fracgood = fltarr(nphoto)
for ip=0L, nphoto-1 do begin
  for i=0L, nfile-1 do begin
    markasbad = (qfinal[ip]) AND (mean(objivar[*,i,ip] GT 0) LT minfracthresh)
    if (markasbad) then begin
      splog, 'Warning: Rejecting std star in fiber = ', $
        iphoto[ip] + 1 + nfiber * (spectroid[0] - 1), $
        ' with too many IVAR=0 pixels'
      qfinal[ip] = 0B
    endif
  endfor
endfor
ifinal = where(qfinal,nfinal) ; This is the list of the good stars
if (nfinal EQ 0) then begin
  splog, 'ABORT: No good fluxing stars!'
  status=-141L
  mlquitmanga3d,status
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Reject any stars with a bad chi^2/DOF either
; in the full spectrum or in just the absorp. line regions.
; Do not reject more than half the stars.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

chi2limit = 2.0 ; ???
chi2list = (kindx.chi2 / (kindx.dof>1)) $
 > (kindx.linechi2 / (kindx.linedof>1))
chi2list = chi2list + 100 * (kindx.linedof LT 10) ; Bad if < 10 pixels
while (max(chi2list) GT chi2limit AND total(qfinal) GT nphoto/2.) do begin
  chi2max = max(chi2list, iworst)
  splog, 'Rejecting std star in fiber = ', $
   iphoto[iworst] + 1 + nfiber * (spectroid[0] - 1), $
   ' with chi2=', chi2max
  chi2list[iworst] = 0
  qfinal[iworst] = 0B
endwhile

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Reject any stars with a very low S/N in the absorp. line regions.
; Do not reject more than half the stars.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

snlimit = 2.0 ; ??? DRL- 2 looks REALLY ugly....
; Do not reject any stars that are already rejected above...
snlist = kindx.linesn_median + (snlimit+1) * (qfinal EQ 0)
while (min(snlist) LT snlimit AND total(qfinal) GT nphoto/2.) do begin
  snmin = min(snlist, iworst)
  splog, 'Rejecting std star in fiber = ', $
   iphoto[iworst] + 1 + nfiber * (spectroid[0] - 1), $
   ' with median line S/N=', snmin
  snlist[iworst] = snlimit + 1
  qfinal[iworst] = 0B
endwhile

if (total(qfinal) EQ 0) then begin
  splog, 'ABORT: No good spectro-photo stars!'
  status=-142L
  mlquitmanga3d,status
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Loop over each exposure, and compute the PCA fit to MRATIO
; using outlier-rejection.
; Iterate, rejecting entire stars if they are terrible fits.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

iblue = where(strmatch(camname,'b*'), nblue)
ired = where(strmatch(camname,'r*'), nred)

qdone = 0L
iiter = 0L
while (qdone EQ 0) do begin
  iiter = iiter + 1
  splog, 'Iteration #', iiter

  ifinal = where(qfinal,nfinal) ; This is the list of the good stars

  ; The MRATIO vectors are the "raw" flux-calib vectors for each expos+CCD
  mmask = modflux GT 0
  mratio = objflux / (modflux*mmask + (1-mmask))
  mrativar = objivar * modflux^2
  flatarr = 0 * mratio

  ; Ignore regions around the stellar features
  mrativar = mrativar * (1 - spflux_masklines(loglam, /stellar))

  ; For each camera (blue or red), divide-out a low-order polynomial from
  ; MRATIO each star to get them all to the same mean flux levels.
  ; This takes out gross throughput differences between exposures.
  ; Also, it will remove the ~5% large-scale spectrophotometry errors
  ; between invidual stars, both from spectrograph throughput variations
  ; and from slight mis-typing of the stars.

  if (nblue GT 0) then $
    flatarr[*,iblue,ifinal] = spflux_mratio_flatten(loglam[*,iblue,ifinal], $
     mratio[*,iblue,ifinal], mrativar[*,iblue,ifinal], pres=pres_b)

  if (nred GT 0) then $
    flatarr[*,ired,ifinal] = spflux_mratio_flatten(loglam[*,ired,ifinal], $
     mratio[*,ired,ifinal], mrativar[*,ired,ifinal], pres=pres_r)

  mratio[*,*,ifinal] = mratio[*,*,ifinal] / flatarr[*,*,ifinal]
  mrativar[*,*,ifinal] = mrativar[*,*,ifinal] * flatarr[*,*,ifinal]^2

  ; Do the B-spline fits for the blue CCDs.
  if (nblue GT 0) then begin
    everyn = nblue * nfinal * 10
    sset_b = spflux_bspline(loglam[*,iblue,ifinal], $
     mratio[*,iblue,ifinal], mrativar[*,iblue,ifinal], $
     everyn=everyn, outmask=mask_b)
  endif

  ; Do the B-spline fits for the red CCDs.
  ; Fit a 2-dimension B-spline using the airmass as the 2nd dimension,
  ; but only if the airmass spans at least 0.10 and there are at
  ; least 3 good stars.

  if (nred GT 0) then begin
    everyn = nred * nfinal * 1.5
    if (max(airmass) - min(airmass) GT 0.10 AND nfinal GE 3) then begin ; ???
      ; Get an airmass value for every *pixel* being fit
      thisair = airmass[*,ired,iphoto[ifinal]]
    endif else begin
      thisair = 0
    endelse
    sset_r = spflux_bspline(loglam[*,ired,ifinal], $
     mratio[*,ired,ifinal], mrativar[*,ired,ifinal], $
     everyn=everyn, outmask=mask_r, airmass=thisair)
  endif

  ; Find which star has the most pixels rejected, and reject
  ; that star if it's bad enough

  fracgood = fltarr(nfinal)
  for k=0L, nfinal-1 do begin
    if (nblue GT 0) then fracgood[k] += 0.5 * mean(mask_b[*,*,k])
    if (nred GT 0) then fracgood[k] += 0.5 * mean(mask_r[*,*,k])
  endfor
  minfrac = min(fracgood, iworst)
  if (minfrac LT minfracthresh) then begin
    if (nfinal LE nphoto/2.) then begin
      splog, 'WARNING: Already rejected ', nphoto-nfinal, ' of ', $
       nphoto, ' std stars'
      qdone = 1B
    endif else begin
      splog, 'Rejecting std star in fiber = ', $
       iphoto[ifinal[iworst]] + 1 + nfiber * (spectroid[0] - 1), $
       ' with fracgood=', minfrac
      qfinal[ifinal[iworst]] = 0B
    endelse
   endif else begin
     qdone = 1B ; No other stars to reject
  endelse
endwhile

kindx[ifinal].qgood = 1
splog, 'Rejected ', nphoto-nfinal, ' of ', nphoto, ' std stars'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Plot fluxing vectors and their polynomial offsets for individual stars
; in individual exposures.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

mratfit = 0 * mratio
if (nblue GT 0) then $
  mratfit[*,iblue,ifinal] = bspline_valu(loglam[*,iblue,ifinal], sset_b)
if (nred GT 0) then $
  if (tag_exist(sset_r,'NPOLY')) then $
    mratfit[*,ired,ifinal] = bspline_valu(loglam[*,ired,ifinal], sset_r, $
      x2=airmass[*,ired,iphoto[ifinal]]) $
  else $
     mratfit[*,ired,ifinal] = bspline_valu(loglam[*,ired,ifinal], sset_r)

!p.multi = [0,1,2]
explist = expnum[uniq(expnum, sort(expnum))]
colorvec = ['default','red','green','blue','cyan','magenta','grey']
xrange = 10^minmax(loglam[*,*,ifinal])
ii = where(mrativar[*,*,ifinal] GT 0, ct)
if (ct GT 1) then yrange = minmax((mratfit[*,*,ifinal])[ii]) $
else yrange = minmax(mratfit[*,*,ifinal])
plottitle = 'PLATE=' + string(plateid[0], format='(i4.4)') $
  + ' MJD=' + string(maxmjd, format='(i5.5)')
for iexp=0, n_elements(explist)-1 do begin
  djs_plot, [0], [0], xrange=xrange, yrange=yrange, /xstyle, /ystyle, $
    /nodata, xtitle='Wavelength [Ang]', ytitle='Flux-calib', $
    title=plottitle+' Exp #' + string(explist[iexp], format='(i8.8)')
  kk = where(expnum EQ explist[iexp], kct) ; blue+red files for this exp
  for j=0, nfinal-1 do begin
    thiscolor = colorvec[j MOD n_elements(colorvec)]
    for k=0, kct-1 do begin
      djs_oplot, 10^loglam[*,kk[k],ifinal[j]], $
;      mratio[*,kk[k],ifinal[j]] * flatarr[*,kk[k],ifinal[j]], $
       mratio[*,kk[k],ifinal[j]], $
       psym=3, color=thiscolor
;      djs_oplot, 10^loglam[*,kk[k],ifinal[j]], $
;;     mratfit[*,kk[k],ifinal[j]] * flatarr[*,kk[k],ifinal[j]], $
;      mratfit[*,kk[k],ifinal[j]], $
;      color=thiscolor
    endfor
    djs_xyouts, 0.9*xrange[0]+0.1*xrange[1], $
      yrange[0] + (j+1)*(yrange[1]-yrange[0])/(nfinal+1), $
      'Fiber '+strtrim(iphoto[ifinal[j]]+(spectroid[0]-1)*nfiber,2), $
      color=thiscolor
  endfor
endfor
!p.multi = 0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Construct the final (B-splined) flux-calibration vectors
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
for ifile=0, nfile-1 do begin
  splog, prelog=fileandpath(objname[ifile])

  ; Is this a blue CCD?
  ii = where(ifile EQ iblue, ct)
  if (ct EQ 1) then begin
    thisloglam = loglam[*,ifile,ifinal]
    thisset = sset_b
    thisflatarr = flatarr[*,iblue[ii],ifinal]
    thispres = pres_b[*,ii,ifinal]
  endif

  ; Is this a red CCD?
  ii = where(ifile EQ ired, ct)
  if (ct EQ 1) then begin
    thisloglam = loglam[*,ifile,ifinal]
    thisset = sset_r
    thisflatarr = flatarr[*,ired[ii],ifinal]
    thispres = pres_r[*,ii,ifinal]
  endif

  thismratio = mratio[*,ifile,ifinal]
  thismrativar = mrativar[*,ifile,ifinal]
  if (tag_exist(thisset,'NPOLY')) then $
    x2 = airmass[*,ifile,iphoto[ifinal]] $
  else $
    x2 = 0

  ; Evaluate the B-spline for the stars at their measured wavelengths
  ; in this exposure, then modulated by the mean FLATARR
  ; for the stars in this exposure.
  ; We re-fit the B-spline to exactly recover what we had before,
  ; just modulated by the lower-order polynomial FLATARR.

  logmin = min(thisloglam[where(thismrativar GT 0)], max=logmax)
  tmploglam = wavevector(logmin, logmax)
  flatarr_mean = 0 * tmploglam
  for i=0L, nfinal-1 do $
    flatarr_mean = flatarr_mean $
    + poly(tmploglam-3.5d0, thispres[*,0,i]) / nfinal

  ; Rather than re-generating the B-spline, I'll simply cheat and
  ; multiply the B-spline coefficients at their wavelengths.
  ; This is a bit of a hack, since there are NORD more values
  ; for FULLBKPT than there are for COEFF, so there's not exactly
  ; a 1-to-1 mapping between the two.
  tmpmult = interpol(flatarr_mean, tmploglam, thisset.fullbkpt)
  tmpmult = tmpmult[1:n_elements(thisset.fullbkpt)-thisset.nord]
  if (keyword_set(x2)) then begin
    for ipoly=0, thisset.npoly-1 do $
      thisset.coeff[ipoly,*] = thisset.coeff[ipoly,*] * tmpmult
    for ipoly=0, thisset.npoly-1 do $
      thisset.icoeff[ipoly,*] = thisset.icoeff[ipoly,*] * tmpmult
  endif else begin
    thisset.coeff = thisset.coeff * tmpmult
    thisset.icoeff = thisset.icoeff * tmpmult
  endelse

  if (keyword_set(x2)) then begin
    x2_min = min(x2, max=x2_max)
     splog, 'Exposure ', objname[ifile], $
       ' spans airmass range ', x2_min, x2_max
    tmpflux1 = bspline_valu(tmploglam, thisset, x2=x2_min+0*tmploglam)
    tmpflux2 = bspline_valu(tmploglam, thisset, x2=x2_max+0*tmploglam)
  endif else begin
    tmpflux1 = bspline_valu(tmploglam, thisset)
    tmpflux2 = 0
  endelse

  ; Make plots of the spectro-photometry data for this exposure only,
  ; overplotting the global fit to all exposures in red.

  ; The following info is just used for the plot title
  plottitle = 'PLATE=' + string(plateid[ifile], format='(i4.4)') $
    + ' MJD=' + string(mjd[ifile], format='(i5.5)') $
    + ' Spectro-Photo Calib for ' + camname[ifile] + '-' $
    + string(expnum[ifile], format='(i8.8)')

  !p.multi = [0,1,2]
  logrange = logmax - logmin
  spflux_plotcalib, $
    thisloglam, thismratio, thismrativar, $
    tmploglam, tmpflux1/flatarr_mean, tmpflux2/flatarr_mean, $
    logrange=(logmin+[0,1]*logrange/2.), plottitle=plottitle
  spflux_plotcalib, $
    thisloglam, thismratio, thismrativar, $
    tmploglam, tmpflux1/flatarr_mean, tmpflux2/flatarr_mean, $
    logrange=(logmin+[1,2]*logrange/2.)
  !p.multi = 0

  ; Create header cards describing the fit range
  hdr = ['']
  sxaddpar, hdr, 'WAVEMIN', 10.^logmin
  sxaddpar, hdr, 'WAVEMAX', 10.^logmax

  ; Generate the pixel map of the flux-calibration for this exposure+CCD
  mlframe_read, objname[ifile], loglam=loglam1
  if (tag_exist(thisset,'NPOLY')) then x2 = airmass[*,ifile,*] $
  else x2 = 0
  calibimg = float( bspline_valu(loglam1, thisset, x2=x2) )

  ; Set to zero any pixels outside the known flux-calibration region
  qbad = loglam1 LT logmin OR loglam1 GT logmax
  ibad = where(qbad, nbad)
  if (nbad GT 0) then calibimg[ibad] = 0

  minval = min(calibimg[where(qbad EQ 0)], max=maxval)
  if (maxval/minval GT 20. OR minval LT 0) then $
    splog, 'WARNING: Min/max fluxcalib = ', minval, maxval $
  else $
    splog, 'Min/max fluxcalib = ', minval, maxval

  ; Write the output file

  ; Put the Kurucz models for this exposure in the output structure
  kindx.modelflux = reform(modflux[*,ifile,*], npix, nphoto)

  calibfile = djs_filepath(string(camname[ifile], expnum[ifile], $
    format='("spFluxcalib-", a2, "-", i8.8, ".fits")'), $
    root_dir=datapath)
  mwrfits, calibimg, calibfile, hdr, /create
  mwrfits, thisset, calibfile
  mwrfits, kindx, calibfile
  spawn, ['gzip','-f',calibfile], /noshell
  ;    splog, prelog=''
endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; BOSS routine then called spcoadd_v5 to actually *apply* the flux
; calibration and combine b+r spectra.  We'll do that here.
; Technically we shouldn't need to read stuff in again... but it's
; easier to adapt from BOSS this way.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;fileb1,filer1,obsparam,slitmap,wave
binsz = 1.0d-4
zeropoint = 3.5D
window = 100

; Sort filenames such that this list contains first the blue then the red
spframes=objname
nfiles = n_elements(spframes)
filenames = spframes[sort(spframes)]

camnames=camname
ncam = N_elements(camnames)
nexpvec = lonarr(ncam)
exptimevec = fltarr(ncam)

; Loop through each 2D output and read in the data

; Start by determining the size of all the files
npixarr = lonarr(nfiles)
for ifile=0, nfiles-1 do begin
  mlframe_read, filenames[ifile], hdr=objhdr
  npixarr[ifile] = sxpar(objhdr,'NAXIS1')
endfor
npixmax = max(npixarr)
nobj = sxpar(objhdr,'NAXIS2') ; Number of fibers per spectrograph
; DRL unhappy with this...  Why?

for ifile=0, nfiles-1 do begin
  ; Read in all data from this input file.
  ; Reading the plug-map structure will fail if its structure is
  ; different between different files.

  splog, 'Reading file #', ifile, ': ', filenames[ifile]
  mlframe_read, filenames[ifile], objflux=tempflux, objivar=tempivar, $
    mask=temppixmask, wset=tempwset, dispset=tempdispset, $
    skyflux=tempsky, ximg=tempximg, superflat=tempsuperflat, $
    hdr=hdr, adderr=adderr

  if (ifile EQ 0) then $
    hdrarr = ptr_new(hdr) $
  else $
    hdrarr = [hdrarr, ptr_new(hdr)]

  ; Read header info
  thismjd = sxpar(hdr, 'MJD')
  if (NOT keyword_set(mjdlist)) then mjdlist = thismjd $
  else mjdlist = [mjdlist, thismjd]
  cameras = strtrim(sxpar(hdr, 'CAMERAS'),2)
  expstr = string(sxpar(hdr, 'EXPOSURE'), format='(i8.8)')

  ; Solve for wavelength and lambda-dispersion at each pixel in the image
  ;need to reset junk since the array lengths change
  junk=0
  traceset2xy, tempwset, junk, tempwave
  traceset2xy, tempdispset, junk, tempdispersion

  ; Here is the correct conversion from pixels to log-lambda dispersion.
  ; We are converting from the dispersion in units of spFrame pixel sizes
  ; to the dispersion in units of the new rebinned pixel size, which is
  ; BINSZ in log-lambda units.
      
  ; this probably should be fixed elsewhere but limit where the fit range
  tempxmax=tempwset.xmax
  tempwset.xmax=(size(tempdispersion,/dimens))[0]-1
  correct_dlam, tempdispersion, 0, tempwset, dlam=binsz, /inverse
  tempwset.xmax=tempxmax

  dims = size(tempflux, /dimens)
  npix = dims[0]
  nfib = dims[1]

  ; Make a map of the size of each pixel in delta-(log10-Angstroms),
  ; and re-normalize the flux to electrons/(dloglam)
  correct_dlam, tempflux, tempivar, tempwset, dlam=binsz
  correct_dlam, tempsky, 0, tempwset, dlam=binsz

  ; Determine if this is a blue or red spectrum
  icam = (where(cameras EQ camnames))[0]
  if (icam EQ -1) then $
    message, 'Unknown camera ' + cameras
  nexpvec[icam] = nexpvec[icam] + 1
  exptimevec[icam] = exptimevec[icam] + sxpar(hdr, 'EXPTIME')

  ; Apply spectro-photometric calibration
  expnum = sxpar(hdr, 'EXPOSURE')
  calibfile = djs_filepath(string(camnames[icam], expnum, $
    format='("spFluxcalib-", a2, "-", i8.8, ".fits")'), $
    root_dir=datapath)
  calibfile = (findfile(calibfile+'*'))[0]
splog,'DRL- calibtest'

  if (keyword_set(calibfile)) then begin
    calibfac = mrdfits(calibfile, 0, calibhdr, /silent)
  endif else begin
    splog, 'WARNING: Reading default flux-calib vectors for camera=' $
       + camnames[icam]
    calibfac = fcalib_default(camnames[icam], tempwave, exptimevec[icam])
  endelse
  minval = 0.05 * mean(calibfac)
  divideflat, tempflux, invvar=tempivar, calibfac, minval=minval
  divideflat, tempsky, calibfac, minval=minval
  temppixmask = temppixmask $
    OR ((calibfac LE minval OR keyword_set(calibfile) EQ 0) $
    * pixelmask_bits('BADFLUXFACTOR'))

  ; This is ordinarily where BOSS would do a flux correction factor to 
  ; account for DAR, among other things.  I am NOT doing this.
  ; It's possible that I need to do part of it though...

  ; Apodize the errors
  ; Do this only for the dichroic overlap region, which are the first
  ; rows in both the blue and red CCD's.

  if (keyword_set(window)) then begin
    swin = window < npix
    indx = lindgen(swin)
    tempivar[indx,*] = $
    tempivar[indx,*] * (indx # replicate(1,nfib)) / swin
  endif

  ; Concatenate data from all images
  if (ifile EQ 0) then begin
    ; Construct the image arrays
    flux = make_array(npixmax,nobj*nfiles,type=size(tempflux,/type))
    fluxivar = make_array(npixmax,nobj*nfiles,type=size(tempivar,/type))
    wave = make_array(npixmax,nobj*nfiles,type=size(tempwave,/type))
    dispersion = make_array(npixmax,nobj*nfiles,type=size(tempdisp,/type))
    pixelmask = make_array(npixmax,nobj*nfiles,type=size(temppixmask,/type))
    skyflux = make_array(npixmax,nobj*nfiles,type=size(tempsky,/type))
    ximg = make_array(npixmax,nobj*nfiles,type=size(tempximg,/type))
    superflat = make_array(npixmax,nobj*nfiles,type=size(tempsuperflat,/type))

    ; Append as vectors...
    camerasvec = cameras
    label = makelabel(hdr)
    filenum = lonarr(nfib) + ifile
  endif else begin
    ; Append as vectors...
    camerasvec = [camerasvec, cameras]
    label = [label, makelabel(hdr)]
    filenum = [filenum, lonarr(nfib) + ifile]
  endelse

  flux[0:npixarr[ifile]-1,nobj*ifile:nobj*(ifile+1)-1] = tempflux
  fluxivar[0:npixarr[ifile]-1,nobj*ifile:nobj*(ifile+1)-1] = tempivar
  wave[0:npixarr[ifile]-1,nobj*ifile:nobj*(ifile+1)-1] = tempwave
  ; Pad the wavelengths with reasonable values
  if (npixarr[ifile] LT npixmax) then begin
    dwave = tempwave[npixarr[ifile]-1,*] - tempwave[npixarr[ifile]-2,*]
    addwave = tempwave[npixarr[ifile]-1,*] $
      ## (1+lonarr(npixmax-npixarr[ifile])) $
      + dwave ## (1+lindgen(npixmax-npixarr[ifile]))
    wave[npixarr[ifile]:npixmax-1,nobj*ifile:nobj*(ifile+1)-1] = addwave
  endif
  dispersion[0:npixarr[ifile]-1,nobj*ifile:nobj*(ifile+1)-1] = tempdispersion
  pixelmask[0:npixarr[ifile]-1,nobj*ifile:nobj*(ifile+1)-1] = temppixmask
  skyflux[0:npixarr[ifile]-1,nobj*ifile:nobj*(ifile+1)-1] = tempsky
  ximg[0:npixarr[ifile]-1,nobj*ifile:nobj*(ifile+1)-1] = tempximg
  superflat[0:npixarr[ifile]-1,nobj*ifile:nobj*(ifile+1)-1] = tempsuperflat
endfor

tempflux = 0
tempivar = 0
tempwave = 0
tempdispersion = 0
temppixmask = 0
tempsky = 0

; Check how many exposures we have in each of the (4) cameras

for icam=0, ncam-1 do begin
  junk = where(camerasvec EQ camnames[icam], nmatch)
  splog, 'Files for camera ' + camnames[icam] + ':', nmatch
  if (icam EQ 0) then nminfile = nmatch $
  else nminfile = nminfile < nmatch
endfor
; ??? Should make this routine robust to fewer files!!!
if (nminfile LT 1) then begin
  splog, 'ABORT: At least 1 file needed for each camera'
  status=-143L
  mlquitmanga3d,status
endif

; Construct output data structures, including the wavelength scale
totalpix = (size(flux, /dimens))[0]
nonzero = where(fluxivar GT 0.0)
minfullwave = min(wave[nonzero])
maxfullwave = max(wave[nonzero])

; Get max and min wavelength from good pixels
if (NOT keyword_set(wavemin)) then begin
  spotmin = long((minfullwave - zeropoint)/binsz) + 1L
  spotmax = long((maxfullwave - zeropoint)/binsz)
  wavemin = spotmin * binsz + zeropoint
  wavemax = spotmax * binsz + zeropoint
endif else begin
  spotmin = 0L
  if (NOT keyword_set(wavemax)) then begin
    spotmax = long((maxfullwave - wavemin)/binsz)
    wavemax = spotmax * binsz + wavemin
  endif else spotmax = long((wavemax - wavemin)/binsz)
endelse

; This is where output wavelengths are set
; kludge to be FIXED to BOSSCAL wavelengths
; so that everything is on same grid without
; requiring multiple resamplings
;nfinalpix = spotmax - spotmin + 1L
;finalwave = dindgen(nfinalpix) * binsz + wavemin
nfinalpix = (size(waveset))[1]
finalwave = dindgen(nfinalpix) * binsz + alog10(waveset[0])

; This seems to be BOSS scaling up from one detector to two??
;nfiber = 2 * nfib

finalflux = fltarr(nfinalpix, nfiber)
finalivar = fltarr(nfinalpix, nfiber)
finalandmask = lonarr(nfinalpix, nfiber)
finalormask = lonarr(nfinalpix, nfiber)
finaldispersion = fltarr(nfinalpix, nfiber)
finalsky = fltarr(nfinalpix, nfiber)

; Note that currently this is combining ALL fibers

; Combine each fiber, one at a time
for ifiber=0, nfiber-1 do begin
  ; Find the first occurance of fiber number IFIBER+1
  indx = (where(slitmap.fiberid EQ ifiber+1))[0]

  if (indx NE -1) then begin
    splog, 'Fiber', ifiber+1, ' ', slitmap[indx].ifuname, $
      slitmap[indx].mag, format = '(a, i5.4, a, a, f6.2, 5f6.2)'

    ; Identify all objects with the same XFOCAL,YFOCAL plate position, and
    ; combine all these objects into a single spectrum.
    ; If all pluggings are identical, then this will always be
    ; the same fiber ID.
    ; Also, insist that the object type is not 'NA', which would
    ; occur for unplugged fibers. <--- Disable this for BOSS ???

    ; DRL- comment this out, do it for ALL fibers.
    ; Ah- this is how we combine b+r, is by finding both locations
    ; where a fiber is in both b and r rows.
    ; WEIRD- why does BOSS do it this v awkward way???
    indx=[where(slitmap.fiberid eq ifiber+1),where(slitmap.fiberid eq ifiber+1)+nfiber]
    if (ifiber eq 0) then splog,indx
;    indx = where(abs(slitmap.xfocal - plugmap[indx].xfocal) LT 0.0001 $
 ;     AND abs(plugmap.yfocal - plugmap[indx].yfocal) LT 0.0001)
  endif

  if (indx[0] NE -1) then begin
    temppixmask = pixelmask[*,indx]
    combine1fiber, wave[*,indx], flux[*,indx], fluxivar[*,indx], $
      finalmask=temppixmask, indisp=dispersion[*,indx], $
      skyflux=skyflux[*,indx], $
      newloglam=finalwave, newflux=bestflux, newivar=bestivar, $
      andmask=bestandmask, ormask=bestormask, newdisp=bestdispersion, $
      newsky=bestsky, $
      nord=nord, binsz=binsz, bkptbin=bkptbin, maxsep=maxsep, $
      maxiter=50, upper=3.0, lower=3.0, maxrej=1

    finalflux[*,ifiber] = bestflux
    finalivar[*,ifiber] = bestivar
    finalandmask[*,ifiber] = bestandmask
    finalormask[*,ifiber] = bestormask
    finaldispersion[*,ifiber] = bestdispersion
    finalsky[*,ifiber] = bestsky

    ; DRL- a really stupid kludge to take out wacky values
    ; in the far blue.  Should do something better.
    for j=0,100 do begin
      if abs(finalflux[j,ifiber]) gt 100. then finalflux[j,ifiber]=0.
    endfor

    ; The following adds the COMBINEREJ bit to the input pixel masks
    pixelmask[*,indx] = temppixmask
  endif else begin
    splog, 'Fiber', ifiber+1, ' NO DATA'
    finalandmask[*,ifiber] = pixelmask_bits('NODATA')
    finalormask[*,ifiber] = pixelmask_bits('NODATA')
  endelse
endfor


   ; Modify the 1st file's header to use for the combined plate header.

   bighdr = *hdrarr[0]

; BOSS does a 'flux distortion' correction here, do I need to???


  ;Indiv calib frames?  DRL should remove


;   for ifile=0, nfiles-1 do begin
;      thisfile = fileandpath(filenames[ifile], path=thispath)
;      thisfile = djs_filepath(repstr(thisfile,'spSFrame','spTFrame'), $
;       root_dir=thispath)
;      splog, 'Writing file #', ifile, ': ', thisfile, prename=filenames[ifile]
;      indx = where(filenum EQ ifile, nthis)
;
;      hdr = *hdrarr[ifile]
;      sxaddpar, hdr, 'BUNIT', '1E-17 erg/cm^2/s/Ang'
;
;      ; Apply the flux-distortion image to each individual frame, by
;      ; interpolating off the full wavelength-scale distortion image
;      ; onto the wavelength mapping of each individual exposure+CCD.
;
;      for i=0L, nthis-1 do begin
;         thisflux1 = flux[*,indx[i]]
;         thisivar1 = fluxivar[*,indx[i]]
;         thissky1 = skyflux[*,indx[i]]
;         ;;j = slitmap[indx[i]].fiberid - 1
;         ;;thisicorr = interpol(invcorrimg[*,j], finalwave, wave[*,indx[i]])
;         ;;divideflat, thisflux1, invvar=thisivar1, thisicorr, minval=minicorrval
;         flux[*,indx[i]] = thisflux1
;         fluxivar[*,indx[i]] = thisivar1
;
;         ;divideflat, thissky1, thisicorr, minval=minicorrval
;         skyflux[*,indx[i]] = thissky1
;      endfor
;
;      mwrfits, flux[*,indx], thisfile, hdr, /create
;      mwrfits, fluxivar[*,indx], thisfile
;      mwrfits, pixelmask[*,indx], thisfile
;      mwrfits, wave[*,indx], thisfile
;      mwrfits, dispersion[*,indx], thisfile
;      mwrfits, slitmap[indx], thisfile
;      mwrfits, skyflux[*,indx], thisfile
;      mwrfits, ximg[*,indx], thisfile
;      mwrfits, superflat[*,indx], thisfile
;   endfor
;   splog, prename=''

   ;----------
   ; Clear memory

   wave = 0
   flux = 0
   fluxivar = 0
   temppixmask = 0
   dispersion = 0
   skyflux = 0
   superflat = 0



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Write the corrected spCFrame files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   ;---------------------------------------------------------------------------
   ; Create the output header
   ;---------------------------------------------------------------------------

   ;----------
   ; Print roll call of bad fibers and bad pixels.

   fiber_rollcall, finalandmask, finalwave

   ;----------
   ; Remove header cards that were specific to this first exposure
   ; (where we got the header).

   ncoeff = sxpar(bighdr, 'NWORDER')
   for i=2, ncoeff-1 do sxdelpar, bighdr, 'COEFF'+strtrim(string(i),2)

   sxdelpar, bighdr, ['SPA', 'IPA', 'IPARATE']
   sxdelpar, bighdr, 'EXPOSURE'
   sxdelpar, bighdr, 'REQTIME'
   sxdelpar, bighdr, 'QUALITY'
   sxdelpar, bighdr, 'FILENAME'
   sxdelpar, bighdr, 'SEQID'
   sxdelpar, bighdr, 'DARKTIME'
   sxdelpar, bighdr, 'CAMERAS'
   sxdelpar, bighdr, 'PLUGMAPO'
   for i=1, 4 do sxdelpar, bighdr, 'GAIN'+strtrim(string(i),2)
   for i=1, 4 do sxdelpar, bighdr, 'RDNOISE'+strtrim(string(i),2)
   sxdelpar, bighdr, ['CAMCOL', 'CAMROW']
   sxdelpar, bighdr, ['AMPLL', 'AMPLR', 'AMPUL', 'AMPUR']
   sxdelpar, bighdr, ['FFS', 'FF', 'NE', 'HGCD']
   sxdelpar, bighdr, ['SPEC1', 'SPEC2']
   sxdelpar, bighdr, 'NBLEAD'
   sxdelpar, bighdr, 'PIXFLAT'
   sxdelpar, bighdr, 'PIXBIAS'
   sxdelpar, bighdr, 'FLATFILE'
   sxdelpar, bighdr, 'ARCFILE'
   sxdelpar, bighdr, 'OBJFILE'
   sxdelpar, bighdr, 'FRAMESN2'
   sxdelpar, bighdr, 'DEREDSN2'

   ;----------
   ; Average together some of the fields from the individual headers.

   cardname = [ 'AZ', 'ALT', 'TAI', 'WTIME', 'AIRTEMP', 'DEWPOINT', $
    'DEWDEP', 'DUSTA', 'DUSTB', 'DUSTC', 'DUSTD', 'GUSTS', 'HUMIDITY', $
    'HUMIDOUT', 'PRESSURE', 'WINDD', 'WINDS', 'TEMP01', 'TEMP02', $
    'TEMP03', 'TEMP04', 'HELIO_RV', 'SEEING20', 'SEEING50', 'SEEING80', $
    'RMSOFF20', 'RMSOFF50', 'RMSOFF80', 'XCHI2', 'SKYCHI2', $
    'WSIGMA', 'XSIGMA' ]
   sxcombinepar, hdrarr, cardname, bighdr, func='average'

   sxcombinepar, hdrarr, 'TAI-BEG', bighdr, func='min'
   sxcombinepar, hdrarr, 'TAI-END', bighdr, func='max'

   sxcombinepar, hdrarr, 'XCHI2', bighdr, func='max', outcard='XCHI2MAX', $
    after='XCHI2'
   sxcombinepar, hdrarr, 'XCHI2', bighdr, func='min', outcard='XCHI2MIN', $
    after='XCHI2'

   sxcombinepar, hdrarr, 'SKYCHI2', bighdr, func='max', outcard='SCHI2MAX', $
    after='SKYCHI2'
   sxcombinepar, hdrarr, 'SKYCHI2', bighdr, func='min', outcard='SCHI2MIN', $
    after='SKYCHI2'

   sxcombinepar, hdrarr, 'WSIGMA', bighdr, func='max', outcard='WSIGMAX', $
    after='WSIGMA'
   sxcombinepar, hdrarr, 'WSIGMA', bighdr, func='min', outcard='WSIGMIN', $
    after='WSIGMA'

   sxcombinepar, hdrarr, 'XSIGMA', bighdr, func='max', outcard='XSIGMAX', $
    after='XSIGMA'
   sxcombinepar, hdrarr, 'XSIGMA', bighdr, func='min', outcard='XSIGMIN', $
    after='XSIGMA'

   ; Add the NGUIDE keywords for all headers of one flavor of CAMERAS
   ; (e.g., for all the 'b1' exposures if the first frame is 'b1'.)
   cardname = 'NGUIDE'
   sxcombinepar, hdrarr[0], cardname, bighdr, func='total'
   cameras0 = sxpar(*(hdrarr[0]), 'CAMERAS')
   for ihdr=1, n_elements(hdrarr)-1 do begin
      if (sxpar(*(hdrarr[ihdr]), 'CAMERAS') EQ cameras0) then $
       sxcombinepar, hdrarr[ihdr], cardname, bighdr, func='total'
   endfor

   ;----------
   ; Use the MJD passed as a keyword, which will typically be for the most
   ; observation, and be consistent with the output file names

   ; DRL for some reason this was vector??? Kludge.
   if (keyword_set(mjd)) then $
    sxaddpar, bighdr, 'MJD', obsparam.mjd

   ; Get the list of MJD's used for these reductions, then convert to a string
   mjdlist = mjdlist[uniq(mjdlist, sort(mjdlist))]
   mjdlist = strtrim(strcompress(string(mjdlist,format='(99a)')),2)
   sxaddpar, bighdr, 'MJDLIST', mjdlist, after='MJD'

   ;----------
   ; Add new header cards

   sxaddpar, bighdr, 'VERSCOMB', idlspec2d_version(), $
    ' Version of idlspec2d for combining multiple spectra', after='VERS2D'
   sxaddpar, bighdr, 'NEXP', nfiles, $
    ' Number of exposures in this file', before='EXPTIME'
   for ifile=0,nfiles-1 do $
    sxaddpar, bighdr, string('EXPID',ifile+1, format='(a5,i2.2)'), label[ifile], $
     ' ID string for exposure '+strtrim(ifile+1,2), before='EXPTIME'
   if (keyword_set(bestexpnum)) then $
    sxaddpar, bighdr, 'BESTEXP', bestexpnum, before='EXPID01'

   sxaddpar, bighdr, 'EXPTIME', min(exptimevec), $
    ' Minimum of exposure times for all cameras'
   for icam=0, ncam-1 do $
    sxaddpar, bighdr, 'NEXP_'+camnames[icam], nexpvec[icam], $
     ' '+camnames[icam]+' camera number of exposures', before='EXPTIME'
   for icam=0, ncam-1 do $
    sxaddpar, bighdr, 'EXPT_'+camnames[icam], exptimevec[icam], $
     ' '+camnames[icam]+' camera exposure time (seconds)', before='EXPTIME'
   sxaddpar, bighdr, 'SPCOADD', systime(), $
    ' SPCOADD finished', after='EXPTIME'

   sxaddpar, bighdr, 'NWORDER', 2, ' Linear-log10 coefficients'
   sxaddpar, bighdr, 'NWORDER', 2, ' Linear-log10 coefficients'
   sxaddpar, bighdr, 'WFITTYPE', 'LOG-LINEAR', ' Linear-log10 dispersion'
   sxaddpar, bighdr, 'COEFF0', wavemin, $
    ' Central wavelength (log10) of first pixel'
   sxaddpar, bighdr, 'COEFF1', binsz, ' Log10 dispersion per pixel'

   sxaddpar, bighdr, 'NAXIS1', n_elements(bestflux)
   sxaddpar, bighdr, 'NAXIS2', nfiber

   spawn, 'uname -n', uname
   sxaddpar, bighdr, 'UNAME', uname[0]

   ;----------
   ; Compute the fraction of bad pixels in total, and on each spectrograph.
   ; Bad pixels are any with SKYMASK(INVVAR)=0, excluding those where
   ; the NODATA bit is set in the pixel mask.

   ifib1 = where(slitmap.spectrographid EQ 1, nfib1)
   ifib2 = where(slitmap.spectrographid EQ 2, nfib2)
   qbadpix = skymask(finalivar, finalandmask, finalormask) EQ 0 $
    AND (finalandmask AND pixelmask_bits('NODATA')) EQ 0
   if (nfib1 GT 0) then $
    fbadpix1 = total(qbadpix[*,ifib1]) / (nfib1 * nfinalpix) $
   else $
    fbadpix1 = 0
   if (nfib2 GT 0) then $
    fbadpix2 = total(qbadpix[*,ifib2]) / (nfib2 * nfinalpix) $
   else $
    fbadpix2 = 0
   if (nfib1 GT 0 AND nfib2 GT 0) then $
    fbadpix = total(qbadpix[*,[ifib1,ifib2]]) / ((nfib1+nfib2) * nfinalpix) $
   else if (nfib1 GT 0) then $
    fbadpix = fbadpix1 $
   else if (nfib2 GT 0) then $
    fbadpix = fbadpix1 $
   else $
    fbadpix = 0

   sxaddpar, bighdr, 'FBADPIX', fbadpix, ' Fraction of bad pixels'
   sxaddpar, bighdr, 'FBADPIX1', fbadpix1, ' Fraction of bad pixels on spectro-1'
   sxaddpar, bighdr, 'FBADPIX2', fbadpix2, ' Fraction of bad pixels on spectro-2'

   ;----------
   ; Add keywords for IRAF-compatability

   add_iraf_keywords, bighdr, wavemin, binsz

   mkhdr, hdrfloat, finalivar, /image, /extend
   add_iraf_keywords, hdrfloat, wavemin, binsz

   mkhdr, hdrlong, finalandmask, /image, /extend
   add_iraf_keywords, hdrlong, wavemin, binsz

   ; First write the file with the flux distortion vectors
   ; DRL- didn't use this...
   ;mwrfits, corrimg, distortfitsfile, bighdr, /create

   fulloutname = str_replace(fileb1,'spSFrame-b1-','spCFrame-')

   ; HDU #0 is flux
   sxaddpar, bighdr, 'BUNIT', '1E-17 erg/cm^2/s/Ang'
   mwrfits, finalflux, fulloutname, bighdr, /create

   ; HDU #1 is inverse variance
   sxaddpar, hdrfloat, 'BUNIT', '1/(1E-17 erg/cm^2/s/Ang)^2'
   mwrfits, finalivar, fulloutname, hdrfloat

   ; HDU #2 is AND-pixelmask
   mwrfits, finalandmask, fulloutname, hdrlong

   ; HDU #3 is OR-pixelmask
   mwrfits, finalormask, fulloutname, hdrlong

   ; HDU #4 is dispersion map
   sxaddpar, hdrfloat, 'BUNIT', 'pixels'
   mwrfits, finaldispersion, fulloutname, hdrfloat

   ; HDU #5 is slitmap
   mwrfits, slitmap, fulloutname

   ; HDU #6 is the sky
   mwrfits, finalsky, fulloutname

   ; HDU #7 is the wavelength solution
   mwrfits, finalwave, fulloutname

return,status
end
