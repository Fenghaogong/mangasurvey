


pro ml_extract_object, outname, objhdr, image, invvar, plugsort, wset, $
 xarc, lambda, xtrace, fflat, fibermask, color=color, proftype=proftype, $
 widthset=widthset, dispset=dispset, skylinefile=skylinefile, $
 plottitle=plottitle, superflatset=superflatset, do_telluric=do_telluric, $
 bbspec=bbspec, splitsky=splitsky, ccdmask=ccdmask, VISUAL=VISUAL, SURVEY=survey, fiberparam=fiberparam
 
   on_error,0
  compile_opt idl2
  
   configuration=obj_new('configuration', sxpar(objhdr, 'MJD'))

   objname = strtrim(sxpar(objhdr,'OBJFILE'),2) 
   flavor  = strtrim(sxpar(objhdr,'FLAVOR'),2) 
   camera  = strtrim(sxpar(objhdr,'CAMERAS'),2) 
   nfiber = fiberparam.nfiber
   radius = fiberparam.radius
   totfiber = fiberparam.totfiber
   nbundle = (long(total(fiberparam.bundlegap NE 0)))[0]
   bundleid = fiberparam.bundleid
 
   ;------------------
   ; Identify very bright objects
   ; Do a boxcar extraction, and look for fibers where the median counts are 10000 ADU per row.

   fextract = extract_boxcar(image*(invvar GT 0), xtrace, nfiber=nfiber, radius=radius)
   scrunch = djs_median(fextract, 1) ; Find median counts/row in each fiber
   whopping = find_whopping(scrunch, 10000.0, whopct)
   scrunch_sort = sort(scrunch)
   i5 = n_elements(scrunch)/20
   i95 = i5 * 19

   splog, 'Whopping fibers: ', whopping
   splog, 'Median counts in all fibers = ', djs_median(scrunch)
   splog, 'Number of bright fibers = ', whopct
   
   if (whopct GT 20) then begin
      splog, 'WARNING: Disable whopping terms ' + objname
      whopping = -1
      whopct = 0
   endif   

   ;------------------------------------------------------------
   ;  Check for bad pixels within 3 pixels of trace
   badcheck = extract_boxcar((invvar LE 0), xtrace, radius=2.5, nfiber=nfiber)
   badplace = where(badcheck GT 0)

   nx = (size(fextract,/dim))[0] 
   ny = (size(fextract,/dim))[1] 
   pixelmask = lonarr(nx,ny)

   badcolumns = where(total(badcheck GT 0,1) GT 0.45 * nx)
   if (badplace[0] NE -1) then pixelmask[badplace] = pixelmask[badplace] OR pixelmask_bits('NEARBADPIXEL')
   if (badcolumns[0] NE -1) then fibermask[badcolumns] = fibermask[badcolumns] OR pixelmask_bits('MANYBADCOLUMNS')

   if (whopping[0] NE -1) then begin
      ; Set the mask bit for whopping fibers themselves
      fibermask[whopping] = fibermask[whopping] OR pixelmask_bits('WHOPPER')

      ; Set the mask bit for fibers near whopping fibers, excluding the whopping fibers themselves.  Note that a fiber could still have both
      ; WHOPPER and NEARWHOPPER set if it is both bright and near another bright fiber.
      wp = [whopping - 2 , whopping -1, whopping+1 , whopping+2]
      wp = wp[ where(wp GE 0 AND wp LT ny) ]
      fibermask[wp] = fibermask[wp] OR pixelmask_bits('NEARWHOPPER')
   endif

   ;-----
   ; Inherit any mask bits from the ccdmask, by setting the pixmask bits for anything that would be hit in a boxcar extraction
   if (keyword_set(ccdmask)) then begin
      for ibit=0, 31 do begin
         thischeck = extract_boxcar((ccdmask AND 2L^ibit) NE 0, xtrace, radius=2.5, nfiber=nfiber)
         pixelmask = pixelmask OR (2L^ibit * (thischeck GT 0))
      endfor
   endif

   ;-----------------------------------------------------------------------
   ;  This is a kludge to fix first and last column ???
   ;-----------------------------------------------------------------------
  if (configuration->extract_object_fixcolumns()) then begin
     image[0,*] = image[0,*]*0.7
     image[nx-1,*] = image[nx-1,*]*0.7
  endif

   ;
   ;  First we should attempt to shift trace to object flexure
   xnow = match_trace(image, invvar, xtrace,radius=radius, nfiber=nfiber)
   bestlag = median(xnow-xtrace)

   splog, 'Shifting traces by match_trace ', bestlag
   if (abs(bestlag) GT 1.0) then begin
      splog, 'WARNING: pixel shift is large!'
   endif

   highrej = 10  ; just for first extraction steps
   lowrej = 10   ; just for first extraction steps
                 ; We need to check npoly with new scattered light backgrounds
   npoly = 16    ; maybe more structure, lots of structure
   nrow = (size(image))[2]
   yrow = lindgen(nrow) 
   nfirst = n_elements(yrow)

   splog, 'Extracting frame '+objname+' with 4 step process'

   traceset2xy, widthset, xx, sigma2
   ntrace = (size(sigma2,/dimens))[1]
   wfixed = [1,1] ; Fit gaussian height + width (fixed center position)
   nterms = n_elements(wfixed)

   ;-----------------------------------------------------------------------
   ;  Now, subtract halo image and do final extraction with all rows
   ;-----------------------------------------------------------------------
   ; (6) Final extraction
   splog, 'Step 6: Final Object extraction'

   highrej = 8
   lowrej = 5
   wfixed = [1] ; Fit to height only (fixed width + center)
   nterms = n_elements(wfixed)
   reject = [0.2,0.2,0.2]
   npoly = 0

   ;TEMP MANGA SOLUTION -  Use extract_bundle_image for red, extract_image for blue
    if (color EQ 'red') then begin
       extract_bundle_image, image, invvar, xnow, sigma2, flux, fluxivar, proftype=proftype, wfixed=wfixed, ansimage=ansimage3, highrej=highrej, lowrej=lowrej, npoly=2L, $ 
        chisq=chisq, ymodel=ymodel, pixelmask=pixelmask, reject=reject, /relative, nperbun=20L, buffsize=8L, visual=visual, survey=survey, nbundle=nbundle, nfiber=nfiber, radius=radius, bundleid=bundleid
    endif
    
    if (color EQ 'blue') then begin
       extract_image, image, invvar, xnow, sigma2, flux, fluxivar, proftype=proftype, wfixed=wfixed, ansimage=ansimage3, highrej=highrej, lowrej=lowrej, npoly=2L, $
        chisq=chisq, ymodel=ymodel, pixelmask=pixelmask, reject=reject, /relative, visual=visual, survey=survey
    endif

; DRL kludge
 ;       testfile=str_replace(outname,'spFrame','spFrameUncal')
 ;       mwrfits, flux, testfile, /create
; end kludge

    ;write out ymodel from science extraction
    if keyword_set(writefiles) then begin
      modelfile = 'YMODEL-'+fileandpath(outname)
      mwrfits, image, modelfile, /create, /lscale
      mwrfits, ymodel, modelfile
      spawn, ['gzip','-f',modelfile], /noshell
    endif
    
   ; Replace the extracted fluxes with bbspec extractions
   if (keyword_set(bbspec)) then begin
      basisfile = 'spBasisPSF-*-'+strmid(sxpar(objhdr,'ARCFILE'),4,11)+'.fits'
      splog, 'Calling BBSPEC with '+basisfile
      tmproot = strmid(sxpar(objhdr,'FILENAME'),4,11)
      bbspec_extract, image, invvar, bbflux, bbfluxivar, basisfile=basisfile, ximg=xnow, ymodel=bb_ymodel, tmproot=tmproot, /batch 

      ; Deal with case of only the first few spectra being re-extracted...
      dims = size(bbflux,/dimens)
      flux[0:dims[0]-1,0:dims[1]-1] = bbflux
      fluxivar[0:dims[0]-1,0:dims[1]-1] = bbfluxivar

      ; Re-apply the rejection from the old extraction code...
      fluxivar *= ((pixelmask AND (pixelmask_bits('PARTIALREJECT') + pixelmask_bits('FULLREJECT'))) EQ 0)

      thisfile = 'ymodel-'+fileandpath(outname)
      mwrfits, bb_ymodel, thisfile, /create
      mwrfits, ymodel, thisfile
      splog, 'BBSPEC output file ymodel-'+thisfile
   endif

   ;----------------------------------------------------------------------
   ; Can we find cosmic rays by looking for outlandish ansimage ratios???
   ; a = where(ansimage[lindgen(ntrace)*nterms, *] LT (-2*ansimage[lindgen(ntrace)*nterms+1, *])

   ;------------------
   ; QA chisq plot for fit calculated in extract image (make QAPLOT ???)
   if keyword_set(VISUAL) then getwindow,/open
   xaxis = lindgen(n_elements(chisq)) + 1
   ymax = 2.*median(chisq)
   djs_plot, xaxis, chisq, xrange=[0,N_elements(chisq)], xstyle=1, yrange=[0,ymax], ystyle=1, xtitle='Row number',  ytitle = '\chi^2', title=plottitle+'Extraction chi^2 for '+objname

   djs_oplot, !x.crange, [1,1]
;x   djs_oplot, yrow, secondchisq[yrow], color='blue'
;x   djs_oplot, yrow, firstchisq[yrow], color='green'

   xyouts, 100, 0.05*!y.crange[0]+0.95*!y.crange[1], 'BLACK = Final chisq extraction'
;x   xyouts, 100, 0.08*!y.crange[0]+0.92*!y.crange[1], 'BLUE = Initial-scattered chisq extraction'
;x   xyouts, 100, 0.08*!y.crange[0]+0.89*!y.crange[1], 'GREEN = Initial chisq extraction'

   ;------------------
   ; Flat-field the extracted object fibers with the global flat
   ofluxivar=fluxivar
   divideflat, flux, invvar=fluxivar, fflat, /quiet, prediv=oldflux
 
   pixelmask = pixelmask OR ((fflat LT 0.5) * pixelmask_bits('LOWFLAT'))

   ;write out fluxes
   thisfile = 'outflux-'+fileandpath(outname)
   mwrfits, flux, thisfile, objhdr, /create
   mwrfits, oldflux, thisfile
   mwrfits, ofluxivar, thisfile

   ;------------------
   ; Apply heliocentric correction
   ; Correct LAMBDA, which is used to shift to vacuum wavelengths.

   helio=0.0
   ra = sxpar(objhdr, 'RA', count=ct_ra)
   dec = sxpar(objhdr, 'DEC', count=ct_dec)
   if (ct_ra NE 1 OR ct_dec NE 1) then splog, 'WARNING: Missing RA and/or DEC from header'

   ;--------------------------------------------------------
   ; Call standard proc to determine time-stamps

   get_tai, objhdr, tai_beg, tai_mid, tai_end

   ; Set TAI equal to the time half-way through the exposure
   ; If all these keywords are present in the header, they will be either type FLOAT or DOUBLE.  Note that SDSS will put NaN in the header for these values if they are unknown.
   if ( size(ra, /tname) NE 'INT' AND size(dec, /tname) NE 'INT' AND size(tai_mid, /tname) NE 'INT' AND finite(ra) AND finite(dec) AND finite(tai_mid) ) then begin
      helio = heliocentric(ra, dec, tai=tai_mid)
      splog, 'Heliocentric correction = ', helio, ' km/s'
      sxaddpar, objhdr, 'HELIO_RV', helio, ' Heliocentric correction (added to velocities)'
   endif else begin
      splog, 'WARNING: Header info not present to compute heliocentric correction'
   endelse
   if (size(tai_mid, /tname) EQ 'INT' OR finite(tai_mid) EQ 0) then begin
      splog, 'WARNING: Header info not present to compute airmass correction to sky level'
      tai_mid = 0
   endif

   ;------------------
   ; Shift to skylines and fit to vacuum wavelengths
   vacset = fitvacset(xarc, lambda, wset, arcshift, helio=helio, airset=airset)
   ; No longer make the following QA plot ???
   ;qaplot_skydev, flux, fluxivar, vacset, plugsort, color, title=plottitle+objname
   sxaddpar, objhdr, 'VACUUM', 'T', ' Wavelengths are in vacuum'

   ;------------------
   ;  If present, reconstruct superflat and normalize
   sxaddpar, objhdr, 'SFLATTEN', 'F', ' Superflat has not been applied'
   superfit = 0

   if keyword_set(superflatset) AND keyword_set(airset) then begin
     superfit = float(smooth_superflat(superflatset, airset, plottitle=plottitle+'Smooth superflat for '+objname))
     if keyword_set(superfit) then begin
       divideflat, flux, invvar=fluxivar, superfit, /quiet
       sxaddpar, objhdr, 'SFLATTEN', 'T', ' Superflat has been applied'
     endif
   endif  

   ; Save current pixelmask for later MaNGA use
   mangamask=pixelmask

   ;----------
   ; Add keywords to object header
   sxaddpar, objhdr, 'VERS2D', idlspec2d_version(), ' Version of idlspec2d for 2D reduction', after='VERSREAD'
   if (keyword_set(osigma)) then sxaddpar, objhdr, 'OSIGMA',  sigma, ' Original guess at spatial sigma in pix'
   sxaddpar, objhdr, 'PREJECT', reject[1], ' Profile area rejection threshold'
   sxaddpar, objhdr, 'LOWREJ', lowrej, ' Extraction: low rejection'
   sxaddpar, objhdr, 'HIGHREJ', highrej, ' Extraction: high rejection'
   sxaddpar, objhdr, 'SCATPOLY', npoly, ' Extraction: Order of scattered light polynomial'
   sxaddpar, objhdr, 'PROFTYPE', proftype, ' Extraction profile: 1=Gaussian'
   sxaddpar, objhdr, 'NFITPOLY', nterms, ' Extraction: Number of parameters in each profile'
   sxaddpar, objhdr, 'XCHI2', mean(chisq), ' Extraction: Mean chi^2'

   snvec = djs_median(flambda*sqrt(flambdaivar),1)
   if (color EQ 'blue') then begin
      filter = 'g'
      mag = plugsort.mag[1] 
   endif else begin
      filter = 'i'
      mag = plugsort.mag[3]
   endelse
   sncoeff = fitsn(mag, snvec, sncode='spreduce', filter=filter, redden=sxpar(objhdr,'REDDEN*'), sn2=sn2, dered_sn2=dered_sn2)

   sxaddpar,objhdr,'FRAMESN2', sn2[0], "(S/N)^2 at fidicial magnitude"
   sxaddpar,objhdr,'DEREDSN2', dered_sn2[0], "Extinction corrected (S/N)^2 (like quick redux)", after='FRAMESN2'
   sxaddpar,objhdr,'EQUINOX',2000.0,after='DEC'
   sxaddpar,objhdr,'RADECSYS', 'FK5', after='EQUINOX'
   sxaddpar,objhdr,'AIRMASS', float(tai2airmass(ra, dec, tai=tai_mid)) * (ct_ra EQ 1) * (ct_dec EQ 1) * (tai_mid GT 0), after='ALT'
   spadd_guiderinfo, objhdr
   sxaddpar, objhdr, 'EXTEND', 'T', after='NAXIS2'

   ;----------
   ; Write extracted, lambda-calibrated, sky-subtracted spectra to disk
   mwrfits, flux, outname, objhdr, /create;non-sky subtracted flux
   mwrfits, fluxivar, outname   ; non-sky subtracted inverse variance
   mwrfits, mangamask, outname ;pixel mask before silly sky issues
   mwrfits, vacset, outname    ;trace-set for wavelength sol; wset
   mwrfits, dispset, outname   ;trace-set for dispersion sol
   mwrfits, plugsort, outname  ;plugmap
   mwrfits, skyimg, outname    ;sky flux (BAD- old BOSS stuff)
   mwrfits, xnow, outname      ;x pos on CCD
   mwrfits, superfit, outname  ;superflat vector from quartz lamps
;   mwrfits, skystruct, outname
;   mwrfits, scatter, outname
;   if (keyword_set(do_telluric)) then mwrfits, telluricfactor, outname ; This array only exists for red frames.

   ;write arrays to a common wavelength grid and add extensions to fits file
   getCommonMangaWave, vacset, flux, fluxivar, mangamask, FILE=outname

   spawn, ['gzip', '-f', outname], /noshell

   heap_gc
   obj_destroy,configuration

   return   
 
 end