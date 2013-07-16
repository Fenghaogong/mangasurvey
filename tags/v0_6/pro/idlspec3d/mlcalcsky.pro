;+
; function mlcalcsky
;
; This is the low-level routine that actually does the sky fitting.
; Based on BOSS skysubtract.pro with a fresh minty MaNGA flavor
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 01/25/2013
;   Last modified: 01/25/2013
;
; REVISION HISTORY:
;   v1: 25-Jan-2013  D. Law
;       Imported code from BOSS skysubtract.pro, start integrating
;       with MaNGA algorithms and calling
;   v2: 23-Jun-2013  D. Law
;       Revised to incorporate a call to use dispersion values
;       in final 2d fit.
;-

function mlcalcsky, objflux, objivar, wave, slitmap, objsub, objsubivar, $
 iskies=iskies, fibermask=fibermask, nord=nord, upper=upper, $
 lower=lower, maxiter=maxiter, pixelmask=pixelmask, thresh=thresh, $
 dispval=dispval, $
 npoly=npoly, relchi2set=relchi2set, $
 tai=tai, nbkpt=nbkpt, newmask=newmask

   if (size(objflux, /n_dimen) NE 2) then message, 'OBJFLUX is not 2-D'
   if (size(objivar, /n_dimen) NE 2) then message, 'OBJIVAR is not 2-D'

   dims = size(objflux, /dimens)
   ncol = dims[0]
   nrow = dims[1]

   ; Set default parameters (can't set any other way currently)
   if (n_elements(maxiter) EQ 0) then maxiter = 5
   if (NOT keyword_set(upper)) then upper = 10.0
   if (NOT keyword_set(lower)) then lower = 10.0
   if (NOT keyword_set(nbkpt)) then nbkpt = ncol
   if (NOT keyword_set(thresh)) then thresh = 4.0

   if (n_elements(fibermask) NE nrow) then fibermask = bytarr(nrow) 

   if ((size(slitmap, /dimens))[0] NE nrow) then $
    message, 'SLITMAP does not have same size as nrow'

   if ( (size(wave, /dimens))[1] NE nrow) then $
    message, 'WSET does not have same size as nrow'

   ;----------
   ; Find sky fibers, if iskies not specified, default is to 
   ; pick everything called 'SKY2' that is properly plugged and working
   iskies = (n_elements(iskies) eq 0) ? where((slitmap.ifuname eq 'SKY2') AND (slitmap.plugstatus eq 1) AND (fibermask EQ 0), nskies) : iskies
   nskies = (n_elements(nskies) eq 0) ? n_elements(iskies) : nskies

   splog, 'Number of sky fibers = ', nskies
   if (nskies EQ 0) then begin
      splog, 'ABORT: No valid sky fibers in SLITMAP'
      return, 0
   endif

   ; At each wavelength, don't reject more than 10% of the sky pixels
   ; per iteration.  No need to set GROUPSIZE=2, because the B-spline code
   ; re-sorts all the data into wavelength order first, which effectively
   ; transposes the data.
   groupsize = nskies
   maxrej = ceil(0.10*nskies)

   if NOT keyword_set(tai) then airmass = replicate(1.0, nrow) $
    else airmass = float(tai2airmass(slitmap.ra, slitmap.dec, tai=tai))

   minairmass = min(airmass, max=maxairmass)
   splog, (maxairmass GT 2.5) ? 'WARNING: ' : '', $
    'Airmass range = ', minairmass, maxairmass

   airmass_correction = replicate(1.0,ncol) # airmass

   splog, 'Warning: Disabling airmass terms in sky-sub!' ; ???
   airmass_correction[*] = 1

   skywave = wave[*,iskies]
   skyflux = objflux[*,iskies]
   skyivar = objivar[*,iskies]

   divideflat, skyflux, invvar=skyivar, airmass_correction[*,iskies]

   ;----------
   ; Mask any sky pixels where LOWFLAT or NEARBADPIXEL are set.
   if (keyword_set(pixelmask)) then begin
      qbad = ((pixelmask[*,iskies] AND pixelmask_bits('LOWFLAT')) NE 0) $
           OR ((pixelmask[*,iskies] AND pixelmask_bits('NEARBADPIXEL')) NE 0)
      ibad = where(qbad, nbad)
      splog, 'Discarding ', float(nbad)/n_elements(skywave), $
       ' (fractional) of the sky pixels as bad'
      if (nbad GT 0) then skyivar[ibad] = 0
   endif

   ;----------
   ; Sort sky points by wavelengths

   isort = sort(skywave)
   skywave = skywave[isort]
   skyflux = skyflux[isort]
   skyivar = skyivar[isort]

   ;----------
   ; Compute "supersky" with a spline fit.
   ; Use the EVERYN parameter to space the spline points according to
   ; the density of data points.

   bkpt = 0
   everyn = floor(2.*nskies/3) > 1
   sset = bspline_iterfit(skywave, double(skyflux), $
    invvar=double(skyivar), $
    nord=nord, everyn=everyn, bkpt=bkpt, $
    upper=upper, lower=lower, maxiter=maxiter, $
    maxrej=maxrej, groupsize=groupsize, $
    outmask=outmask, yfit=skyfit, requiren=2)

   if (NOT keyword_set(sset)) then begin
      splog, 'ABORT: Fit sky is all zeros'
      return, 0
   endif

   ;----------
   ; This code has been added to slightly reduce the number of
   ; breakpoints, from 3100 to NCOL, where we now space them such
   ; that there is an equal amount of S/N between each.
   ; This also gives better behavior near the boundaries.

   igood = where(skyivar GT 0 AND outmask NE 0, ngoodpix)
   minwave  = min(skywave[igood],max=maxwave)
   snsqrt = sqrt((skyfit * sqrt(skyivar) > 0))
   ipos = where(snsqrt GT 0, npos)
   gkern = gauss_kernel(2.0*nskies)

   if (npos GT n_elements(gkern) AND nbkpt GE 16) then begin

      ;----------
      ; Construct a vector with the summed (and smoothed) S/N

      snsqrt[ipos] = convol(snsqrt[ipos], gkern)
      snsum = snsqrt
      for i=1L, n_elements(snsqrt)-1 do snsum[i] = snsum[i] + snsum[i-1]

      ;----------
      ; Select break points with the same amount of S/N between each.
      ; We select specific, tabulated wavelengths, but then smooth the
      ; ones that we select to prevent digitization problems.

      iplace = long(snsum/max(snsum)*nbkpt) < (nbkpt-2)
      iplace = uniq(iplace)
      bkpt = skywave[iplace]
      bkpt[0] = minwave
      newnbk = n_elements(bkpt)
      i1 = 2*nord - 1
      i2 = newnbk - 2*nord
      width = 7 < (i2-i1+1)
      if (i1 LT newnbk AND i2 GT 0 AND width GT 1) then $
       bkpt[i1:i2] = (smooth(bkpt,width))[i1:i2]
      ii = where(bkpt GE minwave AND bkpt LE maxwave, newnbk)
      bkpt = bkpt[ii]

      ;----------
      ; Pad with (NORD-1) break points far to the left, and that
      ; many far to the right.  The factor of 10 in spacing is arbitrary!

      lowdiff = 10.0 * abs(skywave[igood[0]+1] - skywave[igood[0]])
      highdiff = lowdiff
      fullbkpt = [ bkpt[0] - lowdiff*(reverse(lindgen(nord-1))+1), bkpt, $
       bkpt[newnbk-1] + highdiff*(lindgen(nord-1)+1) ]

      ;----------
      ; Re-do the super-fit with the new break points.

      sset = bspline_iterfit(skywave, double(skyflux), $
       invvar=double(skyivar), $
       nord=nord, fullbkpt=fullbkpt, $
       upper=upper*1.5, lower=lower*1.5, maxiter=maxiter, $
       maxrej=maxrej, groupsize=groupsize, $
       outmask=outmask, yfit=skyfit, requiren=2)
 
      if (NOT keyword_set(sset)) then begin
         splog, 'ABORT: Fit sky is all zeros'
         return, 0
      endif

   endif

   fullfit = bspline_valu(wave, sset) 

   ;----------
   ; Re-do the super-sky fit using the variable PSF, if DISPSET is set.

   if (keyword_set(npoly) AND keyword_set(fullbkpt)and(keyword_set(dispval))) then begin
     splog,'DRL adding in dispval'
     
;      fullx2 = replicate(1.0,ncol) # findgen(nrow)
      fullx2=dispval
      x2 = (fullx2[*,iskies])[isort]

      sset2d = bspline_iterfit(skywave, double(skyflux), $
       invvar=double(skyivar*outmask), $
       nord=nord, npoly=npoly, fullbkpt=fullbkpt, $
       upper=1.5*upper, lower=1.5*lower, maxiter=maxiter, $
       maxrej=maxrej, groupsize=groupsize, $
       yfit=skyfit, x2=double(x2), xmin=0., xmax=nrow, requiren=2)

      if (keyword_set(sset2d)) then begin
         if (total(sset2d.coeff) EQ 0) then begin
            splog, 'WARNING: 2-D b-spline failed!'
            sset2d = 0
         endif
      endif

      if (keyword_set(sset2d)) then begin
         sset = sset2d
         fullfit = bspline_valu(wave, sset, x2=fullx2) 
      endif

      if (keyword_set(sset2d)) then begin
         if (max(abs(sset2d.coeff)) GT 0) then qgood2d = 1B
      endif 
      if (keyword_set(qgood2d)) then begin
         splog, 'Applying 2-dimensional sky b-spline'
         sset = sset2d
         fullfit = bspline_valu(wave, sset, x2=fullx2)
      endif else begin
         splog,' WARNING: 2-dimensional sky b-spline failed'
      endelse
   endif

   ;----------
   ; Sky-subtract the entire image

   objsub = objflux - float(fullfit) * airmass_correction

   ;----------
   ; Fit to sky variance (not inverse variance)

   posvar = where(skyivar GT 0)
   if (posvar[0] NE -1) then begin
      skyvariance = 1.0 / skyivar[posvar]
      skyvarset = bspline_iterfit(skywave[posvar], skyvariance, $
       invvar=skyivar[posvar], nord=nord, bkpt=bkpt, $
       upper=upper, lower=lower, maxiter=maxiter, $
       maxrej=maxrej, groupsize=groupsize, requiren=2)

      skyvarfit = bspline_valu(wave, skyvarset) * airmass_correction^2

      skyvarfit = (skyvarfit>0) * (objivar GT 0) ;???
   endif

   ;----------
   ; Store "super" sky information in a structure
   ; We can't name it, because it could change size each time

   skystruct = create_struct( $
    'ISKIES', iskies, $
    'WAVE', skywave, $
    'FLUX', skyflux, $
    'INVVAR', skyivar, $
    'FULLBKPT', sset.fullbkpt, $
    'COEFFS', sset.coeff)

   ;----------
   ; Now attempt to model variance with residuals on sky fibers.
   ; This is difficult since variance has noise, so only do this if there
   ; are at least 3 sky fibers.

   if (nskies GE 3 AND ngoodpix GT ncol $
    AND NOT keyword_set(novariance)) then begin

      ; We don't need to re-evaluate SKYFIT below, because we already
      ; have it in memory.
;      skyfit = (fullfit[*,iskies])[isort]

      ; Note that SKYFLUX and SKYFIT below are the sky and the fit
      ; where we've divided out the airmass correction first.
      skychi2 = (skyflux - skyfit)^2 * skyivar

      ; Bin according to the break points used in the supersky fit.

      if (keyword_set(npoly)) then npoly1 = npoly $
       else npoly1 = 1L
      nbin = N_elements(bkpt) - 1
      relwave = fltarr(nbin)
      relchi2 = fltarr(nbin)
      for ibin=0L, nbin-1L do begin
         ; Locate data points in this bin, excluding masked pixels
         ii = where(skywave GE bkpt[ibin] AND skywave LT bkpt[ibin+1] $
          AND skyivar GT 0, nn)

         if (nn GT 2 AND nn GT (npoly1+1)) then begin
            ; Find the mean wavelength for these points
            relwave[ibin] = total(skywave[ii]) / nn

            ; Find the mean relative chi^2, assuming gaussian statistics.
            ; But this evaluation is wrecked by any outliers.
;            relchi2[ibin] = total(skychi2[ii]) / (nn-1)

            ; The following evaluation looks at the 67th percentile of the
            ; points, which is much more robust.
            pos67 = ceil(2.*nn/3.) - 1
            tmpchi2 = skychi2[ii] * (1.0*nn / (nn - npoly1))
            relchi2[ibin] = tmpchi2[ (sort(tmpchi2))[pos67] ] 

            ; Burles counter of bin number...
;            print, format='("Bin ",i4," of ",i4,a1,$)', $
;             ibin, nbin, string(13b)

         endif
      endfor

      ; Trim to only those bins where we set RELCHI2

      ii = where(relwave NE 0, nbin)
      if (nbin GT 0) then begin
         relwave = relwave[ii]
         relchi2 = relchi2[ii]
 
         ;----------
         ; Spline fit RELCHI2, only for the benefit of getting a smooth function
         ; Also, force the fit to always be >= 1, such that we never reduce the
         ; formal errors.

         relchi2set = bspline_iterfit(relwave, relchi2, nord=3, $
           upper=30, lower=30, maxiter=maxiter, everyn=2, requiren=1)
         relchi2fit = bspline_valu(wave, relchi2set) > 1
         splog, 'Median sky-residual chi2 = ', median(relchi2)
         splog, 'Max sky-residual chi2 = ', max(relchi2)
      endif
   endif

   if (NOT keyword_set(relchi2fit)) then begin
      splog, 'WARNING: Too few sky fibers or pixels to model sky-sub variance!!'
      relchi2fit = 1
   endif

   ;----------
   ; Modify OBJSUBIVAR with the relative variance

   objsubivar = objivar
   if (keyword_set(skyvarfit) AND n_elements(relchi2fit) GT 1) then begin
      posvar = where(objivar GT 0)
      if (posvar[0] NE -1) then begin
        objvar = 1.0 / objivar[posvar]
        objsubivar[posvar] = $
         1.0 / (objvar + ((relchi2fit[posvar] > 1)-1.0)*skyvarfit[posvar])
      endif
   endif 

   ; Reselect the values of SKYIVAR from OBJSUBIVAR
   ; (Comment this out, since it's not actually needed below)
;   skyivar = (objsubivar[*,iskies])[*]
;   skyivar = skyivar[isort]

   ;----------
   ; Create the output pixel mask.

   if (keyword_set(pixelmask)) then newmask = pixelmask $
    else newmask = make_array(size=size(objflux), /long)

   ;----------
   ; If any pixels on the image are outside of the wavelength range
   ; covered by the "supersky", then the formal errors are infinite
   ; for those pixels after skysubtraction.  Set the mask bit 'NOSKY'
   ; and set SKYSUBIVAR=0.

   ii = where(skyivar GT 0, ni) ; Note that SKYWAVE is already sorted
   iout = where(wave LT skywave[ii[0]] OR wave GT skywave[ii[ni-1]])
   if (iout[0] NE -1) then objsubivar[iout] = 0.0

   if (iout[0] NE -1 AND keyword_set(newmask)) then $
    newmask[iout] = newmask[iout] OR pixelmask_bits('NOSKY')

   ;----------
   ; Set the BADSKYCHI mask bit at any wavelength where the relative chi^2
   ; for the sky fibers is greater than THRESH (and insist that the NOSKY
   ; bit is not set for that pixel).
   ; Also, look for the Red Monster (any bad region contiguous in wavelength).

   if (keyword_set(relchi2) AND keyword_set(newmask)) then begin
      newmask = newmask OR pixelmask_bits('BADSKYCHI') $
       * (relchi2fit GT thresh) $
       * ((newmask AND pixelmask_bits('NOSKY')) EQ 0)

      if (arg_present(relchi2set)) then $
       redmonster, relwave, relchi2, wave, thresh=thresh, pixelmask=newmask
   endif

   ;----------
   ; Set the mask bit 'BRIGHTSKY' for any pixels where the sky
   ; level is more than the (sky-subtracted) object flux + 10 * error,
   ; and where the sky is greater than 2.0 times a median sky.
   ; Grow this mask by 1 neighboring pixel in each direction.

   if (keyword_set(newmask)) then begin
      ; Compute a median sky vector for each fiber
      medsky = 0 * fullfit
      for irow=0L, nrow-1L do $
       medsky[*,irow] = djs_median(fullfit[*,irow], width=99, $
        boundary='reflect')

      qbright = (objsubivar NE 0) $
       AND (fullfit GT objsub + 10.0 / sqrt(objsubivar + (objsubivar EQ 0))) $
       AND (fullfit GT 2.0 * medsky)
      qbright = convol(qbright, [1,1,1], /center, /edge_truncate)
      ibright = where(qbright)
      if (ibright[0] NE -1) then $
       newmask[ibright] = newmask[ibright] OR pixelmask_bits('BRIGHTSKY')
   endif

   return, skystruct

end
