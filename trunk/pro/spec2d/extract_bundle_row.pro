;+
; NAME:
;   extract_bundle_row
;
; PURPOSE:
;   Fit the fiber profiles and background in a single row with least
;   squares, chunking it up by bundle, using Gaussians plus polynomial
;   background.
;
;   Based upon extract_row.pro, which was a call to extract_row.c
;
; NOTES ON CONVERSION TO EXTRACT_BUNDLE_ROW (A. Bolton, Utah, 2011Aug):
;   The extraction guts of this program are entirely new,
;   written in pure IDL rather than in IDL-wrapped C.
;   However, I have attempted (1) to preserve most of the same
;   outward functionality in terms of inputs and outputs, at least
;   those that are relevant to the BOSS case, and (2) to retain
;   the same rejection-iteration behavior, so that the same sorts
;   of cases trigger rejection as before.  I did not parse all
;   the details of the rejection logic.  One thing I know is that
;   the rejection in current form would be more efficient if it were
;   done on each bundle separately, rather than on the entire row.
;   But that would require a fine-grained refactoring that would risk
;   breaking the rejection logic, so I didn't attempt it.
;
;   Some deprecated (and even retained) keyword input/output variables
;   may have unexpected behavior is used in a context other than the
;   current BOSS-implementation calls in EXTRACT_IMAGE.
;
;   I have retained a fair bit of the original code in commented-out
;   form, for forensic purposes.  Anything removed during the
;   conversion to bundle form is commented out with a ";#".
;
; CALLING SEQUENCE:
;; First four positional arguments inherited from extract_row:
;   ans = extract_row( fimage, invvar, xcen, sigma,
;; Next keyword arguments new to extract_bundle_row:
;; (skew and kurt not currently enabled)
;              [nperbun=, buffsize=, skew=, kurt=,
;; Next keyword retained but reinterpreted:
;              npoly=, 
;; Next keywords retained in more or less the same role (I think):
;              maxiter=, lowrej=, highrej=, niter=, reducedChi=, xvar=,
;              mask=, relative=, diagonal=, pixelmask=, reject=,, ymodel=
;; Next keyword arguments deprecated:
;              fscat=, proftype=, wfixed=, inputans=, iback=, bfixarr=,
;              fullcovar=, wfixarr=, squashprofile=, whopping=,
;              wsigma=, nBand=  ])
;
; INPUTS:
;   fimage     - Vector [nCol]
;   invvar     - Inverse variance [nCol]
;   xcen       - Initial guesses for X centers [nFiber]
;   sigma      - Sigma of gaussian profile; (scalar or [nFiber])
;
; OPTIONAL KEYWORDS (new):
;   nperbun    - Number of fibers per bundle (default is 20).
;   buffsize   - Pixel buffer size on the very outside of the image
;                (default is 8)
;   npoly      - order of polynomial background in bundle, default=2 (linear).
;
; OPTIONAL KEYWORDS (retained):
;   maxiter    - Maximum number of profile fitting iterations; default to 10
;   lowrej     - Negative sigma deviation to be rejected; default to 5
;   highrej    - Positive sigma deviation to be rejected; default to 5
;   reject     - Three-element array setting partial and full rejection
;                thresholds for profiles; default [0.2, 0.6, 0.6].
;                When there is less than REJECT[2] of the area is left,
;                  then drop fitting of all higher-order terms.
;                When there is less than REJECT[1] of the area is left,
;                  then the pixel is rejected (inverse variance is set to 0).
;                When there is less than REJECT[0] of the area is left,
;                  then assume that there's no fiber there, and don't fit
;                  for that fiber at all.
;   relative   - Set to use reduced chi-square to scale rejection threshold
;
; DEPRECATED KEYWORDS:
;   inputans   - Input fit, excluding background and whopping terms
;                [ncoeff*nFiber]
;                The array is sorted as follows:
;                  [ncoeff] values for fiber #0
;                   ...
;                  [ncoeff] values for fiber #(nFiber-1)
;   squashprofile - ???
;   nband      - Band-width of covariance matrix of fiber profiles: default 1
;   whopping   - X locations to center additional "whopping" terms to describe
;                the exponentail tails of flux near bright fibers; default
;                to -1, which means not to use any such terms.
;   wsigma     - Sigma width for exponential whopping profiles; default to 25
;
; MODIFIED INPUTS (OPTIONAL, retained):
;   xvar       - X values of fimage and invvar; default is findgen(NX).
;   mask       - Image mask: 1=good, 0=bad [NX]
;   pixelmask  - Bits set for each fiber due to extraction rejection
;                [nFiber]
;
; MODIFIED INPUTS (OPTIONAL, deprecated):
;   wfixed     - Array to describe which parameters to fix in the profile;
;                0=fixed, 1=float; default to [1].
;                The number of parameters to fit per fiber is determined
;                this way; e.g. nCoeff = n_elements(wfixed), so the default
;                is to fit only 1 parameter per fiber.  For the (default)
;                Gaussian profile, this is the height of the Gaussian.
;                Note that WFIXED is used to build the array WFIXARR.
;   iback      - 1D array of input background coeff 
;                (needed if fixed parameters are non-zero)
;   bfixarr    - 1D integer array to specify which terms of the background
;                coefficients to fix; 0=fixed, 1=float.
;   wfixarr    - 1D integer array to specify which parameters in the full fit
;                to fix; 0=fixed, 1=float.
;                The array is sorted as follows:
;                  [ncoeff] values for fiber #0
;                   ...
;                  [ncoeff] values for fiber #(nFiber-1)
;                  [npoly] values for the background polynomial terms
;                  [whoppingct] values for the whopping terms
;
; OUTPUTS (reinterpreted):
;   ans        - Output fit [nFiber]
;                (extracted Gaussian profile amplitudes)
;
; OPTIONAL OUTPUTS (retained):
;   ymodel     - Evaluation of best fit [nCol]
;   diagonal   - 1D diagonal of covariance matrix
;   niter      - Number of rejection iterations performed
;   reducedChi - Reduced chi ???
;
; OPTIONAL OUTPUTS (deprecated):
;   fscat      - Scattered light contribution in each fiber [nFiber]
;   fullcovar  - 2D covariance matrix.  This is a symmetric matrix, and we
;                only fill the lower triangle.  Computing this increases CPU
;                time by a factor of 2 or 3.
;
; REVISION HISTORY:
;    8-Aug-1999  extract_row Written by Scott Burles, Chicago 
;    Sep 2010    converted to extract_bundle_row by Adam S. Bolton, Utah
;    Aug 2011    attempted documentation cleanup, Adam S. Bolton, Utah
;-
;------------------------------------------------------------------------------
function extract_bundle_row, fimage, invvar, xcen, sigma, ymodel=ymodel, $
 fscat=fscat, proftype=proftype, wfixed=wfixed, inputans=inputans, $
 iback=iback, bfixarr=bfixarr, xvar=xvar, mask=mask, relative=relative, $
 squashprofile=squashprofile, diagonal=p, fullcovar=fullcovar, $
 wfixarr=wfixarr, npoly=npoly, maxiter=maxiter, $
 lowrej=lowrej, highrej=highrej, niter=niter, reducedChi=reducedChi, $
 whopping=whopping, wsigma=wsigma, pixelmask=pixelmask, reject=reject, $
 oldreject=oldreject, nband = nband, contribution=contribution, $
 nperbun=nperbun, buffsize=buffsize, skew=skew, kurt=kurt, VISUAL=VISUAL, SURVEY=survey,$
 NBUNDLE=nbundle, NFIBER=nfiber, RADIUS=radius, BUNDLEID=bundleid

   on_error, 0
   compile_opt idl2
   
   ; Need 4 parameters
   if (N_params() LT 4) then message, 'Wrong number of parameters'

   ntrace = n_elements(xcen)
   nx = n_elements(fimage)

   if (n_elements(sigma) NE ntrace) then begin
      sigma1 = sigma[0]
      sigma = xcen*0.0 + sigma1
   endif 

   if (n_elements(npoly) EQ 0) then npoly = 2L ; order of background is now per bundle!
   if (~keyword_set(nband)) then nband = 1L
   if (~keyword_set(maxiter)) then maxiter = 10
   if (~keyword_set(highrej)) then highrej = 15.0
   if (~keyword_set(lowrej)) then lowrej = 20.0
   if (~keyword_set(wfixed)) then wfixed = [1]
   if (~keyword_set(proftype)) then proftype = 1 ; Gaussian
   relative = keyword_set(relative) 
   squashprofile = keyword_set(squashprofile) 
   if (~keyword_set(wsigma)) then wsigma = 25.0
   if (~keyword_set(nperbun)) then nperbun = 20L
   if (~keyword_set(buffsize)) then buffsize = 8L

   ; Here we want a three element array where both are between 0 and 1, and the first is larger than the second.  
   ; The first threshold sets the minimum area required to perform single parameter profile fitting
   ; The second threshold is the minimum area required not to reject the pixel in the final extracted spectrum.
   ; The third parameter is the area required in the profile fit containing good pixels to do a full fit.

   if (n_elements(reject) NE 3) then reject = [0.2, 0.6, 0.6]

   if (n_elements(pixelmask) NE ntrace OR size(pixelmask,/tname) NE 'LONG') then pixelmask = lonarr(ntrace)

   if (~keyword_set(whopping)) then whopping = -1
   if (whopping[0] EQ -1) then whoppingct = 0L else whoppingct = n_elements(whopping)

   if (~keyword_set(xvar)) then xvar = findgen(nx) $
    else if (nx NE n_elements(xvar)) then message, 'Number of elements in FIMAGE and XVAR must be equal'

   if (~keyword_set(mask)) then mask = bytarr(nx) + 1b $
    else if (nx NE n_elements(mask)) then message, 'Number of elements in FIMAGE and MASK must be equal'

   ncoeff = n_elements(wfixed)

   if (nx NE n_elements(invvar)) then message, 'Number of elements in FIMAGE and INVVAR must be equal'

   ;----------
   ; Allocate memory for the C subroutine.
   ymodel = fltarr(nx)

   ;----------
   ; Test which points are good
   qgood = invvar GT 0.0 AND mask NE 0
   igood = where(qgood, ngood)

   ;----------
   ; Set the following variables before any possible RETURN statement.
   reducedChi = 0.0
   niter = 0
   ans = fltarr(ntrace)
   p = fltarr(ntrace)

   if (ngood EQ 0) then return, ans

   ;----------
   ; Check that XCEN is sorted in increasing order with separations of at least 3 pixels.

   junk = where(xcen[0:ntrace-2] GE xcen[1:ntrace-1] - 3, ct)
   if (ct GT 0) then splog, 'XCEN is not sorted or not separated by greater than 3 pixels.'

   ;----------
   ; Build the fixed parameter array if it was not passed.
      wfixarr = lonarr(ntrace) + 1

   ;-----------
   ; PUTTING IN LOOP OVER BUNDLES:

   ; Find number of bundles:
   nbun = nbundle     

    ; Find breaks between the bundles and set limits of pixels to be associated with each bundle:
    t_hi = total(nfiber,/cumul)-1  
    t_lo = [0,t_hi[0:nbun-2]+1]

    xc_lo = xcen[t_lo]
    xc_hi = xcen[t_hi]
    if (nbun gt 1) then begin
       midpoints = 0.5 * (xc_lo[1:nbun-1] + xc_hi[0:nbun-2])
       jmax = floor(midpoints)    ;jmax, min are xcen max and mins for the different blocks, midpoints is midpoint xcen between neighboring blocks
       jmin = jmax + 1
       jmin = [(floor(xc_lo[0]) - buffsize) > 0, jmin]
       jmax = [jmax, (ceil(xc_hi[nbun-1]) + buffsize) < (nx - 1)]
    endif else begin
       jmin = (floor(xc_lo) - buffsize) > 0
       jmax = (ceil(xc_hi) + buffsize) < (nx - 1)
    endelse

    ;print, 'xc_lo, xc_hi, mid, jmin, jmax for block 1,2'
    ;print, columnize(xc_lo[17:20], xc_hi[17:20], [midpoints[17:19],0], jmin[17:20], jmax[17:20])

    ; Mask outside the block-fitting zone:
    mask = mask * (xvar ge min(jmin)) * (xvar le max(jmax))

  ; The loop over blocks:
     totalreject = 0
     partial = lonarr(ntrace)
     fullreject = lonarr(ntrace)
     finished = 0

    if keyword_set(VISUAL) then begin
      blocks = ['127_2', '19_5', '5" cal', '127_2', '3" cal', '61_1', '19_4', '127_2', '61_1', '19_3', '127_2', '3" cal', '2" cal', '127_1', $
      '19_2', '127_1', '5" cal', '127_1', '19_1', '2" cal', '127_1'] 
      h=fimage#(fltarr(2)+1)
      h=bytscl(h,min=min(fimage),max=max(fimage))
      getwindow,/open, next=rowwin
      setcolors, /system_variables, /silent
      ml_tvimage, h[jmin[0]:jmax[nbun-1],*], /axis, xra=[jmin[0], jmax[nbun-1]]
      xyouts, jmin, fltarr(nbun)+1.5, strtrim(sindgen(nbun)+1,2), charsize=1,color=!red
      for line=0,nbun-1 do oplot, fltarr(2)+jmin[line], [0,2], color=!blue, thick=2
      p1=!p & x1=!x & y1=!y
      ;window, /free, /pixmap, xsize=!d.x_size, ysize=!d.y_size
      ;pixid = 0
      ;device, copy=[0,0,!d.x_size,!d.y_size,0,0,rowwin]
      getwindow,/open, next=subrowwin
    endif

   while(finished NE 1) do begin 

      ; Apply current mask to invvar:
      workinvvar = (FLOAT(invvar * mask))

      for ibun = 0L, nbun-1 do begin

         nperbun = nfiber[ibun]            ; new MANGA code  
         fiberbase = findgen(nfiber[ibun]) ; new MANGA code

        ; Make the extracting basis:
         workxcen = xcen[t_lo[ibun]:t_hi[ibun]]                                             ;working xcenters for block i
         worksigma = sigma[t_lo[ibun]:t_hi[ibun]]                                           ;working sigmas of Gaussian profile for block i  [size nperbun]
         workpar = (transpose([[workxcen], [worksigma]]))[*]                                ;working parameters combined into 1-vector  -- xcen, sig, xcen, sig, xcen, sig, etc.....
         workx = xvar[jmin[ibun]:jmax[ibun]]                                                ;working x values of block i , findgen index
         workxpoly = 2. * (workx - min(workx)) / (max(workx) - min(workx)) - 1.             ;working polynomial range -1,1 for block i
         workbasis = gausspix(workx, workpar)                                               ;   size [xsize of block i, nperbun+npoly_background]     
         if (npoly gt 0) then workbasis = [[workbasis], [flegendre(workxpoly, npoly)]]      ;working basis for the polynomial  npoly is the order of the background per block 
         workimage = fimage[jmin[ibun]:jmax[ibun]]                                          ;working image for block i 

         if keyword_set(VISUAL) then begin
            wset, rowwin
            !p=p1 & !x=x1 & !y=y1
            ;device, copy=[0,0,!d.x_size,!d.y_size,0,0,pixid]
            if ibun ge 1 then oplot, fltarr(2)+jmin[ibun-1], [0,2], color=!blue, thick=2
            oplot, fltarr(2)+jmin[ibun],[0,2], color=!yellow, thick=2
            oplot, fltarr(2)+jmax[ibun]+1, [0,2], color=!yellow, thick=2
            wset, subrowwin
            erase
            ml_tvimage, workimage#(fltarr(2)+1), /axis, xra=minmax(workx), axkeywords={title:'Block '+strtrim(ibun+1,2)+': '+strtrim(blocks[ibun],2)+'; Fibers '+string(t_lo[ibun]+1,f='(I3)')+'-'+string(t_hi[ibun]+1,f='(I3)'),charsize:2,xtitle:'x /pixel'}
            ;oplot, workx, workbasis[*,0], psym=10, color=!red
            ;oplot, workx, workbasis[*,1], psym=10, color=!green
            oplot, workxcen, fltarr(n_elements(workxcen))+1.8, psym=2, color=!red   ;x-centers of the fibers in block i
         endif

        ; Associate pixels with fibers:
         pixelfiber = round(interpol(fiberbase, workxcen, workx))

        ; Determine the good area fraction for each fiber, limiting to pixels "associated with" that fiber:
         bworkinvvar = workinvvar[jmin[ibun]:jmax[ibun]]
         goodarea = 0. * fltarr(nperbun)
         for ijk = 0L, nperbun-1 do begin
            totalarea = total((pixelfiber eq ijk) * (workbasis[*,ijk]))
            totalgood = total((pixelfiber eq ijk) * (bworkinvvar gt 0.) * (workbasis[*,ijk]))
            if (totalarea gt 0.) then goodarea[ijk] = totalgood / totalarea
         endfor

        ; Populate the rejection masks:
         partial[t_lo[ibun]:t_hi[ibun]] = goodarea lt reject[2]
         fullreject[t_lo[ibun]:t_hi[ibun]] = goodarea lt reject[1]
        
        ; Remove fitting for any fiber with area less than reject[0]:
         use_component = where(goodarea ge reject[0], n_use)
        ; Don't bother fitting at all if there are no unmasked fibers:
        if (n_use gt 0) then begin
        
                 if (npoly gt 0) then use_component = [use_component, nperbun + lindgen(npoly)] ;1 for each fiber and 2 for the background component fit
                 nfit_this = n_elements(use_component)
                 workbasis = workbasis[*,use_component]
        
        ; The extraction steps for block i:
                 itworkbasis = transpose(workbasis * (bworkinvvar # replicate(1., nfit_this)))
                 icovar = itworkbasis # workbasis  ;coefficients to the linear equations
                 beta = itworkbasis # workimage    ; right-hand side of equations
                 workidx = lindgen(nfit_this)
                 extractivar = icovar[workidx,workidx]
                 la_choldc, icovar, status=cholstat     ;cholesky decomposition
                 if (cholstat eq 0) then begin
                    workcoeffs = la_cholsol(icovar, beta)   ;solutions to the system of linear equations [u,v,w etc...]
                 endif else begin
                    workcoeffs = replicate(0., nfit_this)    ;array not positive definite
                    extractivar[*] = 0.
                 endelse
                 workmodel = workbasis # workcoeffs        ; model profile fits to the fibers in block i
                 if keyword_set(VISUAL) then begin
                    mm = minmax(workmodel)
                    iwm = interpol([0.0,1.0], float([floor(mm[0]),ceil(mm[1])]), workmodel)
                    oplot, workx, iwm, psym=10, color=!blue    ;model fits
                    iwc = interpol([1.0,1.5], float([floor(min(workcoeffs[0:nperbun-1])),ceil(max(workcoeffs))]), workcoeffs[0:nperbun-1])
                    oplot, workxcen, iwc, psym=symcat(17), color=!orange  ;relative flux between fibers within block i
                 endif
                 
        ; Unpack results, accounting for masked fibers:
                 fullcoeffs = replicate(0., nperbun + npoly)
                 fullexivar = replicate(0., nperbun + npoly)
                 fullcoeffs[use_component] = workcoeffs
                 fullexivar[use_component] = extractivar
        
        ; Populate the model:
                 ymodel[jmin[ibun]:jmax[ibun]] = workmodel
                 ans[t_lo[ibun]:t_hi[ibun]] = fullcoeffs[0:nperbun-1]
                 p[t_lo[ibun]:t_hi[ibun]] = sqrt(fullexivar[0:nperbun-1])
                 
        ;         if ibun eq 18 then begin
        ;          print, 'First Fiber: sigma, wcoeffs, owc, % change:', worksigma[0], workcoeffs[0], owc[0], (workcoeffs[0]-owc[0])/owc[0]*100.
        ;          print, 'Mid   Fiber: sigma, wcoeffs, owc, % change:', worksigma[nperbun/2], workcoeffs[nperbun/2], owc[nperbun/2], (workcoeffs[nperbun/2]-owc[nperbun/2])/owc[nperbun/2]*100.
        ;          print, 'Last  Fiber: sigma, wcoeffs, owc, % change:', worksigma[nperbun-1], workcoeffs[nperbun-1], owc[nperbun-1], (workcoeffs[nperbun-1]-owc[nperbun-1])/owc[nperbun-1]*100.
        ;          owc=workcoeffs
        ;         endif
                 
        endif
      endfor

      ; Apply fullreject mask to output flux and inverse-error:
      p = p * (fullreject eq 0)
      ans = ans * (fullreject eq 0)
      
      if keyword_set(VISUAL) then begin
        wset, rowwin
        !p=p1 & !x=x1 & !y=y1
        mm=minmax(ans)
        ians = interpol([0.0,1.0],float([floor(mm[0]),ceil(mm[1])]), ans)
        oplot, xcen,ians, color=!green, psym=10   ;relative flux across the slit for the given row on CCD
      endif


      diffs = (fimage - ymodel) * sqrt(workinvvar) 
      chisq = total(diffs*diffs)
      countthese = total(wfixarr)
      reducedChi = chisq / ((total(mask) - countthese) > 1)
      errscale = 1.0
      if (relative) then errscale = sqrt((chisq/total(mask)) > 1.0)

      finished = 1
      ;the rejection algorithm
      if (~keyword_set(oldreject)) then begin
         
         goodi = where(workinvvar GT 0)
         if (goodi[0] NE -1) then begin 
           indchi = diffs[goodi] / (highrej*errscale)
           neg = where(indchi LT 0)
           if (neg[0] NE -1) then indchi[neg] = -indchi[neg] * highrej / lowrej

           badcandidate = where(indchi GT 1.0, badct)
           if (badct GT 0) then begin
             finished = 0
             ; find groups
             padbad = [-2, badcandidate, nx + 2]
             startgroup = where(padbad[1:badct] - padbad[0:badct-1] GT 1, ngrp)
             endgroup = where(padbad[2:badct+1] - padbad[1:badct] GT 1)
             for i=0, ngrp - 1 do begin
               worstdiff = max(indchi[startgroup[i]:endgroup[i]],worst)
               mask[goodi[badcandidate[worst+startgroup[i]]]] = 0
             endfor
             totalreject = totalreject + ngrp
           endif
         endif
      endif else begin
        badhigh = where(diffs GT highrej*errscale, badhighct)
        badlow = where(diffs LT -lowrej*errscale, badlowct)
        if (badhighct GT 0) then begin
           mask[badhigh] = 0
           totalreject = totalreject + badhighct
           finished = 0
        endif 
        if (badlowct GT 0) then begin
           mask[badlow] = 0
           totalreject = totalreject + badlowct
           finished = 0
        endif

        diffs = (fimage - ymodel) * sqrt(invvar) 
        if (finished EQ 0) then begin
           ii = where(diffs GE -lowrej*errscale AND diffs LE highrej*errscale)
           if (ii[0] NE -1) then mask[ii] = 1 ; These points still good
        endif
      endelse

      niter = niter + 1
      if (niter EQ maxiter) then finished = 1
      
   endwhile

   ;----------
   ; Add bits to PIXELMASK
   pixelmask = pixelmask OR (pixelmask_bits('PARTIALREJECT') * fix(partial))
   pixelmask = pixelmask OR (pixelmask_bits('FULLREJECT') * fix(fullreject))

   if keyword_set(VISUAL) then begin
    getwindow,kill=rowwin
    getwindow,kill=subrowwin
   endif

   return, ans
   
end

