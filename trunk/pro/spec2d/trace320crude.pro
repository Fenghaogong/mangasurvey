;+
; NAME:
;   trace320crude
;
; PURPOSE:
;   Calling script to return 320 full traces using TRACE_CRUDE.
;
; CALLING SEQUENCE:
;   xset = trace320crude( fimage, invvar, [ ystart=, nmed=, $
;    xmask=, yset=, maxerr=, maxshifte=, maxshift0=, xerr=, maxdev=, ngrow=, $
;    fibermask=, cartid=, flathdr=, padding=, plottitle= ] )
;
; INPUTS:
;   fimage     - Image
;   cartid     - Cartridge ID from plugmap
;
; OPTIONAL INPUTS FOR TRACE320CEN:
;   ystart     - Y position in image to search for initial X centers; default
;                to the central row
;   nmed       - Number of rows to median filter around YSTART; default to 21
;   plottitle  -
;
; OPTIONAL INPUTS FOR TRACE_CRUDE:
;   flathdr    - FITS header for determining CARTID and MJD
;   invvar     - Inverse variance (weight) image
;   radius     - Radius for centroiding; default to 3.0
;   maxerr     - Maximum error in centroid allowed for valid recentering;
;                default to 0.2
;   maxshifte  - Maximum shift in centroid allowed for valid recentering;
;                default to 0.1
;   maxshift0  - Maximum shift in centroid allowed for initial row;
;                default to 0.5
;   padding    - number added to ndegree in loop to fix deviant centroids
;
; OPTIONAL INPUTS:
;   maxdev     - Maximum deviation of X in pixels; default to rejecting any
;                XPOS positions that deviate by more than 1.0 pixels from
;                a polynomial remapping of the centroids from other rows.
;   ngrow      - For each trace, replace all centroids within NGROW rows
;                of a bad centroid with the predicted centroid locations.
;                Default to 5.
;   fibermask  - Fiber status bits, set nonzero for bad status [NFIBER]
;
; OUTPUTS:
;   xset       - X centers for all traces
;
; OPTIONAL OUTPUTS:
;   yset       - Y centers for all traces
;   xerr       - Errors for XSET
;   xmask      - Mask set to 1 for good fiber centers, 0 for bad;
;                same dimensions as XSET.
;   fibermask  - (Modified.)
;
; COMMENTS:
;   Without djs_maskinterp, hot columns skew traces 
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   djs_maskinterp()
;   fibermask_bits()
;   ml_trace_crude()
;   ml_trace_fweight()
;   trace320cen()
;
; REVISION HISTORY:
;   13-Sep-1999  Written by David Schlegel, Princeton.
;    8-Jul-2001  Added djs_maskinterp call
;   05-Oct-2010  ASB added masking of rows with invvar all zero
;-
;------------------------------------------------------------------------------
function trace320crude, image, invvar, ystart=ystart, nmed=nmed, $
 xmask=xmask, radius=radius, yset=yset, maxerr=maxerr, maxshifte=maxshifte, $
 maxshift0=maxshift0, xerr=xerr, maxdev=maxdev, ngrow=ngrow, $
 fibermask=fibermask, cartid=cartid, flathdr=flathdr, padding=padding, $
 plottitle=plottitle, VISUAL=VISUAL, fiberparam=fiberparam

   if (~keyword_set(maxdev)) then maxdev = 1.0
   if (~keyword_set(ngrow)) then ngrow = 5
   if (~keyword_set(padding)) then padding=0
   if n_elements(survey) eq 0 then survey='MANGA'

   if (~keyword_set(radius)) then begin
     case strupcase(survey) of
        'MANGA': radius = 2.0
        'BOSS':  radius = 3.0
        else:  radius = 3.0
     endcase
   endif

   camname = strtrim(sxpar(flathdr, 'CAMERAS'),2)
   mjd = sxpar(flathdr, 'MJD')
   if (keyword_set(cartid) * keyword_set(camname) * keyword_set(mjd) EQ 0) then message, 'Must set CARTID, CAMERAS, MJD in flat header!'

   ;testing purposes - sets name of file to load
   data = (n_elements(data) eq 0) ? '' : data
   fname='TEST_FIBERPARAM'
   case strupcase(data) of 
      'ORIG': parname = 'manga_test.par'
      'VARY': parname = 'manga_test_vary.par'
      'VMISS': parname = 'manga_test_vary.par'
      'REAL': begin
          parname = 'manga_back_opFibers.par'
          fname='FIBERPARAM'
          end
      else:begin
        parname = 'manga_back_opFibers.par'
        fname='FIBERPARAM'
        end
   endcase

   case strupcase(survey) of 
      'MANGA': fiberparam = yanny_readone(djs_filepath(parname, root_dir=getenv('IDLSPEC2D_DIR'), subdir='opfiles'), fname)
      'BOSS': fiberparam = yanny_readone(djs_filepath('test.par', root_dir=getenv('IDLSPEC2D_DIR'), subdir='opfiles'), 'TEST_FIBERPARAM')
      
      ;assume default is BOSS
      else: fiberparam = yanny_readone(djs_filepath('opFibers.par', root_dir=getenv('IDLSPEC2D_DIR'), subdir='opfiles'), 'FIBERPARAM') 
   endcase

   if (~keyword_set(fiberparam)) then message, 'opFibers.par file not found!'
   i = where(fiberparam.cartid EQ cartid AND fiberparam.camname EQ camname AND fiberparam.mjd LE mjd, ct)
   if (ct EQ 0) then message, 'No match for this CARTID + MJD in opFibers.par!'
   isort = i[reverse(sort(fiberparam[i].mjd))]
   fiberparam = fiberparam[isort[0]]
   ; Assume the fiber bundles used are the first NBUNDLE ones...  ;change from .fiberspace to .bundlegap
   nbundle = (long(total(fiberparam.bundlegap NE 0)))[0]
   nfibspace = (long(total(fiberparam.fiberspace NE 0)))[0]
   nfiber = fiberparam.nfiber
   radius = fiberparam.radius
   if (total(fiberparam.bundlegap[0:nbundle-1] EQ 0) GT 0) then message, 'Some BUNDLEGAP parameters are zero!'
   if (total(fiberparam.fiberspace[0:nfibspace-1] EQ 0) GT 0) then message, 'Some FIBERSPACE parameters are zero!'

   ;compute bundle indicies for each fiber
   d=lindgen(fiberparam.totfiber)
   h=histogram(total(nfiber,/cumul)-1,min=0,/binsize,reverse_indices=ri)
   ind=ri[0:n_elements(h)-1]-ri[0]
   bundleid = d[ind]   
   fiberparam = jjadd_tag(fiberparam, 'nbundle',nbundle)
   fiberparam = jjadd_tag(fiberparam, 'nfibspace',nfibspace)
   fiberparam = jjadd_tag(fiberparam, 'bundleid',bundleid,/array)
   fiberparam = jjadd_tag(fiberparam, 'nend', total(nfiber,/cumul)-1, /array)
   nend = total(nfiber,/cumul)-1
   fiberparam = jjadd_tag(fiberparam, 'nstart', [0,nend[0:nbundle-2]+1], /array)   

   print, "SURVEY SAYS: "+survey+"!"
   print, "       Total Fibers: ", fiberparam.totfiber
   print, "       Number of Bundles: ", nbundle
   print, "       Number of Fiber Spacings (1 per bundle) : ",nfibspace
   print, "       Number of Fibers per Bundle: ", nfiber
   print, "       Radius of Extraction for each Bundle (in pixels): ", radius

;----------------------------------
   ; If INVVAR is set, then start by interpolating over bad pixels

   if (keyword_set(invvar)) then fimage = djs_maskinterp(image, (invvar LE 0), iaxis=0) $
      else fimage = image
  
;---------------------------------
   ; Find the X-centers in the row specified by YSTART (default is central row on CCD)
    ystart=2056
    xposition = getpeaks(fimage[*,ystart], numpeak=numpeak, xgood=xgood,visual=visual,survey=survey,ystart=ystart, nfiber=fiberparam.nfiber,$
        totfiber=fiberparam.totfiber,bundleid=bundleid,nbundle=nbundle,miss=miss,fsep=fiberparam.fiberspace)
   
    if n_elements(miss) ne 0 then print, 'MISSING FIBERS: ', miss.fiberid

  ;FLAG if xposition null
  if n_elements(xposition) eq 0 then begin
    splog, "WARNING: XPOSITION IS NULL - TRACE_CEN FAILED!"
      flag = 'No xposition array'
      return,0
  endif

 ;FLAG if xposition > ntrace
 if n_elements(xposition) gt fiberparam.totfiber then begin
    splog, 'WARNING: XPOSITION IS LARGER THAN TOTAL FIBER NUMBER!'
    return,0
 endif

  if keyword_set(VISUAL) then begin
    getwindow, next=fftrace
    window,fftrace
    !p.multi=[0,1,2]
    setcolors, /system_variables, /silent
    plot, fimage[*,2056], xra=[250,600],xstyle=1,/nodata,ytickformat='(A1)',tit='Survey: '+strupcase(survey)+'; CCD - Flat'
    ml_tvimage, fimage[250:600,1950:2150],/axis,yra=[1950,2150],xra=[250,600],/overplot
    oplot, xposition, fltarr(fiberparam.totfiber)+ystart, psym=2, color=!red
    p1=!p & x1=!x & y1=!y
    plot, fimage[*,2056], xra=[250,600],xstyle=1,psym=10,tit='Row '+strtrim(ystart,2)
    oplot, xposition, fltarr(fiberparam.totfiber)+round(max(fimage[*,2056])*0.1), psym=2, color=!red
    !p.multi=0
  endif

   ntrace = n_elements(xposition)
   if (NOT keyword_set(fibermask)) then fibermask = bytarr(ntrace)

;---------------------------------
   ; Trace out the centers of the rest of the rows (fix this radius)
   print, 'TRACING FIBERS ALONG THE REST OF THE CCD'
   splog, 'TRACING FIBERS ALONG THE REST OF THE CCD'
   
   xset = ml_trace_crude(fimage, invvar, xstart=xposition, ystart=ystart, radius=radius, yset=yset, maxerr=maxerr, maxshifte=maxshifte, maxshift0=maxshift0, xerr=xerr, nfiber=nfiber)
   xmask = xerr LT 990  ; =1 for good centers, =0 for bad

   ; Mask contributions from completely bad rows (ticket #1025)
   nullrow = total(invvar gt 0., 1) eq 0.
   wh_null = where(nullrow, n_null)
   if (n_null gt 0) then xmask[wh_null,*] = 0B

;--------------------------------------------------------------------
   ; Mark this trace as potentially bad (xgood[itrace] = 0)
   ; if either the initial extraction row had bad pixels,
   ; or the initial extraction row was off the left or right
   ; edge of the CCD.

   ncol = (size(invvar,/dimen))[0]
   if ncol GT 0 then begin
    ;loop over bundles
     nend = total(nfiber,/cumul)-1
     nstart = [0,nend[0:n_elements(nend)-2]+1]
     for nbund=0,nbundle-1 do begin
          ;loop over fibers within bundle
         for itrace=0, nfiber[nbund]-1 do begin
           ix1 = floor((xposition[nstart[nbund]:nend[nbund]])[itrace]-radius[nbund]) > 0L
           ix2 =  ceil((xposition[nstart[nbund]:nend[nbund]])[itrace]+radius[nbund]) < (ncol-1)
           if (ix1 GT ncol-1 OR ix2 LT 0) then begin
              ; Case where the initial extraction position was off the left or right edge of the CCD.
              (xgood[nstart[nbund]:nend[nbund]])[itrace] = 0
           endif else begin
              ; Check if any bad pixels at initial centroiding position.
              junk = where(invvar[ix1:ix2,ystart] LE 0, nbad)
              if (nbad GT 0) then begin
                tmp=xgood[nstart[nbund]:nend[nbund]]
                tmp[itrace] = 0
                xgood[nstart[nbund]:nend[nbund]]=tmp
              endif  
           endelse
         endfor
     endfor
   endif

;--------------------------------------------------
   ; Compare the traces in each row to those in row YSTART.
   ; Our assumption is that those centers should be a polynomial mapping
   ; of the centers from row YSTART.  Centers that are deviant from this
   ; mapping are replaced with the position predicted by this mapping.

   ny = (size(fimage, /dimens))[1]
   ndegree = 4 ; Five terms

   ;----------
   ; Loop to find all deviant centroids, and add these to the mask XMASK.
   ; XMASK=1 for good.

   for iy=0, ny-1 do begin
      xcheck = xgood AND xmask[iy,*] ; Test for good fiber & good centroid
      if (total(xcheck) GT ((ndegree+2) > (0.2 * ntrace))) then begin
        ; perform a polynomial to row iy, across the xpositions
        res = djs_polyfit(xposition, reform(xset[iy,*]), ndegree, variance=xcheck, yfit=xfit)
        xdiff = xfit - xset[iy,*]
        ibad = where(abs(xdiff) GT maxdev)
        if (ibad[0] NE -1) then xmask[iy,ibad] = 0

        ;overplot some traces and the fits to the xpositions
        if keyword_set(VISUAL) then begin
            index = findgen(6)*30+1970
            if ifany(index mod iy eq 0) then begin
               getwindow, set=fftrace
               !p=p1 & !x=x1 & !y=y1
               oplot, xset[findgen(6)*30+1970,*],yset[findgen(6)*30+1970,*], psym=2, color=!blue
               lineind = where(index mod iy eq 0,clineind)
               if clineind ne 0 then oplot, xfit, fltarr(fiberparam.totfiber)+iy, color=!green
               !p.multi=0
               ;getwindow,/open
               ;window,2
               ;plot, xdiff, psym=10, ytit='Position Fit-Xset(crude)', xtit='Fiber Number', tit='Row '+strtrim(iy,2)+': Residual of x-position fit with crude positions',charsize=2
             endif
        endif

      endif else begin
        xmask[iy,*] = 0 ; Too few good centroids in this row; mark all as bad
      endelse
   endfor

;--------------------------------------
   ; Smooth the bad centroids to NGROW adjacent rows (of the same trace)

   for itrace=0, ntrace-1 do begin
      xmask[*,itrace] = smooth( xmask[*,itrace]+0.0, 2*ngrow+1) EQ 1
   endfor

;----------------------------------------
   ; Loop to fix deviant centroids

   for iy=0, ny-1 do begin
      ixbad = where(xmask[iy,*] EQ 0, nbad)
      if (nbad GT 0 AND nbad LT ntrace-1-(ndegree+padding)) then begin
         ixgood = where(xmask[iy,*] EQ 1)

         if (!version.release LT '5.4') then begin
            coeff = polyfitw(xposition, xset[iy,*], xmask[iy,*], ndegree, xfit)
         endif else begin
            coeff = polyfitw(xposition, xset[iy,*], xmask[iy,*], ndegree, xfit, /double)
         endelse

         xset[iy,ixbad] = xfit[ixbad]
      endif
   endfor

;----------------------------------------------------------------------
   ;  This is new code to replace traces which have been 
   ;  affected by bad columns (or masked pixels in general) at the 
   ;  starting row.  Any start position which is offset from the true
   ;  position will produce systematic errors in all of the centroids
   ;  which are corrected above and based on xposition.
   ;  
   ;  This algorithm works as follows:
   ;  1) Select from the near eight neighbors (-4 to +4) traces
   ;       which have fewer than 100 bad pixels. (Checktrace)
   ;  2) Calculate the mean trace of the selected neighbors = meantrace
   ;  3) Calculate the offset of the problem trace with the selected neighbors
   ;           a) use only rows (goodrows) which have good centroids in 
   ;                        all selected neighbors AND in the problem trace.
   ;           b) calculate median offset w.r.t. neighbors in goodrows only
   ;           c) Use this mean offset to correct the meantrace to the correct 
   ;                zero point position

   problemtraces = where(xgood EQ 0, ct)
   nrow = (size(xset,/dimen))[1]
   tracenum = lindgen(ntrace)
   tmp_xpos = ml_trace_fweight(fimage, xset, yset, radius=radius, xerr=tmp_xerr, invvar=invvar, nfiber=nfiber)   

   xorig = xset
   badpix = total(tmp_xerr EQ 999,1) 

   if ct GT 0 then splog, 'Warning: Fixing traces: ', fix(problemtraces)
   for ii=0, ct-1 do begin
      itrace = problemtraces[ii] 

;----------------------------------------------------------------------------
;	Really simple minded loop to check for nearest 8 neighbors who might be suitable for substitution
      checktrace = -1
      for icheck = itrace-4 > 0, (itrace+4) < (ntrace-1) do begin

      ;-------------------------------------------------------------------
      ; accept only good traces (xgood), which are not itself, and have less than 100 bad centroids returned from ml_trace_fweight
         if xgood[icheck] AND icheck NE itrace AND badpix[icheck] LT 100 then begin
           if checktrace[0] EQ -1 then checktrace = icheck $
                else checktrace = [checktrace, icheck]
         endif
      endfor

      ;-------------------------------------------------------------------
      ;  Checktrace contains the selected good neighbors
      ncheck = n_elements(checktrace)

      ;---------------------------------------------------------------------
      ;  Need at least two good neighbors to perform correction
      if ncheck GE 2 then begin
        clean = total(tmp_xerr[*,checktrace] EQ 999,2) EQ 0
        meantrace = total(xset[*,checktrace],2)/ncheck 

        goodrows = where(clean AND tmp_xerr[*,itrace] NE 999, ngoodrows)

      ;---------------------------------------------------------------------
      ;  Require at least 100 common good rows to carry-on
        if ngoodrows GE 100 then begin
          xmask[goodrows,itrace] = 1
          xset_good = xset[goodrows,*]
          offset = tmp_xpos[goodrows,itrace] # replicate(1,ncheck) - xset_good[*,checktrace]
          shift =  mean(djs_median(offset,1))
          xset[*,itrace] = meantrace + shift
        endif else begin
           splog, 'Fiber ', fix(itrace), ' Only ', fix(ngoodrows), ' rows to adjust trace, skipping'
             xmask[itrace,*] = 0
        endelse

      endif else begin
           splog, 'Fiber ', fix(itrace), ' Only ', fix(ncheck), ' neighboring good fibers, skipping'
           xmask[itrace,*] = 0
      endelse

    endfor

;------------------------------------------
   ; Set FIBERMASK bit for any fiber with more than 20% of its positions
   ; --->let's change it to 45%, there are far too many traces flagged as
   ;      BADTRACE with so many bad columns and a 3 pixel radius!) ???
   ; masked, which includes any positions off the CCD.  Do not pay any
   ; attention to XGOOD, since that may indicate that a fiber is only
   ; bad near YSTART.

;in sims - fibermask filled with 2's here
   ibad = where(total(1-xmask, 1) GT 0.45*ny)
   if (ibad[0] NE -1) then fibermask[ibad] = fibermask[ibad] OR fibermask_bits('BADTRACE')

   ;plot differences between starting fiber positions and final fiber positions
   if keyword_set(VISUAL) then begin
      getwindow,/open
      med = median(xset[2056,*]-xposition,/even)
      plot, xset[2056,*]-xposition, charsize=2, xtit='Fiber Number', ytit='(Final - Initial) Fiber Position - Row 2056', tit='Survey: '+strupcase(survey),psym=2, yra=[0.0,-0.0]+med
      oplot, [0,1000], [0,0], color=!red
   endif 


   return, xset
end
;------------------------------------------------------------------------------
