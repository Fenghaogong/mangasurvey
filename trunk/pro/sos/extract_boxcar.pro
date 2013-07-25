;+
; NAME:
;   extract_boxcar
;
; PURPOSE:
;   Extract the total flux within a boxcar window at many positions.
;
; CALLING SEQUENCE:
;   fextract = extract_boxcar( fimage, xcen, ycen, [radius=radius] )
;
; INPUTS:
;   fimage     - Image
;   xcen       - Initial guesses for X centers
;   ycen       - Y positions corresponding to "xcen" (expected as integers)
;
; OPTIONAL KEYWORDS:
;   radius     - Radius of extraction; default to 3.0
;
; OUTPUTS:
;   fextract   - Extracted flux at each position specified by (xcen, ycen)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;   Dynamic link to extract_boxcar.c
;
; REVISION HISTORY:
;   24-Mar-1999  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
function extract_boxcar, fimage, xcen, ycen, radius=radius, idl=idl, VISUAL=VISUAL, SURVEY=survey, DEBUG=DEBUG,nbundle=nbundle,nfiber=nfiber

on_error,0
compile_opt idl2

   ; Need 2 parameters
   if (N_params() LT 2) then begin
      print, 'Syntax - fextract = extract_boxcar( fimage, xcen, [ycen, radius=radius] )'
      return, -1
   endif
   if (~keyword_set(radius)) then radius = 3.0

   ;build y centers if not input
   if (~keyword_set(ycen)) then begin
	      ycen=xcen
        ndim = (size(xcen))[0]
        if (ndim EQ 1) then ycen = findgen(N_elements(xcen)) $
          else if (ndim EQ 2) then begin
            npix = (size(xcen))[1]
            nTrace = (size(xcen))[2]
            for i=0,nTrace-1 do ycen[*,i] = findgen(npix)
          endif else message, 'xcen is not 1 or 2 dimensional'
   endif
         
   if (N_elements(xcen) NE N_elements(ycen)) then message, 'Number of elements in XCEN and YCEN must be equal'

   nx = (size(fimage))[1]
   ny = (size(fimage))[2]
   ncen = N_elements(xcen)

;   if (min(ycen) LT 0 OR max(ycen) GT ny-y) then message, 'YCEN contains values out of range'
   fextract = float(0 * xcen)

    ;expand radius to full array for variable extraction
    d=lindgen(total(nfiber))
    h=histogram(total(nfiber,/cumul)-1,min=0,/binsize,reverse_indices=ri)
    ind=ri[0:n_elements(h)-1]-ri[0]
    radius = radius[ind]##(fltarr(ny)+1) 

   ;does the extraction, either inside IDL or in an external library call ;  extracts the flux for a fiber at each xcenter +- radius pixels (3, 6 pixel span)
   if keyword_set(idl) then begin
      left = xcen - radius
      right = xcen + radius
      fextract = extract_asymbox2(fimage, left, right, ycen)
   endif else begin
     soname = filepath('libspec2d.'+idlutils_so_ext(), root_dir=getenv2('MANGAPIPE_DIR'), subdirectory='lib')
     result = call_external(soname, 'ml_extract_boxcar', nx, ny, float(fimage), float(radius), ncen, float(xcen), long(ycen), fextract)
   endelse

  ;visualize the approx extraction
  if keyword_set(VISUAL) then begin
    left = xcen - radius
    right = xcen + radius
    getwindow,/open, next=boxcar
    plot, fimage[400:600,2000:2500], xra=[400,600], xstyle=1, /nodata, yra=[2000,2500], ystyle=1, xtitle='X pixels (Fiber)',$
      ytitle='Y pixels (Wavelength)',title='Survey: '+strupcase(survey)+'; Boxcar flux extraction ',charsize=2, xtickformat='(A1)', ytickformat='(A1)'
    ml_tvimage, fimage[400:600,2000:2500], /axis, xra=[400,600], yra=[2000,2500], /overplot, axkeywords={charsize:2}
    index=[24,25,30,31,35,36,50]
    for i=0,n_elements(index)-1 do begin
      oplot, xcen[findgen(140)*30,index[i]], ycen[findgen(140)*30,index[i]], color=!red, psym=2
      oplot, left[*,index[i]], ycen[*,index[i]], color=!green
      oplot, right[*,index[i]], ycen[*,index[i]], color=!green
    endfor

    ;double size of extraction from 3 to 6 pixels (12 pixel span)
;    bleft = xcen-(radius*2.0)
;    bright= xcen+(radius*2.0)
;    ftemp = extract_asymbox2(fimage,bleft,bright,ycen)
;    getwindow,set=boxcar
;    oplot, bleft[*,30], ycen[*,30], color=!blue
;    oplot, bright[*,30], ycen[*,30], color=!blue
 
     ;TOP of CCD
    getwindow,/open
    plot, fimage[400:600,3000:3500], xra=[400,600], xstyle=1, /nodata, yra=[3000,3500], ystyle=1, xtitle='X pixels (Fiber)',$
      title='Survey: '+strupcase(survey)+'; Top of CCD ',charsize=2, xtickformat='(A1)', ytickformat='(A1)'
    ml_tvimage, fimage[400:600,3000:3500], /axis, xra=[400,600], yra=[3000,3500], /overplot, axkeywords={charsize:2}
    for i=0,n_elements(index)-1 do begin
    oplot, xcen[findgen(140)*30,index[i]], ycen[findgen(140)*30,index[i]], psym=2, color=!red
    oplot, left[*,index[i]], ycen[*,index[i]], color=!green
    oplot, right[*,index[i]], ycen[*,index[i]], color=!green
    endfor
    
    ;END of CCD
    getwindow,/open
    plot, fimage[3500:3800,2000:3500], xra=[400,600], xstyle=1, /nodata, yra=[2000,3500], ystyle=1, xtitle='X pixels (Fiber)',$
      title='Survey: '+strupcase(survey)+'; End of CCD ',charsize=2, xtickformat='(A1)', ytickformat='(A1)'
    ml_tvimage, fimage[3500:3800,2000:3500], /axis, xra=[3500,3800], yra=[2000,3500], /overplot, axkeywords={charsize:2}                              
    if survey eq 'manga' then index=[ntrace-25,ntrace-24,ntrace-2,ntrace-1] else if survey eq 'boss' then index=[482,483,498,499]
    for i=0,n_elements(index)-1 do begin
    oplot, xcen[findgen(140)*30,index[i]], ycen[findgen(140)*30,index[i]], psym=2, color=!red
    oplot, left[*,index[i]], ycen[*,index[i]], color=!green
    oplot, right[*,index[i]], ycen[*,index[i]], color=!green
    endfor    
    
    ;plot example extracted spectra
    getwindow,/open
    plot, fextract[*,30], psym=10, charsize=2,_extra=gang_plot_pos(2,1,0,0), tit='Survey: '+strupcase(survey)+'; Boxcar Flux Extractions from Flat Field',xtit='Wavelength'
    xyouts, 0.7, 0.8, 'Fiber 30', charsize=2, /normal
    plot, fextract[*,31], psym=10, charsize=2,_extra=gang_plot_pos(2,1,0,1)
    xyouts, 0.7, 0.4, 'Fiber 31', charsize=2, /normal
    
    ;plot approx extracted spectra for over-extraction (radius doubled)
;    getwindow,/open
;    plot, ftemp[*,30], tit='Survey: '+strupcase(survey)+'; Boxcar extraction of +-'+strtrim(radius*2.0,2)+' pixels around xcen', charsize=2,xtit='Wavelength'
;    oplot, fextract[*,30], color=!red
;    legend, ['+- '+strtrim(radius*2.0,2)+' pixels', '+- '+strtrim(radius,2)+' pixels'], /right, color=[!white,!red], charsize=2, psym=0
;    xyouts, 4000,1e5, 'Fiber 30', charsize=2
 
  endif

  if keyword_set(DEBUG) then begin
    
  endif

   return, fextract
end
;------------------------------------------------------------------------------
