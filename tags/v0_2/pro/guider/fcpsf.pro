;+
; Name: fcpsf
; Purpose: construct the FCPSF cube (x,y, lambda) for a given PSF profile, 
;          observing ha, dec, and position on plate (xfocal,yfocal).
;          The calculation takes into account of the following:
;          (a) the seeing variation with wavelength due to atmosphere
;          (b) PSF variation due to focal plane shape changing with wavelength.
;          (c) Differential Atmosphere Refraction
;          (d) integration of flux within fiber aperture
; Usage:
;    fcpsf, psfpara, ha, dec, xfocal, yfocal,cube, wavearr=wavearr,
;           xarr=xarr,yarr=yarr,pixelscale=pixelscale,magnify=magnify,
;           aperture=aperture,debug=debug
; Input: 
;     psfpara:  4-element array [amp1, amp2, sig1, sig2] for the double gaussian
;              profile measured from the guider images.
;               F = amp1 * exp(-r^2/(2*sig1^2)) + amp2 * exp(-r^2/(2*sig2^2))
;              sig1 and sig2 are in arcsec. 
;              This should integrate to 1.0 over a 2-d surface.
;     ha : observing mid-point Hour angle in degree!
;     dec:  Declination of the object, not the plate
;     xfocal, yfocal: position on the plate in millimeter
; Optional Input: 
;     wavearr: a 1-D array giving the wavelength grid for the cube in Angstroms.
;     pixelscale: the pixelscale in arcsec. The default is 0.1". 
;                 Smaller pixelscale makes the code very slow. 
;                 Instead, one should use the magnify keyword. 
;     cubewidth:  the cubewidth in arcsec. the default is 8".
;     magnify:  the magnification for the convolution with the 
;               focus offset kernel. Default magnify=1
;               For pixelscale=0.1", magnify=1, the result is good to 2%
;                              with magnify=2, the result is better than 1%
;                                  magnify=3 or 5 makes it nearly perfect.
;                 
;     aperture: the fiber diameter in arcsec. The default is 2".
; Keyword: 
;     nodar : the default is to apply an DAR offset. This can be turned off 
;             by specifying nodar=1.
;  
; Output:  
;     cube:  the fiber-convolved PSF cube (x, y, wavelength)
; Optional output:
;     xarr:  the x-coordinate for the cube, 0 is the center of the star 
;            at 5300A. +x = +RA
;     yarr:  the y-coordinate for the cube, 0 is the center of the star 
;            at 5300A. +y = +Dec
;
;   Example:
;IDL>  fcpsf, param, ha, dec, xfocal,yfocal, cube, wavearr=wavearr, xarr=xarr,yarr=yarr, magnify=2
;IDL>  dRA=-0.417 ; on-sky fiber offset relative to star center at 5300A, East is +x
;IDL>  dDec=0.721 ; fiber offset relative to star center at 5300A, North is +y
;IDL>  xpos= interpol(findgen(n_elements(xarr)), xarr, dRA)
;IDL>  ypos= interpol(findgen(n_elements(yarr)), yarr, dDec)
;IDL>  thrupt = interpolate(cube, xpos, ypos, wavearr, /grid)
;IDL>  plot, wavearr, thrupt
;   
; Revision history:
;    Written on Feb 5, 2013 by Renbin Yan 
;          (Univ. of Kentucky, yanrenbin@gmail.com)
;     
;-
pro fcpsf, psfpara, ha, dec, xfocal, yfocal,cube,wavearr=wavearr,xarr=xarr,yarr=yarr,pixelscale=pixelscale,magnify=magnify,aperture=aperture,debug=debug,nodar=nodar

;  read in PSF at 5400A
;  double gaussian PSF
   amp1 = psfpara[0]
   amp2 = psfpara[1]
   sig1 = psfpara[2]
   sig2 = psfpara[3]
   
   if n_elements(wavearr) eq 0 then begin
      maxwave = 10600. & minwave=3500.
      dwave = 50.
      nwave = (maxwave-minwave)/dwave
      wavearr = minwave+findgen(nwave)*dwave
   endif else begin
      nwave = n_elements(wavearr)
   endelse
;  adjust the width of the PSF for different wavelengths
   wave0 = 5400.
   factor = (wavearr/wave0)^(-0.2)
   factorsq = factor*factor

   sig1arr = sig1*factor
   sig2arr = sig2*factor
   sig1sq = sig1arr*sig1arr
   sig2sq = sig2arr*sig2arr

   amp1arr = amp1/factorsq ; to keep amp*sig^2 fixed.
   amp2arr = amp2/factorsq

   platescale = 3.62730 ; (mm/arcmin) as given in Gunn et al. 2006

   if not keyword_set(pixelscale) then pixelscale = 0.1 ; arcsec/pixel --- pixel scale for the constructed cube
   if not keyword_set(cubewidth) then cubewidth = 8.0 ; arcsec ---- size of the cube
   npixel = ceil(cubewidth/pixelscale)

;  build the cube of PSF
   psfcube = fltarr(npixel,npixel,nwave)
;---------
   x1 = (findgen(npixel)+0.5-(npixel/2.))*pixelscale ; center of each pixel
   y1 = (findgen(npixel)+0.5-(npixel/2.))*pixelscale ; center of each pixel

   xb = [x1-0.5*pixelscale,x1[npixel-1]+0.5*pixelscale]; boundaries of pixels
   yb = [y1-0.5*pixelscale,y1[npixel-1]+0.5*pixelscale]

   for i=0,nwave-1 do begin
      tx1 = xb/sqrt(2.*sig1sq[i])
      ty1 = yb/sqrt(2.*sig1sq[i])
      ex1 = sqrt(!pi)/2.*sqrt(2*sig1sq[i])*erf(tx1)
      ey1 = sqrt(!pi)/2.*sqrt(2*sig1sq[i])*erf(ty1)
      exp1_x = ex1[1:npixel]-ex1[0:npixel-1]
      exp1_y = ey1[1:npixel]-ey1[0:npixel-1]

      tx2 = xb/sqrt(2.*sig2sq[i])
      ty2 = yb/sqrt(2.*sig2sq[i])
      ex2 = sqrt(!pi)/2.*sqrt(2*sig2sq[i])*erf(tx2)
      ey2 = sqrt(!pi)/2.*sqrt(2*sig2sq[i])*erf(ty2)
      exp2_x = ex2[1:npixel]-ex2[0:npixel-1]
      exp2_y = ey2[1:npixel]-ey2[0:npixel-1]

      psfcube[*,*,i] = amp1arr[i]*(exp1_x # exp1_y) + amp2arr[i]*(exp2_x # exp2_y)
   endfor

;  figure out the focus offset as a function of wavelength given xfocal, yfocal
   focusoffset,xfocal,yfocal,wavearr,focusarr
   ; focusarr returned in unit of millimeter 

;  built a kernel cube according to the focusoffset at each wavelength.
;  convolve the PSF at each wave with the kernel 
   if NOT keyword_set(magnify) then magnify = 1
;   magcube = rebin(psfcube,npixel*magnify,npixel*magnify,nwave)
   if magnify gt 1 then begin
       magcube = fltarr(npixel*magnify,npixel*magnify,nwave)
       for i=0,nwave-1 do begin
         magcube[*,*,i] = interpolate(psfcube[*,*,i],findgen(npixel*magnify)/magnify-mean(findgen(magnify)/magnify),findgen(npixel*magnify)/magnify-mean(findgen(magnify)/magnify),/grid,/cubic)
       endfor
   endif else magcube = psfcube
   concube = magcube*0.0
   for i=0,nwave-1 do begin
     kernel=focusoff_kernel(focusarr[i]*1000.,pixel=pixelscale*platescale*1000./60/magnify)
     concube[*,*,i] = convol(magcube[*,*,i],kernel,total(kernel),/center,/edge_zero)
   endfor
   cube2 = rebin(concube,npixel,npixel,nwave)

   if NOT keyword_set(nodar) then begin
;  figure out the direction and magnitude of DAR shift
     dar = mldar(ha/15.,dec,wavearr,parangle,offsetx,offsety,waveref=5300)
     shiftx = offsetx/pixelscale
     shifty = offsety/pixelscale
   endif else begin
     shiftx = fltarr(nwave)
     shifty = fltarr(nwave)
   endelse

;  Pad the cube with zeros so that when shifting, the wrapped-around part are zero.
   leftmargin = ceil(abs((min(shiftx)) < 0))
   rightmargin = ceil(max(shiftx)) > 0
   bottommargin = ceil(abs(min(shifty) < 0))
   topmargin = ceil(max(shifty) > 0)
   padcube = fltarr(npixel+leftmargin+rightmargin,npixel+bottommargin+topmargin,nwave)
   padcube[leftmargin:leftmargin+npixel-1,bottommargin:bottommargin+npixel-1,*] = cube2
   if leftmargin gt 0 then x1=[min(x1)+(findgen(leftmargin)-leftmargin)*pixelscale,x1]
   if rightmargin gt 0 then x1=[x1,max(x1)+(findgen(rightmargin)+1)*pixelscale]
;   x1 = [min(x1)+(findgen(leftmargin)-leftmargin)*pixelscale,x1,max(x1)+(findgen(rightmargin)+1)*pixelscale]
   if bottommargin gt 0 then y1 = [min(y1)+(findgen(bottommargin)-bottommargin)*pixelscale,y1]
   if topmargin gt 0 then y1 = [y1,max(y1)+(findgen(topmargin)+1)*pixelscale]
;   y1 = [min(y1)+(findgen(bottommargin)-bottommargin)*pixelscale,y1,max(y1)+(findgen(topmargin)+1)*pixelscale]

;  shift the PSF cube according to the DAR shift
;  and then convolve with fiber size 2"
   if NOT keyword_set(aperture) then aperture = 2 ;(arcsec) aperture diameter
   kernelsize = ceil(aperture/pixelscale)+2
   kernelsize = 2*(ceil(kernelsize/2.0))+1
   aperker = fltarr(kernelsize,kernelsize)

   magnification = 5
   msiz = fix(kernelsize*magnification)
   mkernel = fltarr(msiz,msiz)
   xm = findgen(msiz)-msiz/2
   x2 = xm#(fltarr(msiz)+1)
   y2 = (fltarr(msiz)+1)#xm
   rr = sqrt(x2*x2+y2*y2)*pixelscale/magnification
   t = where(rr le aperture/2.)
   mkernel[t] = 1.0
   fiberkernel = rebin(mkernel,kernelsize,kernelsize)

   shcube = padcube*0.0
   outcube= padcube*0.0
   for i=0,nwave-1 do begin
     xint = findgen(npixel+leftmargin+rightmargin)-shiftx[i]
     yint = findgen(npixel+bottommargin+topmargin)-shifty[i]
     shcube[*,*,i] = interpolate(padcube[*,*,i],xint,yint,/grid)
     outcube[*,*,i] = convol(shcube[*,*,i],fiberkernel,total(fiberkernel),/center,/edge_zero)
   endfor

if keyword_set(debug) then begin
   save,psfcube,cube2,shcube,outcube,filename='fcpsfcube.sav'      
endif
   xarr = x1
   yarr = y1
   cube = outcube
;  convolve with fiber aperture to get FCPSF
end
