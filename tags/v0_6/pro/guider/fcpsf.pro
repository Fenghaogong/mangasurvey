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
;           aperture=aperture,raobj=raobj,racen=racen,deccen=deccen,nodar=nodar,debug=debug,
; Input: 
;     psfpara:  4-element array [amp1, amp2, sig1, sig2] for the double gaussian
;              profile measured from the guider images.
;               F = amp1 * exp(-r^2/(2*sig1^2)) + amp2 * exp(-r^2/(2*sig2^2))
;              sig1 and sig2 are in arcsec. 
;              This should integrate to 1.0 over a 2-d surface.
;     ha : observing mid-point Hour angle in degree!
;     dec:  Declination of the object, not the plate
;     xfocal, yfocal: position on the plate in millimeter
;     raobj: RA of the object --- required unless nodar is set
;     racen: RA of the plate center --- required unless nodar is set
;     deccen: Dec of the plate center --- required unless nodar is set.
;     
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

function convolpix, array, kernel, scalefactor, index
   siz = size(array,/dimen)
   nx = siz[0] & ny = siz[1]
   ksiz = size(kernel,/dimen)
   nkx = fix(ksiz[0]) & nky = fix(ksiz[1])
   if (nkx mod 2) eq 0 or (nky mod 2) eq 0 then message,'The size of kernel should be odd numbered.'

   if nx lt nkx or ny lt nky then message,'The kernel needs to be smaller than the array.'
   output = fltarr(n_elements(index))
   for i=0,n_elements(index)-1 do begin
     tx = index[i] mod nx
     ty = index[i]/nx
     cutout = array[((tx-nkx/2) > 0):((tx+nkx/2) < (nx-1)),((ty-nky/2) > 0):((ty+nky/2) < (ny-1))]
     cutsiz = size(cutout,/dimen)
     if cutsiz[0] lt nkx or cutsiz[1] lt nky then begin 
;       print,'Warning: kernel reaches array boundary, pad it with zero.'
       if ty+nky/2 gt ny-1 then begin
         cutout=[[cutout],[fltarr(cutsiz[0])#fltarr(nky-cutsiz[1])]]
         cutsiz = size(cutout,/dimen)
       endif
       if ty-nky/2 lt 0 then begin
         cutout=[[fltarr(cutsiz[0])#fltarr(nky-cutsiz[1])],[cutout]]
         cutsiz = size(cutout,/dimen)
       endif
       if tx+nkx/2 gt nx-1 then begin
         tr_cutout=transpose(cutout)
         tr_cutout = [[tr_cutout],[fltarr(cutsiz[1])#fltarr(nkx-cutsiz[0])]]
         cutout = transpose(tr_cutout)
         cutsiz = size(cutout,/dimen)
       endif
       if tx-nkx/2 lt 0 then begin
         tr_cutout=transpose(cutout)
         tr_cutout = [[fltarr(cutsiz[1])#fltarr(nkx-cutsiz[0])],[tr_cutout]]
         cutout = transpose(tr_cutout)
         cutsiz = size(cutout,/dimen)
       endif
     endif
     if cutsiz[0] ne nkx or cutsiz[1] ne nky then begin
         message,'How come the size still not match?'
     endif
     output[i] = total(cutout*kernel)/scalefactor
   endfor 
   return,output
end

pro fcpsf, psfpara, ha, dec, xfocal, yfocal,cube,wavearr=wavearr,xarr=xarr,yarr=yarr,pixelscale=pixelscale,magnify=magnify,aperture=aperture,debug=debug,nodar=nodar,resetwave=resetwave,raobj=raobj,racen=racen,deccen=deccen

  time0=systime(/sec)

;  read in PSF at 5400A
;  double gaussian PSF
   amp1 = psfpara[0]
   amp2 = psfpara[1]
   sig1 = psfpara[2]
   sig2 = psfpara[3]
   
   if n_elements(wavearr) eq 0 or keyword_set(resetwave) then begin
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
   if not keyword_set(cubewidth) then cubewidth = 12.0 ; arcsec ---- size of the cube
   npixel = ceil(cubewidth/pixelscale)
   npixel = (npixel/2)*2+1 ; make npixel an odd number

;  build the cube of PSF
   psfcube = fltarr(npixel,npixel,nwave)
;---------
   x1 = (findgen(npixel)-fix(npixel/2))*pixelscale ; center of each pixel
   y1 = (findgen(npixel)-fix(npixel/2))*pixelscale ; center of each pixel

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
   npixmag = npixel*magnify
   xmag = interpolate(x1,findgen(npixel*magnify)/magnify-mean(findgen(magnify)/magnify))
   ymag = interpolate(y1,findgen(npixel*magnify)/magnify-mean(findgen(magnify)/magnify))
   dist = sqrt((xmag#(fltarr(npixmag)+1))^2+((fltarr(npixmag)+1)#ymag)^2)
   for i=0,nwave-1 do begin
     kernel=focusoff_kernel(focusarr[i]*1000.,pixel=pixelscale*platescale*1000./60/magnify)

;     kernelsize=max(size(kernel,/dimen))
;     partcube = magcube[*,npixmag/2-kernelsize/2:npixmag/2+kernelsize/2,i]
     index = (npixmag/2)*npixmag+indgen(npixmag/2)+npixmag/2
     convalue = convolpix(magcube[*,*,i],kernel,total(kernel),index)
;     tmpcube = convol(partcube,kernel,total(kernel),/center,/edge_zero)
;     concube[*,*,i] = interpol(tmpcube[npixmag/2:npixmag-1,kernelsize/2],dist[npixmag/2:npixmag-1,npixmag/2],dist)
     concube[*,*,i] = interpol(convalue,dist[index],dist)
;     concube[*,*,i] = convol(magcube[*,*,i],kernel,total(kernel),/center,/edge_zero)
   endfor
   cube2 = rebin(concube,npixel,npixel,nwave)


;  convolve with fiber size 2"
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

   outcube= cube2*0.0
   
;   partcube = cube2[*,npixel/2-kernelsize/2:npixel/2+kernelsize/2,*]
   dist = sqrt((x1#(fltarr(npixel)+1))^2+((fltarr(npixel)+1)#y1)^2)
;   outindex = npixel/2+indgen(npixel/2+1)
   for i=0,nwave-1 do begin
     index = (npixel/2)*npixel+npixel/2+indgen(npixel/2)
     convalue = convolpix(cube2[*,*,i],fiberkernel,1.0,index)
     outcube[*,*,i] = interpol(convalue,dist[index],dist)
;     psfplane[*,i] = interpol(convalue,dist[index],x1[outindex])
;     tmpcube= convol(partcube[*,*,i],fiberkernel,/center,/edge_zero)
;     outcube[*,*,i] = interpol(tmpcube[npixel/2:npixel-1,kernelsize/2],dist[npixel/2:npixel-1,npixel/2],dist)
   endfor
;stop

   shiftx = fltarr(nwave)
   shifty = fltarr(nwave)
   if NOT keyword_set(nodar) then begin
;  figure out the direction and magnitude of DAR shift
     for i=0,nwave-1 do begin
       dar = mldar(ha/15.,dec,wavearr[i],parangle,offsetx,offsety,waveref=5300,RAOBJ=raobj,RACEN=racen,DECCEN=deccen,/distort)
       shiftx[i] = offsetx/pixelscale
       shifty[i] = offsety/pixelscale
     endfor
   endif 

;  Pad the cube with zeros so that when shifting, the wrapped-around part are zero.
   leftmargin = ceil(abs((min(shiftx)) < 0))
   rightmargin = ceil(max(shiftx)) > 0
   bottommargin = ceil(abs(min(shifty) < 0))
   topmargin = ceil(max(shifty) > 0)
   padcube = fltarr(npixel+leftmargin+rightmargin,npixel+bottommargin+topmargin,nwave)
   padcube[leftmargin:leftmargin+npixel-1,bottommargin:bottommargin+npixel-1,*] = outcube
   if leftmargin gt 0 then x1=[min(x1)+(findgen(leftmargin)-leftmargin)*pixelscale,x1]
   if rightmargin gt 0 then x1=[x1,max(x1)+(findgen(rightmargin)+1)*pixelscale]
;   x1 = [min(x1)+(findgen(leftmargin)-leftmargin)*pixelscale,x1,max(x1)+(findgen(rightmargin)+1)*pixelscale]
   if bottommargin gt 0 then y1 = [min(y1)+(findgen(bottommargin)-bottommargin)*pixelscale,y1]
   if topmargin gt 0 then y1 = [y1,max(y1)+(findgen(topmargin)+1)*pixelscale]
;   y1 = [min(y1)+(findgen(bottommargin)-bottommargin)*pixelscale,y1,max(y1)+(findgen(topmargin)+1)*pixelscale]

;  shift the PSF cube according to the DAR shift
   shcube = padcube*0.0
   for i=0,nwave-1 do begin
     xint = findgen(npixel+leftmargin+rightmargin)-shiftx[i]
     yint = findgen(npixel+bottommargin+topmargin)-shifty[i]
     shcube[*,*,i] = interpolate(padcube[*,*,i],xint,yint,/grid)
   endfor

if keyword_set(debug) then begin
   save,psfcube,cube2,shcube,outcube,filename='fcpsfcube.sav'      
endif
   xarr = x1
   yarr = y1
   cube = shcube 
;  convolve with fiber aperture to get FCPSF
   print,'Time used:',systime(/sec)-time0,' sec'
end
