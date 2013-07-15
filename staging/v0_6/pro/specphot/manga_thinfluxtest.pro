; function manga_thinfluxtest
;
; Simulate flux losses from 2''/3'' calibration fibers
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 07/24/2012
;   Last modified: 10/17/2012
;
; Modification history:
;   v1: Implemented with commentary.
;   v2: Updated 10/17/2012 to work with revised library filenames
;
forward_function mldar

pro manga_thinfluxtest
radpdeg=3.14159/180.
PixScale=0.05; arcsec/pixel for the test image

BoxSize=6./PixScale; Boxsize for test image
Image=dblarr(BoxSize,BoxSize)
fluxtot=1000.D; Total flux of standard star (arb units)
Image[BoxSize/2.,BoxSize/2.]=fluxtot

seeing0=1.3; Seeing, in arcseconds

nrand=500; Number of random draws

; Loop over wavelength
lam=findgen(65)*100.+3600.
valA=fltarr(65,nrand)
valAmean=fltarr(65)
valAsig=fltarr(65)
valB=fltarr(65,nrand)
valBmean=fltarr(65)
valBsig=fltarr(65)

; Loop over random position error of 0.15'' RMS
randpar1=randomn(systime_seed,nrand)*0.15; Offset error
randpar2=randomu((systime_seed+1),nrand)*360.; Direction of offset
randparX=randpar1*cos(randpar2*3.14159/180.); Offset along X
randparY=randpar1*sin(randpar2*3.14159/180.); Offset along Y

OpenW,UnitW,'test.dat',/get_lun;christy_am1p2_1p5see.dat',/Get_Lun
printf,unitw,'# lambda		valAmean	valAsig		valBmean	valBsig		ratio	rsig'

; Map of radii from image center
radmap=dblarr(BoxSize,BoxSize)

;for i=64,64 do begin
for i=0,64 do begin
  print,i,lam[i]

  seeing=seeing0*((lam[i]/5500.)^(-0.2))
  seeing=seeing/PixScale
  twoseeing=2*seeing

  BlurredImage1=mlfilterimg(Image,fwhm_gaussian=seeing,/all_pixels)*9./13.
  BlurredImage2=mlfilterimg(Image,fwhm_gaussian=twoseeing,/all_pixels)*4./13.
  SimImage=BlurredImage1+BlurredImage2
 ; writefits,'SimImage.fits',SimImage

  ; Define fiber center
  ; Fiducial
  xcen0=BoxSize/2.
  ycen0=BoxSize/2.

  ; Move for DAR
  darHA=2.09; Gives airmass=1.2
  darTIME=35./60.
  darDEC=15.0
  darWAVE=lam[i]
  shift=mldar(darHA,darDEC,darWAVE,parangA,shiftX,shiftY,waveref=5000.)

  radmap[*,*]=0.D

  for q=0,nrand-1 do begin
  ; Image center shift is sum of DAR and random component
  shiftX=-shift*sin(parangA*radpdeg)+randparX[q]
  shiftY=shift*cos(parangA*radpdeg)+randparY[q]
  xcen=xcen0-shiftX/PixScale
  ycen=ycen0-shiftY/PixScale

  ; Define 2d arrays describing the fiber
  ; fibermap is -1 everywhere except where a fiber covers, where
  ; it is equal to the fiber number
  fiberradA=1.0;arcsec
  fiberradA=fiberradA/PixScale;pixels
  fiberradB=1.5;arcsec
  fiberradB=fiberradB/PixScale;pixels

  for j=0,BoxSize-1 do begin
     for k=0,BoxSize-1 do begin
        radmap[j,k]=sqrt((j-xcen)*(j-xcen)+(k-ycen)*(k-ycen))
     endfor
  endfor

; Amount of flux down the fibers is equal to sum of values at all
; radii less than fiber radius
  valA[i,q]=total(SimImage[where(radmap lt fiberradA)])/fluxtot
  valB[i,q]=total(SimImage[where(radmap lt fiberradB)])/fluxtot
  endfor

  temp1=0.
  temp2=0.
  temp3=0.
  temp3=mlmeanclip(valA[i,*],temp1,temp2,CLIPSIG=3.0)
  valAmean[i]=temp1
  valAsig[i]=temp2
  temp3=mlmeanclip(valB[i,*],temp1,temp2,CLIPSIG=3.0)
  valBmean[i]=temp1
  valBsig[i]=temp2

  ratio=valAmean[i]/valBmean[i]
; propagate uncertainty to the ratio
  runc=sqrt(ratio*ratio*(valAsig[i]*valAsig[i]/valAmean[i]/valAmean[i]+valBsig[i]*valBsig[i]/valBmean[i]/valBmean[i]))

  printf,unitw,lam[i],valAmean[i],valAsig[i],valBmean[i],valBsig[i],ratio,runc,FORMAT='(7(7x,d12.6))'
endfor

close,unitw
free_lun,unitw

return
end

