;+
; function mlfibersn
;
; Make a high resolution S/N ratio map for each fiber that went into
; the stacked cube.  Do it for blue and red wavelengths
; (g band and i-band)
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 03/07/2013
;   Last modified: 03/07/2013
;
; REVISION HISTORY:
;   v1: 07-Mar-2013  D. Law
;       First written based on BOSS S/N calculations in
;       quickextract.pro
;-
function mlfibersn,xrel,yrel,fiber_status,wave,flux,ivar,SNgMap,SNrMap,SNiMap,SNTable,flipstat=flipstat

; Shared use variables
common MANGA_SHARE, pixscale, platescale, tstart

; Status starts off OK
status=0

nfiber=(size(fiber_status))[1]

; g-band wavelength coverage (Blue ccd)
wrange_gband = [4000,5500]
gselect=where((wave gt wrange_gband[0])and(wave le wrange_gband[1]),ngsel)
gmid=gselect[(size(gselect))[1]/2]; Midpoint of the band

; r-band wavelength coverage (Blue/Red ccd)
wrange_rband = [5611,6969]
rselect=where((wave gt wrange_rband[0])and(wave le wrange_rband[1]),nrsel)
rmid=rselect[(size(rselect))[1]/2]; Midpoint of the band

; i-band wavelength coverage (Red ccd)
wrange_iband = [6910,8500]
iselect=where((wave gt wrange_iband[0])and(wave le wrange_iband[1]),nisel)
imid=iselect[(size(iselect))[1]/2]; Midpoint of the band

; Calculate S/N values in these ranges
meangsn=fltarr(nfiber)
meanrsn=fltarr(nfiber)
meanisn=fltarr(nfiber)
for i=0,nfiber-1 do begin
  snvec=djs_median(flux[gselect,i]*sqrt(ivar[gselect,i]),width=ngsel,boundary='reflect')
  meangsn[i]=djs_mean(snvec)

  snvec=djs_median(flux[rselect,i]*sqrt(ivar[rselect,i]),width=nrsel,boundary='reflect')
  meanrsn[i]=djs_mean(snvec)

  snvec=djs_median(flux[iselect,i]*sqrt(ivar[iselect,i]),width=nisel,boundary='reflect')
  meanisn[i]=djs_mean(snvec)
endfor

; Put values into the table
SNTable=fltarr(nfiber,9)
SNTable[*,0]=-xrel[gmid,*]*pixscale/3600.D; g X posn in degrees (flip b/c E left)
SNTable[*,1]=yrel[gmid,*]*pixscale/3600.D; g Y posn in degrees
SNTable[*,2]=meangsn[*]*meangsn[*]; g S/N ratio^2

SNTable[*,3]=-xrel[rmid,*]*pixscale/3600.D; r X posn in degrees (flip b/c E left)
SNTable[*,4]=yrel[rmid,*]*pixscale/3600.D; r Y posn in degrees
SNTable[*,5]=meanrsn[*]*meanrsn[*]; r S/N ratio^2

SNTable[*,6]=-xrel[imid,*]*pixscale/3600.D; i X posn in degrees (flip b/c E left)
SNTable[*,7]=yrel[imid,*]*pixscale/3600.D; i Y posn in degrees
SNTable[*,8]=meanisn[*]*meanisn[*]; i S/N ratio^2

; Define sizes
resln=0.1; Arcsec per pixel for this map
Xsize=2*(max(abs(xrel*pixscale/resln)) > max(abs(yrel*pixscale/resln)))+25
Ysize=Xsize
SNgMap=fltarr(Xsize,Ysize)
SNrMap=fltarr(Xsize,Ysize)
SNiMap=fltarr(Xsize,Ysize)

; Set up vectors
corediam=2.0/resln
coordX = rebin(dindgen(Xsize), [Xsize, Ysize])
coordY = rebin(transpose(dindgen(Ysize)), [Xsize, Ysize])
radmap = replicate(0.,Xsize,Ysize)
xmin=fix(Xsize/2.-corediam)
xmax=fix(Xsize/2.+corediam)
ymin=fix(Ysize/2.-corediam)
ymax=fix(Ysize/2.+corediam)
workingregion=[n_elements(radmap[xmin:xmax]),n_elements(radmap[ymin:ymax])]

; Allow what is 'good' fiber status to flip with a keyword to
; accomodate different conventions for 0/1 good/bad
goodstatus=1
if (keyword_set(flipstat)) then goodstatus=0

; Make the map for g wavelengths
for j=0,nfiber-1 do begin
  x=xrel[gmid,j]*pixscale/resln+Xsize/2.
  y=yrel[gmid,j]*pixscale/resln+Ysize/2.

  xmin=x-corediam > 0
  xmax=x+corediam < Xsize-1
  ymin=y-corediam > 0
  ymax=y+corediam < Ysize-1

  radmap[*,*]=corediam; radmap is corediam by default (big enough to ignore)
  radmap[xmin:xmax,ymin:ymax]=sqrt((coordX[xmin:xmax,ymin:ymax]-x)^2 +(coordY[xmin:xmax,ymin:ymax]-y)^2)

  ; Assign values
  if (fiber_status[j] eq goodstatus) then SNgMap[where(radmap lt corediam/2.)]+=meangsn[j]*meangsn[j]
endfor

; Make the map for r wavelengths
for j=0,nfiber-1 do begin
  x=xrel[rmid,j]*pixscale/resln+Xsize/2.
  y=yrel[rmid,j]*pixscale/resln+Ysize/2.

  xmin=x-corediam > 0
  xmax=x+corediam < Xsize-1
  ymin=y-corediam > 0
  ymax=y+corediam < Ysize-1

  radmap[*,*]=corediam; radmap is corediam by default (big enough to ignore)
  radmap[xmin:xmax,ymin:ymax]=sqrt((coordX[xmin:xmax,ymin:ymax]-x)^2 +(coordY[xmin:xmax,ymin:ymax]-y)^2)

  ; Assign values
  if (fiber_status[j] eq goodstatus) then SNrMap[where(radmap lt corediam/2.)]+=meanrsn[j]*meanrsn[j]
endfor

; Make the map for i wavelengths
for j=0,nfiber-1 do begin
  x=xrel[imid,j]*pixscale/resln+Xsize/2.
  y=yrel[imid,j]*pixscale/resln+Ysize/2.

  xmin=x-corediam > 0
  xmax=x+corediam < Xsize-1
  ymin=y-corediam > 0
  ymax=y+corediam < Ysize-1

  radmap[*,*]=corediam; radmap is corediam by default (big enough to ignore)
  radmap[xmin:xmax,ymin:ymax]=sqrt((coordX[xmin:xmax,ymin:ymax]-x)^2 +(coordY[xmin:xmax,ymin:ymax]-y)^2)

  ; Assign values
  if (fiber_status[j] eq goodstatus) then SNiMap[where(radmap lt corediam/2.)]+=meanisn[j]*meanisn[j]
endfor

return, status
end
