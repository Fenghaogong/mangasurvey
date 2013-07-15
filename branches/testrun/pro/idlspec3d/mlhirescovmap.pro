;+
; function mlhirescovmap
;
; Make a high resolution coverage map for the final data cube.
; This is for illustrative purposes, so to save space and time only make
; it for every 1000 Angstroms from 3600 to 9600.
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 11/23/2012
;   Last modified: 01/15/2013
;
; REVISION HISTORY:
;   v1: 23-Nov-2012  D. Law
;       First written.
;   v1.1: 15-Jan-2013 D. Law
;       Tweaked logging
;-
function mlhirescovmap,CovMap,xrel,yrel,fiber_status,wave

; Shared use variables
common MANGA_SHARE, pixscale, platescale, tstart

; Status starts off OK
status=0

; Sampling points, 3600-9600
sample=findgen(7)*1000.+3600.
sampleloc=value_locate(wave,sample)
nsample=(size(sample))[1]

; Define sizes
nftotal=(size(fiber_status))[1]
resln=0.1; Arcsec per pixel
Xsize=2*(max(abs(xrel*pixscale/resln)) > max(abs(yrel*pixscale/resln)))+25
Ysize=Xsize
CovMap=intarr(Xsize,Ysize,nsample)

; Set up vectors
temp=intarr(Xsize,Ysize)
corediam=2.0/resln
coordX = rebin(dindgen(Xsize), [Xsize, Ysize])
coordY = rebin(transpose(dindgen(Ysize)), [Xsize, Ysize])
radmap = replicate(0.,Xsize,Ysize)
xmin=fix(Xsize/2.-corediam)
xmax=fix(Xsize/2.+corediam)
ymin=fix(Ysize/2.-corediam)
ymax=fix(Ysize/2.+corediam)
workingregion=[n_elements(radmap[xmin:xmax]),n_elements(radmap[ymin:ymax])]

for i=0,nsample-1 do begin
  temp[*,*]=0
  for j=0,nftotal-1 do begin
     x=xrel[sampleloc[i],j]*pixscale/resln+Xsize/2.
     y=yrel[sampleloc[i],j]*pixscale/resln+Ysize/2.

    xmin=x-corediam > 0
    xmax=x+corediam < Xsize-1
    ymin=y-corediam > 0
    ymax=y+corediam < Ysize-1


    radmap[*,*]=corediam; radmap is corediam by default (big enough to ignore)
    radmap[xmin:xmax,ymin:ymax]=sqrt((coordX[xmin:xmax,ymin:ymax]-x)^2 +(coordY[xmin:xmax,ymin:ymax]-y)^2)
    ; Assign coverage
    if (fiber_status[j] eq 0) then temp[where(radmap lt corediam/2.)]+=1
  endfor
  CovMap[*,*,i]=temp
endfor

return,status
end
