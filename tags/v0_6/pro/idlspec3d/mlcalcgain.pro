;+
;
; Quick program to calculate gain from a flatfield flux image,
; need to improve this routine and integrate into the main pipeline,
; have it put the typical gain in the image header?
;
;-

pro mlcalcgain

blue='spFlat1flux-b1-00153995.fits'
red='spFlat1flux-r1-00153995.fits'

bim=readfits(blue)
rim=readfits(red)

nb=(size(bim))[1]
nr=(size(rim))[1]

xb=indgen(nb)
xr=indgen(nr)

status=mlreadslm('6653','56298',slitmap)
goodfib=where(slitmap.fsize eq 2,ngood)

bim=bim[*,goodfib]
rim=rim[*,goodfib]

bgain=fltarr(ngood)
rgain=fltarr(ngood)
for i=0,ngood-1 do begin
  bvec=bim[*,i]
  rvec=rim[*,i]

  nord=3
  everyn=20

  ssetb=bspline_iterfit(xb,double(bvec),nord=nord,everyn=everyn,bkpt=0,maxrej=0,yfit=yfitb)
  ssetr=bspline_iterfit(xr,double(rvec),nord=nord,everyn=everyn,bkpt=0,maxrej=0,yfit=yfitr)

  ratiob=(bvec-yfitb)/sqrt(yfitb)
  ratior=(rvec-yfitr)/sqrt(yfitr)
  b=mlmeanclip(ratiob[1000:3000],b1,b2)
  r=mlmeanclip(ratior[1000:3000],r1,r2)
  bgain[i]=b2*b2
  rgain[i]=r2*r2
endfor

bavggain=mlmeanclip(bgain,temp1,temp2)
ravggain=mlmeanclip(rgain,temp1,temp2)
print,'Average blue gain:',1./bavggain,' e-/ADU'
print,'Average red gain:',1./ravggain,' e-/ADU'
print,'Multiply fluxes by this number to get e- units'

return
end
