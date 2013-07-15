; program mlplotbm
;
; Plots a bundle map in units of arcsec
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 11/19/2012
;   Last modified: 01/16/2013
;
; REVISION HISTORY:
;   v1: 19-Nov-2012  D. Law
;       First written.
;   v1.1: 16-Jan-2013 D. Law
;       Revised for new metrology database format
pro mlplotbm,bmap,bhdr

device,decomposed=0
loadct,39; Rainbow+white

fiberdiam=2.5;arcsec
corediam=2.0;arcsec
scale=60.;microns per arcsec
scalemm=scale/1000.;mm per arcsec

nfiber=long(yanny_par(bhdr,'nfiber'))
xwin=sqrt(nfiber/!PI)*fiberdiam
ywin=sqrt(nfiber/!PI)*fiberdiam
name=string(yanny_par(bhdr, 'ifuname'))

x=fltarr(2)
y=x
plot,x,y,xrange=[xwin,-xwin],yrange=[-ywin,ywin],xtitle='Delta RA (arcsec)',ytitle='Delta Dec (arcsec)',charsize=2
for i=0,nfiber-1 do begin
  xcen=bmap[i].xpmm/scalemm
  ycen=bmap[i].ypmm/scalemm
  fnum=strcompress(string(bmap[i].ise),/remove_all)
  gbu=bmap[i].gbu
  label=fnum

  if gbu eq 1 then mldrawcirc,xcen,ycen,corediam,label
  ;else mldrawcirc,xcen,ycen,corediam,label,/dead
endfor

; Label N and E
xyouts,xwin+1,0, 'E',charsize=2,alignment=0.5
xyouts,0,ywin+1, 'N',charsize=2,alignment=0.5

; Show the pin
if nfiber lt 60 then labdist=1.0 $
else labdist=2.0
xyouts,-(xwin+labdist),0, 'pin',charsize=2,alignment=0.5

; Label bundle names
xyouts,xwin,ywin,name,charsize=2,alignment=0.5

return
end
