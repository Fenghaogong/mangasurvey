; program mldrawcirc
;
; Draws fiber circles with fiber number overtop
; (or address if so provided)
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 11/19/2012
;   Last modified: 12/17/2012
;
; REVISION HISTORY:
;   v1: 19-Nov-2012  D. Law
;       First written.

pro mldrawcirc,xcen,ycen,corediam,fnum,dead=dead,color=color

; White is default color
if (keyword_set(color)) then color=color $
else color=255
; Note that /dead flag overrules color flag

theta=findgen(360)
x=corediam/2.*cos(theta*!PI/180.)+xcen
y=corediam/2.*sin(theta*!PI/180.)+ycen

if (keyword_set(dead)) then begin
  oplot,x,y,color=250
  xyouts,xcen,ycen,fnum,charsize=2,alignment=0.5,color=250
endif else begin
  oplot,x,y,color=color
  xyouts,xcen,ycen,fnum,charsize=2,alignment=0.5,color=color
endelse

return
end
