; program mlplotplate
;
; Plots the arrangement of fibers on the plugged plate
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 11/19/2012
;   Last modified: 01/17/2013
;
; REVISION HISTORY:
;   v1: 19-Nov-2012  D. Law
;       First written.
;   v1.2: 09-Jan-2013  D. Law
;       Added point key, color coding for fiber size
;   v1.3: 17-Jan-213  D. Law
;       Tweaked for new slitmap style and ma-numbering,
;       improved generality to arbitrary number of IFUs per plate

pro mlplotplate,slitmap,shdr

platelabel=strcompress('plate:'+string(yanny_par(shdr,'plate')),/remove_all)
mjdlabel=strcompress('mjd:'+string(yanny_par(shdr,'mjd')),/remove_all)

cenra=slitmap[0].cenra
cendec=slitmap[0].cendec

device,decomposed=0
loadct,39; Rainbow+white

platerad=1.5; degrees

theta=findgen(360)
x=platerad*cos(theta*!PI/180.)/cos(cendec*!PI/180.)+cenra
y=platerad*sin(theta*!PI/180.)+cendec

plot,x,y,xrange=[cenra+platerad/cos(cendec*!PI/180.),cenra-platerad/cos(cendec*!PI/180.)],yrange=[cendec-platerad,cendec+platerad],xtitle='RA',ytitle='DEC',charsize=2,xstyle=1,ystyle=1

; Plot IFUs
; IFUs are where names start with 'ma'.  Pick only the first entry for each.
ifuindex=where((strmid(slitmap.ifuname,0,2) eq 'ma')and(slitmap.fnum eq 1))
nifu=(size(ifuindex))[1]
ifu=slitmap[ifuindex]
oplot,ifu.ra,ifu.dec,psym=6,color=100
xyouts,ifu.ra,ifu.dec+0.1,ifu.ifuname,alignment=0.5,charsize=2,color=100

; Plot STD2 (everything with STD2 in the name)
index=where(strpos(slitmap.ifuname,'STD2') ne -1)
if (size(index))[0] eq 1 then begin
  std2=slitmap[index]
  oplot,std2.ra,std2.dec,psym=2
endif

; Plot STD3 (everything with STD3 in the name)
index=where(strpos(slitmap.ifuname,'STD3') ne -1)
if (size(index))[0] eq 1 then begin
  std3=slitmap[index]
  oplot,std3.ra,std3.dec,psym=2,color=150
endif

; Plot STD5 (everything with STD5 in the name)
index=where(strpos(slitmap.ifuname,'STD5') ne -1)
if (size(index))[0] eq 1 then begin
  std5=slitmap[index]
  oplot,std5.ra,std5.dec,psym=2,color=210
endif

; Plot SKY2 (everything with SKY2 in the name)
index=where(strpos(slitmap.ifuname,'SKY2') ne -1)
if (size(index))[0] eq 1 then begin
  sky2=slitmap[index]
  oplot,sky2.ra,sky2.dec,psym=1
endif

; Plot SKY3 (everything with SKY3 in the name)
index=where(strpos(slitmap.ifuname,'SKY3') ne -1)
if (size(index))[0] eq 1 then begin
  sky3=slitmap[index]
  oplot,sky3.ra,sky3.dec,psym=1,color=150
endif

; Plot SKY5 (everything with SKY5 in the name)
index=where(strpos(slitmap.ifuname,'SKY5') ne -1)
if (size(index))[0] eq 1 then begin
  sky5=slitmap[index]
  oplot,sky5.ra,sky5.dec,psym=1,color=210
endif

; Show a key
xyouts,(cenra+platerad/cos(cendec*!PI/180.))-0.1,cendec+platerad-0.2,'* = STD',charsize=2
xyouts,(cenra+platerad/cos(cendec*!PI/180.))-0.1,cendec+platerad-0.4,'+ = SKY',charsize=2
xyouts,(cenra+platerad/cos(cendec*!PI/180.))-0.1,cendec-platerad+0.1,'2"',charsize=2
xyouts,(cenra+platerad/cos(cendec*!PI/180.))-0.1,cendec-platerad+0.3,'3"',charsize=2,color=150
xyouts,(cenra+platerad/cos(cendec*!PI/180.))-0.1,cendec-platerad+0.5,'5"',charsize=2,color=210
xyouts,(cenra-platerad/cos(cendec*!PI/180.))+0.7,cendec+platerad-0.2,platelabel,charsize=2
xyouts,(cenra-platerad/cos(cendec*!PI/180.))+0.6,cendec+platerad-0.4,mjdlabel,charsize=2
;http://www.astro.virginia.edu/class/oconnell/astr511/idl_5.1_html/idl22c.htm

return
end
