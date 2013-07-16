;+
; function mlreaddrillfiles
;
; Reads information from the drill files (plDrillPos and the CMM file)
; about where the holes in a plate were actually drilled vs where
; they were meant to be drilled.  Unfortunately these files aren't in
; useful order, so much sorting needs to be done.
;
; Note that the code requires the CMM data from the UWash shop first
; be converted to Yanny-style format using mlconvcmm.pro
;
; Uses this information to calculate x,y offsets *in microns* in the focal plane
; of actual from fiducial holes for a set of input RA/DEC coordinates.
;
; Note that if a hole is not found for a given RA/DEC the
; corresponding xdrill, ydrill will simply be zero.  This is because
; the CMM cannot measure holes near the +/- y edge of the plate and
; there will therefore always be some holes with no CMM data.
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 11/27/2012
;   Last modified: 01/15/2013
;
; REVISION HISTORY:
;   v1: 27-Nov-2012  D. Law
;       First written.
;   v1.1: 10-Dec-2012  D. Law
;       Modified to read in CMM data from Yanny-style file.
;   v1.1: 15-Jan-2013 D. Law
;       Tweaked logging
;-
function mlreaddrillfiles, drillfile, drillcmmfile,objra,objdec,xdrill,ydrill,xresid,yresid,rresid,nominalx,nominaly

; Status starts off OK
status=0

; Drill offsets are zero by default unless an appropriate hole and
; associated offset is found to make them otherwise.
xdrill=0.
ydrill=0.

splog,strcompress('Drillfile: '+drillfile)
; Check whether drillfile exists, return error code -60 if not
junk=findfile(drillfile,count=ct)
if (ct EQ 0) then return,-60L

splog,strcompress('DrillCMMfile: '+drillcmmfile)
; Check whether drillcmmfile exists, return error code -61 if not
junk=findfile(drillcmmfile,count=ct)
if (ct EQ 0) then return,-61L

; Read drillfile structure
drillpos=yanny_readone(drillfile,'DRILLPOS',hdr=hdr)

; Read drillcmmfile.
drillcmm=yanny_readone(drillcmmfile,'CMMPOS',hdr=cmmhdr)
xglobal=double(yanny_par(cmmhdr, 'FitOffsetX'))
yglobal=double(yanny_par(cmmhdr, 'FitOffsetY'))
scaleglobal=double(yanny_par(cmmhdr, 'FitScale'))
rotglobal=double(yanny_par(cmmhdr, 'FitRotAngle'))

nomx=drillcmm.nomx
nomy=drillcmm.nomy
measx=drillcmm.measx
measy=drillcmm.measy

; Transform measured x,y locations by the global rotation, offset
; and scale parameters.  Based on matrix code provided by Conor
; Sayres.  I think the sign on coefficients is correct.
; NEED TO CHECK WITH CONOR.
c0=xglobal
c1=yglobal
c3=scaleglobal*cos(rotglobal*!PI/180.)
c2=-scaleglobal*sin(rotglobal*!PI/180.)
;c3=sqrt(scaleglobal*scaleglobal/(1+(tan(rotglobal*!PI/180.))^2))
;c2=sqrt(scaleglobal*scaleglobal-c3*c3)

measy=(measy-c1+c2/c3*measx-c2/c3*c0)/(c3+c2*c2/c3)
measx=(measx-c0-c2*measy)/c3

  ; Figure out which (if any) of the holes is a good match
  temp=sqrt(((drillpos.ra-objra)*cos(objdec*!PI/180.))^2 + ((drillpos.dec-objdec)^2))
  ; If offset of the hole from the coordinate is < 1 arcsec, it's a match
  match=where(temp*3600. lt 1.)

  ; If exactly 1 matching hole was found in the drillfile, process it
  if (size(match))[1] eq 1 then begin
    xfocal=drillpos[match].xFocal
    yfocal=drillpos[match].yFocal
    xflat=drillpos[match].xFlat
    yflat=drillpos[match].yFlat
;print,'CMM info:',xfocal,yfocal,xflat,yflat
;stop
    ; Find if there is a corresponding hole in the CMM data by
    ; matching xflat to nomx and yflat to nomy
    temp2=sqrt((xflat-nomx)^2+(yflat-nomy)^2)
    match2=where(temp2 lt 0.01)

    ; If exactly 1 matching hole was found in the drillCMMfile, process it
    if (size(match2))[1] eq 1 then begin
      xdrill=(measx[match2]-nomx[match2])*1000.
      ydrill=(measy[match2]-nomy[match2])*1000.

      xresid=drillcmm[match2].qpresidx
      yresid=drillcmm[match2].qpresidy
      rresid=drillcmm[match2].qpresidr
      nominalx=nomx[match2]
      nominaly=nomy[match2]
;      print,nomx[match2],nomy[match2],drillcmm[match2].qpresidx,drillcmm[match2].qpresidy,drillcmm[match2].qpresidr

    endif
  endif

; Strictly we've computed the offset in the *flat* plate plane, not in
; the curved focal plane.  However, the scales of the coordinates are
; sufficiently similar that small shifts in one translate quite
; accurately to small shifts in the other.

; If any of the offsets are bigger than 60 microns, there's a problem!
if (sqrt(xdrill*xdrill+ydrill*ydrill) gt 60.) then begin
  splog,strcompress('Problem in drill file '+drillcmmfile)
  splog,'One or more holes over 60 microns from intended position!'
  return,-62L
endif

return,status
end

