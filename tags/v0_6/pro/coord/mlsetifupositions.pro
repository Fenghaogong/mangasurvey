;+
; function mlsetifupositions
;
; Figure out what the relative position of every fiber in an IFU was.
;
; Note that this code assumes that objra, objdec are single numbers and
; xinbundle, and yinbundle are vectors.
; If supplying just single values for x/yinbundle make sure to put them in brackets
; (e.g., [0.],[0.])
;
; Incorporates position within bundle, error of drilled holes from
; fiducial locations, base shift at guide wavelength due to spatial
; DAR arising from observation at a specific time, offset due to
; dithering, and shift due to DAR.  Incorporates effects of changing
; scale and rotation from guiding optimized over entire plate.
;
; Note that xpos, ypos returned have dimension [nwave, nfiber]
; Note that objra, objdec should be in decimal degrees.
; Input wavelength vector 'wave' must be in Angstroms
; Input argument 'obsparam' is from the .exp file and contains
; information about the exposure.
;
; Output is in units of arcsec in the -xfocal,yfocal coordinate frame
; (i.e., North up, East left)
; (+x = -RA, +y = +DEC)
; Assumes that xinbundle, yinbundle are
; (+x = +RA, +y = +DEC)
;
; Returns 0 if everything is ok, returns an error code if there was a problem.
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 11/22/2012
;   Last modified: 05/29/2013
;
; REVISION HISTORY:
;   v1: 22-Nov-2012  D. Law
;       First written.
;   v1.1: 27-Nov-2012  D. Law
;       Modified to be useful for single fibers, return in units of arcsec
;   v1.2: 30-Dec-2012  D. Law
;       Adding dither positions 4-9, tweak signs to be correct
;   v1.3: 21-Jan-2013 D. Law
;       Tweaked logging, changed units on bundle metrology
;   v1.4: 29-May-2013 D. Law
;       Modified calculation of base offsets due to scale changes.
;       Now better accounts for anticipated changes in plate scale,
;       rotation, and shifts to track a given field.
;-
function mlsetifupositions,obsparam,objra,objdec,xinbundle,yinbundle,wave,xpos,ypos

; Shared use variables
common MANGA_SHARE, pixscale, platescale, tstart

; Status starts off OK
status=0

; platescale in mm/arcsec
platescalemm=platescale/1000.

nwave=(size(wave))[1]
nfiber=(size(xinbundle))[1]

; The zero-point for everything is the center of the bundle, if
; observed at the drill HA, and the guide wavelength.

;;;;;;;;;
; Hack May 9 2013 to force x,y positions in bundle to be wrong
; by a specified rotation and x,y shift

;randno=randomn(systime_seed,3)
;print,randno
;rotang=3.*randno[0];degrees
;xbad=0.1*platescalemm*randno[1];0.1'' x error
;ybad=0.2*platescalemm*randno[2];0.1'' y error
;xnew=xinbundle*cos(rotang*!PI/180.)+yinbundle*sin(rotang*!PI/180.)+xbad
;ynew=-xinbundle*sin(rotang*!PI/180.)+yinbundle*cos(rotang*!PI/180.)+ybad
;xinbundle=xnew
;yinbundle=ynew
;splog,'Deliberately screwing up the fiber positions!'
; Output offsets in arcsec and degrees
;splog,xbad/platescalemm,ybad/platescalemm,rotang
; End hack
;;;;;;;;;;



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Account for offset of drilled holes from the fiducial positions
; using the plDrillPos and CMM file
; See http://www.apo.nmsu.edu/Telescopes/SDSS/eng.papers/19990112_PlugPlateDistortion/19990112.html
; for some background.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Read info from plDrillPos file (fiducial hole locations)
platenum=obsparam.plate
path=strcompress(getenv('MANGACORE_DIR')+'/'+getenv('MANGACORE_VER')+'/drill/',/remove_all)

drillfile=strcompress(path+'plDrillPos-'+platenum+'.par'); Fiducial positions
drillcmmfile=strcompress(path+'plate'+platenum+'_CMM.par'); Measured positions
; For the set of hole ra,dec specified read the drill files and
; calculate the offsets xdrill and ydrill from where the hole was
; meant to be.
status=mlreaddrillfiles(drillfile,drillcmmfile,objra,objdec,xdrill,ydrill)
;print,xdrill,ydrill
; If something failed, exit
if status ne 0 then mlquitmanga3d,status
; xdrill, ydrill are the shifts in xfocal, yfocal in microns of the real hole
; from the ideal hole.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Calculate the base shift of bundle center at guide wavelength due to spatial DAR
; arising from offset of the target from the drilled position,
; given the actual ra,dec of target and the drill hour angle
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Define the rough scale of the focal plane
focalscale=16.5338; arcsec/mm
; This is based on the SDSS platedesign code, which specifies
; 217.7358 mm/degree)

; Get drill hour angle from obsparam.  Might not ultimately want this
; from here as it will cause problems combining data from multiple
; plates???? Why?

 ; Convert from drill HA to effective drill LST in decimal hours
lstdrill=mlcalclstdrill(obsparam.drillha,obsparam.cenra) 

; Actual LST of observation in decimal hours
lst=mlcalclst(obsparam.LSTstart,obsparam.exptime)

; Set up 150 fake locations across plate to mimic having holes
; drilled in lots of places that the guider tries to optimize for
mlrandplateloc,obsparam.cenra,obsparam.cendec,150,rarand,decrand

; Coordinates of reference grid of random locations at drill HA
ad2xyfocal,rarand,decrand,xguide0,yguide0,racen=obsparam.cenra,deccen=obsparam.cendec,lst=lstdrill*15.
; Coordinates of reference grid at observed HA
ad2xyfocal,rarand,decrand,xguide1,yguide1,racen=obsparam.cenra,deccen=obsparam.cendec,lst=lst*15.
; Fit rotation, scale, shift parameters in guide targets
ha_fit, xguide0, yguide0, xguide1, yguide1, xnew=xtmp1, ynew=ytmp1, rot=grot, scale=gscale, xshift=gxshift, yshift=gyshift

; Drill hour angle defines the base position of the hole center
; (xfocal0,yfocal0) on the focal plane in mm
ad2xyfocal,objra,objdec,xfocal0,yfocal0,racen=obsparam.cenra,deccen=obsparam.cendec,lst=lstdrill*15.,lambda=5500.

; Actual position of object (xfocal,yfocal) on the focal plane
ad2xyfocal,objra,objdec,xfocal,yfocal,racen=obsparam.cenra,deccen=obsparam.cendec,lst=lst*15.,lambda=5500.
; Apply guider scaling, rotation, and shift at this hour angle
ha_apply, xfocal, yfocal, xnew=gxfocal, ynew=gyfocal, rot=grot, scale=gscale, xshift=gxshift, yshift=gyshift

; The difference between the hole center (where object would fall if
; observed at drill HA) and the actual object location at the observed
; HA after accounting for guide corrections sets the base
; position offset
; Define +x is West, +y is North
; Convert from mm to arcsec
xbase=-(gxfocal-xfocal0)*focalscale
ybase=(gyfocal-yfocal0)*focalscale
; Flip BOTH xbase and ybase, 'coz if the GALAXY effectively moved NE,
; that's like saying the IFU moved SW, relatively speaking.
xbase=-xbase
ybase=-ybase

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Calculate the base shift of bundle center at guide wavelength
; due to dithering
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Determine dither offset
fiberdiam=150.; Fiber outer diameter in microns
ditherrad=fiberdiam/3./platescale; Dither offset radius from zero-position in arcsec
  ; (Note that this is the distance assuming all 3 positions are offset
  ; from a central point.  Division by 3.0 is EXACT: 3.0=(2*cos(30))**2)
dangle=0.
; Dither offset (see diagram in obs strategy)
; This corresponds to moving SW for posn 1, E for posn 2, NW for posn3
if (obsparam.dposn eq 1) then dangle=240.*!PI/180.
if (obsparam.dposn eq 2) then dangle=0.*!PI/180.
if (obsparam.dposn eq 3) then dangle=120.*!PI/180.
if (obsparam.dposn eq 4) then dangle=300.*!PI/180.
if (obsparam.dposn eq 5) then dangle=180.*!PI/180.
if (obsparam.dposn eq 6) then dangle=60.*!PI/180.
if (obsparam.dposn eq 7) then dangle=0.*!PI/180.
if (obsparam.dposn eq 8) then dangle=210.*!PI/180.
if (obsparam.dposn eq 9) then dangle=270.*!PI/180.
; APOGEE dithers
if (obsparam.dposn eq 10) then dangle=240.*!PI/180.
if (obsparam.dposn eq 11) then dangle=0.*!PI/180.
if (obsparam.dposn eq 12) then dangle=120.*!PI/180.
; Coordinates where N up, E left
if (obsparam.dposn le 6) then begin
  xdither=-ditherrad*cos(dangle)
  ydither=ditherrad*sin(dangle)
endif
if (obsparam.dposn eq 7) then begin
  xdither=0.
  ydither=0.
endif
if ((obsparam.dposn eq 8)or(obsparam.dposn eq 9)) then begin
  xdither=-ditherrad*2*cos(30.*!PI/180.)*cos(dangle)
  ydither=ditherrad*2*cos(30.*!PI/180.)*sin(dangle)
endif
if (obsparam.dposn ge 10) then begin
  xdither=-ditherrad*cos(dangle)*0.694
  ydither=ditherrad*sin(dangle)*0.694
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Determine the offset due to DAR as a function of wavelength
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Actual hour angle of the observation in decimal hours (HA=LST-RA)
ha=lst-objra/15.
;print,objra,objdec,ha,lst,' ',obsparam.taistart

for j=0,nwave-1 do begin
  ; Figure out DAR and include effects of chromatic distortion
  ; in the telescope optics using optical model.
  dar=mldar(ha,objdec,wave[j],parang,xdar,ydar,alt,waveREF=5500.,TC=obsparam.tcel,RH=obsparam.rhumid,P=obsparam.pres,RAOBJ=objra,RACEN=obsparam.cenra,DECCEN=obsparam.cendec,/distort)
  ; Now the xdar and ydar numbers that we got were in xfocal, yfocal coordinates
  ; Flip the xcoordinate so that positive to the West
  xdar=-xdar
 ; if ((j eq 70)or(j eq 900)or(j eq 1900)or(j eq 2900)or(j eq 3900)or(j eq 4502)) then print,wave[j],j,-xdar,ydar
  ; Assign total offsets as a function of fiber number and wavelength
  ; using all of the above calculated numbers.
  for k=0,nfiber-1 do begin
    xpos[j,k]=-xinbundle[k]/platescalemm+xdither+xbase-xdar+xdrill/platescale
    ypos[j,k]=yinbundle[k]/platescalemm+ydither+ybase-ydar+ydrill/platescale
    ;if (j eq 1000) then print,xdither+xbase-xdar+xdrill/platescale
    ;if (j eq 1000) then print,ydither+ybase-ydar+ydrill/platescale
  endfor
endfor

return,status
end
