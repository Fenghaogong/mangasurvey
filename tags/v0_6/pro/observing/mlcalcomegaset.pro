;+
;
; function mlcalcomegaset
;
; Calculate the Omega value for a set of 3 observations at
; a given declination taken starting at a given HA, for a
; given IFU location on the plate.
;
; Required input:
;   platedec: Plate central declination in degrees
;   HAstart: Starting hour angle in decimal hours
;
; Optional input:
;   HAdrill: HA for which plate is drilled (default is 0.0)
;   theta: Position angle of IFU on plate in degrees (default=90 E of N)
;   objdist: Distance of IFU from center of plate in degrees (default
;      1.5 degrees, edge of plate)
;   lambda: Wavelength to consider in Angstroms (default 3600 Angstroms)
;   setlength: Amount of time to complete 1 set of 3 dithers (default 1 hr)
;
; Optional output:
;   baseshift: Max movement of bundle centroid away from drilled
;      position in arcsec
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 06/10/2013
;   Last modified: 06/10/2013
;
;   v1: 10-Jun-2013 D. Law
;       Implemented function based on previous non-svn Omega calculators
;-

function mlcalcomegaset,platedec,HAstart,HAdrill=HAdrill,theta=theta,objdist=objdist,lambda=lambda,setlength=setlength,baseshift=baseshift

; Default constants
; Default drill HA is 0
if (keyword_set(HAdrill)) then HAdrill=HAdrill $
else HAdrill=0.
; Default PA of IFU on plate is 90 degrees (E)
if (keyword_set(theta)) then theta=theta $
else theta=90.
; Default distance of IFU from center of plate is 1.5 degrees
if (keyword_set(objdist)) then objdist=objdist $
else objdist=1.5
; Default lambda is 3600 Angstroms
if (keyword_set(lambda)) then lambda=lambda $
else lambda=3600.
; Default observing set length is 1 hr
if (keyword_set(setlength)) then setlength=setlength $
else setlength=1.0

; RA of plate center in hours.
; Technically irrelevant, but easiest if it's big enough that
; LST tested do not have to wrap to 24-XX.
platera=9.0

; Define observation hour angles
HAmid=HAstart+setlength/2.
HAend=HAstart+setlength

; Define some LST
lstdrill=HAdrill+platera
lststart=HAstart+platera
lstmid=HAmid+platera
lstend=HAend+platera

objdec=platedec+objdist*cos(theta/!RADEG); IFU declination in degrees
objra=platera+objdist/15./cos(platedec/!RADEG)*sin(theta/!RADEG); IFU RA in hours

; Set up 150 fake locations across plate to mimic having holes
; drilled in lots of places that the guider tries to optimize for
; Fix the rngseed to aid in possible debugging of higher routines
mlrandplateloc,platera*15.,platedec,150,rarand,decrand,rngseed=10.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Coordinates at the drill time
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Coordinates of reference grid of random locations at guide wave
ad2xyfocal,rarand,decrand,xguide0,yguide0,racen=platera*15.,deccen=platedec,lst=lstdrill*15.
; Coordinates of object of interest at wave of interest
; (This is because here we're interested in change AT that wavelength,
; not wrt the guide wavelength)
ad2xyfocal,objra*15.,objdec,xfocal0,yfocal0,racen=platera*15.,deccen=platedec,lst=lstdrill*15.,lambda=lambda

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Coordinates at start of window
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Coordinates of reference grid
ad2xyfocal,rarand,decrand,xguide1,yguide1,racen=platera*15.,deccen=platedec,lst=lststart*15.
; Coordinates of object of interest
ad2xyfocal,objra*15.,objdec,xfocal1,yfocal1,racen=platera*15.,deccen=platedec,lst=lststart*15.,lambda=lambda
;; Fit rotation, scale, shift parameters in guide targets
ha_fit, xguide0, yguide0, xguide1, yguide1, xnew=xtmp1, ynew=ytmp1, rot=rot1, scale=scale1, xshift=xshift1, yshift=yshift1
; Apply guider scaling, rotation, and shift at this hour angle
ha_apply, xfocal1, yfocal1, xnew=xfocal1g, ynew=yfocal1g, rot=rot1, scale=scale1, xshift=xshift1, yshift=yshift1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Coordinates at middle of window
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Coordinates of reference grid
ad2xyfocal,rarand,decrand,xguide2,yguide2,racen=platera*15.,deccen=platedec,lst=lstmid*15.
; Coordinates of object of interest
ad2xyfocal,objra*15.,objdec,xfocal2,yfocal2,racen=platera*15.,deccen=platedec,lst=lstmid*15.,lambda=lambda
;; Fit rotation, scale, shift parameters in guide targets
ha_fit, xguide0, yguide0, xguide2, yguide2, xnew=xtmp2, ynew=ytmp2, rot=rot2, scale=scale2, xshift=xshift2, yshift=yshift2
; Apply guider scaling, rotation, and shift at this hour angle
ha_apply, xfocal2, yfocal2, xnew=xfocal2g, ynew=yfocal2g, rot=rot2, scale=scale2, xshift=xshift2, yshift=yshift2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Coordinates at end of window
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Coordinates of reference grid
ad2xyfocal,rarand,decrand,xguide3,yguide3,racen=platera*15.,deccen=platedec,lst=lstend*15.
; Coordinates of object of interest
ad2xyfocal,objra*15.,objdec,xfocal3,yfocal3,racen=platera*15.,deccen=platedec,lst=lstend*15.,lambda=lambda
;; Fit rotation, scale, shift parameters in guide targets
ha_fit, xguide0, yguide0, xguide3, yguide3, xnew=xtmp3, ynew=ytmp3, rot=rot3, scale=scale3, xshift=xshift3, yshift=yshift3
; Apply guider scaling, rotation, and shift at this hour angle
ha_apply, xfocal3, yfocal3, xnew=xfocal3g, ynew=yfocal3g, rot=rot3, scale=scale3, xshift=xshift3, yshift=yshift3

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Calculate Omega
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Approx scale is just 
scale=16.5338 ; arcsec/mm (based on SDSS platedesign code, which
; specifies 217.7358 mm/degree)

; Omega is the maximum difference between x,y positions
; This is probably between start and end of window, but not
; necessarily since object paths describe strange loops
omega12=sqrt((xfocal2g-xfocal1g)^2+(yfocal2g-yfocal1g)^2)*scale
omega13=sqrt((xfocal3g-xfocal1g)^2+(yfocal3g-yfocal1g)^2)*scale
omega23=sqrt((xfocal3g-xfocal2g)^2+(yfocal3g-yfocal2g)^2)*scale

; Pick the largest one (probably omega13)
omega=omega12>omega13>omega23

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Calculate maximum base shift
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Baseshift is the movement of bundle centroid away from drilled position
base1=sqrt((xfocal1g-xfocal0)^2+(yfocal1g-yfocal0)^2)*scale
base2=sqrt((xfocal2g-xfocal0)^2+(yfocal2g-yfocal0)^2)*scale
base3=sqrt((xfocal3g-xfocal0)^2+(yfocal3g-yfocal0)^2)*scale

baseshift=base1>base2>base3

return,omega
end
