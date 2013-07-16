;+
;
; function mlcalcomegaplate
;
; Front-end to mlcalcomegaset that attempts to calculate the maximum
; Omega value experienced by the worst-located IFU on the plate.
; Quite simplistic in this version- just tests NSEW edges of plate
; (E/W are usually about the worst).
;
; Required input:
;   platedec: Plate central declination in degrees
;   HAstart: Starting hour angle in decimal hours
;
; Optional input:
;   HAdrill: HA for which plate is drilled (default is 0.0)
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

function mlcalcomegaplate,platedec,HAstart,HAdrill=HAdrill,lambda=lambda,setlength=setlength,baseshift=baseshift

; Default constants
; Default drill HA is 0
if (keyword_set(HAdrill)) then HAdrill=HAdrill $
else HAdrill=0.
; Default lambda is 3600 Angstroms
if (keyword_set(lambda)) then lambda=lambda $
else lambda=3600.
; Default observing set length is 1 hr
if (keyword_set(setlength)) then setlength=setlength $
else setlength=1.0

; Run mlcalcomegaset for IFUs at NESW edges of plate
omegaN=mlcalcomegaset(platedec,HAstart,theta=0.,objdist=1.5,HAdrill=HAdrill,lambda=lambda,setlength=setlength,baseshift=baseshiftN)
omegaE=mlcalcomegaset(platedec,HAstart,theta=90.,objdist=1.5,HAdrill=HAdrill,lambda=lambda,setlength=setlength,baseshift=baseshiftE)
omegaS=mlcalcomegaset(platedec,HAstart,theta=180.,objdist=1.5,HAdrill=HAdrill,lambda=lambda,setlength=setlength,baseshift=baseshiftS)
omegaW=mlcalcomegaset(platedec,HAstart,theta=270.,objdist=1.5,HAdrill=HAdrill,lambda=lambda,setlength=setlength,baseshift=baseshiftW)

; Pick maximum values
omega=omegaN>omegaE>omegaS>omegaW
baseshift=baseshiftN>baseshiftE>baseshiftS>baseshiftW

return,omega
end
