;+
; function mlhalimit
;
; Returns the hour angle limit for a given declination based
; on Omega calculations and the requirement Omega<0.5 for
; 1 total hour of observing.
; Values are currently hard-coded based on analysis by DRL
; as of May 31 2013.  Based on a 5th order polynomial
; fit.
;
; Required input:
;   decl: Plate declination in degrees
;
; Returned value: Maximum HA to start/finish an exposure
;   in decimal hours.
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 04/31/2013
;   Last modified: 04/31/2013
;
;   v1: 31-May-2013 D. Law
;       First written
;-

function mlhalimit,decl

; Boundary cases
if ((decl lt -10)or(decl gt 80)) then begin
halimit=0.
endif else begin

funcfit=[$
1.59349D,$
0.109658D,$
-0.00607871D,$
0.000185393D,$
-2.54646d-06,$
1.16686d-08$
]

halimit=funcfit[0]+funcfit[1]*decl+funcfit[2]*decl*decl+funcfit[3]*decl*decl*decl+funcfit[4]*decl*decl*decl*decl+funcfit[5]*decl*decl*decl*decl*decl

if(halimit lt 0.) then halimit=0.
endelse

return,halimit
end
