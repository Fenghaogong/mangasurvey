; function mlcalclstdrill
;
; Calculates effective LST in decimal hours for the drill hour angle
; given the plate central RA.  cenra is in decimal degrees, hadrill
; is in sexagesimal hours.
; Assume + hour angle is West of meridian (HA=LST-RA)
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 11/22/2012
;   Last modified: 11/22/2012
;
; REVISION HISTORY:
;   v1: 22-Nov-2012  D. Law
;       First written.
function mlcalclstdrill,hadrill,cenra

temp = strsplit(hadrill,':',/extract)

tsign=1.
if (STRPOS(hadrill,'-') ne -1) then tsign=-1.

lstdrill=tsign*(abs(temp[0]) + temp[1]/60.D + temp[2]/3600.D) + cenra/15.

if (lstdrill lt 0.) then lstdrill=lstdrill+24.
if (lstdrill gt 24.0) then lstdrill=lstdrill-24.

return,lstdrill
end
