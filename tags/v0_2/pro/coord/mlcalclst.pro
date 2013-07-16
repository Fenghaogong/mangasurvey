; function mlcalclst
;
; Calculates LST in decimal hours at midpoint of a given exposure, provided
; a string-formatted LST at the start and an exposure time
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 11/22/2012
;   Last modified: 11/22/2012
;
; REVISION HISTORY:
;   v1: 22-Nov-2012  D. Law
;       First written.
function mlcalclst,lststart,exptime

temp = strsplit(lststart,':',/extract)
lst=temp[0] + temp[1]/60.D + temp[2]/3600.D + exptime/2./3600.D

return,lst
end
