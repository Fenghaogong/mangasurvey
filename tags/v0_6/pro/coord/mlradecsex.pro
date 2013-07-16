; program mlradecsex
;
; Converts between sexagesimal and decimal RA/declination values
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 10/03/2012
;   Last modified: 12/30/2012
;
; REVISION HISTORY:
;   v1: 03-Oct-2012  D. Law
;       Adapted old code to manga pipeline edition
;   v1.1: 30-DEC-2012 D. Law
;       Tweaked trim routine to mltrim, 2 decimal places

pro mlradecsex,rainput,decinput,raoutput,decoutput,inpsex=inpsex

; If /inpsex was specified, then treat input as sexagesimal
; and make decimal output
if keyword_set(inpsex) then begin

RA = strsplit(rainput,':',/extract)
raoutput=0.D + RA[0] + RA[1]/60.D + RA[2]/3600.D
raoutput=raoutput*15.D

DEC=strsplit(decinput,':',/extract)
DECsign=1.
if (STRPOS(DEC[0],'-') ne -1) then DECsign=-1.
decoutput=DECsign*(0.D +abs(DEC[0]) + DEC[1]/60.D + DEC[2]/3600.D)

; If /inpsex was not specified, then treat input as decimal
; and make sexagesimal output
endif else begin

RA=rainput/15.D
RA1=FIX(RA)
RA1b=(RA-RA1)*60.D
RA2=FIX(RA1b)
RA3=(RA1b-RA2)*60.D

DEC=abs(decinput)
DEC1=FIX(DEC)
DEC1b=(DEC-DEC1)*60.D
DEC2=FIX(DEC1b)
DEC3=(DEC1b-DEC2)*60.D

raoutput=string(RA1)+':'+string(RA2)+':'+MLTRIM(round(RA3*100)/100.)

decoutput=''
if (decinput lt 0.D) then decoutput='-'
decoutput+=string(DEC1)+':'+string(DEC2)+':'+MLTRIM(round(DEC3*100)/100.)

raoutput=strcompress(raoutput,/remove_all)
decoutput=strcompress(decoutput,/remove_all)

endelse

return
end
