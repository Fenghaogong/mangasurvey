; program mlhasex
;
; Converts between sexagesimal and decimal hour angle values
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
;       Tweaked trim routine to mltrim, 2 decimal places, leading zeroes

pro mlhasex,hainput,haoutput,inpsex=inpsex

; If /inpsex was specified, then treat input as sexagesimal
; and make decimal output
if keyword_set(inpsex) then begin

HA = strsplit(hainput,':',/extract)
HAsign=1.
if (STRPOS(HA[0],'-') ne -1) then HAsign=-1.
haoutput=HAsign*(0.D + abs(HA[0]) + HA[1]/60.D + HA[2]/3600.D)

; If /inpsex was not specified, then treat input as decimal
; and make sexagesimal output
endif else begin

HA=abs(hainput)
HA1=FIX(HA)
HA1b=(HA-HA1)*60.D
HA2=FIX(HA1b)
HA3=(HA1b-HA2)*60.D

; Add leading zeroes if necessary
if (HA1 lt 10.) then HA1=strcompress('0'+string(HA1),/remove_all) $
else HA1=strcompress(string(HA1),/remove_all)
if (HA2 lt 10.) then HA2=strcompress('0'+string(HA2),/remove_all) $
else HA2=strcompress(string(HA2),/remove_all)
if (HA3 lt 10.) then HA3=strcompress('0'+MLTRIM(round(HA3*100)/100.),/remove_all) $
else HA3=strcompress(MLTRIM(round(HA3*100)/100.),/remove_all)

haoutput=''
if (hainput lt 0.D) then haoutput='-'
haoutput+=HA1+':'+HA2+':'+HA3

haoutput=strcompress(haoutput,/remove_all)

endelse

return
end
