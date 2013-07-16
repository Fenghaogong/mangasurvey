;+
; function mlmatchmjd
;
; Given a list of filenames with MJD in them, find which one
; is the best match to a given MJD (require filename mjd earlier
; than or equal to target MJD)
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 01/16/2013
;   Last modified: 01/17/2013
;
; REVISION HISTORY:
;   v1: 17-JAN-2013  D. Law
;       First written.
;-
function mlmatchmjd,filenames,basestring,mjdtarget

if (size(filenames))[0] eq 0 then begin
  splog,'WARNING- No match for ',basestring,' with mjd < ',mjdtarget
  exit
endif

nnames=(size(filenames))[1]
mjd=lonarr(nnames)

basestring = file_search(basestring)+'/'  ; This ensures any double '/' in the path is removed.

for i=0,nnames-1 do begin
  index=strsplit(filenames[i],basestring,/regex)
  mjd[i]=long(strmid(filenames[i],index,5))
endfor

; MJD is a match if less than or equal to target
; MJD<50000 is automatical fail (we weren't taking data before 1982)
okmjd=where((mjd le mjdtarget) and (mjd ge 50000L))
temp=max(mjd[okmjd],maxat)

if nnames gt 1 then splog,'Multiple matches found,'
splog,'Using ',filenames[maxat]

return,filenames[maxat]
end
