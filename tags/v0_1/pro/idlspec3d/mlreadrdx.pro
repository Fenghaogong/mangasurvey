; function mlreadrdx
;
; Reads a reduction setup (.rdx) file
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 11/19/2012
;   Last modified: 01/15/2013
;
; REVISION HISTORY:
;   v1: 19-Nov-2012  D. Law
;       First written.
;   v1.1: 21-Jan-2013 D. Law
;       Tweaked logging, format
function mlreadrdx,rdxfile,redx,rhdr

; Status starts off OK
status=0

; Check whether file exists, return error code -40 if not
junk=findfile(rdxfile,count=ct)
if (ct EQ 0) then return,-40L

; Read reduction info
redx=yanny_readone(rdxfile,'REDUX',hdr=rhdr)

; If yanny_readone failed, return error code -41
if (status ne 0) then return,-41L

; Check that the reduction file is ok
; Check that all 
;status=mlcheckrdx()

return,status
end
