; function mlreadexp
;
; Reads an exposure information (exp*.par) file
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 11/19/2012
;   Last modified: 01/21/2013
;
; REVISION HISTORY:
;   v1: 19-Nov-2012  D. Law
;       First written.
;   v1.1: 21-Jan-2013 D. Law
;       Tweaked logging
function mlreadexp,expfile,obsparam

; Status starts off OK
status=0

splog,'Using EXPfile: ',expfile

; Check whether file exists, return error code -30 if not
junk=findfile(expfile,count=ct)
if (ct EQ 0) then return,-30L

; Read observing parameter block
obsparam=yanny_readone(expfile,'OPARAM',hdr=ehdr)

; If yanny_readone failed, return error code -31
if (status ne 0) then return,-31L

; Check that the exposure file is ok
;status=mlcheckexp(expfile,obsparam,targets)

return,status
end
