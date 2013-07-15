; function mlreadslm
;
; Reads a slitmap information file.  You need one of these
; per plate, per mjd.
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 11/13/2012
;   Last modified: 01/21/2013
;
; REVISION HISTORY:
;   v1: 13-Nov-2012  D. Law
;       First written.
;   v1.1: 15-Jan-2013 D. Law
;       Tweaked logging, changed format
function mlreadslm,plate,mjd,slitmap

; Status starts off OK
status=0

slmfile=strcompress(getenv('MANGAROOT')+'/mangadb/slitmaps/slitmap-'+string(plate)+'-'+string(mjd)+'.par',/remove_all)

splog,slmfile

; Check whether file exists, return error code -20 if not
junk=findfile(slmfile,count=ct)
if (ct EQ 0) then return,-20L

; Read first slitmap parameter structure
blocks=yanny_readone(slmfile,'BLOCKS',hdr=shdr,errcode=status)

; If read failed, return error code -21
if (status ne 0) then return,-21L

; Read second slitmap parameter structure
slitmap=yanny_readone(slmfile,'SMAP',hdr=shdr,errcode=status)

; If read failed, return error code -22
if (status ne 0) then return,-22L

; Check that the slitmap is ok
status=mlcheckslm(slitmap,blocks,shdr)

return,status
end
