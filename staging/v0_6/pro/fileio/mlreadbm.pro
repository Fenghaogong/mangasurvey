;+
; function mlreadbm
;
; Reads a bundlemap information (.bm) file.  You need one of these
; for each IFU bundle, but these should be in standard reference 
; files and not change once built.
;
; Returns 0 if everything is ok, returns an error code if there was a problem.
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 11/15/2012
;   Last modified: 01/18/2013
;
; REVISION HISTORY:
;   v1: 15-Nov-2012  D. Law
;       First written.
;   v1.1: 18-Jan-2013 D. Law
;       Tweaked logging, modified to Matt's metrology and mjd archiving
;-
function mlreadbm,ifuname,mjd,bmap

; Status starts off OK
status=0

; File root for metrology data
bmtarget=strcompress(getenv('MANGAROOT')+'/mangadb/metrology/',/remove_all)
; Find the most recent mjd of metrology data that is before mjd of observation
bmfile=file_search(bmtarget,'[0-9]????')
bmtarget=strcompress(mlmatchmjd(bmfile,bmtarget,mjd)+'/',/remove_all)
; Find the most recent version of that metrology data
bmfile=file_search(bmtarget+'*')
; If no versions, fail
if (size(bmfile))[0] eq 0 then begin
  splog,'FAILED TO FIND METROLOGY DATA'
  return, -10L
endif
; If more than one version, pick the latest
if (size(bmfile))[1] gt 1 then bmtarget=bmfile[(size(bmfile))[1]-1]
; Find the metrology file in this folder that matches ifuname
bmfile=file_search(bmtarget+'/'+ifuname+'*.par')
; If no match, or more than one match, fail
if (size(bmfile))[0] eq 0 then begin
  splog,'FAILED TO FIND METROLOGY DATA'
  return, -10L
endif
if (size(bmfile))[1] ne 1 then begin
  splog,'ERROR, MULTIPLE METROLOGY DATA'
  return, -10L
endif

splog,bmfile

; Read bundlemap parameter structure
bmap=yanny_readone(bmfile,'BUNDLEMAP',hdr=bhdr,errcode=status)

; If yanny_readone failed, return error code -11
if (status ne 0) then return,-11L

; Check that the bundlemap is ok
status=mlcheckbm(ifuname,bmap,bhdr)

return,status
end
