; function mlsetwcalib
;
; Set the master wavelength grid from calibration file
;
; Returns 0 if everything is ok, returns an error code if there was a problem.
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 11/21/2012
;   Last modified: 01/17/2013
;
; REVISION HISTORY:
;   v1: 21-Nov-2012  D. Law
;       First written.
;   v1.1: 17-Jan-2013 D. Law
;       Tweaked logging, filepath
function mlsetwcalib,wave

; Status starts off OK
status=0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Calibration file listing standard wavelength grid
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

calibfile=strcompress(getenv('MANGAROOT')+'/mangadb/bosscal/56280/calibmatrix.fits',/remove_all)
splog,calibfile

; Check whether calibration file exists, return error code -120 if not
junk=findfile(calibfile,count=ct)
if (ct EQ 0) then return,-120L

calibmatrix=readfits(calibfile)
nwave=(size(calibmatrix))[2]
; wave is the master wavelength vector for the analysis, log spacing
; from 3500 to 10000 Angstroms
wave=dblarr(nwave)
wave[*]=calibmatrix[0,*]

return,status
end
