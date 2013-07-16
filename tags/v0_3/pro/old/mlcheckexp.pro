; function mlcheckexp
;
; Purpose is to check an exposure (.exp) file to make sure that it's ok.
; I.e., 
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 11/23/2012
;   Last modified: 01/15/2015
;
; REVISION HISTORY:
;   v1: 23-Nov-2012  D. Law
;       First written.
;   v1.1: 15-Jan-2013 D. Law
;       Tweaked logging
function mlcheckexp,expfile,obsparam,targets

; Status starts off OK
status=0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Check that only 1 entry in OBSPARAM
; If fails, error code -34
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if (size(obsparam))[1] ne 1 then begin
  splog,strcompress('More than one entry for OPARAM in file '+expfile)
  status=-34L
  return,status
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Check that dither position is either 1, 2, or 3
; If fails, error code -35
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

dith=obsparam[0].dposn
if ((dith ne 1)and(dith ne 2)and(dith ne 3)) then begin
  splog,strcompress('Invalid entry ('+string(dith)+') for dither position (dposn) in file '+expfile)
  status=-35L
  return,status
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Check that no duplicated ifu or galaxy names in BUNDLES
; If fails, error code -36 or -37
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

nbund=(size(targets))[1]
for i=0,nbund-2 do begin
  for j=i+1,nbund-1 do begin
    if targets[i].ifuname eq targets[j].ifuname then begin
      splog,strcompress('Duplicate of ifuname '+targets[i].ifuname+' in exposure file '+expfile)
      status=-36L
      return,status
    endif
    if targets[i].objname eq targets[j].objname then begin
      splog,strcompress('Duplicate of objname '+targets[i].ifuname+' in exposure file '+expfile)
      status=-37L
      return,status
    endif
  endfor
endfor

return,status
end
