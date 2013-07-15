; function mlcounttotalfibers
;
; Counts the number of total fibers to combine together considering
; all fibers in each exposure.  Looks up how many fibers are in
; an ifu of the given name.
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 11/20/2012
;   Last modified: 11/20/2012
;
; REVISION HISTORY:
;   v1: 20-Nov-2012  D. Law
;       First written.
function mlcounttotalfibers,names,nexp

total=0L
for i=0,nexp-1 do begin
  nfiber=mlgetbundlesize(names[i])
  if nfiber lt 0 then mlquitmanga3d,-42L
  total+=nfiber
endfor

return,total
end
