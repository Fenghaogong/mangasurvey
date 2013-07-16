; function mlfluxcor
;
; Placeholder function for flux correction
;
; Returns 0 if everything ok, returns an error code if there was a problem.
;
; Author:
;   David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   First written: 01/23/2013
;   Last modified: 01/23/2013
;
; REVISION HISTORY:
;   v1: 23-Jan-2013  D. Law
;       Placeholder function written

function mlfluxcor, Cfile1


; Write out sky subtracted files
; In placeholder function, simply copy them to
; new names
Ffile1=str_replace(Cfile1,'spCFrame','spFFrame')
spawn,strcompress('cp '+Cfile1+' '+Ffile1)

return,0
end
